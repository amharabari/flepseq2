#!/usr/bin/env python3



import os
import csv
import pysam
import edlib
import logging
import argparse
import regex
import re
import gffutils
import pandas as pd

from tqdm import tqdm
from pathlib import Path
from itertools import groupby
from typing import Iterator, Tuple
from collections import defaultdict
from rich.logging import RichHandler

__description__ = "Process nanopore reads from a minimap2 alignment file (.bam, .sam or .cram)"
__author__      = "David Pflieger _ Helene Zuber"
__email__       = "david.pflieger@ibmp-cnrs.unistra.fr"
__credits__     = [""]
__license__     = "MIT"
__version__     = "3.0.0" # major, minor, patch
__maintainer__  = ["David Pflieger Helene Zuber"]
__status__      = "Development"

_SCRIPT_DIR = Path(__file__).resolve().parent

# DORADO ALIGNEMENT TAGS DESCRIPTION
# https://github.com/nanoporetech/dorado/blob/release-v0.8/documentation/SAM.md

def set_logger(args) -> logging.Logger:
    """
    Set up a logger to print to stdout. By default, logging level is set to INFO.
    If --debug is used, it prints every level: ERROR, WARNING, INFO & DEBUG.
    If --color is used, the output is colorized using RichHandler.
    """
    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO, 
        format= "%(message)s", #if args.color else "%(asctime)s - %(message)s", 
        datefmt="[%X]", 
        handlers=[RichHandler()] if args.color else [logging.StreamHandler()]
    ) 
    return(logging.getLogger(__name__))

def arguments():
    # Definition of arguments using argparse
    parser = argparse.ArgumentParser(description="Process Illumina paired-end alignment files")
    parser.add_argument("input", type=Path, help="Input alignment file (.bam, .sam, .cram)")
    parser.add_argument("-a", "--annotation", type = str, required=True, help="Input annotation file (must be in GFF or GTF format). \
                        If GTF format, it should contains only exon information, transcript will be deduced automatically, use then 'transcript' as feature.")
    parser.add_argument("-g", "--genome", required=True, type = str, help = "Reference genome file in fasta format.")
    parser.add_argument("-o", "--output", required=True, type=Path, help= "Output CSV file for results.")
    parser.add_argument("-f", "--feature", type=str, nargs="+", default="mRNA", help="Feature(s) to extract from the annotation file. Example: '-f exon transcript'. \
                        Use 'transcript' if annotation contains only exons. Default: ['mRNA'].")
    parser.add_argument("-n", "--name", type=str, default="gene_id", help="Annotation attribute for gene/transcript IDs (GFF/GTF). Default: 'gene_id'.")
    parser.add_argument("-m", "--min_softclip", type=int, default=3, help="Minimum softclip length (bp) to consider. Default: 3bp (inclusive).")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of CPU threads to use. Default is 4.") # TODO: use joblib or multiprocess to speed up thing
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s (v{__version__})")
    parser.add_argument("--progress", action="store_true", help="Enable progress bar estimation.")
    parser.add_argument("--color", action="store_true", help="Enable color mode.")
    parser.add_argument("--debug", action="store_true", help="Enable debugging mode.")
    return parser.parse_args()

def load_annotations(args, aln_chromosomes: list[str]) -> list[gffutils.Feature]:
    log.info(f"Loading annotation file...")
    # Gffutils create a database of the annotation file
    annot_db_path = args.annotation + ".db" # We just add .db to the given annotation file
    # if "transcript" in args.feature:
    #     args.name = ["transcript_id"] + args.name 
    if not os.path.exists(annot_db_path):
        # Create annotation database if it does not exist
        log.info(f"Creating database {annot_db_path}...")
        gff_db = gffutils.create_db(args.annotation, dbfn=annot_db_path, keep_order=True, gtf_gene_key = args.name, verbose = True, merge_strategy="create_unique")

    gff_db = gffutils.FeatureDB(annot_db_path, keep_order=True)
 
    # Retrieve only lines matching the alignment's chromosomes and user's input features
    features_dict = defaultdict(list)
    features_list = list(gff_db.features_of_type(args.feature, order_by='start'))

    # Put features in a dict indexed by chromosome, to speed up search
    for chr in aln_chromosomes:
        for feature in features_list:
            if feature.seqid == chr:
                feature.sequence = feature.sequence(args.genome, use_strand=True)
                features_dict[chr].append(feature)
                #log.debug(feature)
    
    log.info("Estimating homopolymers content...")
    log.debug(features_list[0:10])
    
    if len(features_list) == 0: 
        log.info(f"No feature detected with {args.feature}. Exiting now...")
        exit()
    # DataFrame with features for chromosomes/seqnames detected in the alignement header
    features_df = pd.DataFrame([
        {
            'Chr': f.seqid, 
            'Name':  f.attributes.get(args.name, default='transcript_id')[0], 
            'Start': f.start, 
            'End': f.end, 
            'Strand': f.strand, 
            'FeatureType': f.featuretype, 
            'Sequence': f.sequence[:10] + "..." + f.sequence[-10:], 
            **homopolymer_estimation(f.sequence)
        }  for _, features in features_dict.items() for f in features]
        ).sort_values("HP_content", ascending=False)

    #features_df.to_csv("Features.stats.csv", index=False, sep='\t')
    #log.info(f"File Features.stats.csv created!")
    log.debug(features_df.head)
    log.info(f"Number of feature(s) found: \n\t {features_df[['Chr', 'FeatureType']].value_counts().to_string()} \n")
    return features_dict, gff_db #$

def is_aln_csorted(aln: pysam.AlignmentFile) -> bool:
    """
    Determine if the alignment file is sorted by coordinates or query name.
    Returns True if the alignment file is sorted by coordinates or False if sorted by query name.
    Raises an error, if the alignment file is not sorted by either coordinates or query name. 
    """
    log.debug(f"Alignment header: \n{aln.header}\n")
    sort_order = aln.header['HD']['SO']

    if sort_order == "coordinate":
        log.info("Alignment file is sorted by coordinates.")
        return True
    elif sort_order == "queryname":
        log.info("Alignment file is sorted by query name.")
        return False
    else:
        raise argparse.ArgumentTypeError(f"⚠️  Error! The alignment file is not sorted by coordinate or query name: {sort_order}")

def reverse_complement(seq: str):
    if seq:
        pairing = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        return "".join(pairing[nucl] for nucl in reversed(seq))
    else:
        return ""

def find_feature(read: pysam.AlignedSegment, features: list[gffutils.Feature], FLANKING_BP: int = 150, WIDTH_THRESHOLD = 0.3) -> list[gffutils.Feature, gffutils.Feature, str]:
    """
    Return a read annotation if any found. The read start position should be between the feature start and end.
    If the read start is inside a feature but match in the end of the feature, we skip this feature to deal overlapping annotations. 
    Threshold is set around 2/3 of the feature width to deal with overlaping features
    """

    read_feature = None # Set default to None

    # Go through the genomic features 
    for feature in features:
        # Skip search if read is already annotated
        if read_feature: break

        # if the read start AND the read end is inside an annotation, no doubt
        if feature.start - FLANKING_BP <= read.reference_start and read.reference_end <= feature.end + FLANKING_BP:
            read_feature = feature
            # # Check feature strand and get the distance from the read end to the feature 3' end 
            # if feature == "-":
            #     dist_from_3prime = read.reference_start - feature.start
            # else: 
            #     dist_from_3prime = feature.end - read.reference_end

            # # Get allowed distance from 3' end to determine annotation. 
            # feature_width_threshold = (feature.end - feature.start) * WIDTH_THRESHOLD
            
            # print(read.reference_start, read.reference_end, feature.start, feature.end)
            # print(feature.attributes.get("ID")[0], feature.strand, (feature.end - feature.start), feature_width_threshold, abs(dist_from_3prime))
            
            # # If end of read is around the 3' end of the feature
            # if abs(dist_from_3prime) < feature_width_threshold:
            #     print("OK")
            #     read_feature = feature
    return read_feature

def extract_softclips(read: pysam.AlignedSegment, feature: gffutils.Feature) -> list[str]:
    """
    Extract softclip sequences from a read CIGAR
    We extract the FIRST or LAST softclip sequence from the CIGAR code.
    The softclip sequence is reverse complemented when the feature is on reverse strand
    """

    # Store left and right softclips sequences, if any
    # Left should be adapter and right should be the tail
    softclips = {"left": None, "right": None}
 
    # To get the first and last element in the cigar tuple: 
    # https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples

    # Last value of CIGAR tuple should be "4" for SOFTCLIP
    # Extract left part of cigar code 5S145M => 5S
    if read.cigartuples[0][0] == 4:
        softclips["left"] = read.query_sequence[:read.cigartuples[0][1]]
    # Extract right part of CIGAR: 145M5S => 5S
    if read.cigartuples[-1][0] == 4:
        softclips["right"] = read.query_sequence[-read.cigartuples[-1][1]:]

    # If strand -, tail is left softclip and if strand +, tail is right softclip
    if feature.strand == "-":
        softclips["left"], softclips["right"] = reverse_complement(softclips["right"]), reverse_complement(softclips["left"])

    return softclips.values()

def trim_adapter(seq: str, adapter: str):
    """ Align and trim an adapter sequence from sequence using edlib alignment infix mode """
    assert len(seq) > 0, "Cannot trim empty sequence!"
    try:
        result = edlib.align(adapter, seq, task="path", mode='HW', k = 4) # HW = infix mode
        start_adapt, end_adapt = result['locations'][0]
        trim_seq, adapter_seq= seq[:start_adapt], seq[start_adapt:]
        # nice = edlib.getNiceAlignment(result, adapter, seq)
        # print(result)
        # print("\n".join(nice.values()))
        # print(f"{seq} {trim_seq} {adapter_seq}")
        # print("----" *10)
        return [trim_seq, adapter_seq]
    except:
        return ["", ""]

def get_composition(seq: str) -> list: # modified by HZ
    """
    Get the composition and the percentage of ATGC in a DNA sequence. 
    Ignore all nucleotides other than ATGC (like N).
    """
    nucl_count = {'A': 0, 'T': 0, 'G': 0, 'C': 0} #mod HZ
    nucl_perc = {'%A': 0, '%T': 0, '%G': 0, '%C': 0} #mod HZ
    tail_tag = {'tag':"none"}
    if seq:
        for nucl in seq:
            if nucl in nucl_count.keys():
                nucl_count[nucl] += 1
    total_nucl = sum(nucl_count.values())
    if total_nucl > 0:
        tail_tag['tag'] = "other"
        for nucl, count in nucl_count.items():
            perc = round(count * 100.0 / total_nucl, 2)
            nucl_perc[f"%{nucl}"] = perc
            if perc > 70:
                tail_tag['tag'] = f"{nucl}-tail"
    return list(nucl_count.values()) + list(nucl_perc.values()) + list(tail_tag.values())

def run_length_encoding(seq: str) -> Iterator[Tuple[str, int]]:
    """
    Returns run length encoded Tuples for string. 
    See https://stackoverflow.com/questions/18948382/run-length-encoding-in-python
    A memory efficient (lazy) and pythonic solution using generators
    """
    return ((x, sum(1 for _ in y)) for x, y in groupby(seq))

def homopolymer_estimation(sequence: str, min_homopolymer_length: int = 5) -> dict:
    """
    Estimate the level of homopolymer in a sequence and its last 50 bp
    """

    rle_full = run_length_encoding(sequence)
    rle_last_50 = run_length_encoding(sequence[-50:])
    
    # Filter to find homopolymers (runs >= threshold)
    homopolymers_full = [(base, length) for base, length in rle_full if length >= min_homopolymer_length]
    homopolymers_last_50 = [(base, length) for base, length in rle_last_50 if length >= min_homopolymer_length]

    # Calculate percentage of homopolymer content
    total_hp_bases = sum(length for _, length in homopolymers_full)
    last_50_hp_bases = sum(length for _, length in homopolymers_last_50)

    total_content_percentage = round((total_hp_bases / len(sequence)) * 100, 2) if sequence else 0
    last_50_content_percentage = round((last_50_hp_bases / min(len(sequence), 50)) * 100, 2)

    return {
        "HP_content": total_content_percentage,
        "50bp_HP_content": last_50_content_percentage,
        "Total_HP": homopolymers_full,
        "50bp_HP": homopolymers_last_50
    }

def realign_softclip(softclip_seq: str, feature_seq: str):
    """
    Use edlib to align softclip sequence to feature sequence
    """

    len_softclip = len(softclip_seq)
    
    # Perform the alignment with edlib
    # Modes are: global(NW), prefix(SHW) and infix(HW)
    # Infix mode is best suited here
    result = edlib.align(softclip_seq, feature_seq, task="path", mode = "HW")  # 'path' gives the CIGAR string

    # print("\n----------------------------")
    # print("Softclip sequence:", softclip_seq)
    # # print("Feature:", feature_seq)
    # print(result)
    # nice = edlib.getNiceAlignment(result, softclip_seq, feature_seq)
    # print("\n".join(nice.values()))
    # print("---------------------------------------")

    # If it aligns, it is not a real softclip
    # result["locations"] --> number of gaps. Here we limit to two since the softclips are small compared to 
    # the features they will be aligned to.
    if result["editDistance"] <= get_max_edits(len_softclip) and len(result["locations"]) <= 2: 
        # print("Found feature overlap")
        return ""
    else:
        # print("No overlap with feature")
        return(softclip_seq)

def realign_tail(softclip_seq: str, feature_seq: str, max_length: int = 30):
    len_softclip = len(softclip_seq)
    # Extract the last five nucleotides 
    softclip_split = softclip_seq.split(feature_seq[-5:], 1)
    softclip_split_1bp = softclip_seq.split(feature_seq[-5:-1], 1)
    # If the feature sequence is exactly found, we just trim it
    if len(softclip_split) > 1:
        tail = softclip_split[1]
        return(tail)
    if len(softclip_split_1bp) > 1:
        tail = softclip_split_1bp[1]
        return(tail)
    
    # Otherwise we match sequence and softclip using regex to allow for some mismatches/errors
    max_errors = get_max_errors(len_softclip)
    pattern = regex.compile(f"({feature_seq}){{e<={max_errors}}}")
    match = pattern.search(softclip_seq)

    # print(f"-----------------------")
    # print(f"Reference: {feature_seq}")
    # print(f"Softclip: {softclip_seq}, {len_softclip} bp")
    
    if match:
        tail = softclip_seq[match.end():]
        # print("Pattern", pattern)
        # print("Match found at position:", match.span())
        # print("Matched sequence:", match.group(1))
    else: 
        tail = softclip_seq

    # print(f"Extracted tail: {tail}")
    # print(f"-----------------------")

    return tail

def robust_get_tag(read: pysam.AlignedSegment, tag_name: str, default_value = None):
    """ 
    Retrieve the specified tag from an AlignedSegment read returning a default value if tag is not present
    """
    try:  
        return read.get_tag(tag_name)
    except KeyError:
        return default_value

def library_size(args):
    log.info('Retrieving flagstats info')
    flag = pysam.flagstat(f"{args.input}", f"-@ {args.threads}")
    count = int(flag.split()[0])
    return count

def flagstat(args):
    output = pysam.flagstat(str(args.input))
    lines = output.split('\n')
    values = [int(line.split()[0]) for line in lines if line]
    index = [
        'total', 'secondary', 'supplementary', 'duplicates', 'mapped',
        'paired', 'read1', 'read2', 'properly_paired', 'both_mapped',
        'singletons', 'mate_diff_chr', 'mate_diff_chr_mapq5'
    ]

    return pd.DataFrame({'count': values}, index = index)

def isolate_addtail(tail, stretch_size=6): # modified by HZ
    add_tail = ""
    a_stretch_found = False
    a_stretch = ""  # variable pour stocker le stretch d'As
    rle = run_length_encoding(tail)

    for nucl, count in rle:
        if a_stretch_found:
            # Après avoir trouvé le stretch, on récupère la queue
            add_tail += nucl * count
        else:
            # On cherche le stretch d'As
            if nucl == 'A' and count > stretch_size:
                a_stretch_found = True
                a_stretch = nucl * count  # stocke le stretch d'As

    return a_stretch, add_tail

def get_RA5(seq: str, adapter: str):
    """ Align and extract the RA5 adapter sequence from sequence using edlib alignment infix mode """
    assert len(seq) > 0, "Cannot trim empty sequence!"
    try:
        result = edlib.align(adapter, seq, task="path", mode='HW', k = 4) # HW = infix mode
        start_adapt, end_adapt = result['locations'][0]
        adapter_seq= seq[start_adapt:]
        return(adapter_seq)
    except:
        return("")


if __name__ == "__main__":

    # Set args, see https://docs.python.org/3/library/argparse.html
    args = arguments()

    # Set log, see https://docs.python.org/3/library/logging.html
    log = set_logger(args)

    # Output file header
    FIELDNAMES = ["Read.name", 
                  "Read.chromosome", "Read.start", "Read.end", "Read.flag", #"Read.cigar", 
                  "Feature.name", "Feature.chromosome", "Feature.start", "Feature.end", "Feature.strand",
                  "CDS.start_codon", "CDS.stop_codon", "CDS_dist_5prime", "CDS_dist_3prime", "Degradation_Tag",#$
                  "Left_softclip", "RA5", "RA5_tag", "Right_trim_softclip", "Right_adapter", "PolyA","Add_tail", "PolyA_length_basecall", "PolyA_length_signal", 
                  "Add_tail_length", "Dist_from_5prime", "Dist_from_3prime", "A", "T", "G", "C", "%A", "%T", "%G", "%C", "Tail_tag"
                  ]
    
    # 5' to 3' adapter _ 3' adapter
    ADAPTER = "CTGTAGGCACCATCAAT"

    # 5' to 3' adapter _ 5' adapter
    ADAPTER_RA5 = "CTTTCCCTACACGACGCTCTTCCGATCT" # remy
    #ADAPTER_RA5= "GTTCAGAGTTCTACAGTCCGA" # CG après primer GUUCAGAGUUCUACAGUCCGACGAUC
    #total_reads = flagstats(args)

    # Go through an alignement file: https://pysam.readthedocs.io/en/latest/api.html
    with pysam.AlignmentFile(args.input, threads = args.threads) as aln, open(args.output, 'w') as csvfile: 
        log.info(f"Parsing alignment file: {args.input.name}")

        # Check if the aln is sorted by coordinate or queryname 
        is_csorted = is_aln_csorted(aln) 

        # Retrieve aln chromosomes from header
        aln_chromosomes = [aln['SN'] for aln in aln.header['SQ']]
        log.info(f"Chromosomes/contigs identified in alignment: {aln_chromosomes}")

        # Load annotations for chromosomes/contigs that matches the alignment file
        # https://daler.github.io/gffutils/index.html
        #features = load_annotations(args, aln_chromosomes)
        features, gff_db = load_annotations(args, aln_chromosomes) #$
        
        # Init the csv writer
        writer = csv.writer(csvfile, delimiter='\t')
        writer.writerow(FIELDNAMES)

        # Init variables
        results, reads_skipped = [], []
        reads_without_annotation, reads_supplementary, reads_unmapped, reads_secondary, final_reads, final_reads_RA5 = 0, 0, 0, 0, 0, 0
        #adapter, tail, tail_adapter, add_tail, a_stretch = "", "", "", "", ""
        

        # Go through each read
        for read in tqdm(aln, desc="Processing", unit=" reads", disable=None, total = library_size(args) if args.progress else None):
            RA5            = ""
            RA5_tag        = "NO_RA5"
            tail           = ""
            tail_adapter   = ""
            add_tail       = ""
            a_stretch      = ""
            softclip_left  = None
            softclip_right = None
            polya_len      = None
            deg_tag        = "No_CDS"
            # Skip reads that are not primary aligned. There is no "is_primary" function in pysam
            # UNMAPPED: Reads that did not map but still present in alignement file
            if read.is_unmapped:
                reads_unmapped += 1
                continue
            # SUPPLEMENTARY: Represent parts of a read that map to different locations on the genome, typically in the context of split or chimeric alignments. 
            # For example, when a read spans a large structural variation or a complex genomic rearrangement, the alignment may be split into multiple segments.
            if read.is_supplementary:
                reads_supplementary += 1 
                continue
            # SECONDARY: Represent alternative alignments for a read that did not score as the primary (best) alignment. 
            # These alignments are typically produced when a read has multiple possible mapping locations, 
            # and one location is chosen as the primary while others are marked as secondary
            if read.is_secondary:
                reads_secondary += 1
                continue

            # Find read annotation/feature
            read_feature = find_feature(read, features[read.reference_name])
            
            # If read has an annotation
            if read_feature:
                final_reads += 1
                # Get dorado polyA length estimation from the "pt" tag in the alignment
                polya_len = robust_get_tag(read, "pt")
                # Extract the softclips sequences by parsing the CIGAR code
                softclip_left, softclip_right = extract_softclips(read, read_feature)
                # If right softclip is present => tail
                if softclip_right:
                    # Check for adapter and separate tail and adapter 
                    tail, tail_adapter = trim_adapter(softclip_right, adapter = ADAPTER)
                    # Check for additional nucleotide in the end of the polyA tail
                    a_stretch, add_tail = isolate_addtail(tail)
                    # HZ addition: extract final 'T' if present
                    match = re.search(r"(T+)$", tail)
                    add_tail = match.group(1) if match else add_tail
                # If left softclip_left => look for RA5 _ HZ modification
                if softclip_left:
                    RA5 = get_RA5(softclip_left, adapter = ADAPTER_RA5)
                    if len(RA5) > 0:
                        RA5_tag = "RA5"
                        final_reads_RA5 += 1
                    else:
                        RA5_tag = "NO_RA5"

                #print(read.query_name, read.reference_name, tail, tail_adapter, add_tail, polya_len, read_feature.strand, read_feature.start, read_feature.end)
        
                if read_feature.strand == "-":
                    dist_from_3prime = read.reference_start - read_feature.start
                    dist_from_5prime = read.reference_end - read_feature.end # HZ addition
                else: 
                    dist_from_3prime = read_feature.end - read.reference_end
                    dist_from_5prime = read_feature.start - read.reference_start # HZ addition
                
                cds_parts = list(gff_db.children(read_feature, featuretype="CDS", order_by="start"))

                if cds_parts:
                    #$ Physical boundaries if the total CDS length
                    cds_min = cds_parts[0].start
                    cds_max = cds_parts[-1].end

                    if read_feature.strand == "-":
                        start_codon_pos = cds_max
                        stop_codon_pos = cds_min

                        #$ 5' is high coord (end) and 3' is low coord (start)
                        cds_dist_5prime = start_codon_pos - read.reference_end     
                        cds_dist_3prime = read.reference_start - stop_codon_pos
                    else:
                        start_codon_pos = cds_min
                        stop_codon_pos = cds_max

                        #$ 5' is low coord (start) and 3' is high coord (end)
                        cds_dist_5prime = read.reference_start - start_codon_pos
                        cds_dist_3prime = stop_codon_pos - read.reference_end
                    
                    #$ Degradation tagging, using a 3bp buffer

                    is_5_in_cds = cds_dist_5prime >= -10
                    is_3_in_cds = cds_dist_3prime >= -3

                    is_5_outside_cds = cds_dist_5prime < -10
                    is_3_outside_cds = cds_dist_3prime < -3

                    if is_5_outside_cds and is_3_outside_cds:
                        deg_tag = "Intact"
                    elif is_5_in_cds and is_3_outside_cds:
                        deg_tag = "5'_truncated"
                    elif is_5_outside_cds and is_3_in_cds:
                        deg_tag = "3'_truncated"
                    elif is_5_in_cds and is_3_in_cds:
                        deg_tag = "Both_truncated"
                else:
                    start_codon_pos = "NA"
                    stop_codon_pos = "NA"
                    cds_dist_5prime = "NA"
                    cds_dist_3prime = "NA"
                    deg_tag = "No_CDS"

                infos = [read.query_name, 
                        read.reference_name, read.reference_start, read.reference_end, read.flag, #read.cigarstring, 
                        read_feature.attributes.get(args.name, default='transcript_id')[0], read_feature.seqid, read_feature.start, read_feature.end, read_feature.strand,
                        start_codon_pos, stop_codon_pos, cds_dist_5prime, cds_dist_3prime, deg_tag,#$
                        softclip_left, RA5, RA5_tag, tail, tail_adapter, a_stretch, add_tail, len(a_stretch) if a_stretch else 0, polya_len, len(add_tail) if add_tail else 0,
                        dist_from_5prime, dist_from_3prime, *get_composition(add_tail)
                        ]

                writer.writerow(infos)

            else: 
                reads_without_annotation += 1 

        # Convert the list to DataFrame 
        log.info(f"Reads unmapped: {reads_unmapped }")
        log.info(f"Reads secondary: {reads_secondary}")
        log.info(f"Reads supplementary: {reads_supplementary}")
        #log.info(f"Total reads skipped: {reads_unmapped + reads_secondary + reads_supplementary}")
        log.info(f"Reads not overlapping annotation: {reads_without_annotation}")
        log.info(f"Final number of reads: {final_reads}")
        log.info(f"Final number of reads containing RA5 sequence: {final_reads_RA5}")
        log.info(f"DONE! 👍")
