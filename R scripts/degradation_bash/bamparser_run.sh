#!/bin/bash
#SBATCH -p fast                              # partition
#SBATCH -n 9                                 # number of parallel tasks (for 9 barcodes: 07–15)
#SBATCH -c 2                                 # CPUs per task
#SBATCH --mem-per-cpu=4000                   # 4 GB per CPU
#SBATCH -o bamparser.%N.%j.out.log           # stdout log
#SBATCH -e bamparser.%N.%j.err.log           # stderr log

# Paths
BAM_DIR="/shared/home/aharabari/Projects/FLEP-seq2/1_Mapping_Flepseq2"
OUT_DIR="/shared/home/aharabari/Projects/FLEP-seq2/2_Output"
SCRIPT="/shared/home/aharabari/Projects/FLEP-seq2/scripts/bamparser_v4.0.py"
GFF="/shared/home/aharabari/Projects/inputs/Araport11_CLEAN.gff"
FASTA="/shared/home/aharabari/Projects/inputs/TAIR10_chr_all.fas"
GFF_DB="${GFF}.db"

# Activate environment once in the job
eval "$(micromamba shell hook --shell bash)"
micromamba activate flepseq

# 1. DATABASE CHECK (single-threaded)
if [ ! -f "$GFF_DB" ]; then
    echo "Library database not found. Building it now..."
    FIRST_BAM=$(ls "$BAM_DIR"/*.bam | head -n 1)
    python3 "$SCRIPT" "$FIRST_BAM" \
        -a "$GFF" -g "$FASTA" \
        -o /tmp/init.tsv \
        -f mRNA -n ID \
        --threads 1
else
    echo "Library database found! Skipping initialization and starting workers..."
fi

echo "Starting parallel processing with srun..."

# 2. Loop over desired barcodes and launch one srun per BAM
for bam in "$BAM_DIR"/*_barcode0[7-9].bam "$BAM_DIR"/*_barcode1[0-5].bam
do
    # If the glob does not match anything, skip
    [ -e "$bam" ] || continue

    name=$(basename "$bam" .bam)

    srun --exclusive -N1 -n1 \
        python3 "$SCRIPT" "$bam" \
            -a "$GFF" \
            -g "$FASTA" \
            -o "$OUT_DIR/${name}_bamparsed.tsv" \
            -f mRNA \
            -n ID \
            --threads "$SLURM_CPUS_PER_TASK" \
            --progress \
            --color &
done

wait
echo "All files processed!"