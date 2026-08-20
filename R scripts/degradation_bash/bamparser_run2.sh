#!/usr/bin/env bash
# Define paths
BAM_DIR="/shared/home/aharabari/Projects/FLEP-seq2/1_Mapping_Flepseq2"
OUT_DIR="/shared/home/aharabari/Projects/FLEP-seq2/2_Output"
SCRIPT="/shared/home/aharabari/Projects/FLEP-seq2/scripts/bamparser_v5.0.py"
GFF="/shared/home/aharabari/Projects/inputs/Araport11_CLEAN.gff"
FASTA="/shared/home/aharabari/Projects/inputs/TAIR10_chr_all.fas"

# The path to the database the script creates
GFF_DB="${GFF}.db"

# 1. DATABASE CHECK: 
if [ ! -f "$GFF_DB" ]; then
    echo "Library database not found. Building it now..."
    FIRST_BAM=$(ls "$BAM_DIR"/*.bam | head -n 1)
    python3 "$SCRIPT" "$FIRST_BAM" -a "$GFF" -g "$FASTA" -o /tmp/init.tsv -f mRNA -n ID --threads 1
else
    echo "Library database found! Skipping initialization and starting workers..."
fi
# 2. PARALLEL PROCESSING:
echo "Starting parallel processing with xargs..."
ls "$BAM_DIR"/*.bam | xargs -I {} -P 4 bash -c '
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate flepseq
    file="{}"
    name=$(basename "$file" .bam)
    python3 "'$SCRIPT'" "$file" \
        -a "'$GFF'" \
        -g "'$FASTA'" \
        -o "'$OUT_DIR'/${name}_bamparsed.tsv" \
        -f mRNA \
        -n ID \
        --threads 2 \
        --progress \
        --color
'
echo "All files processed!"