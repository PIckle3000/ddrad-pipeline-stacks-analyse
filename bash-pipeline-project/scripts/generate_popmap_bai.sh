#!/bin/bash
source "$PROJECT_ROOT/config/config.env"

INPUT_DIR="$BAM_DIR"         
POPMAP_OUT="$POPMAP"           
POP_NAME="$SPECIES"           

echo "Generating popmap at: $POPMAP_OUT"

#make the popmap with the PEaln_sorted.bam files
find "$INPUT_DIR" -maxdepth 1 -type f -name "*_aln_sorted.bam" | while read -r bamfile; do
    filename=$(basename "$bamfile" .bam)
    echo -e "${filename}\t${POP_NAME}"
done | sort | uniq > "$POPMAP_OUT"

touch "$POPMAP_OUT.done"