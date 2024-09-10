#!/bin/bash

# Set default values for parameters
FLANK_SIZE=30
GENOME_OLD="genome.old.fasta"
GENOME_NEW="genome.new.fasta"
VARIANTS_FILE="variants.bed"
OUTPUT_BED="lifted_over.bed"

# Help message
usage() {
    echo "Usage: $0 [-f flank_size] [-o old_genome] [-n new_genome] [-v variants_file] [-b output_bed]"
    echo "  -f : Number of bases up/downstream to extract (default: 30)"
    echo "  -o : Fasta file of the old genome"
    echo "  -n : Fasta file of the new genome"
    echo "  -v : Bed file with variants (chromosome, position)"
    echo "  -b : Output bed file with lifted-over positions"
    exit 1
}

# Parse command line arguments
while getopts "f:o:n:v:b:" opt; do
    case ${opt} in
        f ) FLANK_SIZE=$OPTARG ;;
        o ) GENOME_OLD=$OPTARG ;;
        n ) GENOME_NEW=$OPTARG ;;
        v ) VARIANTS_FILE=$OPTARG ;;
        b ) OUTPUT_BED=$OPTARG ;;
        * ) usage ;;
    esac
done

# Check required arguments
if [[ -z "$GENOME_OLD" || -z "$GENOME_NEW" || -z "$VARIANTS_FILE" ]]; then
    usage
fi

# Temporary files for sequences and alignments
FLANKS_FILE="flanks.fasta"
ALIGN_OLD="alignment_old.sam"
ALIGN_NEW="alignment_new.sam"

# Step 1: Extract flanking sequences from the old genome using samtools
echo "Extracting flanking sequences..."
while read -r CHR POS _; do
    START=$((POS - FLANK_SIZE))
    END=$((POS + FLANK_SIZE))
    
    # Extract the sequence using samtools faidx
    samtools faidx "$GENOME_OLD" "$CHR:$START-$END" >> "$FLANKS_FILE"
done < "$VARIANTS_FILE"

# Step 2: Align the flanking sequences to the old genome
echo "Aligning flanking sequences to the old genome..."
bwa mem "$GENOME_OLD" "$FLANKS_FILE" > "$ALIGN_OLD"

# Step 3: Align the flanking sequences to the new genome
echo "Aligning flanking sequences to the new genome..."
bwa mem "$GENOME_NEW" "$FLANKS_FILE" > "$ALIGN_NEW"

# Step 4: Parse the alignments and generate the output BED file
echo "Generating output BED file..."
echo -e "#chrom\tstart\tend\tnew_chrom\tnew_start\tnew_end" > "$OUTPUT_BED"

awk 'BEGIN{OFS="\t"} 
    NR==FNR {if($3=="=") {old_pos[$1]=$4}} 
    NR!=FNR {if($3=="=" && $1 in old_pos) {print $1, old_pos[$1], old_pos[$1]+60, $1, $4, $4+60}}' \
    "$ALIGN_OLD" "$ALIGN_NEW" >> "$OUTPUT_BED"

# Clean up
rm "$FLANKS_FILE" "$ALIGN_OLD" "$ALIGN_NEW"

echo "Lift-over process complete. Results saved in $OUTPUT_BED"
