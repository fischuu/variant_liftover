#!/bin/bash

# Set default values for parameters
FLANK_SIZE=30
GENOME_OLD="genome.old.fasta"
GENOME_NEW="genome.new.fasta"
VARIANTS_FILE="variants.bed"
OUTPUT_BED="lifted_over.bed"
OUTPUT_VCF="lifted_over.vcf"
VALIDATION_FILE="validation_results.txt"

# Help message
usage() {
    echo "Usage: $0 [-f flank_size] [-o old_genome] [-n new_genome] [-v variants_file] [-b output_bed] [-c output_vcf]"
    echo "  -f : Number of bases up/downstream to extract (default: 30)"
    echo "  -o : Fasta file of the old genome"
    echo "  -n : Fasta file of the new genome"
    echo "  -v : Bed or vcf file with variants (chromosome, position)"
    echo "  -b : Output bed file with lifted-over positions"
    echo "  -c : Output VCF file with lifted-over positions (if VCF input is provided)"
    exit 1
}

# Parse command line arguments
while getopts "f:o:n:v:b:c:" opt; do
    case ${opt} in
        f ) FLANK_SIZE=$OPTARG ;;
        o ) GENOME_OLD=$OPTARG ;;
        n ) GENOME_NEW=$OPTARG ;;
        v ) VARIANTS_FILE=$OPTARG ;;
        b ) OUTPUT_BED=$OPTARG ;;
        c ) OUTPUT_VCF=$OPTARG ;;
        * ) usage ;;
    esac
done

# Check required arguments
if [[ -z "$GENOME_OLD" || -z "$GENOME_NEW" || -z "$VARIANTS_FILE" ]]; then
    usage
fi

# Index the genomes for samtools
if [[ ! -f "$GENOME_OLD.fai" ]]; then
    echo "Indexing old genome..."
    samtools faidx "$GENOME_OLD"
fi

if [[ ! -f "$GENOME_NEW.fai" ]]; then
    echo "Indexing new genome..."
    samtools faidx "$GENOME_NEW"
fi

# Index the genomes for Bowtie2
if [[ ! -f "${GENOME_OLD}.1.bt2" ]]; then
    echo "Indexing old genome with Bowtie2..."
    bowtie2-build "$GENOME_OLD" "${GENOME_OLD%.fasta}"
fi

if [[ ! -f "${GENOME_NEW}.1.bt2" ]]; then
    echo "Indexing new genome with Bowtie2..."
    bowtie2-build "$GENOME_NEW" "${GENOME_NEW%.fasta}"
fi

# Check if the input file is a VCF and convert to BED if necessary
if [[ "$VARIANTS_FILE" == *.vcf || "$VARIANTS_FILE" == *.vcf.gz ]]; then
    echo "Converting VCF to BED..."
    
    # If VCF is gzipped, use zcat, else cat
    if [[ "$VARIANTS_FILE" == *.vcf.gz ]]; then
        VCF_TOOL="zcat"
    else
        VCF_TOOL="cat"
    fi

    # Convert VCF to BED format
    VARIANTS_BED="variants_converted.bed"
    
    #$VCF_TOOL "$VARIANTS_FILE" | \
    #grep -v "^#" | \
    #awk -v FS="\t" -v OFS="\t" '{print $1, $2-1, $2}' > "$VARIANTS_BED"

    $VCF_TOOL "$VARIANTS_FILE" | \
    grep -v "^#" | \
    awk -v FS="\t" -v OFS="\t" '{print $1, $2-1, $2, $4, $5}' > "$VARIANTS_BED"

    VARIANTS_FILE_ORIG="$VARIANTS_FILE"
    VARIANTS_FILE="$VARIANTS_BED"
    INPUT_WAS_VCF=true
else
    INPUT_WAS_VCF=false
fi

# Temporary files for sequences and alignments
FLANKS_FILE="flanks.fasta"
ALIGN_OLD="alignment_old.sam"
ALIGN_NEW="alignment_new.sam"

# Step 1: Extract flanking sequences from the old genome using samtools
# Temporary files for sequences and alignments
OUTPUT_FILE="flanking_and_variant.txt"
: > "$OUTPUT_FILE"

# Step 1: Extract left, right flanking sequences and variant itself from the old genome using samtools
echo "Extracting left, right flanking sequences and the variant itself..."
while read -r CHR POS _; do
    START=$((POS - FLANK_SIZE))
    END=$((POS + FLANK_SIZE))
    
    # Extract the sequence using samtools faidx
    samtools faidx "$GENOME_OLD" "$CHR:$START-$END" >> "$FLANKS_FILE"
done < "$VARIANTS_FILE"


while read -r CHR POS END REF ALT; do
    # Calculate the positions for left flank, variant, and right flank
    LEFT_START=$((POS - FLANK_SIZE + 1))
    LEFT_END=$((POS))
    
    RIGHT_START=$((POS + 2))
    RIGHT_END=$((POS + FLANK_SIZE + 1))
    
    VARIANT_START=$((POS + 1))
    VARIANT_END=$((POS + 1))

    LEFT_FLANK=$(samtools faidx "$GENOME_OLD" "$CHR:$LEFT_START-$LEFT_END" | grep -v "^>" | tr -d '\n')
    VARIANT_SEQ=$(samtools faidx "$GENOME_OLD" "$CHR:$VARIANT_START-$VARIANT_END" | grep -v "^>" | tr -d '\n')
    RIGHT_FLANK=$(samtools faidx "$GENOME_OLD" "$CHR:$RIGHT_START-$RIGHT_END" | grep -v "^>" | tr -d '\n')
    
    VARIANT_NAME="${CHR}_${VARIANT_END}";

    # Write the extracted sequences and REF/ALT into the output file
    echo -e "$LEFT_FLANK\t$RIGHT_FLANK\t$VARIANT_SEQ\t$VARIANT_NAME\t$REF\t$ALT" >> "$OUTPUT_FILE"
done < "$VARIANTS_FILE"



echo "Flanking sequences and variant extraction complete. Results saved to $OUTPUT_FILE."


# Clean up extra newlines in the FASTA file
echo "Cleaning up the FASTA file..."
awk '/^>/ {if (NR!=1) print ""; print; next} {printf "%s", $0}' "$FLANKS_FILE" > "${FLANKS_FILE}_cleaned"
mv "${FLANKS_FILE}_cleaned" "$FLANKS_FILE"

# Step 2: Align the flanking sequences to the old genome using bowtie2
echo "Aligning flanking sequences to the old genome..."
bowtie2 -x "${GENOME_OLD%.fasta}" -f "$FLANKS_FILE" -S "$ALIGN_OLD"

# Step 3: Align the flanking sequences to the new genome using bowtie2
echo "Aligning flanking sequences to the new genome..."
bowtie2 -x "${GENOME_NEW%.fasta}" -f "$FLANKS_FILE" -S "$ALIGN_NEW"

# Step 4: Parse the alignments and generate the output BED file
echo "Generating output files..."
if [ "$INPUT_WAS_VCF" = true ]; then
# Generate VCF output
# Generate VCF output
echo -e "##fileformat=VCFv4.2" > "$OUTPUT_VCF"
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> "$OUTPUT_VCF"

# Initialize counters directly in awk
awk -v FLANK_SIZE=$FLANK_SIZE '
    BEGIN {
        FS = "\t"; 
        old_entries_count = 0;
        processed_entries_count = 0;
        processed_vcf_entries_count = 0;
    }

    # Process old VCF
    FILENAME == ARGV[1] {
        if ($1 ~ /^#/) next;

        pos_key = $1 ":" $2;
        old_ref[pos_key] = $4;
        old_alt[pos_key] = $5;
        old_id[pos_key] = $3;
        old_filter[pos_key] = $7;
        old_info[pos_key] = $8;

        old_entries_count++;

        next;
    }

    # Process old alignment SAM
    FILENAME == ARGV[2] {
        if ($1 ~ /^@/) next;

        if ($2 == 0 && !($1 in processed_old_reads)) {
            old_chrom[$1] = $3;  # Store old chromosome information
            old_pos[$1] = $4 + FLANK_SIZE + 1;  # Adjusted to FLANK_SIZE + 1
            old_ref[$1] = substr($10, FLANK_SIZE + 2, 1);  # Adjusted to capture correct base

            processed_old_reads[$1] = 1;
            processed_entries_count++;
        }
        next;
    }

    # Process new alignment SAM
    FILENAME == ARGV[3] {
        if ($1 ~ /^@/) next;

        # Check if this alignment is mapped (SAM flag field $2 != 4 means mapped)
        if ($2 != 4 && $1 in old_pos) {
            new_pos = $4 + FLANK_SIZE + 1;  # Adjusted to FLANK_SIZE + 1
            ref_base = substr($10, FLANK_SIZE + 2, 1);  # Adjusted to get correct base

            # Construct keys for both old and new alignments
            old_pos_key = old_chrom[$1] ":" old_pos[$1];  # Use the old chromosome
            new_pos_key = $3 ":" new_pos;

            if (old_pos_key in old_alt) {
                alt_value = old_alt[old_pos_key];
                id_value = old_id[old_pos_key];
                filter_value = old_filter[old_pos_key];
                info_value = old_info[old_pos_key];
                old_ref_value = old_ref[old_pos_key];
                old_chrom_value = old_chrom[$1];  # Store old chromosome for INFO field

                # Print the tab-separated VCF line
                print $3 "\t" new_pos "\t" id_value "\t" ref_base "\t" alt_value "\t.\t" filter_value "\tOLD=" old_pos[$1] ";OLD_REF=" old_ref_value ";OLD_CHROM=" old_chrom_value ";INFO=" info_value;
                processed_vcf_entries_count++;

                # Mark the read as processed so that only the first valid alignment is used
                processed_new_reads[new_pos_key] = 1;
            } else {
                # Print debug information for missing ALT entries
                print "Warning: No ALT entry found for " old_pos_key " - REF: " old_ref[old_pos_key] > "/dev/stderr";
            }
        }
        
    }

    END {
        print "Summary:" > "/dev/stderr";
        print "Total old VCF entries: " old_entries_count > "/dev/stderr";
        print "Total old alignment entries processed: " processed_entries_count > "/dev/stderr";
        print "Total new alignment entries processed: " length(processed_old_reads) > "/dev/stderr";  # Count entries in processed_old_reads
        print "Total VCF entries written: " processed_vcf_entries_count > "/dev/stderr";
    }
' "$VARIANTS_FILE_ORIG" "$ALIGN_OLD" "$ALIGN_NEW" >> "$OUTPUT_VCF"

# Debug output message
echo "Debugging complete. Check stderr for detailed debug information."

else
    # Generate BED output
    echo -e "#chrom\tstart\tend\tnew_chrom\tnew_start\tnew_end" > "$OUTPUT_BED"

    awk -v FLANK_SIZE=$FLANK_SIZE 'BEGIN{OFS="\t"} 
        NR==FNR {if($3=="=") {old_pos[$1]=$4; old_start[$1]=$2}} 
        NR!=FNR {if($3=="=" && $1 in old_pos) {print $1, old_start[$1], old_start[$1]+2*FLANK_SIZE, $1, $4, $4+2*FLANK_SIZE}}' \
        "$ALIGN_OLD" "$ALIGN_NEW" >> "$OUTPUT_BED"
fi

# Clean up
#rm "$ALIGN_OLD" "$ALIGN_NEW"

echo "Lift-over process complete."
if [ "$INPUT_WAS_VCF" = true ]; then
    echo "VCF results saved in $OUTPUT_VCF"
else
    echo "BED results saved in $OUTPUT_BED"
fi
echo "Validation results saved in $VALIDATION_FILE"
