#!/bin/bash

echo "loading modules"
module load nixpkgs/16.09 StdEnv/2020 plink/2.00-10252019-avx2
module load r/4.0.2

plink --bfile ../AMP_PD_FILTERED_ALL_CHR --keep ./inter_files/PRS_top_25_fids.txt --make-bed --recode vcf --out ./inter_files/PRS_top_25
plink --bfile ../AMP_PD_FILTERED_ALL_CHR --keep ./inter_files/PRS_bottom_25_fids.txt --make-bed --recode vcf --out ./inter_files/PRS_bottom_25

mkdir -p SKAT_AMP_PD

ANNOVAR_DIR=..//annovar
HUMANDB=${ANNOVAR_DIR}/humandb
SPLIT_DIR="./inter_files"
FILTER_DIR="./filtered_snp"
VCF_FILES=("PRS_bottom_25" "PRS_top_25")

CHROMS=($(seq 1 22))

# Ensure output directories exist
mkdir -p "$FILTER_DIR"

for PREFIX in "${VCF_FILES[@]}"; do
    SPLIT_PATH="${SPLIT_DIR}/${PREFIX}_split"
    mkdir -p $SPLIT_PATH
    
    for CHR in "${CHROMS[@]}"; do
        CHR_VCF="${SPLIT_PATH}/${PREFIX}_chr${CHR}.vcf"
        ANN_OUTPUT="${SPLIT_PATH}/${PREFIX}_chr${CHR}.out.annovar.hg38_multianno.txt"
	
        # Extract header and chromosome-specific variants
        grep "^#" ./inter_files/${PREFIX}.vcf > "$CHR_VCF"
        grep -w "^$CHR" ./inter_files/${PREFIX}.vcf >> "$CHR_VCF"
        
        if [ -s "$CHR_VCF" ]; then
            echo "Processing: $CHR_VCF"

            if [ ! -f "$ANN_OUTPUT" ]; then
                echo "Running ANNOVAR for chromosome ${CHR}..."
                perl ${ANNOVAR_DIR}/table_annovar.pl \
                    "$CHR_VCF" "$HUMANDB" \
                    --buildver hg38 \
                    --out "${SPLIT_PATH}/${PREFIX}_chr${CHR}.out.annovar" \
                    --remove \
                    --protocol refGene,ljb26_all,dbnsfp41c \
                    --operation g,f,f \
                    --nastring . \
                    -vcfinput
            else
                echo "Annotation already exists for chromosome ${CHR}, skipping ANNOVAR."
            fi

            FILTER_NON_SYN="${FILTER_DIR}/${PREFIX}_chr${CHR}.nonsyn.txt"
            FILTER_LOF="${FILTER_DIR}/${PREFIX}_chr${CHR}.LOF.txt"
            FILTER_ENCODE="${FILTER_DIR}/${PREFIX}_chr${CHR}.ENCODE.txt"
            FILTER_CADD="${FILTER_DIR}/${PREFIX}_chr${CHR}.CADD.txt"

            awk -F "\t" '{ if (NR==1 || $9 ~ /nonsyn/ ) print $84"\t"$1"\t"$2"\t"$4"\t"$5"\t"$7}' "$ANN_OUTPUT" > "$FILTER_NON_SYN"
            awk -F "\t" '{ if (NR==1 || $9 ~ /stop/ || $9 ~ /frame/ || ($9 ~ /intronic_splicing/ && $9 ~ /[+-][1-2]/) ) print $84"\t"$1"\t"$2"\t"$4"\t"$5"\t"$7}' "$ANN_OUTPUT" > "$FILTER_LOF"
            awk -F "\t" '{ if (NR==1 || $6 == "exonic") print $84"\t"$1"\t"$2"\t"$4"\t"$5"\t"$7}' "$ANN_OUTPUT" > "$FILTER_ENCODE"
            awk -F "\t" '{ if (NR==1 || $31 >= 20) print $84"\t"$1"\t"$2"\t"$4"\t"$5"\t"$7}' "$ANN_OUTPUT" > "$FILTER_CADD"

            echo "filtered $ANN_OUTPUT annotation"
            # Run PLINK filtering for each category
            for CATEGORY in "nonsyn" "LOF" "ENCODE" "CADD"; do
                PLINK_INPUT="${SPLIT_DIR}/${PREFIX}"
                PLINK_OUTPUT="${FILTER_DIR}/${PREFIX}_chr${CHR}.${CATEGORY}"
                FILTER_FILE="${FILTER_DIR}/${PREFIX}_chr${CHR}.${CATEGORY}.txt"

                if [ -s "$FILTER_FILE" ]; then
                    plink --bfile "$PLINK_INPUT" --extract "$FILTER_FILE" --make-bed --out "$PLINK_OUTPUT"
                    echo "$PLINK_OUTPUT made"

                    echo "filtering covariances"
                    Rscript filter_cov.r "${PLINK_OUTPUT}.fam"
                    
                    echo "making SetID files for SKAT-O"
                    Rscript make_SetID.r "$FILTER_FILE"
                    
                    echo "Running SKAT-O"
                    Rscript SKAT_TMEM16F.R "$(basename $PLINK_OUTPUT)"
                else
                    echo "Skipping PLINK for ${CATEGORY}, no variants found."
                fi
            done

        else
            echo "Skipping missing or empty VCF: $CHR_VCF"
        fi
    done
done

echo "Processing complete."

