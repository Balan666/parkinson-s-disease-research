#!/bin/bash

echo "loading modules"
module load nixpkgs/16.09 StdEnv/2020 plink/2.00-10252019-avx2
module load r/4.0.2

plink --bfile ../AMP_PD_FILTERED_ALL_CHR --keep ./inter_files/PRS_top_25_fids.txt --make-bed --recode vcf --max-maf 0.01 --out ./inter_files/PRS_top_25
plink --bfile ../AMP_PD_FILTERED_ALL_CHR --keep ./inter_files/PRS_bottom_25_fids.txt --make-bed --recode vcf --max-maf 0.01 --out ./inter_files/PRS_bottom_25

plink --bfile ../AMP_PD_FILTERED_ALL_CHR --keep ./inter_files/PRS_bottom_25_case_fids.txt --make-bed --recode vcf --max-maf 0.01 --out ./inter_files/PRS_bottom_25_case

plink --bfile ../AMP_PD_FILTERED_ALL_CHR --keep ./inter_files/PRS_top_25_control_fids.txt --make-bed --recode vcf --max-maf 0.01 --out ./inter_files/PRS_top_25_control

mkdir -p SKAT_AMP_PD

ANNOVAR_DIR=../ANNOVAR/annovar
HUMANDB=${ANNOVAR_DIR}/humandb
SPLIT_DIR="./inter_files"
FILTER_DIR="./filtered_snp"
VCF_FILES=("PRS_bottom_25" "PRS_top_25" "PRS_bottom_25_case" "PRS_top_25_control")
GENERAL_FILES=("PRS_bottom_25" "PRS_top_25") 
RESULTING_GROUPS=("PRS_bottom_25" "PRS_top_25" "PRS_b_case_top_control") 

CHROMS=($(seq 1 22))

# Ensure output directories exist
mkdir -p "$FILTER_DIR"

# SPLITTING INTO CHROMOSOMES
for PREFIX in "${VCF_FILES[@]}"; do
    SPLIT_PATH="${SPLIT_DIR}/${PREFIX}_split"
    mkdir -p $SPLIT_PATH
    
    for CHR in "${CHROMS[@]}"; do
        CHR_VCF="${SPLIT_PATH}/${PREFIX}_chr${CHR}.vcf"
	
        # Extract header and chromosome-specific variants
        grep "^#" ./inter_files/${PREFIX}.vcf > "$CHR_VCF"
        grep -w "^$CHR" ./inter_files/${PREFIX}.vcf >> "$CHR_VCF"
    
    done
done        
        

# ANNOTATING GENERAL GROUPS SNPS
for PREFIX in "${GENERAL_FILES[@]}"; do
    SPLIT_PATH="${SPLIT_DIR}/${PREFIX}_split"
    mkdir -p $SPLIT_PATH

    for CHR in "${CHROMS[@]}"; do
        CHR_VCF="${SPLIT_PATH}/${PREFIX}_chr${CHR}.vcf"
        ANN_OUTPUT="${SPLIT_PATH}/${PREFIX}_chr${CHR}.out.annovar.hg38_multianno.txt"
	        
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
    done
done


for PREFIX in "${VCF_FILES[@]}"; do    
    for CHR in "${CHROMS[@]}"; do

        # Determine which general group this PREFIX belongs to
        if [[ "$PREFIX" == PRS_bottom_25* ]]; then
            BASE_GROUP="PRS_bottom_25"
        elif [[ "$PREFIX" == PRS_top_25* ]]; then
            BASE_GROUP="PRS_top_25"
        else
            echo "Unrecognized group in $PREFIX"
            continue
        fi

        # construct annotation path from the general group
        ANNOTATION="./inter_files/${BASE_GROUP}_split/${BASE_GROUP}_chr${CHR}.out.annovar.hg38_multianno.txt"

        # Define output paths
        FILTER_NON_SYN="${FILTER_DIR}/${PREFIX}_chr${CHR}.nonsyn.txt"
        FILTER_LOF="${FILTER_DIR}/${PREFIX}_chr${CHR}.LOF.txt"
        FILTER_ENCODE="${FILTER_DIR}/${PREFIX}_chr${CHR}.ENCODE.txt"
        FILTER_CADD="${FILTER_DIR}/${PREFIX}_chr${CHR}.CADD.txt"
        FILTER_ALL="${FILTER_DIR}/${PREFIX}_chr${CHR}.ALL.txt"

        # Run filtering commands
        awk -F "\t" '{ if (NR==1 || $9 ~ /nonsyn/ ) print $84"\t"$1"\t"$2"\t"$4"\t"$5"\t"$7}' "$ANNOTATION" > "$FILTER_NON_SYN"
        awk -F "\t" '{ if (NR==1 || $9 ~ /stop/ || $9 ~ /frame/ || ($9 ~ /intronic_splicing/ && $9 ~ /[+-][1-2]/) ) print $84"\t"$1"\t"$2"\t"$4"\t"$5"\t"$7}' "$ANNOTATION" > "$FILTER_LOF"
        awk -F "\t" '{ if (NR==1 || $6 == "exonic") print $84"\t"$1"\t"$2"\t"$4"\t"$5"\t"$7}' "$ANNOTATION" > "$FILTER_ENCODE"
        awk -F "\t" '{ if (NR==1 || $31 >= 20) print $84"\t"$1"\t"$2"\t"$4"\t"$5"\t"$7}' "$ANNOTATION" > "$FILTER_CADD"
        awk -F "\t" '{print $84"\t"$1"\t"$2"\t"$4"\t"$5"\t"$7}' "$ANNOTATION" > "$FILTER_ALL"
        
        echo "filtered $ANNOTATION annotation for ${FILTER_DIR}/${PREFIX}_chr${CHR}"
    done
done


for CHR in "${CHROMS[@]}"; do
    for CATEGORY in "nonsyn" "LOF" "ENCODE" "CADD" "ALL"; do
        files=(${FILTER_DIR}/PRS_*_25_c*_chr${CHR}.${CATEGORY}.txt)

        awk 'NR == 1 || FNR > 1' ${files[@]} > "${FILTER_DIR}/PRS_b_case_top_control_chr${CHR}.${CATEGORY}.txt"

    done
done

for PREFIX in "${RESULTING_GROUPS[@]}"; do
    for CHR in "${CHROMS[@]}"; do    
        # Run PLINK filtering for each category
        for CATEGORY in "nonsyn" "LOF" "ENCODE" "CADD" "ALL"; do
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
    done
done

echo "Processing complete."
