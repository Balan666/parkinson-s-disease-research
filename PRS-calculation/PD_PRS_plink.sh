#!/bin/bash
## For PD PRS calculation using Plink
## No clumping, covar or p-value assumption yet

# how to use
## bfile_prefix=~/scratch/AMP_PD/all_chrs_file
## out=~/scratch/AMP_PD
## name=AMP_PD
## bash PD_PRS.sh $base $bfile_prefix $out $name

base=$1
bfile_prefix=$2
out=$3
name=$4

## 1. no pruning in target data
## 2. no clumping
## the snp in gwas file is renamed into chr:pos format
## gwas file is cut, only 3 columns left
## base=/lustre03/project/6004655/COMMUN/runs/eyu8/data/PRS/PD/UKB/meta5_hg38_chr_pos.tab
if [ -f $base ];then 
echo "the summary stat file for PD $base is present"
else
echo "please re-specify the path for your PD GWAS summary stat which should be meta5_hg38_chr_pos.tab"
fi 

## modify bim file into chr:pos format
module load scipy-stack/2020a python/3.8.10
python scripts/convert_bim_idto_chr_pos.py ${bfile_prefix} ${out}

## run plink2
module load nixpkgs/16.09 StdEnv/2020 plink/2.00-10252019-avx2
plink2 \
    --bim ${out}/${bfile_prefix}_modified.bim \
    --fam ${out}/${bfile_prefix}.fam \
    --bed ${out}/${bfile_prefix}.bed \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump $base \
    --clump-snp-field SNP \
    --clump-field p \
    --out $name

awk 'NR!=1{print $3}' ${name}.clumped >  ${name}.valid.snp
awk '{print $1,$7}' $base > SNP.pvalue

echo "0.001 0 0.001" > range_list 
echo "0.05 0 0.05" >> range_list
echo "0.1 0 0.1" >> range_list
echo "0.2 0 0.2" >> range_list
echo "0.3 0 0.3" >> range_list
echo "0.4 0 0.4" >> range_list
echo "0.5 0 0.5" >> range_list

plink \
    --bim ${out}/${bfile_prefix}_modified.bim \
    --fam ${out}/${bfile_prefix}.fam \
    --bed ${out}/${bfile_prefix}.bed \
    --score $base 1 2 5 header \
    --q-score-range range_list SNP.pvalue \
    --extract ${name}.valid.snp \
    --out ${name}

##prunning
plink \
    --bim ${out}/${bfile_prefix}_modified.bim \
    --fam ${out}/${bfile_prefix}.fam \
    --bed ${out}/${bfile_prefix}.bed \
    --indep-pairwise 200 50 0.25 \
    --out ${name}

##incorporating principal components (PCs) as covariates
plink \
    --bim ${out}/${bfile_prefix}_modified.bim \
    --fam ${out}/${bfile_prefix}.fam \
    --bed ${out}/${bfile_prefix}.bed \
    --extract ${name}.prune.in \
    --pca 6 \
    --out ${name}_6PC

plink \
    --bim ${out}/${bfile_prefix}_modified.bim \
    --fam ${out}/${bfile_prefix}.fam \
    --bed ${out}/${bfile_prefix}.bed \
    --extract ${name}.prune.in \
    --pca 10 \
    --out ${name}_10PC

plink \
    --bim ${out}/${bfile_prefix}_modified.bim \
    --fam ${out}/${bfile_prefix}.fam \
    --bed ${out}/${bfile_prefix}.bed \
    --extract ${name}.prune.in \
    --pca 15 \
    --out ${name}_15PC


# calculate zscore
module load scipy-stack/2020a python/3.8.10
python calculate_zscore_PRScs_PLINK2.py ${out} ${name}
