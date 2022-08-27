#!/bin/bash -l
#SBATCH -A sedaghat
#SBATCH --time=10:00:00
#SBATCH -p ram1t,amd2tb,amd512
#SBATCH --mem=100G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lin00374@umn.edu
#SBATCH -o /home/panwei/lin00374/snpnet/pbs/slurm-%j.out

cd $SLURM_SUBMIT_DIR
module load plink
for i in {1..22}
do
plink --bfile /home/panwei/shared/UKBiobankIndiv/genetic_data/ukb_cal_chr${i}_v2 --keep keep.id --geno 0.05 --maf 0.01 --hwe 1e-6 --make-bed --out chr${i} 
done

plink --bfile chr1 --merge-list  /home/panwei/lin00374/metabolites/genotype/mergelist.txt --make-bed --out 110K_QCed

module load plink/2.00-alpha-091019

plink2 --bfile 110K_QCed --pca 10 approx --out 110k_PC


