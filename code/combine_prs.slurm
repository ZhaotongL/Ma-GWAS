#!/bin/bash
#SBATCH -A sedaghat
#SBATCH --time=4:00:00
#SBATCH -p ram1t,amd2tb,amd512
#SBATCH --mem=100G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lin00374@umn.edu
#SBATCH -o /home/panwei/lin00374/snpnet/pbs/slurm-%j.out

module load plink

for i in {1..20}
do
    cat $(ls *PC${i}_*) > prs_PC${i}_wgt.txt
    rm $(ls *PC${i}_*chr*)
    plink --bfile ../../110K_QCed --score prs_PC${i}_wgt.txt 2 4 6 sum --out prs_PC$i
done

