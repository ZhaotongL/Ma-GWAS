#!/bin/bash -l
#SBATCH -A sedaghat
#SBATCH --time=24:00:00
#SBATCH -p ram1t,amd2tb,amd512
#SBATCH --mem=100G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lin00374@umn.edu
#SBATCH -o /home/panwei/lin00374/snpnet/pbs/slurm-%j.out


#module load plink
module load plink/2.00-alpha-091019
plink2 --bfile 110K_QCed --pheno ukb_int.txt --covar noMET_20PC.txt  --covar-variance-standardize --glm hide-covar --out no_metabolites
plink2 --bfile 110K_QCed --pheno ukb_int.txt --covar MET_20PC.txt  --covar-variance-standardize --glm hide-covar  --out with_metabolites_PC
#plink2 --bfile 110K_QCed --pheno ukb_int.txt --covar MET_20PC.txt  --covar-col-nums 3-19 --covar-variance-standardize --glm hide-covar  --out with_metabolites_5PC
#plink2 --bfile 110K_QCed --pheno ukb_int.txt --covar MET_20PC.txt  --covar-col-nums 3-24 --covar-variance-standardize --glm hide-covar  --out with_metabolites_10PC
#plink2 --bfile 110K_QCed --pheno ukb_int.txt --covar MET_20PC_PRScs.txt  --covar-variance-standardize --glm hide-covar  --out with_metabolites_PC_PRScs




