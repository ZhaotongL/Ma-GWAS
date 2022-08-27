**prepare_data.R**: Prepare UKB individual phenotypes and metabolite PCs (mPCs) via PCA. More explorations on mPCs.

**QC.slurm**: QC UKB genotype data and obtain top 10 genotype principal components.

**gwas.slurm**: Perform GWAS with different sets of covariates.

**compare_result.R**: Compare results between different GWAS (See Table 2, 5 in the main text)

**null_mod.R**: Perform analysis with only covariates (no SNP). (See Table 4 in the main text)

**mpc_gwas_prscs.slurm**: Perform GWAS on each of the 20 mPCs, and apply PRS-CS with the mPC GWAS results.

**combine_prs.slurm**: Combine results from PRS-CS and obtain the genetic components of the mPCs.