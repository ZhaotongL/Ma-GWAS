library(data.table)
library(tidyr)
library(dplyr)

ld_block = fread('/home/panwei/lin00374/metabolites/fourier_ls-all.bed')
ld_block$chr = substr(ld_block$chr,4,5)

a = fread(list.files(pattern='no.*glm'))
b = fread(list.files(pattern='with_metabolites_PC.PHENO1.*glm'))
#b = fread(list.files(pattern='with.*PRScs.*glm'))
#b = fread(list.files(pattern='5PC.*glm'))
#b = fread('with_metabolites_PC_PRScs.PHENO1.glm.linear')
c = fread(list.files(pattern='tsv'))
c %>% separate(variant,c("CHR","POS","A2","A1")) -> c
c$POS = as.numeric(c$POS)
c$CHR = as.numeric(c$CHR)


a$P = as.numeric(a$P)
b$P = as.numeric(b$P)
a.sig = a %>% filter(P<5e-8)
b.sig = b %>% filter(P<5e-8)
c.sig = c %>% filter(pval<5e-8)
a.sig$locus = 0
b.sig$locus = 0
for(i in 1:nrow(a.sig)){
    temp = ld_block %>% filter(chr==a.sig$`#CHROM`[i])
    id = findInterval(a.sig$POS[i],temp$start)
    a.sig$locus[i]=which(ld_block$chr==temp$chr[id] & ld_block$start == temp$start[id])
}

for(i in 1:nrow(b.sig)){
    temp = ld_block %>% filter(chr==b.sig$`#CHROM`[i])
    id = findInterval(b.sig$POS[i],temp$start)
    b.sig$locus[i]=which(ld_block$chr==temp$chr[id] & ld_block$start == temp$start[id])
}
c.sig = na.omit(c.sig)
c.sig$locus = 0
for(i in 1:nrow(c.sig)){
    temp = ld_block %>% filter(chr==c.sig$`CHR`[i])
    id = findInterval(c.sig$POS[i],temp$start)
    c.sig$locus[i]=which(ld_block$chr==temp$chr[id] & ld_block$start == temp$start[id])
}

a.sig_locus = unique(a.sig$locus)
b.sig_locus = unique(b.sig$locus)
c.sig_locus = unique(c.sig$locus)
new_locus = setdiff(b.sig_locus,a.sig_locus)
new_snps = setdiff(b.sig$ID,a.sig$ID)
replicated_locus = new_locus[new_locus %in% c.sig_locus]
print(length(new_locus))
print(length(replicated_locus))

c1 = fread('Meta-analysis_Wood_et_al+UKBiobank_2018.txt.gz')
c1.sig = c1 %>% filter(P<5e-8) ## GIANT, ICBP
c1.sig$locus = 0
for(i in 1:nrow(c1.sig)){
    temp = ld_block %>% filter(chr==c1.sig$`CHR`[i])
    id = findInterval(c1.sig$POS[i],temp$start)
    c1.sig$locus[i]=which(ld_block$chr==temp$chr[id] & ld_block$start == temp$start[id])
}
c1.sig_locus = unique(c1.sig$locus)
replicated_locus_1 = new_locus[new_locus %in% c1.sig_locus]
print(length(replicated_locus_1))



Res = list(a.sig = a.sig, b.sig = b.sig, c.sig = c.sig, c1.sig = c1.sig,
            new_locus = new_locus,
            replicated_locus = replicated_locus,
            replicated_locus_1 = replicated_locus_1)
save(Res,file='Res_ldblock.RData')


