library(dplyr)
library(data.table)
#y =  fread(list.files(pattern='ukb'));y = na.omit(y)
met = fread('/home/panwei/lin00374/metabolites/redo/dat.phe')
#id = intersect(y$f.eid,met$IID);fwrite(cbind(id,id),'keep.id',sep='\t',col.names=F)
id = read.table('keep.id')

met = merge(met,id,by.x=c('FID','IID'),by.y=c('V1','V2'))
met = met[,-c(5:14)]

met = as.data.frame(met)

### mean imputation ###
for(i in 1:ncol(met)) {
  met[ , i][is.na(met[ , i])] <- mean(met[ , i], na.rm = TRUE)
}

met.X = met[,-(1:4)]
### PCA ###
met.pca = prcomp(met.X,scale=TRUE)
met.20PC = met.pca$x[,1:20]
met.20PC = cbind(met[,1:4],met.20PC)

gen.PC = fread('110k_PC.eigenvec')
met.dat = merge(gen.PC,met.20PC,by.x=c('#FID','IID'),by.y=c('FID','IID'))
fwrite(met.dat,'MET_20PC.txt',sep='\t')
fwrite(met.dat[,1:14],'noMET_20PC.txt',sep='\t')


y =  fread(list.files(pattern='ukb.*txt$'))
y = merge(id,y,by.x='V1',by.y='f.eid')
y[,3] = qnorm((rank(y[,3],na.last="keep")-0.5)/sum(!is.na(y[,3]))) # inverse norm transformtaion
fwrite(y,'ukb_int.txt',sep='\t')

### Detail of PCA ###
library(factoextra)
eig.val <- get_eigenvalue(met.pca)
fwrite(data.frame(PC=paste0('PC',1:249),round(eig.val[,-1],2)),'PCvar_explained.txt',sep='\t')
var.pca <- get_pca_var(met.pca)
save(var.pca,file='res_pca.RData')
pca.con = as.data.frame(var.pca$contrib)
pca.con = pca.con[,1:20]
cn = colnames(pca.con)
pca.con$field = rownames(pca.con)
met_dic = fread('/home/panwei/shared/UKBiobankIndiv/metabolites/Metabolite_names_and_classes.txt') %>% select(-V1)
pca.con = merge(pca.con,met_dic, by.x='field',by.y='field.tab')
list_of_datasets = list()
for(i in 1:20){
    pca.con %>% select(-cn[-i]) %>% arrange(desc(.[2])) %>% head(n=10) -> tmp_df
    list_of_datasets[[paste0("PC",i)]] <- tmp_df
}
openxlsx::write.xlsx(list_of_datasets, file = "top10var_top20PC.xlsx")

