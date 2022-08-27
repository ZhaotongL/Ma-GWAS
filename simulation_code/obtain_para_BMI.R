library(dplyr)
library(data.table)
library(plink2R)
out_name = 'BMI_rs4639527'
  plink_command = paste0("module load plink; \n","plink --bfile /home/panwei/lin00374/metabolites/redo/BMI/110K_QCed",
                        " --extract rs4639527 ",
                        " --make-bed --out ",
                        out_name,sep="")
  plink_msg=system(plink_command,intern=TRUE)

  snp = tryCatch(read_plink(paste(out_name,sep="")),
                 error=function(e){cat("ERROR :",
                                       conditionMessage(e),
                                       "\n")})
met = fread('/home/panwei/lin00374/metabolites/redo/BMI/MET_20PC_PRScs.txt')
y = fread('/home/panwei/lin00374/metabolites/redo/BMI/ukb_int.txt')
M1.lm.df = cbind(y[,3],met[,3:34])
M1.lm.df = scale(M1.lm.df) %>% as.data.frame()
colnames(M1.lm.df)[1] = 'Y'
M1.null.mod = lm(Y~.,data=M1.lm.df)
M1.null.mod.summ = summary(M1.null.mod)

bed_mat = na.omit(snp$bed)
keep_id = unlist(strsplit(rownames(bed_mat),':'))[seq(1,by=2,length.out=nrow(bed_mat))]
met = met %>% filter(IID %in% keep_id) 
y = y %>% filter(V1 %in% keep_id)
M1.lm.df = cbind(y[,3],bed_mat,met[,3:34])
M1.lm.df = scale(M1.lm.df) %>% as.data.frame()
colnames(M1.lm.df)[1] = 'Y'
M1.mod = lm(Y~.,data=M1.lm.df)
M1.mod.summ = summary(M1.mod)

save.image(file='BMI_rs4639527.RData')

