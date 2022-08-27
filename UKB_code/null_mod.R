library(data.table)
library(dplyr)
y = fread('ukb_int.txt')
y = na.omit(y)

med = fread('MET_20PC.txt')
#y = y[which(y$V1 %in% med$IID),]
#med = fread('MET_20PC_PRScs.txt') ## environmental effect only
M1.lm.df = cbind(y[,3],med[,3:34])
M0.lm.df = cbind(y[,3],med[,3:14])
colnames(M0.lm.df)[1] = colnames(M1.lm.df)[1] = 'Y'


M0.null.mod = lm(Y~.,data=M0.lm.df)
M0.null.mod.summ = summary(M0.null.mod)
M1.null.mod = lm(Y~.,data=M1.lm.df)
M1.null.mod.summ = summary(M1.null.mod)
print(M0.null.mod.summ$adj.r.squared)
print(M1.null.mod.summ$adj.r.squared)
#cor_Y = cor(lm.df)[1,]

#Res = list(null.mod = null.mod,cor_Y= cor_Y)
save(Res,file='Null_res.RData')

# For binary trait Hypertension
library(rms)
M0.null.mod = lrm(Y~.,data=M0.lm.df)
print(M0.null.mod)
M1.null.mod = lrm(Y~.,data=M1.lm.df)
print(M1.null.mod)

# likelihood ratio test
library(lmtest)
lrt_M1_M0 = lrtest(M1.null.mod, M0.null.mod) 
save(lrt_M1_M0,file='lrt_M1_M0.RData')
