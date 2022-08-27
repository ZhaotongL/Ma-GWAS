rm(list=ls())
library(parallel)
library(tidyr)
ncore <- detectCores() -2
print(ncore)
load('HPT_rs6482423.RData')
n = 90000
beta_vec = c(-0.05,0,0.05)
gamma_vec = 1
cov1_eff = M1.null.mod.summ$coefficients[2:13,1]
cov2_eff = M1.null.mod.summ$coefficients[14:33,1]
paras = crossing(beta_vec,gamma_vec)

for(i in 1:nrow(paras)){
    beta = paras$beta_vec[i]
    gamma = paras$gamma_vec[i]
    out1 = mclapply(1:1000,function(x){
        set.seed(x)
        sample_id = sample(1:nrow(bed_mat),n,replace=F)
        sim_bed = scale(bed_mat[sample_id,])
        sim_cov2 = scale(met[sample_id,15:34])
        sim_cov1 = scale(met[sample_id,3:14])
        sige = 1 - var(beta*sim_bed + sim_cov1 %*% cov1_eff + M1.null.mod.summ$coefficients[1,1] + sim_cov2 %*% (gamma * cov2_eff))
        sim_lp = beta*sim_bed + sim_cov1 %*% cov1_eff + M1.null.mod.summ$coefficients[1,1] + sim_cov2 %*% (gamma * cov2_eff) 
        sim_p = exp(sim_lp)/(1+exp(sim_lp))
        sim_y = rbinom(prob=sim_p,n=n,size=1)
        m1 = summary(glm(sim_y~sim_cov1+sim_cov2+sim_bed,family='binomial'))
        m1.5 = summary(glm(sim_y~sim_cov1 + sim_cov2[,1:5]+sim_bed,family='binomial'))
        m1.10 = summary(glm(sim_y~sim_cov1 + sim_cov2[,1:10]+sim_bed,family='binomial'))
        m0 = summary(glm(sim_y~sim_cov1+sim_bed,family='binomial'))
        list(PC20= m1$coefficients[nrow(m1$coefficients),],PC0=m0$coefficients[nrow(m0$coefficients),],
             PC5= m1.5$coefficients[nrow(m1.5$coefficients),],PC10=m1.10$coefficients[nrow(m1.10$coefficients),])
        },mc.cores=ncore)
    error_ind = which(unlist(lapply(out1,function(x){class(x)=='try-error'})))
    if(length(error_ind)>0){
    print(error_ind)
    out1 = out1[-error_ind]}
    var_name = names(out1[[1]])
    for(vn in var_name){
        t = paste0(vn,'<-unlist(lapply(out1,function(x){x$',vn,'[1]}))')
        eval(parse(text=t))
        t1 = paste0("write.table(t(c(paras[i,],",vn,")),'",vn,"_b.txt',quote=F,row.names=F,append=T,col.names=F)")
        eval(parse(text=t1))
        t = paste0(vn,'<-unlist(lapply(out1,function(x){x$',vn,'[2]}))')
        eval(parse(text=t))
        t1 = paste0("write.table(t(c(paras[i,],",vn,")),'",vn,"_se.txt',quote=F,row.names=F,append=T,col.names=F)")
        eval(parse(text=t1))
        t = paste0(vn,'<-unlist(lapply(out1,function(x){x$',vn,'[4]}))')
        eval(parse(text=t))
        t1 = paste0("write.table(t(c(paras[i,],",vn,")),'",vn,"_p.txt',quote=F,row.names=F,append=T,col.names=F)")
        eval(parse(text=t1))
    }
}
