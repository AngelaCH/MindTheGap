rm(list=ls(all=TRUE)) 

library(doParallel)
library(mice)
library(BBmisc)
library(pwr)

procs <- as.numeric(Sys.getenv("MOAB_PROCCOUNT"))
registerDoParallel(cores=procs)

args=(commandArgs(TRUE))
print(args)
if(length(args)==0){
  print("j not there")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

allgmms<-readRDS("allgmms.rds")

filler_pmm<-function(i,j){
  library(mice)
  missingdata<-allgmms[[i]][[1]]
  completedata<-allgmms[[i]][[2]]
  mm<-c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100)
  pred<- quickpred(missingdata)
  m<-mm[j]
  imp <- mice(missingdata, m=m, predictorMatrix=pred)
  ########
  #mean/mu
  ########
  Q_mu <- vector("numeric",length=m)
  U_mu <- vector("numeric",length=m)
  for (l in 1:m) {
    Q_mu[l] <- mean(complete(imp,l)$V2)
    U_mu[l] <- var(complete(imp,l)$V2)/(nrow(missingdata))
  }
  result_mu<-pool.scalar(Q_mu,U_mu)
  mu_Qbar_pmm<-result_mu$qbar
  mu_Ubar_pmm<-result_mu$ubar
  mu_B_pmm<-result_mu$b
  mu_T_pmm<-result_mu$t 
  mu_lo95_pmm<-result_mu$qbar-qt(.95,df=result_mu$df)*sqrt(result_mu$t)
  mu_hi95_pmm<-result_mu$qbar+qt(.95,df=result_mu$df)*sqrt(result_mu$t)
  mu_CI_pmm<-mu_hi95_pmm-mu_lo95_pmm
  mu_r_pmm<-result_mu$r
  mu_lambda_pmm<-(result_mu$b+result_mu$b/result_mu$m)/result_mu$t
  mu_fmi_pmm<-result_mu$f
  mu_sb_pmm<-(mu_Qbar_pmm-5.2)/(1.44/sqrt(nrow(missingdata)))  #difference between an average estimate and the pop parameter, divide by standard deviation (empirical standard error)
  ########
  #sd
  ########
  Q_sd <- vector("numeric",length=m)
  U_sd <- vector("numeric",length=m)
  for (o in 1:m) {
    Q_sd[o] <- sd(complete(imp,o)$V2)
    U_sd[o] <- var(complete(imp,o)$V2)/(2*(nrow(missingdata)))
  }
  result_sd<-pool.scalar(Q_sd,U_sd)
  sd_Qbar_pmm<-result_sd$qbar
  sd_Ubar_pmm<-result_sd$ubar
  sd_B_pmm<-result_sd$b
  sd_T_pmm<-result_sd$t
  sd_lo95_pmm<-result_sd$qbar-qt(.95,df=result_sd$df)*sqrt(result_sd$t)
  sd_hi95_pmm<-result_sd$qbar+qt(.95,df=result_sd$df)*sqrt(result_sd$t)
  sd_CI_pmm<-sd_hi95_pmm-sd_lo95_pmm
  sd_r_pmm<-result_sd$r
  sd_lambda_pmm<-(result_sd$b+result_sd$b/result_sd$m)/result_sd$t
  sd_fmi_pmm<-result_sd$f
  sd_sb_pmm<-(sd_Qbar_pmm-1.44)/(1.44/sqrt(2*nrow(missingdata)))
  ########
  #byx
  ########
  fit_byx <- with(data=imp, exp=lm(V2~V1))
  est_byx <- pool(fit_byx)
  byx_Qbar_pmm<-est_byx$qbar[[2]]
  byx_Ubar_pmm<-est_byx$ubar[2,2]
  byx_B_pmm<-est_byx$b[2,2]
  byx_T_pmm<-summary(est_byx)[2,3]
  byx_lo95_pmm<-summary(est_byx)[2,6]
  byx_hi95_pmm<-summary(est_byx)[2,7]
  byx_CI_pmm<-summary(est_byx)[2,7]-summary(est_byx)[2,6]
  byx_r_pmm<-est_byx$r[[2]]
  byx_lambda_pmm<-est_byx$lambda[[2]]
  byx_fmi_pmm<-est_byx$fmi[[2]]
  byx_sd<-summary(lm(V2~V1,as.data.frame(completedata)))[[4]][2,2]
  byx<-summary(lm(V2~V1,as.data.frame(completedata)))[[4]][2,1]
  byx_sb_pmm<-(byx_Qbar_pmm-byx)/byx_sd  
  ########
  #r2byx & pwr
  ######## 
  est_r2byx <- pool.r.squared(fit_byx)
  rsquarebyx <- matrix(NA, nrow = m, ncol = 3, dimnames = list(1:m, c("R^2", "Fisher trans F^2", "se()")))
  for (v in 1:m) {
    fit <- fit_byx$analyses[[v]]
    rsquarebyx[v, 1] <- sqrt(summary(fit)$r.squared)
    rsquarebyx[v, 2] <- 0.5 * log((rsquarebyx[v, 1] + 1)/(1 - rsquarebyx[v, 1]))
    rsquarebyx[v, 3] <- 1/(length(summary(fit)$residuals) - 3)
  }
  rsquarebyxfit <- pool.scalar(rsquarebyx[, 2], rsquarebyx[, 3])
  r2byx_Qbar<-est_r2byx[1]
  byx_pwr_pmm<-pwr.f2.test(u=1,v=(nrow(missingdata)-2),f2=r2byx_Qbar/(1-r2byx_Qbar),sig.level=.05)[[5]]
  ########
  #bxy
  ########
  fit_bxy <- with(data=imp, exp=lm(V1~V2))
  est_bxy <- pool(fit_bxy)
  bxy_Qbar_pmm<-est_bxy$qbar[[2]]
  bxy_Ubar_pmm<-est_bxy$ubar[2,2]
  bxy_B_pmm<-est_bxy$b[2,2]
  bxy_T_pmm<-summary(est_bxy)[2,3]
  bxy_lo95_pmm<-summary(est_bxy)[2,6]
  bxy_hi95_pmm<-summary(est_bxy)[2,7]
  bxy_CI_pmm<-summary(est_bxy)[2,7]-summary(est_bxy)[2,6]
  bxy_r_pmm<-est_bxy$r[[2]]
  bxy_lambda_pmm<-est_bxy$lambda[[2]]
  bxy_fmi_pmm<-est_bxy$fmi[[2]]
  bxy_sd<-summary(lm(V1~V2,as.data.frame(completedata)))[[4]][2,2]
  bxy<-summary(lm(V1~V2,as.data.frame(completedata)))[[4]][2,1]
  bxy_sb_pmm<-(bxy_Qbar_pmm-bxy)/bxy_sd 
  ########
  #r2bxy & pwr
  ######## 
  est_r2bxy <- pool.r.squared(fit_bxy)
  rsquarebxy <- matrix(NA, nrow = m, ncol = 3, dimnames = list(1:m, c("R^2", "Fisher trans F^2", "se()")))
  for (v in 1:m) {
    fit <- fit_bxy$analyses[[v]]
    rsquarebxy[v, 1] <- sqrt(summary(fit)$r.squared)
    rsquarebxy[v, 2] <- 0.5 * log((rsquarebyx[v, 1] + 1)/(1 - rsquarebyx[v, 1]))
    rsquarebxy[v, 3] <- 1/(length(summary(fit)$residuals) - 3)
  }
  rsquarebxyfit <- pool.scalar(rsquarebxy[, 2], rsquarebxy[, 3])
  r2bxy_Qbar<-est_r2bxy[1]
  bxy_pwr_pmm<-pwr.f2.test(u=1,v=(nrow(missingdata)-2),f2=r2bxy_Qbar/(1-r2bxy_Qbar),sig.level=.05)[[5]]
  ########
  #rxy
  ########
  Q_rxy <- vector("numeric",length=m)
  U_rxy <- vector("numeric",length=m)
  rho <- vector("numeric",length=m)
  for (p in 1:m) {
    rho[p] <- cor(complete(imp,p)$V1,complete(imp,p)$V2)
    Q_rxy[p] <- 1/2*log((1+rho[p])/(1-rho[p]))
    U_rxy[p] <- 1/(nrow(missingdata)-3)
  }
  result_rxy<-pool.scalar(Q_rxy,U_rxy)
  rxy_Qbar_pmm<-(exp(1)^(2*result_rxy$qbar)-1)/(exp(1)^(2*result_rxy$qbar)+1)
  rxy_Ubar_pmm<-result_rxy$ubar
  rxy_B_pmm<-result_rxy$b
  rxy_T_pmm<-result_rxy$t
  rxy_lo95_pmm<-(exp(1)^(2*(result_rxy$qbar-1.96*sqrt(1/(nrow(missingdata)-3))))-1)/(exp(1)^(2*(result_rxy$qbar-1.96*sqrt(1/(nrow(missingdata)-3))))+1)
  rxy_hi95_pmm<-(exp(1)^(2*(result_rxy$qbar+1.96*sqrt(1/(nrow(missingdata)-3))))-1)/(exp(1)^(2*(result_rxy$qbar+1.96*sqrt(1/(nrow(missingdata)-3))))+1)
  rxy_CI_pmm<-rxy_hi95_pmm-rxy_lo95_pmm
  rxy_r_pmm<-result_rxy$r
  rxy_lambda_pmm<-(result_rxy$b+result_rxy$b/result_rxy$m)/result_rxy$t
  rxy_fmi_pmm<-result_rxy$f 
  rxy_sd<-1/sqrt(nrow(missingdata)-3)
  rxycor<-cor(completedata[,1],completedata[,2])
  rxy_z<-.5*log((1+rxycor)/(1-rxycor))
  rxy_sb_pmm<-(result_rxy$qbar-rxy_z)/rxy_sd
  rxy_pwr_pmm<-pwr.r.test(n=nrow(missingdata),r=rxy_Qbar_pmm,sig.level=.05)[[4]]
  #in Liste speichern
  zwischenliste<-c(mu_Qbar_pmm,mu_Ubar_pmm,mu_B_pmm,mu_T_pmm,mu_lo95_pmm,mu_hi95_pmm,mu_CI_pmm,mu_r_pmm,mu_lambda_pmm,mu_fmi_pmm,mu_sb_pmm,
                   sd_Qbar_pmm,sd_Ubar_pmm,sd_B_pmm,sd_T_pmm,sd_lo95_pmm,sd_hi95_pmm,sd_CI_pmm,sd_r_pmm,sd_lambda_pmm,sd_fmi_pmm,sd_sb_pmm,
                   byx_Qbar_pmm,byx_Ubar_pmm,byx_B_pmm,byx_T_pmm,byx_lo95_pmm,byx_hi95_pmm,byx_CI_pmm,byx_r_pmm,byx_lambda_pmm,byx_fmi_pmm,byx_sb_pmm,byx_pwr_pmm,
                   bxy_Qbar_pmm,bxy_Ubar_pmm,bxy_B_pmm,bxy_T_pmm,bxy_lo95_pmm,bxy_hi95_pmm,bxy_CI_pmm,bxy_r_pmm,bxy_lambda_pmm,bxy_fmi_pmm,bxy_sb_pmm,bxy_pwr_pmm,
                   rxy_Qbar_pmm,rxy_Ubar_pmm,rxy_B_pmm,rxy_T_pmm,rxy_lo95_pmm,rxy_hi95_pmm,rxy_CI_pmm,rxy_r_pmm,rxy_lambda_pmm,rxy_fmi_pmm,rxy_sb_pmm,rxy_pwr_pmm)
  return(zwischenliste)
}

filler_norm<-function(i,j){
  library(mice)
  missingdata<-allgmms[[i]][[1]]
  completedata<-allgmms[[i]][[2]]
  mm<-c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100)
  pred<- quickpred(missingdata)
  m<-mm[j]
  imp <- mice(missingdata, meth="norm", m=m, predictorMatrix=pred)
  ########
  #mean/mu
  ########
  Q_mu <- vector("numeric",length=m)
  U_mu <- vector("numeric",length=m)
  for (l in 1:m) {
    Q_mu[l] <- mean(complete(imp,l)$V2)
    U_mu[l] <- var(complete(imp,l)$V2)/(nrow(missingdata))
  }
  result_mu<-pool.scalar(Q_mu,U_mu)
  mu_Qbar_norm<-result_mu$qbar
  mu_Ubar_norm<-result_mu$ubar
  mu_B_norm<-result_mu$b
  mu_T_norm<-result_mu$t 
  mu_lo95_norm<-result_mu$qbar-qt(.95,df=result_mu$df)*sqrt(result_mu$t)
  mu_hi95_norm<-result_mu$qbar+qt(.95,df=result_mu$df)*sqrt(result_mu$t)
  mu_CI_norm<-mu_hi95_norm-mu_lo95_norm
  mu_r_norm<-result_mu$r
  mu_lambda_norm<-(result_mu$b+result_mu$b/result_mu$m)/result_mu$t
  mu_fmi_norm<-result_mu$f
  mu_sb_norm<-(mu_Qbar_norm-5.2)/(1.44/sqrt(nrow(missingdata)))  #difference between an average estimate and the pop parameter, divide by standard deviation (empirical standard error)
  ########
  #sd
  ########
  Q_sd <- vector("numeric",length=m)
  U_sd <- vector("numeric",length=m)
  for (o in 1:m) {
    Q_sd[o] <- sd(complete(imp,o)$V2)
    U_sd[o] <- var(complete(imp,o)$V2)/(2*(nrow(missingdata)))
  }
  result_sd<-pool.scalar(Q_sd,U_sd)
  sd_Qbar_norm<-result_sd$qbar
  sd_Ubar_norm<-result_sd$ubar
  sd_B_norm<-result_sd$b
  sd_T_norm<-result_sd$t
  sd_lo95_norm<-result_sd$qbar-qt(.95,df=result_sd$df)*sqrt(result_sd$t)
  sd_hi95_norm<-result_sd$qbar+qt(.95,df=result_sd$df)*sqrt(result_sd$t)
  sd_CI_norm<-sd_hi95_norm-sd_lo95_norm
  sd_r_norm<-result_sd$r
  sd_lambda_norm<-(result_sd$b+result_sd$b/result_sd$m)/result_sd$t
  sd_fmi_norm<-result_sd$f
  sd_sb_norm<-(sd_Qbar_norm-1.44)/(1.44/sqrt(2*nrow(missingdata)))
  ########
  #byx
  ########
  fit_byx <- with(data=imp, exp=lm(V2~V1))
  est_byx <- pool(fit_byx)
  byx_Qbar_norm<-est_byx$qbar[[2]]
  byx_Ubar_norm<-est_byx$ubar[2,2]
  byx_B_norm<-est_byx$b[2,2]
  byx_T_norm<-summary(est_byx)[2,3]
  byx_lo95_norm<-summary(est_byx)[2,6]
  byx_hi95_norm<-summary(est_byx)[2,7]
  byx_CI_norm<-summary(est_byx)[2,7]-summary(est_byx)[2,6]
  byx_r_norm<-est_byx$r[[2]]
  byx_lambda_norm<-est_byx$lambda[[2]]
  byx_fmi_norm<-est_byx$fmi[[2]]
  byx_sd<-summary(lm(V2~V1,as.data.frame(completedata)))[[4]][2,2]
  byx<-summary(lm(V2~V1,as.data.frame(completedata)))[[4]][2,1]
  byx_sb_norm<-(byx_Qbar_norm-byx)/byx_sd  
  ########
  #r2byx & pwr
  ######## 
  est_r2byx <- pool.r.squared(fit_byx)
  rsquarebyx <- matrix(NA, nrow = m, ncol = 3, dimnames = list(1:m, c("R^2", "Fisher trans F^2", "se()")))
  for (v in 1:m) {
    fit <- fit_byx$analyses[[v]]
    rsquarebyx[v, 1] <- sqrt(summary(fit)$r.squared)
    rsquarebyx[v, 2] <- 0.5 * log((rsquarebyx[v, 1] + 1)/(1 - rsquarebyx[v, 1]))
    rsquarebyx[v, 3] <- 1/(length(summary(fit)$residuals) - 3)
  }
  rsquarebyxfit <- pool.scalar(rsquarebyx[, 2], rsquarebyx[, 3])
  r2byx_Qbar<-est_r2byx[1]
  byx_pwr_norm<-pwr.f2.test(u=1,v=(nrow(missingdata)-2),f2=r2byx_Qbar/(1-r2byx_Qbar),sig.level=.05)[[5]]
                      
 ########
 #bxy
 ########
 fit_bxy <- with(data=imp, exp=lm(V1~V2))
 est_bxy <- pool(fit_bxy)
 bxy_Qbar_norm<-est_bxy$qbar[[2]]
 bxy_Ubar_norm<-est_bxy$ubar[2,2]
 bxy_B_norm<-est_bxy$b[2,2]
 bxy_T_norm<-summary(est_bxy)[2,3]
 bxy_lo95_norm<-summary(est_bxy)[2,6]
 bxy_hi95_norm<-summary(est_bxy)[2,7]
 bxy_CI_norm<-summary(est_bxy)[2,7]-summary(est_bxy)[2,6]
 bxy_r_norm<-est_bxy$r[[2]]
 bxy_lambda_norm<-est_bxy$lambda[[2]]
 bxy_fmi_norm<-est_bxy$fmi[[2]]
 bxy_sd<-summary(lm(V1~V2,as.data.frame(completedata)))[[4]][2,2]
 bxy<-summary(lm(V1~V2,as.data.frame(completedata)))[[4]][2,1]
 bxy_sb_norm<-(bxy_Qbar_norm-bxy)/bxy_sd 
 ########
 #r2bxy & pwr
 ######## 
 est_r2bxy <- pool.r.squared(fit_bxy)
 rsquarebxy <- matrix(NA, nrow = m, ncol = 3, dimnames = list(1:m, c("R^2", "Fisher trans F^2", "se()")))
 for (v in 1:m) {
   fit <- fit_bxy$analyses[[v]]
   rsquarebxy[v, 1] <- sqrt(summary(fit)$r.squared)
   rsquarebxy[v, 2] <- 0.5 * log((rsquarebyx[v, 1] + 1)/(1 - rsquarebyx[v, 1]))
 }
 rsquarebxyfit <- pool.scalar(rsquarebxy[, 2], rsquarebxy[, 3])
 r2bxy_Qbar<-est_r2bxy[1]
 bxy_pwr_norm<-pwr.f2.test(u=1,v=(nrow(missingdata)-2),f2=r2bxy_Qbar/(1-r2bxy_Qbar),sig.level=.05)[[5]]
 ########
 #rxy
 ########
 Q_rxy <- vector("numeric",length=m)
 U_rxy <- vector("numeric",length=m)
 rho <- vector("numeric",length=m)
 for (p in 1:m) {
  rho[p] <- cor(complete(imp,p)$V1,complete(imp,p)$V2)
  Q_rxy[p] <- 1/2*log((1+rho[p])/(1-rho[p]))
  U_rxy[p] <- 1/(nrow(missingdata)-3)
 }
 result_rxy<-pool.scalar(Q_rxy,U_rxy)
 rxy_Qbar_norm<-(exp(1)^(2*result_rxy$qbar)-1)/(exp(1)^(2*result_rxy$qbar)+1)
 rxy_Ubar_norm<-result_rxy$ubar
 rxy_B_norm<-result_rxy$b
 rxy_T_norm<-result_rxy$t
 rxy_lo95_norm<-(exp(1)^(2*(result_rxy$qbar-1.96*sqrt(1/(nrow(missingdata)-3))))-1)/(exp(1)^(2*(result_rxy$qbar-1.96*sqrt(1/(nrow(missingdata)-3))))+1)
 rxy_hi95_norm<-(exp(1)^(2*(result_rxy$qbar+1.96*sqrt(1/(nrow(missingdata)-3))))-1)/(exp(1)^(2*(result_rxy$qbar+1.96*sqrt(1/(nrow(missingdata)-3))))+1)
 rxy_CI_norm<-rxy_hi95_norm-rxy_lo95_norm
 rxy_r_norm<-result_rxy$r
 rxy_lambda_norm<-(result_rxy$b+result_rxy$b/result_rxy$m)/result_rxy$t
 rxy_fmi_norm<-result_rxy$f 
 rxy_sd<-1/sqrt(nrow(missingdata)-3)
 rxycor<-cor(completedata[,1],completedata[,2])
 rxy_z<-.5*log((1+rxycor)/(1-rxycor))
 rxy_sb_norm<-(result_rxy$qbar-rxy_z)/rxy_sd
 rxy_pwr_norm<-pwr.r.test(n=nrow(missingdata),r=rxy_Qbar_norm,sig.level=.05)[[4]]
 #in Liste speichern
 zwischenliste<-c(mu_Qbar_norm,mu_Ubar_norm,mu_B_norm,mu_T_norm,mu_lo95_norm,mu_hi95_norm,mu_CI_norm,mu_r_norm,mu_lambda_norm,mu_fmi_norm,mu_sb_norm,
            sd_Qbar_norm,sd_Ubar_norm,sd_B_norm,sd_T_norm,sd_lo95_norm,sd_hi95_norm,sd_CI_norm,sd_r_norm,sd_lambda_norm,sd_fmi_norm,sd_sb_norm,
            byx_Qbar_norm,byx_Ubar_norm,byx_B_norm,byx_T_norm,byx_lo95_norm,byx_hi95_norm,byx_CI_norm,byx_r_norm,byx_lambda_norm,byx_fmi_norm,byx_sb_norm,byx_pwr_norm,
            bxy_Qbar_norm,bxy_Ubar_norm,bxy_B_norm,bxy_T_norm,bxy_lo95_norm,bxy_hi95_norm,bxy_CI_norm,bxy_r_norm,bxy_lambda_norm,bxy_fmi_norm,bxy_sb_norm,bxy_pwr_norm,
            rxy_Qbar_norm,rxy_Ubar_norm,rxy_B_norm,rxy_T_norm,rxy_lo95_norm,rxy_hi95_norm,rxy_CI_norm,rxy_r_norm,rxy_lambda_norm,rxy_fmi_norm,rxy_sb_norm,rxy_pwr_norm)
 return(zwischenliste)
}

multifun <- function(i,j){
  c(filler_pmm(i,j),filler_norm(i,j))
}

#Save resulting data set for j imputations
Egwene <- foreach(i=1:120000,.combine=rbind) %dopar% {
  multifun(i,j)
}

#Save Egwene_j (with index j)
Egwenename <- sprintf("Egwene_%02d",j)

myfun <- function(var, env = globalenv()) {
  assign(eval(substitute(var)), Egwene, envir = env)
  print(get(var))
}

myfun(Egwenename)

rm(Egwene)

#save wd gmm_j (with index j)
rm(allgmms,args,filler_pmm,filler_norm,Egwenename,i,myfun,procs)

save.image(sprintf("Egwene_%02d.RData", j))

rm(list=ls(all=TRUE)) 