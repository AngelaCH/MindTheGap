#gapper 4.0


args=(commandArgs(TRUE))

library(MBESS)
library(MASS)
library(mice)
library(snowfall)
library(BBmisc)


sfInit(parallel=T,cpus=16)


generator <- function(i){
  
  library(MBESS)
  library(MASS)
  library(mice)
  library(snowfall)
  library(BBmisc)
  library(random)
  

  set.seed(randomNumbers(1,-999999999,999999999,col=1))
  ##################################################################################################################################
  ##################################################################################################################################
  ##################################################################################################################################
  ##################################################################################################################################
  #Funktion, die aus Vorgaben wie Anzahl, mean, sd und cor einen multinormalverteilten Datensatz erstellt
  
  datmat <- function(N, mvec, sdvec, cor){
    covarmat <- cor2cov(cor, sdvec)
    datamat <- mvrnorm(n=N, mu=mvec, Sigma=covarmat, tol=1e-6, empirical=TRUE)#die erzeugten Datensätze haben die gewünschte Verteilung bis auf 1e-6 genau
    return(datamat)
  }
  
  ##################################################################################################################################
  #Funktion, die MCAR löscht
  
  mcar <- function(N, P, y){
    #generalmat<-readRDS("generalmat.rds")#da später innerhalb der generator() Funktion auf mcar zugegriffen wird, muss generalmat, das vorher in generator() erstellt wird, gesondert gespeichert werden, um in mcar() abrufbar zu sein
    mcarmat <- generalmat
    x <- round (N*(P/100))
    r <- sample(N,N)
    delete <- r[1:x]
    delete <- sort(delete)
    for(i in delete){
      mcarmat[i,y] <- NA
    }
    return(mcarmat)
  }
  
  ##################################################################################################################################
  #Funktion, die MAR löscht
  
  mar <- function(N, P, y, x, v){
    #add a new column to marmat which indicates the quantile x is in
    quant <- quantile(generalmat[,x])
    quant.x.generalmat <- ifelse(generalmat[,x]<=quant[2],1,
                                 ifelse(generalmat[,x]<=quant[3],2,
                                        ifelse(generalmat[,x]<=quant[4],3,
                                               4)))
    marmat <- c(generalmat[],quant.x.generalmat)
    dim(marmat) <- c(N,v + 1)
    while(sum(is.na(marmat[,y]))!=round((N*P/100))){
      marmat <- c(generalmat[],quant.x.generalmat)
      dim(marmat) <- c(N,v + 1)
      #how many missings are there in total? And how many cases are in each quantile?
      missing <- (N*P/100)
      z1 <- sum(marmat[,v+1]==1)
      z2 <- sum(marmat[,v+1]==2)
      z3 <- sum(marmat[,v+1]==3)
      z4 <- sum(marmat[,v+1]==4)
      #get the random Numbers for the seeds from random.org and determine % of missing in each quantile (based on the 25%/50% linear missings in Collins et al. (2001))
      #r <- randomNumbers((N+5),1,100,1)
      r <- sample(100,(N+5),replace=T)
      m1 <- (P/250)
      m2 <- (P*2/250)
      m3 <- (P*3/250)
      m4 <- (P*4/250)
      #determine number of missings in each quantile, by working with random.org draws (see above)
      r1 <- 0
      for(i in c(1:z1)){
        if(r[i]<=(m1*100)){
          r1 <- r1+1
        }
      }
      r2 <- 0
      for(i in c((z1+1):(z1+z2))){
        if(r[i]<=(m2*100)){
          r2 <- r2+1
        }
      }
      r3 <- 0
      for(i in c((z1+z2+1):(z1+z2+z3))){
        if(r[i]<=(m3*100)){
          r3 <- r3+1
        }
      }
      r4 <- 0
      for(i in c((z1+z2+z3+1):N)){
        if(r[i]<=(m4*100)){
          r4 <- r4+1
        }
      }
      #we want (r1+r2+r3+r4)==missing. Now ri can't be bigger than zi, but when sum(zi)<missing, additional ri have to be random.org-ly drawn.
      #because we have always a 10-20-30-40 ratio, the rest is distributed to the quantiles in this fashion
      #if in this additional rest-draw ri>zi, than ri==zi and the rest is taken care of in the end by looking where there is still space and filling it (see below)
      rest1 <- 0
      rest2 <- 0
      rest3 <- 0
      rest4 <- 0
      rest <- round(missing-(r1+r2+r3+r4))
      while(rest>0){
        #restr <- randomNumbers(rest,1,100,1)
        restr <- sample(100,rest)
        for(i in c(1:rest)){
          if((60<restr[i])&(r4<z4)){
            r4 <- r4+1
          }
          if(r4>z4){
            rest4 <- (r4-z4)
            r4 <- z4
          }
          if((30<restr[i])&(restr[i]<=60)&(r3<z3)){
            r3 <- r3+1
          }
          if(r3>z3){
            rest3 <- (r3-z3)
            r3 <- z3
          }
          if((10<restr[i])&(restr[i]<=30)&(r2<z2)){
            r2 <- r2+1
          }
          if(r2>z2){
            rest2 <- (r2-z2)
            r2 <- z2
          }
          if((restr[i]<=10)&(r1<z1)){
            r1 <- r1+1
          }
          if(r1>z1){
            rest1 <- (r1-z1)
            r1 <- z1
          }
          rest <- round(missing-(r1+r2+r3+r4))
        }
      }
      #as mentioned above, ri can never be bigger than zi, but sum(ri) can be bigger than missing.
      #If this occures, the ri will be reduced by the amount of rest in the ratio 40-30-20-10.
      #if so much is reduced that any ri<0, then this ri==0 and any other ri>0 will be reduced.
      rest1n <- 0
      rest2n <- 0
      rest3n <- 0
      rest4n <- 0
      if(rest<0){
        restr <- sample(100,-rest)
        for(i in c(1:-rest)){
          if((60<restr[i])&(r1>0)){
            r1 <- r1-1
          }
          if(r1<0){
            rest1n <- r1
            r1 <- 0
          }
          if((30<restr[i])&(restr[i]<=60)&(r2>0)){
            r2 <- r2-1
          }
          if(r2<0){
            rest2n <- r2
            r2 <- 0
          }
          if((10<restr[i])&(restr[i]<=30)&(r3>0)){
            r3 <- r3-1
          }
          if(r3<0){
            rest3n <- r3
            r3 <- 0
          }
          if((restr[i]<=10)&(r4>0)){
            r4 <- r4-1
          }
          if(r4<0){
            rest4n <- r4
            r4 <- 0
          }
        }
      }
      if((rest1n+rest2n+rest3n+rest4n)<0){
        for(i in c(1:-(rest1n+rest2n+rest3n+rest4n))){
          ifelse(r1>0,r1<-r1-1,
                 ifelse(r2>0,r2<-r2-1,
                        ifelse(r3>0,r3<-r3-1,r4<-r4-1)))
        }
      }
      #for each quantile, generate a list of row-numbers in which a datapoint of the quantile is to be found, and based on that generate a list for the rows of the to be deleted items in y
      g2.row.quant1 <- vector("numeric")
      j <- 1
      for(i in 1:N) {
        if(marmat[i,v + 1]==1) {
          g2.row.quant1[j] <- i
          j <- j+1
        }
      }
      set.seed(r[N+1])
      deletelistquant1 <- sample(g2.row.quant1,r1)
      for(i in deletelistquant1){
        marmat[i,y]<- NA
      } 
      g2.row.quant2 <- vector("numeric")
      j <- 1
      for(i in 1:N) {
        if(marmat[i,v + 1]==2) {
          g2.row.quant2[j] <- i
          j <- j+1
        }
      }
      set.seed(r[N+2])
      deletelistquant2 <- sample(g2.row.quant2,r2)
      for(i in deletelistquant2){
        marmat[i,y]<- NA
      }
      g2.row.quant3 <- vector("numeric")
      j <- 1
      for(i in 1:N) {
        if(marmat[i,v + 1]==3) {
          g2.row.quant3[j] <- i
          j <- j+1
        }
      } 
      set.seed(r[N+3])
      deletelistquant3 <- sample(g2.row.quant3,r3)
      for(i in deletelistquant3){
        marmat[i,y]<- NA
      }
      g2.row.quant4 <- vector("numeric")
      j <- 1
      for(i in 1:N) {
        if(marmat[i,v + 1]==4) {
          g2.row.quant4[j] <- i
          j <- j+1
        }
      }
      set.seed(r[N+4])
      deletelistquant4 <- sample(g2.row.quant4,r4)
      for(i in deletelistquant4){
        marmat[i,y]<- NA
      }
      
      if((rest1+rest2+rest3+rest4)>0){
        resteloeschen <- vector("numeric")
        j <- 1
        for(i in c(1:N)){
          if(marmat[i,y]!=NA){
            resteloeschen[j] <- i
            j <- j+1
          }
        }
        set.seed(r[N+5])
        resteloeschenlist <- sample(resteloeschen,(rest1+rest2+rest3+rest4))
        for(i in resteloeschenlist){
          marmat[i,y]<-NA
        }
      }
      #delete the extra quantile column in marmat
      marmat <- marmat[,-(v+1)]
    }
    return(marmat)
  }
  
  ##################################################################################################################################
  #Funktion, die MARCONVEX löscht
  
  marconvex <- function(N, P, y, x, v){
    #add a new column to marmat which indicates the quantile x is in
    quant <- quantile(generalmat[,x])
    quant.x.generalmat <- ifelse(generalmat[,x]<=quant[2],1,
                                 ifelse(generalmat[,x]<=quant[3],2,
                                        ifelse(generalmat[,x]<=quant[4],3,4)))
    marmat <- c(generalmat[],quant.x.generalmat)
    dim(marmat) <- c(N,v + 1)
    while(sum(is.na(marmat[,y]))!=round((N*P/100))){
      marmat <- c(generalmat[],quant.x.generalmat)
      dim(marmat) <- c(N,v + 1)
      #how many missings are there in total? And how many cases are in each quantile?
      missing <- (N*P/100)
      z1 <- sum(marmat[,v+1]==1)
      z2 <- sum(marmat[,v+1]==2)
      z3 <- sum(marmat[,v+1]==3)
      z4 <- sum(marmat[,v+1]==4)
      #get the random Numbers for the seeds from random.org and determine % of missing in each quantile (based on the 25%/50% linear missings in Collins et al. (2001))
      #r <- randomNumbers((N+5),1,100,1)
      r <- sample(100,(N+5),replace=T)
      m1 <- (P*4/250)
      m2 <- (P/250)
      m3 <- (P/250)
      m4 <- (P*4/250)
      #determine number of missings in each quantile, by working with random.org draws (see above)
      r1 <- 0
      for(i in c(1:z1)){
        if(r[i]<=(m1*100)){
          r1 <- r1+1
        }
      }
      r2 <- 0
      for(i in c((z1+1):(z1+z2))){
        if(r[i]<=(m2*100)){
          r2 <- r2+1
        }
      }
      r3 <- 0
      for(i in c((z1+z2+1):(z1+z2+z3))){
        if(r[i]<=(m3*100)){
          r3 <- r3+1
        }
      }
      r4 <- 0
      for(i in c((z1+z2+z3+1):N)){
        if(r[i]<=(m4*100)){
          r4 <- r4+1
        }
      }
      #we want (r1+r2+r3+r4)==missing. Now ri can't be bigger than zi, but when sum(zi)<missing, additional ri have to be random.org-ly drawn.
      #because we have always a 40-10-10-40 ratio, the rest is distributed to the quantiles in this fashion
      #if in this additional rest-draw ri>zi, than ri==zi and the rest is taken care of in the end by looking where there is still space and filling it (see below)
      rest1 <- 0
      rest2 <- 0
      rest3 <- 0
      rest4 <- 0
      rest <- round(missing-(r1+r2+r3+r4))
      while(rest>0){
        restr <- sample(100, rest)
        for(i in c(1:rest)){
          if((60<restr[i])&(r4<z4)){
            r4 <- r4+1
          }
          if(r4>z4){
            rest4 <- (r4-z4)
            r4 <- z4
          }
          if((50<restr[i])&(restr[i]<=60)&(r3<z3)){
            r3 <- r3+1
          }
          if(r3>z3){
            rest3 <- (r3-z3)
            r3 <- z3
          }
          if((40<restr[i])&(restr[i]<=50)&(r2<z2)){
            r2 <- r2+1
          }
          if(r2>z2){
            rest2 <- (r2-z2)
            r2 <- z2
          }
          if((restr[i]<=40)&(r1<z1)){
            r1 <- r1+1
          }
          if(r1>z1){
            rest1 <- (r1-z1)
            r1 <- z1
          }
          rest <- round(missing-(r1+r2+r3+r4))
        }
      }
      #as mentioned above, ri can never be bigger than zi, but sum(ri) can be bigger than missing.
      #If this occures, the ri will be reduced by the amount of rest in the ratio 10-40-40-10.
      #if so much is reduced that any ri<0, then this ri==0 and any other ri>0 will be reduced.
      rest1n <- 0
      rest2n <- 0
      rest3n <- 0
      rest4n <- 0
      if(rest<0){
        restr <- sample(100,-rest)
        for(i in c(1:-rest)){
          if((90<restr[i])&(r1>0)){
            r1 <- r1-1
          }
          if(r1<0){
            rest1n <- r1
            r1 <- 0
          }
          if((50<restr[i])&(restr[i]<=90)&(r2>0)){
            r2 <- r2-1
          }
          if(r2<0){
            rest2n <- r2
            r2 <- 0
          }
          if((10<restr[i])&(restr[i]<=50)&(r3>0)){
            r3 <- r3-1
          }
          if(r3<0){
            rest3n <- r3
            r3 <- 0
          }
          if((restr[i]<=10)&(r4>0)){
            r4 <- r4-1
          }
          if(r4<0){
            rest4n <- r4
            r4 <- 0
          }
        }
      }
      if((rest1n+rest2n+rest3n+rest4n)<0){
        for(i in c(1:-(rest1n+rest2n+rest3n+rest4n))){
          ifelse(r2>0,r2<-r2-1,
                 ifelse(r3>0,r3<-r3-1,
                        ifelse(r4>0,r4<-r4-1,r1<-r1-1)))
        }
      }
      #for each quantile, generate a list of row-numbers in which a datapoint of the quantile is to be found, and based on that generate a list for the rows of the to be deleted items in y
      g2.row.quant1 <- vector("numeric")
      j <- 1
      for(i in 1:N) {
        if(marmat[i,v + 1]==1) {
          g2.row.quant1[j] <- i
          j <- j+1
        }
      }
      set.seed(r[N+1])
      deletelistquant1 <- sample(g2.row.quant1,r1)
      for(i in deletelistquant1){
        marmat[i,y]<- NA
      } 
      g2.row.quant2 <- vector("numeric")
      j <- 1
      for(i in 1:N) {
        if(marmat[i,v + 1]==2) {
          g2.row.quant2[j] <- i
          j <- j+1
        }
      }
      set.seed(r[N+2])
      deletelistquant2 <- sample(g2.row.quant2,r2)
      for(i in deletelistquant2){
        marmat[i,y]<- NA
      }
      g2.row.quant3 <- vector("numeric")
      j <- 1
      for(i in 1:N) {
        if(marmat[i,v + 1]==3) {
          g2.row.quant3[j] <- i
          j <- j+1
        }
      } 
      set.seed(r[N+3])
      deletelistquant3 <- sample(g2.row.quant3,r3)
      for(i in deletelistquant3){
        marmat[i,y]<- NA
      }
      g2.row.quant4 <- vector("numeric")
      j <- 1
      for(i in 1:N) {
        if(marmat[i,v + 1]==4) {
          g2.row.quant4[j] <- i
          j <- j+1
        }
      }
      set.seed(r[N+4])
      deletelistquant4 <- sample(g2.row.quant4,r4)
      for(i in deletelistquant4){
        marmat[i,y]<- NA
      }
      
      if((rest1+rest2+rest3+rest4)>0){
        resteloeschen <- vector("numeric")
        j <- 1
        for(i in c(1:N)){
          if(marmat[i,y]!=NA){
            resteloeschen[j] <- i
            j <- j+1
          }
        }
        set.seed(r[N+5])
        resteloeschenlist <- sample(resteloeschen,(rest1+rest2+rest3+rest4))
        for(i in resteloeschenlist){
          marmat[i,y]<-NA
        }
      }
      #delete the extra quantile column in marmat
      marmat <- marmat[,-(v+1)]
    }
    return(marmat)
  }
  
  ##################################################################################################################################
  #Funktion, die MARSINISTER löscht
  
  marsinister <- function(N, P, y, x, v){   #notice: the last column/variable v is always considered as z!
    z_x <- scale(generalmat[,x])
    z_z <- scale(generalmat[,v])
    marmat <- c(generalmat[],-abs(z_x-z_z))
    dim(marmat) <- c(N,v + 1)
    while(sum(is.na(marmat[,y]))!=round((N*P/100))){
      marmat <- c(generalmat[],-abs(z_x-z_z))
      dim(marmat) <- c(N,v + 1)
      #add a new column to marmat which indicates the quantile x is in
      quant <- quantile(marmat[,v+1])
      quant.x.generalmat <- ifelse(marmat[,v+1]<=quant[2],1,
                                   ifelse(marmat[,v+1]<=quant[3],2,
                                          ifelse(marmat[,v+1]<=quant[4],3,4)))
      marmat <- c(marmat[],quant.x.generalmat)
      dim(marmat) <- c(N,v + 2)
      #how many missings are there in total? And how many cases are in each quantile?
      missing <- (N*P/100)
      z1 <- sum(marmat[,v+2]==1)
      z2 <- sum(marmat[,v+2]==2)
      z3 <- sum(marmat[,v+2]==3)
      z4 <- sum(marmat[,v+2]==4)
      #get the random Numbers for the seeds from random.org and determine % of missing in each quantile (based on the 25%/50% linear missings in Collins et al. (2001))
      r <- sample(100,(N+5),replace=T)
      m4 <- (P/250)
      m3 <- (P*2/250)
      m2 <- (P*3/250)
      m1 <- (P*4/250)
      #determine number of missings in each quantile, by working with random.org draws (see above)
      r1 <- 0
      for(i in c(1:z1)){
        if(r[i]<=(m1*100)){
          r1 <- r1+1
        }
      }
      r2 <- 0
      for(i in c((z1+1):(z1+z2))){
        if(r[i]<=(m2*100)){
          r2 <- r2+1
        }
      }
      r3 <- 0
      for(i in c((z1+z2+1):(z1+z2+z3))){
        if(r[i]<=(m3*100)){
          r3 <- r3+1
        }
      }
      r4 <- 0
      for(i in c((z1+z2+z3+1):N)){
        if(r[i]<=(m4*100)){
          r4 <- r4+1
        }
      }
      #we want (r1+r2+r3+r4)==missing. Now ri can't be bigger than zi, but when sum(zi)<missing, additional ri have to be random.org-ly drawn.
      #because we have always a 10-20-30-40 ratio, the rest is distributed to the quantiles in this fashion
      #if in this additional rest-draw ri>zi, than ri==zi and the rest is taken care of in the end by looking where there is still space and filling it (see below)
      rest1 <- 0
      rest2 <- 0
      rest3 <- 0
      rest4 <- 0
      rest <- round(missing-(r1+r2+r3+r4))
      while(rest>0){
        restr <- sample(100,rest)
        for(i in c(1:rest)){
          if((60<restr[i])&(r4<z4)){
            r4 <- r4+1
          }
          if(r4>z4){
            rest4 <- (r4-z4)
            r4 <- z4
          }
          if((30<restr[i])&(restr[i]<=60)&(r3<z3)){
            r3 <- r3+1
          }
          if(r3>z3){
            rest3 <- (r3-z3)
            r3 <- z3
          }
          if((10<restr[i])&(restr[i]<=30)&(r2<z2)){
            r2 <- r2+1
          }
          if(r2>z2){
            rest2 <- (r2-z2)
            r2 <- z2
          }
          if((restr[i]<=10)&(r1<z1)){
            r1 <- r1+1
          }
          if(r1>z1){
            rest1 <- (r1-z1)
            r1 <- z1
          }
          rest <- round(missing-(r1+r2+r3+r4))
        }
      }
      #as mentioned above, ri can never be bigger than zi, but sum(ri) can be bigger than missing.
      #If this occures, the ri will be reduced by the amount of rest in the ratio 40-30-20-10.
      #if so much is reduced that any ri<0, then this ri==0 and any other ri>0 will be reduced.
      rest1n <- 0
      rest2n <- 0
      rest3n <- 0
      rest4n <- 0
      if(rest<0){
        restr <- sample(100,-rest)
        for(i in c(1:-rest)){
          if((60<restr[i])&(r1>0)){
            r1 <- r1-1
          }
          if(r1<0){
            rest1n <- r1
            r1 <- 0
          }
          if((30<restr[i])&(restr[i]<=60)&(r2>0)){
            r2 <- r2-1
          }
          if(r2<0){
            rest2n <- r2
            r2 <- 0
          }
          if((10<restr[i])&(restr[i]<=30)&(r3>0)){
            r3 <- r3-1
          }
          if(r3<0){
            rest3n <- r3
            r3 <- 0
          }
          if((restr[i]<=10)&(r4>0)){
            r4 <- r4-1
          }
          if(r4<0){
            rest4n <- r4
            r4 <- 0
          }
        }
      }
      if((rest1n+rest2n+rest3n+rest4n)<0){
        for(i in c(1:-(rest1n+rest2n+rest3n+rest4n))){
          ifelse(r1>0,r1<-r1-1,
                 ifelse(r2>0,r2<-r2-1,
                        ifelse(r3>0,r3<-r3-1,r4<-r4-1)))
        }
      }
      #for each quantile, generate a list of row-numbers in which a datapoint of the quantile is to be found, and based on that generate a list for the rows of the to be deleted items in y
      g2.row.quant1 <- vector("numeric")
      j <- 1
      for(i in 1:N) {
        if(marmat[i,v + 2]==1) {
          g2.row.quant1[j] <- i
          j <- j+1
        }
      }
      set.seed(r[N+1])
      deletelistquant1 <- sample(g2.row.quant1,r1)
      for(i in deletelistquant1){
        marmat[i,y]<- NA
      } 
      g2.row.quant2 <- vector("numeric")
      j <- 1
      for(i in 1:N) {
        if(marmat[i,v + 2]==2) {
          g2.row.quant2[j] <- i
          j <- j+1
        }
      }
      set.seed(r[N+2])
      deletelistquant2 <- sample(g2.row.quant2,r2)
      for(i in deletelistquant2){
        marmat[i,y]<- NA
      }
      g2.row.quant3 <- vector("numeric")
      j <- 1
      for(i in 1:N) {
        if(marmat[i,v + 2]==3) {
          g2.row.quant3[j] <- i
          j <- j+1
        }
      } 
      set.seed(r[N+3])
      deletelistquant3 <- sample(g2.row.quant3,r3)
      for(i in deletelistquant3){
        marmat[i,y]<- NA
      }
      g2.row.quant4 <- vector("numeric")
      j <- 1
      for(i in 1:N) {
        if(marmat[i,v + 2]==4) {
          g2.row.quant4[j] <- i
          j <- j+1
        }
      }
      set.seed(r[N+4])
      deletelistquant4 <- sample(g2.row.quant4,r4)
      for(i in deletelistquant4){
        marmat[i,y]<- NA
      }
      
      if((rest1+rest2+rest3+rest4)>0){
        resteloeschen <- vector("numeric")
        j <- 1
        for(i in c(1:N)){
          if(marmat[i,y]!=NA){
            resteloeschen[j] <- i
            j <- j+1
          }
        }
        set.seed(r[N+5])
        resteloeschenlist <- sample(resteloeschen,(rest1+rest2+rest3+rest4))
        for(i in resteloeschenlist){
          marmat[i,y]<-NA
        }
      }
      #delete the extra quantile column in marmat
      marmat <- marmat[,-(v+1)]
      marmat <- marmat[,-(v+1)]
    }
    return(marmat)
  }
  
  ##################################################################################################################################
  #Funktion, die MNAR löscht
  
  mnar <- function(N, P, y, v){
    #add a new column to marmat which indicates the quantile x is in
    quant <- quantile(generalmat[,y])
    quant.y.generalmat <- ifelse(generalmat[,y]<=quant[2],1,
                                 ifelse(generalmat[,y]<=quant[3],2,
                                        ifelse(generalmat[,y]<=quant[4],3,4)))
    marmat <- c(generalmat[],quant.y.generalmat)
    dim(marmat) <- c(N,v + 1)
    while(sum(is.na(marmat[,y]))!=round((N*P/100))){
      marmat <- c(generalmat[],quant.y.generalmat)
      dim(marmat) <- c(N,v + 1)
      #how many missings are there in total? And how many cases are in each quantile?
      missing <- (N*P/100)
      z1 <- sum(marmat[,v+1]==1)
      z2 <- sum(marmat[,v+1]==2)
      z3 <- sum(marmat[,v+1]==3)
      z4 <- sum(marmat[,v+1]==4)
      #get the random Numbers for the seeds from random.org and determine % of missing in each quantile (based on the 25%/50% linear missings in Collins et al. (2001))
      r <- sample(100,(N+5),replace=T)
      m1 <- (P/250)
      m2 <- (P*2/250)
      m3 <- (P*3/250)
      m4 <- (P*4/250)
      #determine number of missings in each quantile, by working with random.org draws (see above)
      r1 <- 0
      for(i in c(1:z1)){
        if(r[i]<=(m1*100)){
          r1 <- r1+1
        }
      }
      r2 <- 0
      for(i in c((z1+1):(z1+z2))){
        if(r[i]<=(m2*100)){
          r2 <- r2+1
        }
      }
      r3 <- 0
      for(i in c((z1+z2+1):(z1+z2+z3))){
        if(r[i]<=(m3*100)){
          r3 <- r3+1
        }
      }
      r4 <- 0
      for(i in c((z1+z2+z3+1):N)){
        if(r[i]<=(m4*100)){
          r4 <- r4+1
        }
      }
      #we want (r1+r2+r3+r4)==missing. Now ri can't be bigger than zi, but when sum(zi)<missing, additional ri have to be random.org-ly drawn.
      #because we have always a 10-20-30-40 ratio, the rest is distributed to the quantiles in this fashion
      #if in this additional rest-draw ri>zi, than ri==zi and the rest is taken care of in the end by looking where there is still space and filling it (see below)
      rest1 <- 0
      rest2 <- 0
      rest3 <- 0
      rest4 <- 0
      rest <- round(missing-(r1+r2+r3+r4))
      while(rest>0){
        restr <- sample(100,rest)
        for(i in c(1:rest)){
          if((60<restr[i])&(r4<z4)){
            r4 <- r4+1
          }
          if(r4>z4){
            rest4 <- (r4-z4)
            r4 <- z4
          }
          if((30<restr[i])&(restr[i]<=60)&(r3<z3)){
            r3 <- r3+1
          }
          if(r3>z3){
            rest3 <- (r3-z3)
            r3 <- z3
          }
          if((10<restr[i])&(restr[i]<=30)&(r2<z2)){
            r2 <- r2+1
          }
          if(r2>z2){
            rest2 <- (r2-z2)
            r2 <- z2
          }
          if((restr[i]<=10)&(r1<z1)){
            r1 <- r1+1
          }
          if(r1>z1){
            rest1 <- (r1-z1)
            r1 <- z1
          }
          rest <- round(missing-(r1+r2+r3+r4))
        }
      }
      #as mentioned above, ri can never be bigger than zi, but sum(ri) can be bigger than missing.
      #If this occures, the ri will be reduced by the amount of rest in the ratio 40-30-20-10.
      #if so much is reduced that any ri<0, then this ri==0 and any other ri>0 will be reduced.
      rest1n <- 0
      rest2n <- 0
      rest3n <- 0
      rest4n <- 0
      if(rest<0){
        restr <- sample(100,-rest)
        for(i in c(1:-rest)){
          if((60<restr[i])&(r1>0)){
            r1 <- r1-1
          }
          if(r1<0){
            rest1n <- r1
            r1 <- 0
          }
          if((30<restr[i])&(restr[i]<=60)&(r2>0)){
            r2 <- r2-1
          }
          if(r2<0){
            rest2n <- r2
            r2 <- 0
          }
          if((10<restr[i])&(restr[i]<=30)&(r3>0)){
            r3 <- r3-1
          }
          if(r3<0){
            rest3n <- r3
            r3 <- 0
          }
          if((restr[i]<=10)&(r4>0)){
            r4 <- r4-1
          }
          if(r4<0){
            rest4n <- r4
            r4 <- 0
          }
        }
      }
      if((rest1n+rest2n+rest3n+rest4n)<0){
        for(i in c(1:-(rest1n+rest2n+rest3n+rest4n))){
          ifelse(r1>0,r1<-r1-1,
                 ifelse(r2>0,r2<-r2-1,
                        ifelse(r3>0,r3<-r3-1,r4<-r4-1)))
        }
      }
      #for each quantile, generate a list of row-numbers in which a datapoint of the quantile is to be found, and based on that generate a list for the rows of the to be deleted items in y
      g2.row.quant1 <- vector("numeric")
      j <- 1
      for(i in 1:N) {
        if(marmat[i,v + 1]==1) {
          g2.row.quant1[j] <- i
          j <- j+1
        }
      }
      set.seed(r[N+1])
      deletelistquant1 <- sample(g2.row.quant1,r1)
      for(i in deletelistquant1){
        marmat[i,y]<- NA
      } 
      g2.row.quant2 <- vector("numeric")
      j <- 1
      for(i in 1:N) {
        if(marmat[i,v + 1]==2) {
          g2.row.quant2[j] <- i
          j <- j+1
        }
      }
      set.seed(r[N+2])
      deletelistquant2 <- sample(g2.row.quant2,r2)
      for(i in deletelistquant2){
        marmat[i,y]<- NA
      }
      g2.row.quant3 <- vector("numeric")
      j <- 1
      for(i in 1:N) {
        if(marmat[i,v + 1]==3) {
          g2.row.quant3[j] <- i
          j <- j+1
        }
      } 
      set.seed(r[N+3])
      deletelistquant3 <- sample(g2.row.quant3,r3)
      for(i in deletelistquant3){
        marmat[i,y]<- NA
      }
      g2.row.quant4 <- vector("numeric")
      j <- 1
      for(i in 1:N) {
        if(marmat[i,v + 1]==4) {
          g2.row.quant4[j] <- i
          j <- j+1
        }
      }
      set.seed(r[N+4])
      deletelistquant4 <- sample(g2.row.quant4,r4)
      for(i in deletelistquant4){
        marmat[i,y]<- NA
      }
      
      if((rest1+rest2+rest3+rest4)>0){
        resteloeschen <- vector("numeric")
        j <- 1
        for(i in c(1:N)){
          if(marmat[i,y]!=NA){
            resteloeschen[j] <- i
            j <- j+1
          }
        }
        set.seed(r[N+5])
        resteloeschenlist <- sample(resteloeschen,(rest1+rest2+rest3+rest4))
        for(i in resteloeschenlist){
          marmat[i,y]<-NA
        }
      }
      #delete the extra quantile column in marmat
      marmat <- marmat[,-(v+1)]
    }
    return(marmat)
  }
  
  ##################################################################################################################################
  ##################################################################################################################################
  ##################################################################################################################################
  
  generalmatsmat <- matrix(,nrow=144000,ncol=5)
  j<-1
  for(corr in (1:6)){
    for(muster in (1:5)){
      for(N in (1:8)){
        for(Nmis in (1:6)){
          for(Faelle in (1:100)){
            generalmatsmat[j,]<-c(corr,muster,N,Nmis,Faelle)
            j<-j+1
          }
        }
      }
    }
  }
  

  
  #Funktion, die die missing-Liste und die vollständige Liste zeilenweise füllen kann
  
generator <- function(i){
  mvec <- c(5,5.2,0)
  sdvec <- c(1,1.44,1)
  corr<-generalmatsmat[i,1]
  muster<-generalmatsmat[i,2]
  Namount<-generalmatsmat[i,3]
  Nmis<-generalmatsmat[i,4]
  ifelse(corr==1,cor<-matrix(c(1,.6,.12,.6,1,.2,.12,.2,1),nrow=3,ncol=3),
         ifelse(corr==2,cor<-matrix(c(1,.6,.2,.6,1,.5,.2,.5,1),nrow=3,ncol=3),
                ifelse(corr==3,cor<-matrix(c(1,.6,.48,.6,1,.8,.48,.8,1),nrow=3,ncol=3),
                       ifelse(corr==4,cor<-matrix(c(1,.3,.12,.3,1,.2,.12,.2,1),nrow=3,ncol=3),
                              ifelse(corr==5,cor<-matrix(c(1,.3,.2,.3,1,.5,.2,.5,1),nrow=3,ncol=3),
                                     cor<-matrix(c(1,.3,.48,.3,1,.8,.48,.8,1),nrow=3,ncol=3))))))
  ifelse(Namount==1,N<-25,
         ifelse(Namount==2,N<-50,
                ifelse(Namount==3,N<-75,
                       ifelse(Namount==4,N<-100,
                              ifelse(Namount==5,N<-125,
                                     ifelse(Namount==6,N<-150,
                                            ifelse(Namount==7,N<-175,
                                                   N<-200)))))))
  ifelse(Nmis==1,P<-10,
         ifelse(Nmis==2,P<-20,
                ifelse(Nmis==3,P<-30,
                       ifelse(Nmis==4,P<-40,ifelse(Nmis==5,P<-50,
                                                   P<-60)))))
  generalmat<-datmat(N,mvec,sdvec,cor)
  #saveRDS(generalmat,"generalmat.rds")
  x<-1
  y<-2
  mcarvar<-3
  v<-3 #z
  if(muster==1){
    return(list(mcar(N,P,y),generalmat))
  }
  if(muster==2){
    return(list(mar(N,P,y,x,v),generalmat))
  }
  if(muster==3){
    return(list(marconvex(N,P,y,x,v),generalmat))
  }
  if(muster==4){
    return(list(marsinister(N,P,y,x,v),generalmat))
  }
  if(muster==5){
    return(list(mnar(N,P,y,v),generalmat))
  }
}


gmm<-sfLapply(c(70990:71000),generator)

sfStop()

save.image("gmm.RData")

