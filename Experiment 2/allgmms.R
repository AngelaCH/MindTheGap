allgmms<-list()

length(allgmms)<-120000

for(i in 1:120){
  name <- c("/home/fr/fr_fr/fr_ac1002/Experiment_2/",as.character(sprintf("gmm_%03d",i)),".RData")
  load(paste(name,collapse=""))
  end<-(i*1000)
  start<-(end-999)
  allgmms[start:end]<-get(sprintf("gmm_final_%03d",i))
}

saveRDS(allgmms,"/home/fr/fr_fr/fr_ac1002/Experiment_2/allgmms.rds")