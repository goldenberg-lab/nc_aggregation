#!/hpf/tools/centos6/R/3.1.1/bin/Rscript

patfile="/hpf/largeprojects/agoldenb/aziz/FHS/processed/LDL_High/patients.txt"
outputfile="/hpf/largeprojects/agoldenb/aziz/FHS/processed/LDL_High/patients1.txt"
indph=2;
indcov=c(5,6);

patfile=commandArgs(TRUE)[1];
outputfile=commandArgs(TRUE)[2];
indph=as.numeric(commandArgs(TRUE)[3]);
indcov=as.numeric(unlist(strsplit(commandArgs(TRUE)[4], ",")))

allpatients=read.table(patfile,sep="\t",stringsAsFactors =FALSE);
cov=as.factor(apply(allpatients,1,function(x)paste(x[indcov],collapse="_")));
phcov=data.frame(allpatients[,indph],cov)
covtab=table(phcov);

sample_from=function(covtab,phcov,i){
mini=min(covtab[,i]);
sampled=c();
for (j in 1:nrow(covtab))sampled=c(sampled,sample(which(phcov[,2]==colnames(covtab)[i] & phcov[,1]==rownames(covtab)[j]),mini));
return(sampled)
}

sampled=c();
for (i in 1: ncol(covtab))sampled=c(sampled,sample_from(covtab,phcov,i))

write.table(allpatients[sampled,],outputfile,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)

