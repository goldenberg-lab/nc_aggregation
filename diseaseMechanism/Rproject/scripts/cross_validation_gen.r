#!/hpf/tools/centos6/R/3.1.1/bin/Rscript
if (length(commandArgs(TRUE))<3){print("Error: incorrect number of arguments");if(NA)print("Error");}
file=commandArgs(TRUE)[1];
outputdir=commandArgs(TRUE)[2];
nfold=as.numeric(commandArgs(TRUE)[3]);
header=FALSE; if(as.numeric(commandArgs(TRUE)[4])) header=TRUE;
paired=FALSE; if(length(commandArgs(TRUE))>4 && as.numeric(commandArgs(TRUE)[5])) paired=TRUE;

clinical=read.table(file,header=header,sep="\t");
indCases=which(clinical[,2]==1);
resampled=sample(1:length(indCases));
indCases=indCases[resampled];#making groups random
indContr=which(clinical[,2]==0);
if (!paired)resampled=sample(1:length(indContr));
indContr=indContr[resampled];#making groups random
groupassigment=rep(0,nrow(clinical));
groupassigment[indCases]=1+((1:length(indCases))%%nfold);
groupassigment[indContr]=1+((1:length(indContr))%%nfold);
for (i in 1:nfold){
write.table(clinical[which(groupassigment!=i),],paste(outputdir,"patients_training",i,".txt",sep=""),row.names=FALSE,col.names=FALSE,sep="\t")
write.table(clinical[which(groupassigment==i),],paste(outputdir,"patients_validation",i,".txt",sep=""),row.names=FALSE,col.names=FALSE,sep="\t")
}

