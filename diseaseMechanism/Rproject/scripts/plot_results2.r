#Same as plot result , but works if only performance files are downloaded in the same directory:path="/home/aziz/Desktop/aziz/diseaseMechanism/";
path="/home/aziz/Desktop/aziz/diseaseMechanism/";
folder="simul/simulrareonly"
output="res/"
library(ggplot2)
#methods=c("Z-Integration","CAST","C-Alpha","VT","SKAT-O","Diff. Expr","SKAT-O+dmGWAS")#with expression
methods=c("Conflux","CAST","C-Alpha","VT","SKAT-O","SKAT-O+dmGWAS")#exome only methods
#methods=c("Conflux","SKAT-O","SKAT-O+dmGWAS")
#methods=c("Conflux","no_network","no_harm","no_network_harm")#robustness to network/harmfulness
#methods=c("m=20","m=5","m=10","m=50","m=100");#robustness to meancgene parameter
#methods=c("r=0.05","r=0.25","r=0.5","r=0.75","r=1");#robustness to ratioSignal parameter
#methods=c("0.1","0.05","0.2","0.4");#robustness to decay parameter
variations=1;
mc=c(10,20,50,100);
nrep=20;nmeasures=3*length(methods);if ("SKAT-O+dmGWAS" %in% methods)nmeasures=nmeasures-1;
toremove=1+((1:variations)-1)*4;#WARNING removing first column which contains convergene status, for each variation of our method
nc=c(100,200,400);
m=mc[2];
dir.create(paste(path,folder,"/",output,sep=""))
  resdata=array(0,dim=c(nrep,nmeasures,length(nc)));
  for(i in 1:length(nc)){
    n=nc[i];
    codepoint=paste(m,"_0_",n,sep="")
    prefix=paste("performance",codepoint,sep="");
    files=list.files(path=paste(path,folder,sep=""), pattern=paste(prefix,".*\\.txt",sep=""),full.names=TRUE)
    u=1;
    for (j in 1:length(files)){
      x=read.table(files[j]);
      resdata[u:(u+dim(x)[1]-1),,i]=as.matrix(x[,-toremove]);#for wrong robustness3 experiment was[,-(1:4)][,-toremove]
      u=u+dim(x)[1];
    }
    if (u!=nrep+1){print(paste("problem in",codepoint));for(v in u:nrep) resdata[v,,i]=resdata[sample(1:(u-1),1),,i]}
  }

  #sensitivity
  d1=c();for (i in 1:length(nc))d1=c(d1,resdata[,3*(1:length(methods))-2,i]);
  #precision
  d2=c();for (i in 1:length(nc))d2=c(d2,resdata[,3*(1:length(methods))-1,i]);
  #N-precision #if dmgwas is there (last position, then nprecision is sensitivity for it)
  d3=c();if (length(methods)>1 & methods[length(methods)]=="SKAT-O+dmGWAS"){
  for (i in 1:length(nc))d3=c(d3,resdata[,3*(1:(length(methods)-1)),i], resdata[,3*length(methods)-2,i]);
  } else {for (i in 1:length(nc))d3=c(d3,resdata[,3*(1:length(methods)),i]);}
  methlab=as.factor(rep(sapply(methods,function(x)rep(x,nrep)),length(nc)));
  nlab=as.factor(as.numeric(t(matrix(nc,length(nc),nrep*length(methods)))));
  methlab=factor(as.character(methlab),levels=c(methods[-1],methods[1]));#NEW 

  pdf(paste(path,folder,"/",output,"res_sensitivity_",m,".pdf",sep=""));
  update_geom_defaults("point", list(colour = NULL));
  theme_set(theme_gray(base_size = 30));
  ggplot(data=data.frame(Sensitivity=d1,Methods=methlab,Sample_size=nlab),aes(x=Sample_size,y=Sensitivity))+ geom_boxplot(aes(fill=Methods,color=Methods),width=0.6,alpha=0.65)+ theme_bw() +theme(axis.text=element_text(size=18),axis.title=element_text(size=22,face="bold"),legend.text=element_text(size=16),legend.title=element_text(size=22))#+coord_fixed(ratio = 8)#+stat_summary(fun.y=median,geom="line",aes(group=method,color=method))# Think about coloring or removing outlier ,outlier.colour=NA removes them
  dev.off();
  pdf(paste(path,folder,"/",output,"res_precision_",m,".pdf",sep=""));
  ggplot(data=data.frame(Precision=d2,Methods=methlab,Sample_size=nlab),aes(x=Sample_size,y=Precision))+ geom_boxplot(aes(fill=Methods,color=Methods),width=0.6,alpha=0.65)+ theme_bw()+theme(axis.text=element_text(size=18),axis.title=element_text(size=22,face="bold"),legend.text=element_text(size=16),legend.title=element_text(size=22))#+coord_fixed(ratio = 8)#+stat_summary(fun.y=median,geom="line",aes(group=method,color=method))#
  dev.off();
  pdf(paste(path,folder,"/",output,"res_nprecision_",m,".pdf",sep=""));
ggplot(data=data.frame(N_Precision=d3,Methods=methlab,Sample_size=nlab),aes(x=Sample_size,y=N_Precision))+ geom_boxplot(aes(fill=Methods,color=Methods),width=0.6,alpha=0.65)+theme_bw()+theme(axis.text=element_text(size=18),axis.title=element_text(size=22,face="bold"),legend.text=element_text(size=16),legend.title=element_text(size=22))#+coord_fixed(ratio = 8)#+stat_summary(fun.y=median,geom="line",aes(group=method,color=method))#
  dev.off();



