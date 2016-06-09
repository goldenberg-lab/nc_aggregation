path="/home/aziz/Desktop/aziz/diseaseMechanism/";
folder="simul_propagate/"
library(ggplot2)
methods=c("Z-Integration","CAST","C-Alpha","VT","SKAT-O","Diff. Expr","SKAT-O+dmGWAS")
methods=c("Conflux","CAST","C-Alpha","VT","SKAT-O","SKAT-O+dmGWAS")
#methods=c("Z-Integration","CAST","C-Alpha","RWAS","SKAT","Diff. Expr","SKAT+dmGWAS")
mc=c(10,20,50,100);
nrep=20;nmeasures=3*length(methods)-1;
nc=c(50,100,200,400);
m=mc[1]
  resdata=array(0,dim=c(nrep,nmeasures,length(nc)));
  for(i in 1:length(nc)){
    n=nc[i];
    codepoint=paste(m,"_0_",n,sep="")
    prefix=paste(path,folder,"simul_",codepoint,"/",sep="");
    files=list.files(path=prefix, pattern=paste("(","performance",").*\\.txt$",sep=""))
    u=1;
    for (j in 1:length(files)){
      x=read.table(paste(prefix,files[j],sep=""));
      resdata[u:(u+dim(x)[1]-1),,i]=as.matrix(x);
      u=u+dim(x)[1];
    }
    if (u!=nrep+1){print(paste("problem in",codepoint));for(v in u:nrep) resdata[v,,i]=resdata[sample(1:(u-1),1),,i]}
  }

  #sensitivity
  d1=c();for (i in 1:length(nc))d1=c(d1,resdata[,3*(1:length(methods))-2,i]);
  #precision
  d2=c();for (i in 1:length(nc))d2=c(d2,resdata[,3*(1:length(methods))-1,i]);
  #N-precision
  d3=c();for (i in 1:length(nc))d3=c(d3,resdata[,3*(1:(length(methods)-1)),i], resdata[,3*length(methods)-2,i]);
  methlab=as.factor(rep(sapply(methods,function(x)rep(x,nrep)),length(nc)));
  nlab=as.factor(as.numeric(t(matrix(nc,length(nc),nrep*length(methods)))));
  methlab=factor(as.character(methlab),levels=c(methods[-1],methods[1]));#NEW 

  pdf(paste(path,folder,"res_sensitivity_",m,".pdf",sep=""));
  update_geom_defaults("point", list(colour = NULL));
  theme_set(theme_gray(base_size = 30));
  ggplot(data=data.frame(Sensitivity=d1,Methods=methlab,Sample_size=nlab),aes(x=Sample_size,y=Sensitivity))+ geom_boxplot(aes(fill=Methods,color=Methods),width=0.6,alpha=0.65)+ theme_bw() +theme(axis.text=element_text(size=18),axis.title=element_text(size=22,face="bold"),legend.text=element_text(size=16),legend.title=element_text(size=22))#+coord_fixed(ratio = 8)#+stat_summary(fun.y=median,geom="line",aes(group=method,color=method))# Think about coloring or removing outlier ,outlier.colour=NA removes them
  dev.off();
  pdf(paste(path,folder,"res_precision_",m,".pdf",sep=""));
  ggplot(data=data.frame(Precision=d2,Methods=methlab,Sample_size=nlab),aes(x=Sample_size,y=Precision))+ geom_boxplot(aes(fill=Methods,color=Methods),width=0.6,alpha=0.65)+ theme_bw()+theme(axis.text=element_text(size=18),axis.title=element_text(size=22,face="bold"),legend.text=element_text(size=16),legend.title=element_text(size=22))#+coord_fixed(ratio = 8)#+stat_summary(fun.y=median,geom="line",aes(group=method,color=method))#
  dev.off();
  pdf(paste(path,folder,"res_nprecision_",m,".pdf",sep=""));
ggplot(data=data.frame(N_Precision=d3,Methods=methlab,Sample_size=nlab),aes(x=Sample_size,y=N_Precision))+ geom_boxplot(aes(fill=Methods,color=Methods),width=0.6,alpha=0.65)+theme_bw()+theme(axis.text=element_text(size=18),axis.title=element_text(size=22,face="bold"),legend.text=element_text(size=16),legend.title=element_text(size=22))#+coord_fixed(ratio = 8)#+stat_summary(fun.y=median,geom="line",aes(group=method,color=method))#
  dev.off();



#plot results from medulloblastoma
path="/home/aziz/Desktop/aziz/diseaseMechanism/";

for (n in 1:3){for (m in 0:3){
x=unlist(read.table(paste(path,"medulloblastoma/results/final/lik",n,"_",m,".txt",sep="")))
y=unlist(read.table(paste(path,"medulloblastoma/results/final/pred",n,"_",m,".txt",sep="")))
indreal=which(1:length(x)%%11==1)
print(paste(n,":",m,":",length(which(x[-indreal] > x[indreal[1]])),"size=",length(x[-indreal])));
print(paste(n,":",m,":",length(which(y[-indreal] < y[indreal[1]])),"size=",length(x[-indreal])));
}
}
#TODO I can plot distribution of permutation statistics and one point for the real analysis to show how far we are from null hypothesis



#plot results for robustness
path="/home/aziz/Desktop/aziz/diseaseMechanism/";
folder="simul_robust/"
library(ggplot2)
methods=c("Z-Integration","No-Net","No-Expr","No-Harm","No-Net-Harm_expr")
methods=c("Conflux-All","No-Net","No-Harm","No-Net-Harm","SKAT-O")
#methods=c("Z-Integration","CAST","C-Alpha","RWAS","SKAT","Diff. Expr","SKAT+dmGWAS")
mc=c(10,20);
nrep=20;nmeasures=3*length(methods);
nc=c(50,100,200,400);
m=mc[2]
  resdata=array(0,dim=c(nrep,nmeasures,length(nc)));
  for(i in 1:length(nc)){
    n=nc[i];
    codepoint=paste(m,"_0_",n,sep="")
    prefix=paste(path,folder,"simul_",codepoint,"/",sep="");
    files=list.files(path=prefix, pattern=paste("(","performance",").*\\.txt$",sep=""))
    u=1;
    for (j in 1:length(files)){
      x=read.table(paste(prefix,files[j],sep=""));
      resdata[u:(u+dim(x)[1]-1),,i]=as.matrix(x);
      u=u+dim(x)[1];
    }
    if (u!=nrep+1){print(paste("problem in",codepoint));for(v in u:nrep) resdata[v,,i]=resdata[sample(1:(u-1),1),,i]}
  }

  #sensitivity
  d1=c();for (i in 1:length(nc))d1=c(d1,resdata[,3*(1:length(methods))-2,i]);
  #precision
  d2=c();for (i in 1:length(nc))d2=c(d2,resdata[,3*(1:length(methods))-1,i]);
  #N-precision
  d3=c();for (i in 1:length(nc))d3=c(d3,resdata[,3*(1:length(methods)),i]);
  methlab=as.factor(rep(sapply(methods,function(x)rep(x,nrep)),length(nc)));
  nlab=as.factor(as.numeric(t(matrix(nc,length(nc),nrep*length(methods)))));
  methlab=factor(as.character(methlab),levels=c(methods[-1],methods[1]));#NEW 

  pdf(paste(path,folder,"res_sensitivity_",m,".pdf",sep=""));
  update_geom_defaults("point", list(colour = NULL));
  theme_set(theme_gray(base_size = 30));
  ggplot(data=data.frame(Sensitivity=d1,Methods=methlab,Sample_size=nlab),aes(x=Sample_size,y=Sensitivity))+ geom_boxplot(aes(fill=Methods,color=Methods),width=0.6,alpha=0.65)+ theme_bw() +theme(axis.text=element_text(size=18),axis.title=element_text(size=22,face="bold"),legend.text=element_text(size=16),legend.title=element_text(size=22))#+coord_fixed(ratio = 8)#+stat_summary(fun.y=median,geom="line",aes(group=method,color=method))# Think about coloring or removing outlier ,outlier.colour=NA removes them
  dev.off();
  pdf(paste(path,folder,"res_precision_",m,".pdf",sep=""));
  ggplot(data=data.frame(Precision=d2,Methods=methlab,Sample_size=nlab),aes(x=Sample_size,y=Precision))+ geom_boxplot(aes(fill=Methods,color=Methods),width=0.6,alpha=0.65)+ theme_bw()+theme(axis.text=element_text(size=18),axis.title=element_text(size=22,face="bold"),legend.text=element_text(size=16),legend.title=element_text(size=22))#+coord_fixed(ratio = 8)#+stat_summary(fun.y=median,geom="line",aes(group=method,color=method))#
  dev.off();
  pdf(paste(path,folder,"res_nprecision_",m,".pdf",sep=""));
ggplot(data=data.frame(N_Precision=d3,Methods=methlab,Sample_size=nlab),aes(x=Sample_size,y=N_Precision))+ geom_boxplot(aes(fill=Methods,color=Methods),width=0.6,alpha=0.65)+theme_bw()+theme(axis.text=element_text(size=18),axis.title=element_text(size=22,face="bold"),legend.text=element_text(size=16),legend.title=element_text(size=22))#+coord_fixed(ratio = 8)#+stat_summary(fun.y=median,geom="line",aes(group=method,color=method))#
  dev.off();



