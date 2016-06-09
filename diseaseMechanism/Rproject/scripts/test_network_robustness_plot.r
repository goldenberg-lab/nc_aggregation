#!/hpf/tools/centos6/R/3.1.1/bin/Rscript
modnet=c(-0.5,-0.2,0.2,0.5)
nrep=20;
library(ggplot2)
netperf=array(NA,dim=c(length(modnet),nrep,13));
for (i in 1:length(modnet)){
temp=read.table(paste("performanceprim_",modnet[i],".txt",sep=""));
netperf[i,1:ncol(temp),]=t(temp)
}
dmres=read.table("performance10_0_100.txt");
toplotconflux=data.frame(netperf[1,,c(11,12)] , netperf[1,,c(6,7)],netperf[2,,c(6,7)], netperf[3,,c(6,7)],netperf[4,,c(6,7)])
toplotdms=data.frame(dmres , netperf[1,,c(9,10)],netperf[2,,c(9,10)], netperf[3,,c(9,10)],netperf[4,,c(9,10)])

methods=c("true_net","remove_50","remove_20","add_20","add_50")
nc=c("Conflux", "dmGWAS")
 methlab=as.factor(rep(sapply(methods,function(x)rep(x,nrep)),length(nc)));
 nlab=as.factor(t(matrix(nc,length(nc),nrep*length(methods))));
 methlab=factor(as.character(methlab),levels=c(methods[-1],methods[1]));#NEW 
 
  #sensitivity
  d1=c( as.numeric(as.matrix(toplotconflux[,2*(1:length(methods))-1])), as.numeric(as.matrix(toplotdms[,2*(1:length(methods))-1])));
  #precision
  d2=c( as.numeric(as.matrix(toplotconflux[,2*(1:length(methods))])), as.numeric(as.matrix(toplotdms[,2*(1:length(methods))])));


  pdf("net_sensitivity.pdf");
  update_geom_defaults("point", list(colour = NULL));
  theme_set(theme_gray(base_size = 30));
  ggplot(data=data.frame(Sensitivity=d1,Methods=methlab,Sample_size=nlab),aes(x=Sample_size,y=Sensitivity))+ geom_boxplot(aes(fill=Methods,color=Methods),width=0.6,alpha=0.65)+ theme_bw() +theme(axis.text=element_text(size=18),axis.title=element_text(size=22,face="bold"),legend.text=element_text(size=16),legend.title=element_text(size=22))#+coord_fixed(ratio = 8)#+stat_summary(fun.y=median,geom="line",aes(group=method,color=method))# Think about coloring or removing outlier ,outlier.colour=NA removes them
  dev.off();

  pdf("net_precision.pdf");
  update_geom_defaults("point", list(colour = NULL));
  theme_set(theme_gray(base_size = 30));
  ggplot(data=data.frame(Precision=d2,Methods=methlab,Sample_size=nlab),aes(x=Sample_size,y=Precision))+ geom_boxplot(aes(fill=Methods,color=Methods),width=0.6,alpha=0.65)+ theme_bw() +theme(axis.text=element_text(size=18),axis.title=element_text(size=22,face="bold"),legend.text=element_text(size=16),legend.title=element_text(size=22))#+coord_fixed(ratio = 8)#+stat_summary(fun.y=median,geom="line",aes(group=method,color=method))# Think about coloring or removing outlier ,outlier.colour=NA removes them
  dev.off();

select=which(methlab %in% c("true_net","remove_50","remove_20"));
select=which(methlab %in% c("true_net","add_50","add_20"))
  pdf("net_sensitivityadd.pdf");
  update_geom_defaults("point", list(colour = NULL));
  theme_set(theme_gray(base_size = 30));
  ggplot(data=data.frame(Sensitivity=d1[select],Methods=methlab[select],Sample_size=nlab[select]),aes(x=Sample_size,y=Sensitivity))+ geom_violin(aes(fill=Methods,color=Methods),width=0.6,alpha=0.65)+ theme_bw() +theme(axis.text=element_text(size=18),axis.title=element_text(size=22,face="bold"),legend.text=element_text(size=16),legend.title=element_text(size=22))+stat_summary(fun.y=median,geom="point", position=position_dodge(width=0.64),aes(group=Methods,color=Methods))#+coord_fixed(ratio = 8)#+stat_summary(fun.y=median,geom="line",aes(group=method,color=method))# Think about coloring or removing outlier ,outlier.colour=NA removes them
  dev.off();
  pdf("net_precisionadd.pdf");
  update_geom_defaults("point", list(colour = NULL));
  theme_set(theme_gray(base_size = 30));
  ggplot(data=data.frame(Precision=d2[select],Methods=methlab[select],Sample_size=nlab[select]),aes(x=Sample_size,y=Precision))+ geom_violin(aes(fill=Methods,color=Methods),width=0.6,alpha=0.65)+ theme_bw() +theme(axis.text=element_text(size=18),axis.title=element_text(size=22,face="bold"),legend.text=element_text(size=16),legend.title=element_text(size=22))+stat_summary(fun.y=median,geom="point", position=position_dodge(width=0.64),aes(group=Methods,color=Methods))#+coord_fixed(ratio = 8)#+stat_summary(fun.y=median,geom="line",aes(group=method,color=method))# Think about coloring or removing outlier ,outlier.colour=NA removes them
  dev.off();


