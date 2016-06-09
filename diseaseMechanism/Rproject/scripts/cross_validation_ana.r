#!/hpf/tools/centos6/R/3.1.1/bin/Rscript
if (length(commandArgs(TRUE))<1){print("Error: incorrect number of arguments");if(NA)print("Error");}
dir=commandArgs(TRUE)[1];
#name=commandArgs(TRUE)[2];
fold=1:5;
ratio=(1:9)*10
complexity=1:3
nmeasures=18
output="result_all"
dir.create(paste(dir,output,sep=""))
  
resdata=array(0,dim=c(length(fold),length(complexity),length(ratio),nmeasures));
  for(i in 1:length(fold)){for (j in 1:length(complexity)){ for (k in 1:length(ratio)){    
    file=paste(dir,"/results_",fold[i],"-",complexity[j],"-",ratio[k],"/cross_validation_auc.txt",sep="");
    res=read.table(file);
    resdata[i,j,k,]=as.numeric(as.matrix(res)[1,1:nmeasures])[1:nmeasures];
  }}}

critera=apply(resdata[,,,2],c(2,3),mean)
critera2=apply(resdata[,,,4],c(2,3),mean)
write.table(format(critera,digits=4),paste(dir,output,"/mean_auc.txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE);
write.table(format(critera2,digits=4),paste(dir,output,"/mean_krus.txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE);




dir="/hpf/largeprojects/agoldenb/aziz/diseaseMechanism/simul/simul_10_0_200/rep1/"
fold=1:5;
ratio=((1:9)*10)
complexitydec=1:3
nmeasures=18
output="result_all"
dir.create(paste(dir,output,sep=""))
  
resdata=array(0,dim=c(length(fold),length(complexitydec),length(ratio),nmeasures));
  for(a in 1:length(fold)){for (b in 1:length(complexitydec)){ for (c in 1:length(ratio)){    
    file=paste(dir,"/results_",fold[a],"-",complexitydec[b],"-",ratio[c],"/data.RData",sep="");
    load(file); 
phenomat_v=read.table(paste(datapath,patfile_validation,sep=""),sep="\t",stringsAsFactors =FALSE);
 patientids_v=1:nrow(phenomat_v); names(patientids_v)<- phenomat_v[,1];
 het_v=list();length(het_v)<- nbgenes;for (i in 1:nbgenes){het_v[[i]]<- list(); length(het_v[[i]])<- nrow(phenomat_v);} 
 hom_v=list();length(hom_v)<- nbgenes;for (i in 1:nbgenes){hom_v[[i]]<- list(); length(hom_v[[i]])<- nrow(phenomat_v);}
 mappat_v=1:length(patientids_v); 
 p_v=mappat_v[match(genotype[,2],names(patientids_v))];
 for (i in 1:nrow(genotype)){
  if(!is.na(g[i]) && !is.na(p_v[i])){
    if (genotype[i,3]==1)het_v[[g[i]]][[p_v[i]]] <- c(het_v[[g[i]]][[p_v[i]]], indsnp[ind[i]])
    if (genotype[i,3]==2)hom_v[[g[i]]][[p_v[i]]] <- c(hom_v[[g[i]]][[p_v[i]]], indsnp[ind[i]])
  }
 }
 pheno_v=sumproduct_predict(x,het_v,hom_v,thresh=0.0);
 pheno_v50=sumproduct_predict(x,het_v,hom_v,thresh=0.5);
 pheno_vfinal=cbind(phenomat_v[,c(1,2)],pheno_v,pheno_v50);   
 pheno_tfinal=cbind(phenomat[,c(1,2)], x$predict);
 library(pROC)
   y1=pheno_vfinal[which(pheno_v>0.5),2];y2=pheno_vfinal[which(pheno_v>0.75),2];
   ord=order(pheno_v,decreasing=TRUE);y3=pheno_vfinal[ord[1:10],2];y4=pheno_vfinal[ord[1:20],2];
   stats=c(length(which(y1==1)),length(which(y1==0)),length(which(y2==1)),length(which(y2==0)), length(which(y3==1)),length(which(y3==0)),length(which(y4==1)),length(which(y4==0)));
   signmult=rep(1,length(pheno_v));signmult2=signmult;signmult[which(pheno_vfinal[ord,2]==0)]=-1;signmult2[which(pheno_vfinal[ord,2]==0)]=-2; #multiplicative
   map=mean(cumsum(as.numeric(pheno_vfinal[ord,2]==1))/(1:length(ord)))#mean average precision
   es1=which.max(cumsum(signmult*pheno_vfinal[ord,3]));es2=which.max(cumsum(signmult));es3=which.max(cumsum(signmult2));
   res=c(auc(pheno_tfinal[,2],pheno_tfinal[,3]),auc(pheno_vfinal[,2],pheno_vfinal[,3]),auc(pheno_vfinal[,2],pheno_vfinal[,4]),kruskal.test(pheno_vfinal[, 3],as.factor(pheno_vfinal[,2]))$p.value, length(which(x$h>0.5)),x$status, stats,es1,es2,es3,map,paste(genenames[which(x$h>0.5)],collapse='+'));
    resdata[a,b,c,]=res[1:nmeasures];
  }}}

matconv=resdata[,,,7];class(matconv)<- "numeric"; print(apply(matconv,c(2,3),mean))
matconv2=resdata[,,,8];class(matconv2)<- "numeric"; print(apply(matconv2,c(2,3),mean))
matconv1p=matconv/(matconv+matconv2);print(apply(matconv1p,c(2,3),mean))
matconv3=resdata[,,,15];class(matconv3)<- "numeric"; print(apply(matconv3,c(2,3),mean))
matconv4=resdata[,,,16];class(matconv4)<- "numeric"; print(apply(matconv4,c(2,3),mean))
matconv5=resdata[,,,17];class(matconv5)<- "numeric"; print(apply(matconv5,c(2,3),mean))

write.table(format(apply(matconv,c(2,3),mean),digits=4),paste(dir,output,"/mean_over50.txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE);
write.table(format(apply(matconv2,c(2,3),mean),digits=4),paste(dir,output,"/mean_over50false.txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE);
write.table(format(apply(matconv3,c(2,3),mean),digits=4),paste(dir,output,"/mean_es1.txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE);


