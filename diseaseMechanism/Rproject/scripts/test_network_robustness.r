#!/hpf/tools/centos6/R/3.1.1/bin/Rscript
#dirnet="/hpf/largeprojects/agoldenb/aziz/diseaseMechanism/simul/simul_10_0_100/rep1/"
dirnet=commandArgs(TRUE)[1];
modnet=as.numeric(commandArgs(TRUE)[2]);

load(paste(dirnet,"data.RData",sep=""))

netprim=net1;
tot=sum(unlist(sapply(net1,length)));
if (modnet>0){
toadd=cbind(sample(1:length(net1),tot* modnet,replace=TRUE), sample(1:length(net1),tot*modnet,replace=TRUE))
for (a in 1: nrow(toadd)){ netprim[[toadd[a,1]]]=c(netprim[[toadd[a,1]]], toadd[a,2]); netprim[[toadd[a,2]]]=c(netprim[[toadd[a,2]]], toadd[a,1]); }
for (a in 1: length(netprim))if(length(netprim[[a]]))netprim[[a]]=unique(netprim[[a]])
}
if (modnet<0){
netmat=data.frame(interactorA=rep("gene",200000),interactorB=rep("gene",200000),stringsAsFactors =FALSE)
k=1;for (i in 1:nbgenes)if(length(net1[[i]])){for(j in net1[[i]]){if (j >i){netmat[k,1]=genenames[i];netmat[k,2]=genenames[j];k=k+1;}}}
netmat=netmat[1:(k-1),];
toremove=sample(1:(k-1),tot*0.5*(-modnet));
netmat=netmat[-toremove,];
netfilename=paste(dirnet,"netrem",-modnet,".txt",sep="")
write.table(netmat,netfilename,sep="\t",row.names=FALSE,col.names=FALSE);
netprim=load_network_genes(netfilename,genenames,maxConnections)$network;
}
sum(unlist(sapply(netprim,length)))

powermeasuresprim=c(powermeasures[1,c(2,3,4,17,18)], 0,0,0,0,0,0,0,0); #measure sensitivity and precision for this method + others

xori <-  grid_search(codeDir,nbgenes,nbpatients,nbsnps,harm,harmgene,meancgenes,complexityDistr,pheno,hom,het,net1,NULL, cores,ratioSignal,decay,alpha,netparams=c(netparam,netmax,netmaxNoise,netrelaxed),removeExpressionOnly,truth,propagate=propagate);#WARNING changed e to NULL

xprim <- grid_search(codeDir,nbgenes,nbpatients,nbsnps,harm,harmgene,meancgenes,complexityDistr,pheno,hom,het,netprim,NULL, cores,ratioSignal,decay,alpha,netparams=c(netparam,netmax,netmaxNoise,netrelaxed),removeExpressionOnly,truth,propagate=propagate);#WARNING changed e to NULL

to=which(order(xori$h,decreasing=TRUE) %in% truth);found=which(xori$h>thresh);
sensori=length(which(found %in% truth))/length(truth);#sensitivity
precori=1;if (length(found)){precori=length(which(found %in% truth))/length(found);}#precision
nprecori=length(which(to<= length(truth)))/length(truth);#precision at rank

to=which(order(xprim$h,decreasing=TRUE) %in% truth);found=which(xprim$h>thresh);
powermeasuresprim[6]=length(which(found %in% truth))/length(truth);#sensitivity
if (length(found)){powermeasuresprim[7]=length(which(found %in% truth))/length(found);}else {powermeasuresprim[7]=1;}#precision
powermeasuresprim[8]=length(which(to<= length(truth)))/length(truth);#precision at rank

print(c(sensori,precori,nprecori));
print(powermeasuresprim)
powermeasuresprim[c(11,12,13)]=c(sensori,precori,nprecori);

if(length(methods)){# computing SKAT-O pvalues
library(AssotesteR)
library(SKAT)
genesToTest=unique(trans$gene);
labels=c( rep(1,length(indPatients)),rep(0,length(indControls)));
w=4 #was for loop 1 to 4
pvals=other_methods(trans,labels,genesToTest,methods[w],1000000,mf,dict);
}

if (testdmgwas){#dmGWAS (done on SKAT-O pvalues)
library(dmGWAS)
skatpvals=rep(0.5,nbgenes);skatpvals[genesToTest]=pvals;
skatpvals[which(skatpvals<=10^(-16))]=10^(-16);skatpvals[which(skatpvals>=1)]=0.5#+runif(length(which(skatpvals>=1)))/2;
d1=data.frame(gene=genenames,weight=skatpvals,stringsAsFactors =FALSE);
netmat=data.frame(interactorA=rep("gene",100000),interactorB=rep("gene",100000),stringsAsFactors =FALSE)
k=1;for (i in 1:nbgenes)if(length(netprim[[i]])){for(j in netprim[[i]]){if (j >i){netmat[k,1]=genenames[i];netmat[k,2]=genenames[j];k=k+1;}}}
netmat=netmat[1:(k-1),]
resdmgwas=dms(netmat,d1,expr1=NULL,expr2=NULL,d=1,r=0.1)
sel=resdmgwas$genesets.clear[[as.numeric(rownames(resdmgwas$zi.ordered)[1])]]
powermeasuresprim[9]=length(which(sel %in% genenames[truth]))/length(truth) #sensitivity
if (length(sel)){powermeasuresprim[10]=length(which(sel %in% genenames[truth]))/length(sel)
} else powermeasuresprim[10]= 1#precision
print("NAA method performance assessment: done");
}

write.table(powermeasuresprim,paste(dirnet,"performanceprim","_",modnet,".txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE);


