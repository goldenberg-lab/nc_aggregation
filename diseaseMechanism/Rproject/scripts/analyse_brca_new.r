#!/hpf/tools/centos6/R/3.1.1/bin/Rscript
#awk -F"\t" '{print $1}'  genotype.txt | sort | uniq -c > genotype_counts.txt
#load("/home/aziz/Desktop/aziz/diseaseMechanism_Results/brca/BRCA_methylation/BRCA_methyl_450__All__Both.rda")
#write.table(Data,"/home/aziz/Desktop/aziz/diseaseMechanism_Results/brca/processed/methyl_all.txt",col.names=substr(colnames(Data),9,12),row.names=Des[,1],sep="\t")
#x=load("/home/aziz/Desktop/aziz/diseaseMechanism_Results/brca/BRCA_methylation/BRCA_methyl_450__TSS1500__Both.rda")
#write.table(Data,"/home/aziz/Desktop/aziz/diseaseMechanism_Results/brca/processed/methyl_1500.txt",col.names=substr(colnames(Data),9,12),row.names=Des[,1],sep="\t")


path=commandArgs(TRUE)[1];
datapath=commandArgs(TRUE)[2];
resultdirparent=commandArgs(TRUE)[3];
phenocode1=as.numeric(commandArgs(TRUE)[4]);
phenocode2=as.numeric(commandArgs(TRUE)[5]);
analysetype=as.numeric(commandArgs(TRUE)[6]);
indpermut=as.numeric(commandArgs(TRUE)[7]);
npermut=as.numeric(commandArgs(TRUE)[8]);

genefile="genes.txt";patfile="patients.txt";varfile="variants.txt";exomefile="genotype.txt";exprfile="gene_expression_hareem.txt";methfile="methyl_body_hareem.txt";methfile2="methyl_200_hareem.txt";

#path="/hpf/largeprojects/agoldenb/aziz/diseaseMechanism/"
#datapath="/hpf/largeprojects/agoldenb/aziz/brca/processed/"
#resultdirparent="/hpf/largeprojects/agoldenb/aziz/brca/PatientsProba/"

local=TRUE;

if (local){
path="/home/aziz/Desktop/aziz/diseaseMechanism/";
indpermut=0;npermut=0;
phenocode1=1;#1 to 4 ; luninal A, B, basal, HER2+
phenocode2=3;#1 to 4 ; luninal A, B, basal, HER2+
analysetype=2;#1 exome only; 2 exome+expr; 3 exome+methyl
datapath="/home/aziz/Desktop/aziz/brca/processed/"
resultdirparent="/home/aziz/Desktop/aziz/brca/results/"
genefile="genes.txt";patfile="patients.txt";varfile="variants.txt";exomefile="genotype.txt";exprfile="gene_expression_hareem.txt";methfile="methyl_body_hareem.txt";methfile2="methyl_200_hareem.txt";
}

#Model parameters
meancgenes=20;
complexityDistr=c(0,1,0)/1;meancgenesppat=0;
ratioSignal=0.95;
alpha=0.5;
usenet2=TRUE;
premethod="None";#can be RUVinv, RUVrinv, RUVg, RUV4, RUV2, None
maxConnections=50;
netparams=c(0.9,0.01,0.01,2);
removeExpressionOnly=FALSE;
e=NULL;cores=1;
auto_balance=TRUE;

networkName=paste(path,"networks/BIOGRID3.2.98.tab2/interactions.txt",sep="") ## Biogrid network
#networkName=paste(path,"networks/HPRD_Release9_062910/interactions.txt",sep="") ## HPRD network
#networkName=paste(path,"networks/humannet/interactions150.txt",sep="");## HumanNet
#TODO change to whatever gene network you want to use. Can use other annotations.

codeDir=paste(path,"Rproject/",sep="");
dir.create(resultdirparent);
resultdir=paste(resultdirparent,"result_",phenocode1,"_",phenocode2,"_",analysetype,"_",indpermut,"/",sep="");dir.create(resultdir);

#load files
phenomat=read.table(paste(datapath,patfile,sep=""),sep="\t",stringsAsFactors =FALSE);
genotype=read.table(paste(datapath,exomefile,sep=""),sep="\t",stringsAsFactors =FALSE,colClasses=c("character","character","numeric"));
rawannot=read.table(paste(datapath,varfile,sep=""),sep="\t",stringsAsFactors =FALSE,colClasses=c("character","character","factor","numeric"));
annot=rawannot[order(rawannot[,2]),];print(table(annot[,3]));
genemat=read.table(paste(datapath,genefile,sep=""),sep="\t",stringsAsFactors =FALSE);

varids=1:nrow(annot);names(varids) <- annot[,1];
geneids=1:nrow(genemat);names(geneids)<- genemat[,1];
phenoname=phenomat[,2];
phenos=lapply(unique(phenoname),function(x)which(phenoname==x));
ph1=phenos[[phenocode1]];ph0=phenos[[phenocode2]];
pheno=c(rep(1,length(ph1)), rep(0,length(ph0)));
patientids=1:length(pheno); names(patientids)<- phenomat[c(ph1,ph0),1];

#Balance them by reducing the bigger(sampling or taking first n)
ph1=which(pheno==1);ph0=which(pheno==0);
if (auto_balance){
if(length(ph1)>length(ph0)){ph1=sample(ph1,length(ph0));
}else if (length(ph0)>length(ph1))ph0=sample(ph0,length(ph1));
pheno=c(rep(1,length(ph1)),rep(0,length(ph0)));
includedpatientids=patientids[c(ph1,ph0)];
ph1=1:length(ph1); ph0=(length(ph1)+1):(length(ph1)+length(ph0));
mappat=match(patientids,includedpatientids);
} else {mappat=1:length(patientids);}

nbgenes=length(geneids); nbpatients=length(pheno); genenames= names(geneids);

indsnp=rep(1,nrow(annot));nb=1;
for (i in 2:nrow(annot)){
if(annot[i,2]!=annot[i-1,2]){nb=1;}else nb=nb+1;
indsnp[i]=nb;
}

harm=list(); length(harm) <- nbgenes;
for(i in 1:nbgenes)harm[[i]]<- annot[which(annot[,2]==genenames[i]),4];
nbsnps=rep(0,nbgenes);for (i in 1:nbgenes)nbsnps[i]=length(harm[[i]]);
harmgene=rep(meancgenes/nbgenes,nbgenes);


library(preprocessCore);
library(Matrix)
library(ruv)
#library(RUVSeq)
source(paste(codeDir,"functions/load_network_functions.r",sep=""));
source(paste(codeDir,"functions/process_expr.r",sep=""));

if(analysetype>=2){#####load expression or methylation
if(analysetype==2){exprmatraw=read.table(paste(datapath,exprfile,sep=""),sep="\t");quantnormalize=TRUE;logfirst=TRUE;ksi=1;premethod="None";nfactor=50;}
if(analysetype==3){exprmatraw=read.table(paste(datapath,methfile,sep=""),sep="\t");quantnormalize=FALSE;logfirst=FALSE;premethod="None";nfactor=50;}
if(analysetype==4){exprmatraw=read.table(paste(datapath,methfile2,sep=""),sep="\t");quantnormalize=FALSE;logfirst=FALSE;premethod="None";nfactor=50;}
mapped=mappat[patientids[colnames(exprmatraw)]];#patient in expression/methylation data mapped to includedpatientids (or patientsids if mappat identity);
neg_controls=read.table(paste(datapath,"negative_controls.txt",sep=""),sep="\t",stringsAsFactors =FALSE)[,1];
e=preprocess_expr_all(exprmatraw,mapped,ph1,ph0,genenames,nbpatients,pheno_expr,quantnormalize,logfirst,neg_controls,premethod,ksi,nfactor);
}

if(length(e))print(paste("Number of mild aberrations: ",length(which(abs(e)>2)), ", strong :",length(which(abs(e)>3))))


#exprmatraw2=read.table(paste(datapath,methfile2,sep=""),sep="\t");quantnormalize=FALSE;logfirst=FALSE;premethod="None";nfactor=40;
#e2=preprocess_expr_all(exprmatraw2,mapped,ph1,ph0,genenames,nbpatients,pheno_expr,quantnormalize,logfirst,neg_controls,premethod,ksi,nfactor);
#print(paste("Number of mild aberrations: ",length(which(abs(e2)>2)), ", strong :",length(which(abs(e2)>3))))
#print(paste( length(which(e>2 & e2>2)), length(which(e>2 & e2< -2)), length(which(e< -2 & e2>2)), length(which(e< -2 & e2< -2)) ))
#print(paste( length(which(e>3 & e2>3)), length(which(e>3 & e2< -3)), length(which(e< -3 & e2>3)), length(which(e< -3 & e2< -3)) ))

#which(rowSums(abs(e)>2)>10)
#e1=NULL;x=NULL;
gc();

#load genotypes (list of variants in each inividual and zygocity)
het=list();length(het)<- nbgenes;for (i in 1:nbgenes){het[[i]]<- list(); length(het[[i]])<- nbpatients;} 
hom=list();length(hom)<- nbgenes;for (i in 1:nbgenes){hom[[i]]<- list(); length(hom[[i]])<- nbpatients;} 
#p=mappat[patientids[genotype[,2]]];ind=varids[genotype[,1]];g=geneids[annot[ind,2]];
p=mappat[match(genotype[,2],names(patientids))]; ind=match(genotype[,1],names(varids)); g=match(annot[ind,2],names(geneids));
pres=which(!is.na(p))
varuniq=unique(ind[pres]); 
geneuniq=geneids[annot[varuniq,2]];
transraw=list(gene=genenames[geneuniq],snps=varuniq ,mat=matrix(0,nbpatients,length(varuniq)));
for (i in 1:nrow(genotype)){
  if(!is.na(g[i]) && !is.na(p[i])){
    if (genotype[i,3]==1)het[[g[i]]][[p[i]]] <- c(het[[g[i]]][[p[i]]], indsnp[ind[i]])
    if (genotype[i,3]==2)hom[[g[i]]][[p[i]]] <- c(hom[[g[i]]][[p[i]]], indsnp[ind[i]])
    transraw$mat[p[i],which(varuniq==ind[i])]=genotype[i,3];
  }
}
prescolumns=which(colSums(transraw$mat)>0);
trans=list(gene=transraw$gene[prescolumns],snps=transraw$snps[prescolumns] ,mat=transraw$mat[,prescolumns])

#load network and optionally include second order neighbours as direct link
net1=load_network_genes(networkName,as.character(names(geneids)),maxConnections)$network;
net2=mapply(surround2,net1,1:length(net1) , MoreArgs=list(net1=net1));
net=net1;if (usenet2)net=net2;

#Load and apply the method
source(paste(codeDir,"functions/analyse_results.r",sep=""));
source(paste(codeDir,"newmethod/sumproductmem.r",sep=""))

acc=NULL;lik=NULL;acc0=NULL;lik0=NULL;
if(indpermut==0){
ptm <- proc.time();#Rprof(filename = "Rprof.out")
x<- grid_search(codeDir,nbgenes,nbpatients,nbsnps,harm,harmgene,meancgenes,complexityDistr,pheno,hom,het,net,e, cores,ratioSignal,alpha,netparams,removeExpressionOnly);
print(proc.time()-ptm);#summaryRprof(filename = "Rprof.out")

#Analyse results
bestgenes=order(x$h,decreasing=TRUE)[1:(2*meancgenes)];
print(genenames[bestgenes[1:meancgenes]]);print(x$h[bestgenes[1:meancgenes]])
write.table(t(x$h[bestgenes]),paste(resultdir,"genes.txt",sep=""),row.names=FALSE,col.names=as.character(genenames[bestgenes]))

if (x$status){
print(x$margC)
lik0=sum(x$likelihood);print(lik0);
acc0=t(x$predict-pheno)%*%(x$predict-pheno);print(acc0);
genestodraw=bestgenes;
plot_graphs(x,pheno,resultdir,genestodraw);
print(exp(x$munet[bestgenes,2]));
#pie_plot_patients(x$causes[[7]],bestgenes,genenames,resultdir,TRUE)
}
if (npermut>0)indpermut=indpermut+1
}

if (length(e)){
  for(i in 1:10){
  g=bestgenes[i];
  plot_expr_gene(paste(resultdir,genenames[g],".png",sep=""),g,genenames,exprmatraw,e,ph1,ph0,includedpatientids)
  }
}

#other methods
library(AssotesteR)
library(dmGWAS)
library(SKAT)
source(paste(codeDir,"functions/analyse_results.r",sep=""));
methods=c("CAST","CALPHA","VT","SKAT-O")#c("CAST","CALPHA","VT","SKAT-O");
if(length(methods)){
methodsfile=paste(resultdir,"methods.txt",sep="")
file.create(methodsfile)
genesToTest=unique(trans$gene);
for (w in 1: length(methods)){
pvals=other_methods(trans,pheno,genesToTest,methods[w],10000);
genesToTest1=genesToTest[which(apply(pvals,1,min)<0.05/length(genesToTest))];#8.10e-6 is the significance threshold after multiple hypothesis
write(methods[w],methodsfile, append=TRUE,sep="\t")
if(length(genesToTest1))write(genesToTest1,methodsfile, append=TRUE,sep="\t",ncolumns=length(genesToTest1));
write(genesToTest[order(apply(pvals,1,min))[1:20]],methodsfile, append=TRUE,sep="\t",ncolumns=20)
write(pvals[order(apply(pvals,1,min))[1:20]],methodsfile, append=TRUE,sep="\t",ncolumns=20)
}
print("Gene based tests performance assessment: done");
}

#Permutations
if (indpermut & npermut){
xp=list();length(xp)=npermut;
lik=rep(0,npermut);acc=rep(0,npermut);
if(npermut){ for (i in indpermut:(indpermut+npermut-1)){
phenop=sample(pheno);
xp[[i]]=grid_search(codeDir,nbgenes,nbpatients,nbsnps,harm,harmgene,meancgenes,complexityDistr,phenop,hom,het,net,e, cores,ratioSignal,alpha,netparams,removeExpressionOnly);
plot_graphs(xp[[i]],phenop,resultdir,NULL,i,reorder=TRUE);
lik[i]=sum(xp[[i]]$likelihood);
acc[i]=(t(xp[[i]]$predict-pheno)%*%(xp[[i]]$predict-pheno));
}}
}

write.table(c(acc0,acc),paste(resultdir,"pred.txt",sep=""),row.names=FALSE,col.names=FALSE)
write.table(c(lik0,lik),paste(resultdir,"lik.txt",sep=""),row.names=FALSE,col.names=FALSE)

