#!/hpf/tools/centos6/R/3.1.1/bin/Rscript
path=commandArgs(TRUE)[1];
phenocodesbin=as.numeric(commandArgs(TRUE)[2]);
analysistype=as.numeric(commandArgs(TRUE)[3]);
indpermut=as.numeric(commandArgs(TRUE)[4]);
npermut=as.numeric(commandArgs(TRUE)[5]);

local=TRUE;

if (local){
path="/home/aziz/Desktop/aziz/diseaseMechanism/";
phenocodesbin=3;#1 to 5 ; which subtypes are cases which are controls
analysistype=2;#0 is germline+somatic, 1 is germline, 2 is just somatic
indpermut=1;npermut=0;
}

#model parameters
meancgenes=10;
complexityDistr=c(0.33,0.33,0.33);meancgenesppat=0;
cores=1;
pcontrol=0.9;
alpha=0.4;
usenet2=TRUE;
maxConnections=50;
netparams=c(0.9,0.01,0.01,3);
removeExpressionOnly=FALSE;
e=NULL;
networkName=paste(path,"networks/BIOGRID3.2.98.tab2/interactions.txt",sep="") ## Biogrid network
#networkName=paste(path,"networks/HPRD_Release9_062910/interactions.txt",sep="") ## HPRD network
networkName=paste(path,"networks/humannet/interactions150.txt",sep="");## HumanNet

if(phenocodesbin==1)phenocodes=c(1,0,0);if(phenocodesbin==2)phenocodes=c(0,1,0);if(phenocodesbin==3)phenocodes=c(0,0,1);if(phenocodesbin==4)phenocodes=c(0,1,1);if(phenocodesbin==5)phenocodes=c(1,1,1);
codeDir=paste(path,"Rproject/",sep="");
datadir=paste(path,"medulloblastoma/subtypes/",sep="");
exprName=paste(path,"medulloblastoma/expression.txt",sep="");
exprName2=paste(path,"medulloblastoma/sampleID_to_subtype.txt",sep="");
resultdir=paste(path,"medulloblastoma/results/result",phenocodesbin,"_",analysistype,"_",indpermut,"/",sep="");dir.create(resultdir);
anno=read.table(paste(datadir,"snpsannot.txt",sep=""),sep="\t")
class(anno[,3])<-"character";class(anno[,4])<-"character";
nbgenes=max(anno[,9]);
nbsnps=c();
for (i in 2:length(anno[,10]))if(anno[i,10]==1)nbsnps=c(nbsnps,anno[i-1,10]);
nbsnps=c(nbsnps,anno[length(anno[,10]),10]);

genenames=unique(anno[,1])
phenonames=c("SHH","Group3","Group4")
phenosizes=c(31,21,24);
pheno=c(rep(phenocodes[1],phenosizes[1]),rep(phenocodes[2],phenosizes[2]),rep(phenocodes[3],phenosizes[3])) # 1 sick ; 0 control
nbpatients=length(pheno);


library(hash)
library(preprocessCore);

harm=list(); length(harm) <- nbgenes;
for(i in 1:nbgenes)harm[[i]]<- anno[which(anno[,9]==i),5];
harmgene=rep(meancgenes/nbgenes,nbgenes)
##

#create hashtable snp to index in annotation
collapser=function(x){return(paste(x,collapse="-",sep=""))}
keys=apply(anno[,c(1,2,3,4)],1,collapser)
map =hash(keys=keys,values=1:length(anno[,1]))
exprmapping=read.table(exprName2,colClasses=c("character","character","character")) 


#Add gene interactions
library(Matrix)
source(paste(codeDir,"functions/load_network_functions.r",sep=""));
source(paste(codeDir,"functions/process_expr.r",sep=""));


net1=load_network_genes(networkName,as.character(genenames),maxConnections)$network;
#cliks=networkToCliques(net)
net2=mapply(surround2,net1,1:length(net1) , MoreArgs=list(net1=net1));
net=net1;if (usenet2)net=net2;

#loading genotype
if (analysistype %in% c(0,2,3)){
het=list();length(het)<- nbgenes;for (i in 1:nbgenes){het[[i]]<- list(); length(het[[i]])<- nbpatients;} 
hom=list();length(hom)<- nbgenes;for (i in 1:nbgenes){hom[[i]]<- list(); length(hom[[i]])<- nbpatients;} #JUST FIXED IT

ind=0;sampleIDs=c();filetype="MB";
for (k in 1:length(phenonames)){
files=list.files(path=paste(datadir,phenonames[k],"/",sep=""), pattern=paste("(",filetype,").*\\.txt$",sep=""))
sampleIDs=c(sampleIDs,sapply(files,function(x)paste(filetype,".",substring(unlist(strsplit(x,paste(filetype,"-",sep="")))[2],1,4),sep="")))
for (i in 1:length(files)){
  gen=read.table(paste(datadir,phenonames[k],"/",files[i],sep=""),colClasses=c("character","NULL","numeric","NULL","character","character","character",rep("NULL",4)),sep='\t',skip=1)
  clp=apply(gen[,c(1,3,4,5)],1,collapser)
  for (j in 1:(dim(gen)[1])){
	x=map[[clp[j]]];
	if (gen[j,2]==1)het[[anno[x,9]]][[i+ind]]=c(het[[anno[x,9]]][[i+ind]],anno[x,10])
	if (gen[j,2]==2)hom[[anno[x,9]]][[i+ind]]=c(hom[[anno[x,9]]][[i+ind]],anno[x,10])
  }#print(i)
}
ind=ind+phenosizes[k];
}
}

get_pID=function(x){y=exprmapping[which(exprmapping[,3]==unlist(strsplit(x,"_"))[1]), 1];if (length(y))return(y);return(paste("BL-",substring(unlist(strsplit(x,paste("BL-",sep="")))[2],1,4),sep=""));}

if (analysistype>0){
het2=list();length(het2)<- nbgenes;for (i in 1:nbgenes){het2[[i]]<- list(); length(het2[[i]])<- nbpatients;} 
hom2=list();length(hom2)<- nbgenes;for (i in 1:nbgenes){hom2[[i]]<- list(); length(hom2[[i]])<- nbpatients;} 
ind=0;sampleIDs2=c();filetype="BL";
for (k in 1:length(phenonames)){
files=list.files(path=paste(datadir,phenonames[k],"/",sep=""), pattern=paste("(",filetype,").*\\.txt$",sep=""))
#sampleIDs2=c(sampleIDs2,sapply(files,function(x)paste("MB",".",substring(unlist(strsplit(x,paste(filetype,"-",sep="")))[2],1,4),sep="")))
sampleIDs2=c(sampleIDs2,sapply(files,function(x)paste("MB",".",substring(get_pID(x),4,7),sep="") ))
for (i in 1:length(files)){
  gen=read.table(paste(datadir,phenonames[k],"/",files[i],sep=""),colClasses=c("character","NULL","numeric","NULL","character","character","character",rep("NULL",4)),sep='\t',skip=1)
  clp=apply(gen[,c(1,3,4,5)],1,collapser)
  for (j in 1:(dim(gen)[1])){
	x=map[[clp[j]]];
	if (gen[j,2]==1)het2[[anno[x,9]]][[i+ind]]=c(het2[[anno[x,9]]][[i+ind]],anno[x,10])
	if (gen[j,2]==2)hom2[[anno[x,9]]][[i+ind]]=c(hom2[[anno[x,9]]][[i+ind]],anno[x,10])
  }
}
ind=ind+phenosizes[k];
}
}

if (analysistype>1){
het3=list();length(het3)<- nbgenes;for (i in 1:nbgenes){het3[[i]]<- list(); length(het3[[i]])<- nbpatients;} 
hom3=list();length(hom3)<- nbgenes;for (i in 1:nbgenes){hom3[[i]]<- list(); length(hom3[[i]])<- nbpatients;} 
for (i in 1:nbgenes){for(j in 1:nbpatients){
if (length(het[[i]][[j]])>0) het3[[i]][[j]]<- het[[i]][[j]][which(!(het[[i]][[j]] %in% het2[[i]][[j]]))];
if (length(hom[[i]][[j]])>0) {
hom3[[i]][[j]]<- hom[[i]][[j]][which(!(hom[[i]][[j]] %in% c(het2[[i]][[j]],hom2[[i]][[j]])))];
if (analysistype==2) het3[[i]][[j]]<- sort(unique(c(het3[[i]][[j]], hom[[i]][[j]][which(!(hom[[i]][[j]] %in% hom2[[i]][[j]]) & (hom[[i]][[j]] %in% het2[[i]][[j]]))])));
}}}
het=het3;hom=hom3;
} else if (analysistype==1){het=het2;hom=hom2;sampleIDs=sampleIDs2}

exprbad=read.table(exprName);#WARNING: Need to be the same ordering as exprmapping else BUG. Correcting now by manually reordering:
ord=match(substring(exprmapping[,1],4),substring(colnames(exprbad),4))
expr=exprbad[,ord];
expr2=expr;expr2=normalize.quantiles(as.matrix(expr));colnames(expr2)<- colnames(expr);rownames(expr2)<- rownames(expr);
uniqsub=unique(exprmapping[which(exprmapping[,2] %in% phenonames),2]);
e1=matrix(0,dim(expr2)[1],dim(expr2)[2])
for (i in 1:length(uniqsub)){
ind=which(exprmapping[,2]==uniqsub[i]);
e1[,ind]=t(apply(expr2[,ind],1,medmad));#e=t(apply(patExpr,1,modmad));
}
#colnames(e)<- NULL;
matchexpr=match(sampleIDs,colnames(expr2));
matchgenes=match(genenames,rownames(expr2))
e=matrix(0,nbgenes,nbpatients)
e[which(!is.na(matchgenes)),which(!is.na(matchexpr))]=e1[matchgenes[which(!is.na(matchgenes))],matchexpr[which(!is.na(matchexpr))]]

#Warning: Sometimes snps are repeated in patients files and they are present multiple times in het , hom
#This is often related to frameshift mutations who are ambiguously annotated (multiple times)

source(paste(codeDir,"functions/analyse_results.r",sep=""));
source(paste(codeDir,"newmethod/sumproductmem.r",sep=""))

#Rprof(filename = "Rprof.out")
ptm <- proc.time()
x<- grid_search(codeDir,nbgenes,nbpatients,nbsnps,harm,harmgene,meancgenes,complexityDistr,pheno,hom,het,net,e, cores,pcontrol,alpha,netparams,removeExpressionOnly);
print(proc.time()-ptm)
#summaryRprof(filename = "Rprof.out")


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
pie_plot_patients(x$causes[[7]],bestgenes,genenames,resultdir,TRUE)
}

#Permutations
xp=list();length(xp)=npermut;
lik=rep(0,npermut);acc=rep(0,npermut);
if(npermut){ for (i in 1:npermut){
phenop=sample(pheno);
xp[[i]]=grid_search(codeDir,nbgenes,nbpatients,nbsnps,harm,harmgene,meancgenes,complexityDistr,phenop,hom,het,net,e, cores,pcontrol,alpha,netparams,removeExpressionOnly);
plot_graphs(xp[[i]],phenop,resultdir,NULL,i,reorder=TRUE);
lik[i]=sum(xp[[i]]$likelihood);
acc[i]=(t(xp[[i]]$predict-pheno)%*%(xp[[i]]$predict-pheno));
}}


write.table(c(acc0,acc),paste(resultdir,"pred.txt",sep=""),row.names=FALSE,col.names=FALSE)
write.table(c(lik0,lik),paste(resultdir,"lik.txt",sep=""),row.names=FALSE,col.names=FALSE)

#check memory footprint
#sort( sapply(ls(),function(x){object.size(get(x))}))
#object.size(x=lapply(ls(), get));print(object.size(x=lapply(ls(), get)), units="Mb")

