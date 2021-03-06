#!/usr/bin/env Rscript
path=commandArgs(TRUE)[1];
nb_genes_causal=as.numeric(commandArgs(TRUE)[2]);
ratioExprSnp=as.numeric(commandArgs(TRUE)[3]);
sample_size=as.numeric(commandArgs(TRUE)[4]);
cores=as.numeric(commandArgs(TRUE)[5]);
indrep=as.numeric(commandArgs(TRUE)[6]);
nbrep=as.numeric(commandArgs(TRUE)[7]);
freeMem=FALSE;

local=FALSE;

if(local){
nb_genes_causal=10;
sample_size=100;
ratioExprSnp=0;#c(0,0.66,1,1.5,4) for 100%, 60%,50%,40%,20% snps
path="/home/aziz/Desktop/aziz/diseaseMechanism/";
#path="/project/lipidext/diseaseMechanism/";
cores=1;freeMem=TRUE;nbrep=1;indrep=1;
}
outputdirSuffix=paste(nb_genes_causal,"_",ratioExprSnp,"_",sample_size,sep="");
dir=paste(path,"simul/simulrarecommon/simul_",outputdirSuffix,"/",sep="");#output directory for the simulations
dir.create(dir, recursive=TRUE);
codeDir=paste(path,"Rproject/",sep="");

#simul param
simgenelengths=FALSE;
simmaxmut=10;#max mutation per patients (kill all individuals above)
nbPatients=3000;
ratioSignalSim=1;
maxConnections=50;
maxneighborhoodsize=max(100,2*nb_genes_causal);#max second order neighberhood for selecting model
#Analysis param

npermut=0;
ratioSignal=ratioSignalSim;
complexityDistr=c(0,1,0);
decay=c(0.05,0.1,0.2,0.4);
removeCommon=FALSE;
alpha=0.5;if(sample_size>=500)alpha=0.35;
usenet2=TRUE;removeExpressionOnly=FALSE;
netparam=0.9; #parameter by default in the MRF
netmax=0.01;#max proba that can be originating from network (modulate strength of gene network messages)
netmaxNoise=0.01;#Correct for high degree nodes where accumulation of small signals explodes; max proba originating for accumulation of priors
netrelaxed=1; #maximal net contribution computed on the basis that 1/1+x of the true genes are in the neighbourhood (this means stronger net signals in general)
propagate=FALSE;
methods=c("CAST","CALPHA","VT","SKAT-O");
testExpr=0;testdmgwas=1; #0 if not doing diff. expr , one if we do. Same for doing or not doing dmGWAS.
testVariations=1;#=5; #try without net, without expression or without harmfulness (add harm first, then add net or add expr)
testrobust=0; #0 nothing, 1 network/harm/expr, 2 meancgenes, 3 ratioSignal, 4 complexitydistr/decay
thresh=c(0.5);

library(preprocessCore);
library(Matrix)
#library(modeest)#For mode computations
#library(pROC)#for ROC curve and AUC measurement
#library(robustbase) #for Sn and Qn
#library(asbio)#for biweght midvariance
source(paste(codeDir,"functions/harmonize_expr.r",sep=""));
source(paste(codeDir,"functions/simul_polyphen.r",sep=""));
source(paste(codeDir,"functions/simul_functions.r",sep=""));
source(paste(codeDir,"functions/load_network_functions.r",sep=""));
source(paste(codeDir,"functions/analyse_results.r",sep=""));

exprName=paste(path,"simul_data/healthyExpr65834.txt",sep="");
#networkName=paste(path,"networks/BIOGRID3.2.98.tab2/interactions.txt",sep="");
#networkName=paste(path,"networks/HPRD_Release9_062910/interactions.txt",sep="");## HPRD network
networkName=paste(path,"networks/BIOGRID3.4-132/Biogrid-Human-34-132.txt",sep="")## Biogrid updated
#networkName=paste(path,"networks/humannet/interactions150.txt",sep="");## HumanNet
net=load_network_genes(networkName,NULL,maxConnections);
healthy_expr=as.matrix(read.table(exprName));
genes=intersect(net$genes, rownames(healthy_expr));
genesToExpr=match(genes,rownames(healthy_expr));
healthy_expr=healthy_expr[genesToExpr,];
nbgenes=length(genes);

#load SNPs
allmap=read.table(paste(path,"simul_data/outvar.var",sep=""),colClasses=rep("numeric",4))
sdamage=0.001;#threshold for harmfulness
codelength=max(allmap[,1])+1;
#Simulate polyphen scores
polyph_distr=read.table(paste(path,"simul_data/polyphen.txt",sep=""));
distr1=polyph_distr[which(!is.na(polyph_distr[,2]) & polyph_distr[,1]=="deleterious"),2];
distr2=polyph_distr[which(!is.na(polyph_distr[,2]) & polyph_distr[,1]=="neutral"),2];
ps=apply(cbind(allmap[,3]>sdamage,allmap[,3]),1,simul_polyphen2,distr1=distr1,distr2=distr2);
mf=allmap[,2];#MAF
if(simgenelengths){
lGene=sample(read.table(paste(path,"simul_data/CDS_lengths.txt",sep=""),header=TRUE)[,1],nbgenes);#I sample randomly from gene lengths. I could take the same genes as expression/net but overlap is not 100%
lGene2=round(codelength*lGene/sum(lGene));lGene2[length(lGene2)]=lGene2[length(lGene2)]+codelength-sum(lGene2);#modify length of last gene for coherence with total seq length
convertToGene=unlist(mapply(function(ig,jg)rep(ig,jg),1:length(lGene2),lGene2));
}else {lGene=codelength/nbgenes;convertToGene=1+((1:codelength) %/%lGene);}
if (indrep==1)write.table( simul_annotMatrix(allmap[,1],convertToGene,ps,genes),paste(dir,"annot.txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=c(1,2));

fc <- file(paste(path,"simul_data/outpat.pat",sep=""))
allmut <- lapply(strsplit(readLines(fc), "\t"), as.numeric)
close(fc)
dict=match(1:codelength, allmap[,1]); #require allmap sorted # Avoid using many "which" later
nbDamage=unlist(lapply(allmut,simul_nbDamage,sdamage, allmap, dict))
allmutAvgLength= mean(nbDamage)
allmutSdLength= sd(nbDamage)


#source("functions/setexpressionratio.r");

#Simulate expression perturbations
corrector=0.9#unbiasRatio(nb_genes_causal);#ratio is 1 if model size 5; 0.5 md10; 3.5 md2; For 66% 1.3 md2/5
nbExprAvg=corrector*ratioExprSnp*allmutAvgLength;nbExprSd=corrector*ratioExprSnp*allmutSdLength
nbindiv=length(allmut)
allexp <-list();length(allexp) <-nbindiv;
allexpMagnitude <- list();length(allexpMagnitude) <-nbindiv;
genesID=1:nbgenes
allexp <-lapply(allexp,simul_sample,genesID,nbExprAvg,nbExprSd)
start=0.5;end=0.8;coef=6; # perturbation is unif(start,end)*coef*sd *sign
allexpMagnitude <- lapply(allexp,simul_unif,start=start,end=end)


#Here I can vary sample_size and maybe even nb_genes_causal without having to repeat what preceded

#reapeat experiment nbtry times with the same model
powermeasures=matrix(0,nbrep,4*testVariations+ 3*(testExpr+length(methods))+2* testdmgwas); #measure sensitivity and precision for this method + others
contributions=matrix(0,2,nbrep);
xres=list();length(xres) <- testVariations;
for (z in 1:nbrep){
dirrep=paste(dir,"rep",indrep+z-1,"/",sep='');dir.create(dirrep);
source(paste(codeDir,"main_Simul_Disease.r",sep=""))
print("Gene contributions");
print(rowSums(indivScores_Qual));print(rowSums(indivScores_Quant));
object.size(x=lapply(ls(), get));print(object.size(x=lapply(ls(), get)), units="Mb")
source(paste(codeDir,"main_Simul_Analyze.r",sep=""))


for (w in 1:testVariations){
to=which(order(xres[[w]]$h,decreasing=TRUE) %in% truth);found=which(xres[[w]]$h>thresh);
powermeasures[z,1+4*(w-1)]=xres[[w]]$status;#convergence status
powermeasures[z,2+4*(w-1)]=length(which(found %in% truth))/length(truth);#sensitivity
if (length(found)){powermeasures[z,3+4*(w-1)]=length(which(found %in% truth))/length(found);}else {powermeasures[z,3+4*(w-1)]=1;}#precision
powermeasures[z,4+4*(w-1)]=length(which(to<= length(truth)))/length(truth);#precision at rank
}
print("Method performance assessment: done");

save.image(paste(dirrep,"data.RData",sep=""))


if(length(methods)){
library(AssotesteR)
library(SKAT)
genesToTest=unique(trans$gene);
labels=c( rep(1,length(indPatients)),rep(0,length(indControls)));
for (w in 1: length(methods)){
pvals=other_methods(trans,labels,genesToTest,methods[w],1000000,mf,dict,removeCommon=removeCommon);
genesToTest1=genesToTest[which(apply(pvals,1,min)<0.05/nbgenes)];#8.10e-6 is the significance threshold after multiple hypothesis
powermeasures[z,(4*testVariations+3*(w-1)+1)]=length(which(genesToTest1 %in% truth))/length(truth) #sensitivity
if (length(genesToTest1)){powermeasures[z,(4*testVariations+3*(w-1)+2)]=length(which(genesToTest1 %in% truth))/length(genesToTest1);} else powermeasures[z,(4*testVariations+3*(w-1)+2)]=1; #precision
powermeasures[z,(4*testVariations+3*(w-1)+3)]=length(which(genesToTest[order(pvals)[1:length(truth)]] %in% truth))/length(truth) #precision at rank
}
print("Gene based tests performance assessment: done");
}

if (testExpr){#Differential expression
pvalsexp=sapply(1:nbgenes,function(x)t.test(expr[,affectedexpr], expr[,-affectedexpr])$p.value)
genesig=which(pvalsexp<0.05/nbgenes);
powermeasures[z,(4*testVariations+3*length(methods)+1)]=length(which(genesig %in% truth))/length(truth) #sensitivity
if (length(genesig))powermeasures[z,(4*testVariations+3*length(methods)+2)]=length(which(genesig %in% truth))/length(genesig) #precision
powermeasures[z,(4*testVariations+3*length(methods)+3)]=length(which(order(pvalsexp)[1:length(truth)] %in% truth))/length(truth) #precision at rank
print("Differential expression performance assessment: done");
}

if (testdmgwas){#dmGWAS (done on SKAT pvalues)
library(dmGWAS)
skatpvals=rep(0.5,nbgenes);skatpvals[genesToTest]=pvals;
skatpvals[which(skatpvals<=10^(-16))]=10^(-16);skatpvals[which(skatpvals>=1)]=0.5#+runif(length(which(skatpvals>=1)))/2;
d1=data.frame(gene=genenames,weight=skatpvals,stringsAsFactors =FALSE);
netmat=data.frame(interactorA=rep("gene",100000),interactorB=rep("gene",100000),stringsAsFactors =FALSE)
k=1;for (i in 1:nbgenes)if(length(net1[[i]])){for(j in net1[[i]]){if (j >i){netmat[k,1]=genenames[i];netmat[k,2]=genenames[j];k=k+1;}}}
netmat=netmat[1:(k-1),]
#resdmgwas2=run_dmgwas(net1,skatpvals,expr[,affectedexpr],expr[,-affectedexpr]);
#sel=chooseModule(resdmgwas,top=0.01,plot=FALSE)
resdmgwas=dms(netmat,d1,expr1=NULL,expr2=NULL,d=1,r=0.1)
sel=resdmgwas$genesets.clear[[as.numeric(rownames(resdmgwas$zi.ordered)[1])]]
powermeasures[z,(4*testVariations+3*(length(methods)+testExpr)+1)]=length(which(sel %in% genenames[truth]))/length(truth) #sensitivity
if (length(sel)){powermeasures[z,(4*testVariations+3*(length(methods)+testExpr)+2)]=length(which(sel %in% genenames[truth]))/length(sel)
} else powermeasures[z,(4*testVariations+3*(length(methods)+testExpr)+2)]= 1#precision
print("NAA method performance assessment: done");
}

}
#end repetitions

write.table(powermeasures,paste(dir,"performance",outputdirSuffix,"_",indrep,".txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE);

if(0){
  print(paste("***Stats:", mean(powermeasures[,1]),mean(powermeasures[,2])  ));
  allresults[[z1]]<- powermeasures;
  if (nbrep>1){
    perc=colSums(contributions); contributionsPercentage=100*contributions/rbind(perc,perc)
    print(paste("contribution averages: ", rowMeans(contributionsPercentage)));
  }
  source("functions/simul_plot.r")
  #results_form=simul_format_results(allresults,param1,param2)
  #simul_box(results_form,param1,param2,"")
  simul_mean_plot(allresults, methods,param1,param2,param3,"")
}

