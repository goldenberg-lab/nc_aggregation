#awk -F"\t" '{print $6"\t"$7"\t"$9"\t"$10"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$40"\t"$42"\t"$46"\t"$35"\t"$36}' Downloads/TOFdata/TOF_ALL_CALLS.txt  |uniq > Downloads/TOFdata/variants0.txt
#awk -F"\t" '{if (NR>1){ if ($9==""){split($8,a,"("); split(a[1],b,";");$9=b[1];} if ($5~/splicing/)$11="0.80001"; if ($7~/frame/ || $7~/stop/)$11=0.99; if($11)print $9"\t"$1"_"$2"_"$3"_"$4"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$7} }' Downloads/TOFdata/variants0.txt | sort | uniq >Downloads/TOFdata/variants1.txt
#annot=read.table("/home/aziz/Downloads/TOFdata/variants1.txt",sep="\t",stringsAsFactors =FALSE);
#genes=unique(annot[,1]); g=rep(0,length(genes));
#u=rep(1,nrow(annot));k=1;j=1; for (i in 2:nrow(annot)){if(annot[i-1,1]==annot[i,1]){k=k+1;} else {g[j]=k;j=j+1;k=1;}; u[i]=k;}; g[j]=k;
#write.table(data.frame( genes, 1:length(genes), g),"/home/aziz/Downloads/TOFdata/genes.txt",sep="\t",col.names=FALSE,row.names=FALSE)
#write.table(data.frame(annot[,c(1,2,8)],annot[,4], rowMeans(annot[,6:7],na.rm=TRUE),u ),"/home/aziz/Downloads/TOFdata/variants.txt",sep="\t",col.names=FALSE,row.names=FALSE)#rowMeans(annot[,3:5],na.rm=TRUE)
#awk -F"\t" '{if (NR>1){a=0;if($3=="CASE")a=1;print $1"\t"a"\t"$2}}' Downloads/TOFdata/TOF_ALL_CALLS.txt  |sort | uniq > Downloads/TOFdata/patients.txt
#awk -F"\t" '{if (NR>1){ if ($21~/splicing/ || $23~/frame/ || $23~/stop/ || ($42)){a=1;if($14=="Hom")a=2; print $6"_"$7"_"$9"_"$10"\t"$1"\t"a"\t"$14"\t"$15   }}}' Downloads/TOFdata/TOF_ALL_CALLS.txt > Downloads/TOFdata/genotype.txt

#awk '{print $1; if ($3==2)print $1}' Downloads/TOFdata/genotype.txt  | sort | uniq -c >  Downloads/TOFdata/genotype_counts.txt
#rmaftab=read.table("Downloads/TOFdata/genotype_counts.txt"); rmaf=rmaftab[,1]/ (2*93); names(rmaf) <- rmaftab[,2];
#annot2=read.table("/home/aziz/Downloads/TOFdata/variants.txt",sep="\t",stringsAsFactors =FALSE);annot2[which(is.na(annot2[,5])),5]=0;
#rmaflocal= rmaf[annot2[,2]];indmafgood=which( annot2[,5]< 0.05 & rmaflocal<0.05 );
#write.table(annot2[indmafgood,],"Downloads/TOFdata/variants3.txt",sep="\t",col.names=FALSE,row.names=FALSE);

#for harmful only
#write.table(data.frame(annot[,c(1,2,8)],annot[,4], rowMeans(annot[,6:7],na.rm=TRUE),u )[which(annot[,4]>0.85 & annot[,3]>0.95),],"/home/aziz/Downloads/TOFdata/variantsharm.txt",sep="\t",col.names=FALSE,row.names=FALSE)#rowMeans(annot[,3:5],na.rm=TRUE)
#repeat last block

path=commandArgs(TRUE)[1];
datapath=as.numeric(commandArgs(TRUE)[2]);
resultdir=as.numeric(commandArgs(TRUE)[3]);
indpermut=as.numeric(commandArgs(TRUE)[4]);
npermut=as.numeric(commandArgs(TRUE)[5]);
genefile="genes.txt";patfile="patients.txt";varfile="variants.txt";exomefile="genotype.txt";exprfile="gene_expression.txt";

local=TRUE;

if (local){
path="/home/aziz/Desktop/aziz/diseaseMechanism/";
indpermut=0;npermut=5;
datapath="/home/aziz/Downloads/TOFdata/"
resultdir=paste(path,"Results/",sep="");
genefile="genes.txt";patfile="patients.txt";varfile="variantsharm4.txt";exomefile="genotype.txt";exprfile="gene_expression.txt";
}

#Model parameters
meancgenes=10;
complexityDistr=c(0.33,0.33,0.33);meancgenesppat=0;
pcontrol=0.95;
alpha=0.4;
usenet2=TRUE;
maxConnections=50;
netparams=c(0.9,0.01,0.01,2);
removeExpressionOnly=TRUE;
ksi=1;#log(ksi+x) to transform expression data
e=NULL;cores=1;
auto_balance=FALSE;

#networkName=paste(path,"networks/BIOGRID3.2.98.tab2/interactions.txt",sep="") ## Biogrid network
#networkName=paste(path,"networks/HPRD_Release9_062910/interactions.txt",sep="") ## HPRD network
networkName=paste(path,"networks/humannet/interactions150.txt",sep="");## HumanNet
#TODO change to whatever gene network you want to use. Can use other annotations.


codeDir=paste(path,"Rproject/",sep="");
dir.create(resultdir);


#load your genes annotation (Description in README file) 
genemat=read.table(paste(datapath,genefile,sep=""),sep="\t",stringsAsFactors =FALSE);

#load your patients file
phenomat=read.table(paste(datapath,patfile,sep=""),sep="\t",stringsAsFactors =FALSE);
phenomat=phenomat[which(phenomat[,3]!="DILATION"),]#CHANGED NEW line

geneids=genemat[,2];names(geneids)<- genemat[,1];
patientids=1:nrow(phenomat);names(patientids)<- phenomat[,1]; #patientids=phenomat[,2]; names(patientids)<- phenomat[,1];#CHANGED
nbsnps=genemat[,3];
pheno=phenomat[,2]#pheno=phenomat[,3];#CHANGED
#dys=which(sapply(as.character(phenomat[,3]),function(x)grepl("DYS",x)));pheno[dys]=1;pheno[-dys]=0;#CHANGED NEW LINE
#pheno=1-pheno;

#load your exome variants annotation
annot=read.table(paste(datapath,varfile,sep=""),sep="\t",stringsAsFactors =FALSE);
varids=1:nrow(annot);names(varids) <- annot[,2];


#Balance them by reducing the bigger(sampling or taking first n)
ph1=which(pheno==1);ph0=which(pheno==0);
if (auto_balance){
if(length(ph1)>length(ph0)){ph1=sample(ph1,length(ph0));
}else if (length(ph0)>length(ph1))ph0=sample(ph0,length(ph1));
pheno=c(rep(1,length(ph1)),rep(0,length(ph0)))
includedpatientids=patientids[c(ph1,ph0)];
ph1=1:length(ph1); ph0=(length(ph1)+1):(length(ph1)+length(ph0));
mappat=match(patientids,includedpatientids);
} else {mappat=1:length(patientids);}

#Inputs
nbgenes=length(geneids); nbpatients=length(pheno); genenames= names(geneids);
harm <- list(); length(harm) <- nbgenes;for (i in 1:nbgenes)harm[[i]] <- rep(0,nbsnps[i]);

gids=geneids[annot[,1]];#gene index
for (i in 1: nrow(annot))harm[[ gids[i] ]][annot[i,6]] =annot[i,4];
harmgene=rep(meancgenes/nbgenes,nbgenes);#CHANGED
harmgene[bestgenesnow]=0.5; #NEW LINE TO REMOVE; TRY TO OVERFIT


library(preprocessCore);
library(Matrix)
source(paste(codeDir,"functions/load_network_functions.r",sep=""));
source(paste(codeDir,"functions/process_expr.r",sep=""));

#load expression
exprmat=read.table(paste(datapath,exprfile,sep=""),sep="\t");
mapped=mappat[patientids[colnames(exprmat)]];#patient in expression data mapped to includedpatientids (or patientsids if mappat identity);
cases=match(ph1,mapped);
controls=match(ph0,mapped);
####log transform and quantile normalize each group
exprcases=log(ksi+normalize.quantiles(as.matrix(exprmat)[,cases]));
exprcontr=log(ksi+normalize.quantiles(as.matrix(exprmat)[,controls]));
####Robust standardization
e1=matrix(0,nrow(exprmat),nbpatients);rownames(e1) <- rownames(exprmat);
e1[,ph1]=t(apply(exprcases,1,medmad));
e1[,ph0]=t(apply(exprcontr,1,medmad));
mapexprToGene=match(rownames(exprmat),genenames);
e=matrix(0,nbgenes,nbpatients);
e[mapexprToGene,]=e1;


#load genotypes (list of variants in each inividual and zygocity)
genotype=read.table(paste(datapath,exomefile,sep=""),sep="\t",stringsAsFactors =FALSE);
het=list();length(het)<- nbgenes;for (i in 1:nbgenes){het[[i]]<- list(); length(het[[i]])<- nbpatients;} 
hom=list();length(hom)<- nbgenes;for (i in 1:nbgenes){hom[[i]]<- list(); length(hom[[i]])<- nbpatients;} 
p=mappat[patientids[genotype[,2]]];ind=varids[genotype[,1]];g=geneids[annot[ind,1]];
for (i in 1:nrow(genotype)){
  if(!is.na(g[i]) && !is.na(p[i])){
    if (genotype[i,3]==1)het[[g[i]]][[p[i]]] <- c(het[[g[i]]][[p[i]]], annot[ind[i],6])
    if (genotype[i,3]==2)hom[[g[i]]][[p[i]]] <- c(hom[[g[i]]][[p[i]]], annot[ind[i],6])
  }
}

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
x<- grid_search(codeDir,nbgenes,nbpatients,nbsnps,harm,harmgene,meancgenes,complexityDistr,pheno,hom,het,net,e, cores,pcontrol,alpha,netparams,removeExpressionOnly);
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
pie_plot_patients(x$causes[[7]],bestgenes,genenames,resultdir,TRUE)
}
if (npermut>0)indpermut=indpermut+1
}


#Permutations
if (indpermut & npermut){
xp=list();length(xp)=npermut;
lik=rep(0,npermut);acc=rep(0,npermut);
if(npermut){ for (i in indpermut:(indpermut+npermut-1)){
phenop=sample(pheno);
xp[[i]]=grid_search(codeDir,nbgenes,nbpatients,nbsnps,harm,harmgene,meancgenes,complexityDistr,phenop,hom,het,net,e, cores,pcontrol,alpha,netparams,removeExpressionOnly);
plot_graphs(xp[[i]],phenop,resultdir,NULL,i,reorder=TRUE);
lik[i]=sum(xp[[i]]$likelihood);
acc[i]=(t(xp[[i]]$predict-pheno)%*%(xp[[i]]$predict-pheno));
}}
}

write.table(c(acc0,acc),paste(resultdir,"pred.txt",sep=""),row.names=FALSE,col.names=FALSE)
write.table(c(lik0,lik),paste(resultdir,"lik.txt",sep=""),row.names=FALSE,col.names=FALSE)

#Imbalanced Analyse:
HSPG2"   "PLEC"    "DNAH7"   "ADAM21"  "COL4A3"  "C2CD3"   "LRP1"    "EXOC3L"  "MYOM2"   "CELA1"   "GUCY2F"  "DCHS2"   "EXOSC6"  "VCAN"   "PYGL"    "LPA"     "ALPK1"   "COL18A1" "APC2"    "MSR1
-3rd gene probably false positive because of imbalance
-1st and second have odd ratio of 2 (if we correct for imbalance)
-ADAM21 have 10 patients with harmful snv (9 have rs72735759), no controls
-COL4A3 have at least 7 patients with potentially harmful mutations, no controls
-C2CD3 at least 10 patients with harmful mutations, 1 control
-LRP1 9 patients, 1 control
-EXOC3L 10 patients, 2 controls
-MYOM2 Odd ratio of 2.4 become less after correction for imbalance (5 controls, 12 patients)
-CELA1 9 patients no controls, (sets of 4 and 4 patients have same snvs : very close (6 nt) not in LD a priori)
-GUCY2F 3 heterozygous patients, 2 homozygous, no controls 
-DCHS2 10 cases, 2 controls
-EXOSC6 at least 6 patients (6 with the same snv), no control
-VCAN   8 patients 3 controls
-PYGL 10 cases, 2 controls
-LPA 8 cases, 2 controls
-ALPK1 10 cases, 2 controls
-COL18A1 10 cases (1 homozygous), 3 controls (at least for both)
-APC2 6 cases , 0 harmful controls
-MSR1 6 cases, 1 control

#imbalanced with MAF 0.05
"VCAN"    "HSPG2"   "GOLGA3"  "USH2A"   "BDP1"    "C2CD3"   "COL4A3"  "ADAM21"  "LRP1"    "SYNE1"   "PODXL"   "COL6A1"  "PLEKHG2" "PCSK4"   "PTPRH"   "ELSPBP1" "ZC3H3"   "PLEC"    "EPHA5"   "DDX58"
-16 vs 3
-22 vs 10
-15 vs 1
-21 vs 6
-12 vs 3
-16 vs 2
-7 vs 0
-9 vs 1
-8 vs 1
-18 vs 6
-10 vs 1
-14(1 homo) vs 3


#GOrilla analysis: extracellular matrix/structure organization, glomerular basement membrane development, lipoprotein metabolic process

#Analysis without Dilation MAF 0.05
"BCLAF1" "PCSK4"  "PTPRH"  "EPHA5"  "DDX58"  "VCAN"   "BCAM"   "AMBN"   "BTN3A2" "PODXL"	"DEFA4"    "SLC22A24" "FGGY"     "TNS1"     "OR52E2"   "TTC21A"   "RIF1"     "EXOSC6"  "C15orf40" "CCDC34" 
-10 vs 2
-8 vs 2

#Analysis without Dilation rare
"EXOSC6"  "RIF1"    "MYOM2"   "GUCY2F"  "COL18A1" "C2CD3"   "DCHS2"  "VCAN"    "ADAM21"  "EXOC3L2" "HSPG2"    "DDX58"   "MYO5B"    "ALMS1"    "MYO1E"    "HTT"      "KIAA0562" "MXRA5"   "APC2"     "LTF" 
-5 to 0
-6 to 0
-VCAN (6 to 3), ADAM21 5 to 1, HSPG2 10 to 7




#harmful only
"FAT4"      "VCAN"      "LAMA5"     "C2CD3"     "DNAH3"     "RP1L1"    "ALMS1"     "PDE6A"     "PHKB"      "BCR"       "ARHGAP33"  "SBF1"    "ROBO4"     "MYO5B"     "DENND2C"   "HSPG2"     "C14orf102" "MYOM2"    "NRP2"      "TNS1" 

"VCAN"     "PLEKHG2"  "DDX58"    "ADAM21"   "BDP1"     "SYNE1"   "ZC3H3"    "SERPINA1" "SHROOM3"  "FGGY"     "PCSK4"    "APBA3"   "C2CD3"    "ZNF19"    "SIGLEC12" "ZSWIM4"   "LRRC56"   "AQP7"     "IL10RA"   "SNCAIP" 

"ZSWIM4"    "ZNF434"    "LMAN1L"    "SELE"      "FREM2"     "PCNT"      "ZC3H3"     "MOS"       "LAMC2"     "FARP2"     "PRR15"     "ADAM21"   "PRKRA"     "OR4L1"     "BHLHE22"   "C20orf151" "VCAN"      "GPR151"    "SAT2"      "INADL"
#harmful only no dilation
"ALMS1"     "FAT4"      "VCAN"      "SBF1"      "C2CD3"     "ONECUT1"  "SHROOM3"   "WT1"       "MYO5B"     "TSHZ1"     "RIF1"      "TEP1"     "BCR"       "DNAH9"     "KIAA1409"  "TNS1"      "C14orf102" "FBN3"     "DLEC1"     "SLC43A2" 

"DDX58"    "VCAN"     "PCSK4"    "SHROOM3"  "C15orf40" "ADAM21"  "FARP1"    "FGGY"     "TEP1"     "EPB41L2"  "AMBN"     "SLC8A3"   "C2CD3"    "ZC3H3"    "ACOXL"    "GUCY2F"   "PTPRH"    "DCAKD"   "ZNF19"    "PIK3C2B"

"USP8"      "ARSH"      "ADAM21"    "OR7C1"     "TRIM22"    "LMAN1L"    "MAST2"     "PRR15"     "ROS1"      "LAMC2"     "VCAN"      "CD6"       "ZNF434"    "ARSA"      "PCSK4"     "LRRC6"     "C20orf151" "C15orf52"  "BHLHE22"   "VPS53"
