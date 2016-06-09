simul_sample= function(x ,vec,n, sd){
nt=rnorm(1,n,sd); #could be normal around n
return(sample(vec,nt))
}

simul_unif= function(x,start=0,end=1){
return(runif(length(x),start,end ))
}

#simulate phenotype
simul_pheno= function(l1,l2,l3,model,sdamage,convertToGene,allmap,dict){
sToPhenoCoef=50; #I can either count damaging mutation or add them as min(50*s,1)
effects1=0;effects2=0;
ngs=convertToGene[l1];
ngsInModel=which(ngs %in% model)
if (length(ngsInModel)>0){
      x=which(allmap[dict[l1[ngsInModel]],3]>sdamage)
      if (length(x)>0){ 
      effects1=sum( apply( cbind(50*allmap[dict[l1[ngsInModel[x]]],3],1), 1, min) ); #Was length(y) #Add snps related causes in model     
      }
}    
expInModel=which(l2 %in% model)
if (length(expInModel)>0) effects2=sum(l3[expInModel]) ;#Was length(expInModel); #Add expression related causes in model 
return(c(effects1,effects2))
}

simul_pheno_number= function(l1,l2,l3,model,sdamage,convertToGene,allmap,dict){
effects1=0;effects2=0;y1=rep(0,length(model));y2=y1;
ngs=convertToGene[l1];
ngsInModel=which(ngs %in% model)
if (length(ngsInModel)>0){
      x=which(allmap[dict[l1[ngsInModel]],3]>sdamage);
      if (length(x)>0){ 
      effects1=length(x);
      }
}    
expInModel=which(l2 %in% model)
if (length(expInModel)>0) {
effects2=length(which(l3[expInModel]<3))+2*length(which(l3[expInModel]>=3) ) ;
}
return(c(effects1,effects2))
}

simul_pheno_number_Qual= function(l1,model,sdamage,convertToGene,allmap,dict){
y1=rep(0,length(model));
ngs=convertToGene[l1];
ngsInModel=which(ngs %in% model)
if (length(ngsInModel)>0){
      x=which(allmap[dict[l1[ngsInModel]],3]>sdamage);
      for (i in 1:length(model))y1[i]=length(which(ngs[ngsInModel[x]]==model[i]));
}    
return(y1)
}

simul_pheno_number_Quant= function(l2,l3,model,sdamage,convertToGene,allmap,dict){
y2=rep(0,length(model));
expInModel=which(l2 %in% model)
if (length(expInModel)>0) {
effects2=length(which(l3[expInModel]<3))+2*length(which(l3[expInModel]>=3) ) ;
for (i in 1:length(model))y2[i]=length(which(l2[expInModel]==model[i]));
}
return(y2)
}

#Number of damaging snps per patient
simul_nbDamage=function(l1,sdamage,allmap,dict){
x=which(allmap[dict[l1],3]>sdamage)
return(length(x))
}

#get zygocyty from list of snps 
#l11 is unique(l1)
getzygocity=function(l1,l11){
zyg=rep(1,length(l11))
if (length(l11)<length(l1)){
x=l1[which(duplicated(l1))]
zyg[which(l11 %in%x )]=2;
}
return(zyg)
}



#From a patient snps list to a matrix describing them
simul_indivToMatrix= function(l1,convertToGene,ps,genes,dict){
l11=unique(l1)
zyg=getzygocity(l1,l11);
mat=data.frame(l11,genes[convertToGene[l11]],ps[dict[l1[match(l11,l1)]]],zyg)
return(mat)
}

#create annotation for all snps
simul_annotMatrix= function(l1,convertToGene,ps,genes){
mat=data.frame(genes[convertToGene[l1]],l1,ps)
return(mat)
}

#Form a 0,1,2 matrix for all patients and controls snps
simul_snpsToMatrix= function(snps,indPatients,indControls,convertToGene,ps,genes,dict){
allindiv=c(indPatients,indControls);
l=c();
for (i in allindiv)l=c(l,snps[[i]]);
l1=sort(unique(l));
dict2=match(1:length(convertToGene), l1);
l1gene=convertToGene[l1];
mat=matrix(0,length(allindiv),length(l1));
k=1;
for (i in allindiv){
l11=unique(snps[[i]])
zyg=getzygocity(snps[[i]],l11);
mat[k, dict2[l11] ]=zyg;
k=k+1;
}
return(list(mat=mat, gene=l1gene,snps=l1));
}


simul_expr_background=function(healthy_expr,nbGenes,nbPatients,nbControls,sampling=TRUE ){
patient_gene_expr=matrix(0,nbGenes,nbPatients);
control_gene_expr=matrix(0,nbGenes,nbControls);
for (i in 1:nbGenes){
bw=sd(healthy_expr[i,]);
m=mean(healthy_expr[i,]);
if(sampling){
bootstrap=sample(healthy_expr[i,],nbPatients+nbControls,replace=TRUE);
noise=rnorm(nbPatients+nbControls,0,bw/10);#I am using noice equal to one 10th of the standard deviation(to avoid adding a lot of noisse))
sig=bootstrap+noise;
}else{
sig=rnorm(nbPatients+nbControls,m,bw);
}
patient_gene_expr[i,] =(sig)[1:nbPatients];
control_gene_expr[i,] =(sig)[(nbPatients+1):(nbPatients+nbControls)];
}
#patient_gene_expr[patient_gene_expr<0] <- 0; control_gene_expr[control_gene_expr<0] <- 0; don't need this since I work in log space
return(list(pat_expr=patient_gene_expr, ctr_expr=control_gene_expr))
}

simul_expr_change=function(ret,indPatients,indControls,allexp,allexpMagnitude,coef){
  pat=ret$pat_expr;ctr=ret$ctr_expr;
  ng=dim(pat)[1];
  x=runif(ng)-0.5;
  perturb_direction=x/abs(x);
  sd_expp=apply(ret$pat_expr,1,sd);
  sd_expc=apply(ret$ctr_expr,1,sd);
  for (i in 1:length(indPatients)){
  ind=indPatients[i];
  pat[allexp[[ind]],i]=pat[allexp[[ind]],i]+ perturb_direction[allexp[[ind]]]*allexpMagnitude[[ind]]*coef*sd_expp[allexp[[ind]]];
  }
  for (i in 1:length(indControls)){
  ind=indControls[i];
  ctr[allexp[[ind]],i]=ctr[allexp[[ind]],i]+ perturb_direction[allexp[[ind]]]*allexpMagnitude[[ind]]*coef*sd_expc[allexp[[ind]]];
  }
  #remove negative values
  #pat[pat<0]=0;ctr[ctr<0]=0;# not needed since work in log
return(list(pat_expr=pat,ctr_expr=ctr))
}

unbiasRatio=function(modelSize){
if (modelSize<3)return(3.5);
if (modelSize==3)return(2);
if (modelSize==4)return(1.4);
if (modelSize==5)return(1.1);
if (modelSize==6)return(0.85);
if (modelSize>6)return(5/modelSize);
}

