#factor SNPs to Quality Generate Matrix
factorQualGen= function(maxsnps){
  factorQual=array(0,dim=c(2,maxsnps+1,3))
  for(i in 1:maxsnps)factorQual[2,i,]=c(0,0,1);
  factorQual[1,1,]=c(1,0,0);
  for (i in 2:maxsnps){
    probhet=(0.5)^(i-2)
    factorQual[1,i,]=c(0,probhet,1-probhet)
  }
  return(factorQual)
}

#factor SNPs to Quality # removable function
factorQualC= function(hett,homt,factorQual){
return(factorQual[1+ (homt>0), 1+hett, ]);
}

#factor variants to Quantity:
factorQuantGen= function(sizeF,ebins,relbins,ecrossbins=ebins){#1 Very Low only; 2 Very high only; 3 Lower only; 4 Higher only; 5 Low only; 6 High only 7 No response; 8 none
#ebin -2 is for the elements smaller than -2; ebin 0 for around 0; ebin 2 is for larger than 2
  mat=array(0,c(sizeF,3,length(ebins),length(relbins),2,length(ebins)));mat[,1,,,,]=1; #by default
  for (k in 1:2)for (l in 1:length(ecrossbins))for (i in which(ebins<= -3))mat[1,,i,,k,l]=factorQuantGenOne(ebins[i],relbins);
  for (k in 1:2)for (l in 1:length(ecrossbins))for (i in which(ebins>= 3))mat[2,,i,,k,l]=factorQuantGenOne(ebins[i],relbins);
 
  for (k in 1:2)for (l in 1:length(ecrossbins))for (i in which(ebins<= -2.5))mat[3,,i,,k,l]=factorQuantGenOne(ebins[i],relbins);
  for (k in 1:2)for (l in 1:length(ecrossbins))for (i in which(ebins>= 2.5))mat[4,,i,,k,l]=factorQuantGenOne(ebins[i],relbins);

  for (k in 1:2)for (l in 1:length(ecrossbins))for (i in which(ebins<= -2))mat[5,,i,,k,l]=factorQuantGenOne(ebins[i],relbins);
  for (k in 1:2)for (l in 1:length(ecrossbins))for (i in which(ebins>= 2))mat[6,,i,,k,l]=factorQuantGenOne(ebins[i],relbins);
  return(mat)
}

factorQuantGenOne= function(e,relbins){#1 none; 2 Very Low only; 3 Very high only; 4 No response; 5 Low only; 6 High only;
  #ebin -2 is for the elements smaller than -2; ebin 0 for around 0; ebin 2 is for larger than 2
  one=matrix(0,3,length(relbins));
  coefE=1/(e*e);one[1,]=coefE; #
  #one[2,]=(1-coefE)*relbins;one[3,]=(1-coefE)*(1-relbins);#WARNING CHANGED TO MATCH PAPER DESCRIPTION
  #one[2,]=(1-coefE)*(1- 0.5*(e %in% c(1,7)));one[3,]=(1-coefE)*0.5*(e%in% c(1,7));#Possible to have A=1 or 2
  one[2,]=(1-coefE);one[3,]=0;#Only Possible to have A=1 #WARNING changed to fit cancer data better
  return(one)
}

factorQuantC= function(factorQuant,f,eharm, e,rel,diff,ecross=1){#e and rel should be rounded to index (multiply, integer divide)
  return(factorQuant[f,eharm,e,rel,diff,ecross]);
}

#factor Quantity/Quality to gene Generate Matrix
factorGeneGen= function(){
  geneFArray=array(0,dim=c(3,3,2,3));
  geneFArray[,,1,1]=1;
  geneFArray[,,1,c(2,3)]=0;
  for (i in 1:3)for(j in 1:3){
    x=min(i+j-2,2);
    geneFArray[i,j,2,]=as.numeric(x==c(0,1,2));
  }
  geneFArray[2,2,2,]=c(0,0.5,0.5);#if half from quality, half from quantity then the damage could be same allele or cumulative
  return(geneFArray)
}

#factor Quantity/Quality to gene
factorGene= function(Q,E,H,G,geneFArray){
  return(geneFArray[Q,E,H,G]);
}

factorPCommonGen=function(nbgenes,mx,possibleComplexity=3,pcase0=0.95,decay=c(0.05,0.1,0.2,0.4)){
  # Depends on the number of active genes and P (phenotype)
  size=min(mx,nbgenes);
  factorP=array(0,dim=c(2,nbgenes+1,possibleComplexity));
  factorP[1,1,]=1-pcase0;
  for (i in 2:size)factorP[1,i,]=decay[1:possibleComplexity]^(i-1);
  factorP[2,,]=1-factorP[1,,];
  return(factorP)
}

factorRegGen=function(nbgenes,mx,meancgenes=4,pnone=0.0,sd=1){
  size=min(mx,nbgenes)
  factorReg=rep(0,nbgenes+1);
  #factorReg[1:(size+1)]=dbinom(0:size,size,meancgenes/size); # also should think about using exp for regularization
  factorReg[1:(size+1)]=dnorm(0:size,meancgenes,meancgenes/sd);
  factorReg[1]=pnone;
  return(factorReg)
}

factorLayRegGen=function(nodesize,meancgenes,layersize,decay=0.5,pnone=0.8){# work only when meancgenes>layersize
  factorLayReg=matrix(0,2,2*nodesize+2);
  factorLayReg[1,1]=1;
  factorLayReg[2,]=1-factorLayReg[1,];
  if (nodesize>1) {
    coef=dbinom(1:(2*nodesize+1),2*nodesize+1,meancgenes/layersize);
  } else{ coef=dbinom(1:layersize,layersize,layersize/meancgenes);}
  #ind=0:(dim(factorLayReg)[2]-2);coef=(1-pnone)^ind;
  factorLayReg[2,-1]=coef;
  return(factorLayReg)
}

factorNetGen=function(netparams,maxConnections,nbgenes,meancgenes){
  factorNet=array(0.5,dim=c(2,2,maxConnections));
  netparam=netparams[1];netmax=netparams[2];netmaxNoise=netparams[3];
  epsi=meancgenes/nbgenes;
  denom=(1:meancgenes);if (netparams[4]) denom=denom/(1+netparams[4]);
  temp=1-1/(1+exp(log(netmax/epsi)/denom));
  maxSignal=c(temp,rep(temp[length(temp)],maxConnections-length(temp)));
  maxNoise=0.5+0.5*log(netmaxNoise/epsi)/(2*epsi*(1:maxConnections));#For less than 600 degree, this is not the bottleneck compared to signal. So for now I won' consider sum
  degparam=mapply(min,netparam,maxSignal,maxNoise)

  factorNet[2,2,]=degparam;factorNet[2,1,]=1-degparam;
  #factorNet[1,]=0.5;#Asymetric #by default 0.5

  return(factorNet)
}

factorLayGen=function(nodesize){
  factorLay=matrix(0,2,2*nodesize+1);
  factorLay[1,1]=1;factorLay[1,2]=0.0;
  #factorLay[2,1]=0.1;
  factorLay[2,]=1-factorLay[1,];
  #factorLay[2,2]=1;factorLay[2,3]=1;
  #if (nodesize>2){factorLay[2,4]=1;factorLay[2,5]=1;factorLay[2,6]=1;factorLay[2,7]=1}
  return(factorLay);
}


factorPGen=function(){ #Removabe function
  # TODO think about changing 0.5 to nbPatients/nbControls
  #Think about prior 0.1/0.9 that a control does actually have the disease genotype
  factorP=matrix(0,2,2);
  factorP[1,]=log(0.5);
  factorP[2,1]=log(0.001);
  factorP[2,2]=log(0.999);
  return(factorP)
}

