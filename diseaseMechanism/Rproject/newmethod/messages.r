
compute_etah= function(muh2,muh){
  muh2[muh2==-Inf] <- -10000;
  all=apply(muh2,c(1,3),sum)+muh;
  all=all-apply(all,1,max);
  x=-muh2;
  for (i in 1:(dim(muh2)[2]))x[,i,]=normalize_etah(x[,i,]+all)
  return(x)
}

compute_etaf= function(i,muf2,all){
  return(normalize_etaf(all-muf2[,i,]));
}

compute_etaC= function(i,muC2,all){
  return(normalize_etaC(all-muC2[[i]]));
}

compute_muC2= function(mugppat,pheno,factorP,mx){
  nonzer=which(exp(mugppat[,1])!=1);
  n=length(nonzer);lp=n-1;b=1; 
  if(n==0){z=1;}else if(n==1){z=exp(mugppat[nonzer,]);}else {z=(compute_mu1_cardinality(exp(mugppat),lp,n,nonzer,b))[[1]];}
  res=z%*%((1-pheno)*factorP[1,1:length(z),]+ pheno* factorP[2,1:length(z),]);
  return(log(res/sum(res)))
}

compute_mue_zero= function(nbgenes,nbpatients,factorQuant){
  x=array(-Inf,dim=c(nbgenes,nbpatients,3));
  x[,,1]=0
  return(x);
}

compute_mue= function(exetaf,ed,factorQuant){#per gene per patient
  x=exetaf%*%factorQuant[,,ed,5,1,1];
  return(log(x/sum(x)));
}

compute_muf2= function(exmue2,ed,factorQuant){#per gene per patient
  x=factorQuant[,,ed,5,1,1]%*%exmue2;
  return(log(x/sum(x)));
}

compute_mug2=function(mu,exetaC,pheno,factorP){
  potential=( (pheno*factorP[2,,]+ (1-pheno)*factorP[1,,])%*%exetaC)
  return(comupute_mu_cardinality(exp(mu),potential,init=c(0.5,0.3,0.2)));
}

compute_mureg=function(etareg,factorReg){
  y=rowMeans(etareg);
  expreg=exp(etareg-cbind(y,y));
  expetareg=expreg/rowSums(expreg);
  mureg=comupute_mu_cardinality(expetareg,factorReg);
  return(mureg);
}

sumallminusone=function(i,mat){if(dim(mat)[1]==2)return(mat[-i,]); return(colSums(mat[-i,]))}
whichis=function(vec,i){return(which(vec==i))}

compute_etanet= function(munet,incoming,net){
  if (length(net)==0)return(NULL);
  if (length(net)==1){ex=exp(incoming);return(log(ex/sum(ex))); }
  allminusone=mapply(sumallminusone,1:length(net),MoreArgs=list(mat=munet));
  etanetraw=t(allminusone+incoming);
  return(logexpnormalized_2d(etanetraw));
}

compute_munet_one= function(etanet,pos,factorNet){
  if (length(etanet)==2)return(log(exp(etanet)%*%factorNet));
  return(log(exp(etanet[pos,])%*%factorNet ));
}

compute_munet= function(igene,net,etanet,factorNet){
  src=net[[igene]];
  if (length(src)==0)return(NULL);
  pos=mapply(whichis,net[src],MoreArgs=list(i=igene));
  munetraw=t(mapply(compute_munet_one,etanet[src],pos,MoreArgs=list(factorNet=factorNet[,,length(src)])));
  return(logexpnormalized_2d(munetraw));
}

compute_muq2= function(mug2, mue,etah,geneFArray,propagate=TRUE){
  res=array(0,dim=c(dim(mug2)[1:2],3));
  for (i in 1:3){for (j in 1:3){for (k in 1:3){
    if (propagate){res[,,i]=res[,,i]+exp(mug2[,,k]+mue[,,j])*geneFArray[i,j,2,k];
    }else{res[,,i]=res[,,i]+ exp(mug2[,,k]+mue[,,j]+etah[,,1])*geneFArray[i,j,1,k]+ exp(mug2[,,k]+mue[,,j]+etah[,,2])*geneFArray[i,j,2,k];}
  }}}
  return(lognormalized(res))
}

compute_mue2= function(mug2, muq, etah,geneFArray,propagate=TRUE){
  res=array(0,dim=c(dim(mug2)[1:2],3));
  for (i in 1:3){for (j in 1:3){for (k in 1:3){
    if (propagate){res[,,j]=res[,,j]+exp(mug2[,,k]+muq[,,i])*geneFArray[i,j,2,k];
    }else{res[,,j]=res[,,j]+ exp(mug2[,,k]+muq[,,i]+etah[,,1])*geneFArray[i,j,1,k]+ exp(mug2[,,k]+muq[,,i]+etah[,,2])*geneFArray[i,j,2,k];}
  }}}
  return(lognormalized(res))
}

#compute_mug= function(muq, mue, etah,geneFArray){
#  res=array(0,dim=c(nbgenes,nbpatients,3));
#  for (i in 1:3){for (j in 1:3){for (k in 1:3){
#    res[,,k]=res[,,k]+ exp(muq[,,i]+mue[,,j]+etah[,,1])*geneFArray[i,j,1,k]+ exp(muq[,,i]+mue[,,j]+etah[,,2])*geneFArray[i,j,2,k];
#  }}}
#  return(lognormalized(res))
#}
#compute_muh2= function(mug2, muq, mue,geneFArray){
#  res=array(0,dim=c(nbgenes,nbpatients,2));
#  for (i in 1:3){for (j in 1:3){for (k in 1:3){
#    res[,,1]=res[,,1]+ exp(mug2[,,k]+muq[,,i]+mue[,,j])*geneFArray[i,j,1,k];
#    res[,,2]=res[,,2]+ exp(mug2[,,k]+muq[,,i]+mue[,,j])*geneFArray[i,j,2,k];
#  }}}
#  return(lognormalized(res))
#}

#Faster versions of compute_muh2 and compute_mug taking advantage of the zeros in the factor
compute_muh2= function(mug2, muq, mue,geneFArray){
  res=array(0,dim=c(dim(mug2)[1:2],2));
  res[,,2]=res[,,2]+ exp(mug2[,,1]+muq[,,1]+mue[,,1])*geneFArray[1,1,2,1];
  for (i in 1:3){for (j in 1:3){
    res[,,1]=res[,,1]+ exp(mug2[,,1]+muq[,,i]+mue[,,j])*geneFArray[i,j,1,1];
    res[,,2]=res[,,2]+ exp(mug2[,,2]+muq[,,i]+mue[,,j])*geneFArray[i,j,2,2];
    res[,,2]=res[,,2]+ exp(mug2[,,3]+muq[,,i]+mue[,,j])*geneFArray[i,j,2,3];
  }}
  return(lognormalized(res))
}

compute_mug= function(muq, mue, etah,geneFArray){
  res=array(0,dim=c(dim(muq)[1:2],3));
  res[,,1]=res[,,1]+ exp(muq[,,1]+mue[,,1]+etah[,,2])*geneFArray[1,1,2,1];
  for (i in 1:3){for (j in 1:3){
    res[,,1]=res[,,1]+ exp(muq[,,i]+mue[,,j]+etah[,,1])*geneFArray[i,j,1,1];
    res[,,2]=res[,,2]+ exp(muq[,,i]+mue[,,j]+etah[,,2])*geneFArray[i,j,2,2];
    res[,,3]=res[,,3]+ exp(muq[,,i]+mue[,,j]+etah[,,2])*geneFArray[i,j,2,3];
  }}
  return(lognormalized(res))
}

compute_muy2= function(muq2, muy, het, hom, factorQual){#for now I am only computing the het messages for het and hom for hom
  #w is the probability of having a damaging hom hit, z is the distribution over sum of damaging het
  #muy2=matrix(0,3, dim(muy)[2]);
  if (length(het)+length(hom) ==0)return(NULL)
  n=length(het);l1=list(1);
  if (length(hom)==0){w=0;muy2hom=c();} else { w=1-exp(sum(muy[n+(1:length(hom)),1]));}
  if(n==0){
    muy2het=c();
  }else if(n==1){
    muy2het=(w*factorQual[2, 1:2, ]+ (1-w)*factorQual[1, 1:2, ]) %*% exp(muq2); ## could be faster when w=0 (often the case)    ##removedc(, factorQual[2, 1, ]%*% exp(muq2)
    l1[[1]]=exp(muy[1,1:2]);
  }else {
    l1=compute_mu1_cardinality_all(exp(muy[1:n,1:2]),n-1,n);
    potential=(w*factorQual[2, 1:(n+1), ]+ (1-w)*factorQual[1, 1:(n+1), ]) %*% exp(muq2);
    l2=compute_mu2_cardinality_all(l1,potential[1:length(l1[[1]])],n-1,n);
    muy2het=simplify2array(l2[n:(2*n-1)]);##removed rbind(,factorQual[2, 1, ]%*% exp(muq2))
  }
  if (length(hom)){
    #default=compute_convolution2(c(potential[1:length(l1[[1]])],0),l1[[1]],2);
    if (length(hom)==1){
      muy2hom= c(l1[[1]] %*% (factorQual[1, 1:length(l1[[1]]), ]%*% exp(muq2)),0, l1[[1]] %*% (factorQual[2, 1:length(l1[[1]]), ]%*% exp(muq2)))
    }else{
      w1=rep(0,length(hom));
      for (i in 1:length(hom))w1[i]=sum(muy[n+((1:length(hom))[-i]),1])## was w1=1-exp(log(1-w)- muy[1,(n+1):(n+length(hom))]);
      w2=1-exp(w1);
      muy2hom=matrix(0,3,length(hom))
      for (i in 1:length(hom)){
        muy2hom[1,i]= l1[[1]] %*% (((1-w2[i])*factorQual[1, 1:length(l1[[1]]), ] + w2[i]*factorQual[2, 1:length(l1[[1]]), ])%*% exp(muq2));
        muy2hom[3,i]= l1[[1]] %*% (factorQual[2, 1:length(l1[[1]]), ]%*% exp(muq2));
      }
    }
  }

  if(length(muy2het))muy2het=rbind(muy2het,0);
  res=cbind(muy2het,muy2hom);
  return(log(t(res)/colSums(res)))
}

compute_mux2= function(muy2,het,hom,mux2){
  x=mux2;
  if (length(mux2)){
    mapply(function(i,het,hom,muy2){ if (length(het))x[i,het,1:2]<<- muy2[1:length(het),1:2]; if(length(hom))x[i,hom,1:2]<<- muy2[length(het)+(1:length(hom)),c(1,3)]; },1:(dim(mux2)[1]),het,hom,muy2);
    #for (i in 1:(dim(mux2)[1])){
    #  if (length(het[[i]])) x[i,het[[i]],1:2]= muy2[[i]][1:length(het[[i]]),1:2];
    #  if (length(hom[[i]])) x[i,hom[[i]],1:2]= muy2[[i]][length(het[[i]])+(1:length(hom[[i]])),c(1,3)];
    #}
  }
  return(x);
}

compute_etax=function(i,het,hom,mux2,all,mux){
  etax=matrix(0,ncol(mux),2);
  indc=c(het,hom);
  if (length(indc)){
    x= exp(all[indc,1]-mux2[i,indc,1]+mux[1,indc]);
    y= exp(all[indc,2]-mux2[i,indc,2]+mux[2,indc]);
    etax_s=x+y+0.0000000000000000001;
    etax[indc,1]= log(x/etax_s);
    etax[indc,2]= log(y/etax_s);
  }
  return(etax)
}

compute_muy= function(i,het,hom,muy,etax){
  x=muy;
  x[,1]= etax[i,c(het,hom),1];
  if (length(het))x[ 1:length(het),2]= etax[i,het,2]
  if (length(hom))x[ length(het)+(1:length(hom)),3]= etax[i,hom,2]
  return(x)
}

compute_dhet=function(muy,het){
  if(length(het)==0)return(1);
  n=length(het);
  if(n==0)return(1);
  if(n==1)return(exp(muy[1,1:2]));
  dhet=(compute_mu1_cardinality_all(exp(muy[1:n,1:2]),n-1,n))[[1]]; 
  return(dhet);
}

compute_muq= function(muy,het,hom,factorQual){
  if ((length(het)+length(hom))==0)return(c(0,-Inf,-Inf))
  if (length(hom)>0){
    phom=1-exp(sum(muy[length(het)+(1:length(hom)),1]));
  } else phom=0;
  if (phom==1)return(c(-Inf,-Inf,0));
  dhet=compute_dhet(muy,het);
  phet=t(dhet)%*%factorQual[1,1:(length(dhet)),];
  return(log((1-phom)*phet+c(0,0,phom)))
}


test_new=function(){
library(Matrix)
het=list(c(1,3),c(2,3),c(2),c())
hom=list(c(),c(),c(),c());
mux2=list(Matrix(0,4,3,sparse=TRUE),Matrix(0,4,3,sparse=TRUE))
mux=log(rbind(c(0.5,0.7,0.3), c(0.5,0.3,0.7)))
muy=list(matrix(0,2,3),matrix(0,2,3),matrix(0,1,3),matrix(0,0,3))
etax=compute_etax(mux2,mux)
muy=mapply(compute_muy,1:4,het,hom,muy,MoreArgs=list(etax=etax),SIMPLIFY=FALSE)
dhet=mapply(compute_dhet,muy,het,SIMPLIFY=FALSE)
muq=mapply(compute_muq,muy,dhet,het,hom,MoreArgs=list(factorQual=factorQual),SIMPLIFY=FALSE)
#muq2=muq2[[1]][c(1,101,2,102)]
muy2=mapply(compute_muy2,muq2, muy, het, hom,MoreArgs=list(factorQual=factorQual),SIMPLIFY=FALSE)
mux2new=compute_mux2(muy2,het,hom,mux2)
}


