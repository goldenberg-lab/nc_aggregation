
initialize= function(nbgenes, nbpatients, nbsnps=c(),ncol=2,init=0){
  x=list(); length(x) <- nbgenes; 
  for (i in 1:nbgenes){
    x[[i]]<- list(); length(x[[i]])<- nbpatients; 
    if (length(nbsnps)>0) x[[i]]<- lapply(x[[i]], function(l,init,ncol,nbs){return(matrix(init,ncol,nbs))} ,init=init,ncol=ncol,nbs=nbsnps[i])
  }
  return(x)
}

initializeMem= function(nbgenes, nbpatients, ncol=2,init=0){
  x=array(init,dim=c(nbgenes,nbpatients,ncol));
  return(x)
}

transformmuq2= function(muq2m,muq2){
  for (i in 1:nbgenes)for(j in 1:nbpatients){
    muq2[[i]][[j]]=muq2m[i,j,];
  }
  return(muq2)
}

eToBins=function(e,ebins){#suppose ebins symmetric with 0 in the middle
  mid=1+length(ebins)%/%2;
  if (abs(e)<ebins[mid+1])return (mid);
  for (i in 1:(mid-1))if (e<=ebins[i])return (i);
  for (i in length(ebins):(mid+1))if (e>=ebins[i])return (i);
}

damping=function(x,y,alpha){return((1-alpha)*x+alpha*y);}

listbydim=function(m,d){
  return( lapply(apply(m,d, list),function(x)x[[1]]) )
}

listbydim13=function(m){
  return(lapply(1:(dim(m)[1]),function(x)m[x,,])  )
}

listbydim12=function(m){
  return(lapply(1:(dim(m)[1]),function(x)m[x,])  )
}

warp_reg= function(m,v){
  x=order(v,decreasing=TRUE);
  ini=m[1];mwarped=m/2;
  mwarped[x[1:20]]=ini*20;mwarped[x[21:100]]=ini*10;mwarped[x[101:200]]=ini*7;mwarped[x[201:350]]=ini*3.5;mwarped[x[351:500]]=ini*1.7;mwarped[x[501:1000]]=ini;
  mwarped=mwarped*sum(m)/sum(mwarped);
  return(mwarped);
}

estimate_causes= function(margg,margq,marge,muy,muy2,het,hom){
  top=10;snpthreshold=0.1;
  causes=list();length(causes) <- dim(margg)[2];
  xq=exp(margg[,,2]+marginal_aggr(margq));
  xe=exp(margg[,,2]+marginal_aggr(marge));
  for (i in 1:(dim(margg)[2])){
    mat=data.frame(stringsAsFactors=FALSE,check.names=FALSE,gene = numeric(0), prob = numeric(0), type = character(0), id=numeric(0));
    x=sort(c(xq[,i],xe[,i]),decreasing=TRUE);
    y=order(c(xq[,i],xe[,i]),decreasing=TRUE);
    for (j in 1:top){
      if (y[j]> dim(margg)[1]){mat=rbind(mat,data.frame(gene=y[j]%%(dim(margg)[1]),prob=x[j],type="Expr",id=0 ) );
      }else {
        mat=rbind(mat,data.frame(gene=y[j],prob=x[j], type="Qual",id=0 ))
        u=1-exp(logexpnormalized_2d(muy[[y[j]]][[i]]+muy2[[y[j]]][[i]]))[,1]
        potentialsnps=which(u> (snpthreshold*max(u,na.rm=TRUE)));
        if(length(potentialsnps)){
          toAdd=data.frame(y[j],x[j]*u[potentialsnps], "Qual",c(het[[y[j]]][[i]],hom[[y[j]]][[i]])[potentialsnps] );names(toAdd)<- c("gene","prob","type","id");
          mat=rbind(mat,toAdd);
        }
      }
      causes[[i]]<- mat;
    }
  }
  return(causes)
}

evaluate= function(gold,g){
  return(sum(abs(gold-g)))
}

evaluate_parallel= function(gold,g){
  return(sum(mapply(evaluate,gold,g)))
}

sum2lists= function(l1,l2){
  return(mapply(function(x,y){return(x+y);},l1,l2 ,SIMPLIFY=FALSE))
}

substract2lists= function(l1,l2){
  return(mapply(function(x,y){return(x-y);},l1,l2 ,SIMPLIFY=FALSE))
}

is_diverging=function(x){#require x of length 4
return ( ((x[2]-x[1])*(x[3]-x[1]) <0) & ((x[2]-x[1])*(x[4]-x[1]) >0))
}

verif=function(x,y)any(abs(x-y)>0.00000000001)
verifexp=function(x,y)any(abs(exp(x)-exp(y))>0.00000000001)

shuffler= function(nbgenes){
  x=sample(1:nbgenes);
  y=order(x);
  return(list(pos=x,posinverse=y))
}

build_layer= function(prevlayersize,layersize){# return start,end for the nodes indexed affected to each node in the new layer
  mat=matrix(1,layersize,2);
  rest=prevlayersize%%layersize;
  avg=prevlayersize%/%layersize;
  k=1;
  if (rest){for (i in 1:rest){mat[i,1]=k; mat[i,2]=k+avg; k=k+ avg+1;  }}
  for (i in (rest+1): layersize){mat[i,1]=k; mat[i,2]=k+avg-1; k=k+ avg;}
  return(mat);
}

colSumsOrzero=function(mat){
  if (length(mat)==1)return(mat);
  if (length(mat))return(colSums(mat))
  return(c(0,0));
}

multiplierMat= function(nb){
  return(matrix(1,nb,nb)-diag(1,nb))
}

compute_likelihood1= function(mugppat,pheno,margC,factorP,mx){
  nonzer=which(exp(mugppat[,1])!=1);
  n=length(nonzer);
  if(n==0){z=1;}else if(n==1){z=exp(mugppat[nonzer,]);}else {z=(compute_mu1_cardinality(exp(mugppat),n-1,n,nonzer,1))[[1]];}
  res=(z%*%(pheno*factorP[2,1:length(z),]+ (1-pheno)*factorP[1,1:length(z),]))%*%t(margC);
  res2=(z%*%(pheno*factorP[1,1:length(z),]+ (1-pheno)*factorP[2,1:length(z),]))%*%t(margC);
  return(log(res/(res+res2)))
}

predict_pheno= function(mugppat,margC,factorP,mx,thresh=0){#threshold is to remove the cumulative small effects of inactive genes. TODO
  nonzer=which(exp(mugppat[,1])!=1);
  n=length(nonzer);
  if(n==0){z=1;}else if(n==1){z=exp(mugppat[nonzer,]);}else {z=(compute_mu1_cardinality(exp(mugppat),n-1,n,nonzer,1))[[1]];}
  res1=(z%*%factorP[1,1:length(z),])%*%t(margC);
  res2=(z%*%factorP[2,1:length(z),])%*%t(margC);
  return(res2/(res2+res1))
}

marginal_x=function(mux2,mux){
if(length(mux)==0)return(NULL);
all=t(apply(mux2,c(2,3),sum))+mux;
all2=exp(all-rbind(colMeans(all),colMeans(all)));
return(t(log(all2/rbind(colSums(all2),colSums(all2)))))
}

nb_mutation_gene=function(het,hom,pheno)length(unlist(het[which(pheno==1)]))+length(unlist(hom[which(pheno==1)]));

