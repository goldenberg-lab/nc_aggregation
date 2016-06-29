
#Belief propagation grid search function definition
grid_search= function(codeDir,nbgenes,nbpatients,nbsnps,harm,harmgene,meancgenes,complexityDistr,pheno,hom,het,net=NULL,e=NULL, cores=1,ratioSignal=0.9,decay=c(0.05,0.1,0.2,0.4),alpha=0.5,netparams=c(0.9,0.01,0.01,0),removeExpressionOnly=FALSE,propagate=TRUE,truth=NULL,internalLoop=TRUE){
  stepmax=1;multiplier=2;
  ratioSignalt=ratioSignal;diverging=TRUE;
  step=0;
  while (diverging && step<stepmax){
    ptm <- proc.time()
    x <- sumproduct(codeDir,nbgenes,nbpatients,nbsnps,harm,harmgene,meancgenes,complexityDistr,pheno,hom,het,net=net,e=e, cores=cores,ratioSignal=ratioSignalt,decay=decay,alpha=alpha,netparams=netparams,removeExpressionOnly=removeExpressionOnly,propagate=propagate,truth=truth,internalLoop=internalLoop);
    print(proc.time()-ptm)
    diverging=(x$status==0) ; step=step+1;
    if (diverging){
      print(genenames[which(x$h>0.2)]);print(x$h[which(x$h>0.2)]);print(x$margC);
      if (step<stepmax){ratioSignalt=multiplier*ratioSignalt;print(paste("Belief propagation did not converge. Now trying with ratioSignal=",ratioSignalt));}
    }
  }
  return(x);
}


sumproduct= function(codeDir,nbgenes,nbpatients,nbsnps,harm,harmgene,meancgenes,complexityDistr,pheno,hom,het,mes,net=NULL,e=NULL,cores=1,ratioSignal=0.9,decay=c(0.05,0.1,0.2,0.4),alpha=0.5,netparams=c(0.9,0.01,0.01,0), removeExpressionOnly=FALSE,propagate=TRUE,truth=NULL,internalLoop=TRUE){
  #pcase0 :uncertainty about phenotype, alpha is the damping coefficient, netparams are the default param, the max proba contrib by real signal and the max proba contrib by noise
  #sum product algorithm
  #We compute messages from leaves to root and top to bottom then messages from root to leaves and bottom to top(indicators).

  #mu messages from factors to nodes, eta from nodes to factor.
  #mux prior to X; mux2 other way
  #etax X to factor (bw Y); # should take all mux2 and one mux to compute it
  #muy factor XY to Y; muy2 factor QY to Y; should use all etay and the etaq2 to compute it
  #etay =mes$muy and etay2=mes$muy2
  #muq YQ factor to Q use all etay and big factor to compute
  #muq2 from GQ factor to Q; is computed using etag2, etah and etae
  #Same for mue2 uses etag2, etah and etaq
  #etaq=mes$muq and etaq2=mes$muq2 and etae=mes$mue and etae2=mes$mue2
  #mug from EQGH to G , uses etae, etaq  and etah
  #muh prior on genes; muh2 the other way uses etae, etaq, etag2
  #etah uses all other muh2 and one muh
  #etag=mes$mug ,etag2=mes$mug2

  verbo=TRUE;
  toomuch=20*meancgenes;#if sum of genes marginals goes over this stop the algorithm (regularization will fail)
  maxiter=30;convergencedelta=0.005*nbpatients;#was =2   
  convergencedeltaInternal=0.1*convergencedelta;
  iteralphareduction=10; #iteration at which alpha is halved to improve convergence
  repmax=1;if (internalLoop)repmax=19;reps2=1;#Internal loops sizes
  possibleComplexity=length(complexityDistr);#n+1 number of mutations a control can have (or number of mutations needed to cause disease)
  lf=6;nullpriorexpr=0.;iterexpr=4;  #after iterexpr iterations I start considering expression levels
  ebins=c(-3,-2.5,-2,0,2,2.5,3);relbins=c((1:10)/10);#releff=exp(-abs(log(relbins)));
  pcase0=max((1-ratioSignal)/(2-ratioSignal),0.1/1.1)#pacase0 can't be more extreme than ratioSignal of 90%
  sd=2;

  #library(parallel)#for sumproduct
  source(paste(codeDir,"newmethod/factors.r",sep=""))
  source(paste(codeDir,"newmethod/normalization_functions.r",sep=""))
  source(paste(codeDir,"newmethod/extra_functions.r",sep=""))
  source(paste(codeDir,"newmethod/messages.r",sep=""))
  source(paste(codeDir,"newmethod/messages_parallel.r",sep=""))
  source(paste(codeDir,"newmethod/messages_cardinality.r",sep=""))

  #going from top to bottom and from leaves (snps) to root
  delta=100;iter=1;iter2=1;
  g=marginal(mes$mug,mes$mug2)

  while((delta>convergencedelta && iter<maxiter) | (iter<iterexpr+2) ){
    ptm1 <- proc.time();#Rprof(filename = "Rprof.out")summaryRprof(filename = "Rprof.out")
    gold=g;
    if((iter%%iteralphareduction)==0){alpha=alpha*0.5;print(paste("Alpha reduced to ",alpha));}
    etax=mapply(compute_etax_parallel,mes$mux2,mes$mux,het,hom,SIMPLIFY=FALSE);
    mes$muy=mapply(compute_muy_parallel,het,hom,mes$muy,etax,SIMPLIFY=FALSE);
    muql=mapply(compute_muq_parallel,mes$muy,het,hom,MoreArgs=list(factorQual=mes$factorQual),SIMPLIFY=FALSE);
    muqm=matrix(unlist(muql), nrow=3);
    for( i in 1:3)mes$muq[,,i]=t(matrix(muqm[i,],ncol=nbgenes));
    #This part is not necessary if mue2 is unchanged
    if (length(e) && iter>iterexpr){
      etaf=compute_etaf_parallel(mes$muf2,mes$muf);
      mes$mue=aperm(mapply(compute_mue_parallel,listbydim13(exp(etaf)),listbydim12(ed),MoreArgs=list(factorQuant),SIMPLIFY="array"),c(3,2,1));
    } 
    if (verbo)print(paste("Higher Tier done",(proc.time()-ptm1)[3]));

    reps=c(0);repn=0;h=0;
    passed=FALSE;
    alpha1=alpha;
    while (!passed && repn<repmax){
      ptm <- proc.time()
      newmuh2=compute_muh2(mes$mug2,mes$muq,mes$mue,mes$geneFArray);
      alphat=1; st=2*meancgenes*(1+4/sd)+1; sti=0;
      while(st> 2*meancgenes*(1+4/sd)){
        if (sti)print(paste("Oscillation too strong ",sti,"sum=",st, ". Damping more",alphat))
        temp=damping(mes$muh2,newmuh2,alphat);
        incomingup=mes$muh+apply(temp,c(1,3),sum);
        incomingup=incomingup-apply(incomingup,1,max);
        alphat=alphat*0.75;st=sum(exp(incomingup[,2]));sti=sti+1;
      }
      if (alphat<0.75)alpha1=alpha1*0.75;
      mes$muh2=temp;
      #Regularization messages
      etareg=incomingup+mes$munetall;
      mes$mureg=compute_mureg(etareg,mes$factorReg);#Warning why 2 #Error dim is 2,2  
      incoming=split(mes$mureg+incomingup, 1:nbgenes);
      #Gene network messages
      if(length(net) ){
        for (r2 in 1:reps2){
  	  etanet=mapply(compute_etanet,mes$munet,incoming,net);
	  mes$munet=mapply(damping,mes$munet,mapply(compute_munet,1:nbgenes,MoreArgs=list(net=net,etanet=etanet,factorNet=factorNet)),MoreArgs=list(alpha=alpha),SIMPLIFY=FALSE);
        }
        mes$munetall=logexpnormalized_2d(t(mapply(colSumsOrzero,mes$munet)));
        if(any(is.na(mes$munetall)))debug(mes$munetall,mes$mureg,incomingup,mes$muh2[,,2],mes$mug2[,,2],mes$mug[,,2]); 
        mes$munetall=t(apply(mes$munetall,1,function(x){if( (x[2]-x[1])> (maxnetcontrib[2]-maxnetcontrib[1]))return(maxnetcontrib);return(x);}))
      }
      #Regularization messages (second time: after network dispersion)
      etareg=incomingup+mes$munetall;
      mes$mureg=compute_mureg(etareg,mes$factorReg);#Warning why 2 #Error dim is 2,2 
 
      etah=compute_etah(mes$muh2,mes$muh+mes$mureg+mes$munetall);
      mes$mug=compute_mug(mes$muq,mes$mue,etah,mes$geneFArray);  
      mes$mugppat=listbydim(mes$mug,2);     
      mes$mug2=damping(mes$mug2,aperm(mapply(compute_mug2,mes$mugppat,lapply(mes$etaC,function(x)t(exp(x))),pheno,MoreArgs=list(factorP=mes$factorP), SIMPLIFY="array"),c(1,3,2)),alpha1);
      hex=exp(incomingup+mes$mureg+mes$munetall);
      h=hex[,2]/rowSums(hex);if (verbo)print(paste("sum of genes=",sum(h),". time=",(proc.time()-ptm)[3]));
      if(sum(h)>toomuch){print("Unlikely to converge. Try changing ratioSignal or alpha");return(list(status=FALSE));}
      mes$hsave[iter2,]=h;iter2=iter2+1;
      if(length(truth)){print(which(order(h,decreasing=TRUE) %in% truth))}else {print(genenames[order(h,decreasing=TRUE)[1:5]]);}
      passed= (abs(reps[length(reps)]-sum(h)) <convergencedeltaInternal);#passed=any(abs(reps-sum(h)) <convergencedeltaInternal) #just check change since last step now
      if (repn>3)if (repn%%4==1 & is_diverging(c(reps[(repn-2):repn],sum(h))) ){alpha1=alpha1*0.75;print(paste("alpha1=",alpha1));} else {if (repn==13)alpha1=alpha1*0.75;}     
      reps=c(reps,sum(h));repn=repn+1;
      #print(h[7759]);#save(list=ls(),file=paste("/hpf/largeprojects/agoldenb/aziz/data",iter2,".RData",sep=""))
    }

    mes$muC2=mapply(damping,mes$muC2,mapply(compute_muC2,mes$mugppat,pheno,MoreArgs=list(factorP=mes$factorP,mx=mes$mx),SIMPLIFY=FALSE),MoreArgs=list(alpha=alpha),SIMPLIFY=FALSE);  
    mes$etaC=compute_etaC_parallel(mes$muC2,mes$muC);
    mes$margC=mes$muC+Reduce('+',mes$muC2);mes$margC=exp(mes$margC-max(mes$margC));margC=mes$margC/sum(mes$margC);
    if (verbo)print("Lower Tier done");

    if (length(e) && iter>=iterexpr){
      mes$mue2=compute_mue2(mes$mug2,mes$muq,etah,mes$geneFArray,propagate);
      #This part is not necessary if mes$mue2 is unchanged
      mes$muf2=aperm(mapply(compute_muf2_parallel,listbydim13(exp(mes$mue2)),listbydim12(ed),MoreArgs=list(factorQuant),SIMPLIFY="array"),c(3,2,1));
    }
    mes$muq2m=compute_muq2(mes$mug2,mes$mue,etah,mes$geneFArray,propagate);
    mes$muq2=lapply(listbydim13(mes$muq2m),listbydim12);
    mes$muy2=mapply(compute_muy2_parallel,mes$muq2,mes$muy,het,hom,MoreArgs=list(factorQual=mes$factorQual),SIMPLIFY=FALSE);
    mes$mux2=mapply(compute_mux2,mes$muy2,het,hom,mes$mux2,SIMPLIFY=FALSE);

    hex=exp(etareg+mes$mureg);
    h=hex[,2]/rowSums(hex);
    g=marginal(mes$mug,mes$mug2);
    delta=sum(abs(exp(gold[,,2])-exp(g[,,2])));
    print(paste("Iteration: ",iter,". delta= ",delta, " time=", (proc.time()-ptm1)[3]));
    #write.table(h,paste(dir,"genes",(iter%%10),".txt",sep=""))
    iter=iter+1; 
  }

  margg=marginal(mes$mug,mes$mug2);
  margq=marginal(mes$muq,mes$muq2m);
  marge=marginal(mes$mue,mes$mue2);
  mes$margC=mes$muC+Reduce('+',mes$muC2);mes$margC=exp(mes$margC-max(mes$margC));mes$margC=mes$margC/sum(mes$margC);
  gppat=listbydim(g,2);
  likelihood=mapply(compute_likelihood1,mes$mugppat,pheno,MoreArgs=list(margC=mes$margC,factorP=mes$factorP,mx=mes$mx))
  mes$predict=mapply(predict_pheno,mes$mugppat,MoreArgs=list(margC=mes$margC,factorP=mes$factorP,mx=mes$mx,thresh=0.0))
  mes$causes=estimate_causes(margg,margq,marge,mes$muy,mes$muy2,het,hom)
  print(paste("finished with alpha=",alpha))
  if (delta>convergencedelta)h=colMeans(mes$hsave[(iter2-20):(iter2-1),])
  mes$status=(delta<=convergencedelta);
  mes$g=g[,,2];
  mes$h=h;
  mes$iter=iter-1;
  mes$likelihood=likelihood;

  return(mes)
  #return(list(status=(delta<=convergencedelta),g=g[,,2],h=h,iter=iter-1,causes=causes,hsave=hsave,margC=margC,munetall=munetall,pcase0=pcase0,mureg=mureg,likelihood=likelihood,predict=predict,mux2=mux2,mux=mux,mug=mug,muq=muq,muq2=muq2,mue=mue,mue2=mue2,muf=muf,muf2=muf2,factorP=factorP,factorQual=factorQual,geneFArray=geneFArray,mx=mx))
}

sumproduct_predict=function(x,het_v,hom_v,thresh){
    nbpatients_v=length(het_v[[1]])
    nbgenes=length(het_v)
    muy=list();length(muy) <- nbgenes;for (j in 1:nbgenes){muy[[j]]=list();length(muy[[j]]) <- nbpatients_v; for (k in 1:nbpatients_v)muy[[j]][[k]]=matrix(-Inf,length(het_v[[j]][[k]])+length(hom_v[[j]][[k]]),3);}
    mue=initializeMem(nbgenes,nbpatients_v,3,-Inf);mue[,,1]=0;
    muq=initializeMem(nbgenes,nbpatients_v,3); 

    etaxall=mapply(marginal_x,x$mux2,x$mux,SIMPLIFY=FALSE);
    etax=list();length(etax)=length(etaxall);
    for (i in 1:nbgenes){
      etax[[i]]=array(0,dim=c(nbpatients_v,dim(etaxall[[i]])));
      for (j in 1:nbpatients_v)etax[[i]][j,,]=etaxall[[i]];
    }
    muy=mapply(compute_muy_parallel,het_v,hom_v,muy,etax,SIMPLIFY=FALSE);
    muql=mapply(compute_muq_parallel,muy,het_v,hom_v,MoreArgs=list(factorQual=x$factorQual),SIMPLIFY=FALSE);
    muqm=matrix(unlist(muql), nrow=3); 
    for( i in 1:3)muq[,,i]=t(matrix(muqm[i,],ncol=nbgenes));
    margH=x$h; margH[which(x$h<thresh)]=0.000000000000001; #which(x$h<thresh)
    etah=array(0,dim=c(nbgenes,nbpatients_v,2));for(j in 1:nbpatients_v)etah[,j,]=log(cbind(1-margH,margH));
    mug=compute_mug(muq,mue,etah,x$geneFArray);  
    mugppat=listbydim(mug,2);
    predict=mapply(predict_pheno,mugppat,MoreArgs=list(margC=x$margC,factorP=x$factorP,mx=x$mx,thresh=thresh))#threshold is not used here
    return(predict);
}

init=function(codeDir,nbgenes,nbpatients,nbsnps,harmgene,meancgenes,complexityDistr,pheno,hom,het,ratioSignal=0.9,netparams=c(0.9,0.01,0.01,0),internalLoop=TRUE) {
  # Returns a named list with the default initializations for all messages
  verbo=TRUE;

  source(paste(codeDir,"newmethod/factors.r",sep=""))
  source(paste(codeDir,"newmethod/extra_functions.r",sep=""))
  if (length(e)){
  if(removeExpressionOnly){e[which(mapply(nb_mutation_gene,het,hom,MoreArgs=list(pheno=pheno))==0),]=0;}
    ed=apply(e,1:2,eToBins,ebins=ebins);
  }
  maxiter=30;convergencedelta=0.005*nbpatients;#was =2   
  repmax=1;if (internalLoop)repmax=19;reps2=1;#Internal loops sizes
  pcase0=max((1-ratioSignal)/(2-ratioSignal),0.1/1.1)#pacase0 can't be more extreme than ratioSignal of 90%
  maxsnps=max(nbsnps);
  ratmaxnetcontrib=netparams[2]*nbgenes/((1-netparams[2]) *meancgenes);#ratmaxnetcontrib should always be >1 : net maximal contribution more than prior. netparams[2] should be set accordingly
  maxnetcontrib=log(c(1/(ratmaxnetcontrib+1),ratmaxnetcontrib/(ratmaxnetcontrib+1)));
  lf=6;nullpriorexpr=0.;iterexpr=4;  #after iterexpr iterations I start considering expression levels
  possibleComplexity=length(complexityDistr);#n+1 number of mutations a control can have (or number of mutations needed to cause disease)
  ebins=c(-3,-2.5,-2,0,2,2.5,3);relbins=c((1:10)/10);#releff=exp(-abs(log(relbins)));

  #parameters nbgenes, nbpatients, nbsnp in a vector nbsnp[gene], harm is a list(genes) of tables of harmf. pred.

  #intialize all messages
  mux=list(); length(mux) <- nbgenes; for (j in 1:nbgenes)mux[[j]]=log(rbind(1-harm[[j]],harm[[j]]));
  initsnp=log(c(0.4,0.6));initsnp2=log(c(0.37,0.63));
  mux2=list(); length(mux2) <- nbgenes; for (j in 1:nbgenes) if(nbsnps[j]){
    mux2[[j]]=array(0,dim=c(nbpatients,nbsnps[j],2));
    for (k in 1:nbpatients)
      if(length(c(het[[j]][[k]],hom[[j]][[k]]))){
        if (pheno[k]){
          mux2[[j]][k,c(het[[j]][[k]],hom[[j]][[k]]),1]=initsnp[1]; mux2[[j]][k,c(het[[j]][[k]],hom[[j]][[k]]),2]=initsnp[2];
        } else{
          mux2[[j]][k,c(het[[j]][[k]],hom[[j]][[k]]),1]=initsnp2[2]; mux2[[j]][k,c(het[[j]][[k]],hom[[j]][[k]]),2]=initsnp2[1];
        }
    }    
  }
  muy=list();length(muy) <- nbgenes;for (j in 1:nbgenes){muy[[j]]=list();length(muy[[j]]) <- nbpatients; for (k in 1:nbpatients)muy[[j]][[k]]=matrix(-Inf,length(het[[j]][[k]])+length(hom[[j]][[k]]),3);}
  #initMat=lapply(1:maxsnps,function(x)matrix(0,3,x));#To Make muy2 computation faster (don't have to intialize a new matrix every time)
  muq=initializeMem(nbgenes,nbpatients,3); 

  muf=matrix(log((1-nullpriorexpr)/lf),nbgenes,lf);
  if (lf>6)
	  muf[,7]=log(nullpriorexpr);
  muf2=initializeMem(nbgenes,nbpatients,lf);

  mue=initializeMem(nbgenes,nbpatients,3,-Inf);
  mue[,,1]=0;

  muh=matrix(0,nbgenes,2);
  muh[,1]=log(1-harmgene);
  muh[,2]=log(harmgene);
  muh2=initializeMem(nbgenes,nbpatients,2);

  munet=list();length(munet) <- nbgenes;
  if(length(net)){for (j in 1:nbgenes)munet[[j]]=matrix(0,length(net[[j]]),2);}
  munetall=matrix(0,nbgenes,2);

  mug=initializeMem(nbgenes,nbpatients,3);

  muC=log(complexityDistr);
  muC2=listbydim12(matrix(0,nbpatients,possibleComplexity));

  etaC=list();
  length(etaC)=nbpatients; 
  for (i in 1:nbpatients)etaC[[i]]=t(muC);

  mue2=initializeMem(nbgenes,nbpatients,3);
  mug2=initializeMem(nbgenes,nbpatients,3,log(1/3));
  #Changed mug2 initialization (optional two lines)
  #mug2[,which(pheno>0.5),1]=log(0.3);mug2[,which(pheno>0.5),2]=log(0.35);mug2[,which(pheno>0.5),3]=log(0.35);
  #mug2[,which(pheno<0.5),1]=log(0.4);mug2[,which(pheno<0.5),2]=log(0.3);mug2[,which(pheno<0.5),3]=log(0.3);
  hsave=matrix(0,maxiter*repmax,nbgenes);


  #Precompute finite states factors
  factorQual=factorQualGen(maxsnps);
  if (length(e))factorQuant=factorQuantGen(dim(muf)[2],ebins,relbins);#standardized expression go in 3
  geneFArray=factorGeneGen();
  mx=max(meancgenes*5,50)
  factorP=factorPCommonGen(nbgenes,mx,possibleComplexity,pcase0,decay);
  sd=2;factorReg=factorRegGen(nbgenes,mx,meancgenes,sd=sd);#WARNING should test this *2 (not the real model size)
  if (length(net))factorNet=factorNetGen(netparams,max(unlist(lapply(net,length))),nbgenes,meancgenes);


  if (verbo)print("Initialization finished");
  gc()
  return(list(status=FALSE,g=NULL,h=NULL,iter=0,causes=NULL,hsave=hsave,margC=NULL,munetall=munetall,pcase0=pcase0,mureg=NULL,likelihood=NULL,predict=NULL,mux2=mux2,mux=mux,mug=mug,mug2=mug2,muq=muq,muq2=NULL,mue=mue,mue2=mue2,muf=muf,muf2=muf2,factorP=factorP,factorQual=factorQual,geneFArray=geneFArray,mx=mx,muy=muy,muh2=muh2,muh=muh,factorReg=factorReg,etaC=etaC,muC2=muC2,muC=muC))

}


addvariant=function(mes,geneid,harmfulness,varhet,varhom,pheno,nbpatients,nbsnps) {
  # Add a variant to the model, adding to dimensionality and initializing messages as appropriate
  # geneid: id of which gene variant belongs to
  # harmfulness: our prior harmfulness score
  # vals: a vector of size nbpatients of 0,1, or 2 for each patient

  # Add to mux
  mes$mux[[geneid]]<-cbind(mes$mux[[geneid]],log(rbind(1-harmfulness, harmfulness)));

  # Add to mux2
  initsnp=log(c(0.4,0.6));initsnp2=log(c(0.37,0.63));
  to_add<-array(0,dim=c(nbpatients,1,2));
  # For patients that have it, make message nonzero
  for(k in 1:nbpatients) {
	if(length(c(varhet[[k]],varhom[[k]]))) {
      if(pheno[[k]]) {
        to_add[k,c(varhet[[k]],varhom[[k]]),1]=initsnp[1];
        to_add[k,c(varhet[[k]],varhom[[k]]),2]=initsnp[2];
      } else {
        to_add[k,c(varhet[[k]],varhom[[k]]),1]=initsnp[2];
        to_add[k,c(varhet[[k]],varhom[[k]]),2]=initsnp[1];
      }
    }
  }

  # Merge in new variant to old mux2
  mes$mux2[[geneid]]<-abind(mes$mux2[[geneid]], to_add, along=2);
  
  # Add to muy
  # Will add an entry to the end of each patients matrix that has the variant
  for(k in 1:nbpatients) {
    if(length(c(varhet[[k]],varhom[[k]]))) {
	  mes$muy[[geneid]][[k]]<-rbind(mes$muy[[geneid]][[k]], matrix(-Inf, 1, 3));
	}
  }

  return mes
}

debug=function(munetall,mureg,incomingup,muh2,mug2,mug){
  print("munetall");print(munetall);
  print("mureg");print(mureg);
  print("incomingup");print(incomingup);
  print("muh2");print(muh2);
  print("mug2");print(mug2);
  print("mug");print(mug);
}


