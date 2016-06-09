
compute_etax_parallel=function(mux2,mux,het,hom){
  if(length(mux)==0)return(NULL);
  mux2=remInf(mux2);
  if(length(mux)==1){all=apply(mux2,c(2,3),sum);}else{ all=apply(mux2,c(2,3),sum);}
  all=all-rowMeans(all);
  #etax=array(0, dim=c(dim(mux2[[1]]),2));#simple_triplet_zero_matrix(dim(mux2[[1]])[1], dim(mux2[[1]])[2], mode = "double")
  etax=mapply(compute_etax,1:(dim(mux2)[1]), het,hom, MoreArgs=list(mux2=mux2,all=all,mux=mux),SIMPLIFY="array",USE.NAMES = FALSE)
  return(aperm(etax,c(3,1,2)))
}

compute_muy_parallel= function(het,hom,muy,etax){
  return(mapply(compute_muy,1:length(het),het,hom,muy,MoreArgs=list(etax=etax),SIMPLIFY=FALSE));
}

compute_muq_parallel= function(muy,het,hom,factorQual){
  return(mapply(compute_muq,muy,het,hom,MoreArgs=list(factorQual=factorQual),SIMPLIFY=FALSE));#### Was compute_muq_general
}

compute_etaf_parallel=function(muf2,muf){
  muf2=remInf(muf2)
  all=muf+apply(muf2,c(1,3),sum);
  all=all-apply(all,1,max);
  res= mapply(compute_etaf,1:(dim(muf2)[2]),MoreArgs=list(muf2=muf2,all=all),SIMPLIFY="array");
  return(aperm(res,c(1,3,2)))
}

compute_etaC_parallel= function(muC2,muC){
  muC2=lapply(muC2,remInf);
  all=muC+Reduce('+',muC2);
  res= lapply(1:length(muC2),compute_etaC,muC2=muC2,all=all);
  return(res);
}

compute_mue_parallel= function(exetaf,ed,factorQuant){#for one gene
  return(mapply(compute_mue,listbydim12(exetaf),ed,MoreArgs=list(factorQuant=factorQuant)));
}

compute_muf2_parallel= function(exmue2,ed,factorQuant){#for one gene
  return(mapply(compute_muf2,listbydim12(exmue2),ed,MoreArgs=list(factorQuant=factorQuant)));
}

compute_mug_parallel= function(muq, mue, etah){
  return(mapply(compute_mug,muq, mue, etah,SIMPLIFY=FALSE));
}

compute_muq2_parallel= function(mug2, mue, etah, propagate){
  return(mapply(compute_muq2,mug2, mue, etah, MoreArgs=list(propagate=propagate),SIMPLIFY=FALSE));
}

compute_mue2_parallel= function(mug2, muq, etah, propagate){
  return(mapply(compute_mue2,mug2, muq, etah, MoreArgs=list(propagate=propagate), SIMPLIFY=FALSE));
}

compute_muh2_parallel= function(mug2, muq, mue){
  return(mapply(compute_muh2,mug2, muq, mue, SIMPLIFY=FALSE));
}

compute_muy2_parallel= function(muq2, muy, het, hom, factorQual){
  return(mapply(compute_muy2,muq2, muy,het, hom, MoreArgs=list(factorQual=factorQual),SIMPLIFY=FALSE));#was compute_muy2_general
}

compute_mux2_parallel= function(muy2,het,hom,mux2){
  return(mapply(compute_mux2,muy2, het, hom, MoreArgs=list(mux2=mux2), SIMPLIFY=FALSE));
}

