normalize_etax= function(m){#dim=state*snps
m2=exp(m);
s=m2[1,]+m2[2,];
return(log(m2/rbind(s,s)));
}

normalize_etaC= function(m){#dim=possiblecomplexity
mx=max(m);
m2=exp(m-mx);
return(log(m2/sum(m2)));
}

normalize_etaf= function(m){#can use this for etax too if we transpose everything(just changed the rbinds by matrix)
m2=exp(m);
return(log(m2/rowSums(m2)));
}

normalize_etah= function(m){
m2=exp(m);
return(log(m2/(m2[,1]+m2[,2])));
}

logexpnormalized= function(ex){
ex2=exp(ex);
tot=apply(ex2,c(1,2),sum);
for (i in 1:(dim(ex)[3]))ex2[,,i]=log(ex2[,,i]/tot);
return(ex2)
}

logexpnormalized_2d= function(mu){
ex=exp(mu);
y=rowSums(ex);
return(log(ex/matrix(y,ncol=dim(ex)[2],nrow=dim(ex)[1])))
}

lognormalized= function(ex){
ex2=ex;
tot=ex[,,1];for (i in 2:(dim(ex)[3]))tot=tot+ex[,,i];
for (i in 1:(dim(ex)[3]))ex2[,,i]=log(ex[,,i]/tot);
return(ex2)
}

logminusmax= function(lg){
return(lg-max(lg));
}

marginal= function(mug,mug2){
x=mug+mug2;
return(logexpnormalized(x))
}

marginal_aggr=function(marg){#return non zeros (from multiple dim to 1)
x=exp(marg);
return(log( 1-x[,,1] ))
}

marginal_parallel= function(mug,mug2){
return(mapply(marginal,mug,mug2))
}

remInf= function(m){
m[m==-Inf] <- -10000;
return(m)
}

remUnstable= function(x,eps){
x[x<eps] <- 0;# 0 or eps/2
return(x)
}

remZeros= function(x,eps=10^(-17)){
last=length(x);
while (x[last]<=eps){last=last-1;}# or <eps
return(x[1:last]);
}

