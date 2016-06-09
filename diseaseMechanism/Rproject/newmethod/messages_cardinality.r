
comupute_mu_cardinality=function(mu,potential,init=rep(1/ncol(mu), ncol(mu)) ){
mu2exp=matrix(0,nbgenes,ncol(mu));
for (i in 1: ncol(mu))mu2exp[,i]=init[i];

nonzer=which(mu[,1]!=1);
n=length(nonzer);
if (n==0){return(log(mu2exp));}
if (n==1){mu2exp[nonzer[1],]=potential[1:ncol(mu2exp)];return(log(mu2exp/rowSums(mu2exp)));}
#counting number of internal messages
lp=n-1;
b=max(2^max(0,ceiling(log2(n) -7) ),1); b=1; #if (mx<501)b=1; #was mx <200 #Avoiding the problem

l1=compute_mu1_cardinality(mu,lp,n,nonzer,b);
l2=compute_mu2_cardinality(l1,potential[1:length(l1[[1]])],lp,n,nonzer,b)
default=compute_convolution2(c(l2[[1]],rep(0,ncol(mu2exp)-1)),l1[[1]],ncol(mu2exp));
for (i in 1:ncol(mu2exp))mu2exp[,i]=default[i];

for (i in 1:n)mu2exp[nonzer[i],]=l2[[lp+i]];#WARNING verify each message is of the correct size 

return(log(mu2exp/rowSums(mu2exp)));
}

compute_mu1_cardinality=function(mu,lp,n,nonzer,b){
l1=list();length(l1)<- lp+n;
#forward messages
for(i in 1:n)l1[[lp+i]]<- mu[nonzer[i],];
for(i in lp:b){l1[[i]]= compute_convolution1(l1[[2*i]],l1[[2*i+1]]);}
if(b>1){for(i in (b-1):1){l1[[i]]= compute_convolution1(l1[[2*i]],l1[[2*i+1]]);}}

return(l1)
}

compute_mu1_cardinality_all=function(mu,lp,n){
l1=list();length(l1)<- lp+n;
#forward messages
for(i in 1:n)l1[[lp+i]]<- mu[i,];
for(i in lp:1){l1[[i]]= compute_convolution1(l1[[2*i]],l1[[2*i+1]]);}
return(l1)
}

compute_mu2_cardinality=function(l1,potential,lp,n,nonzer,b){
l2=list();length(l2)<- lp+n;
#backward messages
l2[[1]]=potential;
ind=1:(lp+n); inddiv=ind%/%2; indpair=ind+1-2*(ind%%2);
if(b>2)for (i in 2:b){l2[[i]]= compute_convolution2(l2[[inddiv[i]]],l1[[indpair[i]]],length(l1[[i]])); }
for (i in (b+1):(lp+n)){l2[[i]]= compute_convolution2(l2[[inddiv[i]]],l1[[indpair[i]]],length(l1[[i]]));} 

return(l2)
}

compute_mu2_cardinality_all=function(l1,potential,lp,n){
l2=list();length(l2)<- lp+n;
#backward messages
l2[[1]]=potential;
ind=1:(lp+n); inddiv=ind%/%2; indpair=ind+1-2*(ind%%2);
for (i in 2:(lp+n)){l2[[i]]= compute_convolution2(l2[[inddiv[i]]],l1[[indpair[i]]],length(l1[[i]]));} 
return(l2)
}


compute_convolution1=function(m1,m2){#m2 should be smaller
s=rep(0,length(m1)+length(m2)-1);
l1=get_last(m1);
for (i in 1:get_last(m2)){ind=i:(i+l1-1); s[ind]=s[ind]+m1[1:l1]*m2[i];}##CHANGED
return(s);
}

compute_convolution1_ftt=function(m1,m2){
return(convolve(m1,rev(m2),type="o"));
}

compute_convolution2=function(m1,m2,l){#m2 should be smaller (m1 is coming from the sum node)
s=rep(0,l);
l2=get_last(m2);l1=min(get_last(m1),l);
s[1:l1]=vapply(1:l1,function(i)sum(m1[i:(i+l2-1)]*m2[1:l2]),c(1),USE.NAMES = FALSE)
return(s);
}

compute_convolution2_fft=function(m1,m2){#m2 should be smaller (m1 is coming from the sum node)
return(convolve(m1,m2,type="o"));
}

get_last=function(x,eps=10^(-18)){
last=length(x);
while (x[last]<=eps){last=last-1;}# or <eps
return(last);
}
