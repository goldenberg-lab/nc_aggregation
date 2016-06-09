medmad= function(x,div=NULL,med=NULL){
if(!length(div)){div=mad(x);med=median(x);}
if ( ((length(which(x==min(x)))/length(x)) >0.2) | ((length(which(x==max(x)))/length(x)) >0.2))return((x-mean(x))/(sd(x)+0.01)) #correct for "many 0s problem"  or /meddis(x,mean(x))
return((x-med)/(div+0.01))
}

modmad= function(x,method="shorth"){
mod=mlv(x,method=method)$M;
return((x-mod)/meddis(x,mod));
}

modmad2= function(x,method="Tsybakov"){
mod=mlv(x,method=method,bw=NULL)$M;
return((x-mod)/meddis(x,mod));
}

meddis= function(x,m){
return(1.4826*median(abs(x-m)));
}

medtau= function(x){
return((x-median(x))/scaleTau2(x))
}

medSn= function(x){
return((x-median(x))/Sn(x))
}

medQn= function(x){
return((x-median(x))/Qn(x))
}

medrbw= function(x){
return((x-median(x))/sqrt(as.numeric(r.bw(x))))
}

less_marginal=function(x,y){if (abs(x)<abs(y)){return(x);}else return(y);}

less_marginal_mat=function(x,y){
u=x;
for (i in 1:nrow(x))for(j in 1:ncol(x))u[i,j]=less_marginal(x[i,j],y[i,j]);
return(u);
}


robust_standardize=function(x){#take an object with two matrices $exprcases $exprcontr
var1=apply(x$exprcases,1,mad);
var2=apply(x$exprcontr,1,mad);
med1=apply(x$exprcases,1,median);
med2=apply(x$exprcontr,1,median);
cases1=mapply(function(i,varj,medj)medmad(x$exprcases[i,],varj,medj),1:length(var1),var1,med1);
cases2=mapply(function(i,varj,medj)medmad(x$exprcases[i,],varj,medj),1:length(var2),var2,med2);
contr1=mapply(function(i,varj,medj)medmad(x$exprcontr[i,],varj,medj),1:length(var1),var1,med1);
contr2=mapply(function(i,varj,medj)medmad(x$exprcontr[i,],varj,medj),1:length(var2),var2,med2);
return(list(exprcases=t(less_marginal_mat(cases1,cases2)),exprcontr= t(less_marginal_mat(contr1,contr2))))
}

clean_rows=function(x,genenames){
g_in= rownames(x) %in% genenames;
na_rows=which(sapply(1:nrow(x),function(i)return(any(is.na(x[i,])) | !g_in[i])));
if(length(na_rows))return(x[-na_rows,]); 
return(x);
}

preprocess_expr_all=function(exprmatraw,mapped,ph1,ph0,genenames,nbpatients,pheno_expr,quantnormalize,logfirst,neg_controls,premethod,ksi,nfactor){
cases=match(ph1,mapped);controls=match(ph0,mapped);
indcasespresent=which(!is.na(cases));indcontrolspresent=which(!is.na(controls));
###remove lines with missing values or not in genenames. exprmat is ordered as follows: all cases then all controls
pheno_expr=c(rep(1,length(indcasespresent)),rep(0,length(indcontrolspresent)))
exprmat=clean_rows(as.matrix(exprmatraw)[,c(cases[indcasespresent],controls[indcontrolspresent])],genenames)
###Remove unwanted variance
goodneg_controls=neg_controls#find_good_negative_controls(exprmat,pheno_expr,quantnormalize,logfirst,neg_controls,premethod,ksi,20,200,genenames)
x=preprocess_expr(exprmat,pheno_expr,quantnormalize,logfirst,goodneg_controls,premethod,ksi,k=nfactor);
####Robust standardization
x2=robust_standardize(x);
e1=matrix(0,nrow(exprmat),nbpatients);rownames(e1) <- rownames(exprmat);
e1[,ph1[indcasespresent]]=x2$exprcases;#was t(apply(x$exprcases,1,medmad));
e1[,ph0[indcontrolspresent]]=x2$exprcontr;#was t(apply(x$exprcontr,1,medmad));
mapexprToGene=match(rownames(exprmat),genenames);
e=matrix(0,length(genenames),nbpatients);
e[mapexprToGene,]=e1;
return(e);
}

preprocess_expr=function(exprmat,pheno_expr,quantnormalize,logfirst,neg_controls,premethod,ksi,k=40,logafter=FALSE){
ncases=length(which(pheno_expr==1));ruv=NULL; exprgenes=rownames(exprmat);
#quantile normalize and/or log transform if needed
if(quantnormalize)exprmat=normalize.quantiles(exprmat);
if(logfirst)exprmat=log(exprmat+ksi);
#Apply another RUV method (or not)
ctl=rep(FALSE,nrow(exprmat));ctl[which(exprgenes %in% neg_controls)]=TRUE;
if (premethod=="None"){cleaned_expr2=exprmat;
}else if (premethod=="PCA"){ruv=prcomp(exprmat,center=TRUE,scale=FALSE);cleaned_expr2=ruv$x[,k:ncol(ruv$x)] %*% t(ruv$rotation[,k:ncol(ruv$x)]);
}else if (premethod=="RUVg"){ruv=RUVg(exprmat, ctl,k=k);cleaned_expr2=ruv$normalizedCounts
}else{
if (premethod=="RUVinv")ruv=RUVinv(t(exprmat), as.matrix(pheno_expr), ctl);
if (premethod=="RUVrinv")ruv=RUVrinv(t(exprmat), as.matrix(pheno_expr), ctl);
if (premethod=="RUV4")ruv=RUV4(t(exprmat), as.matrix(pheno_expr), ctl,k=k);
if (premethod=="RUV2")ruv=RUV2(t(exprmat), as.matrix(pheno_expr), ctl,k=k);
cleaned_expr=exprmat- t(ruv$W%*%ruv$alpha)
cleaned_expr2=cleaned_expr; #indneg=which(cleaned_expr<0); if (length(indneg)){cleaned_expr2[cleaned_expr2<0]<- runif(length(indneg))*ksi;print(paste("negative:",length(indneg)))}
}
exprcases=cleaned_expr2[,1:ncases];
exprcontr=cleaned_expr2[,(ncases+1):ncol(cleaned_expr2)];
####log transform after
if (logafter){exprcases=log(ksi+exprcases);exprcontr=log(ksi+exprcontr);}
return(list(exprcases=exprcases,exprcontr=exprcontr,ruv=ruv))
}

find_good_negative_controls_old=function(exprmat,pheno_expr,quantnormalize,logfirst,neg_controls,premethod,ksi,div,genenames){
u=rep(0,div);v=rep(0,div);
for (i in 1:div){
ind=which(1:length(neg_controls) %% div ==(i-1));
x=preprocess_expr(exprmat,pheno_expr,quantnormalize,logfirst,neg_controls[-ind],premethod,ksi);
e1=matrix(0,nrow(exprmat),ncol(exprmat));rownames(e1) <- rownames(exprmat);
e1[,1:ncol(x$exprcase)]=t(apply(x$exprcases,1,medmad));
e1[,(ncol(x$exprcase)+1):ncol(exprmat)]=t(apply(x$exprcontr,1,medmad));
mapexprToGene=match(rownames(exprmat),genenames);
e=matrix(0,length(genenames),ncol(exprmat));
e[mapexprToGene,]=e1;
u[i]=length(which(abs(e)>3));
v[i]=length(which(abs(e)>2));
}
good_bins=order(4*u+v)[1:(div/2)];
return(neg_controls[which( ((1:length(neg_controls) %% div) +1) %in% good_bins)])
}

find_good_negative_controls=function(exprmat,pheno_expr,quantnormalize,logfirst,neg_controls,premethod,ksi,div,n,genenames){
sampling=matrix(0,div,length(neg_controls));
u=rep(0,div);v=rep(0,div);
for (i in 1:div){
ind=sample(1:length(neg_controls),n);
sampling[i, ind]=1;
x=preprocess_expr(exprmat,pheno_expr,quantnormalize,logfirst,neg_controls[ind],premethod,ksi);
e1=matrix(0,nrow(exprmat),ncol(exprmat));rownames(e1) <- rownames(exprmat);
e1[,1:ncol(x$exprcase)]=t(apply(x$exprcases,1,medmad));
e1[,(ncol(x$exprcase)+1):ncol(exprmat)]=t(apply(x$exprcontr,1,medmad));
mapexprToGene=match(rownames(exprmat),genenames);
e=matrix(0,length(genenames),ncol(exprmat));
e[mapexprToGene,]=e1;
u[i]=length(which(abs(e)>3));
v[i]=length(which(abs(e)>2));
}
y=sampling*(4*u+v); y[y==0]<- NA;
good_ctr=order(colMeans(y,na.rm=TRUE))[1:n];
return(neg_controls[good_ctr])
}

plot_expr_gene=function(file,g,genenames,exprmatraw,e,ph1,ph0,includedpatientids){
png(file)
x=exprmatraw[which(rownames(exprmatraw)==genenames[g]),]
par(mfrow=c(2,2));
hist(log(as.numeric(x[which(names(x)%in% names(includedpatientids[ph1]))])+1),main="Cases");hist(log(as.numeric(x[which(names(x)%in% names(includedpatientids[ph0]))])+1),main="Controls");
hist(e[g,ph1],main="Cases after robust standardization");hist(e[g,ph0],main="Controls after robust standardization");
dev.off();
}


