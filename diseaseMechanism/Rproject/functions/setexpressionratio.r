nb_genes_causal=5;
nbPatients=3000;
ratioExprSnps=c(0.66,1,1.5,4);#(0.5,1,1.5,2,3)for 10genes give (0.33 0.52 0.66 0.73 0.71) expr 50 genes0.35 0.54  20genes  0.36 0.526 0.65 0.708 0.773   5genes 0.35 0.52 0.647 0.66 0.756
vals=matrix(0,2,length(ratioExprSnps));k=1;
nbrep=20;
for (ratioExprSnp in ratioExprSnps){
contributions=matrix(0,2,nbrep);
for (z in 1:nbrep){

#Simulate expression perturbations
corrector=0.9#unbiasRatio(nb_genes_causal);#ratio is 1 if model size 5; 0.5 md10; 3.5 md2; For 66% 1.3 md2/5
nbExprAvg=corrector*ratioExprSnp*allmutAvgLength;nbExprSd=corrector*ratioExprSnp*allmutSdLength
nbindiv=length(allmut)
allexp <-list();length(allexp) <-nbindiv;
allexpMagnitude <- list();length(allexpMagnitude) <-nbindiv;
genesID=1:nbGenes
allexp <-lapply(allexp,simul_sample,genesID,nbExprAvg,nbExprSd)
allexpMagnitude <- lapply(allexp,simul_unif)

maxConnections=100;
candidates=c();
n=0;
index_causal_one=1;
while(length(candidates)<nb_genes_causal){
	index_causal_one=sample(genesToNet,1,replace=F);
	surround=surrounding(net$network,index_causal_one,c(),2,maxConnections);
	candidates=surround[which(surround %in% genesToNet)];
	n=n+1;
}
index_causal=c(index_causal_one,sample(candidates,nb_genes_causal-1,replace=F));
index_causal_genes=netToGenes[index_causal];#index in expr / harm data

#print(paste("n=",n,". Causal genes:",paste(genes[index_causal_genes], collapse=" "),"indexes:",paste(index_causal_genes, collapse=" ")));
model= index_causal_genes# Select a group of causal genes

indivScores=mapply(simul_pheno_number,allmut,allexp,allexpMagnitude,MoreArgs=list(model,sdamage,lGene,allmap, dict))
indPatients=order(colSums(indivScores),decreasing=TRUE)[1:nbPatients]
contrib=rowSums(indivScores[,indPatients])
contributions[,z]=contrib/sum(contrib);
}
vals[,k]=rowMeans(contributions);
k=k+1;
}

#Trying a few standardizations
e=t(apply(patExpr,1,medmad));e1=t(apply(ctrExpr,1,medmad));
#e=t(apply(patExpr,1,medQn));e1=t(apply(ctrExpr,1,medQn));
#e=t(apply(patExpr,1,medrbw));e1=t(apply(ctrExpr,1,medrbw));
#e=t(apply(patExpr,1,modmad));e1=t(apply(ctrExpr,1,modmad));
#e=t(apply(patExpr,1,modmad2));e1=t(apply(ctrExpr,1,modmad2));
#check a few things about e
thres=2.;
abExpr=apply(abs(e),1,function(e){return(length(which(abs(e) >thres) ))});
abExpr1=apply(abs(e1),1,function(e1){return(length(which(abs(e1) >thres) ))});
print(sum(abExpr[truth]));print(sum(abExpr1[truth]));



#test the Standardization approaches

rat=0.1;coef=12;
nbExprAvg2=rat*nbgenes;nbExprSd2=nbExprSd*nbExprAvg2/nbExprAvg;
allexp <-list();length(allexp) <- sample_size*2;
allexpMagnitude <- list();length(allexpMagnitude) <- sample_size*2;
genesID=1:nbGenes
allexp <-lapply(allexp,simul_sample,genesID,nbExprAvg2,nbExprSd2)
start=0.5;end=0.8; # perturbation is unif(start,end)*coef*sd *sign
allexpMagnitude <- lapply(allexp,simul_unif,start=start,end=end)

res=simul_expr_background(expr_bg,nbGenes,sample_size,sample_size,FALSE);
res2=simul_expr_change(res,1:sample_size,(sample_size+1):(sample_size*2),allexp,allexpMagnitude,coef);
m1=rowMeans(res$pat_expr);m11=rowMeans(res$ctr_expr);
s1=apply(res$pat_expr,1,sd);s11=apply(res$ctr_expr,1,sd);
m2=rowMeans(res2$pat_expr);m21=rowMeans(res2$ctr_expr);
s2=apply(res2$pat_expr,1,sd);s21=apply(res2$ctr_expr,1,sd);
m2med=apply(res2$pat_expr,1,median);m21med=apply(res2$ctr_expr,1,median);
s2med=apply(res2$pat_expr,1,mad);s21med=apply(res2$ctr_expr,1,mad);
m2mod=apply(res2$pat_expr,1,mod);m21mod=apply(res2$ctr_expr,1,mod);
s2mod=apply(res2$pat_expr,1,modsd);s21mod=apply(res2$ctr_expr,1,modsd);
#m2mod2=apply(res2$pat_expr,1,mod2);m21mod2=apply(res2$ctr_expr,1,mod2);
#s2mod2=apply(res2$pat_expr,1,modsd2);s21mod2=apply(res2$ctr_expr,1,modsd2);

indx=1:200#length(m1);#truth
mserm= mean(((m1-m11)^2)[indx]);msers= mean(((s1-s11)^2)[indx]);
print(c(mserm,msers))
msem1= mean(((m1-m2)^2)[indx]);mses1= mean(((s1-s2)^2)[indx]);
print(c(msem1,mses1));
msemmed=mean(((m1-m2med)^2)[indx]);msesmed=mean(((s1-s2med)^2)[indx]);
print(c(msemmed,msesmed));
msemmod=mean(((m1-m2mod)^2)[indx]);msesmod=mean(((s1-s2mod)^2)[indx]);
print(c(msemmod,msesmod));
msemmod2=mean(((m1-m2mod2)^2)[indx]);msesmod2=mean(((s1-s2mod2)^2)[indx]);
print(c(msemmod2,msesmod2));

palette=c("black","red","blue")
png(paste("sd_estimation_",rat,"_",coef,".png",sep=""))
plot(s1[indx],s2[indx],col=palette[1],xlab="sd before perturbation",ylab="estimated sd after perturbation");abline(coef=c(0,1));
points(s1[indx],s2med[indx],col=palette[2]);
points(s1[indx],s2mod[indx],col=palette[3]);points(s1[indx],s2mod2[indx],col=palette[4]);
legend(0.8,legend=c("sd","med_sd","mod_sd"),col=palette)
dev.off()

mod= function(x,method="Venter"){
return(mlv(x,method=method)$M);
}

modsd= function(x){
mod=mod(x);
return(meddis(x,mod));
}

mod2= function(x,method="Asselin"){
return(mlv(x,method=method)$M);
}

modsd2= function(x){
mod=mod2(x);
return(meddis(x,mod));
}

meddis= function(x,m){
return(1.4826*median(abs(x-m)));
}truth


