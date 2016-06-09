#veriables needed from mainSimulNew
#dirrep,healthy_expr,networkName (not codeDir)
#genes,nbgenes,maxConnections,maxneighborhoodsize,simmaxmut
#nb_genes_causal,sample_size,nbPatients, ratioSignalSim
#freeMem,allmut,allexp,allexpMagnitude,sdamage,convertToGene,allmap, dict,ps, coef


candidates=c();
net=load_network_genes(networkName,genes,maxConnections)$network;

surrounding=mapply(surround2,net,1:length(net) , MoreArgs=list(net1=net));
neighborhoodsizes=unlist(lapply(surrounding,length))
acceptable=which(neighborhoodsizes>=nb_genes_causal & neighborhoodsizes<= maxneighborhoodsize)
index_causal_one=sample(acceptable,1,replace=F);
candidates=surrounding[[index_causal_one]];

index_causal=c(index_causal_one,sample(candidates,nb_genes_causal-1,replace=F));
index_causal_genes=index_causal;#index in expr / harm data

print(paste("Causal genes:",paste(genes[index_causal_genes], collapse=" "),"indexes:",paste(index_causal_genes, collapse=" ")));
write.table(cbind(c(index_causal_genes),c(genes[index_causal_genes])), paste(dirrep,"truth.txt",sep=""),col.names=c("Causal-genes-index", "genes"),row.names=FALSE);
length(candidates);
truth= index_causal_genes# Select a group of causal genes

indivScores=mapply(simul_pheno_number,allmut,allexp,allexpMagnitude,MoreArgs=list(truth,sdamage,convertToGene,allmap, dict))

dead=length(which(sort(colSums(indivScores),decreasing=TRUE)>simmaxmut));
nbPatientsSignal=round(ratioSignalSim*nbPatients); 
nbPatientsNoSignal=nbPatients-nbPatientsSignal;
indPatients=order(colSums(indivScores),decreasing=TRUE)[(dead+1):(dead+nbPatientsSignal)];
if (nbPatientsNoSignal)indPatients=c(indPatients, sample(((dead+1):nbindiv)[-indPatients],nbPatientsNoSignal))
print(paste("Contributions from SNPs/Expression",paste(rowSums(indivScores[,indPatients]),collapse=" ")))
contributions[,z]=rowSums(indivScores[,indPatients]);
indControls=sample(((dead+1):nbindiv)[-indPatients],nbPatients);
indPatients=sample(indPatients,sample_size);
indControls=sample(indControls,sample_size);

#patient/control snp matrix
indivScores_Qual=mapply(simul_pheno_number_Qual,allmut[c(indControls,indPatients)],MoreArgs=list(truth,sdamage,convertToGene,allmap, dict))
indivScores_Quant=mapply(simul_pheno_number_Quant,allexp[c(indControls,indPatients)],allexpMagnitude[c(indControls,indPatients)],MoreArgs=list(truth,sdamage,convertToGene,allmap, dict))
trans=simul_snpsToMatrix(allmut,indPatients,indControls,convertToGene,ps,genes,dict);
patientnames=paste("Indiv_",c(indPatients,indControls),sep="");#same patient order as rows in trans

#cleaning previous content just in case
rem=file.remove(paste(dirrep,list.files(path=dirrep,pattern="^snps"),sep=""))

#preprocess Expression background
x=normalize.quantiles(healthy_expr);
expr_bg=log(x);

res=simul_expr_background(expr_bg,nbgenes,sample_size,sample_size,FALSE);#WARNING the last argument means I am simulating gaussian. Ramove if I want sampling from real
res2=simul_expr_change(res,indPatients,indControls,allexp,allexpMagnitude,coef);
write.table(cbind(res2$pat_expr,res2$ctr_expr),paste(dirrep,"gene_expression.txt",sep=""),sep="\t",row.names=genes,col.names=patientnames)


exomefile=paste(dirrep,"genotype.txt",sep="");file.create(exomefile);
for (i in 1:ncol(trans$mat)){
patind=which(trans$mat[,i]>0);
write(paste(trans$snps[i],patientnames[patind],trans$mat[patind,i],sep="\t"),exomefile, append=TRUE,sep="\t",ncolumns=1);
}

annot=read.table(paste(dir,"annot.txt",sep=""),sep="\t",colClasses=c("character","character","numeric"));
presentsnps=which(annot[,2]%in% trans$snps);
write.table( data.frame(patientnames,c(rep(1,length(indPatients)),rep(0,length(indPatients)))),paste(dirrep,"patients.txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE);
write.table( data.frame(annot[presentsnps,2:1],"missense",annot[presentsnps,3]),paste(dirrep,"variants.txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=c(1,2,3));
write.table( genes, paste(dirrep,"genes.txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE);
write.table( trans$mat,paste(dirrep,"trans_snp_matrix.txt",sep="") ,sep="\t",row.names=FALSE,col.names=paste(trans$gene,trans$snps,sep="__"));


#freeing memory
if (freeMem){allmut=NULL;allexpMagnitude=NULL; allexp=NULL;}
gc();#print(sum( sapply(mget(ls(),envir = parent.frame()),object.size) ));

