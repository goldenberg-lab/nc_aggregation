load_network= function(fileName){

data= read.table(fileName,stringsAsFactors =FALSE, sep="\t",na.strings = "null",header=TRUE);
print("file loaded.");
#build a list of gene names and their ranks and put them in a hashtable
hashenv<-new.env()
j=1;
genes=c();
for(i in 1:dim(data)[1]){
	if ( is.null(hashenv[[ data[i,1] ]]) ){  
	hashenv[[ data[i,1] ]]<- j;
	genes=c(genes,data[i,1]);
	j=j+1;
	}
        if ( is.null(hashenv[[ data[i,2] ]]) ){  
	hashenv[[ data[i,2] ]]<- j;
	genes=c(genes,data[i,2]);
	j=j+1;
	}
}
print(paste("There are", length(genes), "genes in the network."));
network <- list();
for (i in 1:length(genes))network[[i]] <- 0;

for(i in 1:dim(data)[1]){
u= hashenv[[data[i,1]]];
v= hashenv[[data[i,2]]];
if(network[[u]][1]==0){network[[u]]<- v;
}else if (!(v %in% network[[u]]))network[[u]]<- c(network[[u]],v);
if(network[[v]][1]==0){network[[v]]<- u;
}else if (!(u %in% network[[v]]))network[[v]]<- c(network[[v]],u);
}

print("Network constructed.");
return (list(network=network,genes=genes));
}

surrounding= function(net, rootg , candidates, K, maxConnections){

    distance_one = net[[rootg]];
    genes_to_check = c(distance_one);
    if (K>1){
    for(i in 1:length(distance_one)){
        print(distance_one[i])
        if (length(net[[distance_one[i]]])<maxConnections)genes_to_check = c(genes_to_check, surrounding(net,distance_one[i],candidates,K-1, maxConnections) );
    }
    }
    genes_to_check = unique(genes_to_check);
#remove rootg
x=which(genes_to_check!=rootg) 

return(genes_to_check[x]);
}
