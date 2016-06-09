load_network_genes= function(fileName,genenames,rem=0){

data= read.table(fileName,stringsAsFactors =FALSE, sep="\t",na.strings = "null",header=TRUE);
print("file loaded.");
data2=unique(data[,c(1,2)]);
if(length(genenames)==0)genenames=unique(sort(c(data[,1],data[,2])));
l=list();length(l)<- length(genenames);
#build a list of gene names and their ranks and put them in a hashtable

u=match(data2[,1],genenames);
v=match(data2[,2],genenames);

for(i in 1:dim(data2)[1]){
if(!is.na(u[i]) && !is.na(v[i]) && u[i]!=v[i]){l[[u[i]]]=c(l[[u[i]]],v[i]); l[[v[i]]]=c(l[[v[i]]],u[i]);}
}
n=0;for(i in 1:length(genenames))if(length(l[[i]])){l[[i]] <- unique(sort(l[[i]]));n=n+length(l[[i]]); }

p=0;for(i in 1:length(genenames))if(length(l[[i]])>0)p=p+1;

if (rem){
 s=unlist(lapply(l,length));
 si=which(s>rem);
 for (i in si){
  for (j in l[[i]]){
    if(length(l[[j]])>1){l[[j]]=l[[j]][-which(l[[j]]==i)];} else {l[j] <- list(NULL);p=p-1;}
  } 
  n=n-2*length(l[[i]]); l[i] <- list(NULL);p=p-1;
 }
}

print(paste("Network constructed. Contains ",n ,"interactions involving ",p," of the genes"));
return (list(network=l,genes=genenames));
}

surround2= function(x,i,net1,nofirstorder=FALSE){
 if (length(x)){
   nett=unique(sort(unlist( net1[x] ))); 
   self=which(nett==i); 
   if (nofirstorder) self=unique(c(self,which(nett %in% net1)));
   if (length(self))nett=nett[-self];
   return(nett);
 } else{ return(NULL);}
}


