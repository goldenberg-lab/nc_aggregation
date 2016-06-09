harmonize_expr =function(data,verbose){
m=rep(0,dim(data)[2])
data2=data;
for(i in 1:(dim(data)[2])){
d=density(data[,i]);
m[i]=d$x[which.max(d$y)];
}
avg=mean(m);
if (verbose)print(m);
for(i in 1:(dim(data)[2])){
data2[,i]=data[,i]*avg/m[i];
}
return(data2);
}
