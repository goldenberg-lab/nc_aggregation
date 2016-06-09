simul_polyphen= function(ent){
harm=ent[1];s=ent[2];
x=0;
if (harm){
  if (s<0.0005) {x=rnorm(1,0.2,0.3);
  }else if (s<0.001) {x=rnorm(1,0.3,0.3);
  }else if (s<0.005) {x=rnorm(1,0.6,0.2);
  }else if (s<0.01) {x=rnorm(1,0.7,0.2);
  }else if (s<0.02) {x=rnorm(1,0.8,0.2);
  }else if (s<0.05) {x=rnorm(1,0.9,0.2);
  }else x=rnorm(1,1,0.2);

}else{
x=rnorm(1,0.15,0.2);
}

return( max(0,min(x,1)) )
}

simul_polyphen2= function(ent,distr1,distr2){
harm=ent[1];s=ent[2];
x=0;
if (harm){
x=sample(distr1,1);
}else{
x=sample(distr2,1);
}

return( max(0,min(x,1)) )
}

