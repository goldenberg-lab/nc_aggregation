
hashed=function(str,val){
e=new.env(size=2*length(str),hash=TRUE);
mapply(function(s,v)e[[s]]<- v,str,val);
return(e)
}

get_hashed=function(str,hasht){
return(sapply(str,function(x){y=hasht[[x]];if(!length(y))return(NA);return(y)}))
}

verifyDependencies=function(phenomat,genotype,annot,genemat,exprmatraw){
#verify that patient ids are the same/intersect in clinical,gene expression and exome
#verify that variant ids are the same/intersect in annotation (variants) and exome
#verify that gene names intersect between the annotation, expression and optional gene names list
dataok=1;
if (length(exprmatraw)){ 
intersectpatient1=intersect(phenomat[,1],colnames(exprmatraw));
if (length(intersectpatient1)==0){print("Problem with Patients Ids: The ids in the expression data are entirely different from the ids in the clinical data.");dataok=0;}
diff=length(phenomat[,1])-length(intersectpatient1);if (diff>0.5*length(phenomat[,1]))print(paste("Warning : more than half of the patients in the clinical data are not found in the expression data."))
}
intersectpatient2=intersect(phenomat[,1],unique(genotype[,2]));
if (length(intersectpatient2)==0){print("Problem with Patients Ids: The ids in the exome data are entirely different from the ids in the clinical data.");dataok=0;}
diff=length(phenomat[,1])-length(intersectpatient2);
if (diff>0.5*length(phenomat[,1])){print(paste("Problem : more than half of the patients in the clinical data are not found in the exome data."));dataok=0;} else{
if (diff>0){print(paste("Problem :",diff,"of the patients in the clinical data are not found in the exome data."));dataok=0;}}

intersectvariant=intersect(unique(annot[,1]),unique(genotype[,1]));
if (length(intersectvariant)==0){print("Problem with Variants Ids: The ids in the exome data are entirely different from the ids in the annotation file");dataok=0;}
diff=length(unique(genotype[,1]))-length(intersectvariant);if (diff)print(paste("Warning :",diff,"variants are found in the exome data but not in the annotation."))

if (length(exprmatraw)){ 
intersectgene1=intersect(genemat[,1],rownames(exprmatraw));
if (length(intersectgene1)==0){print("Problem with Gene names: The names in the expression data are entirely different from the ones in the exome/given list of genes");dataok=0;}
}
intersectgene2=intersect(genemat[,1],unique(annot[,2]));
if (length(intersectgene2)==0){print("Problem with Gene names: The names in the variants annotation are entirely different from the ones in the given list of genes");dataok=0;}

if (any(is.na(annot[,4]))){print("Problem with Variants annotation: Some harmfulness predictions have missing values");dataok=0;}

return(dataok)
}
