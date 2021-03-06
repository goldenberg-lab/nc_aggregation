# Create our macro for merging 2 regions
# Works similarly to the one for adding a variant, but just takes in 2 genes
mergevariants.db <- defmacro(gene1, gene2, expr={
		                     harm[[gene1]] <- c(harm[[gene1]], harm[[gene2]]);
							 harm <- harm[-gene2];
							 h3k[[gene1]] <- c(h3k[[gene1]], h3k[[gene2]]);
							 h3k <- h3k[-gene2];
							 new.het <- vector("list", nbpatients);
							 new.hom <- vector("list", nbpatients);
               for(i in 1:nbpatients) {
				 if(length(het[[gene1]][[i]]) + length(het[[gene2]][[i]]) > 0) {
					 het[[gene1]][[i]] <- c(het[[gene1]][[i]], (het[[gene2]][[i]] + nbsnps[[gene1]]));
				 }
			     if(length(hom[[gene1]][[i]]) + length(hom[[gene2]][[i]]) > 0) {
					 hmo[[gene1]][[i]] <- c(hom[[gene1]][[i]], (hom[[gene2]][[i]] + nbsnps[[gene1]]));
				 }
               }
			   het <- het[-gene2];
			   hom <- hom[-gene2];
			   nbsnps[[gene1]] <- nbsnps[[gene1]] + nbsnps[[gene2]];
			   nbsnps <- nbsnps [-gene2];
			   nbgenes <- nbgenes - 1;
               mes <- merge_regions(mes,nbpatients,nbgenes,nbsnps,gene1,gene2);
});


# Create our macro for adding variants
# It takes in a geneid we are adding a variant to, the harmfulness score for it, the acetylation value, and the status of that variant for each patient (0/1/2)
# It does the actual addition in the model, but also handles adding everything into our database
addvariant.db <- defmacro(geneid, new_harm, new_h3k, vals, expr={
						  nbsnps[[geneid]] <- nbsnps[[geneid]] + 1;
						  harm[[geneid]] <- c(harm[[geneid]], new_harm);
						  h3k[[geneid]] <- c(h3k[[geneid]], new_h3k);
						  for(i in 1:length(vals)) {
							  if(vals[[i]] == 1) {
								  het[[geneid]][[i]] <- c(het[[geneid]][[i]], nbsnps[[geneid]]);
							  } else if(vals[[i]] == 2) {
								  hom[[geneid]][[i]] <- c(hom[[geneid]][[i]], nbsnps[[geneid]]);
							  }
						  }
						  mes<-addvariant(mes,geneid,new_harm,vals,pheno,nbpatients);
});


addvariant=function(mes,geneid,harmfulness,vals,pheno,nbpatients) {
  # Add a variant to the model, adding to dimensionality and initializing messages as appropriate
  # geneid: id of which gene variant belongs to
  # harmfulness: our prior harmfulness score
  # vals: numerical of length nbpatients with 0,1,2
  # pheno: The phenotypes of all patients
  

  library("abind");

  # Add to mux
  mes$mux[[geneid]]<-cbind(mes$mux[[geneid]],log(rbind(1-harmfulness, harmfulness)));

  # Add to mux2
  initsnp=log(c(0.4,0.6));initsnp2=log(c(0.37,0.63));
  to_add<-array(0,dim=c(nbpatients,1,2));
  # For patients that have it, make message nonzero
  for(k in 1:nbpatients) {
	if(vals[[k]] > 0) {
      if(pheno[[k]]) {
        to_add[k,1,1]=initsnp[1];
        to_add[k,1,2]=initsnp[2];
      } else {
        to_add[k,1,1]=initsnp[2];
        to_add[k,1,2]=initsnp[1];
      }
    }
  }

  # Merge in new variant to old mux2
  mes$mux2[[geneid]]<-abind(mes$mux2[[geneid]], to_add, along=2);
  
  # Add to muy
  # Will add an entry to the end of each patients matrix that has the variant
  for(k in 1:nbpatients) {
    if(vals[[k]]) {
	  mes$muy[[geneid]][[k]]<-rbind(mes$muy[[geneid]][[k]], matrix(-Inf, 1, 3));
	}
  }

  return(mes);
}


merge_regions=function(mes,nbpatients,nbgenes,nbsnps,reg1,reg2) {
  library("abind");
  # Merge two regions in the model, as given by indices reg1 and reg2

  # Need to handle the follow messages changing properly: mux, muy, mux2, muq, muq2, mug, mug, mug2, muh, muh2, factorQual, mue, mue2
  

  # Have to reconstruct some factors to adjust the size
  mes$factorP<-factorPCommonGen(nbgenes,mes$mx,mes$possibleComplexity,mes$pcase0,mes$decay);
  sd=2;
  mes$factorReg<-factorRegGen(nbgenes,mes$mx,mes$meancgenes,sd=sd);#WARNING should test this *2 (not the real model size)

  # We're not using expression, so we're just going to squash mue, mue2
  mes$mue <- mes$mue[1:nbgenes,,];
  mes$mue2 <- mes$mue2[1:nbgenes,,];

  # Create new factorQual
  mes$factorQual <- factorQualGen(max(nbsnps));

  # Not using the network either, so we will fix those dimensionalities
  mes$munet=list();length(mes$munet) <- nbgenes;
  mes$munetall=matrix(0,nbgenes,2);

  # mux: just need to collapse harmfulness scores, remove second region
  mes$mux[[reg1]]<-cbind(mes$mux[[reg1]],mes$mux[[reg2]]);
  mes$mux<-mes$mux[-reg2];

  # muy: just need to concatenate hits for each patient, remove second region
  for(i in 1:nbpatients) {
    mes$muy[[reg1]][[i]] <- rbind(mes$muy[[reg1]][[i]],mes$muy[[reg2]][[i]]);
  }
  mes$muy <- mes$muy[-reg2];

  # mux2: need to concatenate variants for our new region, remove second region
  mes$mux2[[reg1]] <- abind(mes$mux2[[reg1]], mes$mux2[[reg2]], along=2);
  mes$mux2 <- mes$mux2[-reg2];

  # muq: As a starting point for our new belief, we'll do this based on previous 2 regions:
  # 1,1 or 0.5,1 -> 1, 0.5, 0.5 or 0,1 -> 0.5, 0,0 or 0.5,0 -> 0
  new_muq <- array(0,dim=c(1, nbpatients,3));
  new_muq[,,1] <- (mes$muq[reg1,,1] * mes$muq[reg2,,1]) +
  				 (mes$muq[reg1,,1] * mes$muq[reg2,,2]) +
				 (mes$muq[reg1,,2] * mes$muq[reg2,,1]);
  new_muq[,,2] <- (mes$muq[reg1,,2] * mes$muq[reg2,,3]) + 
  				 (mes$muq[reg1,,3] * mes$muq[reg2,,1]) +
	 			 (mes$muq[reg1,,2] * mes$muq[reg2,,2]);
  new_muq[,,3] <- (mes$muq[reg1,,2] * mes$muq[reg2,,3]) + 
  				 (mes$muq[reg1,,3] * mes$muq[reg2,,2]) + 
				 (mes$muq[reg1,,2] * mes$muq[reg2,,2]);
  mes$muq <- abind(mes$muq[0:(reg1-1),,],new_muq,mes$muq[(reg1+1):(nbgenes+1),,],along=1);
  print(dim(mes$muq));
  mes$muq <- mes$muq[-reg2,,];

  # muq2: We do the same kind of thing as in muq
  new_muq2 <- array(0,dim=c(nbpatients,3));
  # We need these in matrices so we can access columns properly
  matrix_muq2_reg1 <- lapply(mes$muq2[[reg1]], as.matrix);
  matrix_muq2_reg1 <- t(do.call(cbind, matrix_muq2_reg1));
  matrix_muq2_reg2 <- lapply(mes$muq2[[reg2]], as.matrix);
  matrix_muq2_reg2 <- t(do.call(cbind, matrix_muq2_reg2));

  new_muq2[,1] <- (matrix_muq2_reg1[,1] * matrix_muq2_reg2[,1]) +
  				  (matrix_muq2_reg1[,1] * matrix_muq2_reg2[,2]) +
				  (matrix_muq2_reg1[,2] * matrix_muq2_reg2[,1]);
  new_muq2[,2] <- (matrix_muq2_reg1[,1] * matrix_muq2_reg2[,3]) + 
  				  (matrix_muq2_reg1[,3] * matrix_muq2_reg2[,1]) +
	 			  (matrix_muq2_reg1[,2] * matrix_muq2_reg2[,2]);
  new_muq2[,3] <- (matrix_muq2_reg1[,2] * matrix_muq2_reg2[,3]) + 
  				  (matrix_muq2_reg1[,3] * matrix_muq2_reg2[,2]) + 
				  (matrix_muq2_reg1[,2] * matrix_muq2_reg2[,2]);
  mes$muq2[[reg1]] <- new_muq2;
  mes$muq2 <- mes$muq2[-reg2];

  # mug: Again, same kind of thing as muq
  new_mug <- array(0,dim=c(1,nbpatients,3));
  new_mug[,,1] <- (mes$mug[reg1,,1] * mes$mug[reg2,,1]) +
  				 (mes$mug[reg1,,1] * mes$mug[reg2,,2]) +
				 (mes$mug[reg1,,2] * mes$mug[reg2,,1]);
  new_mug[,,2] <- (mes$mug[reg1,,1] * mes$mug[reg2,,3]) + 
  				 (mes$mug[reg1,,3] * mes$mug[reg2,,1]) +
	 			 (mes$mug[reg1,,2] * mes$mug[reg2,,2]);
  new_mug[,,3] <- (mes$mug[reg1,,2] * mes$mug[reg2,,3]) + 
  				 (mes$mug[reg1,,3] * mes$mug[reg2,,2]) + 
				 (mes$mug[reg1,,2] * mes$mug[reg2,,2]);
  mes$mug <- abind(mes$mug[0:(reg1-1),,],new_mug,mes$mug[(reg1+1):(nbgenes+1),,],along=1);
  mes$mug <- mes$mug[-reg2,,];

  # mug2: Again, same kind of thing as muq2
  new_mug2 <- array(0,dim=c(1,nbpatients,3));
  new_mug2[,,1] <- (mes$mug2[reg1,,1] * mes$mug2[reg2,,1]) +
  				  (mes$mug2[reg1,,1] * mes$mug2[reg2,,2]) +
				  (mes$mug2[reg1,,2] * mes$mug2[reg2,,1]);
  new_mug2[,,2] <- (mes$mug2[reg1,,1] * mes$mug2[reg2,,3]) + 
  				  (mes$mug2[reg1,,3] * mes$mug2[reg2,,1]) +
	 			  (mes$mug2[reg1,,2] * mes$mug2[reg2,,2]);
  new_mug2[,,3] <- (mes$mug2[reg1,,2] * mes$mug2[reg2,,3]) + 
  				  (mes$mug2[reg1,,3] * mes$mug2[reg2,,2]) + 
				  (mes$mug2[reg1,,2] * mes$mug2[reg2,,2]);
  mes$mug2 <- abind(mes$mug2[0:(reg1-1),,],new_mug2,mes$mug2[(reg1+1):(nbgenes+1),,],along=1);
  mes$mug2 <- mes$mug2[-reg2,,];

  # muh: just need to remove second region, we're not giving our merged one any more prior weight than others
  mes$muh <- mes$muh[-reg2,];

  # muh2: We'll be taking an OR of the beliefs about the relevance of the two genes
  new_muh2 <- array(0,dim=c(1,nbpatients,2));
  # new value is 0 iff previous were both 0
  new_muh2[,,1] <- mes$muh2[reg1,,1] * mes$muh2[reg2,,1];
  # ohterwise it should be 1
  new_muh2[,,2] <- (mes$muh2[reg1,,1] * mes$muh2[reg2,,2]) + 
				  (mes$muh2[reg1,,2] * mes$muh2[reg2,,1]) +
				  (mes$muh2[reg1,,2] * mes$muh2[reg2,,2]);
  mes$muh2 <- abind(mes$muh2[0:(reg1-1),,],new_muh2,mes$muh2[(reg1+1):(nbgenes+1),,],along=1);
  mes$muh2 <- mes$muh2[-reg2,,];

  return(mes);
}


get.new.objective=function(func,args,obj,info) {
  # Given a function func (basically either addvariant_db or mergevariants_db macros), and the arguments given by arg, it will call the function, rerun message passing until convergence, and then evaluate the objective function.
  # obj is our objective function, which must be a function with only two arguments: info and a single float representing the likelihood of the model
  # Info MUST contain all arguments needed for sumproduct, as well as any other information needed to compute the objective function
  # Returns info, with all modified fields updated as well as a new field obj.val containing the resulting value of the objective function
  
  # Unload the contents of info into our current environment
  list2env(info,envir=environment());
  ## print(mes$mux);
  print(args);
  # Call the function to get our modification
  mes <- do.call(func, args);
  #print(nbgenes);
  #print(mes$mux);
  # Rerun message passing until convergence
  mes<- sumproduct(codeDir,nbgenes,nbpatients,nbsnps,harm,harmgene,meancgenes,complexityDistr,pheno,hom,het,mes,net=net,e=e, cores=corse,ratioSignal=ratioSignal,decay=decay,alpha=alpha,netparams=netparams,removeExpressionOnly=removeExpressionOnly,propagate=propagate);

  # Cram all our new stuff back into info so we can pass it into the objective function
  info <- as.list.environment(environment());

  # Calculate result of running objective function and store it
  info$obj.val <- obj(info,sum(info$mes$likelihood));

  # Return our new info, so the updated model and new objective is all available
  return(info);
}


  # Create our macro for testing and setting. This will take in a function (again, basically add.variant.db or merge.variants.db), and arguments, will get the new objective function obj we get from doing the operation and rerunning sumproduct, and setting it to the new result if it's beneficial according to the objective function.
  # Updates the resulting info, the same if it wasn't beneficial, otherwise the updated version
test.and.set <- defmacro(func,args,obj,info,expr={
  obj.old.val <- obj(info,sum(info$mes$likelihood));
  new.info <- get.new.objective(func,args,obj,info);
  if(new.info$obj.val > obj.old.val) {
    info <- new.info;
    print(paste(obj.old.val, " -> ", new.info$obj.val));
  }
});

merge.collapse=function(obj,info) {
  # Starting from the first and second regions, try to collapse each pair of adjacent regions until we get to then end, using obj as the objective function criteria for whether or not we merge.
  # Returns the updated info
  merge.index <- 1;
  while(merge.index < length(info$nbsnps)) {
    prev.length <- length(info$nbsnps);
    test.and.set(mergevariants.db,list(merge.index,merge.index+1),obj,info);
    if(prev.length == length(info$nbsnps)) {
      merge.index <- merge.index + 1;
    } else {
      print(paste("Merging ", merge.index, " and ", merge.index+1));
	}
  }
  return(info);
}

merge.collapse.check.all=function(obj,info) {
  # Works similarly to merge.collapse, but goes through checking all pairs first to see which is best, then does collapse and repeat
  # Returns the updated info
  while(1) {
	  obj.old.val <- obj(info,sum(info$mes$likelihood));
	  best.obj <- obj.old.val;
	  best.ind <- -1;
	  best.info <- NULL;
	  merge.index <- 1;
	  while(merge.index < length(info$nbsnps)) {
		new.info <- get.new.objective(mergevariants.db,list(merge.index,merge.index+1),obj,info);
	    if(new.info$obj.val > best.obj) {
			best.obj <- new.info$obj.val;
			best.ind <- merge.index;
			best.info <- new.info;
		}
		merge.index <- merge.index+1;
	  }
	  if(best.ind > -1) {
		  print(paste("Merging ", best.ind, " and ", best.ind+1));
		  info <- best.info;
	  } else {
		  break;
	  }
  }
  return(info);
}
likelihood.objective=function(info,lik) {
  # The simplest objective function, just the likelihood
  return(lik);
}

sqrt.acet.sqrt.penalty.objective=function(info,lik) {
  # This objective function is likelihood - sqrt(sum(size of regions/2)^2) + sum_regions(sqrt(sum_varinregion(h3k for variant + 1)))
  penalty <- sum(as.numeric(lapply(info$nbsnps, function(x) (x^2/4))));
  #penalty <- sqrt(sum(as.numeric(lapply(info$nbsnps, function(x) (x^2/4)))));
  acet <- sum(as.numeric(lapply(info$h3k,function(x) sqrt(sum(as.numeric(x))+1))));
  return(lik - penalty + acet);
}

sqrt.penalty.objective=function(info,lik) {
  # This objective function is likelihood - sqrt(sum(size of regions/2)^2)
  penalty <- sqrt(sum(as.numeric(lapply(info$nbsnps, function(x) (x^2/4)))));
  return(lik - penalty);
}
