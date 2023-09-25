coding <- function(chr,gene_name,genofile,obj_nullmodel,genes,
                   rare_maf_cutoff=0.01,rv_num_cutoff=2,
                   QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                   Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                   Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,silent=FALSE){

	## evaluate choices
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	phenotype.id <- as.character(obj_nullmodel$id_include)

	## get SNV id, position, REF, ALT (whole genome)
	filter <- seqGetData(genofile, QC_label)
	if(variant_type=="variant")
	{
		SNVlist <- filter == "PASS"
	}

	if(variant_type=="SNV")
	{
		SNVlist <- (filter == "PASS") & isSNV(genofile)
	}

	if(variant_type=="Indel")
	{
		SNVlist <- (filter == "PASS") & (!isSNV(genofile))
	}

	position <- as.numeric(seqGetData(genofile, "position"))
	variant.id <- seqGetData(genofile, "variant.id")

	rm(filter)
	gc()

	### Gene
	kk <- which(genes[,1]==gene_name)

	sub_start_loc <- genes[kk,3]
	sub_end_loc <- genes[kk,4]

	is.in <- (SNVlist)&(position>=sub_start_loc)&(position<=sub_end_loc)
	variant.id.gene <- variant.id[is.in]

	rm(position)
	gc()

	seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

	## Gencode_Exonic
	GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
	## Gencode
	GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
	## Meta.SVM.Pred
	MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))

	################################################
	#           Coding
	################################################
	variant.id.gene <- seqGetData(genofile, "variant.id")
	lof.in.coding <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|(GENCODE.EXONIC.Category=="nonsynonymous SNV")|(GENCODE.EXONIC.Category=="synonymous SNV")
	variant.id.gene <- variant.id.gene[lof.in.coding]

	seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

	## Gencode_Exonic
	GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
	## Gencode
	GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
	## Meta.SVM.Pred
	MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))

	## Annotation
	Anno.Int.PHRED.sub <- NULL
	Anno.Int.PHRED.sub.name <- NULL

	if(variant_type=="SNV")
	{
		if(Use_annotation_weights)
		{
			for(k in 1:length(Annotation_name))
			{
				if(Annotation_name[k]%in%Annotation_name_catalog$name)
				{
					Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
					Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))

					if(Annotation_name[k]=="CADD")
					{
						Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
					}

					if(Annotation_name[k]=="aPC.LocalDiversity")
					{
						Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
						Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
						Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
					}
					Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
				}
			}

			Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
			colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
		}
	}

	################################################
	#                  plof_ds
	################################################
	variant.id.gene <- seqGetData(genofile, "variant.id")
	lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|((GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D"))
	variant.id.gene.category <- variant.id.gene[lof.in.plof]

	seqSetFilter(genofile,variant.id=variant.id.gene.category,sample.id=phenotype.id)

	## genotype id
	id.genotype <- seqGetData(genofile,"sample.id")
	# id.genotype.match <- rep(0,length(id.genotype))

	id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
	phenotype.id.merge <- data.frame(phenotype.id)
	phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
	id.genotype.match <- phenotype.id.merge$index

	## Genotype
	Geno <- seqGetData(genofile, "$dosage")
	Geno <- Geno[id.genotype.match,,drop=FALSE]
	
	# get the rsID 
	rsIDs <- seqGetData(genofile, "annotation/id")  
	
	## impute missing
	if(!is.null(dim(Geno)))
	{
		if(dim(Geno)[2]>0)
		{
			if(geno_missing_imputation=="mean")
			{
				Geno <- matrix_flip_mean(Geno)$Geno
			}
			if(geno_missing_imputation=="minor")
			{
				Geno <- matrix_flip_minor(Geno)$Geno
			}
		}
	}

	## Annotation
	Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.plof,]

	pvalues <- 0
	try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)

	results_plof_ds <- list()
	if(class(pvalues)=="list")
	{
	  # change the original code to output list
		results_temp <- list()
		results_temp$Gene_name <- as.character(genes[kk,1])
		results_temp$Chr <- chr
		results_temp$Category <- "plof_ds"
		results_temp$'#SNV' <- pvalues$num_variant
		
		# add the two kinds of IDs to the results
		results_temp$rsIDs <- rsIDs[pvalues$RV_label]
		results_temp$variantIDs <- variant.id.gene.category[pvalues$RV_label]

		results_temp <- c(results_temp,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
		pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
		pvalues$results_STAAR_A_1_1)
		
		results_temp$'ACAT-O' <- pvalues$results_ACAT_O
		results_temp$'STAAR-O' <- pvalues$results_STAAR_O

		results_plof_ds <- c(results_plof_ds, results_temp)
	}
	

	#####################################################
	#                      plof
	#####################################################
	# variant.id.gene <- seqGetData(genofile, "variant.id")
	lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")
	variant.id.gene.category <- variant.id.gene[lof.in.plof]

	seqSetFilter(genofile,variant.id=variant.id.gene.category,sample.id=phenotype.id)

	## genotype id
	id.genotype <- seqGetData(genofile,"sample.id")
	# id.genotype.match <- rep(0,length(id.genotype))

	id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
	phenotype.id.merge <- data.frame(phenotype.id)
	phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
	id.genotype.match <- phenotype.id.merge$index

	## Genotype
	Geno <- seqGetData(genofile, "$dosage")
	Geno <- Geno[id.genotype.match,,drop=FALSE]

	# get the rsID 
	rsIDs <- seqGetData(genofile, "annotation/id") 
	
	## impute missing
	if(!is.null(dim(Geno)))
	{
		if(dim(Geno)[2]>0)
		{
			if(geno_missing_imputation=="mean")
			{
				Geno <- matrix_flip_mean(Geno)$Geno
			}
			if(geno_missing_imputation=="minor")
			{
				Geno <- matrix_flip_minor(Geno)$Geno
			}
		}
	}

	## Annotation
	Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.plof,]

	pvalues <- 0
	try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)

	results_plof <- list()
	if(class(pvalues)=="list")
	{
		# change the original code to output list
		results_temp <- list()
		results_temp$Gene_name <- as.character(genes[kk,1])
		results_temp$Chr <- chr
		results_temp$Category <- "plof"
		results_temp$'#SNV' <- pvalues$num_variant
		
		# add the two kinds of IDs to the results
		results_temp$rsIDs <- rsIDs[pvalues$RV_label]
		results_temp$variantIDs <- variant.id.gene.category[pvalues$RV_label]
		
		results_temp <- c(results_temp,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
		                  pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
		                  pvalues$results_STAAR_A_1_1)
		
		results_temp$'ACAT-O' <- pvalues$results_ACAT_O
		results_temp$'STAAR-O' <- pvalues$results_STAAR_O
		
		results_plof <- c(results_plof,results_temp)
	}

	#############################################
	#             synonymous
	#############################################
	lof.in.synonymous <- (GENCODE.EXONIC.Category=="synonymous SNV")
	variant.id.gene.category <- variant.id.gene[lof.in.synonymous]

	seqSetFilter(genofile,variant.id=variant.id.gene.category,sample.id=phenotype.id)

	## genotype id
	id.genotype <- seqGetData(genofile,"sample.id")
	# id.genotype.match <- rep(0,length(id.genotype))

	id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
	phenotype.id.merge <- data.frame(phenotype.id)
	phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
	id.genotype.match <- phenotype.id.merge$index

	## Genotype
	Geno <- seqGetData(genofile, "$dosage")
	Geno <- Geno[id.genotype.match,,drop=FALSE]

	# get the rsID 
	rsIDs <- seqGetData(genofile, "annotation/id") 
	
	## impute missing
	if(!is.null(dim(Geno)))
	{
		if(dim(Geno)[2]>0)
		{
			if(geno_missing_imputation=="mean")
			{
				Geno <- matrix_flip_mean(Geno)$Geno
			}
			if(geno_missing_imputation=="minor")
			{
				Geno <- matrix_flip_minor(Geno)$Geno
			}
		}
	}

	## Annotation
	Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.synonymous,]

	pvalues <- 0
	try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)

	results_synonymous <- list()
	if(class(pvalues)=="list")
	{
		# change the original code to output list
		results_temp <- list()
		results_temp$Gene_name <- as.character(genes[kk,1])
		results_temp$Chr <- chr
		results_temp$Category <- "synonymous"
		results_temp$'#SNV' <- pvalues$num_variant
		
		# add the two kinds of IDs to the results
		results_temp$rsIDs <- rsIDs[pvalues$RV_label]
		results_temp$variantIDs <- variant.id.gene.category[pvalues$RV_label]
		
		results_temp <- c(results_temp,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
		                  pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
		                  pvalues$results_STAAR_A_1_1)
		
		results_temp$'ACAT-O' <- pvalues$results_ACAT_O
		results_temp$'STAAR-O' <- pvalues$results_STAAR_O
		
		results_synonymous <- c(results_synonymous,results_temp)
	}


	#################################################
	#        missense
	#################################################
	lof.in.missense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")
	variant.id.gene.category <- variant.id.gene[lof.in.missense]

	seqSetFilter(genofile,variant.id=variant.id.gene.category,sample.id=phenotype.id)

	## genotype id
	id.genotype <- seqGetData(genofile,"sample.id")
	# id.genotype.match <- rep(0,length(id.genotype))

	id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
	phenotype.id.merge <- data.frame(phenotype.id)
	phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
	id.genotype.match <- phenotype.id.merge$index

	## Genotype
	Geno <- seqGetData(genofile, "$dosage")
	Geno <- Geno[id.genotype.match,,drop=FALSE]

	# get the rsID 
	rsIDs <- seqGetData(genofile, "annotation/id") 
	
	## impute missing
	if(!is.null(dim(Geno)))
	{
		if(dim(Geno)[2]>0)
		{
			if(geno_missing_imputation=="mean")
			{
				Geno <- matrix_flip_mean(Geno)$Geno
			}
			if(geno_missing_imputation=="minor")
			{
				Geno <- matrix_flip_minor(Geno)$Geno
			}
		}
	}

	## Annotation
	Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.missense,]

	pvalues <- 0
	try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)

	results <- list()
	if(class(pvalues)=="list")
	{
		# change the original code to output list
		results_temp <- list()
		results_temp$Gene_name <- as.character(genes[kk,1])
		results_temp$Chr <- chr
		results_temp$Category <- "missense"
		results_temp$'#SNV' <- pvalues$num_variant
		
		# add the two kinds of IDs to the results
		results_temp$rsIDs <- rsIDs[pvalues$RV_label]
		results_temp$variantIDs <- variant.id.gene.category[pvalues$RV_label]
		
		results_temp <- c(results_temp,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
		                  pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
		                  pvalues$results_STAAR_A_1_1)
		
		results_temp$'ACAT-O' <- pvalues$results_ACAT_O
		results_temp$'STAAR-O' <- pvalues$results_STAAR_O
		
		results <- c(results,list(results_temp))
	}

	#################################################
	#         disruptive missense
	#################################################
	lof.in.dmissense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D")
	variant.id.gene.category <- variant.id.gene[lof.in.dmissense]

	seqSetFilter(genofile,variant.id=variant.id.gene.category,sample.id=phenotype.id)

	## genotype id
	id.genotype <- seqGetData(genofile,"sample.id")
	# id.genotype.match <- rep(0,length(id.genotype))

	id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
	phenotype.id.merge <- data.frame(phenotype.id)
	phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
	id.genotype.match <- phenotype.id.merge$index

	## Genotype
	Geno <- seqGetData(genofile, "$dosage")
	Geno <- Geno[id.genotype.match,,drop=FALSE]

	# get the rsID 
	rsIDs <- seqGetData(genofile, "annotation/id") 
	
	## impute missing
	if(!is.null(dim(Geno)))
	{
		if(dim(Geno)[2]>0)
		{
			if(geno_missing_imputation=="mean")
			{
				Geno <- matrix_flip_mean(Geno)$Geno
			}
			if(geno_missing_imputation=="minor")
			{
				Geno <- matrix_flip_minor(Geno)$Geno
			}
		}
	}

	## Annotation
	Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.dmissense,]

	pvalues <- 0
	try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)

	if(class(pvalues)=="list")
	{
		# change the original code to output list
		results_temp <- list()
		results_temp$Gene_name <- as.character(genes[kk,1])
		results_temp$Chr <- chr
		results_temp$Category <- "disruptive_missense"
		results_temp$'#SNV' <- pvalues$num_variant
		
		# add the two kinds of IDs to the results
		results_temp$rsIDs <- rsIDs[pvalues$RV_label]
		results_temp$variantIDs <- variant.id.gene.category[pvalues$RV_label]
		
		results_temp <- c(results_temp,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
		                  pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
		                  pvalues$results_STAAR_A_1_1)
		
		results_temp$'ACAT-O' <- pvalues$results_ACAT_O
		results_temp$'STAAR-O' <- pvalues$results_STAAR_O
		
		results <- c(results,list(results_temp))		
		
	}

	if(length(results)!=0)
	{
	  
		if(length(results)==1)
		{
			if(results[[1]]$Category!="disruptive_missense")
			{
			  # modify the original code to adapt to list
			  new_elements <- list("SKAT(1,25)-Disruptive" = 1,
			                       "SKAT(1,1)-Disruptive" = 1,
			                       "Burden(1,25)-Disruptive" = 1,
			                       "Burden(1,1)-Disruptive" = 1,
			                       "ACAT-V(1,25)-Disruptive" = 1,
			                       "ACAT-V(1,1)-Disruptive" = 1)
			  
			  # Append the new elements to the results list
			  results <- c(results[[1]], new_elements)
				results_missense <- results
				results_ds <- list()
			}else
			{
				results_missense <- list()
				results_ds <- results[[1]]
				results <- list()
			}
		}

		if(length(results)!=0)
		{
			if(length(results)==2)
			{
			  # modify the original code to adapt to list
			  new_elements <- list("SKAT(1,25)-Disruptive" = results[[2]][["SKAT(1,25)"]],
			                       "SKAT(1,1)-Disruptive" = results[[2]][["SKAT(1,1)"]],
			                       "Burden(1,25)-Disruptive" = results[[2]][["Burden(1,25)"]],
			                       "Burden(1,1)-Disruptive" = results[[2]][["Burden(1,1)"]],
			                       "ACAT-V(1,25)-Disruptive" = results[[2]][["ACAT-V(1,25)"]],
			                       "ACAT-V(1,1)-Disruptive" = results[[2]][["ACAT-V(1,1)"]])
			  results_m <- c(results[[1]], new_elements)
			  
				apc_num <- (length(results_m)-20)/6
				p_seq <- c(1:apc_num,1:apc_num+(apc_num+1),1:apc_num+2*(apc_num+1),1:apc_num+3*(apc_num+1),1:apc_num+4*(apc_num+1),1:apc_num+5*(apc_num+1),(6*apc_num+9):(6*apc_num+14))
				results_m[["STAAR-O"]] <- CCT(as.numeric(results_m[7:length(results_m)][p_seq]))
				results_m[["STAAR-S(1,25)"]] <- CCT(as.numeric(results_m[7:length(results_m)][c(1:apc_num,6*apc_num+9)]))
				results_m[["STAAR-S(1,1)"]] <- CCT(as.numeric(results_m[7:length(results_m)][c(1:apc_num+(apc_num+1),6*apc_num+10)]))
				results_m[["STAAR-B(1,25)"]] <- CCT(as.numeric(results_m[7:length(results_m)][c(1:apc_num+2*(apc_num+1),6*apc_num+11)]))
				results_m[["STAAR-B(1,1)"]] <- CCT(as.numeric(results_m[7:length(results_m)][c(1:apc_num+3*(apc_num+1),6*apc_num+12)]))
				results_m[["STAAR-A(1,25)"]] <- CCT(as.numeric(results_m[7:length(results_m)][c(1:apc_num+4*(apc_num+1),6*apc_num+13)]))
				results_m[["STAAR-A(1,1)"]] <- CCT(as.numeric(results_m[7:length(results_m)][c(1:apc_num+5*(apc_num+1),6*apc_num+14)]))

				results_ds <- list()
				results_ds <- c(results_ds,results[[2]])

				results <- list()
				results <- c(results,results_m)
			}
		}
	}else
	{
		results <- list()
		results_ds <- list()
	}

	results_coding <- list(plof=results_plof,plof_ds=results_plof_ds,missense=results,disruptive_missense=results_ds,synonymous=results_synonymous)

	seqResetFilter(genofile)

	return(results_coding)
}

