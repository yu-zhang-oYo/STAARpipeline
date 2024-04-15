# setwd("/Users/yzh10/Library/CloudStorage/OneDrive-IndianaUniversity/research/nuMoM2b/whole_sequencing/GA")
# Check if 'igraph' package is installed; install if it's not, then load it
if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}
library(igraph)

library(MendelianRandomization)

# 
# Step1: check STAAR-O/STAAR-B. 
  # 1. ?check P values of burden in each annotation 
  # 2. ?for one gene, maybe a kind of annotation is used as instrument while the other one is not.
# Step2:
#   
#   1. comon1 + common2  + gene1_miss + gene1_plof

# question: 
# 1. how to choose suitable genes and variants?
# 2. considering the independence?


# calculate the correlation of two rare gene set to check if they are indpendent
# Weighted_gene_by_gene is a list of Geno_Rare_Weighted of each gene
gene_cor <- function(Weighted_gene_by_gene){
  # Initialize a list to hold combined data for each annotation
  Weighted_gene_by_annotation <- list()
  # Determine the number of annotations (assuming the number is consistent across all genes)
  num_annotations <- ncol(Weighted_gene_by_gene[[1]])
  annotation_names <- colnames(Weighted_gene_by_gene[[1]])
  # Loop over each annotation
  for (i in 1:num_annotations) {
    Weighted_gene_by_annotation[[i]] <- do.call(cbind, lapply(Weighted_gene_by_gene, function(x) x[, i]))
    # Set the name of each list element to correspond to the annotation
    names(Weighted_gene_by_annotation)[i] <- annotation_names[i]
  }
  # Calculate the correlation matrix between genes for each annotation
  correlation_between_genes <- lapply(Weighted_gene_by_annotation, cor, method = "pearson")
  return(correlation_between_genes)
}

# input the orginal results
select_LE_sig_gene <- function(results, le_threshold = 0.3, sig_threshold=5e-2) {
  
  # Subset the list to keep only non-empty elements
  results <- Filter(length, results)
  names(results) <- paste0(sapply(results, function(x){x$Gene_name}), "_", names(results))
  # get the Geno_Rare_Weighted variables
  Weighted_gene_by_gene <- lapply(results, function(x) x$Geno_Rare_Weighted)
  gene_cors <- gene_cor(Weighted_gene_by_gene)
  num_gene <- length(Weighted_gene_by_gene)
  
  annotation_names <- colnames(Weighted_gene_by_gene[[1]])
  
  # get the P values of burden so that we can select the minimum P values of a group of LD genes
  pvalues_burden <- lapply(results,function(x) x[annotation_names])
  # re-arrange the P values so that it is P values of different gene in each annotation
  pvalues_burden_by_annotation <- list()
  for (i in 1:length(annotation_names)) {
    pvalues_burden_by_annotation[[i]] <- lapply(pvalues_burden, function(x) x[[annotation_names[i]]])
    names(pvalues_burden_by_annotation)[i] <- annotation_names[i]
  }
  
  # Initialize list to store genes with minimum p-value under linkage equilibrium for each annotation
  selected_genes_by_annotation <- list()
  
  for (ann in 1:length(annotation_names)) {
    # Get the correlation matrix for the current annotation
    cor_matrix <- gene_cors[[ann]]
    # To replace row and column names with numeric indices, so that we can save the gene position as the results
    rownames(cor_matrix) <- seq_len(num_gene)
    colnames(cor_matrix) <- seq_len(num_gene)
    # Get the p-values for the current annotation
    pvalues <- unlist(lapply(pvalues_burden, `[[`, annotation_names[ann]))
    
    # Convert the correlation matrix to an adjacency matrix for LD
    adj_matrix <- cor_matrix > le_threshold
    diag(adj_matrix) <- FALSE  # Remove self-loops

    # Create a graph from the adjacency matrix
    g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
    
    # Find connected components of the graph, each represents a group of genes in LD
    # find the the index of the minimum P values in each groups
    components <- components(g)
    min_p_gene_index <- NULL
    # Find the gene with the smallest p-value in each component
    for (comp_idx in seq_len(components$no)) {
      genes_in_component <- V(g)[components$membership == comp_idx]
      gene_indices <- as.numeric(names(genes_in_component))
      min_p_gene_index <- c(min_p_gene_index, gene_indices[which.min(pvalues[gene_indices])])
    }
    # only keep the significant P values
    # ? do we need to select the gene according to the same values, or different annotation
    min_p_sig_gene_index <- min_p_gene_index[pvalues[min_p_gene_index] < sig_threshold]
    selected_genes_by_annotation[[annotation_names[ann]]] <- min_p_sig_gene_index  # this is the index of the selected genes
  }
  
  # select the Burden_Effect_Size of the selected genes for each annotation
  Burden_Effect_Size_by_gene <- lapply(results, function(x) x$Burden_Effect_Size)
  # Initialize a list to hold combined data for each annotation
  Burden_Effect_Size_by_annotation <- list()
  # Determine the number of annotations (assuming the number is consistent across all genes)
  num_annotations <- nrow(Burden_Effect_Size_by_gene[[1]])
  annotation_names <- rownames(Burden_Effect_Size_by_gene[[1]])
  # Loop over each annotation
  for (ith_anno in 1:num_annotations) {
    Burden_Effect_Size_by_annotation[[ith_anno]] <- do.call(rbind, lapply(Burden_Effect_Size_by_gene, function(x) x[ith_anno, ]))
    Burden_Effect_Size_by_annotation[[ith_anno]] <- Burden_Effect_Size_by_annotation[[ith_anno]][selected_genes_by_annotation[[annotation_names[ith_anno]]], ]
    
    # Set the name of each list element to correspond to the annotation
    names(Burden_Effect_Size_by_annotation)[ith_anno] <- annotation_names[ith_anno]
  }
  
  return(Burden_Effect_Size_by_annotation)
}



# example
results <- results_coding
Burden_Effect_Size_by_annotation <- select_LE_sig_gene(results, le_threshold = 0.3, sig_threshold=1e-2)

# input of rare_X is the result for each annotation
MR_gene <- function(common_X, rare_X, common_Y, rare_Y){
  
  num_annotation <- length(rare_X)  # the number of annotation
  MR_estimates <- matrix(nrow = num_annotation, ncol = 3)  # Prepare a matrix to store estimates, standard errors, and p-values
  
  for (ith_anno in 1:num_annotation) {
    
    # Create an MRInput object
    MR_object <- mr_input(
      bx = c(common_X$Est, rare_X[[ith_anno]][, "Burden_Est"]),
      bxse = c(common_X$Est_se, rare_X[[ith_anno]][, "Burden_SE_Est"]),
      by = c(common_Y$Est, rare_Y[[ith_anno]][, "Burden_Est"]),
      byse = c(common_Y$Est_se, rare_Y[[ith_anno]][, "Burden_SE_Est"])
    )
    
    # Perform MR analysis using IVW
    MR_result <- mr_ivw(MR_object)
    MR_estimates[ith_anno,] <- c(MR_result$Estimate, MR_result$StdError, MR_result$Pvalue)
  }
  colnames(MR_estimates) <- c("Estimate", "StdError", "Pvalue")
  rownames(MR_estimates) <- names(rare_X)
  MR_O <- CCT(MR_estimates[, "Pvalue"])
  return(list(MR_anno=MR_estimates, MR_O=MR_O))
}

# # common is the estimates and standard error of estimates of common variants,
# # from the output of individual_analysis.
# # input of rare_X is the result for each genes
# MR_gene <- function(common_X, rare_X, common_Y, rare_Y){
# 
#   num_annotation <- nrow(rare_X[[1]])  # the number of annotation
#   MR_estimates <- matrix(NA, nrow = num_annotation, ncol = 3)  # Prepare a matrix to store estimates, standard errors, and p-values
# 
#   for (anno in 1:num_annotation) {
# 
#     # Aggregate the data for the current annotation across all genes
#     bx = unlist(lapply(rare_X, function(x) x[anno, "Burden_Est"]))
#     bxse = unlist(lapply(rare_X, function(x) x[anno, "Burden_SE_Est"]))
#     by = unlist(lapply(rare_Y, function(y) y[anno, "Burden_Est"]))
#     byse = unlist(lapply(rare_Y, function(y) y[anno, "Burden_SE_Est"]))
# 
#     # Create an MRInput object
#     MR_object <- mr_input(
#       bx = c(common_X$Est, bx),
#       bxse = c(common_X$Est_se, bxse),
#       by = c(common_Y$Est, by),
#       byse = c(common_Y$Est_se, byse)
#     )
# 
#     # Perform MR analysis using IVW
#     MR_result <- mr_ivw(MR_object)
#     MR_estimates[anno,] <- c(MR_result$Estimate, MR_result$StdError, MR_result$Pvalue)
#   }
#   colnames(MR_estimates) <- c("Estimate", "StdError", "Pvalue")
#   rownames(MR_estimates) <- rownames(rare_X[[1]])
#   MR_O <- CCT(MR_estimates[, "Pvalue"])
#   return(list(MR_anno=MR_estimates, MR_O=MR_O))
# }

################## get the results for all the significant genes ##################
# when we only have the outcome individual data
genes_name <- c("ABCA1", "LCAT", "ABCA1", "SCARB1", "APOC3", "CD36", "LPL", "LIPC", "CETP", "CD300LG")
# genes on chromosome 19
genes_name <- c("LDLR", "APOE", "TM6SF2")

results_coding1 <- c()
for(gene_name in genes_name)
{
  print(gene_name)
  results <- Gene_Centric_Coding(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,rare_maf_cutoff=0.05,rv_num_cutoff=2,
                                 QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                 Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                 Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  results_coding1 <- append(results_coding1,results)
}
# save(results_coding,file=paste0(output_path,output_file_name,"_all_",chr,".Rdata"))
# seqClose(genofile)


################## only select the significant gene category on each genes for TC##################
results <- results_coding1
gene_category <- c("LDLR_plof_ds", "LDLR_disruptive_missense", "LDLR_missense", "APOE_plof_ds", "APOE_disruptive_missense", "APOE_missense")
burden11 <- data.frame(Burden_Est=c(12.819, 12.081, 7.112, -7.674, -7.555, -5.927), Burden_pvalue=c(7.67E-30, 5.52E-26, 7.19E-18, 1.18E-08, 2.61E-08, 6.25E-08))
Z_Score <- qnorm(burden11$Burden_pvalue / 2, lower.tail = FALSE)
# Calculate Standard Error
burden11$Burden_SE_Est <- abs(burden11$Burden_Est / Z_Score)
rownames(burden11) <- gene_category
# the input is the buden effect of each gene for each annotation
# NOTE: revise the coding.R so that we can calculate the burden effect for certain gene_category, and annotation
Burden_Effect_by_annotation_X <- list("Burden(1,1)" = burden11)
pvalues_X <- list("Burden(1,1)" = c(7.67E-30, 5.52E-26, 7.19E-18, 1.18E-08, 2.61E-08, 6.25E-08))
annotation_names_use <- c("Burden(1,1)")
annotation_names_use <- annotation_names 


# input the orginal results
select_LE_sig_gene <- function(results, gene_category, annotation_names_use, Burden_Effect_Size_by_annotation_X, le_threshold = 0.3, sig_threshold=5e-2) {
  
  # Subset the list to keep only non-empty elements
  results <- Filter(length, results)
  # ? each annotation has diffferent significant genes, how to consider this into the code
  names(results) <- paste0(sapply(results, function(x){x$Gene_name}), "_", names(results))
  results <- results[gene_category]
  # get the Geno_Rare_Weighted variables
  Weighted_gene_by_gene <- lapply(results, function(x) x$Geno_Rare_Weighted)
  gene_cors <- gene_cor(Weighted_gene_by_gene)
  num_gene <- length(Weighted_gene_by_gene)
  
  annotation_names <- colnames(Weighted_gene_by_gene[[1]])
  
  # get the P values of burden so that we can select the minimum P values of a group of LD genes
  pvalues_burden <- lapply(results,function(x) x[annotation_names])
  # re-arrange the P values so that it is P values of different gene in each annotation
  pvalues_burden_by_annotation <- list()
  for (i in 1:length(annotation_names)) {
    pvalues_burden_by_annotation[[i]] <- lapply(pvalues_burden, function(x) x[[annotation_names[i]]])
    names(pvalues_burden_by_annotation)[i] <- annotation_names[i]
  }
  
  # Initialize list to store genes with minimum p-value under linkage equilibrium for each annotation
  selected_genes_by_annotation_use <- list()
  
  for (ann in 1:length(annotation_names_use)) {
    # Get the correlation matrix for the current annotation
    cor_matrix <- gene_cors[[ann]]
    # To replace row and column names with numeric indices, so that we can save the gene position as the results
    rownames(cor_matrix) <- seq_len(num_gene)
    colnames(cor_matrix) <- seq_len(num_gene)

    # Convert the correlation matrix to an adjacency matrix for LD
    adj_matrix <- abs(cor_matrix) > le_threshold
    diag(adj_matrix) <- FALSE  # Remove self-loops
    
    # Create a graph from the adjacency matrix
    g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
    # plot(g, main = "Graph of Genes in LD")
    
    # Find a maximal independent set
    # independent_sets is a list, randomly select one if don't have Pvalues, or select the set with minimum P values
    independent_sets <- largest_ivs(g)
    if(is.null(pvalues_X)){
      selected_genes_by_annotation_use[[annotation_names_use[ann]]]  <- gene_category[unlist(sample(independent_sets, 1))]
    } else {
      pvalues_X_each_anno <- pvalues_X[[annotation_names_use[ann]]]
      # Using lapply to calculate the sum of P-values in each set
      sums_pvalues <- lapply(independent_sets, function(x) {
        sum(pvalues_X_each_anno[x])  # Sum the P-values at these indices
      })
      # Find the index of the set with the minimum sum of P-values
      min_index <- which.min(unlist(sums_pvalues))
      selected_genes_by_annotation_use[[annotation_names_use[ann]]] <- gene_category[independent_sets[[min_index]]]
    }
  }
   
  # select the Burden_Effect_Size of the selected genes for each annotation
  Burden_Effect_Size_by_gene <- lapply(results, function(x) x$Burden_Effect_Size)
  # Initialize a list to hold combined data for each annotation
  Burden_Effect_Size_by_annotation <- list()
  # Determine the number of annotations (assuming the number is consistent across all genes)
  num_annotations <- nrow(Burden_Effect_Size_by_gene[[1]])
  annotation_names <- rownames(Burden_Effect_Size_by_gene[[1]])
  # Loop over each annotation
  for (ith_anno in annotation_names_use) {
    Burden_Effect_Size_by_annotation[[ith_anno]] <- do.call(rbind, lapply(Burden_Effect_Size_by_gene, function(x) x[ith_anno, ]))
    Burden_Effect_Size_by_annotation[[ith_anno]] <- Burden_Effect_Size_by_annotation[[ith_anno]][selected_genes_by_annotation_use[[ith_anno]], ]
  }
  return(selected_genes_by_annotation_use = selected_genes_by_annotation_use, Burden_Effect_Size_by_annotation = Burden_Effect_Size_by_annotation)
}





