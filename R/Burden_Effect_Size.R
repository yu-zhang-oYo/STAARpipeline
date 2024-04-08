# refer to the effect size calculation in metaSTAAR https://github.com/xihaoli/MetaSTAAR/blob/main/R/Burden_Burden_Effect_Size_meta.R
# Note: only check the function BurdenEffectSize_O_SMMAT because we have the data for this, still need to check the other two functions.


# the effect size under 3 situations
BurdenEffectSize_O <- function(G, X, working, sigma, fam, residuals, weights_B) {
  if (fam == 0) {
    tX_G <- t(X) %*% G
    Cov <- t(G) %*% G - t(tX_G) %*% solve(t(X) %*% X) %*% tX_G
  } else {
    tX_G <- t(X) %*% diag(working) %*% G
    Cov <- t(diag(working) %*% G) %*% G - t(tX_G) %*% solve(t(X) %*% diag(working) %*% X) %*% tX_G
  }
  
  x <- t(residuals) %*% G
  
  # Initialize matrix to store results for each annotation
  burden_Results <- matrix(0, nrow = ncol(weights_B), ncol = 5)
  # colnames(burden_Results) <- c("Burden_Est", "Burden_SE_Est", "Burden_Score_Stat", "Burden_SE_Score", "Burden_pvalue")
  
  for (i in 1:ncol(weights_B)) {
    weight_vector <- weights_B[, i]
    sum0 <- sum(x * weight_vector)
    sumw <- sum((weight_vector %*% Cov) * weight_vector)
    
    # Storing the burden estimate and its SE in the respective columns of the matrix
    burden_Results[i, 1] <- sum0 / sumw  # Burden effect size estimate
    burden_Results[i, 2] <- 1 / sqrt(sumw)  # Standard error estimate
    burden_Results[i, 3] <- sum0  # Burden score sta
    burden_Results[i, 4] <- sqrt(sumw)  # Standard error of burden score
    burden_Results[i, 5] <- pchisq(sum0^2/sumw, 1, lower.tail=FALSE)  # P value of burden score
  }
  
  return(burden_Results)
}

BurdenEffectSize_O_SMMAT <- function(G, P, residuals, weights_B) {
  Cov <- t(as.matrix(P %*% G)) %*% G  # P is not a matrix
  # Cov <- t(P %*% G) %*% G
  
  x <- t(residuals) %*% G
  
  # Initialize matrix to store results for each annotation
  burden_Results <- matrix(0, nrow = ncol(weights_B), ncol = 5)
  # colnames(burden_Results) <- c("Burden_Est", "Burden_SE_Est", "Burden_Score_Stat", "Burden_SE_Score", "Burden_pvalue")
  
  for (i in 1:ncol(weights_B)) {
    weight_vector <- weights_B[, i]
    sum0 <- sum(x * weight_vector)
    sumw <- sum((weight_vector %*% Cov) * weight_vector)
    
    # Storing the burden estimate and its SE in the respective columns of the matrix
    burden_Results[i, 1] <- sum0 / sumw  # Burden effect size estimate
    burden_Results[i, 2] <- 1 / sqrt(sumw)  # Standard error estimate
    burden_Results[i, 3] <- sum0  # Burden score sta
    burden_Results[i, 4] <- sqrt(sumw)  # Standard error of burden score
    burden_Results[i, 5] <- pchisq(sum0^2/sumw, 1, lower.tail=FALSE)  # P value of burden score
  }
  
  return(burden_Results)
}


BurdenEffectSize_O_SMMAT_sparse <- function(G, Sigma_i, Sigma_iX, cov, residuals, weights_B) {
  # Calculate the covariance matrix
  tSigma_iX_G <- t(Sigma_iX) %*% G
  Cov <- t(Sigma_i %*% G) %*% G - tSigma_iX_G %*% cov %*% tSigma_iX_G
  
  x <- t(residuals) %*% G
  
  # Initialize matrix to store burden effect size and SE for each annotation
  burden_Results <- matrix(0, nrow = ncol(weights_B), ncol = 5)
  # colnames(burden_Results) <- c("Burden_Est", "Burden_SE_Est", "Burden_Score_Stat", "Burden_SE_Score", "Burden_pvalue")
  
  for (i in 1:ncol(weights_B)) {
    weight_vector <- weights_B[, i]
    sum0 <- sum(x * weight_vector)
    sumw <- sum((weight_vector %*% Cov) * weight_vector)
    
    # Calculate and store burden estimate and its SE
    burden_Results[i, 1] <- sum0 / sumw  # Burden effect size estimate
    burden_Results[i, 2] <- 1 / sqrt(sumw)  # Standard error estimate
    burden_Results[i, 3] <- sum0  # Burden score sta
    burden_Results[i, 4] <- sqrt(sumw)  # Standard error of burden score
    burden_Results[i, 5] <- pchisq(sum0^2/sumw, 1, lower.tail=FALSE)  # P value of burden score
  }
  
  return(burden_Results)
}


# final calculation of effect size
# or we add this into the STAAR.R code so that it can output all the burden effect size with other results under different settings, 
# we need to revise the code in STAAR_pipeline to output the effect sizes
Burden_EffectSize <- function(genotype, obj_nullmodel, annotation_phred = NULL, rare_maf_cutoff = 0.01, rv_num_cutoff = 2) {
  if(!inherits(genotype, "matrix") && !inherits(genotype, "Matrix")){
    stop("genotype is not a matrix!")
  }

  if(dim(genotype)[2] == 1){
    stop("Number of rare variant in the set is less than 2!")
  }

  annotation_phred <- as.data.frame(annotation_phred)
  if(dim(annotation_phred)[1] != 0 && dim(genotype)[2] != dim(annotation_phred)[1]){
    stop("Dimensions don't match for genotype and annotation!")
  }

  if(inherits(genotype, "sparseMatrix")){
    genotype <- as.matrix(genotype)
  }

  genotype <- matrix_flip(genotype)
  MAF <- genotype$MAF
  RV_label <- as.vector((MAF<rare_maf_cutoff)&(MAF>0))
  Geno_rare <- genotype$Geno[, RV_label]

  rm(genotype)
  gc()
  annotation_phred <- annotation_phred[RV_label, , drop = FALSE]

  if(sum(RV_label) >= rv_num_cutoff){
    G <- as(Geno_rare, "dgCMatrix")
    MAF <- MAF[RV_label]
    rm(Geno_rare)
    gc()

    ## Calculate weights for burden test
    annotation_rank <- 1 - 10^(-annotation_phred / 10)
    ## beta(1,25)
    w_1 <- dbeta(MAF, 1, 25)
    ## beta(1,1)
    w_2 <- dbeta(MAF, 1, 1)

    ## Burden weights
    if(dim(annotation_phred)[2] == 0){
      w_B <- as.matrix(cbind(w_1, w_2))
    } else {
      w_B <- as.matrix(cbind(w_1, annotation_rank * w_1, w_2, annotation_rank * w_2))
    }

    if(obj_nullmodel$relatedness){
      if(!obj_nullmodel$sparse_kins){
        P <- obj_nullmodel$P
        P_scalar <- sqrt(dim(P)[1])
        P <- P * P_scalar

        residuals.phenotype <- obj_nullmodel$scaled.residuals * sqrt(P_scalar)

        Burden_Effect_Size <- BurdenEffectSize_O_SMMAT(G, P, residuals.phenotype, weights_B=w_B)
        
      } else {
        Sigma_i <- obj_nullmodel$Sigma_i
        Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
        cov <- obj_nullmodel$cov

        residuals.phenotype <- obj_nullmodel$scaled.residuals

        Burden_Effect_Size <- BurdenEffectSize_O_SMMAT_sparse(G, Sigma_i, Sigma_iX, cov, residuals.phenotype, weights_B=w_B)
      }
    } else {
      X <- model.matrix(obj_nullmodel)
      working <- obj_nullmodel$weights
      sigma <- sqrt(summary(obj_nullmodel)$dispersion)
      fam <- if(obj_nullmodel$family[1] == "binomial") 1 else 0

      residuals.phenotype <- obj_nullmodel$y - obj_nullmodel$fitted.values

      Burden_Effect_Size <- BurdenEffectSize_O(G, X, working, sigma, fam, residuals.phenotype, weights_B=w_B)
    }
    
    colnames(Burden_Effect_Size) <- c("Burden_Est", "Burden_SE_Est", "Burden_Score_Stat", "Burden_SE_Score", "Burden_pvalue")
    

    if(dim(annotation_phred)[2] == 0){
      rownames(Burden_Effect_Size) <- c("Burden(1,25)","STAAR-B(1,25)","Burden(1,1)","STAAR-B(1,1)")

    }else{
      rownames(Burden_Effect_Size) <- c("Burden(1,25)",
                                          paste0("Burden(1,25)-",colnames(annotation_phred)),
                                         "Burden(1,1)",
                                          paste0("Burden(1,1)-",colnames(annotation_phred))
                                          )
    }

    return(list(
      num_variant = sum(RV_label),
      cMAC = sum(G),
      RV_label = RV_label,
      Burden_Effect_Size = Burden_Effect_Size
    ))
  } else {
    stop(paste0("Number of rare variant in the set is less than ", rv_num_cutoff, "!"))
  }
}

# # test by running the function
# effectsize <- Burden_EffectSize(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff)


