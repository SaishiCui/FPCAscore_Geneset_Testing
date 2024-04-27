if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("slingshot")



rm(list=ls())

library(slingshot)
library(mclust, quietly = TRUE)
library(RColorBrewer)
library(grDevices)
library(fda)
library(Hotelling)
library(FNN)
library(car)
library(compiler)


### Functions 



simulation_data <- function(seed, sigma1, sigma2){
  
  
  set.seed(seed)  
  ## non-DE genes
  non_DE_matrix <- matrix(rep(rnorm(100, 1.25, sigma1) , each  = 1000), nrow = 100, byrow = T)
  
  
  # deactivation
  deactivation_matrix = matrix(NA, nrow = 100, ncol = 1000)
  vector <- exp(atan((c(1000:1)-650)/100))+1 
  for (i in 1:1000) {
    deactivation_matrix[,i] = rnorm(100, vector[i], sigma2)
  }
  
  
  
  
  # activation
  activation_matrix = matrix(NA, nrow = 100, ncol = 1000)
  vector <- exp(atan((c(1:1000)-650)/100))+1 
  for (i in 1:1000) {
    activation_matrix[,i] = rnorm(100, vector[i], sigma2)
  }
  
  
  
  ## cross up
  
  cross_up_matrix = matrix(NA, nrow = 100, ncol = 1000)
  vector <- exp(atan((c(1:1000)-500)/200))+1 
  for (i in 1:1000) {
    cross_up_matrix[,i] = rnorm(100, vector[i], sigma2)
  }
  
  
  ## cross cos
  
  cross_down_matrix = matrix(NA, nrow = 100, ncol = 1000)
  vector <-  exp(atan((c(1000:1)-500)/200))+1 
  for (i in 1:1000) {
    cross_down_matrix[,i] = rnorm(100, vector[i], sigma2)
  }
  
  
  ## parallel up
  
  parallel_up_matrix = matrix(NA, nrow = 100, ncol = 1000)
  vector <- exp(atan((c(1:1000)-500)/300))+2
  for (i in 1:1000) {
    parallel_up_matrix[,i] = rnorm(100, vector[i], sigma2)
  }
  
  
  ## parallel down
  
  parallel_down_matrix = matrix(NA, nrow = 100, ncol = 1000)
  vector <-  exp(atan( (c(1:1000)-500)/300))+1 
  for (i in 1:1000) {
    parallel_down_matrix[,i] = rnorm(100, vector[i], sigma2)
  }
  
  
  
  
  
  means <- rbind(
    non_DE_matrix,
    deactivation_matrix,
    activation_matrix,
    cross_up_matrix,
    cross_down_matrix,
    parallel_up_matrix,
    parallel_down_matrix
  )
  
  means[means < 0] = 0
  
  counts <- apply(means,2,function(cell_means){
    total <- rnbinom(1, mu = 7000, size = 4)
    rmultinom(1, total, cell_means)
  })
  
  rownames(counts) <- c(paste0('G1_',1:100), paste0('G2_',1:100), paste0('G3_',1:100),
                        paste0('G4_',1:100), paste0('G5_',1:100), paste0('G6_',1:100), 
                        paste0('G7_',1:100))
  colnames(counts) <- paste0('c',1:1000)
  sce <- SingleCellExperiment(assays = List(counts = counts))
  return(sce)
} ## Simulation setting function

FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}  ## Normalization function

pseudotime_estimate <- function(data){
  
  par(mfrow = c(2,2))
  pca <- prcomp(t(log1p(data)), scale. = FALSE)
  rd1 <- pca$x[,1:2]
  plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)
  reducedDims(sce) <- SimpleList(PCA = rd1)
  
  
  cl1 <- Mclust(rd1)$classification
  colData(sce)$GMM <- cl1
  plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)
  
  
  sce <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'PCA')
  
  plot(reducedDims(sce)$PCA, col = brewer.pal(9,'Set1')[sce$GMM], pch=16, asp = 1)
  lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')
  
  colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
  plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
  plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
  lines(SlingshotDataSet(sce), lwd=2, col='black')
  pseudotime <- order(sce@colData$slingPseudotime_1, decreasing = F)
  par(mfrow = c(1,1))
  outlist = list(pseudotime = pseudotime, sce= sce)
  return(outlist)
} ## Slingshot (pseudotime inference) function

gaussian_kernel <- function(x, h) {
  return(1 / (sqrt(2 * pi) * h) * exp(-x^2 / (2 * h^2)))
} ## gaussian kernel function

nw_estimator_weights <- function(x, X, Y, h) {
  weights <- sapply(X, function(xi) gaussian_kernel(x - xi, h))
  weights_normalized <- weights / sum(weights)
  estimate <- sum(weights_normalized * Y) / sum(weights_normalized)
  return(estimate)
}  ## NW smoothing function

function_estimate <- function(realization, largest_nbasis = 15, largest_norder = 10, h = 100){
  
  
  
  
  estimate_mat = matrix(0, nrow = nrow(realization), ncol = ncol(realization))
  estimate = c()
  for (line in 1:ncol(realization)) {
    for (x_eval in 1:nrow(realization)){
      estimate[x_eval] <- nw_estimator_weights(x_eval, X = c(1:nrow(realization)), Y = realization[,line], h = h)
    }
    estimate_mat[,line] = estimate
  }
  
  realization_new = estimate_mat
  
  
  npseudo <- dim(realization_new)[1]
  argvals = seq(1, npseudo, len= npseudo)
  gcv_list = c()
  norder_list = c()
  nbasis_list = c()
  lambda_list  = c()
  for (norder in 4:largest_norder) {
    for (nbasis in norder:largest_nbasis) {
      for (lambda in seq(0, 0.2, 0.02)){
        basisobj = create.bspline.basis(c(1, npseudo), nbasis, norder = norder)
        Par = fdPar(fdobj = basisobj, lambda = lambda)
        ys = smooth.basis(argvals, realization_new, Par)
        xfd = ys$fd
        gcv = mean(ys$gcv)
        gcv_list = append(gcv_list, gcv)
        norder_list = append(norder_list, norder)
        nbasis_list = append(nbasis_list, nbasis)
        lambda_list = append(lambda_list, lambda)
        
      }
    }    
  }  
  
  final_norder <- norder_list[which.min(gcv_list)]
  final_nbasis <- nbasis_list[which.min(gcv_list)]
  final_lambda <- lambda_list[which.min(gcv_list)]
  
  
  basisobj_final = create.bspline.basis(c(1,npseudo), final_nbasis, norder = final_norder)
  Par_final  = fdPar(fdobj = basisobj_final, lambda = final_lambda)
  ys = smooth.basis(argvals, realization,  Par_final)
  xfd = ys$fd
  return(xfd)
} ## Functional curves estimate function

cox_lee <- function(d1, d2, n1, n2, nbreaks, permutation, svd = NULL){
  
  
  y1_bar_inner <- apply(d1, 2, mean)
  y2_bar_inner <- apply(d2, 2, mean)
  
  if(!permutation){
    
    covariance_matrix = matrix(0,ncol = nbreaks, nrow = nbreaks)
    y_matrix = rbind(d1,d2)
    y_bar    = apply(y_matrix,2,mean)
    
    for (N in 1:(n1+n2)) {
      covariance_matrix = covariance_matrix + (y_matrix[N,]-y_bar)%*%t((y_matrix[N,]-y_bar))
    }
    estimated_covariance_matrix <- 1/(n1+n2-1)*covariance_matrix
    svd <- svd(estimated_covariance_matrix)
  }
  
  
  
  
  p0<-max(which( (svd$d[2:nbreaks]/svd$d[1] > 10^{-16}) == T))
  
  mid_matrix = matrix(0,nrow = nbreaks, ncol = nbreaks)
  tild_Tksquare_list = c()
  for (k in 1:p0) {
    mid_matrix = mid_matrix + svd$d[k]^{-1} * (svd$u[,k] %*% t(svd$u[,k]))
    tild_Tksquare_list  = append(tild_Tksquare_list, (n1*n2/(n1+n2))*t(y1_bar_inner-y2_bar_inner)%*%mid_matrix %*% (y1_bar_inner-y2_bar_inner))
  }
  
  Tstar_list <- (tild_Tksquare_list-c(1:p0))*(1/sqrt(2*c(1:p0)))
  Sstar_list <- (n - c(1:p0) -1)/((n-2)*c(1:p0))*tild_Tksquare_list
  
  
  Tstar = Tstar_list[which.max(Tstar_list)]
  Sstar = Sstar_list[which.max(Sstar_list)]
  
  return(list(Tstar=Tstar,Sstar = Sstar, svd = svd))
  
} ### Cox and Lee's test (both chisq transform and F transform)



### Compiling

simulation_data  <- cmpfun(simulation_data )
FQnorm  <- cmpfun(FQnorm )
pseudotime_estimate  <- cmpfun(pseudotime_estimate)
gaussian_kernel  <- cmpfun(gaussian_kernel)
nw_estimator_weights  <- cmpfun(nw_estimator_weights)
function_estimate  <- cmpfun(function_estimate)
cox_lee  <- cmpfun(cox_lee)




