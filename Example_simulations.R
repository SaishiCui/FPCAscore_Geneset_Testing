
source("C:/Users/css22/Desktop/Thesis3/All_required_functions.R")



p_val_chisq = c()
p_val_F = c()
p_val_manova = c()
p_val_t = c()
for (seed in 1:1010) {
  
  sce = simulation_data(seed = seed, sigma1 = 0, sigma2= 0)
  assays(sce)$norm <- FQnorm(assays(sce)$counts)
  slingshot_obj <- pseudotime_estimate(data = assays(sce)$norm)
  pseudotime<- slingshot_obj$pseudotime
  sce <- slingshot_obj$sce
  
  
  n1 = 50
  n2 = 50
  n = n1+n2
  nbreaks = 600
  start = 1
  end = 300
  set.seed(seed)
  geneset1_idx = sample(c(1:100), n1)
  geneset2_idx = sample(c(201:300), n1)
  
  eva <- seq(1, 300, length.out = nbreaks)
  
  para_up_data <- t(assays(sce)$norm[geneset1_idx , pseudotime[start:end]])
  para_down_data <- t(assays(sce)$norm[geneset2_idx, pseudotime[start:end]])
  
  
  functional_estimate_1 <- function_estimate(realization = para_up_data, 
                                             largest_nbasis = 15, largest_norder = 10, h = 100)
  functional_estimate_2 <- function_estimate(realization = para_down_data, 
                                             largest_nbasis = 15, largest_norder = 10, h = 100)
  
  fit1 =  t(eval.fd(eva, functional_estimate_1))
  fit2 =  t(eval.fd(eva, functional_estimate_2))
  
  
  y1_bar <- apply(fit1, 2, mean)
  y2_bar <- apply(fit2, 2, mean)
  
  
  plot(eva, y1_bar, type = "l", lwd = 5, col= "red", ylim = c(0,20))
  lines(eva, y2_bar, type = "l", lwd = 5, col ="blue")
  
  
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
    
  }
  
  test_statistic <- cox_lee(d1=fit1, d2=fit2, n1=n1, n2=n2, nbreaks = nbreaks, permutation =F, svd =NULL)
  Tstar <- test_statistic$Tstar
  Sstar <- test_statistic$Sstar
  
  
  
  
  nharm = 3
  
  
  FPCA_result <- try({
    pca.fd(functional_estimate_1, nharm = nharm)
    pca.fd(functional_estimate_2, nharm = nharm)
  }, silent = TRUE)
  
  if (class(FPCA_result) == "try-error") {
    next
  }else{
    
    fdpca_score_1 = pca.fd(functional_estimate_1, nharm = nharm, centerfns = F)$scores
    fdpca_score_2 = pca.fd(functional_estimate_2, nharm = nharm, centerfns = F)$scores
  }
  
  
  fpca_test_df <- data.frame( rbind(fdpca_score_1, fdpca_score_2), "group" = as.factor(rep(c("G1","G2"), each = n1)))
  manova_test <- summary(manova(cbind(X1,X2,X3) ~ group  , data = fpca_test_df), test = "Wilks")
  p_manova <- manova_test$stats[1,6]
  p_val_manova <- append(p_val_manova, p_manova)
  
  
  p_t <- hotelling.test(.~group, data = fpca_test_df)$pval
  p_val_t<- append(p_val_t, p_t)
  
  
  
  
  Tstar_perm = rep(0,1000)
  Sstar_perm = rep(0,1000)
  
  for (it in 1:1000) {
    
    set.seed(it^2)
    fit_total =  rbind(fit1, fit2)
    perm_idx1 <- sample(c(1:n),   n1)
    perm_idx2 <- setdiff(c(1:n), perm_idx1)
    
    fit1_perm = fit_total[perm_idx1,]
    fit2_perm = fit_total[perm_idx2,]
    
    test_statistic_perm <- cox_lee(d1=fit1_perm, d2=fit2_perm, n1=n1, n2=n2, nbreaks = nbreaks, permutation = T, svd = test_statistic$svd)
    
    Tstar_perm[it] <- test_statistic_perm$Tstar
    Sstar_perm[it] <- test_statistic_perm$Sstar
    
    
  }
  
  p_val_chisq <- append(p_val_chisq, 1-sum(Tstar>Tstar_perm)/1000)
  p_val_F     <- append(p_val_F, 1-sum(Sstar>Sstar_perm)/1000)
  
  
  
  
  print(seed)
  print(paste("Cox-lee-Chisq:", sum(p_val_chisq<0.05)/length(p_val_chisq)))
  print(paste("Cox-lee-F:", sum(p_val_F<0.05)/length(p_val_F)))
  print(paste("FPCA Manova:", sum(p_val_manova<0.05)/length(p_val_manova)))
  print(paste("FPCA Hotelling:", sum(p_val_t<0.05)/length(p_val_t)))
  
}



sum(p_val_chisq[1:1000] < 0.05)/1000*100
sum(p_val_F[1:1000] < 0.05)/1000*100
sum(p_val_manova[1:1000] < 0.05)/1000*100
sum(p_val_t[1:1000] < 0.05)/1000*100








