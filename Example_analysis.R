
source("C:/Users/css22/Desktop/Thesis3/All_required_functions.R")



### Analysis

sce = simulation_data(seed = 123, sigma1 = 0, sigma2= 0)
assays(sce)$norm <- FQnorm(assays(sce)$counts)
slingshot_obj <- pseudotime_estimate(data = assays(sce)$norm)
pseudotime<- slingshot_obj$pseudotime
sce <- slingshot_obj$sce



### Draw plot


plot(c(1:1000),  exp(atan((c(1:1000)-500)/300))+2, type = "l", lwd =5, xlab = "Pseudotime", ylab = "Mean expression" , ylim = c(0,5) )
lines(c(1:1000),  exp(atan((c(1:1000)-500)/300))+1,  type = "l", lwd =5, xlab = "Pseudotime", ylab = "Mean expression"  )
abline(v = 350, lwd = 2, col = "red", lty = "dashed")
abline(v = 700, lwd = 2, col = "red", lty = "dashed")
arrows(x0 = 500, y0 = 2.1, x1 = 500, y1 = 2.7, col = "blue", lwd = 2, length = 0.1)
arrows(x0 = 600, y0 = 2.5, x1 = 600, y1 = 3.1, col = "blue", lwd = 2, length = 0.1)


functional_estimate <- function_estimate(realization = t(assays(sce)$norm[1:100, pseudotime]), 
                                         largest_nbasis = 20, largest_norder =10)
fit = t(eval.fd(c(1:1000), functional_estimate))
plot(x=c(1:1000), fit[1,], type = "l", ylim = c(0, 20), col = "red", ylab = "Expression", xlab = "Pseudotime", 
     main = "Pattern 1 (large variability)")
for (i in 2:100) {

  lines(x=c(1:1000), fit[i,], type ="l", col = "red")

}
functional_estimate <- function_estimate(realization = t(assays(sce)$norm[101:200, pseudotime]), 
                                         largest_nbasis = 20, largest_norder =10)
fit = t(eval.fd(c(1:1000), functional_estimate))
for (i in 1:100) {
  lines(x=c(1:1000), fit[i,], type ="l", col = "blue")
}



functional_estimate <- function_estimate(realization = t(assays(sce)$norm[1:100, pseudotime]), 
                                         largest_nbasis = 50, largest_norder =10)

fit = t(eval.fd(c(1:1000), functional_estimate))
plot(x=c(1:1000), fit[1,], type = "l", ylim = c(0, 20), col = "red", ylab = "Expression", xlab = "Pseudotime", 
     main = "Pattern 2 (large variability)")
for (i in 2:100) {
  
  lines(x=c(1:1000), fit[i,], type ="l", col = "red")
  
}
functional_estimate <- function_estimate(realization = t(assays(sce)$norm[201:300, pseudotime]), 
                                         largest_nbasis = 50, largest_norder =10)
fit = t(eval.fd(c(1:1000), functional_estimate))
for (i in 1:100) {
  lines(x=c(1:1000), fit[i,], type ="l", col = "blue")
}







functional_estimate <- function_estimate(realization = t(assays(sce)$norm[301:400, pseudotime]), 
                                         largest_nbasis = 50, largest_norder =10)

fit = t(eval.fd(c(1:1000), functional_estimate))
plot(x=c(1:1000), fit[1,], type = "l", ylim = c(0,20), col = "red",ylab = "Expression", xlab = "Pseudotime", 
     main = "Pattern 3 (large variability)")
for (i in 2:100) {
  
  lines(x=c(1:1000), fit[i,], type ="l", col = "red")
  
}
functional_estimate <- function_estimate(realization = t(assays(sce)$norm[401:500, pseudotime]), 
                                         largest_nbasis = 50, largest_norder =10)
fit = t(eval.fd(c(1:1000), functional_estimate))
for (i in 1:100) {
  lines(x=c(1:1000), fit[i,], type ="l", col = "blue")
}




functional_estimate <- function_estimate(realization = t(assays(sce)$norm[501:600, pseudotime]), 
                                         largest_nbasis = 50, largest_norder =10)

fit = t(eval.fd(c(1:1000), functional_estimate))
plot(x=c(1:1000), fit[1,], type = "l", ylim = c(0,20), col = "red", ylab = "Expression", xlab = "Pseudotime", 
     main = "Pattern 4 (large variability)")
for (i in 2:100) {
  
  lines(x=c(1:1000), fit[i,], type ="l", col = "red")
  
}
functional_estimate <- function_estimate(realization = t(assays(sce)$norm[601:700, pseudotime]), 
                                         largest_nbasis = 50, largest_norder =10)
fit = t(eval.fd(c(1:1000), functional_estimate))
for (i in 1:100) {
  lines(x=c(1:1000), fit[i,], type ="l", col = "blue")
}

