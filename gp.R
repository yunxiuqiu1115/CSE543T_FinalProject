### gaussian process Figure 2.2
library(MASS)
library(pracma)
n = 1000
x = seq(from = -5, to = 5,length.out = n)
kernel = function(x,y) exp(-(x-y)^2/2)
mu = rep(0,n)
var.noise = 1

Kernel = function(x,x.star){
  nrows = length(x)
  ncols = length(x.star)
  K = matrix(rep(0, nrows*ncols), nrows, ncols)
  for (i in 1:nrows) {
    for (j in 1:ncols) {
      K[i,j] = kernel(x[i],x.star[j])
    }
  }
  return(K)
}

covariance = Kernel(x, x)

gp = mvrnorm(n = 3, mu = mu, Sigma = covariance)
plot(x, gp[1,], col="red", type="l", lty=1) 
lines(x, gp[2,], col="blue",lty=1)
lines(x, gp[3,], col="yellow",lty=1)

# y.low = mu - 1.96*diag(covariance)
# y.high = mu + 1.96*diag(covariance)
# polygon(c(x, rev(x)), c(y.high, rev(y.low)), col = "grey", border = NA)

data = c(-4,-3,-1,0,2)
y = c(-2,0,1,2,-1)
# plot(data, y, pch = 4)
# par(new=TRUE)

K.x <- Kernel(data, data)
K.star <- Kernel(data, x)
L <- chol(K.x)
alpha <- mldivide(t(L),mldivide(L,y))
y.pred <- t(K.star)%*%alpha
v <- mldivide(L,K.star)
covariance.pred <- covariance - t(v)%*%v

posterior = mvrnorm(n = 3, mu = y.pred, Sigma = covariance.pred)
plot(x, posterior[1,], col="red", type="l", lty=1) 
lines(x, posterior[2,], col="blue",lty=1)
lines(x, posterior[3,], col="yellow",lty=1)


