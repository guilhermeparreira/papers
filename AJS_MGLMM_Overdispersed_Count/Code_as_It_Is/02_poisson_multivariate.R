library(TMB)
library(mvtnorm)
library(tictoc)
setwd("~/Google Drive/Mestrado/dissertacao/TMB/MultivariateModels")
model <- "02_poisson_multivariate" #1st try
compile(paste0(model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(model))
# Data
# set.seed(1234)
set.seed(2652)
beta1 <- c(2, 0.8)
beta2 <- c(0.5, -1)
beta3 <- c(-.5, 2)
n <- 100
x1 <- rnorm(n)
X <- model.matrix(~x1)
Sigma <- matrix(c(.5, .25, .25,
                  .25, .2, .18,
                  .25, .18, .5), 
                nrow = 3,
                byrow = T)
cov2cor(Sigma)
U <- rmvnorm(n, mean = c(0,0,0), sigma = Sigma)
mu1 <- exp(X%*%beta1 + U[,1])
mu2 <- exp(X%*%beta2 + U[,2])
mu3 <- exp(X%*%beta3 + U[,3])
Y1 <- rpois(n, lambda = mu1)
Y2 <- rpois(n, lambda = mu2)
Y3 <- rpois(n, lambda = mu3)
## Data set
data <- list(X = X,
             Y = cbind(Y1, Y2, Y3))
## Parameter
params <- list(beta = matrix(c(2, 0.8, .5, -1, -.5, 2), ncol = 3),
               U = matrix(0, ncol = 3, nrow = n), 
               rho = c(0.8,0.5,0.5),
               sigma = c(.5,.2,.5))
# Rinterface(paste0(model, ".cpp"))
openmp(n=3)
obj <- MakeADFun(data = data, 
                 parameters = params,
                 DLL = model,
                 method = "BFGS",
                 hessian = T,
                 silent = T
                 ,random = c("U"))
?optim
# print(obj$report())
# length(obj$report()$Y1)
tic()
opt <- do.call("optim", obj)
toc()
rep <- sdreport(obj)
summary(rep, "fixed", p.value = T)
