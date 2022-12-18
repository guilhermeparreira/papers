library(TMB)
setwd("~/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels")
source("/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/unify_packages.R")

# Simulating Data GLM ------------------------------------------------------------------------
seed <- 1123
beta1 <- c(-3.01, .15)
n <- 50000
set.seed(seed)
x1 <- rnorm(n)
X <- model.matrix(~x1)
mu1 <- exp(X%*%beta1)
## Response variables
phi1 <- .5
l_phi1 <- log(phi1)
True <- c(beta1, phi1)
set.seed(seed)
Y1 <- rnbinom(n, size = exp(l_phi1), mu = mu1)
# Input 
data <- list(Y1 = Y1,
             X = X)
var(Y1);mean(Y1)
parameters <- list(beta1 = beta1,
                   phi = c(l_phi1))
model <- "03_nb_glm" #1st try
compile(paste0(model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(model))
# openmp(n=3)
obj <- MakeADFun(data = data, 
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T)
obj$fn()
# opt <- do.call("optim", obj)
(opt <- nlminb(obj$par, obj$fn, obj$gr,
               control = list(eval.max = 1e9, iter.max = 1e9,
                              abs.tol = 1e-4, rel.tol = 1e-4)))
rep <- sdreport(obj)
fixef <- summary(rep, "fixed")
cor <- summary(rep, "report")
# Joining all together
(pars_used <- TMB_summary(rep_fe, rep_cor, 1))

# install.packages("devtools")
library(mpcmp)
mpcmp::rcomp(10, mu = 500, nu=1)
microbenchmark::microbenchmark(mpcmp::rcomp(10, mu = 500, nu=1), times=10)

model <- "03_nb_glm_II" #1st try
compile(paste0(model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(model))
# openmp(n=3)
# Rinterface(paste0(model,".cpp"))
FreeADFun(obj)
obj <- MakeADFun(data = list(Y = 5), 
                 parameters = list(phi=.5,
                                   mu=.5),
                 DLL = model,
                 hessian = T,
                 silent = T)
### TMB  (AGREE)
obj$fn() 
### R    (AGREE)
dnbinom(5, size = .5, mu = .5)

### Mine (AGREE)
dens.mine <- function(alfa, lambda, x){
  p1 <- (gamma(x+alfa)/(gamma(x+1)*gamma(alfa)))
  p2 <- ((lambda/(lambda+alfa))^alfa)*(alfa/(lambda+alfa))^x
  return(p1*p2)
}
dens.mine(.5,.5,5)
y <- 0:10
plot(dens.mine(.5, .5, y)~y, type = "h")
sum(dens.mine(.5, .5, 0:1000), na.rm = T)
### Book (AGREE!)
dens.book <- function(alfa, lambda, x){
  p1 <- (gamma(x+alfa)/(gamma(x+1)*gamma(alfa)))
  p2 <- ((alfa/(alfa+lambda))^alfa)*(lambda/(lambda+alfa))^x
  return(p1*p2)
}
dens.book(.5, .5, 5)
y <- 0:10
plot(dens.book(0.5, .5, y)~y, type = "h")
sum(dens.book(0.5, .5, 0:1000), na.rm = T)
