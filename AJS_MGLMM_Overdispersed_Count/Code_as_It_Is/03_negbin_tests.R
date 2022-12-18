library(TMB)
library(mvtnorm)
library(tictoc)
library(gamlss.dist)
var_BN <- function(r,p){ #Var of a Neg.Binom ~ (r, p)
  r*(1-p)/p^2
}
mu_BN <- function(r, p){ #Mu of a Neg.Binom ~ (r, p)
  r*(1-p)/p
}
setwd("~/Google Drive/Mestrado/dissertacao/TMB/MultivariateModels")
model <- "06_negbin_tests"
compile(paste0(model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(model))
## Test 2 - Building the Model ---------------------------------------
r = 3
p = .7
n <- 1e5
mu <- mu_BN(r,p)
var <- var_BN(r,p)
Y1 <- rnbinom(n, r, p)
p_cpp <- mu/var
r_cpp <- mu*p_cpp/(1-p_cpp)
obj <- MakeADFun(data = list(Y1 = Y1), 
                 parameters = list(mu = mu,
                                   var = log(var)),
                 DLL = model,
                 hessian = T,
                 method = "BFGS",
                 silent = F)
obj$fn()
opt <- do.call("optim", obj)
opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he)
opt$par
exp(opt$par)
c(mean(Y1), var(Y1))
c(mu, var)

















## Test 1 - Evaluating the density -----------------------------------
## Y ~ NB(prob (p), size (r)) ----------------------
## mu = r*(1-p)/p
## s2 = r*(1-p)/p^2
r = 3
p = .7
dnbinom(x = 2, prob = p, size = r, log = T)
## Y ~ NB (mu, r) ----------------------------------
## p = r/(r+mu)
## s2 = mu + (mu^2)/size
(mu = r*(1-p)/p)
(s2 = r*(1-p)/(p^2))
dnbinom(x = 2, mu = r*(1-p)/p, size = r, log = T)
obj <- MakeADFun(data = list(), 
                 parameters = list(),
                 DLL = model,
                 type = "Fun",
                 silent = F)
print(obj$report())

## ----------------------
# Só para resumir o problema:
  # pelo dnbinom, eu sempre parto do size (r) e do prob (p). Com esses, eu consigo obter o mu e var fácil. Só que sair a partir do mu e do var, e aí gerar a amostra, não tem como. É o caso de implementar (a do GAMLSS faz outra coisa)
## ----------------------
# GAMLSS versus dnbinom() :> They both talk.
size = 1 #k
pi = .5  #pi
mu <- (size*(1-pi))/pi
sigma <- 1/size
gamlss.dist::dNBI(2, mu, sigma, log = T)

# GAMLSS versus TMB() :> They both talk.
