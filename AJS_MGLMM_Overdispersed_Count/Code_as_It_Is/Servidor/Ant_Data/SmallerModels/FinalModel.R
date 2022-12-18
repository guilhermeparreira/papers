######################
# Compoisson Multivariate ------------------------------------------------------------------------
######################
rm(list=ls())
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "Ant_Data/SmallerModels"))
source(paste0(mainpath, "unify_packages.R"))
library(TMB)
model <- "compoisson_multi_nu_fixed" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
gc()
############################################
### Esse foi o melhor modelo!!! - Compoisson NU Fixed
############################################
load("compoisson_nu_fixed.RData")
model <- "compoisson_multi_nu_fixed_41" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
parameters <- list(beta1 = estimates_compoisson_sd$ests[c(1,2,3,4,6),1],
                   beta2 = estimates_compoisson_sd$ests[c(1,2,3,4)+6*1,1],
                   beta3 = estimates_compoisson_sd$ests[c(1,2,3,4)+6*2,1],
                   beta4 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*3,1],
                   beta5 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*4,1],
                   beta6 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*5,1],
                   beta7 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*6,1],
                   beta8 = estimates_compoisson_sd$ests[c(1,2,4,6)+6*7,1],
                   beta9 = estimates_compoisson_sd$ests[c(1,2,3,4)+6*8,1],
                   beta10 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*9,1],
                   beta11 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*10,1],
                   beta12 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*11,1],
                   beta13 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*12,1],
                   beta14 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*13,1],
                   beta15 = estimates_compoisson_sd$ests[c(1,2,3,4,5,6)+6*14,1],
                   beta16 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*15,1],
                   beta17 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*16,1],
                   beta18 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*17,1],
                   beta19 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*18,1],
                   beta20 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*19,1],
                   beta21 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*20,1],
                   beta22 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*21,1],
                   beta23 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*22,1],
                   beta24 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*23,1],
                   beta25 = estimates_compoisson_sd$ests[c(1,2,4,6)+6*24,1],
                   beta26 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*25,1],
                   beta27 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*26,1],
                   beta28 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*27,1],
                   beta29 = estimates_compoisson_sd$ests[c(1,2,4,6)+6*28,1],
                   beta30 = estimates_compoisson_sd$ests[c(1,2,3,4)+6*29,1],
                   beta31 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*30,1],
                   beta32 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*31,1],
                   beta33 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*32,1],
                   beta34 = estimates_compoisson_sd$ests[c(1,2,3,4)+6*33,1],
                   beta35 = estimates_compoisson_sd$ests[c(1,2,3,4)+6*34,1],
                   beta36 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*35,1],
                   beta37 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*36,1],
                   beta38 = estimates_compoisson_sd$ests[c(1,2,3,4)+6*37,1],
                   beta39 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*38,1],
                   beta40 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*39,1],
                   beta41 = estimates_compoisson_sd$ests[c(1,2,3,4,6)+6*40,1],
                   U = U,
                   rho = estimates_compoisson_sd$rho,
                   sigma = estimates_compoisson_sd$sigma)

f1 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f2 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover")
f3 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover")
f4 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f5 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f6 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f7 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f8 <- as.formula("Amblyopone.australis ~ Bare.ground + Shrub.cover + Feral.mammal.dung")
f9 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover")
f10 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f11 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f12 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f13 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f14 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f15 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Volume.lying.CWD + Feral.mammal.dung")
f16 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f17 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f18 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f19 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f20 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f21 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f22 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f23 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f24 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f25 <- as.formula("Amblyopone.australis ~ Bare.ground + Shrub.cover + Feral.mammal.dung")
f26 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f27 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f28 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f29 <- as.formula("Amblyopone.australis ~ Bare.ground + Shrub.cover + Feral.mammal.dung")
f30 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover")
f31 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f32 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f33 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f34 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover")
f35 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover")
f36 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f37 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f38 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover")
f39 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f40 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
f41 <- as.formula("Amblyopone.australis ~ Bare.ground + Canopy.cover + Shrub.cover + Feral.mammal.dung")
#Creatrs cov"i" names
for (i in 1:nr){
  assign(paste0("cov", i), c("Intercepto", trimws(strsplit(gsub(".*,.*, ", "", toString(get(paste0("f",i)))), split = " \\+")[[1]])))
}
a <- list()
for (i in 1:nr){
  a[[i]] <- get(paste0("cov", i))
}

dataX <- as.data.frame(cbind(data$X[,-1],                       # Essas são todas as covs
                             Amblyopone.australis = data$Y[,1]) # Esse é para escrever a resposta
                       )
data <- list(Y = data$Y,
             X1 = model.matrix(f1, dataX),
             X2 = model.matrix(f2, dataX),
             X3 = model.matrix(f3, dataX),
             X4 = model.matrix(f4, dataX),
             X5 = model.matrix(f5, dataX),
             X6 = model.matrix(f6, dataX),
             X7 = model.matrix(f7, dataX),
             X8 = model.matrix(f8, dataX),
             X9 = model.matrix(f9, dataX),
             X10 = model.matrix(f10, dataX),
             X11 = model.matrix(f11, dataX),
             X12 = model.matrix(f12, dataX),
             X13 = model.matrix(f13, dataX),
             X14 = model.matrix(f14, dataX),
             X15 = model.matrix(f15, dataX),
             X16 = model.matrix(f16, dataX),
             X17 = model.matrix(f17, dataX),
             X18 = model.matrix(f18, dataX),
             X19 = model.matrix(f19, dataX),
             X20 = model.matrix(f20, dataX),
             X21 = model.matrix(f21, dataX),
             X22 = model.matrix(f22, dataX),
             X23 = model.matrix(f23, dataX),
             X24 = model.matrix(f24, dataX),
             X25 = model.matrix(f25, dataX),
             X26 = model.matrix(f26, dataX),
             X27 = model.matrix(f27, dataX),
             X28 = model.matrix(f28, dataX),
             X29 = model.matrix(f29, dataX),
             X30 = model.matrix(f30, dataX),
             X31 = model.matrix(f31, dataX),
             X32 = model.matrix(f32, dataX),
             X33 = model.matrix(f33, dataX),
             X34 = model.matrix(f34, dataX),
             X35 = model.matrix(f35, dataX),
             X36 = model.matrix(f36, dataX),
             X37 = model.matrix(f37, dataX),
             X38 = model.matrix(f38, dataX),
             X39 = model.matrix(f39, dataX),
             X40 = model.matrix(f40, dataX),
             X41 = model.matrix(f41, dataX),
             nu = rep(log(1.5), nr))
dyn.load(dynlib(paste0(mainpath, model)))
gc()

# APARENTEMENTE TODA A PARTE MANUAL ESTÁ CORRETA!
# dd <- as.data.frame(cbind(Params = sapply(parameters[-c(42:44)], length),
#                           Cov = sapply(a, length),
#                           Data = sapply(data[-c(1,43)], ncol)))
# 
# sum(dd$Params==dd$Cov)
# sum(dd$Params==dd$Data)
# sum(dd$Cov==dd$Data)

# First Round
# data$nu <- rep(log(1.5), nr)
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
# gc()
# FreeADFun(obj)
# openmp(12)
obj$fn()
a <- system.time(fit_compoisson_f <- nlminb(obj$par, obj$fn, obj$gr,
                                            control = list(eval.max = 1e9, iter.max = 1e9,
                                                           abs.tol = 1e-04, rel.tol = 1e-04)))
save.image("compoisson_nu_fixed.RData")
# rep <- sdreport(obj)
# est_com_res <- sv(fit = fit_compoisson, time = a, n_betas = 11, nr = nr, 
# model = model, method = "nlimnb", threads = threads)
# From std
namebetas <- character(0)
lengthresps <- numeric(0)
for(i in 1:nr){
  namebetas <- c(namebetas, get(paste0("cov",i)))
  lengthresps <- c(lengthresps, rep(i, length(get(paste0("cov",i)))))
}
est_com_res_sd <- sv_sd(obj = obj, fit = fit_compoisson_f, time = a, nr = nr, 
                        model = model, method = "nlimnb", threads = threads,
                        namebeta = namebetas,
                        lengthresp = lengthresps)
save.image("compoisson_nu_fixed.RData")
# load("compoisson_nu_fixed.RData")
############################################
### COMPOISSON FULL FIXED
############################################
model <- "compoisson_multi_41" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
data <- data[-43]
parameters <- list(beta1 = est_com_res_sd$ests[1:5,1],
                   beta2 = est_com_res_sd$ests[6:9,1],
                   beta3 = est_com_res_sd$ests[10:13,1],
                   beta4 = est_com_res_sd$ests[14:18,1],
                   beta5 = est_com_res_sd$ests[19:23,1],
                   beta6 = est_com_res_sd$ests[24:28,1],
                   beta7 = est_com_res_sd$ests[29:33,1],
                   beta8 = est_com_res_sd$ests[34:37,1],
                   beta9 = est_com_res_sd$ests[38:41,1],
                   beta10 = est_com_res_sd$ests[grep("Y10", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta11 = est_com_res_sd$ests[grep("Y11", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta12 = est_com_res_sd$ests[grep("Y12", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta13 = est_com_res_sd$ests[grep("Y13", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta14 = est_com_res_sd$ests[grep("Y14", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta15 = est_com_res_sd$ests[grep("Y15", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta16 = est_com_res_sd$ests[grep("Y16", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta17 = est_com_res_sd$ests[grep("Y17", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta18 = est_com_res_sd$ests[grep("Y18", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta19 = est_com_res_sd$ests[grep("Y19", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta20 = est_com_res_sd$ests[grep("Y20", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta21 = est_com_res_sd$ests[grep("Y21", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta22 = est_com_res_sd$ests[grep("Y22", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta23 = est_com_res_sd$ests[grep("Y23", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta24 = est_com_res_sd$ests[grep("Y24", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta25 = est_com_res_sd$ests[grep("Y25", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta26 = est_com_res_sd$ests[grep("Y26", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta27 = est_com_res_sd$ests[grep("Y27", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta28 = est_com_res_sd$ests[grep("Y28", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta29 = est_com_res_sd$ests[grep("Y29", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta30 = est_com_res_sd$ests[grep("Y30", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta31 = est_com_res_sd$ests[grep("Y31", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta32 = est_com_res_sd$ests[grep("Y32", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta33 = est_com_res_sd$ests[grep("Y33", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta34 = est_com_res_sd$ests[grep("Y34", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta35 = est_com_res_sd$ests[grep("Y35", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta36 = est_com_res_sd$ests[grep("Y36", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta37 = est_com_res_sd$ests[grep("Y37", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta38 = est_com_res_sd$ests[grep("Y38", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta39 = est_com_res_sd$ests[grep("Y39", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta40 = est_com_res_sd$ests[grep("Y40", row.names(est_com_res_sd$ests[1:196,])),1],
                   beta41 = est_com_res_sd$ests[grep("Y41", row.names(est_com_res_sd$ests[1:196,])),1],
                   U = U,
                   rho = est_com_res_sd$rho,
                   sigma = est_com_res_sd$sigma,
                   nu = rep(log(1), nr))

# APARENTEMENTE TODA A PARTE MANUAL ESTÁ CORRETA!
# dd <- as.data.frame(cbind(Params = sapply(parameters[c(1:41)], length),
#                           Data = sapply(data[-c(1)], ncol)))
# 
# sum(dd$Params==dd$Cov)
# sum(dd$Cov==dd$Data)


obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
# gc()
# FreeADFun(obj)
# openmp(12)
obj$fn()
a <- system.time(fit_compoisson_full <- nlminb(obj$par, obj$fn, obj$gr,
                                               control = list(eval.max = 1e9, iter.max = 1e9,
                                                              abs.tol = 1e-04, rel.tol = 1e-04)))
# fit_compoisson_f$objective
save.image("compoisson_nu_fixed.RData")
est_com_res_sd_full <- sv_sd(obj = obj, fit = fit_compoisson_full, time = a, nr = nr, 
                        model = model, method = "nlimnb", threads = threads,
                        namebeta = namebetas, 
                        lengthresp = lengthresps)
save.image("compoisson_nu_fixed.RData")
# Funcionou! :D
est_com_res_sd_full_easy <- sv2(fit = fit_compoisson_full, time = a, nr = nr, rows = nrow(data$Y),
                                model = model, method = "nlminb", threads = threads)
#####################################
#### NEGBIN
#####################################
load("compoisson_nu_fixed.RData")
model <- "negbin_multivariate_41"
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
parameters <- est_com_res_sd_full_easy
parameters <- parameters[-length(parameters)]
names(parameters)[nr+4] <- "phi"
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
obj$fn()
a <- system.time(fit_negbin_full_1 <- nlminb(obj$par, obj$fn, obj$gr,
                                               control = list(eval.max = 1e9, iter.max = 1e9,
                                                              abs.tol = 1e-04, rel.tol = 1e-04)))
save.image("compoisson_nu_fixed.RData")
est_nb_res_sd_full <- sv_sd(obj = obj, fit = fit_negbin_full_1, time = a, nr = nr, 
                             model = model, method = "nlimnb", threads = threads,
                             namebeta = namebetas, 
                             lengthresp = lengthresps)
est_nb_res_sd_full_easy <- sv2(fit = fit_negbin_full_1, time = a, nr = nr, rows = nrow(data$Y),
                                model = model, method = "nlminb", threads = threads)
save.image("compoisson_nu_fixed.RData")
#####################################
#### POISSON
#####################################
load("compoisson_nu_fixed.RData")
model <- "02_poisson_multivariate_41"
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
parameters <- est_nb_res_sd_full_easy
# REFAZER
parameters <- parameters[-c((length(parameters)-1):(length(parameters)))]
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
obj$fn()
a <- system.time(fit_poisson_1 <- nlminb(obj$par, obj$fn, obj$gr,
                                             control = list(eval.max = 1e9, iter.max = 1e9,
                                                            abs.tol = 1e-04, rel.tol = 1e-04)))
save.image("compoisson_nu_fixed.RData")
est_poisson_res_sd_full <- sv_sd(obj = obj, fit = fit_poisson_1, time = a, nr = nr, 
                            model = model, method = "nlimnb", threads = threads,
                            namebeta = namebetas, 
                            lengthresp = lengthresps)
est_poisson_res_sd_full_easy <- sv2(fit = fit_poisson_1, time = a, nr = nr, rows = nrow(data$Y),
                               model = model, method = "nlminb", threads = threads)
save.image("compoisson_nu_fixed.RData")

#################################################################
# Likelihood summary for all
#################################################################
est_poisson_res_sd_full$summary$Model <- "Poisson"
est_poisson_res_sd_full$summary$np <- nrow(est_poisson_res_sd_full$ests)
est_poisson_res_sd_full$summary$SE <- TRUE
m1 <- est_poisson_res_sd_full$summary
est_nb_res_sd_full$summary$Model <- "NB"
est_nb_res_sd_full$summary$np <- nrow(est_nb_res_sd_full$ests)
est_nb_res_sd_full$summary$SE <- TRUE
m2 <- est_nb_res_sd_full$summary
est_com_res_sd_full$summary$Model <- "COM-Poisson Full"
est_com_res_sd_full$summary$np <- nrow(est_com_res_sd_full$ests)
est_com_res_sd_full$summary$SE <- TRUE
m3 <- est_com_res_sd_full$summary
est_com_res_sd$summary$Model <- "COM-Poisson Fixed nu"
est_com_res_sd$summary$np <- nrow(est_com_res_sd$ests)
est_com_res_sd$summary$SE <- TRUE
m4 <- est_com_res_sd$summary
mm <- rbind(m1, m2, m3, m4)
mm$Loglik <- mm$Loglik*-1
mm$AIC <- -2*mm$Loglik + 2*mm$np
mm$BIC <- -2*mm$Loglik + log(nrow(data$Y))*mm$np
mm <- mm[, c(6,7,9,10,1,8,5)]
row.names(mm) <- NULL
mm$Optimizer <- with(mm, ifelse(Optimizer=="bfgs", "BFGS", 
                                ifelse(Optimizer%in%c("nlminb","nlimnb"), "PORT", Optimizer)))
mm <- mm[, -c(7)]
library(knitr)
kable(mm, 
      caption = "Model fit measures for ANT data from the best parametrization for each distribution",
      label = "antfit2",
      align = c("l", rep("c", ncol(mm)-1)),
      booktabs = T,
      digits = 2,
      row.names = F,
      format = "latex")


#################################################################
# SECOND ROUND, BECAUSE COMPOISSON_FULL HAD NaN std.error for some parameters
#################################################################
#################################################################
# FULL COMPOISSON
#################################################################
load("compoisson_nu_fixed.RData")
model <- "compoisson_multi_41" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
only_est <- as.data.frame(est_com_res_sd_full$ests[1:196, ])
lists_betas <- list()
lists_mat <- list()
namesbetas2 <- character()
lengthresps2 <- numeric()
# Cria os parametros e os betas
for (i in 1:nr){
  daux <- only_est[grep(paste0("_Y",i,"$"), row.names(only_est)), ]
  daux <- rbind(daux[1, ], daux[2:nrow(daux), ][!is.nan(daux$`Std. Error`[-1]), ])
  lists_betas[[i]] <- daux[,1]
  covnames <- gsub("_.*", "", row.names(daux))
  if (length(covnames)==1 & covnames[1] == "Intercepto"){
    ff <- as.formula("Amblyopone.australis ~ 1")
  } else {
    ff <- as.formula(paste0("Amblyopone.australis ~ ", paste0(covnames[-1], collapse = " + ")))
  } 
  lists_mat[[i]] <- model.matrix(ff, dataX)
  namesbetas2 <- c(namesbetas2, covnames)
  lengthresps2 <- c(lengthresps2, rep(i, length(covnames)))
}
lists_betas[[nr+1]] <- U
lists_betas[[nr+2]] <- est_com_res_sd_full$rho
lists_betas[[nr+3]] <- est_com_res_sd_full$sigma
lists_betas[[nr+4]] <- est_com_res_sd_full$disp
names(lists_betas) <- c(paste0("beta",1:nr), "U", "rho", "sigma", "nu")
lists_mat[[nr+1]] <- data$Y
names(lists_mat) <- c(paste0("X", 1:nr), "Y")
# Correctin too big values of beta (greater than 5)
# Com os valores iniciais anteriores não iniciou, aí tive que zerar todos os betas.
for (i in 1:nr){
  ll <- length(lists_betas[[i]])
  for (j in 1:ll){
    lists_betas[[i]][j] <- 0
  }
}
obj <- MakeADFun(data = lists_mat,
                 parameters = lists_betas,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
obj$fn()
a <- system.time(fit_compoisson_f2 <- nlminb(obj$par, obj$fn, obj$gr,
                                            control = list(eval.max = 1e9, iter.max = 1e9,
                                                           abs.tol = 1e-04, rel.tol = 1e-04)))
####################################################
#### Redo, because the likelihood value was terrible (2635)
####################################################

ests_firsts <- sv2(fit_compoisson_f2, a, nr, model, method = "nlminb", 
                   threads = threads, rows = nrow(data$Y))
ests_firsts <- ests_firsts[-length(ests_firsts)]
names(ests_firsts)[length(ests_firsts)] <- "nu"
obj <- MakeADFun(data = lists_mat,
                 parameters = ests_firsts,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
# obj$fn()
a <- system.time(fit_compoisson_f2 <- nlminb(obj$par, obj$fn, obj$gr,
                                             control = list(eval.max = 1e9, iter.max = 1e9,
                                                            abs.tol = 1e-04, rel.tol = 1e-04)))
# Na mão dps, pois sobrescrevi objetos
m5 <- data.frame(Model = "COM-Poisson Full after filter",
                  np = sum(sapply(ests_firsts, length))-1230,
                 AIC = -2*-2635 + 2*sum(sapply(ests_firsts, length))-1230,
                 BIC = -2*-2635 + (log(nrow(data$Y))*sum(sapply(ests_firsts, length))-1230),
                 Loglik = -2635,
                 SE = "Não calculado")
####################################################
#### Redo, AGAIN, LAST TRY (2845) - Even worse
####################################################
ests_firsts2 <- sv2(fit_compoisson_f2, a, nr, model, method = "nlminb", 
                   threads = threads, rows = nrow(data$Y))
ests_firsts2 <- ests_firsts2[-length(ests_firsts2)]
names(ests_firsts2)[length(ests_firsts2)] <- "nu"
obj <- MakeADFun(data = lists_mat,
                 parameters = ests_firsts2,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
# obj$fn()
a <- system.time(fit_compoisson_f22 <- nlminb(obj$par, obj$fn, obj$gr,
                                             control = list(eval.max = 1e9, iter.max = 1e9,
                                                            abs.tol = 1e-04, rel.tol = 1e-04)))
# Refazer com esse valor inicial, de cima.

#################################################################
## Compoisson fixed nu (1469)
#################################################################
model <- "compoisson_multi_nu_fixed_41" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
ests_firsts_fixed <- ests_firsts
ests_firsts_fixed <- ests_firsts_fixed[-length(ests_firsts_fixed)]
lists_mat[[43]] <- rep(log(1.5), nr)
names(lists_mat)[[43]] <- "nu"
obj <- MakeADFun(data = lists_mat,
                 parameters = ests_firsts_fixed,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
# obj$fn()
a <- system.time(fit_compoisson_f2 <- nlminb(obj$par, obj$fn, obj$gr,
                                             control = list(eval.max = 1e9, iter.max = 1e9,
                                                            abs.tol = 1e-04, rel.tol = 1e-04)))
# Did not run
est_cmp_res_sd_nu_2 <- sv_sd(obj = obj, fit = fit_compoisson_f2, time = a, nr = nr, 
                                 model = model, method = "nlimnb", threads = threads,
                                 namebeta = namesbetas2, 
                                 lengthresp = lengthresps2)
save.image("compoisson_nu_fixed.RData")



#################################################################
# Likelihood summary for all
#################################################################

est_poisson_res_sd_full$summary$Model <- "Poisson"
est_poisson_res_sd_full$summary$np <- nrow(est_poisson_res_sd_full$ests)
est_poisson_res_sd_full$summary$SE <- TRUE
m21 <- est_poisson_res_sd_full$summary
est_nb_res_sd_full$summary$Model <- "NB"
est_nb_res_sd_full$summary$np <- nrow(est_nb_res_sd_full$ests)
est_nb_res_sd_full$summary$SE <- TRUE
m22 <- est_nb_res_sd_full$summary
est_com_res_sd_full$summary$Model <- "COM-Poisson Full"
est_com_res_sd_full$summary$np <- nrow(est_com_res_sd_full$ests)
est_com_res_sd_full$summary$SE <- TRUE
m23 <- est_com_res_sd_full$summary
est_com_res_sd$summary$Model <- "COM-Poisson Fixed nu"
est_com_res_sd$summary$np <- nrow(est_com_res_sd$ests)
est_com_res_sd$summary$SE <- TRUE
m24 <- est_com_res_sd$summary
### Models with extra filter (only com ṕoisson so far)
m25 <- m5 # add straight to the end

# nu fixed
m26 <- data.frame(Model = "COM-Poisson Fixed nu after filter",
                  np = length(fit_compoisson_f2$par),
                  AIC = -2*-fit_compoisson_f2$objective + 2*length(fit_compoisson_f2$par),
                  BIC = -2*-fit_compoisson_f2$objective + log(nrow(data$Y))*length(fit_compoisson_f2$par),
                  Loglik = -fit_compoisson_f2$objective,
                  SE = "Demorou e parei")


mm2 <- rbind(m21, m22, m23, m24)
mm2$Loglik <- mm2$Loglik*-1
mm2$AIC <- -2*mm2$Loglik + 2*mm2$np
mm2$BIC <- -2*mm2$Loglik + log(nrow(data$Y))*mm2$np
mm2 <- mm2[, c(6,7,9,10,1,8,5)]
row.names(mm2) <- NULL
mm2$Optimizer <- with(mm2, ifelse(Optimizer=="bfgs", "BFGS", 
                                ifelse(Optimizer%in%c("nlminb","nlimnb"), "PORT", Optimizer)))
mm2 <- mm2[, -c(7)]
mm2 <- rbind(mm2, m25, m26)
library(knitr)
mm2 <- mm2[-c(5,6), ]


uma.maiuscula <- function(x) {
  paste(toupper(substring(x, 1,1)), substring(x, 2, nchar(x)),
        sep="", collapse=" ")
}
uma.maiuscula <- Vectorize(uma.maiuscula, "x")
mm2$Model <- as.character(uma.maiuscula(gsub("Full ","", gsub("NB", "negative binominal", gsub("Disp ", "dispersion ", gsub("Var ", "variance ", mm2$Model))))))
mm2$SE <- ifelse(mm2$SE=="TRUE", "\\checkmark", "\\xmark")
mm2[3:4, 1] <- c("COM-Poisson", "Fixed dispersion COM-Poisson")

kable(mm2, 
      caption = toupper("Model fit measures for ANT data from the best parametrization for each distribution"),
      label = "antfit2",
      align = c("l", rep("c", ncol(mm)-1)),
      booktabs = T,
      digits = 2,
      row.names = F,
      escape = F,
      format = "latex")
#### TRV between COM-Poisson Full and COM-Poisson Fixed nu###
LR <- 2*(mm[3,5]-mm[4,5])
df <- mm[3,2]-mm[4,2]
pchisq(LR, df, lower.tail = F)
# library(epicalc)
# model0 <- glm(case ~ induced + spontaneous, family=binomial, data=infert)
# model1 <- glm(case ~ induced, family=binomial, data=infert)
# pchisq(2*(logLik(model0)-logLik(model1)), 1, lower.tail = F)
# 
# lrtest (model0, model1)

load("compoisson_nu_fixed.RData")
library(ggplot2)
library(dplyr)
library(latex2exp)
#################################################################
############# Primary graphic
#################################################################
df <- data.frame(Sigma = c(est_com_res_sd_full$sigma, est_nb_res_sd_full$sigma),
                 Nu = c(exp(est_com_res_sd_full$disp), exp(est_nb_res_sd_full$disp)),
                 Dist = c(rep(paste0("COM-Poisson (", expression(hat(nu)), ")"), 41), 
                          rep(paste0("NB (", expression(hat(phi)), ")"), 41)))
df %>% 
  group_by(Dist) %>% 
  summarise(corre = cor(Sigma, Nu))
df[which.max(df$Nu), 2] <- NA # 46072846421 Excluded
# df
ggplot(data = df, aes(x = Nu, y = Sigma)) +
  geom_point() +
  geom_smooth(se = F) +
  # geom_vline(xintercept = 1) +
  facet_wrap(~Dist, scales = "free_x",
             labeller = label_parsed) +
  labs(x = "Dispersion parameter",
       y = TeX("$\\hat{\\sigma}$")) +
  theme_bw() +
  theme(text = element_text(size = 13))
pathfig <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/Texto/LaTeX/Figuras"
ggsave(paste0(pathfig, "/ant_disp_var.pdf"), width = 20, height = 10, units = "cm")


#################################################################
############# Graphic for all 
#################################################################

np_poisson <- nrow(est_poisson_res_sd_full$ests)
np_cmp <- nrow(est_com_res_sd_full$ests)
dg <- data.frame(Parametros = c(row.names(est_poisson_res_sd_full$ests), row.names(est_nb_res_sd_full$ests), row.names(est_com_res_sd_full$ests)),
                 Estimativas = c(est_poisson_res_sd_full$ests[,1], est_nb_res_sd_full$ests[,1], est_com_res_sd_full$ests[,1]),
                 ErroPadrao = c(est_poisson_res_sd_full$ests[,2], est_nb_res_sd_full$ests[,2], est_com_res_sd_full$ests[,2]),
                 Distribuição = c(rep("Poisson", np_poisson), rep("NB", np_cmp), rep("COM-Poisson", np_cmp)),
                 Tipo = c(c(rep("Média", 196), rep("Variância", 41), rep("Correlação", 820)), #Poisson
                          c(rep("Média", 196), rep("Variância", 41), rep("Correlação", 820), rep("Dispersão", 41)), #NB
                          c(rep("Média", 196), rep("Variância", 41), rep("Correlação", 820), rep("Dispersão", 41)))) #COM

# (STARTING HERE!!!!)
# head(dg,20)
# dg[c(33:35, 52:54), 1] <- c(paste0("exp(phi)", 1:3), paste0("exp(nu)", 1:3))
dg[dg$Tipo=="Dispersão", 1] <- c(paste0("exp(phi)", 1:41), paste0("exp(nu)", 1:41))
dg$Response <- ifelse(dg$Tipo=="Média", gsub(".*_Y", "", dg$Parametros),
                             ifelse(dg$Tipo=="Variância", gsub("sig", "", dg$Parametros),
                                    ifelse(dg$Tipo=="Correlação", gsub("r_", "", dg$Parametros),
                                           ifelse(dg$Tipo=="Dispersão", gsub("[^0-9.-]", "", dg$Parametros), "BU"))))


varsnames <- c(" 1.Amblyopone", " 2.Aphaenogaster", " 3.Camponotus.Ci", " 4.Camponotus.Cl", 
               " 5.Camponotus.Co", " 6.Camponotus.Ni", " 7.Camponotus.Ns", " 8.Cardiocondyla", 
               " 9.Crematogaster", "10.Heteroponera", "11.Iridomyrmex.Bi", "12.Iridomyrmex.Dr", 
               "13.Iridomyrmex.Mj", "14.Iridomyrmex.Pu", "15.Iridomyrmex.Ru", 
               "16.Iridomyrmex.Su", "17.Iridomyrmex.Ss", "18.Melophorus.E", 
               "19.Melophorus.F", "20.Melophorus.H", "21.Meranoplus.A", "22.Monomorium.L", 
               "23.Monomorium.Ro", "24.Monomorium.Sy", "25.Myrmecia", "26.Notoncus.Ca", 
               "27.Notoncus.Ec", "28.Nylanderia", "29.Ochetellus", "30.Paraparatrechina", 
               "31.Pheidole.A", "32.Pheidole.B", "33.Pheidole.E", "34.Pheidole.J", 
               "35.Polyrhachis", "36.Rhytidoponera.A", "37.Rhytidoponera.B", 
               "38.Solenopsis", "39.Stigmacros", "40.Tapinoma", "41.Tetramorium")
respnames <- data.frame(Number = as.character(1:41),
                        Names = factor(varsnames, 
                                       levels = varsnames))
# Real names
head(dg)
# dg <- dg[, -7]
dg <- left_join(dg, respnames, by = c("Response" = "Number"))
## Only beta (exp() or only beta()?)
dg$Response <- as.numeric(dg$Response)
dg$NAN <- factor(as.character(ifelse(is.nan(dg$ErroPadrao), 1, 16)))
dgmedia <- dg %>% 
  filter(Tipo=="Média") %>% 
  mutate(Parametros = factor(gsub("_Y\\d+", "", Parametros), 
                             levels = c("Intercepto", "Bare.ground", "Canopy.cover", "Shrub.cover", "Feral.mammal.dung", "Volume.lying.CWD"),
                             labels = c("Intercept", "Bare.ground", "Canopy.cover", "Shrub.cover", "Feral.mammal.dung", "Volume.lying.CWD")))

dgmedia$ShapeFinal <- factor(with(dgmedia, ifelse(Distribuição == 'Poisson' & NAN == 16, 'Poisson',
                                                  ifelse(Distribuição == 'Poisson' & NAN == 1, 'Poisson without SE',
                                                         ifelse(Distribuição == 'NB' & NAN == 16, 'NB',
                                                                ifelse(Distribuição == 'NB' & NAN == 1, 'NB without SE',
                                                                       ifelse(Distribuição == 'COM-Poisson' & NAN == 16, 'COM-Poisson',
                                                                              ifelse(Distribuição == 'COM-Poisson' & NAN == 1, 'COM-Poisson without SE', 25
                                                                              ))))))))
# shapes <- c(15, 0, 16, 1, 17, 2)
# names(shapes) <- c('Poisson', 'Poisson without SE', 'NB', 'NB without SE', 'COM-Poisson', 'COM-Poisson without SE')
shapes <- c(15, 16, 1, 17, 2)
names(shapes) <- c('Poisson', 'NB', 'NB without SE', 'COM-Poisson', 'COM-Poisson without SE')
########################################
# Graphic for average only!!!
########################################
plot.mean <- function(data, inf, sup, shapes){
  data %>% 
    filter(between(Response, inf, sup)) %>% 
    filter(Parametros != 'Intercept') %>% 
    ggplot() +
    geom_pointrange(aes(x = Parametros, 
                        y = Estimativas,
                        ymin = Estimativas-1.96*ErroPadrao,
                        ymax = Estimativas+1.96*ErroPadrao,
                        shape = ShapeFinal),
                    # It controls the vertical distance between models
                    position = position_dodge(width = .6),
                    size = .5) +
    scale_shape_manual(values = shapes) +
    guides(shape = guide_legend(override.aes = list(linetype = "blank"))) +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          text = element_text(size = 17),
          legend.position = "bottom",
          legend.title = element_blank()) +
    facet_wrap(~Names, scales = "free_x") +
    coord_flip() +
    labs(y = TeX('$\\hat{\\beta}\\pm {1.96 SE}$'))
}

pathfig <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/Texto/LaTeX/Figuras"
plot.mean(dgmedia, 1, 12, shapes)
ggsave(paste0(pathfig, "/ant_beta_final_model_1_12_v_artigo.pdf"), width = 30, height = 24, units = "cm")
plot.mean(dgmedia, 13, 24, shapes)
ggsave(paste0(pathfig, "/ant_beta_final_model_13_24_v_artigo.pdf"), width = 30, height = 24, units = "cm")
plot.mean(dgmedia, 25, 36, shapes)
ggsave(paste0(pathfig, "/ant_beta_final_model_25_36_v_artigo.pdf"), width = 30, height = 24, units = "cm")
plot.mean(dgmedia, 37, 42, shapes)
ggsave(paste0(pathfig, "/ant_beta_final_model_37_42_v_artigo.pdf"), width = 30, height = 24, units = "cm")

plot.mean.old <- function(data, inf, sup){
  data %>% 
    filter(between(Response, inf, sup)) %>% 
    ggplot() +
    geom_pointrange(aes(x = Parametros, 
                        y = Estimativas,
                        ymin = Estimativas-1.96*ErroPadrao,
                        ymax = Estimativas+1.96*ErroPadrao,
                        linetype = Distribuição,
                        shape = NAN),
                    # It controls the vertical distance between models
                    position = position_dodge(width = .6),
                    
                    size = .5) +
    scale_shape_manual(values = c(1,20),
                       labels = c("Without SE", "With SE")) +
    scale_linetype_manual(values=c("solid", "dashed", "dotted")) +
    guides(linetype = guide_legend(override.aes = list(shape = NA)),
           shape = guide_legend(override.aes = list(linetype = "blank"))) +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          text = element_text(size = 17),
          legend.position = "bottom") +
    facet_wrap(~Names, scales = "free_x") +
    coord_flip() +
    labs(y = TeX('$\\hat{\\beta}\\pm {1.96 SE}$'),
         linetype = "Model",
         shape = "SE")
}

#####################
# Only standard deviation
#####################

dg %>% 
  filter(Tipo %in% c("Variância")) %>% 
  ggplot() +
  geom_pointrange(aes(x = Names, 
                      y = Estimativas,
                      ymin = Estimativas-1.96*ErroPadrao,
                      ymax = Estimativas+1.96*ErroPadrao,
                      linetype = Distribuição,
                      shape = NAN),
                  # It controls the vertical distance between models
                  position = position_dodge(width = .6),
                  
                  size = .6) +
  coord_flip() +
  facet_wrap(~Response<21, scales = "free") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        text = element_text(size = 15),
        legend.position = "bottom",
        strip.text = element_blank()) +
  scale_shape_manual(values = c(1,20),
                     labels = c("Without SE", "With SE")) +
  scale_linetype_manual(values=c("solid", "dashed", "dotted")) +
  guides(linetype = guide_legend(override.aes = list(shape = NA)),
         shape = guide_legend(override.aes = list(linetype = "blank"))) +
  labs(y = TeX('$\\hat{\\sigma}\\pm {1.96 SE}$'),
       linetype = "Model",
       shape = "SE")
pathfig <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/Texto/LaTeX/Figuras"
ggsave(paste0(pathfig, "/ant_beta_final_model_var.pdf"), width = 25, height = 20, units = "cm")


########################################
# Dispersion and standard error printing
########################################
# Except Dispersion
library(tidyr)
library(kableExtra)
myHeader1 <- c(1,rep(2,2))
names(myHeader1) <- c(" ", "Negative binomial($\\\\phi$)", "COM-Poisson($\\\\nu$)")
myHeader2 <- c(1,rep(2,3))
names(myHeader2) <- c(" ", "Poisson($\\\\sigma$)", "NB($\\\\sigma$)", "COM-Poisson($\\\\sigma$)")
dgdisp_var <- dg %>% 
  filter(Tipo %in% c("Dispersão", "Variância")) %>% 
  select(-Response, -Parametros, -NAN) %>% 
  # mutate(Response = rep(c("Nmsp", "Nmosp", "Nspfy"), times = 2)) %>% 
  pivot_wider(names_from = c(Distribuição, Tipo),
              values_from = c("Estimativas", "ErroPadrao")) %>%
  select(Names, "Estimativas_NB_Dispersão", "ErroPadrao_NB_Dispersão", 
         "Estimativas_COM-Poisson_Dispersão", "ErroPadrao_COM-Poisson_Dispersão",
         "Estimativas_Poisson_Variância", "ErroPadrao_Poisson_Variância",
         "Estimativas_NB_Variância", "ErroPadrao_NB_Variância",
         "Estimativas_COM-Poisson_Variância", "ErroPadrao_COM-Poisson_Variância") %>% 
  mutate(Estimativas_NB_Dispersão = format(Estimativas_NB_Dispersão, scientific = T, digits = 2),
         ErroPadrao_NB_Dispersão = format(ErroPadrao_NB_Dispersão, scientific = T, digits = 2))
dgdisp <- dgdisp_var[, 1:5]
dgvar <- dgdisp_var[, c(1,6:ncol(dgdisp_var))]
  # mutate_if(is.numeric, funs(as.character(signif(., 2)))) %>% 
kable(dgdisp,
      caption = toupper("Dispersion of parameter estimates and standard errors (SE) for each model and outcome of ANT data"),
      label = "antdisp",
      booktabs = T,
      digits = 3,
      row.names = F,
      escape = F,
      format = "latex",
      col.names = c("Outcome",rep(c("Estimate", "SE"), 2)),
      align = "c") %>% 
  add_header_above(header = myHeader1,
                   escape = F)

kable(dgvar,
      caption = toupper("Standard deviation of random effect estimates and standard errors (SE) for each model and outcome of ANT data"),
      label = "antvar",
      booktabs = T,
      digits = 3,
      row.names = F,
      escape = F,
      format = "latex",
      col.names = c("Outcome",rep(c("Estimate", "SE"), 3)),
      align = "c") %>% 
  add_header_above(header = myHeader2,
                   escape = F)


########################################
# Correlation
########################################

dgcorr <- dg %>% 
  filter(Tipo == "Correlação")
dgcorr$Resp <- gsub("r_", "", dgcorr$Parametros)
dgcorr$Resp <- gsub("Y", "", dgcorr$Resp)
dgcorr <- dgcorr[, -c(6,7)]
resps2 <- stringr::str_split_fixed(dgcorr$Resp, "_", 2)
colnames(resps2) <- c("Var1", "Var2")
dgcorr <- cbind(dgcorr, resps2)
dgcorr$Var1 <- as.numeric(dgcorr$Var1)
dgcorr$Var2 <- as.numeric(dgcorr$Var2)
dgcorr$IcInf <- dgcorr$Estimativas-1.96*dgcorr$ErroPadrao
dgcorr$IcSup <- dgcorr$Estimativas+1.96*dgcorr$ErroPadrao
dgcorr$Valorp <- pnorm(abs(dgcorr$Estimativas/dgcorr$ErroPadrao), lower.tail = F)*2


# Agora operação para colocar em formato de matrix
df_to_mat <- function(dd){
  mat_est <- matrix(numeric(), nrow = 41, ncol = 41, dimnames = list(1:41, 1:41))
  mat_se <- matrix(numeric(), nrow = 41, ncol = 41, dimnames = list(1:41, 1:41))
  mat_inf <- matrix(numeric(), nrow = 41, ncol = 41, dimnames = list(1:41, 1:41))
  mat_sup <- matrix(numeric(), nrow = 41, ncol = 41, dimnames = list(1:41, 1:41))
  mat_p <- matrix(numeric(), nrow = 41, ncol = 41, dimnames = list(1:41, 1:41))
  for (i in 1:41){
    for (j in i:41){
      if (i == j){
        mat_est[i,j] <- 1
        mat_se[i,j] <- 0
        mat_p[i, j] <- 0
      } else{
        dgcorrfilt <- dd[dd$Var1 %in% i & dd$Var2 %in% j, , drop = F]
        mat_est[i, j] <- dgcorrfilt[, 2]
        mat_est[j, i] <- dgcorrfilt[, 2]
        
        mat_se[i, j] <- dgcorrfilt[, 3]
        mat_inf[i, j] <- as.numeric(ifelse(dgcorrfilt[, 10]>1, 1, ifelse(dgcorrfilt[, 10]<(-1),-1, dgcorrfilt[, 10])))
        mat_sup[i, j] <- as.numeric(ifelse(dgcorrfilt[, 11]>1, 1, ifelse(dgcorrfilt[, 11]<(-1),-1, dgcorrfilt[, 11])))
        
        mat_p[i, j] <- 1-dgcorrfilt[, 12]
        mat_p[j, i] <- 1-dgcorrfilt[, 12]
      }
    }
  }
  ll <- list(mat_est, mat_se, mat_inf, mat_sup, mat_p)
  names(ll) <- c("est", "se", "inf", "sup", "p")
  return(ll)
}
library(corrplot) # REVER DAQUI!!!!
gera_grafico <- function(data, namefile, pathfig){
  dgcorrP <- data
  corrp <- df_to_mat(dgcorrP)
  # pdf(file = paste0(pathfig, namefile), width = 6, height = 5.5)
  graph1 <- corrplot(corrp$est, 
                     # mar = c(0,0,0,0), 
                     col = COL1('Greys'),
                     tl.cex = 0.6, 
                     cl.cex = .7,
                     tl.col = 'black',
                     tl.srt = 90,
                     order = 'hclust',
                     pch = "*",
                     type = "upper",
                     diag = F,
                     outline = F,
                     pch.cex = 1.7,       # Controls the size of the star
                     p.mat = 1-corrp$p, #Opposite
                     sig.level = .95)   #Opposite
  print(graph1)
  dev.off()
}

pathfig = "/home/guilherme/Dropbox/SuperDisperso_AJS/Supplementary material for AJS/"
gera_grafico(dgcorr[dgcorr$Distribuição=="Poisson", ], "/cor_ant_poisson.pdf", pathfig)
gera_grafico(dgcorr[dgcorr$Distribuição=="NB", ], "/cor_ant_NB.pdf", pathfig)
gera_grafico(dgcorr[dgcorr$Distribuição=="COM-Poisson", ], "/cor_ant_cmp.pdf", pathfig)
# gera_grafico(dgcorrP, "/cor_ant_cmp.pdf", pathfig)
View(0.2547279144)
sum(dgcorr[dgcorr$Distribuição=="COM-Poisson", ]$Valorp<.05, na.rm = T)
sum(dgcorr[dgcorr$Distribuição=="NB", ]$Valorp<.05, na.rm = T)
sum(dgcorr[dgcorr$Distribuição=="Poisson", ]$Valorp<.05, na.rm = T)



# Para artigo
# Agora operação para colocar em formato de matrix
df_to_mat <- function(dd){
  mat_est <- matrix(numeric(), nrow = 41, ncol = 41, dimnames = list(1:41, 1:41))
  mat_se <- matrix(numeric(), nrow = 41, ncol = 41, dimnames = list(1:41, 1:41))
  mat_inf <- matrix(numeric(), nrow = 41, ncol = 41, dimnames = list(1:41, 1:41))
  mat_sup <- matrix(numeric(), nrow = 41, ncol = 41, dimnames = list(1:41, 1:41))
  mat_p <- matrix(numeric(), nrow = 41, ncol = 41, dimnames = list(1:41, 1:41))
  for (i in 1:41){
    for (j in i:41){
      if (i == j){
        mat_est[i,j] <- 1
        mat_se[i,j] <- 0
        mat_p[i, j] <- 0
      } else{
        dgcorrfilt <- dd[dd$Var1 %in% i & dd$Var2 %in% j, , drop = F]
        mat_est[i, j] <- dgcorrfilt[, 2]
        mat_est[j, i] <- dgcorrfilt[, 2]
        
        mat_se[i, j] <- dgcorrfilt[, 3]
        mat_inf[i, j] <- as.numeric(ifelse(dgcorrfilt[, 10]>1, 1, ifelse(dgcorrfilt[, 10]<(-1),-1, dgcorrfilt[, 10])))
        mat_sup[i, j] <- as.numeric(ifelse(dgcorrfilt[, 11]>1, 1, ifelse(dgcorrfilt[, 11]<(-1),-1, dgcorrfilt[, 11])))
        
        mat_p[i, j] <- 1-dgcorrfilt[, 12]
        mat_p[j, i] <- 1-dgcorrfilt[, 12]
      }
    }
  }
  ll <- list(mat_est, mat_se, mat_inf, mat_sup, mat_p)
  names(ll) <- c("est", "se", "inf", "sup", "p")
  return(ll)
}
dgcorrP <- dgcorr[dgcorr$Distribuição=="COM-Poisson", ]
dgcorrP$Valorp[is.na(dgcorrP$Valorp)] <- 1
corrp <- df_to_mat(dgcorrP)
pathfig = "/home/guilherme/Dropbox/SuperDisperso_AJS/Supplementary material for AJS/"
namefile = "cor_ant_cmp.pdf"
pdf(file = paste0(pathfig, namefile), width = 6, height = 5.5)
graph1 <- corrplot(corrp$est, 
         mar = c(0,0,0,0),
         tl.cex = 0.6,
         cl.cex = .7,
         pch = "*",
         type = "upper",
         order = "hclust",
         diag = F,
         outline = F,
         pch.cex = 1.7,       # Controls the size of the star
         p.mat = corrp$p, #Opposite
         sig.level = .95)
# print(graph1)
graph1$corr
dev.off()
