# Package --------------------------------------------------------------------
rm(list = ls())
library(TMB)
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "NHANES_DATA/Full/"))
source(paste0(mainpath, "unify_packages.R"))
# load("negbin.RData")
gc()
#####################
# Sample port------------------------------------------------------------------------
#####################
# model <- "compoisson_multi" #1st try
# compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
# dyn.load(dynlib(paste0(mainpath, model)))
# #### Sample
# threads <- 12
# TMB::openmp(n=threads)
# parameters <- list(beta = estimates_negbin$beta,
#                    U = U_am,
#                    rho = estimates_negbin$rho,
#                    sigma = estimates_negbin$sigma,
#                    nu = log(rep(1, nr)))
# obj <- MakeADFun(data = data_am,
#                  parameters = parameters,
#                  DLL = model,
#                  hessian = T,
#                  silent = F,
#                  random = "U")
# b <- system.time(fit_compoisson_s <- nlminb(obj$par, obj$fn, obj$gr,
#                                            control = list(eval.max = 1e7, iter.max = 1e7,
#                                                           abs.tol = 1e-04, rel.tol = 1e-04)))
# # Rapidão mesmo.
# save.image("comp_sample.RData")
# load("comp_sample.RData")
# estimates_comp_am <- sv(fit = fit_compoisson_s, time = b, n_betas = n_params,
# nr = nr, model = model, method = "nlimnb", threads = threads)
# dyn.load(dynlib(paste0(mainpath, model)))
#####################
# Sample bfgs (too much memory consuming)------------------------------------------------------------------------
#####################
# threads <- 2
# TMB::openmp(n=threads)
# parameters <- list(beta = estimates_comp_am$beta,
#                    U = U_am,
#                    rho = estimates_comp_am$rho,
#                    sigma = estimates_comp_am$sigma,
#                    nu = estimates_comp_am$disp)
# obj <- MakeADFun(data = data_am,
#                  parameters = parameters,
#                  DLL = model,
#                  hessian = T,
#                  silent = F,
#                  random = "U")
# gc()
# b <- system.time(fit_compoisson_s2 <- optim(obj$par, obj$fn, obj$gr,
#                                             method = "BFGS",
#                                             control = list(maxit = 1e8,
#                                             abstol = 1e-04, reltol = 1e-04)))
# # Rapidão mesmo.
# estimates_comp_am <- sv(fit = fit_compoisson_s2, time = b, n_betas = n_params,
#                         nr = nr, model = model, method = "optim", threads = threads)
# save.image("comp_sample.RData")

#####################
# Full - PORT ------------------------------------------------------------------------
#####################
# load("comp_sample.RData")
# threads <- 8
# openmp(n=threads)
# parameters <- list(beta = estimates_comp_am$beta,
#                    U = U,
#                    rho = estimates_comp_am$rho,
#                    sigma = estimates_comp_am$sigma,
#                    nu = estimates_comp_am$disp)
# dyn.load(dynlib(paste0(mainpath, model)))
# gc()
# # First Round
# obj <- MakeADFun(data = data,
#                  parameters = parameters,
#                  DLL = model,
#                  hessian = T,
#                  silent = F
#                  ,random = c("U"))
# gc()
# a <- system.time(fit_compoisson <- nlminb(obj$par, obj$fn, obj$gr,
#                                    control = list(eval.max = 1e8, iter.max = 1e8,
#                                                   abs.tol = 1e-04, rel.tol = 1e-04)))
# gc()
# estimates_comp_s <- sv(fit = fit_compoisson, time = a, n_betas = n_params, nr = nr, model = model, method = "nlimnb", threads = threads)
# # From std
# estimates_comp_am_2 <- sv_sd(obj = obj, fit = fit_compoisson, time = a, nr = nr, model = model, method = "nlimnb", threads = threads)
# save.image("comp.RData")
# load("comp.RData")
#####################
# Full - BFGS ------------------------------------------------------------------------
#####################
# threads <- 6
# openmp(n=threads)
# parameters <- list(beta = estimates_comp_s$beta,
#                    U = U,
#                    rho = estimates_comp_s$rho,
#                    sigma = estimates_comp_s$sigma,
#                    nu = estimates_comp_s$disp)
# dyn.load(dynlib(paste0(mainpath, model)))
# gc()
# obj <- MakeADFun(data = data,
#                  parameters = parameters,
#                  DLL = model,
#                  hessian = T,
#                  silent = F
#                  ,random = c("U"))
# gc()
# a2 <- system.time(fit_compoisson2 <- optim(obj$par, obj$fn, obj$gr,
#                                     method = "BFGS",
#                                     control = list(maxit = 1e8,
#                                                    abstol = 1e-04, reltol = 1e-04)))
# gc()
# estimates_comp <- sv_sd(obj = obj, fit = fit_compoisson2, time = a2, nr = nr, model = model, method = "bfgs", threads = threads)
# save.image("comp.RData")
# load("comp.RData")
#####################
# Full - PORT ------------------------------------------------------------------------
#####################
# threads <- 5
# openmp(n=threads)
# parameters <- list(beta = estimates_comp$beta,
#                    U = U,
#                    rho = estimates_comp$rho,
#                    sigma = estimates_comp$sigma,
#                    nu = estimates_comp$disp)
# dyn.load(dynlib(paste0(mainpath, model)))
# gc()
# # First Round
# obj <- MakeADFun(data = data,
#                  parameters = parameters,
#                  DLL = model,
#                  hessian = T,
#                  silent = T
#                  ,random = c("U"))
# gc()
# a <- system.time(fit_compoisson <- nlminb(obj$par, obj$fn, obj$gr,
#                                    control = list(eval.max = 1e9, iter.max = 1e9,
#                                                   abs.tol = 1e-04, rel.tol = 1e-04)))
# gc()
# # estimates_comp é o escolhido
# estimates_comp <- sv(fit = fit_compoisson, time = a, n_betas = n_params, nr = nr, 
#                      model = model, method = "nlimnb", threads = threads)
# # From std
# estimates_comp_sd <- sv_sd(obj = obj, fit = fit_compoisson, time = a, nr = nr, 
#                            model = model, method = "nlimnb", threads = threads)
# save.image("comp.RData")
load("comp.RData")

#####################
# Resultado do modelo retirando os betas com NaN no std error (COMPOISSON) ------------------------------------------------------------------------
#####################
model <- "compoisson_tri_fixed_effects" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
threads <- 12
openmp(n=threads)
parameters <- list(beta1 = ests_new$beta[c(1,2,4),1],
                   beta2 = ests_new$beta[c(1,2,4,5),2],
                   beta3 = ests_new$beta[c(1,2,3,4),3],
                   U = U,
                   rho = ests_new$rho,
                   sigma = ests_new$sigma,
                   nu = ests_new$disp)
f1 <- as.formula("Nmsp ~ Race + Marital")
f2 <- as.formula("Nmsp ~ Race + Marital + Age")
f3 <- as.formula("Nmsp ~ Race + Education + Marital")
cov1 <- c("Intercepto", strsplit(gsub(".*~ ", "", format(f1)), split = " \\+ ")[[1]])
cov2 <- c("Intercepto", strsplit(gsub(".*~ ", "", format(f2)), split = " \\+ ")[[1]])
cov3 <- c("Intercepto", strsplit(gsub(".*~ ", "", format(f3)), split = " \\+ ")[[1]])

data <- list(Y = as.matrix(nhanes[, vars_resp]),
             X1 = model.matrix(f1, nhanes),
             X2 = model.matrix(f2, nhanes),
             X3 = model.matrix(f3, nhanes))
dyn.load(dynlib(paste0(mainpath, model)))
gc()
# First Round
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
gc()
a <- system.time(fit_compoisson <- nlminb(obj$par, obj$fn, obj$gr,
                                          control = list(eval.max = 1e9, iter.max = 1e9,
                                                         abs.tol = 1e-04, rel.tol = 1e-04)))
save.image("comp.RData")
load("comp.RData")
rep <- sdreport(obj)
# est_com_res <- sv(fit = fit_compoisson, time = a, n_betas = 11, nr = nr, 
                     # model = model, method = "nlimnb", threads = threads)
# From std
est_com_res_sd <- sv_sd(obj = obj, fit = fit_compoisson, time = a, nr = nr, 
                        model = model, method = "nlimnb", threads = threads,
                        namebeta = c(cov1, cov2, cov3), 
                        lengthresp = c(rep(1, length(cov1)), rep(2, length(cov2)), rep(3, length(cov3))))
save.image("comp.RData")
#####################
### Retirando aqueles com NaN
#####################
model <- "compoisson_tri_fixed_effects" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
threads <- 12
openmp(n=threads)
parameters <- list(beta1 = as.numeric(est_com_res_sd$ests[c(1,2,3),1]),
                   beta2 = as.numeric(est_com_res_sd$ests[c(4,5,6),1]),
                   beta3 = as.numeric(est_com_res_sd$ests[c(8,9,10,11),1]),
                   U = U,
                   rho = est_com_res_sd$rho,
                   sigma = est_com_res_sd$sigma,
                   nu = est_com_res_sd$disp)
f1 <- as.formula("Nmsp ~ Race + Marital")
f2 <- as.formula("Nmsp ~ Race + Marital")
f3 <- as.formula("Nmsp ~ Race + Education + Marital")
cov1 <- c("Intercepto", strsplit(gsub(".*~ ", "", format(f1)), split = " \\+ ")[[1]])
cov2 <- c("Intercepto", strsplit(gsub(".*~ ", "", format(f2)), split = " \\+ ")[[1]])
cov3 <- c("Intercepto", strsplit(gsub(".*~ ", "", format(f3)), split = " \\+ ")[[1]])
data <- list(Y = as.matrix(nhanes[, vars_resp]),
             X1 = model.matrix(f1, nhanes),
             X2 = model.matrix(f2, nhanes),
             X3 = model.matrix(f3, nhanes))
dyn.load(dynlib(paste0(mainpath, model)))
gc()
# data$X1%*%parameters$beta1
# data$X2%*%parameters$beta2
# data$X3%*%parameters$beta3
# First Round
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
gc()
FreeADFun(obj)
obj$fn()
a <- system.time(fit_compoisson_f <- nlminb(obj$par, obj$fn, obj$gr,
                                          control = list(eval.max = 1e9, iter.max = 1e9,
                                                         abs.tol = 1e-04, rel.tol = 1e-04)))
rep <- sdreport(obj)
# From std
# est_com_res_sd$summary
# fit_compoisson_f$objective
est_com_res_sd <- sv_sd(obj = obj, fit = fit_compoisson_f, time = a, nr = nr, 
                        model = model, method = "nlimnb", threads = threads,
                        namebeta = c(cov1, cov2, cov3), 
                        lengthresp = c(rep(1, length(cov1)), rep(2, length(cov2)), rep(3, length(cov3))))
save.image("comp.RData")
#####################
# NEGBIN ------------------------------------------------------------------------
#####################
model <- "negbin_tri_fixed_effects" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
parameters <- list(beta1 = as.numeric(est_com_res_sd$ests[c(1,2,3),1]),
                   beta2 = as.numeric(est_com_res_sd$ests[c(4,5,6),1]),
                   beta3 = as.numeric(est_com_res_sd$ests[c(7,8,9,10),1]),
                   U = U,
                   rho = est_com_res_sd$rho,
                   sigma = est_com_res_sd$sigma,
                   phi = rep(log(1), nr))
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
gc()
openmp(12)
a <- system.time(negbin_f <- nlminb(obj$par, obj$fn, obj$gr,
                                    control = list(eval.max = 1e9, iter.max = 1e9,
                                                   abs.tol = 1e-04, rel.tol = 1e-04)))
est_neg_sd <- sv_sd(obj = obj, fit = negbin_f, time = a, nr = nr, 
                        model = model, method = "nlimnb", threads = threads,
                        namebeta = c(cov1, cov2, cov3), 
                        lengthresp = c(rep(1, length(cov1)), rep(2, length(cov2)), rep(3, length(cov3))))
gc()
#####################
# POISSON ------------------------------------------------------------------------
#####################
model <- "02_poisson_tri_fixed_effects" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
parameters <- list(beta1 = as.numeric(est_neg_sd$ests[c(1,2,3),1]),
                   beta2 = as.numeric(est_neg_sd$ests[c(4,5,6),1]),
                   beta3 = as.numeric(est_neg_sd$ests[c(7,8,9,10),1]),
                   U = U,
                   rho = rep(0, nr),
                   sigma = est_com_res_sd$sigma)
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
gc()
openmp(12)
a <- system.time(poisson_f <- nlminb(obj$par, obj$fn, obj$gr,
                                    control = list(eval.max = 1e9, iter.max = 1e9,
                                                   abs.tol = 1e-04, rel.tol = 1e-04)))
est_poisson_sd <- sv_sd(obj = obj, fit = poisson_f, time = a, nr = nr, 
                    model = model, method = "nlimnb", threads = threads,
                    namebeta = c(cov1, cov2, cov3), 
                    lengthresp = c(rep(1, length(cov1)), rep(2, length(cov2)), rep(3, length(cov3))))
gc()
save.image("comp.RData")

#################################################################
# Likelihood summary for all
#################################################################
load("comp.RData")
est_poisson_sd$summary$Model <- "Poisson"
est_poisson_sd$summary$np <- nrow(est_poisson_sd$ests)
est_poisson_sd$summary$SE <- TRUE
m1 <- est_poisson_sd$summary
est_neg_sd$summary$Model <- "NB"
est_neg_sd$summary$np <- nrow(est_neg_sd$ests)
est_neg_sd$summary$SE <- TRUE
m2 <- est_neg_sd$summary
est_com_res_sd$summary$Model <- "COM-Poisson Full"
est_com_res_sd$summary$np <- nrow(est_com_res_sd$ests)
est_com_res_sd$summary$SE <- TRUE
m3 <- est_com_res_sd$summary
mm <- rbind(m1, m2, m3)
mm$Loglik <- mm$Loglik*-1
mm$AIC <- -2*mm$Loglik + 2*mm$np
mm$BIC <- -2*mm$Loglik + log(1281)*mm$np
mm <- mm[, c(6,7,9,10,1,8,5)]
row.names(mm) <- NULL
mm$Optimizer <- with(mm, ifelse(Optimizer=="bfgs", "BFGS", 
                                ifelse(Optimizer%in%c("nlminb","nlimnb"), "PORT", Optimizer)))
mm <- mm[, -c(7)]
mm$Model <- gsub(" Full", "", mm$Model)
mm$SE <- "\\checkmark"
library(knitr)
library(ggplot2)
library(latex2exp)
library(dplyr)
kable(mm, 
      caption = "Model fit measures for NHANES data from the best parametrization for each distribution",
      label = "nhanesfit2",
      align = c("l", rep("c", ncol(mm)-1)),
      booktabs = T,
      digits = 2,
      row.names = F,
      format = "latex",
      escape = F)
#################################################################
############# Graphic for all 
#################################################################
np_poisson <- nrow(est_poisson_sd$ests)
np_cmp <- nrow(est_com_res_sd$ests)
dg <- data.frame(Parametros = c(row.names(est_poisson_sd$ests), row.names(est_neg_sd$ests), row.names(est_com_res_sd$ests)),
                 Estimativas = c(est_poisson_sd$ests[,1], est_neg_sd$ests[,1], est_com_res_sd$ests[,1]),
                 ErroPadrao = c(est_poisson_sd$ests[,2], est_neg_sd$ests[,2], est_com_res_sd$ests[,2]),
                 Distribuição = c(rep("Poisson", np_poisson), rep("NB", np_cmp), rep("COM-Poisson", np_cmp)),
                 Tipo = c(c(rep("Média", 10), rep("Variância", 3), rep("Correlação", 3)), #Poisson
                          c(rep("Média", 10), rep("Variância", 3), rep("Correlação", 3), rep("Dispersão", 3)), #NB
                          c(rep("Média", 10), rep("Variância", 3), rep("Correlação", 3), rep("Dispersão", 3)))) #COM
# head(dg,20)
dg[c(33:35, 52:54), 1] <- c(paste0("exp(phi)", 1:3), paste0("exp(nu)", 1:3))
save.image("comp.RData")
dg$Response <- factor(substr(dg$Parametros, nchar(dg$Parametros), nchar(dg$Parametros)),
                      labels = c("Nmsp", "Nmosp", "Nspfy"))
load("comp.RData")
library(latex2exp)
library(ggplot2)
library(tidyr)
library(dplyr)
## Only beta (exp() or only beta()?)
dg %>% 
  filter(Tipo=="Média") %>% 
  mutate(Parametros = factor(gsub("_Y\\d", "", Parametros), 
                             levels = c("Intercepto", "Race", "Marital", "Education"),
                             labels = c("Intercept", "Race", "Marital", "Education"))) %>% 
  ggplot() +
  # ggplot(aes(x = Parametros, y = Estimativas, shape = Distribuição)) +
  geom_pointrange(aes(x = Parametros, y = Estimativas, shape = Distribuição,
                      ymin = Estimativas-1.96*ErroPadrao,
                      ymax = Estimativas+1.96*ErroPadrao),
                  position = position_dodge(width = 0.8)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        text = element_text(size = 12),
        legend.position = "top",
        legend.title = element_blank()) +
  facet_wrap(~Response, scales = "free_x") +
  coord_flip() +
  labs(y = TeX('$\\hat{\\beta}\\pm 1.96 SE}$'),
       linetype = "Model") +
  guides(ymin = guide_legend(override.aes = list(ymin = 0, ymax = 0)))
pathfig <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/Texto/LaTeX/Figuras"
# ggsave(paste0(pathfig, "/nhanes_beta_final_model.pdf"), units = "cm")
ggsave(paste0(pathfig, "/nhanes_beta_final_model.pdf"), width = 18, height = 8, units = "cm")

# Except Dispersion
myHeader <- c(1,2,2)

names(myHeader) <- c(" ", "Negative binomial ($\\\\phi$)", "COM-Poisson ($\\\\nu$)")
dg %>% 
  filter(Tipo == "Dispersão") %>% 
  select(-Tipo, -Response, -Parametros) %>% 
  mutate(Response = rep(c("Nmsp", "Nmosp", "Nspfy"), times = 2)) %>% 
  pivot_wider(names_from = Distribuição,
              values_from = c("Estimativas", "ErroPadrao")) %>% 
  select(Response, Estimativas_NB, ErroPadrao_NB, everything()) %>% 
  kable(caption = toupper("Dispersion parameter estimates and standard errors (SE) for each model and outcome of NHANES data"),
        label = "nhanesdisp",
        booktabs = T,
        digits = 3,
        row.names = F,
        escape = F,
        format = "latex",
        col.names = c("Outcome","Estimate", "SE", "Estimate", "SE"),
        align = "c") %>% 
  add_header_above(header = myHeader,
                   escape = F)
 

est_poisson_sd$corr_matrix
est_poisson_sd$ests

est_neg_sd$corr_matrix
est_neg_sd$ests

est_com_res_sd$corr_matrix
est_com_res_sd$ests
