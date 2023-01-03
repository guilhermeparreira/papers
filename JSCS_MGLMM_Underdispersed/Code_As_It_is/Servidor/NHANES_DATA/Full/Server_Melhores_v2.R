# Package --------------------------------------------------------------------
rm(list = ls())
library(TMB)
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "NHANES_DATA/Full/"))
source(paste0(mainpath, "unify_packages.R"))
load("double_v2.RData")

#####################
# Resultado do modelo retirando os betas com NaN no std error (COMPOISSON) ------------------------------------------------------------------------
#####################
model <- "doublepoisson_tri_fixed_effects" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
threads <- 12
openmp(n=threads)
parameters <- list(beta1 = estimates_double_5$beta[c(1,2,3), 1],
                   beta2 = estimates_double_5$beta[, 2],
                   beta3 = estimates_double_5$beta[, 3],
                   U = U,
                   rho = estimates_double_5$rho,
                   sigma = estimates_double_5$sigma,
                   nu = estimates_double_5$disp)
f1 <- as.formula("Nmsp ~ Race + Education")
f2 <- as.formula("Nmsp ~ Race + Education + Marital + Age")
f3 <- as.formula("Nmsp ~ Race + Education + Marital + Age")
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
a <- system.time(fit_double_v9 <- nlminb(obj$par, obj$fn, obj$gr,
                                          control = list(eval.max = 1e9, iter.max = 1e9,
                                                         abs.tol = 1e-04, rel.tol = 1e-04)))
save.image("double_v2.RData")
estimates_double_6 <- sv(fit = fit_double_v9, time = a, n_betas = length(c(cov1, cov2, cov3)), nr = nr,
                            model = model, method = "nlimnb", threads = threads)
# From std
estimates_double_6_sd <- sv_sd(obj = obj, fit = fit_double_v9, time = a, nr = nr, 
                        model = model, method = "nlimnb", threads = threads,
                        namebeta = c(cov1, cov2, cov3), 
                        lengthresp = c(rep(1, length(cov1)), rep(2, length(cov2)), rep(3, length(cov3))))

# Second
parameters <- list(beta1 = estimates_double_6_sd$ests[1:3, 1],
                   beta2 = estimates_double_6_sd$ests[4:8, 1],
                   beta3 = estimates_double_6_sd$ests[9:13, 1],
                   U = U,
                   rho = estimates_double_6$rho,
                   sigma = estimates_double_6$sigma,
                   nu = estimates_double_6$disp)
dyn.load(dynlib(paste0(mainpath, model)))
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
gc()
obj$fn()
b <- system.time(fit_double_v10 <- optim(obj$par, obj$fn, obj$gr,
                                          method = 'BFGS'))

estimates_double_7 <- sv(fit = fit_double_v10, time = a, 
                         n_betas = length(c(cov1, cov2, cov3)), 
                         nr = nr, model = model, method = "nlimnb", threads = threads)

# From std
estimates_double_7_sd <- sv_sd(obj = obj, fit = fit_double_v10, time = a, nr = nr, 
                               model = model, method = "nlimnb", threads = threads,
                               namebeta = c(cov1, cov2, cov3), 
                               lengthresp = c(rep(1, length(cov1)), rep(2, length(cov2)), rep(3, length(cov3))))
save.image("double_v2.RData")

#####################
### COM-Poisson - Retirando aqueles com NaN
#####################
model <- "compoisson_tri_fixed_effects" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))

# Com o double poisson deu ruim. peguei o full do compoisson
parameters <- list(beta1 = c(-0.212830376847286, -0.000423444838758394, 0.000857163844311496
), beta2 = c(-1.41616295420935, 0.00425046307053224, 0.0185903326750339, 
             -0.00921678961348591, -0.00366189767118222), beta3 = c(-1.47120174715083, 
                                                                    0.0221582281718454, -0.0223736760762393, -0.0115880253266172, 
                                                                    0.00147821118325441), rho = c(r_Y1_Y2 = -0.00340009440463812, 
                                                                                                  r_Y1_Y3 = -0.0410213503237788, r_Y2_Y3 = 0.0150502411969709), 
sigma = c(sig1 = 0.459751305287482, sig2 = 1.18117844136332, 
          sig3 = 1.17486153830211), nu = c(`exp(nu)` = 3.40289876412992, 
                                           `exp(nu)` = 2.92925827119835, `exp(nu)` = 3.02669988235857
          ))
parameters$U <- U

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
obj$fn()
a <- system.time(fit_compoisson_f2 <- nlminb(obj$par, obj$fn, obj$gr,
                                          control = list(eval.max = 1e9, iter.max = 1e9,
                                                         abs.tol = 1e-04, rel.tol = 1e-04)))
# From std
# est_com_res_sd$summary
# fit_compoisson_f$objective
est_com_res <- sv(fit = fit_compoisson_f2, time = a, 
                  n_betas = length(c(cov1, cov2, cov3)), 
                  nr = nr, model = model, method = "nlimnb", threads = threads)

est_com_res_sd <- sv_sd(obj = obj, fit = fit_compoisson_f2, time = a, nr = nr, 
                        model = model, method = "nlimnb", threads = threads,
                        namebeta = c(cov1, cov2, cov3), 
                        lengthresp = c(rep(1, length(cov1)), rep(2, length(cov2)), rep(3, length(cov3))))
save.image("double_v2.RData")


#####################
### COM-Poisson 2
#####################
parameters <- list(beta1 = est_com_res_sd$ests[1:3, 1],
                   beta2 = est_com_res_sd$ests[4:8, 1],
                   beta3 = est_com_res_sd$ests[9:13, 1],
                   U = U,
                   rho = est_com_res$rho,
                   sigma = est_com_res$sigma,
                   nu = est_com_res$disp)

dyn.load(dynlib(paste0(mainpath, model)))
threads <- 12
openmp(n=threads)
gc()

# First Round
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
gc()
obj$fn()
a <- system.time(fit_compoisson_f3 <- nlminb(obj$par, obj$fn, obj$gr,
                                             control = list(eval.max = 1e9, iter.max = 1e9,
                                                            abs.tol = 1e-04, rel.tol = 1e-04)))
# From std
# est_com_res_sd$summary
# fit_compoisson_f$objective
est_com_res_sd_2 <- sv_sd(obj = obj, fit = fit_compoisson_f3, time = a, nr = nr, 
                        model = model, method = "nlimnb", threads = threads,
                        namebeta = c(cov1, cov2, cov3), 
                        lengthresp = c(rep(1, length(cov1)), rep(2, length(cov2)), rep(3, length(cov3))))
est_com_res_2 <- sv(fit = fit_compoisson_f3, time = a, 
                  n_betas = length(c(cov1, cov2, cov3)), 
                  nr = nr, model = model, method = "nlimnb", threads = threads)

save.image("double_v2.RData")

#####################
### COM-Poisson 3 (Didn't work)
#####################
# parameters <- list(beta1 = est_com_res_sd_2$ests[1:3, 1],
#                    beta2 = est_com_res_sd_2$ests[4:8, 1],
#                    beta3 = est_com_res_sd_2$ests[9:13, 1],
#                    U = U,
#                    rho = est_com_res_2$rho,
#                    sigma = est_com_res_2$sigma,
#                    nu = est_com_res_2$disp)
# 
# dyn.load(dynlib(paste0(mainpath, model)))
# threads <- 4
# openmp(n=threads)
# gc()
# 
# # First Round
# obj <- MakeADFun(data = data,
#                  parameters = parameters,
#                  DLL = model,
#                  hessian = T,
#                  silent = T
#                  ,random = c("U"))
# gc()
# obj$fn()
# b <- system.time(fit_compoisson_f4 <- optim(obj$par, obj$fn, obj$gr,
#                                          method = 'BFGS'))
# # From std
# # est_com_res_sd$summary
# # fit_compoisson_f$objective
# est_com_res_sd_3 <- sv_sd(obj = obj, fit = fit_compoisson_f4, time = a, nr = nr, 
#                         model = model, method = "nlimnb", threads = threads,
#                         namebeta = c(cov1, cov2, cov3), 
#                         lengthresp = c(rep(1, length(cov1)), rep(2, length(cov2)), rep(3, length(cov3))))
# est_com_res_3 <- sv(fit = fit_compoisson_f3, time = a, 
#                     n_betas = length(c(cov1, cov2, cov3)), 
#                     nr = nr, model = model, method = "nlimnb", threads = threads)
# 
# save.image("double_v2.RData")




#####################
# NEGBIN ------------------------------------------------------------------------
#####################
model <- "negbin_tri_fixed_effects" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
parameters <- list(beta1 = est_com_res_sd_2$ests[1:3, 1],
                   beta2 = est_com_res_sd_2$ests[4:8, 1],
                   beta3 = est_com_res_sd_2$ests[9:13, 1],
                   U = U,
                   rho = est_com_res_2$rho,
                   sigma = est_com_res_2$sigma,
                   phi = rep(log(300), nr)) # I know this from previou code
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL = model,
                 hessian = T,
                 silent = T
                 ,random = c("U"))
obj$fn()
gc()
openmp(12)
a <- system.time(negbin_f <- nlminb(obj$par, obj$fn, obj$gr,
                                    control = list(eval.max = 1e9, iter.max = 1e9,
                                                   abs.tol = 1e-04, rel.tol = 1e-04)))
est_neg_sd <- sv_sd(obj = obj, fit = negbin_f, time = a, nr = nr, 
                        model = model, method = "nlimnb", threads = threads,
                        namebeta = c(cov1, cov2, cov3), 
                        lengthresp = c(rep(1, length(cov1)), rep(2, length(cov2)), rep(3, length(cov3))))
est_neg_1 <- sv(fit = negbin_f, time = a,
                    n_betas = length(c(cov1, cov2, cov3)),
                    nr = nr, model = model, method = "nlimnb", threads = threads)
gc()
#####################
# POISSON ------------------------------------------------------------------------
#####################
model <- "02_poisson_tri_fixed_effects" #1st try
compile(paste0(mainpath, model, ".cpp"), flags = "-O0 -ggdb")
dyn.load(dynlib(paste0(mainpath, model)))
parameters <- list(beta1 = as.numeric(est_neg_sd$ests[c(1:3),1]),
                   beta2 = as.numeric(est_neg_sd$ests[c(4:8),1]),
                   beta3 = as.numeric(est_neg_sd$ests[c(9:13),1]),
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
openmp(2, DLL = model)
a <- system.time(poisson_f <- nlminb(obj$par, obj$fn, obj$gr,
                                    control = list(eval.max = 1e9, iter.max = 1e9,
                                                   abs.tol = 1e-04, rel.tol = 1e-04)))
est_poisson_sd <- sv_sd(obj = obj, fit = poisson_f, time = a, nr = nr, 
                    model = model, method = "nlimnb", threads = threads,
                    namebeta = c(cov1, cov2, cov3), 
                    lengthresp = c(rep(1, length(cov1)), rep(2, length(cov2)), rep(3, length(cov3))))
est_poisson_1 <- sv(fit = poisson_f, time = a,
                n_betas = length(c(cov1, cov2, cov3)),
                nr = nr, model = model, method = "nlimnb", threads = threads)
gc()
save.image("double_v2.RData")
#################################################################
# Likelihood summary for all
#################################################################
load("double_v2.RData")
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
estimates_double_7_sd$summary$Model <- "Double Poisson"
estimates_double_7_sd$summary$np <- nrow(estimates_double_7_sd$ests)
estimates_double_7_sd$summary$SE <- TRUE
m4 <- estimates_double_7_sd$summary

mm <- rbind(m1, m2, m3, m4)
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
      caption = "Goodness-of-fit measures for NHANES data from the best parametrization for each distribution",
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
np_cmp <- nrow(est_com_res_sd$ests)
np_double <- nrow(estimates_double_7_sd$ests)

dg <- data.frame(Parametros = c(row.names(est_poisson_sd$ests), row.names(est_neg_sd$ests),
                                row.names(est_com_res_sd$ests), row.names(estimates_double_7_sd$ests)),
                 Estimativas = c(est_poisson_sd$ests[,1], est_neg_sd$ests[,1], est_com_res_sd$ests[,1], estimates_double_7_sd$ests[,1]),
                 ErroPadrao = c(est_poisson_sd$ests[,2], est_neg_sd$ests[,2], est_com_res_sd$ests[,2], estimates_double_7_sd$ests[,2]),
                 Valorp = c(est_poisson_sd$ests[,4], est_neg_sd$ests[,4], est_com_res_sd$ests[,4], estimates_double_7_sd$ests[,4]),
                 Distribuição = c(rep("Poisson", np_poisson), rep("NB", np_cmp), rep("COM-Poisson", np_cmp), rep("Double Poisson", np_cmp)),
                 Tipo = c(c(rep("Média", 13), rep("Variância", 3), rep("Correlação", 3)), #Poisson
                          c(rep("Média", 13), rep("Variância", 3), rep("Correlação", 3), rep("Dispersão", 3)), #NB
                          c(rep("Média", 13), rep("Variância", 3), rep("Correlação", 3), rep("Dispersão", 3)), #COM
                          c(rep("Média", 13), rep("Variância", 3), rep("Correlação", 3), rep("Dispersão", 3))  #Double
                        )
                 ) 
# head(dg,20)
dg[c(39:41, 61:63, 83:85), 1] <- c(paste0("exp(phi)", 1:3), paste0("exp(nu)", 1:3), paste0("exp(nu)", 1:3))
save.image("double_v2.RData")
dg$Response <- factor(substr(dg$Parametros, nchar(dg$Parametros), nchar(dg$Parametros)),
                      labels = c("Nmsp", "Nmosp", "Nspfy"))
load("double_v2.RData")
library(latex2exp)
library(ggplot2)
library(tidyr)
library(dplyr)
## Only beta (exp() or only beta()?)
dg$NAN <- factor(as.character(ifelse(is.nan(dg$ErroPadrao), 1, 20)))
dg$Distribuição_NAN <- factor(ifelse(is.nan(dg$ErroPadrao), paste0(dg$Distribuição, ' - Without SE'), dg$Distribuição))

dg %>% 
  filter(Tipo=="Média") %>% 
  mutate(Parametros = factor(gsub("_Y\\d", "", Parametros))) %>% 
  ggplot() +
  # ggplot(aes(x = Parametros, y = Estimativas, shape = Distribuição)) +
  geom_pointrange(aes(x = Parametros, 
                      y = Estimativas, 
                      shape = Distribuição_NAN,
                      ymin = Estimativas-1.96*ErroPadrao,
                      ymax = Estimativas+1.96*ErroPadrao),
                  size = .3,
                  linetype = 3,
                  position = position_dodge(width = 0.8)) +
  scale_shape_manual(values = c(20,1,2,3,4)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        text = element_text(size = 10),
        legend.position = "top",
        legend.title = element_blank()) +
  facet_wrap(~Response, scales = "free_x") +
  coord_flip() +
  labs(y = TeX('$\\hat{\\beta}\\pm 1.96 SE$'),
       linetype = "Model") +
  guides(ymin = guide_legend(override.aes = list(ymin = 0, ymax = 0)))

# pathfig <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/Texto/LaTeX/Figuras"
pathfig <- "/home/guilherme/Dropbox/Underdispersed_Count"
# ggsave(paste0(pathfig, "/nhanes_beta_final_model.pdf"), units = "cm")
ggsave(paste0(pathfig, "/nhanes_beta_final_model.pdf"), width = 18, height = 11, units = "cm")

# Except Dispersion
myHeader <- c(1,2,2, 2)
library(kableExtra)
names(myHeader) <- c(" ", "Negative binomial ($\hat\\\\phi$)", "COM-Poisson ($\\\\nu$)", "Double Poisson ($\\\\nu$)")

dg %>% 
  filter(Tipo == "Dispersão") %>% 
  select(-Tipo, -Response, -Parametros, -NAN, -Distribuição_NAN) %>% 
  mutate(Response = rep(c("Nmsp", "Nmosp", "Nspfy"), times = 3)) %>% 
  pivot_wider(names_from = Distribuição,
              values_from = c("Estimativas", "ErroPadrao")) %>% 
  select(Response, Estimativas_NB, ErroPadrao_NB, `Estimativas_COM-Poisson`, `ErroPadrao_COM-Poisson`, everything()) %>% 
  kable(caption = toupper("Dispersion parameter estimates and standard errors (SE) for each model and outcome of NHANES data"),
        label = "nhanesdisp",
        booktabs = T,
        digits = 3,
        row.names = F,
        escape = F,
        format = "latex",
        col.names = c("Outcome",rep(c("Estimate", "SE"), 3)),
        align = "c") %>% 
  add_header_above(header = myHeader,
                   escape = F)
 

est_poisson_sd$corr_matrix
est_poisson_sd$ests

est_neg_sd$corr_matrix
est_neg_sd$ests

est_com_res_sd$corr_matrix
est_com_res_sd$ests

estimates_double_7_sd$ests
estimates_double_7_sd$corr_matrix

dg %>% filter(Tipo=="Média") %>% filter(str_detect(Parametros, "Race")) %>% arrange(Response)
dg %>% filter(Tipo=="Média") %>% filter(str_detect(Parametros, "Marital")) %>% arrange(Response)
dg %>% filter(Tipo=="Média") %>% filter(str_detect(Parametros, "Education")) %>% arrange(Response)
dg %>% filter(Tipo=="Média") %>% filter(str_detect(Parametros, "Age")) %>% arrange(Response)
