rm()
path <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor"
localp <- "/NHANES_DATA/Full"
locals <- "/NHANES_DATA/SmallerModels"
setwd(paste0(path, localp))
library(reshape2)
library(latex2exp)
library(ggplot2)
library(tidyr)
library(dplyr)
library(kableExtra)
library(knitr)
options(knitr.kable.NA = '')
##########################
# Model 1 - Poisson
##########################
load("Poisson.RData")
# ls()
m01 <- as.data.frame(estimates_poisson_sd$ests)
m01$Model <- "Poisson"
m01$MainModel <- "Poisson"
rm(list = ls()[!ls()%in%c("m01", "path", "localp", "locals")])
##########################
# Model 2 - Negbin
##########################
load("negbin.RData")
m02 <- as.data.frame(estimates_negbin_sd$ests)
m02$Model <- "NB"
m02$MainModel <- "NB"
rm(list = ls()[!ls()%in%c("m01", "m02", "path", "localp", "locals")])
##########################
# Model 3 - Negbin - Disp Fixed
##########################
setwd(paste0(path, locals))
load("Negbin_phi_fixed.RData")
#################################################################################
######################### Só preciso fazer exp(disp) - tanto no NB como no COMPOISSON
#################################################################################
b1p <- melt(estimates_negbin$beta)
b1p <- b1p[,-c(1:2), drop = FALSE]
s1p <- melt(estimates_negbin$sigma)
r1p <- melt(estimates_negbin$rho)
m03 <- data.frame(Estimate = unname(rbind(b1p, s1p, r1p)),
           `Std. Error` = NA,
           Z = NA,
           `p value` = NA,
           check.names = F,
           row.names = row.names(m01))
m03$Model <- 'Fixed Disp NB'
m03$MainModel <- "NB"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "path", "localp", "locals")])
##########################
# Model 4 - Negbin - Comum Variance
##########################
load("Negbin_var_comum.RData")
b1p <- melt(estimates_negbin$beta)
b1p <- b1p[,-c(1:2), drop = FALSE]
s1p <- melt(estimates_negbin$sigma)
r1p <- melt(estimates_negbin$rho)
d1p <- exp(melt(estimates_negbin$disp)) # Transformação aqui.
m04 <- data.frame(Estimate = unname(rbind(b1p, s1p, r1p, d1p)),
                 `Std. Error` = NA,
                 Z = NA,
                 `p value` = NA,
                 check.names = F,
                 row.names = c(row.names(m01)[1:15], 
                               'sig1', 
                               row.names(m01)[19:21], 
                               c('exp(phi)1', 'exp(phi)2', 'exp(phi)3')))
m04$Model <- 'Common Var NB'
m04$MainModel <- "NB"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04" ,"path", "localp", "locals")])
##########################
# Model 5 - Negbin - Fixed Var
##########################
load("Negbin_var_fixed.RData")
m05 <- as.data.frame(estimates_negbin_sd$ests)
m05$Model <- "Fixed Var NB"
m05$MainModel <- "NB"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05" ,"path", "localp", "locals")])
##########################
# Model 6 - Compoisson FULL
##########################
setwd(paste0(path, localp))
load("comp.RData")
m06 <- as.data.frame(estimates_comp$ests)
m06$Model <- "COM-Poisson"
m06$MainModel <- "COM-Poisson"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06" ,"path", "localp", "locals")])
##########################
# Model 7 - Compoisson Disp Fixed 
##########################
setwd(paste0(path, locals))
load("compoisson_nu_fixed_1.5.RData")
# estimates_compoisson_s
b1p <- melt(estimates_compoisson$beta)
b1p <- b1p[,-c(1:2), drop = FALSE]
s1p <- melt(estimates_compoisson$sigma)
r1p <- melt(estimates_compoisson$rho)
m07 <- data.frame(Estimate = unname(rbind(b1p, s1p, r1p)),
                 `Std. Error` = NA,
                 Z = NA,
                 `p value` = NA,
                 check.names = F,
                 row.names = row.names(m03))
m07$Model <- 'Fixed Disp COM-Poisson'
m07$MainModel <- "COM-Poisson"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06", "m07" ,"path", "localp", "locals")])
##########################
# Model 8 - Compoisson Comum Var
##########################
load("compoisson_var_comum.RData")
m08 <- as.data.frame(estimates_compoisson_sd$ests)
m08$Model <- "Common Var COM-Poisson"
m08$MainModel <- "COM-Poisson"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06", "m07", "m08","path", "localp", "locals")])
##########################
# Model 9 - Compoisson Fixed Var
##########################
load("compoisson_var_fixed.RData")
m09 <- as.data.frame(estimates_compoisson_sd_2$ests)
m09$Model <- "Fixed Var COM-Poisson"
m09$MainModel <- "COM-Poisson"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06", "m07", "m08", "m09" ,"path", "localp", "locals")])
##########################
# Model 10 - Poisson Rho fixed
##########################
setwd(paste0(path, locals))
load("Poisson_rho_fixed.RData")
# estimates_compoisson_s
b1p <- melt(estimates_poisson$beta)
b1p <- b1p[,-c(1:2), drop = FALSE]
s1p <- melt(estimates_poisson$sigma)
m10 <- data.frame(Estimate = unname(rbind(b1p, s1p)),
                  `Std. Error` = NA,
                  Z = NA,
                  `p value` = NA,
                  check.names = F,
                  row.names = row.names(m03)[1:18])
m10$Model <- "Rho Zero Poisson"
m10$MainModel <- "Poisson"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06", "m07", "m08", "m09", "m10" ,"path", "localp", "locals")])
##########################
# Model 11 - NB Rho fixed
##########################
load("negbin_rho_fixed.RData")
# estimates_compoisson_s
b1p <- melt(estimates_negbin0$beta)
b1p <- b1p[,-c(1:2), drop = FALSE]
s1p <- melt(estimates_negbin0$sigma)
d1p <- exp(melt(estimates_negbin0$disp))
m11 <- data.frame(Estimate = unname(rbind(b1p, s1p, d1p)),
                  `Std. Error` = NA,
                  Z = NA,
                  `p value` = NA,
                  check.names = F,
                  row.names = c(row.names(m10)[1:18], 'exp(phi)1', 'exp(phi)2', 'exp(phi)3'))
m11$Model <- "Rho Zero NB"
m11$MainModel <- "NB"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06", "m07", "m08", "m09", "m10", "m11" ,"path", "localp", "locals")])
##########################
# Model 12 - COM-Poisson rho fixed
##########################
load("comp_rho_fixed.RData")
# estimates_compoisson_s
m12 <- as.data.frame(estimates_comp$ests)
m12 <- m12[!grepl("r_Y", rownames(m12)),]
m12$Model <- "Rho Zero COM-Poisson"
m12$MainModel <- "COM-Poisson"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06", "m07", "m08", "m09", "m10", "m11", "m12" ,"path", "localp", "locals")])
#####################################################3
####### Joining all results
#####################################################3
# sample size 1281
m01$Parametros <- row.names(m01)
m02$Parametros <- row.names(m02)
m03$Parametros <- row.names(m03)
m04$Parametros <- row.names(m04)
m05$Parametros <- row.names(m05)
m06$Parametros <- row.names(m06)
m07$Parametros <- row.names(m07)
m08$Parametros <- row.names(m08)
m09$Parametros <- row.names(m09)
m10$Parametros <- row.names(m10)
m11$Parametros <- row.names(m11)
m12$Parametros <- row.names(m12)
mm <- rbind(m01,m02,m03,m04,m05,m06, m07, m08, m09, m10, m11, m12)
row.names(mm) <- NULL
# Creating variables
mm$Tipo <- ifelse(grepl("b[0-9]", mm$Parametros), "Média",
                  ifelse(grepl("sig[0-9]", mm$Parametros), "Variância",
                         ifelse(grepl("r_Y[0-9]", mm$Parametros), "Correlação",
                                ifelse(grepl("nu", mm$Parametros), "Dispersão",
                                       ifelse(grepl("phi", mm$Parametros), "Dispersão", "Deu Ruim")))))
mm <- mm[order(mm$MainModel),]
mm[mm$Tipo=="Dispersão","Parametros"] <- c(rep(paste0("exp(nu)", 1:3),4), rep(paste0("exp(phi)", 1:3),4))
mm$Response <- factor(substr(mm$Parametros, nchar(mm$Parametros), nchar(mm$Parametros)),
                      labels = c("Nmsp", "Nmosp", "Nspfy"))
names(mm)[1:4] <- c("Estimativas", "ErroPadrao", "Z", "p value")
mm$NAN <- factor(as.character(ifelse(is.nan(mm$ErroPadrao), 1, 20)))
mm$ModelSpec <- ifelse(mm$Model%in%c("COM-Poisson", "NB", "Poisson"), "Full",
                       ifelse(mm$Model%in%c("Fixed Disp COM-Poisson", "Fixed Disp NB"), "Fixed Disp",
                              ifelse(mm$Model%in%c("Common Var COM-Poisson", "Common Var NB"), "Common Var",
                                     ifelse(mm$Model%in%c("Fixed Var COM-Poisson", "Fixed Var NB"), "Fixed Var",
                                            ifelse(mm$Model%in%c("Rho Zero Poisson", "Rho Zero NB", "Rho Zero COM-Poisson"), "Rho 0", "Deu Ruim")))))
#####################################################
####### Beta Graphic
#####################################################
mm %>% 
  filter(Tipo=="Média") %>% 
  mutate(Parametros = factor(gsub("_Y\\d", "", Parametros), 
                             levels = c("b0", "b1", "b2", "b3", "b4"),
                             labels = c("Intercepto", "Race", "Education", "Marital", "Age"))) %>% 
  ggplot(aes(x = Parametros, y = Estimativas)) +
  geom_pointrange(aes(ymin = Estimativas-1.96*ErroPadrao,
                      ymax = Estimativas+1.96*ErroPadrao,
                      linetype = Model,
                      shape = NAN,
                      color = ModelSpec),
                  position = position_dodge(width = 0.65),
                  alpha = .8) +
  scale_shape_manual(values = c(1,20),
                     labels = c("Without SE", "With SE")) +
  scale_linetype_manual(values=1:12) +
  guides(linetype = guide_legend(override.aes = list(shape = NA),
                                 ncol = 3),
         shape = guide_legend(override.aes = list(linetype = "blank"),
                              ncol = 1),
         color = guide_legend(ncol = 1)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        text = element_text(size = 17),
        legend.position = "bottom",
        legend.title = element_blank()) +
  facet_wrap(~Response, scales = "free_x") +
  coord_flip() +
  labs(y = TeX('$\\hat{\\beta}\\pm 1.96 SE}$'),
       linetype = "Model")
pathfig <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/Texto/LaTeX/Figuras"
# ggsave(paste0(pathfig, "/nhanes_beta_final_model.pdf"), units = "cm")
ggsave(paste0(pathfig, "/nhanes_beta_initial_model.pdf"), width = 33, height = 22, units = "cm")
#####################################################
####### Dispersion
#####################################################
myHeader <- c(1,2,2,2)

# names(myHeader) <- c(" ", "Negative binomial ($\\\\phi$)", "COM-Poisson ($\\\\nu$)")
names(myHeader) <- c(" ", "Nmsp", "Nmosp", "Nspfy")
mm %>% 
  filter(Tipo%in%c("Dispersão")) %>% 
  select(-Z, -`p value`, -MainModel, -Tipo, -Parametros, -NAN, -ModelSpec) %>% 
  pivot_wider(names_from = Response,
              values_from = c("Estimativas", "ErroPadrao")) %>% 
  select(Model, Estimativas_Nmsp, ErroPadrao_Nmsp,
         Estimativas_Nmosp, ErroPadrao_Nmosp, everything()) %>% 
  kable(caption = toupper("Dispersion parameter estimates and standard errors (SE) for each initial model and outcome of NHANES data"),
        label = "nhanesdispinit",
        booktabs = T,
        digits = 2,
        row.names = F,
        escape = F,
        format = "latex",
        col.names = c("Model", rep(c("Est", "SE"),3)),
        align = "c") %>% 
  add_header_above(header = myHeader,
                   escape = F) %>% 
  add_footnote(label = 'Empty cells comprises of non returned SE')


#####################################################
####### Variância
#####################################################
mm %>% 
  filter(Tipo=="Variância") %>% 
  filter(!ModelSpec=='Common Var') %>% 
  ggplot(aes(x = Model, y = Estimativas)) +
  geom_pointrange(aes(ymin = Estimativas-1.96*ErroPadrao,
                      ymax = Estimativas+1.96*ErroPadrao,
                      # linetype = Model,
                      shape = NAN,
                      color = ModelSpec),
                  position = position_dodge(width = 0.3),
                  alpha = .8) +
  scale_shape_manual(values = c(1,20),
                     labels = c("Without SE", "With SE")) +
  guides(shape = guide_legend(override.aes = list(linetype = "blank"),
                              nrow = 1),
         color = guide_legend(nrow = 1)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        text = element_text(size = 16),
        legend.position = "bottom",
        legend.title = element_blank()) +
  facet_wrap(~Response, scales = "free_x") +
  coord_flip() +
  labs(y = TeX('$\\hat{\\sigma}\\pm 1.96 SE}$'),
       linetype = "Model")
pathfig <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/Texto/LaTeX/Figuras"
# ggsave(paste0(pathfig, "/nhanes_beta_final_model.pdf"), units = "cm")
ggsave(paste0(pathfig, "/nhanes_var_initial_model.pdf"), width = 28, height = 19, units = "cm")

mm %>% 
  filter(Tipo=="Variância") %>% 
  filter(ModelSpec=='Common Var') %>% 
  select(Model, Estimativas, ErroPadrao) %>% 
  mutate(Ic.Inf = Estimativas-1.96*ErroPadrao,
         Ic.Sup = Estimativas+1.96*ErroPadrao) %>% 
  kable(digits = 4, align = c('l', rep('c', 4)),
        format = 'latex',
        label = 'nhanesvarcommon',
        escape = F,
        booktabs = T,
        col.names = c('Model', 'Estimates', 'SE', 'Lower CI', 'Upper CI'),
        caption = "STANDARD DEVIATION ESTIMATES OF RANDOM EFFECTS AND 95% CONFIDENCE INTERVALS FOR MODELS WITH COMMON VARIANCE") %>% 
  add_footnote(label = 'Empty cells comprises of non returned SE')

#####################################################
####### Correlation
#####################################################

myHeader <- c(1,2,2,2)
names(myHeader) <- c(" ", "Nmsp X Nmosp", "Nmsp X Nspfy", "Nmosp X Nspfy")
mm %>% 
  filter(Tipo%in%c("Correlação")) %>% 
  select(-Z, -`p value`, -MainModel, -Tipo, -Response, -NAN, -ModelSpec) %>% 
  pivot_wider(names_from = Parametros,
              values_from = c("Estimativas", "ErroPadrao")) %>% 
  select(Model, 
         Estimativas_r_Y1_Y2, ErroPadrao_r_Y1_Y2,
         Estimativas_r_Y1_Y3, ErroPadrao_r_Y1_Y3,
         everything()) %>% 
  kable(caption = toupper("Correlation parameter estimates and standard errors (SE) for each initial model and outcome of NHANES data"),
        label = "nhanescorinit",
        booktabs = T,
        digits = 3,
        row.names = F,
        escape = F,
        format = "latex",
        col.names = c("Model", rep(c("Est", "SE"),3)),
        align = "c") %>% 
  add_header_above(header = myHeader,
                   escape = F) %>% 
  add_footnote(label = 'Empty cells comprises of non returned SE')


