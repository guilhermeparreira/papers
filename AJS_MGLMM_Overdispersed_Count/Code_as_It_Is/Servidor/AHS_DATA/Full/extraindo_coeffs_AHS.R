rm()
path <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor"
localp <- "/AHS_DATA/Full"
locals <- "/AHS_DATA/SmallerModels"
setwd(paste0(path, localp))
library(reshape2)
library(latex2exp)
library(ggplot2)
library(tidyr)
library(dplyr)
library(kableExtra)
library(knitr)
options(knitr.kable.NA = '')
#################################################################################
######################### Só preciso fazer exp(disp) - tanto no NB como no COMPOISSON
#################################################################################
##########################
# Model 1 - Poisson
##########################
load("Poisson_final.RData")
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
m03 <- as.data.frame(estimates_negbin_sd$ests)
m03$Model <- 'Fixed Disp NB'
m03$MainModel <- "NB"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "path", "localp", "locals")])
##########################
# Model 4 - Negbin - Comum Variance
##########################
load("Negbin_var_comum.RData")
m04 <- as.data.frame(estimates_negbin_sd$ests)
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
load("compoisson_1.RData")
m06 <- as.data.frame(estimates_compoisson_sd$ests)
m06$Model <- "COM-Poisson"
m06$MainModel <- "COM-Poisson"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06" ,"path", "localp", "locals")])
##########################
# Model 7 - Compoisson Disp Fixed 
##########################
setwd(paste0(path, locals))
load("compoisson_nu_fixed.RData")
# estimates_compoisson_s
m07 <- as.data.frame(estimates_compoisson_sd$ests)
m07$Model <- 'Fixed Disp COM-Poisson'
m07$MainModel <- "COM-Poisson"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06", "m07" ,"path", "localp", "locals")])
##########################
# Model 8 - Compoisson Comum Var
##########################
load("compoisson_var_comum_1.RData")
b1p <- melt(estimates_compoisson$beta)
b1p <- b1p[,-c(1:2), drop = FALSE]
s1p <- melt(estimates_compoisson$sigma)
r1p <- melt(estimates_compoisson$rho)
d1p <- exp(melt(estimates_compoisson$disp)) # Transformação aqui.
m08 <- data.frame(Estimate = unname(rbind(b1p, s1p, r1p, d1p)),
                  `Std. Error` = NA,
                  Z = NA,
                  `p value` = NA,
                  check.names = F,
                  row.names = c(row.names(m01)[1:55], 
                                'sig1', 
                                row.names(m01)[61:70], 
                                paste0('exp(nu)',1:5)))
m08$Model <- "Common Var COM-Poisson"
m08$MainModel <- "COM-Poisson"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06", "m07", "m08","path", "localp", "locals")])
##########################
# Model 9 - Compoisson Fixed Var
##########################
load("compoisson_var_fixed_1.RData")
b1p <- melt(estimates_compoisson$beta)
b1p <- b1p[,-c(1:2), drop = FALSE]
r1p <- melt(estimates_compoisson$rho)
d1p <- exp(melt(estimates_compoisson$disp)) # Transformação aqui.
m09 <- data.frame(Estimate = unname(rbind(b1p, r1p, d1p)),
                  `Std. Error` = NA,
                  Z = NA,
                  `p value` = NA,
                  check.names = F,
                  row.names = c(row.names(m01)[1:55], 
                                row.names(m01)[61:70], 
                                paste0('exp(nu)',1:5)))
m09$Model <- "Fixed Var COM-Poisson"
m09$MainModel <- "COM-Poisson"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06", "m07", "m08", "m09" ,"path", "localp", "locals")])
##########################
# Model 10 - Poisson Rho fixed
##########################
setwd(paste0(path, locals))
load("Poisson_rho_fixed.RData")
m10 <- as.data.frame(estimates_poisson_sd$ests)
m10$Model <- "Rho Zero Poisson"
m10$MainModel <- "Poisson"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06", "m07", "m08", "m09", "m10" ,"path", "localp", "locals")])
##########################
# Model 11 - NB Rho fixed
##########################
load("negbin_rho_fixed.RData")
# estimates_compoisson_s
m11 <- as.data.frame(estimates_negbin_sd1$ests)
m11$Model <- "Rho Zero NB"
m11$MainModel <- "NB"
rm(list = ls()[!ls()%in%c("m01", "m02", "m03", "m04", "m05", "m06", "m07", "m08", "m09", "m10", "m11" ,"path", "localp", "locals")])
##########################
# Model 12 - COM-Poisson rho fixed
##########################
load("compoisson_rho_fixed.RData")
b1p <- melt(estimates_compoisson_s$beta)
b1p <- b1p[,-c(1:2), drop = FALSE]
s1p <- melt(estimates_compoisson_s$sigma)
d1p <- exp(melt(estimates_compoisson_s$disp))
m12 <- data.frame(Estimate = unname(rbind(b1p, s1p, d1p)),
                  `Std. Error` = NA,
                  Z = NA,
                  `p value` = NA,
                  check.names = F,
                  row.names = c(row.names(m10), paste0('exp(nu)',1:nr)))
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
mm[mm$Tipo=="Dispersão","Parametros"] <- c(rep(paste0("exp(nu)", 1:5),4), rep(paste0("exp(phi)", 1:5),4))
varsnames <- c("Ndoc", "Nndoc", "Nmed", "Nhosp", "Nadm")
mm$Response <- factor(substr(mm$Parametros, nchar(mm$Parametros), nchar(mm$Parametros)),
                      labels = varsnames)
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
mmedia <- mm %>% 
  filter(Tipo=="Média") %>% 
  mutate(Parametros = factor(gsub("_Y\\d+", "", Parametros), 
                             labels = c("intercept", "sex(fem)", "age", "income", "levyplus", "freepoor", "freerepa", 
                                        "illnes", "actdays", "chcond not limited", "chcond limited"),
                             levels = paste0("b",0:10)),
         Parametros2 = factor(ifelse(Parametros=="intercept", "Intercept", "Covariates")))
mmedia <- droplevels(mmedia)
dg1 <- mmedia %>% 
  filter(Parametros2 == "Covariates") %>% 
  ggplot(aes(x = Parametros, y = Estimativas)) +
  geom_pointrange(aes(ymin = Estimativas-1.96*ErroPadrao,
                      ymax = Estimativas+1.96*ErroPadrao,
                      linetype = Model,
                      shape = NAN,
                      color = ModelSpec),
                  position = position_dodge(width = 0.7),
                  alpha = .8) +
  theme_bw() +
  theme(axis.title = element_blank(),
        text = element_text(size = 15),
        legend.position = "bottom",
        legend.title = element_blank()) +
  scale_shape_manual(values = c(20, 1),
                     labels = c("With SE", "Without SE")) +
  scale_linetype_manual(values=1:12) +
  guides(linetype = guide_legend(override.aes = list(shape = NA),
                                 ncol = 3),
         shape = guide_legend(override.aes = list(linetype = "blank"),
                              ncol = 1),
         color = guide_legend(ncol = 1)) +
  facet_grid(Parametros2~Response, scales = "free_x") +
  coord_flip() +
  labs(y = TeX('$\\hat{\\beta}\\pm 1.96 SE}$'),
       linetype = "Model")
dg2 <- mmedia %>% 
  filter(Parametros2 == "Intercept") %>% 
  ggplot(aes(x = Parametros, y = Estimativas)) +
  geom_pointrange(aes(ymin = Estimativas-1.96*ErroPadrao,
                      ymax = Estimativas+1.96*ErroPadrao,
                      linetype = Model,
                      shape = NAN,
                      color = ModelSpec),
                  position = position_dodge(width = 0.7),
                  fatten = 2,
                  alpha = .8) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        text = element_text(size = 15),
        legend.position = "bottom",
        legend.title = element_blank(),
        strip.text.x = element_blank()) +
  scale_shape_manual(values = c(20, 1),
                     labels = c("With SE", "Without SE")) +
  scale_linetype_manual(values=1:12) +
  guides(linetype = guide_legend(override.aes = list(shape = NA),
                                 ncol = 3),
         shape = guide_legend(override.aes = list(linetype = "blank"),
                              ncol = 1),
         color = guide_legend(ncol = 1)) +
  facet_grid(Parametros2~Response, scales = "free_x") +
  coord_flip() +
  labs(y = TeX('$\\hat{\\beta}\\pm 1.96 SE}$'),
       linetype = "Model")
# dg2
library(ggpubr)
ggarrange(dg1, dg2, nrow = 2,
          common.legend = T, legend = "bottom",
          align = "v",
          heights = c(4,1.3))
pathfig <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/Texto/LaTeX/Figuras"
# ggsave(paste0(pathfig, "/nhanes_beta_final_model.pdf"), units = "cm")
ggsave(paste0(pathfig, "/ahs_beta_initial_model.pdf"), width = 28, height = 26, units = "cm")
#####################################################
####### Dispersion
#####################################################
myHeader <- c(1,2,2,2,2,2)

# names(myHeader) <- c(" ", "Negative binomial ($\\\\phi$)", "COM-Poisson ($\\\\nu$)")
names(myHeader) <- c(" ", "Ndoc", "Nndoc", "Nmed", "Nhosp", "Nadm")
mm %>% 
  filter(Tipo%in%c("Dispersão")) %>% 
  select(-Z, -`p value`, -MainModel, -Tipo, -Parametros, -NAN, -ModelSpec) %>% 
  pivot_wider(names_from = Response,
              values_from = c("Estimativas", "ErroPadrao")) %>% 
  select(Model, 
         Estimativas_Ndoc, ErroPadrao_Ndoc,
         Estimativas_Nndoc, ErroPadrao_Nndoc,
         Estimativas_Nmed, ErroPadrao_Nmed,
         Estimativas_Nhosp, ErroPadrao_Nhosp,
         Estimativas_Nadm, ErroPadrao_Nadm) %>% 
  kable(caption = toupper("Dispersion parameter estimates and standard errors (SE) for each initial model and outcome of ANT data"),
        label = "ahsdispinit",
        booktabs = T,
        digits = 2,
        row.names = F,
        escape = F,
        format = "latex",
        col.names = c("Model", rep(c("Est", "SE"),5)),
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
  facet_wrap(~Response, scales = "free_x", nrow = 1) +
  coord_flip() +
  labs(y = TeX('$\\hat{\\sigma}\\pm 1.96 SE}$'),
       linetype = "Model")
pathfig <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/Texto/LaTeX/Figuras"
# ggsave(paste0(pathfig, "/nhanes_beta_final_model.pdf"), units = "cm")
ggsave(paste0(pathfig, "/ahs_var_initial_model.pdf"), width = 28, height = 19, units = "cm")

mm %>% 
  filter(Tipo=="Variância") %>% 
  filter(ModelSpec=='Common Var') %>% 
  select(Model, Estimativas, ErroPadrao) %>% 
  mutate(Ic.Inf = Estimativas-1.96*ErroPadrao,
         Ic.Sup = Estimativas+1.96*ErroPadrao) %>% 
  kable(digits = 4, align = c('l', rep('c', 4)),
        format = 'latex',
        label = 'ahsvarcommon',
        escape = F,
        booktabs = T,
        col.names = c('Model', 'Estimates', 'SE', 'Lower CI', 'Upper CI'),
        caption = "STANDARD DEVIATION ESTIMATES OF RANDOM EFFECTS AND 95% CONFIDENCE INTERVALS FOR MODELS WITH COMMON VARIANCE") %>% 
  add_footnote(label = 'Empty cells comprises of non returned SE')

#####################################################
####### Correlation
#####################################################
dgcorr <- mm %>% 
  filter(Tipo == "Correlação")
dgcorr$Resp <- gsub("r_", "", dgcorr$Parametros)
dgcorr$Resp <- gsub("Y", "", dgcorr$Resp)
dgcorr <- dgcorr[, -c(6,7,8)]
resps2 <- stringr::str_split_fixed(dgcorr$Resp, "_", 2)
colnames(resps2) <- c("Var1", "Var2")
dgcorr <- cbind(dgcorr, resps2)
dgcorr$Var1 <- as.numeric(dgcorr$Var1)
dgcorr$Var2 <- as.numeric(dgcorr$Var2)
dgcorr$IcInf <- dgcorr$Estimativas-1.96*dgcorr$ErroPadrao
dgcorr$IcSup <- dgcorr$Estimativas+1.96*dgcorr$ErroPadrao
dgcorr$Valorp <- pnorm(abs(dgcorr$Estimativas/dgcorr$ErroPadrao), lower.tail = F)*2


# Agora operação para colocar em formato de matrix
df_to_mat <- function(dd, nr = 5){
  mat_est <- matrix(numeric(), nrow = nr, ncol = nr)
  mat_se <- matrix(numeric(), nrow = nr, ncol = nr)
  mat_inf <- matrix(numeric(), nrow = nr, ncol = nr)
  mat_sup <- matrix(numeric(), nrow = nr, ncol = nr)
  mat_p <- matrix(numeric(), nrow = nr, ncol = nr)
  for (i in 1:nr){
    for (j in i:nr){
      if (i == j){
        mat_est[i,j] <- 1
        mat_se[i,j] <- 0
      } else{
        dgcorrfilt <- dd[dd$Var1 %in% i & dd$Var2 %in% j, , drop = F]
        mat_est[i, j] <- dgcorrfilt[, 1]
        mat_se[i, j] <- dgcorrfilt[, 2]
        mat_inf[i, j] <- as.numeric(ifelse(dgcorrfilt[, 12]>1, 1, ifelse(dgcorrfilt[, 12]<(-1),-1, dgcorrfilt[, 12])))
        mat_sup[i, j] <- as.numeric(ifelse(dgcorrfilt[, 13]>1, 1, ifelse(dgcorrfilt[, 13]<(-1),-1, dgcorrfilt[, 13])))
        mat_p[i, j] <- dgcorrfilt[, 4]
      }
    }
  }
  ll <- list(mat_est, mat_se, mat_inf, mat_sup, mat_p)
  names(ll) <- c("est", "se", "inf", "sup", "p")
  return(ll)
}
library(corrplot)
gera_grafico <- function(data, namefile, pathfig){
  dgcorrP <- data
  corrp <- df_to_mat(dgcorrP)
  pdf(file = paste0(pathfig, namefile), width = 4, height = 2.5)
  graph1 <- corrplot(
    corrp$est
    # mar = c(0,0,0,0), 
    ,tl.cex = 0.6 
    ,cl.cex = .7
    # pch = "*",
    ,type = "upper"
    ,diag = F
    ,outline = F
    # ,pch.cex = 1.7       # Controls the size of the star
    # p.mat = 1-corrp$p, #Opposite
    # ,sig.level = .95
  )   #Opposite
  # print(graph1)
  dev.off()
}

gera_grafico2 <- function(data, namefile, pathfig){
  dgcorrP <- data
  corrp <- df_to_mat(dgcorrP)
  pdf(file = paste0(pathfig, namefile), width = 4, height = 2.5)
  graph1 <- corrplot(
    corrp$est
    # mar = c(0,0,0,0), 
    ,tl.cex = 0.6 
    ,cl.cex = .7
    ,pch = "*"
    ,type = "upper"
    ,diag = F
    ,outline = F
    ,pch.cex = 1.7       # Controls the size of the star
    ,p.mat = 1-corrp$p #Opposite
    ,sig.level = .95
  )   #Opposite
  # print(graph1)
  dev.off()
}

# Decidi excluir as estrelas de significância, pois apenas 2 modelos tinham significância estatística.
dgcorr$Valorp2Check <- ifelse(dgcorr$`p value`<.05,1,0)

modelszero = c("Common Var COM-Poisson", "Fixed Var COM-Poisson")
for (i in 1:length(modelszero)) {
  nm <- modelszero[i]
  gera_grafico(dgcorr[dgcorr$Model==nm, ], paste0("/cor_ahs_",nm,"_initial.pdf"), pathfig)
}

# Refaz só para aquelas que tem coef significativo (warning é normal)
modelssig = with(dgcorr[dgcorr$Valorp2Check==1,], unique(Model))
modelssig = modelssig[!is.na(modelssig)]
for (i in 1:length(modelssig)) {
  nm <- modelssig[i]
  gera_grafico2(dgcorr[dgcorr$Model==nm, ], paste0("/cor_ahs_",nm,"_initial.pdf"), pathfig)
}
mm %>% 
  filter(Tipo == "Correlação") %>% 
  filter(Model == "Common Var NB")
