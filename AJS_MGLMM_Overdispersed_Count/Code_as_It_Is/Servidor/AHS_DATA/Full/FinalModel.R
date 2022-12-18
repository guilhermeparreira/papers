######################
# Compoisson Multivariate ------------------------------------------------------------------------
######################
rm(list=ls())
mainpath <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/TMB/MultivariateModels/Servidor/"
setwd(paste0(mainpath, "AHS_DATA/Full"))
source(paste0(mainpath, "unify_packages.R"))
library(TMB)
load("compoisson_1.RData")
load("negbin.RData")
load("Poisson_final.RData")
library(ggplot2)
library(dplyr)
library(latex2exp)
# The order of the response variables is different in the results session, compared to the application. 
# If I will change it in the future, I will keep the order from the results
#################################################################
############# Graphic for a ll 
#################################################################

np_poisson <- nrow(estimates_poisson_sd$ests)
np_cmp <- nrow(estimates_compoisson_sd$ests)
dg <- data.frame(Parametros = c(row.names(estimates_poisson_sd$ests), row.names(estimates_negbin_sd$ests), row.names(estimates_compoisson_sd$ests)),
                 Estimativas = c(estimates_poisson_sd$ests[,1], estimates_negbin_sd$ests[,1], estimates_compoisson_sd$ests[,1]),
                 ErroPadrao = c(estimates_poisson_sd$ests[,2], estimates_negbin_sd$ests[,2], estimates_compoisson_sd$ests[,2]),
                 Distribuição = c(rep("Poisson", np_poisson), rep("NB", np_cmp), rep("COM-Poisson", np_cmp)),
                 Tipo = c(c(rep("Média", 55), rep("Variância", 5), rep("Correlação", 10)), #Poisson
                          c(rep("Média", 55), rep("Variância", 5), rep("Correlação", 10), rep("Dispersão", 5)), #NB
                          c(rep("Média", 55), rep("Variância", 5), rep("Correlação", 10), rep("Dispersão", 5)))) #COM

# (STARTING HERE!!!!)
# head(dg,20)
# dg[c(33:35, 52:54), 1] <- c(paste0("exp(phi)", 1:3), paste0("exp(nu)", 1:3))
dg[dg$Tipo=="Dispersão", 1] <- c(paste0("exp(phi)", 1:5), paste0("exp(nu)", 1:5))
dg$Response <- ifelse(dg$Tipo=="Média", gsub(".*_Y", "", dg$Parametros),
                             ifelse(dg$Tipo=="Variância", gsub("sig", "", dg$Parametros),
                                    ifelse(dg$Tipo=="Correlação", gsub("r_", "", dg$Parametros),
                                           ifelse(dg$Tipo=="Dispersão", gsub("[^0-9.-]", "", dg$Parametros), "BU"))))
varsnames <- c("Ndoc", "Nndoc", "Nmed", "Nhosp", "Nadm")
# names(ahs)
# colnames(data$Y)
respnames <- data.frame(Number = as.character(1:5),
                        Names = factor(varsnames, 
                                       levels = varsnames))
# Real names
head(dg)
# dg <- dg[, -7]
dg <- left_join(dg, respnames, by = c("Response" = "Number"))
# save.image("compoisson_nu_fixed.RData")
# load("compoisson_nu_fixed.RData")
## Only beta - Professor does not suggest
data$X
dg$Response <- as.numeric(dg$Response)
dg$NAN <- factor(as.character(ifelse(is.nan(dg$ErroPadrao), 1, 16)))
# Esqueci de incluir o hscore no estudo, mas no artigo do prof, ainda tem 1 parâmetro faltando
# eu vi que no seu artigo da Royal society, vc estima um total de 10 parêmetros de média.
# eu estimei 11 hehe
# o certo era ter estimado 12
# podemos ver isso dps
# 
# colnames(data$X)
dgmedia <- dg %>% 
  filter(Tipo=="Média") %>% 
  mutate(Parametros = factor(gsub("_Y\\d+", "", Parametros), 
                             labels = c("intercept", "sex(fem)", "age", "income", "levyplus", "freepoor", "freerepa", 
                                        "illnes", "actdays", "chcond not limited", "chcond limited"),
                             levels = paste0("b",0:10)))

dgmedia$Parametros2 <- factor(ifelse(dgmedia$Parametros=="intercept", "Intercept", "Covariates"))
########################################
# Graphic for average only!!!
########################################

dg1 <- dgmedia %>% 
  filter(Parametros2 == "Covariates") %>% 
  ggplot() +
  geom_pointrange(aes(x = Parametros, 
                      y = Estimativas,
                      ymin = Estimativas-1.96*ErroPadrao,
                      ymax = Estimativas+1.96*ErroPadrao,
                      linetype = Distribuição),
                  # It controls the vertical distance between models
                  position = position_dodge(width = .6),
                  # alpha = .5,
                  fatten = 2,
                  size = 1,
                  pch = 16
  ) +
  theme_bw() +
  theme(axis.title = element_blank(),
        text = element_text(size = 17),
        legend.position = "bottom") +
  scale_linetype_manual(values=c("solid", "dashed", "dotted")) +
  # guides(linetype = guide_legend(override.aes = list(shape = NA))) +
  guides(linetype = guide_legend(keyheight = 2.9)) +
  # guides(linetype = guide_legend(override.aes = list(shape = NA),keywidth = 3.5)) +
  # guides(linetype = guide_legend(keyheight = 3.5)) +
  facet_grid(Parametros2~Names, scales = "free") +
  coord_flip() +
  labs(y = TeX('$\\hat{\\beta}\\pm {1.96 SE}$'),
       linetype = "Model")
# dg1
dg2 <- dgmedia %>% 
  filter(Parametros2 == "Intercept") %>% 
  ggplot() +
  geom_pointrange(aes(x = Parametros, 
                      y = Estimativas,
                      ymin = Estimativas-1.96*ErroPadrao,
                      ymax = Estimativas+1.96*ErroPadrao,
                      linetype = Distribuição),
                  # It controls the vertical distance between models
                  position = position_dodge(width = .6),
                  # alpha = .5,
                  fatten = 2,
                  size = 1,
                  pch = 16
  ) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        text = element_text(size = 17),
        legend.position = "bottom") +
  scale_linetype_manual(values=c("solid", "dashed", "dotted")) +
  guides(linetype = guide_legend(override.aes = list(shape = NA))) +
  guides(linetype = guide_legend(keyheight = 2.9)) +
  facet_grid(Parametros2~Names, scales = "free") +
  coord_flip() +
  labs(y = TeX('$\\hat{\\beta}\\pm {1.96 SE}$'),
       linetype = "Model")
# dg2
library(ggpubr)
ggarrange(dg1, dg2, nrow = 2,
          common.legend = T, legend = "bottom",
          align = "v",
          heights = c(4,1.2))
pathfig <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/Texto/LaTeX/Figuras"
ggsave(paste0(pathfig, "/ahs_beta_final_model.pdf"), width = 35, height = 25, units = "cm")

### FOR SMJ JOURNAL

dg3 <- dgmedia %>% 
  filter(Parametros2 == "Covariates") %>% 
  ggplot() +
  geom_pointrange(aes(x = Parametros, 
                      y = Estimativas,
                      ymin = Estimativas-1.96*ErroPadrao,
                      ymax = Estimativas+1.96*ErroPadrao,
                      shape = Distribuição),
                  # It controls the vertical distance between models
                  position = position_dodge(width = .6),
                  # alpha = .5,
                  fatten = 2,
                  size = .75
  ) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        text = element_text(size = 17),
        legend.position = "bottom",
        strip.text.y = element_blank()) +
  guides(shape = guide_legend(keyheight = 2.9)) +
  facet_grid(Parametros2~Names, scales = "free") +
  coord_flip() +
  labs(y = TeX('$\\hat{\\beta}\\pm {1.96 SE}$'),
       shape = "Model")
dg3
pathfig2 <- '/home/guilherme/GoogleDrive/Mestrado/dissertacao/Texto/LaTeX/ArtigoRevista/SuperDispersao_SMj'
ggsave(paste0(pathfig2, "/ahs_beta_final_model_covariates.pdf"), width = 35, height = 20, units = "cm")

#####################
# Only standard deviation
#####################

dg %>% 
  filter(Tipo %in% c("Variância")) %>% 
  ggplot() +
  geom_pointrange(aes(x = Distribuição, 
                      y = Estimativas,
                      ymin = Estimativas-1.96*ErroPadrao,
                      ymax = Estimativas+1.96*ErroPadrao,
                      linetype = Distribuição),
                  # It controls the vertical distance between models
                  position = position_dodge(width = .6),
                  size = .6) +
  coord_flip() +
  facet_wrap(~Names, scales = "free", nrow = 1) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 15),
        legend.position = "bottom") +
  scale_linetype_manual(values=c("solid", "dashed", "dotted")) +
    guides(linetype = guide_legend(keyheight = 2.5)) +
  labs(y = TeX('$\\hat{\\sigma}\\pm 1.96 SE}$'),
       linetype = "Model")

pathfig <- "/home/guilherme/GoogleDrive/Mestrado/dissertacao/Texto/LaTeX/Figuras"
ggsave(paste0(pathfig, "/ahs_beta_final_model_var.pdf"), width = 28, height = 8.5, units = "cm")


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
      caption = toupper("Dispersion of parameter estimates and standard errors (SE) for each model and outcome of AHS data"),
      label = "ahsdisp",
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
      caption = toupper("Standard deviation of random effect estimates and standard errors (SE) for each model and outcome of AHS data"),
      label = "ahsvar",
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
# Correlation v2
########################################
## Poisson
aa <- round(estimates_poisson_sd$ests,2)[61:70,]
aa <- round(estimates_negbin_sd$ests,2)[61:70,]
aa <- round(estimates_compoisson_sd$ests,2)[61:70,]
matcrazy <- matrix(numeric(0), ncol = 5, nrow = 5)
diag(matcrazy) <- 1

# c <- 1
# i <- 1
for (i in 1:10){
  nm <- row.names(aa)[i]
  row <- as.numeric(substr(nm, 4,4))
  col <- as.numeric(substr(nm, 7,7))
  fim <- ifelse(aa[i, 4]<.05, ")^{*}", ")")
  matcrazy[row, col] <- paste0(aa[i, 1],"(", aa[i,2],fim)
}
dfcrazy <- as.data.frame(matcrazy)
row.names(dfcrazy) <- c("Ndoc", "Nndoc", "Nmed", "Nhosp", "Nadm")
colnames(dfcrazy) <- c("Ndoc", "Nndoc", "Nmed", "Nhosp", "Nadm")
# dfcrazy <- dfcrazy[-5, -1]
kable(dfcrazy, 
      row.names = F,
      escape = F, 
      format = "latex",
      booktabs = T)