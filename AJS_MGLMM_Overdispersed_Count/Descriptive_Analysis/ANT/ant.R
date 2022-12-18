rm(list=ls())
install.packages('mvabund', dependencies = T)
install.packages('facetscales', dependencies = T)
install.packages('GGally', dependencies = T)
library(mvabund);library(ggplot2)
library(dplyr);library(reshape2)
library(facetscales);library(knitr)
library(kableExtra);library(GGally)
setwd("/home/guilherme/GoogleDrive/Mestrado/dissertacao/Texto/LaTeX/Figuras/")

# Vars Description #############################################################

data.des <- data.frame(Name = c("Bare ground", "Canopy cover", "Shrub.cover", "Volume.lying.CWD", "Feral.mammal.dung"),
                       Description = c("Percent cover of bare ground, as estimated from ten 1x1 metre quadrats",
                                       "Percent canopy cover, as estimated from two 20x20m transects",
                                       "Percent canopy cover, as estimated from two 20x20m transects",
                                       "Estimated volume of Coarse Woody Debris in two 20x20m transects, including all debris >5cm diameter",
                                       "Proportion of quadrats including mammal dung, out of ten 1x1m quadrats"))
kable(data.des, format = "latex", booktabs = T, digits = 3, escape = F,
      row.names = F,
      align = c("l"),
      caption = "COVARIATES COLLECTED IN THE ANT STUDY IN SOUTH-EASTERN AUSTRALIA") %>% 
  kable_styling() %>% 
  column_spec(2, width = "11.5cm")

# Code #########################################################################
data("antTraits")
y <- antTraits$abund # Selecting response variables
X <- antTraits$env # Selecting covariates
head(y[,1:5])
head(X[,1:5])
# Pilosity will be considered as numeric
# Polymorphism I think that needs to be considered as factor
data(antTraits)
# ft = traitglm(antTraits$abund,antTraits$env,antTraits$traits) #to do a fourth corner analysis

# Bare.ground
# Percent cover of bare ground, as estimated from ten 1x1 metre quadrats
# Percentual de solo descoberto. Apenas terra aberta

# Canopy.cover
# Percent canopy cover, as estimated from two 20x20m transects
# Percentual de cobertura do dossel (porcentagem da ocupação das copas das árvores dentro do transecto)
# Deve servir como um grau de ocupação das copas
# Dossel é a parte superior das árvores
# Percentual de cobertura do dossel (copa da árvore)

# Shrub.cover
# Percent canopy cover, as estimated from two 20x20m transects
# Percentual de cobertura de arbusto?

# A terceira é um pouco estranha, mas parece ser um volume de madeira com tamanho maior que 5 cm de diâmetro. 
# Eu imagino que eles mediram os detritos de madeira na florestal e se havia um galho com diâmetro maior que 5 cm, ele tinha o volume contabilizado

# Volume.lying.CWD
# Estimated volume of Coarse Woody Debris in two 20x20m transects, including all debris >5cm diameter.
# "Volume deitado": Volume estimado de detritos lenhosos na área.

# Feral.mammal.dung
# Proportion of quadrats including mammal dung, out of ten 1x1m quadrats.
# Percentual de insetos?

# Names ###########################################################################
dy = data.frame(ny = names(y))
dy <- dy %>% 
  separate(ny, into = c('first', 'sec'), extra = 'merge') %>% 
  mutate(y = paste0(abbreviate(first), '.', abbreviate(sec)))
dy
# Barplot ###########################################################################
nyabrev <- c("Amblyopone", "Aphaenogaster", "Camponotus.Ci", "Camponotus.Cl", 
  "Camponotus.Co", "Camponotus.Ni", "Camponotus.Ns", "Cardiocondyla", "Crematogaster", 
  "Heteroponera", "Iridomyrmex.Bi", "Iridomyrmex.Dr", "Iridomyrmex.Mj", 
  "Iridomyrmex.Pu", "Iridomyrmex.Ru", "Iridomyrmex.Su", "Iridomyrmex.Ss", "Melophorus.E", 
  "Melophorus.F", "Melophorus.H", "Meranoplus.A", "Monomorium.L", "Monomorium.Ro", 
  "Monomorium.Sy", "Myrmecia", "Notoncus.Ca", "Notoncus.Ec", "Nylanderia", 
  "Ochetellus", "Paraparatrechina", "Pheidole.A", "Pheidole.B", "Pheidole.E", 
  "Pheidole.J", "Polyrhachis", "Rhytidoponera.A", "Rhytidoponera.B", 
  "Solenopsis", "Stigmacros", "Tapinoma", "Tetramorium")

names(y) <- paste0(format(1:ncol(y), width = 2), ".", nyabrev)
long <- melt(y)
ggplot(data = long,
       aes(x = value)) +
  geom_bar(stat = "count") +
  theme_bw() +
  labs(x = "Occurrence",
       y = "Frequency") +
  facet_wrap(~variable, scales = "free", nrow = 8) +
  theme(strip.text.x = element_text(size = 6),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8)) +
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
  scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))
ggsave("desc_ant.pdf", width = 15.5, height = 19.5, units = "cm")

# Slide

long <- melt(y)
ord <- names(sort(sapply(sapply(y, table), max), decreasing = T))
y2 <- y[, ord]
long <- melt(y2)
pdf("desc_ant_2.pdf", width=13)
ggplot(data = long,
       aes(x = value)) +
  geom_bar(stat = "count") +
  theme_bw() +
  labs(x = "Occurrence",
       y = "Frequency") +
  facet_wrap(~variable, scales = "free_x", nrow = 4) +
  theme(strip.text.x = element_text(size = 7),
        axis.text = element_text(size = 7.5),
        axis.title = element_text(size = 11)) +
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
  scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))
dev.off()

# Correlation Plot ###########################################################################
y3 <- y
names(y3) <- 1:41
corre <- cor(y3, method = "spearman")
library(corrplot)
whiteblack <- c("white", "black")
pdf(file = "cor_ant.pdf", width = 6, height = 5.5)
corrplot::corrplot(corre, 
                   # mar = c(0,0,0,0), 
                   tl.cex = 0.6, 
                   cl.cex = .7,
                   type = "upper",
                   diag = F,
                   outline = F)
# corrplot(corre, col = whiteblack, bg = "gold2")
dev.off()

# COrrelation text
# diag(corre) <- NA
# min(apply(corre, 2, min, na.rm = T))
# max(apply(corre, 2, max, na.rm = T))

# Generalized Dispersion Index ###########################################################################
source("~/Dropbox/SuperDisperso_AJS/Code/PaperCode.R")
GDI(y) # = 11.54345 Almost equi dispersed
gdi <- 11.54
# Standard Error = 0.9214711
sqrt(sigma_GDI(y)/nrow(ahs))
gdi_se <- .92
# Mean/Var ###########################################################################
# corre <- cor(y, method = "spearman")
# corre[lower.tri(corre, diag = T)] <- NA
mat <- cbind("Mean" = t(t(sapply(y, mean))),
             "Var" = t(t(sapply(y, var))))
mat <- cbind(mat, 
             "DI" = mat[,2]/mat[,1])
colnames(mat) <- c("Mean", "Variance", "DI")
mat #All variables are overdispersed
mat <- as.data.frame(mat)
mat$`GDI(SE)` <- c(rep(NA, (nrow(mat)/2)-1), "11.543(.92)",
                   rep(NA, (nrow(mat)/2)+1))
kable(mat, format = "latex", digits = 3, booktabs = T, escape = F,
      align = rep("c", ncol(mat)),
      label = "tabant",
      caption = "DESCRIPTIVE MEASUREMENTS FOR ANT DATA") %>% 
  kable_styling()

###################
# Se fosse para eu incluir Heavy tailed and Zero Inflation indices, how do I determine nu and phi? From the results of the model?
###################

# Descriptive statistics (Mean + Var + Correlation) ###########################################################################
corre <- cor(y, method = "spearman")
corre[lower.tri(corre, diag = T)] <- NA
mat <- cbind(corre, 
             "Mean" = t(t(sapply(y, mean))),
             "Var" = t(t(sapply(y, var))))
# colnames(mat)[6:7] <- c("Mean", "Variance")
# colnames(mat)[6:8]
options(knitr.kable.NA = '')
kable(mat, format = "latex", booktabs = T, digits = 3, escape = F,
      align = c("c"),
      caption = "Mean, variance and correlation for AHS response variables") %>% 
  kable_styling() %>%
  add_header_above(c(" "=1,"Spearman Correlation $\\\\rho$"=5,"Mean"=1,"Variance"=1),
                   escape = F)
