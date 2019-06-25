require(Hmisc)

# Importing data ###
# setwd("/home/guilherme/Google Drive/TCC/Simulacao")
setwd("/home/guilherme/Google Drive/TCC/ChJS")
data <- read.delim("/home/guilherme/Google Drive/TCC/Simulacao/concrete.txt", dec=".") # Available at http://archive.ics.uci.edu/ml/datasets/concrete+compressive+strength
names(data)[1] <- c("cement")
names(data)[9] <- c("concrete")
data <- data[,c(1,9)]
data$concrete <- sqrt(data$concrete)
# Descriptive Table ####
a <- with(data, rbind(cbind(min(cement), min(concrete)),
                      cbind(mean(cement), mean(concrete)),
                      cbind(median(cement), median(concrete)),
                      cbind(max(cement), max(concrete)),
                      cbind(sd(cement), sd(concrete)),
                      cor(cement,concrete)))
a <- round(a,2)
a <- as.data.frame(a)
a[6,2] <- ""
# Not used in the article ####
# shapiro.test(concrete)
# latex(a, file  = "", 
#       booktabs = TRUE,
#       title    = "Medidas descritivas",
#       rowname  = c("Mínimo","Média","Mediana","Máximo","Desvio padrão","Coeficiente de Correlação"),
#       where    = "H",
#       colheads = c("Cimento", "Concreto"), 
#       digits   = 2,
#       caption  = "Medidas resumo para os Domínios do SF36",
#       label    = "tab:descritivo"
# )
# x11()
# Hist + density ####
setEPS()
postscript("descritivo.eps",width=16)
par(mfrow=c(1,2), mar = c(5, 5, 4, 0.5))
hist(data$concrete, freq=F
     ,xlab = expression(sqrt("Concrete strength (MPa)"))
     ,ylab = "Density"
     ,cex.lab = 1.4
     ,cex.axis = 1
     ,las = 1
     ,main=""
     ,breaks = 10
     ,ylim = c(0,0.3)
)
box()
rug(data$concrete)
lines(density(data$concrete), lty = 1, lwd = 2)
curve(dnorm(x,mean(data$concrete),sd(data$concrete)), from = 0, to = 10, lty=2, add=T, lwd = 2)
legend(80,0.025
       , title = "Adjust"
       , legend = c("Nonparametric", "Parametric")
       , border = "n", col="black"
       , lty = c(1,2), lwd = 2, bty="n"
       , xjust = 1
)
# Scatterplot ####
plot(data$concrete~data$cement, data=data
     , las = 1
     , cex = 1
     , col = "black"
     , xlab = "Quantity of cement(kg)"
     , ylab = expression(sqrt("Concrete strength (MPa)"))
     , cex.lab = 1.5)
abline(coefficients(lm(data$concrete~data$cement)), lty=3, lwd=2)
loess_fit_pre <- predict(loess(data$concrete~data$cement))
lines(data$cement[order(data$cement)], loess_fit_pre[order(loess_fit_pre)], lwd = 2, lty = 7, type="l")

legend(80,85, legend = c("Nonparametric","Parametric")
       , border = "n", col="black"
       , lty = c(7,3), lwd = 2, bty="n"
       , xjust = 0
       , title = "Regression")
dev.off()
# Selecting a sample from NRSS and RSS ####
select <- function(aas){
  # Input:
  # Data from a SRS of size k
  # Description:
  # Select function return the position of the elements to be drawn from a SRS
  len <- length(aas)
  k <- sqrt(len)
  aas <- sort(aas)
  elementos <- rep(0,ncol=k)
  posicao <- rep(0,k)
  i <- 1
  if (len%%2 == 0){
    for(i in 1:k){
      if (i%%2 == 0){l  <- k/2}
      else{l <-  (k+2)/2 }
      posicao[i] <- (l+(i-1)*k)
      elementos[i] <- aas[posicao[i]]
    }
  }
  else{  
    for(i in 1:k){
      posicao[i] <- ((k+1)/2)+(i-1)*k
      elementos[i] <- aas[posicao[i]]
    }
  }
  return(posicao)
}
# Input for testing new code!!!!
# k2 <- 9        ## Tamanho da amostra inicial de um delineamento que fazem uso de ranqueamento
# l <- 25
# data <- data[1:100, ]
urss_process <- function(data, l, k2, po = T){
  dy <- matrix(data[,2], ncol=k2)
  
  if (po){ #Perfect Ordering
    dys <- t(apply(dy, 1, sort))
  } else{  #Imperfect Ordering
    dx <- matrix(data[,1], ncol=k2)
    dys <- matrix(rep(0,l*k2),ncol=k2)
    pos <- t(apply(dx,1,order))         # Every line is a sample
    for (linha in 1:nrow(dy)){                 # Order y by x
      dys[linha,] <- dy[linha,pos[linha,]]
    }
  }
  dys_amostras <- dys[,select(1:k2)]         
  return(dys_amostras)
}
#Input for testing
# m <- 1
# n <- 9
aco_process <- function(data, m, n, po = T){
  # data <- samples
  # m <- sl
  # k <- n
  # -----------------------------------------------------------------------
  # Split the data into two matrices
  # -----------------------------------------------------------------------
  
  x <- matrix(data[,1], n, m*n) # 
  y <- matrix(data[,2], n, m*n) # For k = 3, m*n = 75 and n = 3
  # dim(y)
  if (po){
    y.order <- apply(y, 2, sort)
    # dim(y.order)
    amostra.sel <- matrix(y.order[cbind(rep(1:n, 25), 1:(m*n))], byrow = T, nrow = m)
    
  } else{
    # Indices of the order of matrix x
    i <- apply(x,2,order, decreasing=F)
    
    # Crete indices to order y elements by x
    l <- i[1:nrow(data)]
    c <- rep(1:(n*m),each=n)
    
    y.order <- matrix(y[cbind(l,c)], n, n*m)
    
    # -----------------------------------------------------------------------
    # select only the diagonals
    # -----------------------------------------------------------------------
    ll <- rep(sequence(nrow(x)), m) # line
    cc <- rep(sequence(ncol(x))) # col
    
    #  diagonal
    amostra.sel <- matrix(y.order[cbind(ll,cc)], ncol = n, byrow = TRUE)
  }
  
  return(amostra.sel)
}
# Function to plot the control charts ####
plot_grafico <- function(k, delta, method, lc, lsc, lic, lsc1, lic1, lsc2, lic2, sl, dados, cor, cex, pch.forma, po){
  plot(dados
       , type ="b"
       , col  = cor
       , main =""
       , xlab =""
       , ylab = ""
       , ylim = c(lic-0.2, lsc+1)
       , las  = 1
       , pch  = pch.forma
       , xaxt = "n"
       , cex = cex
       , cex.axis = cex/1.25
  )
  # Axis
  axis(1, at = seq(1, sl+75,2), cex.axis = cex/1.25)
  #Divisão de amostras
  abline(v=sl, lty = "dotted")
  #Method
  text(31,lsc+0.5,paste0(method, ifelse(method == "SRS", "", ifelse(po, "-PR", "-IR"))), cex = cex, font = 2)
  ## Cima
  #LSC
  abline(h=lsc, col = "black", lty="longdash", lwd = cex*1.25)
  mtext("UCL", at=lsc, side=4,line=0.2, adj=0, las=1, cex = cex/1.6)
  #2sigma
  abline(h=lsc2, col = "black", lty ="dashed", lwd = cex)
  mtext(substitute(paste(2, hat(sigma)[bar(X)[method]])
                   , list(method = method))
        , at=lsc2, side=4,line=0.2, adj=0, las=1, cex = cex/1.2)
  #1sigma
  abline(h=lsc1, col = "black", lty ="twodash", lwd = 1.75)
  mtext(substitute(paste(1, hat(sigma)[bar(X)[method]])
                   , list(method = method))
        , at=lsc1, side=4,line=0.2, adj=0, las=1, ann=F, cex = cex/1.2)
  #Central
  abline(h=lc)
  mtext("CL", at=lc, side=4,line=0.2, adj=0, las=1, cex = cex/1.6)
  ## Baixo
  #LCL
  abline(h=lic, col = "black", lty="longdash", lwd = cex*1.25)
  mtext("LCL", at=lic, side=4,line=0.2, adj=0, las=1, outer=F, cex = cex/1.6)
  #2sigma
  abline(h=lic2, col = "black", lty ="dashed", lwd = cex)
  mtext(substitute(paste(2, hat(sigma)[bar(X)[method]])
                   , list(method=method))
        , at=lic2, side=4,line=0.2, adj=0, las=1, cex = cex/1.2)
  #1sigma
  abline(h=lic1, col = "black", lty ="twodash", lwd = 1.75)
  mtext(substitute(paste(1, hat(sigma)[bar(X)[method]])
                   , list(method=method))
        , at=lic1, side=4,line=0.2, adj=0, las=1, cex = cex/1.2)
}
# Função que engloba tudo ####
graficos_controle <- function(data, sl, k, method, snl, A, desc, delta, sig0, pert, cex, po){
  if (method == "NRSS" | method == "RSS"){
    k2 <- k^2
  } else if (method == "SRS"){
    k2 <- k
  } else {stop("A função suporta apenas os delineamentos NRSS, RSS e AAS")}
  
  samples <- data[sample(1:nrow(data), size = sl*k2, replace = TRUE), ]
  
  # nova média do processo, caso entre em descontrole
  mu <- (delta*sig0)/sqrt(k)
  
  if (method == "NRSS"){
    # Faz todo o processo de ordenação e seleção de amostra
    amostras <- urss_process(samples, sl, k2, po)
    amostras_mean <- rowMeans(amostras)
    # Linha Central
    mean_mean <- mean(amostras)
    mcov <- cov(amostras)
    vari <- diag(mcov)
    mcov[lower.tri(mcov, diag=TRUE)] <- 0
    # Variância da média estimada URSS
    var_mean <- (1/k2)*sum(vari)+(2/k2)*sum(colSums(mcov))
    # Limites de controle
    lic <- mean_mean - A*sqrt(var_mean)
    lc <- mean_mean
    lsc <- mean_mean + A*sqrt(var_mean)
    # 1 sigma e 2 sigmas
    lic1 <- mean_mean - 1*sqrt(var_mean)
    lsc1 <- mean_mean + 1*sqrt(var_mean)
    lic2 <- mean_mean - 2*sqrt(var_mean)
    lsc2 <- mean_mean + 2*sqrt(var_mean)
    
    # Coletando novas amostras do processo 
    if (desc){ # Sob descontrole
      data$concrete <- data$concrete + rnorm(nrow(data), mu, sd = pert)
      newdata <- data[sample(1:nrow(data), size = snl*k2, replace = TRUE), ]  # Desloca a média do processo com a pertubação
      # Amostra URSS
      amostras_new <- urss_process(newdata, snl, k2, po)
      mean_amostras_new <- rowMeans(amostras_new)  
    }
    else{ # Sob controle
      newdata <- data[sample(1:nrow(data), size = snl*k2, replace = TRUE), ] 
      # Amostra URSS
      amostras_new <- urss_process(newdata, snl, k2, po)
      mean_amostras_new <- rowMeans(amostras_new)
    }
  } 
  else if (method == "RSS"){
    amostras <- aco_process(samples, sl, k, po)
    amostras_mean <- rowMeans(amostras)
    # Linha Central
    mean_mean <- mean(amostras)
    # Variancia
    variance <- var(as.numeric(amostras))
    parte1 <- variance/k
    parte2 <- sum((colMeans(amostras) - mean_mean)^2)/(k^2)
    # Variância da média estimada ACO
    var_mean <- parte1 - parte2
    # Limites de controle
    lic <- mean_mean - A*sqrt(var_mean)
    lc <- mean_mean
    lsc <- mean_mean + A*sqrt(var_mean)
    # 1 e 2 sigmas
    lic1 <- mean_mean - 1*sqrt(var_mean)
    lsc1 <- mean_mean + 1*sqrt(var_mean)
    lic2 <- mean_mean - 2*sqrt(var_mean)
    lsc2 <- mean_mean + 2*sqrt(var_mean)
    
    # Coletando novas amostras do processo com descontrole
    if (desc) { # Sob descontrole
      data$concrete <- data$concrete + rnorm(nrow(data), mean = mu, sd = pert)
      newdata <- data[sample(1:nrow(data), size = snl*k2, replace = TRUE), ]
      # Amostra ACO
      amostras_new <- aco_process(newdata, snl, k, po)
      mean_amostras_new <- rowMeans(amostras_new)
    }
    else { # Sob controle
      # Coletando novas amostras do processo sem descontrole
      newdata <- data[sample(1:nrow(data), size = snl*k2, replace = TRUE), ]
      # Amostra ACO
      amostras_new <- aco_process(newdata, snl, k, po)
      mean_amostras_new <- rowMeans(amostras_new)
    }
    
  } else{ # AAS
    # Estimando os Limites de controle
    samples <- matrix(samples[,2], ncol = k2)
    amostras_mean <- rowMeans(samples)
    # Linha Central estimada
    mean_mean <- mean(samples)
    # Variancia
    sd <- mean(apply(samples, 1, sd))
    if (k == 3){
      c4 <- 0.8862
    }else{ #k = 5
      c4 <- 0.94
    }
    # Variância da média estimada AAS
    var_mean <- (A*sd)/(c4*sqrt(k2))
    var_mean1 <- (1*sd)/(c4*sqrt(k2))
    var_mean2 <- (2*sd)/(c4*sqrt(k2))
    # Limites de controle
    lic <- mean_mean - var_mean
    lc <- mean_mean
    lsc <- mean_mean + var_mean
    # 1 e 2 sigmas
    lic1 <- mean_mean - var_mean1
    lsc1 <- mean_mean + var_mean1
    lic2 <- mean_mean - var_mean2
    lsc2 <- mean_mean + var_mean2
    
    # Coletando novas amostras do processo
    if (desc){ # Sob descontrole
      data$concrete <- data$concrete + rnorm(nrow(data), mean = mu, sd = pert)
      newdata <- data[sample(1:nrow(data), size = snl*k2, replace = TRUE),2]
      amostras_new <- matrix(newdata, ncol = k2)
      # Amostra AAS
      mean_amostras_new <- rowMeans(amostras_new)  
    }else{    # Sob controle
      newdata <- data[sample(1:nrow(data), size = snl*k2, replace = TRUE),2]
      amostras_new <- matrix(newdata, ncol = k2)
      # Amostra AAS
      mean_amostras_new <- rowMeans(amostras_new) }}
  data_ploting <- as.data.frame(c(amostras_mean,mean_amostras_new))
  data_ploting$pch.forma <- ifelse(as.numeric(data_ploting<lic | data_ploting>lsc),8,
                                   ifelse(as.numeric(data_ploting<lic2 | data_ploting>lsc2),1,16))
  # if(method=="SRS"){
  #   print(var_mean)
  # }else{
  #   print(A*sqrt(var_mean))
  # }
  plot_grafico(k, delta, method, lc, lsc, lic, lsc1, lic1, lsc2, lic2, sl, data_ploting[,1],"black", cex, data_ploting[,2], po)
  # return(data_ploting) # Utilizado para plotar os dados do gráfico e fazer o teste
}


# COLOCANDO ORDENAÇÃO PERFEITA

### sqrt
#3k3delta0   -> j = 4
#5k5delta0   -> j = 4
#3k3delta1.2 -> j = 17  
#5k5delta1.2 -> j = 7

# Inputs ####
k2 <- 25       ## Tamanho da amostra inicial de um delineamento que fazem uso de ranqueamento
k <- 5         ## Tamanho da amostra final
sl <- 25       ## Amostras utilizadas para calcular os limites de controle
snl <- 75      ## Amostras que não serão utilizadas para calcular o lic e lsc
A <- 3         ## Amplitude
desc <- T      ## Descontrole
delta <- 1.2   ## Delta de Descontrole
mu0 <- 5.81    ## Média mu0
sig0 <- 1.45   ## Sigma sob controle
pert <- 0.17   ## Pertubação nos dados
cex <- 1.8     ## Tamanho da fonte
# po <- T        ## Ordenção perfeita

counter <- 0


for(j in 7:7){
  setEPS()
  postscript(paste0("k",k,"d1.2j",j,"_sqrt",".eps"),width=16, height = 15)
  par(mfrow=c(5,1), oma=c(5,0.5,2,0.5),mai=c(0.4,1,0.15,1))
  method <- c("SRS", "RSS", "NRSS")
  for (i in 1:3){
    # counter <- counter + 1
    if (method[i] == "SRS"){
      RNGkind(sample.kind = "Rounding")
      set.seed(j)
      graficos_controle(data = data, sl = sl, k = k, method = method[i],snl = snl, po = F,
                        A = A, desc = desc, delta = delta, sig0 = sig0, pert = pert, cex = cex)
    } else{
      RNGkind(sample.kind = "Rounding")
      set.seed(j)
      graficos_controle(data = data, sl = sl, k = k, method = method[i],snl = snl, po = T, 
                        A = A, desc = desc, delta = delta, sig0 = sig0, pert = pert, cex = cex)
      RNGkind(sample.kind = "Rounding")
      set.seed(j)
      graficos_controle(data = data, sl = sl, k = k, method = method[i],snl = snl, po = F,
                        A = A, desc = desc, delta = delta, sig0 = sig0, pert = pert, cex = cex)
    }
    if(i == 2){ # Trick
      mtext(expression(sqrt("Concrete strength (MPa)")),side=2,line=3, adj=1, cex = cex/1.25)
    }
      
    # print(counter)
  }
  mtext("Samples to estimate LCL and UCL",side=1,line=3, adj=0, cex = cex/1.25)
  mtext("Samples to monitorate the process",side=1,line=3, adj=1, cex = cex/1.25)
  mtext("Samples order",side=1, line=5, cex = cex/1.25)
  dev.off()
}


# Sementes Utilizadas:

#3k3delta0   -> j = 5
#5k5delta0   -> j = 3
#3k3delta1.2 -> j = 5468
#5k5delta1.2 -> j = 5466

# Determinig the pattern of the chart for imperfect ordering

# # Testes
# for(j in 1:200){
#   method <- c("RSS")
#     set.seed(j)
#     a <- graficos_controle(data = data, sl = sl, k = k, method = method,snl = snl, A = A, desc = desc, delta = delta, sig0 = sig0, pert = pert, cex = cex, po = T)
#     if (j==1){
#       a1 <- a
#     } else{
#       a1 <- rbind(a1, a)
#     }
# }
# nrow(a1)
# srs <- a1
# rss <- a1
# nrss <- a1
# delin.k3 <- list(srs, rss, nrss)
# lapply(delin.k3, function(x) nrow(x[x$pch.forma==8, ]))
# 
# srs <- a1
# rss <- a1
# nrss <- a1
# delin.k5 <- list(srs, rss, nrss)
# rf <- cbind(do.call(rbind.data.frame, lapply(delin.k3, function(x) nrow(x[x$pch.forma==8, ]))),
# do.call(rbind.data.frame, lapply(delin.k5, function(x) nrow(x[x$pch.forma==8, ]))))
# names(rf) <- c("k=3","k=5")
# rownames(rf) <- c("srs","rss","nrss")
