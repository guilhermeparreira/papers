#Testando opcoes graficas####
require(Hmisc)
require(lattice)
require(qcc)

# Importando os dados ####
setwd("/home/guilherme/Google Drive/TCC/Simulacao")
data <- read.delim("concrete.txt", dec=".")
str(data)
names(data)[1] <- c("cement")
names(data)[9] <- c("concrete")
data <- data[,c(1,9)]
attach(data)
# Tabela Descritiva ####
a <- rbind(cbind(min(cement), min(concrete)),
           cbind(mean(cement), mean(concrete)),
           cbind(median(cement), median(concrete)),
           cbind(max(cement), max(concrete)),
           cbind(sd(cement), sd(concrete)),
           cor(cement,concrete))
a <- round(a,2)
a <- as.data.frame(a)
a[6,2] <- ""
shapiro.test(concrete)
latex(a, file  = "", 
      booktabs = TRUE,
      title    = "Medidas descritivas",
      rowname  = c("Mínimo","Média","Mediana","Máximo","Desvio padrão","Coeficiente de Correlação"),
      where    = "H",
      colheads = c("Cimento", "Concreto"), 
      digits   = 2,
      caption  = "Medidas resumo para os Domínios do SF36",
      label    = "tab:descritivo"
)
x11()
#Gráficos ####
# Primeiro o Histograma
par(mfrow=c(1,2))
hist(concrete, freq=F
     ,main = "Histograma para a variável resistência de concreto"
     ,xlab = "Resistência de concreto(MPa)"
     ,ylab = "Densidade"
     ,cex.lab = 1.4
     ,cex.axis = 1
     ,las = 1
)
box()
rug(concrete)
lines(density(concrete),col = "darkblue", lwd = 2)
curve(dnorm(x,mean(concrete),sd(concrete)), from = -10, to = 100, col="red", add=T, lwd = 2)
legend(103,0.025
       , title = "Ajuste"
       , legend = c("Não Paramétrico", "Paramétrico")
       , border = "n", col=c("darkblue","red")
       , lty = 1, lwd = 2, bty="n"
       , xjust = 1
)
# Gráfico de correlação
plot(concrete~cement, data=data
     , las = 1
     , cex = 1
     , col = "blue"
     , xlab = "Quantidade de cimento(kg)"
     , ylab = "Resistência de concreto(MPa)"
     , main = "Gráfico de dispersão para as variáveis concreto e cimento"
     , cex.lab = 1.5)
abline(coefficients(lm(concrete~cement)), col="darkorange", lwd=2)
loess_fit_pre <- predict(loess(concrete~cement))
lines(cement[order(cement)], loess_fit_pre[order(loess_fit_pre)], lwd = 2, col = "slategray4", type="l")

legend(80,85, legend = c("Não-Paramétrica", "Paramétrica")
       , border = "n", col=c("slategray4","darkorange")
       , lty = 1, lwd = 2, bty="n"
       , xjust = 0
       , title = "Regressão")
detach(data)
# Dados teste ####
data(pistonrings)
attach(pistonrings)
par(bg = 'white')
q1 <- qcc(pistonrings[1:100,1], type="b", col="blue", main=paste("Gráficos de controle para a média via",method, ",k=",k, "delta=",delta))
abline(h=60)
text(1,50,"1sig")
par(mfrow=c(3,1))
q1<-qcc(pistonrings[1:100,1]
        ,chart.all = F
        ,text=" "
        , type="xbar.one"
        ,add.stats=F
        ,par.restore=F
        ,title=""
)
x11()
oldpar <- par(no.readonly = TRUE)
par(mfrow=c(2,1))
plot(q1,chart.all = F
     ,text=" "
     , type="xbar.one"
     ,add.stats=F
     ,par.restore=F
     ,title="")
par(oldpar)

par(bg="white")
oldpar <- par(no.readonly = TRUE)
par(mfrow=c(3,1), mar=c(0,0.1,0,0.1), oma=c(0, 1, 0, 1),
    mai = c(0,0,0,0))

#c(bottom, left, top, right)
#oma é o de fora, mai o de dentro

plot(pistonrings[1:100,1], title="First samples", add.stats=FALSE, restore.par=FALSE, ylab="", xaxt=NULL, xaxt="n")
par(mfrow=c(1,1))
mtext(expression(paste(2, hat(sigma)[bar(X)])), at=74, side=4,line=0.2, adj=0, las=1)

?par


plot(q1, title="Second samples", add.stats=FALSE, restore.par=FALSE, ylab="", xlab=NULL, xaxt="n")
plot(q3, title="Second samples", add.stats=FALSE, restore.par=FALSE, ylab="", xlab=NULL)
?par
par(oldpar)
mar
(bottom, left, top, right)


mtext("Amostras para estimar os LIC e LSC",side=1,line=0, adj=0)
mtext("Amostras para monitorar o processo",side=1,line=0, adj=1)
mtext("Identificação da amostra",side=1, line=1)

mtext("Amostras para estimar os LIC e LSC",side=3,line=-1.8, adj=0)
mtext("Amostras para monitorar o processo",side=3,line=-1.8, adj=1)

par(mfrow=c(3,1), mai= c(0.2, 1, 0.4, 1), oma=c(1.4, 1, 1, 1))
c(bottom, left, top, right)
par(bg="white")  
method <- "AAS"
delta <- 0.4
k <- 3

sl <- 25
k <- 3
snl <- 75
A <- 3
method <- "AAS"
desc <- T
delta <- 0.4
mu0 <- 35.8
sig0 <- 16.7
pert <- 2
dados <- pistonrings[1:100,1]
lsc <- 74.03
lic <- 73.97
lsc1 <- 74.02
lic1 <- 74.00
lsc2 <- 74.01
lic2 <- 74.005
lc <- 76
# Processos de seleção de amostras da URSS e ACO ####
urss_process <- function(data, l, k2){
  dx <- matrix(data[,1], ncol=k2)
  dy <- matrix(data[,2], ncol=k2)
  
  pos <- t(apply(dx,1,order))                # Posição dos elementos da matriz x (var. auxiliar)
  dys <- matrix(rep(0,l*k2),ncol=k2)         # Matriz vazia para receber as amostras ordenadas
  for (linha in 1:nrow(dy)){                 # Ordena a matriz y por x
    dys[linha,] <- dy[linha,pos[linha,]]
  }
  dys_amostras <- dys[,select(1:k2)]                    # Matriz ordenada com as observações necessárias
  return(dys_amostras)
}
aco_process <- function(base, m, n){
  
  # -----------------------------------------------------------------------
  # separando os dados em duas matrizes
  # -----------------------------------------------------------------------
  x <- matrix(base[,1], n, m*n) # primeiro vetor
  y <- matrix(base[,2], n, m*n) # segundo vetor
  
  # indices da ordem da amostra x
  i <- apply(x,2,order, decreasing=F)
  
  # criando indices para ordenar os elementos de y com base em x
  l <- i[1:nrow(base)]
  c <- rep(1:(n*m),each=n)
  
  y.order <- matrix(y[cbind(l,c)], n, n*m)
  
  # -----------------------------------------------------------------------
  # selecao das diagonais
  # -----------------------------------------------------------------------
  ll <- rep(sequence(nrow(x)), m) # linha selecionada
  cc <- rep(sequence(ncol(x))) # coluna selecionada
  
  # selecao do elementos na diagonal
  amostra.sel <- matrix(y.order[cbind(ll,cc)], ncol = n, byrow = TRUE)
  return(amostra.sel)
}
#Função para plotar o gráfico com legendas e tudo mais ####
plot_grafico <- function(k, delta, method, lc, lsc, lic, lsc1, lic1, lsc2, lic2, sl, dados, cor, cex){
  plot(dados
       , type ="b"
       , col  = cor
       , main =""
       , xlab =""
       , ylab = ""
       , ylim = c(lic-13, lsc+13)
       , las  = 1
       , pch  = 16
       , xaxt = "n"
       , cex = cex
       , cex.axis = cex/1.25
  )
  # Axis
  axis(1, at = seq(1, sl+75,2), cex.axis = cex/1.25)
  #Divisão de amostras
  abline(v=sl, lty = "dotted")
  #Method
  text(31,lsc+10,method, cex = cex, font = 2)
  ## Cima
  #LSC
  abline(h=lsc, col = "red", lty="longdash", lwd = cex*1.25)
  mtext("LSC", at=lsc, side=4,line=0.2, adj=0, las=1, cex = cex/1.1)
  #2sigma
  abline(h=lsc2, col = "black", lty ="dashed", lwd = cex)
  mtext(substitute(paste(2, hat(sigma)[bar(X)[method]])
                   , list(method = method))
        , at=lsc2, side=4,line=0.2, adj=0, las=1, cex = cex/1.25)
  #1sigma
  abline(h=lsc1, col = "grey", lty ="twodash", lwd = 1.75)
  mtext(substitute(paste(1, hat(sigma)[bar(X)[method]])
                   , list(method = method))
        , at=lsc1, side=4,line=0.2, adj=0, las=1, ann=F, cex = cex/1.25)
  #Central
  abline(h=lc)
  mtext("LC", at=lc, side=4,line=0.2, adj=0, las=1, cex = cex/1.1)
  ## Baixo
  #LIC
  abline(h=lic, col = "red", lty="longdash", lwd = cex*1.25)
  mtext("LIC", at=lic, side=4,line=0.2, adj=0, las=1, outer=F, cex = cex/1.1)
  #2sigma
  abline(h=lic2, col = "black", lty ="dashed", lwd = cex)
  mtext(substitute(paste(2, hat(sigma)[bar(X)[method]])
                   , list(method=method))
        , at=lic2, side=4,line=0.2, adj=0, las=1, cex = cex/1.25)
  #1sigma
  abline(h=lic1, col = "grey", lty ="twodash", lwd = 1.75)
  mtext(substitute(paste(1, hat(sigma)[bar(X)[method]])
                   , list(method=method))
        , at=lic1, side=4,line=0.2, adj=0, las=1, cex = cex/1.25)
}
# Função que engloba tudo ####
graficos_controle <- function(data, sl, k, method, snl, A, desc, delta, sig0, pert, cex){
  if (method == "URSS" | method == "ACO"){
    k2 <- k^2
  } else if (method == "AAS"){
    k2 <- k
  } else {stop("A função suporta apenas os delineamentos URSS, ACO e AAS")}
  
  samples <- data[sample(1:nrow(data), size = sl*k2, replace = TRUE), ]

  # nova média do processo, caso entre em descontrole
  mu <- (delta*sig0)/sqrt(k)
  
  if (method == "URSS"){
    # Faz todo o processo de ordenação e seleção de amostra
    amostras <- urss_process(samples, sl, k2)
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
      amostras_new <- urss_process(newdata, snl, k2)
      mean_amostras_new <- rowMeans(amostras_new)  
    }
    else{ # Sob controle
      newdata <- data[sample(1:nrow(data), size = snl*k2, replace = TRUE), ] 
      # Amostra URSS
      amostras_new <- urss_process(newdata, snl, k2)
      mean_amostras_new <- rowMeans(amostras_new)
    }
  } 
  else if (method == "ACO"){
    amostras <- aco_process(samples, sl, k)
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
      amostras_new <- aco_process(newdata, snl, k)
      mean_amostras_new <- rowMeans(amostras_new)
    }
    else { # Sob controle
      # Coletando novas amostras do processo sem descontrole
      newdata <- data[sample(1:nrow(data), size = snl*k2, replace = TRUE), ]
      # Amostra ACO
      amostras_new <- aco_process(newdata, snl, k)
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
  data_ploting$color <- ifelse(as.numeric(data_ploting<lic | data_ploting>lsc),"red",
                               ifelse(as.numeric(data_ploting<lic2 | data_ploting>lsc2),"gold","blue"))
  if(method=="AAS"){
    print(var_mean)
  }else{
    print(A*sqrt(var_mean))
  }
  plot_grafico(k, delta, method, lc, lsc, lic, lsc1, lic1, lsc2, lic2, sl, data_ploting[,1],data_ploting[,2], cex)
}
# Inputs ####
k2 <- 25
sl <- 25
k <- 5
snl <- 75
A <- 3
desc <- T
delta <- 1.2
mu0 <- 35.8
sig0 <- 16.7
pert <- 2
cex <- 2.6
#c(bottom, left, top, right)
#oma é o de fora, mai o de dentro
par(mfrow=c(3,1), oma=c(5,1,2,2),mai=c(0.5,1,0.5,1))
method <- c("AAS", "ACO", "URSS")
for (i in 1:3){
  graficos_controle(data = data, sl = sl, k = k, method = method[i],snl = snl, A = A, desc = desc, delta = delta, sig0 = sig0, pert = pert, cex = cex)
  if(i == 2){
    mtext("Resistência do concreto(MPa)",side=2,line=3, adj=1, cex = cex)
  }else if (i == 3){
    title(substitute(paste("Gráficos de controle para a média com k=",k," e ",delta,"=",delta_value)
                     ,list(k = k, delta_value=delta))
          ,outer=T,line=-1, cex.main=cex*1.5)
    mtext("Amostras para estimar os LIC e LSC",side=1,line=3, adj=0, cex = cex/1.35)
    mtext("Amostras para monitorar o processo",side=1,line=3, adj=1, cex = cex/1.35)
    mtext("Identificação da amostra",side=1, line=5, cex = cex/1.25)
  }
}

dados_ex <- as.data.frame(c(rnorm(25,sd=0.8),rnorm(75,1.5,sd=0.8)))
lf <- mean(dados_ex[26:100,1])
par(mfrow=c(1,1), mai=c(0.5,2,0.5,1))
plot_grafico_2(lf,0, 3, -3, 1, -1, 2, -2, 25, dados_ex, 2)
plot_grafico_2 <- function(lf,lc, lsc, lic, lsc1, lic1, lsc2, lic2, sl, dados, cex){
  dados$color <- ifelse(as.numeric(dados[,1]<lic | dados[,1]>lsc),"red",
                               ifelse(as.numeric(dados[,1]<lic2 | dados[,1]>lsc2),"gold","blue"))
  plot(dados[,1]
       , type ="b"
       , col  = dados[,2]
       , main =""
       , xlab =""
       , ylab = ""
       , ylim = c(lic-1, lsc+1)
       , las  = 1
       , pch  = 16
       , xaxt = "n"
       , cex = cex
       , cex.axis = cex/1.25
       , yaxt = "n"
  )
  # Lado esquerdo
  mtext(expression(mu[0] + frac(3*sigma[0],sqrt(k))), at = lsc, side = 2, las = 1, cex = cex/1.1, line = 2)
  mtext(expression(mu[0] - frac(3*sigma[0],sqrt(k))), at = lic, side = 2, las = 1, cex = cex/1.1, line = 2)
  mtext(expression(mu[0]), at = lc, side = 2, las = 1, cex = cex/1.1, line=2)
  mtext(expression(mu[1]), at = lf, side = 2, las = 1, cex = cex/1.1, line=2)
  abline(h=lf, lwd=2*cex)
  
  mtext("Amostras para estimar os LIC e LSC",side=1,line=3, adj=0, cex = cex/1.35)
  mtext("Amostras para monitorar o processo",side=1,line=3, adj=1, cex = cex/1.35)
  mtext("Identificação da amostra",side=1, line=5, cex = cex/1.25)
  
  # Axis
  axis(1, at = seq(1, sl+75,2), cex.axis = cex/1.25)
  #Divisão de amostras
  abline(v=sl, lty = "dotted")
  ## Cima
  #LSC
  abline(h=lsc, col = "red", lty="longdash", lwd = cex*1.25)
  mtext("LSC", at=lsc, side=4,line=0.2, adj=0, las=1, cex = cex/1.1)
  #2sigma
  abline(h=lsc2, col = "black", lty ="dashed", lwd = cex)
  mtext(substitute(paste(2, sigma[bar(Y)[method]])
                   , list(method = method))
        , at=lsc2, side=4,line=0.2, adj=0, las=1, cex = cex/1.25)
  #1sigma
  abline(h=lsc1, col = "grey", lty ="twodash", lwd = 1.75)
  mtext(substitute(paste(1, sigma[bar(Y)[method]])
                   , list(method = method))
        , at=lsc1, side=4,line=0.2, adj=0, las=1, ann=F, cex = cex/1.25)
  #Central
  abline(h=lc)
  mtext("LC", at=lc, side=4,line=0.2, adj=0, las=1, cex = cex/1.1)
  ## Baixo
  #LIC
  abline(h=lic, col = "red", lty="longdash", lwd = cex*1.25)
  mtext("LIC", at=lic, side=4,line=0.2, adj=0, las=1, outer=F, cex = cex/1.1)
  #2sigma
  abline(h=lic2, col = "black", lty ="dashed", lwd = cex)
  mtext(substitute(paste(2, sigma[bar(Y)[method]])
                   , list(method=method))
        , at=lic2, side=4,line=0.2, adj=0, las=1, cex = cex/1.25)
  #1sigma
  abline(h=lic1, col = "grey", lty ="twodash", lwd = 1.75)
  mtext(substitute(paste(1, sigma[bar(Y)[method]])
                   , list(method=method))
        , at=lic1, side=4,line=0.2, adj=0, las=1, cex = cex/1.25)
}
