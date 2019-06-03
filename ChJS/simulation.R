select <- function(aas){
  # Arguments:
    # Data from a SRS (Sample Size k2)
  # Description
    # This function return the sample positions of a NRSS to be selected from a SRS
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
#Example
a <- rnorm(16)
select(a)

# Probability Density function(PDF) from a ordered sample ####
FdpOrdem <- function(k2, i, x){
  # Arguments:
    # k2 is the size of your initial sample k2
    # i is the ordered position of your sample (result from select function)
    # x is the variable to be integrated (using integrate function)
  # Description:
    # Probability density function from a ordered sample
  (factorial(k2)*(pnorm(x)^(i-1))*((1-pnorm(x))^(k2-i))*dnorm(x))/(factorial(i-1)*factorial(k2-i))
}

#fdpx (x*PDF - expected value)
FdpOrdemX <- function(k2, i, x){
  # Arguments:
    # k2 is the size of your initial sample k2
    # i is the ordered position of your sample (result from select function)
    # x is the variable to be integrated (using integrate function)
  # Description:
    # PDF from a ordered sample multiplied by "x", to be used in the var calculation
  ((x*factorial(k2))*(pnorm(x)^(i-1))*((1-pnorm(x))^(k2-i))*dnorm(x))/(factorial(i-1)*factorial(k2-i))
  
}

# ((x^2)*PDF )
FdpOrdemX2 <- function(x,i,k2){
  # Arguments:
    # k2 is the size of your initial sample k2
    # i is the ordered position of your sample (result from select function)
    # x is the variable to be integrated (using integrate function)
  # Description:
    # PDF from a ordered sample multiplied by "x^2", to be used in the var calculation
  (((x^2)*factorial(k2))*(pnorm(x)^(i-1))*((1-pnorm(x))^(k2-i))*dnorm(x))/(factorial(i-1)*factorial(k2-i))
  
}

# Calculating the Variance
k2 <- 25
i <- 8
ex2 <- integrate(FdpOrdemX2,lower=-Inf,upper=+Inf,k2=k2,i=i)$value
ex  <- (integrate(FdpOrdemX,lower=-Inf,upper=+Inf,k2=k2,i=i)$value)^2
var_i_8 <- ex2 - ex

k2 <- 25
i <- 3
ex2 <- integrate(FdpOrdemX2,lower=-Inf,upper=+Inf,k2=k2,i=i)$value
ex  <- (integrate(FdpOrdemX,lower=-Inf,upper=+Inf,k2=k2,i=i)$value)^2
var_i_3 <- ex2 - ex


# Covariance ####
cov.ord <- function(k2, i, j, x, y){
  # Arguments:
    # k2 is the size of your initial sample k2
    # i is the position of the first sample element
    # j is the position of the second sample element
    # x is the variable to be integrated corresponding to the element i
    # y is the variable to be integrated corresponding to the element j
  # Description:
    # It calculates the covariance between two elements from a ordered sample
  exy <- integrate(function(y) { 
    sapply(y, function(y) {
      integrate(function(x) (x*y*factorial(k2)*(pnorm(x)^(i-1))*((pnorm(y)-pnorm(x))^(j-i-1))*((1-pnorm(y))^(k2-j))*dnorm(x)*dnorm(y))
                /(factorial(i-1)*factorial(j-i-1)*factorial(k2-j)), lower =  -Inf, upper = y)$value
    })
  }, lower = -Inf, upper = +Inf)$value
  ex <- (integrate(FdpOrdemX,lower=-Inf,upper=+Inf,k2=k2,i=i)$value)
  ey <- (integrate(FdpOrdemX,lower=-Inf,upper=+Inf,k2=k2,i=j)$value)
  return(exy-(ex*ey))
}

# Example
correlacao = cov.ord(k2=25,i=3,j=8,x,y)/(sqrt(var_i_8)*sqrt(var_i_3))

# Variance of mean from URSS - under perfect ranking ####
# It is necessary to input the elements positions that you will sample
var_urss2 <- function(urss){
  # Arguments:
    # Sample positions of your sample (TIP: use select())
  # Description:
    # Based on elements positions - this function return the variance of mean from URSS
  k <- length(urss)
  k2 <- k^2
  cov <- matrix(0,nrow=k, ncol=k)
  for(i in 1:k){
    for(j in 1:k){
      if (j > i){
        cov[i,j] <- cov.ord(k=k2, i=urss[i], j=urss[j], x,y) # Covariance
      } else if(i == j){
          cov[i,j] <- (integrate(FdpOrdemX2,lower=-Inf,upper=+Inf,k=k2,i=urss[i])$value)-(integrate(FdpOrdemX ,lower=-Inf,upper=+Inf,k=k2,i=urss[i])$value)^2 #Variance
      }
          else {
            cov[i,j] <- 0
        }
      }
  }
  var <- diag(cov)
  cov[lower.tri(cov, diag=TRUE)] <- 0
  return((1/k2)*sum(var)+(2/k2)*sum(colSums(cov)))
}

# Example
a <- select(rnorm(9))
var_urss2(a)

# Estimating the covariance function ####
# Covariance example
cov325 <- cov.ord(k2=9,i=2,j=5,x,y)
cov328 <- cov.ord(k2=9,i=2,j=8,x,y)
cov358 <- cov.ord(k2=9,i=5,j=8,x,y)

# Variance example
k2 <- 9
i <- 8
ex2 <- integrate(FdpOrdemX2,lower=-Inf,upper=+Inf,k2=k2,i=i)$value
ex  <- (integrate(FdpOrdemX,lower=-Inf,upper=+Inf,k2=k2,i=i)$value)^2
var38 <- ex2 - ex

k2 <- 9
i <- 5
ex2 <- integrate(FdpOrdemX2,lower=-Inf,upper=+Inf,k2=k2,i=i)$value
ex  <- (integrate(FdpOrdemX,lower=-Inf,upper=+Inf,k2=k2,i=i)$value)^2
var35 <- ex2 - ex

k2 <- 9
i <- 2
ex2 <- integrate(FdpOrdemX2,lower=-Inf,upper=+Inf,k2=k2,i=i)$value
ex  <- (integrate(FdpOrdemX,lower=-Inf,upper=+Inf,k2=k2,i=i)$value)^2
var32 <- ex2 - ex

# Result for Covariance
((1/9)*sum(var32,var38,var35))+(2/9)*sum(cov358,cov328,cov325)

#Plot
# values <- c(2,3,4,5,6,7,8,9,10,11,12)
# values <- values^2
# var_mean <- rep(0,length(values))
# for(i in 1: length(values)){
#   var_mean[i] <- var_urss2(select(rnorm(values[i])))
# }
# plot(values,var_mean,xaxt = 'n',type="b", xlab="Valores correspondentes ? AAS", ylab="Vari?ncia da m?dia da URSS",main="Valor da vari?ncia da m?dia URSS versus tamanho da amostra AAS")
# axis(side = 1, at = values,labels = T)
# memory.limit(size=8000)
# memory.size()

# ARL - Perfect ranking ####
aas <- c(9, 16, 25, 36) 
delta <- c(0, 0.1, 0.2, 0.3, 0.4, 0.8, 1.2, 1.6, 2, 2.4, 3.2) 

cms <- function(aas, delta){
  # Argument:
    # aaa: corresponds to the size of your SRS (k2)
  # Description:
    # It calculates ARL for perfect ranking based on NRSS using 3sigma limits for 1,000,000 samples

la <- length(aas) 
ls <- length(delta) 
# Matrix to receive CMS's from different shifts in the process and sample sizes
cms_matrix <- matrix(0, nrow = ls, ncol = la, dimname=list(c(delta),c(sqrt(aas)))) 

for (a in 1:length(aas)){                                             # For every sample size
  k          <- select(rnorm(aas[a])) 
  lk         <- length(k) 
  lk2        <- lk^2 
  var_mean_u <- var_urss2(k)                                          # calculacalcula the variance of the mean which relies on k, i, j, x e y 
  
  # Control limits when the process is under control
  lsc <- 3*sqrt(var_mean_u) 
  lic <- -3*sqrt(var_mean_u) 
  
  # For every shift in the process
      for (i in 1:ls){ 
        set.seed(i*809) 
        dados <- rnorm((lk2*1000000), mean= delta[i]/sqrt(lk), sd= 1) # Generate data with shift
        m <- matrix(dados, ncol=lk2)                                  
        mrow <- t(apply(m,1,sort.int, method="quick"))                
        mrow <- mrow[,k]                                              
        mrow <- cbind(mrow, rowMeans(mrow))                           # It calculates the average for every NRSS sample
        ncolu <- ncol(mrow)
        mrow <- cbind(mrow, mrow[,ncolu]<lic | mrow[,ncolu]>lsc)      # Check if the average of every sample from NRSS is inside the control limits
        p <- sum(mrow[,ncolu+1])/dim(mrow)[1]
        cms_matrix[i,a] <- 1/p 
        print(cms_matrix)
        } 
}
return(cms_matrix)
}
cms_perf <- cms(aas, delta)

# Amostras a serem selecionadas de uma AAS para constituir uma URSS####
# a <- function(){for (k in 3:6)
# {
#   print(select(rnorm(k^2)))
# }}
# b <- a()
# 
# # Te?rico ####
# plot(log(pbinom(00:100,50,0.8)),type="l")


# Big Simulation #####  

# Obtaining the variance of the NRSS mean for imperfect case ####
aas <- c(9, 16, 25, 36)           # Sample Size
pho <- c(0,0.25,0.50,0.75,0.9,1)  # Correlation Coefficient
var <- 1                          # Variance of x and y
l <- 1000000                      # 1m of sample

require(MASS)
var_mean_imp <- function(aas, pho, l, var){
  # Argument:
    # aas: Sample size k^2
    # pho: Correlation coefficient between the variable of interest and the auxiliary variable
    # l  : Amount of samples
    # var: Variance of x and y
  # Description:
    # It calculates the variance of the mean from NRSS in a imperfect ranking scenario
  
  la <- length(aas)
  lp <- length(pho)
  var_mean_imp_values <- matrix(0,nrow=la,ncol=lp,dimnames=list(sqrt(aas),pho))
  seed <- 0
  
  for (i in 1:la){                   #For every sample size
    k2 <- aas[i]
    ktot <- k2*l
    for (p in 1:lp){                 #For every correlation
          seed <- seed + 1
          set.seed(seed)
          d1 <- mvrnorm(n = ktot,
                        mu = c(0,0),
                        Sigma = matrix(c(var,pho[p],pho[p],var),ncol=2,byrow=TRUE))
          dx <- matrix(d1[,1], ncol=k2)
          dy <- matrix(d1[,2], ncol=k2)
          rm(d1)
          pos <- t(apply(dx,1,order))             
          rm(dx)
          dys <- matrix(rep(0,ktot),ncol=k2)      
          for (linha in 1:nrow(dy)){              # Order y by x
            dys[linha,] <- dy[linha,pos[linha,]]
          }
          rm(dy);rm(pos)
          dys_amostras <- dys[,select(1:k2)]      
          rm(dys)
          mcov <- cov(dys_amostras)
          rm(dys_amostras)
          vari <- diag(mcov)
          mcov[lower.tri(mcov, diag=TRUE)] <- 0
          var_mean_imp_values[i,p] <- (1/k2)*sum(vari)+(2/k2)*sum(colSums(mcov))
          print(var_mean_imp_values)              
    }
  }
  return(var_mean_imp_values)
}
require(MASS)
memory.limit(size=10000) # Increase the memory size of computer to run this code (4GB of RAM was not sufficient)

# Teste de Qualidade: OK ####
c <- var_mean_imp(aas,pho,l,var)
str(c)
# var_mean_imp <- read.csv("C:/Users/Guilherme/Google Drive/TCC/Simulacao/var_mean_imp.txt", row.names=1, sep="",
#                          col.names=paste("pho ",c(0,0.25,0.50,0.75,0.90,1),sep=""))

var_urss2(select(rnorm(9)))
var_urss2(select(rnorm(16)))
var_urss2(select(rnorm(25)))
var_urss2(select(rnorm(36)))
var_value <- var_mean_imp[3,4]
var_da_media_matrix_fake <- matrix(var_value,nrow=4, ncol=6)

# ARL for imperfect ranking ####
cms_imp <- function(aas, pho, var, l, var_da_media_matrix, delta){
  # Argument:
    # aas: Sample size k^2
    # pho: Correlation coefficient between the variable of interest and the auxiliary variable
    # l  : Amount of samples
    # var: Variance of x and y
    # var_da_media_matrix: matrix from var_mean_imp(), which contains all values from variance of the mean (imperfect ranking)
  # Description:
    # It calculates the ARL under imperfect ranking scenario using 3sigma limits
  
  la <- length(aas)
  lp <- length(pho)
  ld <- length(delta)
  cms_imp <- list('n= 3' = 1) # It keeps CMS's values
  for (i in 1:la){                 # For every sample size
    k2 <- aas[i]
    k <- select(1:k2)
    lk <- length(k)
    ktot <- k2*l
    cms_matrix <- matrix(0, nrow = ld, ncol = lp, dimnames=list(delta,pho))
    for (p in 1:lp){                               # For every correlation level
      var_mean_u_imp <- var_da_media_matrix[i, p]
      lsc <-  3*sqrt(var_mean_u_imp)                # Controle limits rely on k and pho
      lic <- -3*sqrt(var_mean_u_imp)
      for (d in 1:ld){                             # For every shift in the process
        set.seed(9087832*p*d)
        d1 <- mvrnorm(n = ktot,
                      mu = c(0,delta[d]/sqrt(lk)),
                      Sigma = matrix(c(var,pho[p],pho[p],var),ncol=2,byrow=TRUE))
        dx <- matrix(d1[,1], ncol=k2)
        dy <- matrix(d1[,2], ncol=k2)
        rm(d1)
        pos <- t(apply(dx,1,order))                
        rm(dx)
        dys <- matrix(rep(0,ktot),ncol=k2)         
        for (linha in 1:nrow(dy)){                 # Order y by x
          dys[linha,] <- dy[linha,pos[linha,]]
        }
        rm(dy);rm(pos)
        dys_amostras <- dys[,k]                    
        rm(dys)
        dys_amostras_mean <- rowMeans(dys_amostras)
        rm(dys_amostras)
        mrow <- as.numeric(dys_amostras_mean<lic | dys_amostras_mean>lsc)
        s <- sum(mrow)/length(mrow)
        cms_matrix[d,p] <- 1/s
        gc()
        print(cms_matrix)
      }
    }
    cms_imp[[i]] <- cms_matrix
    print(cms_imp)
  }
  names(cms_imp) <- paste('n=',sqrt(aas),sep=" ")
  return(cms_imp)
}

aas <- c(9,16,25,36)              # Sample Size
pho <- c(0,0.25,0.50,0.75,0.90,1) # Correlation Coefficient
var <- 1                          # Variance of x and y
l <- 1000000                      # 1M of samples
delta <- c(0, 0.1, 0.2, 0.3, 0.4, 0.8, 1.2, 1.6, 2, 2.4, 3.2)
urss_pho_07_imp <- cms_imp(aas, pho, var, l, var_mean_imp, delta)

