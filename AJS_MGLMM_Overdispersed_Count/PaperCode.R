# Taken from Fisher dispersion index for multivariate count distributions: A review and a new proposal
# Supplementar Material
##############

# Data no.1 (Holgate, 1966)
# y=c(rep(0,34),rep(0,12),rep(0,4),rep(0,5),rep(0,2),rep(1,8),rep(1,13),
#     rep(1,3),rep(1,3),rep(2,3),rep(2,6),rep(2,1),rep(2,2),rep(3,1),rep(3,1),
#     rep(3,1),0)
# x=c(rep(0,34),rep(1,12),rep(2,4),rep(3,5),rep(4,2),rep(0,8),rep(1,13),
#     rep(2,3),rep(3,3),rep(0,3),rep(1,6),rep(2,1),rep(3,2),rep(0,1),rep(1,1),
#     rep(3,1),6)
# data=cbind(y,x)
# m=c(mean(data[,1]),mean(data[,2]))
# # GDI
# (t(m^.5)%*%var(data)%*%(m^.5))/(t(m)%*%m)
# # MDI
# (t(m^.5)%*%diag(diag(var(data)))%*%(m^.5))/(t(m)%*%m)
# # bootstrap gdi and mdi (\textbf{10000 replicates})
# n=length(y)
# nb=10000; z=seq(1,n);gdi=numeric(nb);mdi=numeric(nb)
# for(i in 1:nb){
#   zb=sample(z,n,replace=T)
#   datab=cbind(data[zb,1],data[zb,2])
#   m=c(mean(datab[,1]),mean(datab[,2]))
#   gdi[i]=(t(m^.5)%*%var(data)%*%(m^.5))/(t(m)%*%m)
#   mdi[i]=(t(m^.5)%*%diag(diag(var(datab)))%*%(m^.5))/(t(m)%*%m)
# }
# hist(gdi)
# hist(mdi)
# # bootstrap mean and sd
# c(mean(gdi),sd(gdi))
# c(mean(mdi),sd(mdi))
# # 95% CI approx.
# quantile(gdi,c(0.025,0.975))
# quantile(mdi,c(0.025,0.975))

#Function that returns the vector composed by the averages of each variable
aux_esp<-function(data){
  res<-c()
  n=dim(data)
  n=n[2]
  for (i in 1:n){
    res[i]=mean(data[,i])
  }
  return(res)
}
#Function that returns the vector composed by the variances of each variable
aux_var<-function(data){
  res<-c()
  n=dim(data)
  n=n[2]
  for (i in 1:n){
    res[i]=var(data[,i])
  }
  return(res)
}
# Function that returns the matrix composed by the variances on its diagonal and 0 else
aux_diag_var<-function(data){
  res<-c()
  n=dim(data)
  n=n[2]
  for (i in 1:n){
    res[i]=var(data[,i])
  }
  return(diag(res))
}
# Function that calculates the MDI
MDI<-function(data){
  tmp_E=aux_esp(data)
  tmp_D=aux_diag_var(data)
  res=sqrt(tmp_E)%*%tmp_D%*%sqrt((tmp_E))/(tmp_E%*%(tmp_E))
  return(res)
}
# Function that calculates the GDI
GDI<-function(data){
  tmp_E=aux_esp(data)
  tmp_C=cov(data)
  res=sqrt(tmp_E)%*%tmp_C%*%sqrt((tmp_E))/(tmp_E%*%(tmp_E))
  return(res)
}
# Fonction that returns the delta vector
vec_delta<-function(data){
  k=dim(data)
  k=k[2]
  res<-c()
  tmp_E=aux_esp(data)
  for (i in 1:k){
    res[i]=0
    for (j in 1:k){
      res[i]=res[i]+tmp_E[j]/tmp_E[i]*cov(data[,i],data[,j])
    }
    res[i]=(res[i]-2*tmp_E[i]*GDI(data))/(tmp_E%*%(tmp_E))
  }
  for (i in 1:k){
    for(j in i:k){
      if (i==j){
        tmp=tmp_E[i]/(tmp_E%*%(tmp_E))
      }
      else{
        tmp=2*sqrt(tmp_E[i]*tmp_E[j])/(tmp_E%*%(tmp_E))
      }
      res<-c(res,tmp)
    }
  }
  return(res)
}
# Function that returns the gamma3 matrix
mat_gamma3<-function(data){
  k=dim(data)
  k=k[2]
  tmp_E=aux_esp(data)
  res<-c()
  for (j in 1:k){
    tmp_vec<-c()
    for (jj in 1:k){
      for (ll in jj:k){
        tmp=cov(data[,j],data[,jj]*data[,ll])
        tmp_vec=c(tmp_vec,tmp)
      }
    }
    res=rbind(res,tmp_vec)
  }
  rownames(res)=c()
  return(res)
}
# Function that returns the gamma4 matrix
mat_gamma4<-function(data){
  k=dim(data)
  k=k[2]
  tmp_E=aux_esp(data)
  res<-c()
  for (j in 1:k){
    for (l in j:k){
      tmp_vec<-c()
      for (jj in 1:k){
        for (ll in jj:k){
          tmp=cov(data[,j]*data[,l],data[,jj]*data[,ll])
          tmp_vec=c(tmp_vec,tmp)
        }
      }
      res=rbind(res,tmp_vec)
    }
  }
  rownames(res)=c()
  return(res)
}
# Function that returns the gamma matrix
mat_gamma<-function(data){
  tmp=cbind(cov(data),mat_gamma3(data))
  tmp2=cbind(t(mat_gamma3(data)),mat_gamma4(data))
  res=rbind(tmp,tmp2)
  return(res)
}
# Function that calculates the variance of GDI
sigma_GDI<-function(data){
  delta=vec_delta(data)
  gamma=mat_gamma(data)
  res=delta%*%gamma%*%delta
  return(res)
}
vec_lambda<-function(data){ # Function that returns the lambda
  k=dim(data)
  k=k[2]
  res<-c()
  tmp_E=aux_esp(data)
  tmp_V=aux_var(data)
  for (i in 1:k){
    res[i]=(tmp_V[i]-2*tmp_E[i]*MDI(data))/(tmp_E%*%tmp_E)
  }
  for (i in 1:k){
    res[i+k]=tmp_E[i]/(tmp_E%*%tmp_E)
  }
  return(res)
}
mat_pi3<-function(data){ # Function that returns the pi3 matrix
  k=dim(data)
  k=k[2]
  tmp_E=aux_esp(data)
  res<-c()
  for (j in 1:k){
    tmp_vec<-c()
    for (jj in 1:k){
      tmp=cov(data[,j],data[,jj]*data[,jj])
      tmp_vec=c(tmp_vec,tmp)
    }
    res=rbind(res,tmp_vec)
  }
  rownames(res)=c()
  return(res)
}
mat_pi4<-function(data){ # Function that returns the pi4 matrix
  k=dim(data)
  k=k[2]
  tmp_E=aux_esp(data)
  res<-c()
  for (j in 1:k){
    tmp_vec<-c()
    for (jj in 1:k){
      tmp=cov(data[,j]*data[,j],data[,jj]*data[,jj])
      tmp_vec=c(tmp_vec,tmp)
    }
    res=rbind(res,tmp_vec)
  }
  rownames(res)=c()
  return(res)
}
mat_pi<-function(data){ # Function that returns the pi matrix
  tmp=cbind(cov(data),mat_pi3(data))
  tmp2=cbind(t(mat_pi3(data)),mat_pi4(data))
  res=rbind(tmp,tmp2)
  return(res)
}
sigma_MDI<-function(data){ # Function that calculates the variance of MDI
  lambda=vec_lambda(data)
  pi=mat_pi(data)
  res=lambda%*%pi%*%lambda
  return(res)
}
ind_tab<-function(data,n){ # Function that returns a table with all the information
  k=length(n)
  res<-c()
  for (i in 1:k){
    tmp_vec<-c()
    tdata=data[1:n[i],]
    tmp_vec=c(det(cor(tdata)),sigma_GDI(tdata),sigma_MDI(tdata),GDI(tdata)
              ,qnorm(0.975)*sqrt(sigma_GDI(tdata)/n[i]),MDI(tdata),
              qnorm(0.975)*sqrt(sigma_MDI(tdata)/n[i]))
    res=rbind(res,tmp_vec)
  }
  rownames(res)=n
  colnames(res)=c("det(rho)","sigma2GDI","sigma2MDI","GDI +-","E.t(int)",
                  "MDI+-","E.t(int)")
  return(res)
}
# # Example 1
# MDI(data)
# GDI(data)
# # Standard error
# sqrt(sigma_GDI(data))/10
# sqrt(sigma_GDI(data)/nrow(data))
# sqrt(sigma_MDI(data))/10
