# Author: Guilherme Parreira da Silva
# fixed effects are in the linear predictor function, standard deviation and correlation
# 25/09/2020
#########################################################################
# Start values and summary from fit/optimizer (without need for sdreport)
#########################################################################
cria_cor_matrix_tmb_sd_from_vector <- function(nresp, cor){
  dpresent <- diag(nresp)
  counter <- 1
  namescor <- NULL
  for(i in 1:nresp){
    for(j in 1:nresp){
      if (i>j){
        dpresent[i,j] <- round(cor[counter],4)
        namescor <- c(namescor, paste0("r_", "Y", j, "_Y", i))
        counter <- counter + 1
      }
    }
  }
  names(cor) <- namescor
  return(list(corre = cor, matrix_corre = dpresent))
}
cor_from_matrix <- function(corre){
  pos <- which(lower.tri(corre, diag=FALSE), arr.ind = T)
  pos <- pos[order(pos[,1], pos[,2]), , drop = F]
  pars_est_rho <- corre[upper.tri(corre, diag=FALSE)]
  names(pars_est_rho) <- paste0("r_", "Y" ,pos[,2], "_", "Y", pos[, 1])
  return(pars_est_rho)
}


sv <- function(fit, time ,n_betas, nr, model, method = "nlminb", threads, units=60*60){
  pars_estimados <- fit$par
  pars_est_beta <- pars_estimados[grep("beta", names(pars_estimados))]
  resp <- gsub("beta", "", names(pars_est_beta))
  # tresp <- table(resp)
  # matt <- matrix(0, nrow = max(tresp), ncol = nr)
  # for (i in 1:length(tresp)){
  #   matt[, i] <- pars_est_beta[grep(i, names(pars_est_beta))]
  # }
  pars_est_beta <- matrix(pars_est_beta, nrow = n_betas, ncol = nr)
  colnames(pars_est_beta) <- paste0("beta",1:nr)
  pars_est_theta <- pars_estimados[names(pars_estimados)=="rho"]
  
  ############## Variance parameters adjustment #######################################
  if(!grepl("var_fixed", model)){
    pars_est_sigma <- pars_estimados[names(pars_estimados)=="sigma"]
    names(pars_est_sigma) <- paste0("sig", 1:length(pars_est_sigma))
    if (grepl("var_comum", model)){
      corre <- getSigma(pars_est_theta, rep(1, nr), varcov = F) # sigma is fake, only need the dimension
    } else if(!grepl("rho_fixed", model)){ ##Added after using rho fixed
      corre <- getSigma(pars_est_theta, pars_est_sigma, varcov = F)
    }
  } else{
    pars_est_sigma <- NULL
    corre <- getSigma(pars_est_theta, rep(1, nr), varcov = F) # sigma is fake, only need the dimension
  }
  
  ############## Correlation parameters adjustment #######################################
  if(!grepl('rho_fixed', model)){
    pars_est_rho <- cor_from_matrix(corre)
  } else{
    pars_est_rho <- rep(0, (nr*(nr-1))/2)
  }
  
  ############## Dispersion parameters adjustment #####################################
  if(grepl("comp", model) && !grepl("nu_fixed", model)){
    pars_est_disp <- pars_estimados[grep("nu", names(pars_estimados))]
    names(pars_est_disp) <- paste0("log(nu)", 1:nr)
  } else if(grepl("negbin", model) && !grepl("phi_fixed", model)){
    pars_est_disp <- pars_estimados[grep("phi", names(pars_estimados))]
    names(pars_est_disp) <- paste0("log(phi)", 1:nr)
  } else{
    pars_est_disp <- NULL}
  summary <- data.frame(Loglik = as.numeric(fit[2]),
                        Convergence = fit$convergence,
                        Threads = threads, 
                        Time=time[1]/units,
                        Optimizer = method)
  return(list(beta = pars_est_beta,
              sigma = pars_est_sigma,
              rho = pars_est_rho,
              disp = pars_est_disp,
              summary = summary))
}

# Utilizado para os modelos finais (gigantes da ANT), para não ter que digitar um a um os parametros
sv2 <- function(fit, time , nr, model, method = "nlminb", threads, rows, units=60*60){
  pars_estimados <- fit$par
  pars_est_beta <- pars_estimados[grep("beta", names(pars_estimados))]
  # i <- 1
  list <- c()
  for (i in 1:nr){
    list[[i]] <- as.numeric(pars_est_beta[grep(paste0("^beta", i,"$"), names(pars_est_beta))])
  }
  names(list) <- paste0("beta",1:nr)
  pars_est_theta <- pars_estimados[names(pars_estimados)=="rho"]
  
  ############## Variance parameters adjustment #######################################
  if(!grepl("var_fixed", model)){
    pars_est_sigma <- pars_estimados[names(pars_estimados)=="sigma"]
    names(pars_est_sigma) <- paste0("sig", 1:length(pars_est_sigma))
    if (grepl("var_comum", model)){
      corre <- getSigma(pars_est_theta, rep(1, nr), varcov = F) # sigma is fake, only need the dimension
    } else{
      corre <- getSigma(pars_est_theta, pars_est_sigma, varcov = F)
    }
  } else{
    pars_est_sigma <- NULL
    corre <- getSigma(pars_est_theta, rep(1, nr), varcov = F) # sigma is fake, only need the dimension
  }
  pars_est_rho <- cor_from_matrix(corre)
  
  ############## Dispersion parameters adjustment #####################################
  if(grepl("comp", model) && !grepl("nu_fixed", model)){
    pars_est_disp <- pars_estimados[grep("nu", names(pars_estimados))]
    names(pars_est_disp) <- paste0("log(nu)", 1:nr)
  } else if(grepl("negbin", model) && !grepl("phi_fixed", model)){
    pars_est_disp <- pars_estimados[grep("phi", names(pars_estimados))]
    names(pars_est_disp) <- paste0("log(phi)", 1:nr)
  } else{
    pars_est_disp <- NULL}
  summary <- data.frame(Loglik = as.numeric(fit[2]),
                        Convergence = fit$convergence,
                        Threads = threads, 
                        Time=time[1]/units,
                        Optimizer = method)
  list[[i+1]] <- matrix(0, nrow = rows, ncol = nr)
  list[[i+2]] <- pars_est_rho
  list[[i+3]] <- pars_est_sigma
  list[[i+4]] <- pars_est_disp
  list[[i+5]] <- summary
  names(list)[(nr+1):(nr+5)] <- c("U", "rho", "sigma", "disp", "summary")
  return(list)
}
#########################################################################
# Start values and summary from sdreport
#########################################################################
# fixef <- rep_fe
# cor <- rep_cor
# nresp <- n_resp
pretty_matrix <- function(est, std, nresp){
  dpresent <- as.data.frame(diag(nresp))
  for(i in 1:nresp){
    for(j in 1:nresp){
      if (i>j){
        dpresent[i,j] <- paste0(round(est[i,j],4), " (",
                                round(std[i,j],4), ")")
      }
    }
  }
  return(dpresent)
}
cria_cor_matrix_tmb_sd <- function(nresp, cor){
  # cor <- unique(cor)
  # cor <- cor[!cor[, 1]==1, ]
  m_est <- matrix(cor[,1], ncol = nresp, byrow = F)
  m_std <- matrix(cor[,2], ncol = nresp, byrow = F)
  vec_est <- cor_from_matrix(m_est)
  vec_std <- cor_from_matrix(m_std)
  cor <- cbind(Estimate = vec_est, `Std. Error` = vec_std)
  dpresent <- pretty_matrix(m_est, m_std, nresp)
  return(list(corre = cor, matrix_corre = dpresent))
}
# theta <- pars_est_theta
# sigma <- pars_est_sigma
getSigma <- function(theta, sigma, byrow = T, varcov = T){
  nr <- length(sigma)
  L <- diag(1, nrow = nr, ncol = nr)
  L <- cria_cor_matrix_tmb_sd_from_vector(nr, theta)$matrix_corre
  recheio <- L %*% t(L)    
  D <- diag(recheio)
  Dinvsqrt <- diag(1/sqrt(D))     
  if (varcov){ # VarCov
    W <- diag(sigma)     
    out <- W %*% Dinvsqrt %*% recheio %*% Dinvsqrt %*% t(W) #retorna a matriz de variância-covariancia
  } else{      # VarCorr
    out <- Dinvsqrt %*% recheio %*% Dinvsqrt
    # diag(out) <- sigma^2
  }
  return(out) 
}
TMB_summary <- function(fixef, adreport, nresp=5, namebeta = NULL, lengthresp = NULL,
                        model = ""){
  beta <- fixef[grep("beta", rownames(fixef), invert = F),]
  ncovs <- nrow(beta)/nresp
  # ncovs <- 4
  if (is.null(namebeta)){
    rownames(beta) <- paste0(paste0("b", 0:(ncovs-1), "_Y"), rep(1:nresp, each = ncovs))
  } else{
    rownames(beta) <- paste0(namebeta, "_Y", lengthresp)
  }
  sigma <- fixef[grep("sigma", rownames(fixef), invert = F),, drop = F]
  if (nrow(sigma)>0){ #It corrects for var fixed
    rownames(sigma) <- paste0("sig", 1:nrow(sigma))
  }
  report_corre <- rownames(adreport)=="Cor"
  if (nresp>1){ # I still need to check whether it works for 2 responses
    cor <- adreport[report_corre, ]
    corre <- cria_cor_matrix_tmb_sd(nresp = nresp, cor = cor)
    # Dispersion parameter (2nd parameter of a distribution -> PHI_NB; NU_COMP)
    report_extra <- adreport[!report_corre, ]
    #Necessary to accommodate the  model with rho fixed = 0
    if (!grepl('rho_fixed', model)){
      out <- rbind(beta, sigma, cor = corre$corre, report_extra)
    } else{
      out <- rbind(beta, sigma, report_extra)
    }
    out <- cbind(out, out[,1]/out[,2])
    out <- cbind(out, round(pnorm(abs(out[,3]), lower.tail = F)*2, 5))
    colnames(out)[3:4] <- c("Z", "p value")
  } else{
    report_extra <- adreport[!report_corre, ]
    out <- rbind(beta, sigma, report_extra)
    out <- cbind(out, out[,1]/out[,2])
    out <- cbind(out, round(pnorm(abs(out[,3]), lower.tail = F)*2, 5))
  }
  out_final <- list(out, corre$matrix_corre)
  names(out_final) <- c("ests", "cor_mat") #List names
  return(out_final)
}
sv_sd <- function(obj, fit, time, units=60*60, nr, model=NULL, method, threads, namebeta = NULL, lengthresp = NULL){
  rep <- sdreport(obj)
  rep_fe <- summary(rep, "fixed")
  rep_cor <- summary(rep, "report")
  all_estimates <- TMB_summary(fixef = rep_fe, adreport = rep_cor, nresp = nr, namebeta = namebeta, lengthresp = lengthresp, model = model)
  # Estimativa dos parâmetros
  estimates <- all_estimates$ests
  # Start_beta está errado quando as covs não forem iguais para todas as resps
  start_beta <- matrix(estimates[grep("b\\d+_Y", rownames(estimates)), 1],ncol = nr)
  colnames(start_beta) <- paste0("beta", 1:nr)
  if(grepl("var_fixed", model)){
    start_sigma <- NULL
  } else{
    start_sigma <- estimates[grep("sig", rownames(estimates)), 1]
  }
  
  if(grepl('rho_fixed', model)){
    start_rho <- NULL
    # start_rho <- rep(0, (nr*(nr-1))/2)
  } else{
    start_rho <- estimates[grep("r_Y\\d+_Y\\d+", rownames(estimates)), 1]
  }
  
  if(grepl("comp", model) && !grepl("nu_fixed", model)){
    start_disp <- log(estimates[grep("nu", rownames(estimates)), 1])
    names(start_disp) <- paste0("log(nu)", 1:nr)
  } else if(grepl("negbin", model) && !grepl("phi_fixed", model)){
    start_disp <- log(estimates[grep("phi", rownames(estimates)), 1])
    names(start_disp) <- paste0("log(phi)", 1:nr)
  } else{
    start_disp <- NULL
  }
  # Quadro resumo
  summary <- data.frame(Loglik=as.numeric(fit[2]), 
                        Convergence=fit$convergence, 
                        Threads = threads, 
                        Time=time[1]/units,
                        Optimizer = method)
  return(list(beta = start_beta,
              sigma = start_sigma,
              rho = start_rho,
              disp = start_disp,
              ests = estimates,
              summary = summary,
              corr_matrix = all_estimates$cor_mat))
}
#########################################################################
# Start values from MCGLM
#########################################################################
betas_mcglm <- function(fit){
  df <- coef(fit, type = "beta")
  nbetas <- nrow(df[df$Response==1,])
  df$betas <- rep(1:nbetas, times = max(df$Response))
  dfnew <- reshape(df[, -c(2,3)], 
                   idvar = "betas",
                   timevar = "Response",
                   direction = "wide")[, -1]
  names(dfnew) <- paste0("beta", 1:ncol(dfnew))
  dfnew <- as.matrix(dfnew)
  
  res_cov <- as.matrix(residuals(fit))
  sigma <- sqrt(diag(cov(res_cov)))
  return(list(beta = dfnew,
              sigma = sigma))
}
##############
# Send mail --------------------------------------------
##############
send_mail <- function(file, code, 
                      sender = "guilhermeparreira.silva@gmail.com", 
                      destinatary = "guilhermeparreira.silva@gmail.com"){
  gm_auth_configure(path = file)
  my_email <- gm_mime() %>% 
    gm_from(sender) %>% 
    gm_to(destinatary) %>% 
    gm_subject(paste0(code, " finished"))
    # gm_text_body("Code finished")
  gm_send_message(my_email)
}
##############
# Merges different packages results --------------------------------------------
##############
brms_summary <- function(summary){
  # It already prints cor and standard deviation
  tab <- rbind(summary$fixed,
               summary$spec_pars,
               summary$random$id)
  return(tab[,1:2])
}
# FE <- mcmc_poisson$Sol
# VCV <- mcmc_poisson$VCV
MCMCglmm_summary <- function(MCMC){
  # Fixed Effect -------------------------------------------------------------------
  FE <- MCMC$Sol
  test <- colnames(FE)
  nvar <- length(unique(sub(":.*", "", sub(".*trait", "", test))))
  cols_intercepto <- grep(":.*", test, invert = T)
  test[cols_intercepto] <- gsub("trait", "Intercept_", test[cols_intercepto])
  test <- gsub("trait", "", test)
  test[cols_intercepto] <- paste0(sub(".*_", "", test[cols_intercepto]), "_Intercept")
  test <- gsub(":", "_", test)
  init <- matrix(c(colMeans(FE),
                   apply(FE, 2, sd)), ncol = 2,
                 dimnames = list(c(test), c("Estimate", "Std.Error")))
  
  # VCV ---------------------------------------------------------------------------
  VCV <- MCMC$VCV
  VCV <- VCV[, !as.vector(upper.tri(diag(nvar), diag=FALSE))] # Get rid off repeated correlations
  colnames(VCV) <- gsub(".units", "", gsub("trait", "", colnames(VCV))) #Get rid off units and trait
  var_names <- gsub(":.*","", colnames(VCV))==gsub(".*:","", colnames(VCV)) #Variance variable names
  colnames(VCV)[var_names] <- gsub(":.*","", colnames(VCV)[var_names])
  VCV[, var_names] <- apply(VCV[, var_names], 2, sqrt) #Standard deviation
  cor_names <- !var_names                            # Covariance variables
  fors <- (1:length(cor_names))[cor_names] # Indices of the covariance
  for (i in fors){
    names <- colnames(VCV)[i]
    v1 <- gsub(":.*", "", names)
    v2 <- gsub(".*:", "", names)
    VCV[, i] <- VCV[, i]/((VCV[, v1])*(VCV[, v2])) #cov/(sd(1)*sd(2))
  }
  # Padroniza o nome das variáveis (sd, cor)
  colnames(VCV)[var_names] <- paste0("sd(", colnames(VCV)[var_names], ")")
  colnames(VCV)[cor_names] <- paste0("cor(", colnames(VCV)[cor_names], ")")
  # Padroniza a ordem das variáveis (sd, cor)
  ordem_cols <- c(grep("sd\\(", colnames(VCV)),grep("cor\\(", colnames(VCV)))
  VCV <- VCV[, ordem_cols]
  init2 <- matrix(c(colMeans(VCV),
                    apply(VCV, 2, sd)), ncol = 2,
                  dimnames = list(colnames(VCV), c("Estimate", "Std.Error")))
  # Shape ---------------------------------------------------------------------------
  
  # Merges everything ---------------------------------------------------------------------------
  init <- rbind(init, init2)
  return(init)
}
all_summaries <- function(brms_summ, MCMC, TMB_fixef, TMB_cor, nresp){
  brms <- brms_summary(brms_summ)
  MCMC <- MCMCglmm_summary(MCMC)
  TMB <- TMB_summary(TMB_fixef, TMB_cor, nresp)
  rownames(TMB) <- rownames(brms)
  out <- cbind(brms[, 1], MCMC[, 1], TMB[, 1],
  brms[, 2], MCMC[, 2], TMB[, 2])
  out <- round(out, 4)
  colnames(out) <- c(paste0(c("E_"), c("BRMS", "MCMC", "TMB")),
                     paste0(c("SE_"), c("BRMS", "MCMC", "TMB")))
  return(out)
}
##############
# Auxiliar and rarely used function --------------------------------------------
##############
z_to_r <- function(z){
  (exp(2*z)-1)/(exp(2*z)+1)
}
r_to_z <- function(r){
  (1/2)*log((1 + r)/(1-r))
}
TMB_summary_wag <- function(fixef, cor, nresp=2){
  # Beta (Fixed Effects)
  beta <- fixef[grep("beta", rownames(fixef), invert = F),]
  
  # Sigma (Standard error of random effect)
  sigma <- cor[grep("sigma", rownames(cor), invert = F),]
  rownames(sigma) <- c("s2_1", "s2_2")
  # Phi Dispersion parameter (2nd parameter of a distribution)
  phi <- cor[grep("phi", rownames(cor), invert = F),]
  rownames(phi) <- c("phi1", "phi2")
  # Correlation of random effects (Now)
  corre <- cor[grep("rho", rownames(cor), invert = F),]
  out <- rbind(beta, sigma, phi = phi, cor = corre)
  return(out)
}

##############
# Master Code ------------------------------------------------------------------------------
##############
TMB_models <- function(data, parameters, modelTMB, n_params, n_resp, 
                       n, n300, modelname, threads, full = T){
  openmp(threads)
  obj <- MakeADFun(data = data,
                   parameters = parameters,
                   DLL = modelTMB,
                   hessian = T,
                   silent = F
                   ,random = c("U"))
  a <- system.time(fit_NLB <<- print(nlminb(obj$par, obj$fn, obj$gr,
                                            control = list(eval.max = 5000, iter.max = 5000,
                                                           abs.tol = 1e-04, rel.tol = 1e-04))))
  
  # From sd report
  rep <- sdreport(obj)
  rep_fe <- summary(rep, "fixed")
  rep_cor <- summary(rep, "report")
  estimates <- TMB_summary(rep_fe, rep_cor, n_resp)
  
  model_summary <- data.frame(model = modelname, 
                              horas = a[1]/(60*60),
                              convergence = fit_NLB$convergence,
                              loglik = fit_NLB$objective)
  
  # Valores iniciais para Amostra Completa
  if (modelname == "poisson"){
    params_start <- list(
      beta = as.matrix(
        data.frame("beta1" = estimates[grep("b[0-9]+_Y1", row.names(estimates)), 1],
                   "beta2" = estimates[grep("b[0-9]+_Y2", row.names(estimates)), 1],
                   "beta3" = estimates[grep("b[0-9]+_Y3", row.names(estimates)), 1],
                   "beta4" = estimates[grep("b[0-9]+_Y4", row.names(estimates)), 1],
                   "beta5" = estimates[grep("b[0-9]+_Y5", row.names(estimates)), 1])),
      U = matrix(0, ncol = n_resp, nrow = ifelse(full, n, n300)),
      rho = estimates[grep("r_Y", row.names(estimates)), 1],
      sigma = estimates[grep("sig", row.names(estimates)), 1])
  } else if (modelname == "negbin"){
    # Phi
    phi <- fit_NLB$par[names(fit_NLB$par)=="phi"]
    params_start <- list("beta1" = estimates[grep("b[0-9]+_Y1", row.names(estimates)), 1],
                         "beta2" = estimates[grep("b[0-9]+_Y2", row.names(estimates)), 1],
                         "beta3" = estimates[grep("b[0-9]+_Y3", row.names(estimates)), 1],
                         "beta4" = estimates[grep("b[0-9]+_Y4", row.names(estimates)), 1],
                         "beta5" = estimates[grep("b[0-9]+_Y5", row.names(estimates)), 1],
                         U = matrix(0, ncol = n_resp, nrow = ifelse(full, n, n300)),
                         rho = estimates[grep("r_Y", row.names(estimates)), 1],
                         sigma = estimates[grep("sig", row.names(estimates)), 1],
                         phi = phi)
    
  } else if (modelname == "compoisson"){
    nu <- fit_NLB$par[names(fit_NLB$par)=="nu"]
    params_start <- list("beta1" = estimates[grep("b[0-9]+_Y1", row.names(estimates)), 1],
                         "beta2" = estimates[grep("b[0-9]+_Y2", row.names(estimates)), 1],
                         "beta3" = estimates[grep("b[0-9]+_Y3", row.names(estimates)), 1],
                         "beta4" = estimates[grep("b[0-9]+_Y4", row.names(estimates)), 1],
                         "beta5" = estimates[grep("b[0-9]+_Y5", row.names(estimates)), 1],
                         U = matrix(0, ncol = n_resp, nrow = ifelse(full, n, n300)),
                         rho = estimates[grep("r_Y", row.names(estimates)), 1],
                         sigma = estimates[grep("sig", row.names(estimates)), 1],
                         nu = nu)
  }
  
  return(list(params_start,
              model_summary,
              estimates))
}