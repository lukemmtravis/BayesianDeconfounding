library(glmnet)
library(sparsevb)
library(tidyverse)
library(ggplot2)
library(kableExtra)
library(pbapply)
library(mcreplicate)
library(DDL)

setwd('/Users/lmt15/Documents/phd/Causal Inference/Codes')
Rcpp::sourceCpp("SASGibbsSampler.cpp")
# make_design_old = function(n, p, Sigma=NA){
#   if(any(is.na(Sigma))){
#     X = matrix(rnorm(n*p), nrow = n, ncol = p)
#   }else{
#     U = chol(Sigma)
#     normal_sample = matrix(rnorm(n*p), nrow = n, ncol = p)
#     X = normal_sample%*%U  
#   }
#   return(X)
# }

#### Data Generation Functions ####
make_design = function(n, p, rho=0.0, sigma_E=1, AR=FALSE){
  if(!AR){
    Sigma = matrix(rho, nrow = p, ncol = p)  
    diag(Sigma) = rep(sigma_E**2,p)
  }else{
    Sigma = matrix(nrow=p, ncol=p)
    for(i in c(1:p)){
      for(j in c(1:p)){
        if(i == j){
          Sigma[i, j] = sigma_E**2
        }else{
          Sigma[i, j] = rho**abs(i - j)  
        }
      }
    }
  }
  U = chol(Sigma)
  normal_sample = matrix(rnorm(n*p), nrow = n, ncol = p)
  X = normal_sample%*%U
  return(X)
}
make_true_param = function(p, s0, signal_size = NA, include_first_param=TRUE){
  beta = rep(0, p)
  if(include_first_param){
    S = sample(c(2:p), size = (s0 - 1), replace=FALSE)
    if(is.na(signal_size)){
      beta[1] = rnorm(1)
      beta[S] = rnorm(s0 - 1)
    }else{
      beta[1] = signal_size
      beta[S] = signal_size
    }
  }else{
    S = sample(c(1:p), size = s0, replace=FALSE)
    if(is.na(signal_size)){
      beta[S] = rnorm(s0)
    }else{
      beta[S] = signal_size
    } 
  }
  return(beta)
}
make_data = function(n, p, s0, q, signal_size=NA, rho=0.0, noise_var=1, AR=FALSE, dataset='generated', index=1, confounding_layers=1){
  if(dataset == 'riboflavin'){
    X = as.matrix(read.csv('data/riboflavin_normalized.csv'))
    # Override n and p to be correct for the riboflavin data
    n = dim(X)[1]
    p = dim(X)[2]
    beta_0 = make_true_param(p, s0, signal_size=signal_size, include_first_param=TRUE)
    # reorder columns so that index is first.
    X = X[,c(index, c(1:p)[-index])]
    e = rnorm(n, sd=sqrt(noise_var))
    if(q == 0){
      Y = X %*% beta_0 + e
    }else{
      Gamma = matrix(rnorm(q*p), nrow = q, ncol = p)
      # H_hat below is the L2-minimiser of X - H %*% Gamma
      H_hat = X %*% t(Gamma) %*% solve(Gamma %*% t(Gamma))
      delta = rep(signal_size, q)  
      Y = X %*% beta_0 + H_hat %*% delta + e
    } 

  }else if(dataset == 'generated'){
    beta_0 = make_true_param(p, s0, signal_size=signal_size, include_first_param=TRUE)
    E = make_design(n, p, rho=rho, AR=AR)
    e = rnorm(n, sd=sqrt(noise_var))
    
    # if(q > 0){
    #   Gamma = matrix(rnorm(q*p), nrow = q, ncol = p)
    #   H = matrix(rnorm(n*q), nrow = n, ncol = q)
    #   delta = rnorm(q)  
    #   X = H %*% Gamma + E
    #   Y = X %*% beta_0 + H %*% delta + e
    # }else{
    #   X = E
    #   Y = X %*% beta_0 + e
    # }  

    if(q > 0){
      # H is the lowest realisation of confounders
      H = matrix(rnorm(n*q), nrow = n, ncol = q)
      # Gamma transforms the highest hidden set of confounders to X
      Gamma = matrix(rnorm(q*p), nrow = q, ncol = p)
      if(confounding_layers > 1){
        H_hidden = H
        for(i in c(1:(confounding_layers-1))){
          Gamma_hidden = matrix(rnorm(q*q), nrow = q, ncol = q)
          H_hidden = H_hidden %*% Gamma_hidden + matrix(rnorm(n*q), nrow=n, ncol=q)
        }  
      }else{
        H_hidden=H
      }
      
      delta = rnorm(q)  
      X = H_hidden %*% Gamma + E
      
      Y = X %*% beta_0 + H %*% delta + e
    }else{
      X = E
      Y = X %*% beta_0 + e
    }
  }else{
    stop('dataset must be in {generated, riboflavin}')
  }
  

  return(list(X = X, Y = Y, beta_0 = beta_0))
}
make_data_CBM = function(n, p, s0, q, rho=0.0, sigma_E=2){
  # print('running cbm sampling scheme')
  beta_0 = c(rep(1, s0), rep(0, p - s0))
  E = make_design(n, p, rho=rho, sigma_E=sigma_E)
  e = rnorm(n)
  if(q > 0){
    Gamma = matrix(rnorm(q*p), nrow = q, ncol = p)
    H = matrix(rnorm(n*q), nrow = n, ncol = q)
    delta = rnorm(q)  
    X = H %*% Gamma + E
    Y = X %*% beta_0 + H %*% delta + e
  }else{
    X = E
    Y = X %*% beta_0 + e
  }
  return(list(X = X, Y = Y, beta_0 = beta_0))
}

####

#### Transforms ####
trim_transform = function(X, tau_quantile=0.5){
  svd_X = svd(X)
  d = svd_X$d
  # if(is.na(tau_quantile)){
  #   tau = median(d)
  # }else{
  #   tau = as.numeric(quantile(d, tau_quantile)
  # }
  tau = as.numeric(quantile(d, tau_quantile))
  d_tilde = pmin(d, tau)
  F_mat = svd_X$u %*% diag(d_tilde/d) %*% t(svd_X$u)
  return(F_mat)
}
puffer_transform = function(X, tau_quantile=0.5){
  svd_X = svd(X)
  d = svd_X$d
  tau = as.numeric(quantile(d, tau_quantile))
  d_tilde = rep(tau, length(d))
  F_mat = svd_X$u %*% diag(d_tilde/d) %*% t(svd_X$u)
  return(F_mat)
}
identity_transform = function(X){
  n = dim(X)[1]
  return(diag(n))
}
lava_transform = function(X, lambda=1, cbm_lambda=FALSE, cbm_quantile=0.5){
  n = dim(X)[1]
  p = dim(X)[2]
  if(cbm_lambda){
    ds = svd(X)$d
    lambda = (1/n)*as.numeric(quantile(ds, cbm_quantile))**2
  }else{
    lambda=lambda
  }
  A = diag(n) - X %*% solve(t(X)%*%X + n*lambda*diag(p)) %*% t(X)
  svd_A = svd(A)
  L = svd_A$u %*% diag(sqrt(svd_A$d)) %*% t(svd_A$u)
  return(L)
}
oracle_removal_transform = function(X, q=5){
  n = dim(X)[1]
  p = dim(X)[2]
  svd_X = svd(X)
  L = svd_X$u %*% diag(c(rep(0, q), rep(1, n-q))) %*% t(svd_X$u)
  return(L)
}
make_L = function(X, Sigma){
  n = dim(X)[1]
  p = dim(X)[2]
  A = diag(n) - X %*% solve(t(X)%*%X + solve(Sigma)) %*% t(X)
  svd_A = svd(A)
  L = svd_A$u %*% diag(sqrt(svd_A$d)) %*% t(svd_A$u)
  return(L)
}

sigma_est = function(X, Y, L){
  X_tilde = L %*% X
  Y_tilde = L %*% Y
  
  cv_fit = cv.glmnet(X_tilde, Y_tilde)
  lambda_min = cv_fit$lambda.min
  best_fit = glmnet(X_tilde, Y_tilde, lambda = lambda_min)
  beta_hat = best_fit$beta 

  RSS = sum( (Y_tilde - X_tilde %*% beta_hat)**2 )
  sigma_est_sq = RSS/sum(svd(L)$d)

  return(sqrt(sigma_est_sq))
}
####

#### Metrics for Analysis ####
l2_dist = function(v1, v2){
  return(sqrt(sum((v1 - v2)**2)))
}
l1_dist = function(v1, v2){
  return(sum(abs(v1 - v2)))
}
precision = function(beta_hat, beta_0){
  pred_active = which(beta_hat != 0)
  if(length(pred_active) == 0){
    return(1)
  }
  return(mean(beta_0[pred_active] != 0))
}
recall = function(beta_hat, beta_0){
  true_active = which(beta_0 != 0)
  return(mean(beta_hat[true_active] != 0))
}
f1 = function(beta_hat, beta_0){
  prec = precision(beta_hat, beta_0)
  rec = recall(beta_hat, beta_0)
  if(prec == 0 & rec == 0){
    return(0)
  }else{
    return(prec*rec/(prec + rec))  
  }
}
active_hit = function(fit, beta_0){
  active_inds = which(beta_0 != 0)
  CIs = fit$CIs[active_inds,]
  ind_vec = (CIs[,1] <= beta_0[active_inds] & beta_0[active_inds] <= CIs[,2])*1
  hit_prop = mean(ind_vec)
  return(hit_prop)
}
inactive_hit = function(fit, beta_0){
  active_inds = which(beta_0 != 0)
  CIs = fit$CIs[-active_inds,]
  ind_vec = (CIs[,1] <= beta_0[-active_inds] & beta_0[-active_inds] <= CIs[,2])*1
  hit_prop = mean(ind_vec)
  return(hit_prop)
}
active_length = function(fit, beta_0){
  active_inds = which(beta_0 != 0)
  CIs = fit$CIs[active_inds,]
  length = mean(CIs[,2] - CIs[,1])
  return(length)
}
inactive_length = function(fit, beta_0){
  active_inds = which(beta_0 != 0)
  CIs = fit$CIs[-active_inds,]
  length = mean(CIs[,2] - CIs[,1])
  return(length)
}
####

#### Metrics for DD Analysis ####
MAE = function(fit, truth){
  est = fit$beta_hat
  return(abs(est-truth))
}

CI_hit = function(fit, truth){
  ind = 1*(fit$CI[1] <= truth && truth <= fit$CI[2])
  return(ind)
}

int_length = function(fit, truth){
  return(fit$CI[2] - fit$CI[1])
}

# k dimensional metrics
L2 = function(fit, truth){
  est = fit$beta_hat
  k = length(est)
  beta_0_k = truth[1:k]
  return( sqrt(sum( (est-beta_0_k)**2) ) ) 
}

kHit = function(fit, truth){
  # Note, can cope with first diagonal elements of the covariance being 0, then just checks beta_hat
  # in those coordinates is equal to the truth, and performs the routine on the submatrix
  est = fit$beta_hat
  k = length(est)
  beta_0_k = truth[1:k]
  chi_zval = qchisq(0.95, df = k)
  mat = fit$cov_hat

  if('add_length' %in% names(fit)){
    add_length = fit$add_length
    cov=cov.adj(mat, add_length, chi_zval)
  }else{
    cov=mat
  }

  zero_inds = which(diag(mat) == 0)
  non_zero_inds = which(diag(mat) != 0)

  if(length(non_zero_inds) == k){
    test_val = t(beta_0_k - est) %*% solve(cov) %*% (beta_0_k - est) 
    hit = test_val <= chi_zval
  }else if(length(zero_inds) == k){
    hit = all(est == beta_0_k)
  }else{
    zero_equality = all(est[zero_inds] == beta_0_k[zero_inds])
    sub_cov = cov[non_zero_inds, non_zero_inds]
    sub_est = est[non_zero_inds]
    sub_beta_0 = beta_0_k[non_zero_inds]
    test_val = t(sub_beta_0 - sub_est) %*% solve(sub_cov) %*% (sub_beta_0 - sub_est) 
    sub_hit = test_val <= chi_zval
    hit = sub_hit && zero_equality
  }
  return(hit)
}

volume = function(fit){
  k = length(fit$beta_hat)
  mat = fit$cov_hat
  chi_zval = qchisq(0.95, df = k)
  if('add_length' %in% names(fit)){
    add_length = fit$add_length
    cov=cov.adj(mat, add_length, chi_zval)
  }else{
    cov=mat
  }
  svd_cov = svd(cov)
  vol = sqrt(prod(svd_cov$d))
  return(vol)
}
####

#### Single Deconfounding Procedures ####
fit_cbm = function(X, Y, transform=trim_transform, oracle_lambda=FALSE){
  '
  Fits spectral deconfounding method from Cevid et. al, with F given by
  transform(X).
  '
  t1 = Sys.time()
  F_mat = transform(X)
  X_tilde = F_mat %*% X
  Y_tilde = F_mat %*% Y
  
  # Estimate beta using lasso on tilde matrices 
  if(oracle_lambda){
    n = dim(X_tilde)[1]
    p = dim(X_tilde)[2]
    lambda_min =  sqrt(log(p)/n)
  }else{
    cv_fit = cv.glmnet(X_tilde, Y_tilde)
    lambda_min = cv_fit$lambda.min
  }
  best_fit = glmnet(X_tilde, Y_tilde, lambda = lambda_min)
  beta_hat = best_fit$beta
  t2 = Sys.time()
  return(list(beta_hat = beta_hat, time = as.numeric(difftime(t2, t1, units = 'secs'))))
}
fit_vb_posterior = function(X, Y, transform=NULL, credible_intervals=TRUE, first_coord=FALSE){
  '
  Computes a variational approximation to the posterior, with L given by the
  transform. If no transform is given, no transform is applied (and) so this
  is the variational approximation to the spike-and-slab with no deconfounding.
  '
  n = dim(X)[1]
  p = dim(X)[2]
  t1 = Sys.time()
  if (!missing(transform)){
    F_mat = transform(X)
    X_tilde = F_mat %*% X
    Y_tilde = F_mat %*% Y
    svb_fit = svb.fit(X_tilde, Y_tilde)  
  }else{
    svb_fit = svb.fit(X, Y)  
  }
  # svb_mean = svb_fit$mu * svb_fit$gamma
  svb_mean = svb_fit$mu * (svb_fit$gamma > 1/2)*1
  if(credible_intervals){
    if(!first_coord){
      n_samples = 10000
      var_samples = sample_from_VB_posterior(n_samples, svb_fit$mu, svb_fit$sigma, svb_fit$gamma)
      CIs = t(apply(var_samples, 2, function(v) quantile(v, probs = c(0.025, 0.975))))
      t2 = Sys.time()
      return(list(beta_hat = svb_mean,
                CIs = CIs,
                time = as.numeric(difftime(t2, t1, units = 'secs'))))  
    }else{
      beta_hat = svb_fit$mu[1]*svb_fit$gamma[1]
      CI = c(beta_hat - 1.96*svb_fit$sigma[1],  beta_hat + 1.96*svb_fit$sigma[1])
      t2 = Sys.time()
      return(list(beta_hat = beta_hat,
                CI = CI,
                fit_time = as.numeric(difftime(t2, t1, units = 'secs'))))  
    }
    
  }
  t2 = Sys.time()
  return(list(beta_hat = svb_mean, time = as.numeric(difftime(t2, t1, units = 'secs'))))
}

fit_mcmc_posterior = function(X, Y, transform=NULL, burnin=1000, n_samples=5000, credible_intervals=FALSE){
  '
  Samples from the posterior using BoomSpikeSlab (MCMC), with L given by the
  transform. If no transform is given, no transform is applied (and) so this
  is the the spike-and-slab posterior with no deconfounding.
  NOTE: This function achieves the same as fit_svb above
  if transform=lava_transform.
  '
  n = dim(X)[1]
  p = dim(X)[2]
  t1 = Sys.time()
  if (!missing(transform)){
    F_mat = transform(X)
    X_tilde = F_mat %*% X
    Y_tilde = F_mat %*% Y
  }else{
    X_tilde = X
    Y_tilde = Y
  }
  
  # Below using our implementation
  z_init = rep(0, p)
  beta_init = rep(0, p)
  mcmc_chain = sample_mc_cpp(burnin + n_samples, beta_init, z_init, X_tilde, Y_tilde)
  beta_samples = mcmc_chain$beta[-c(1:burnin),]
  z_samples = mcmc_chain$z[-c(1:burnin),]
  
  # use this vector as the inclusion probabilities 
  z_prop = apply(z_samples, 2, mean)
  inclusion = (z_prop > 1/2)*1
  
  beta_hat_samples = t(apply(beta_samples, 1, function(v) v*inclusion))
  beta_hat = apply(beta_hat_samples, 2, mean)
  if(credible_intervals){

    CIs = t(apply(beta_samples*z_samples, 2, function(v) quantile(v, probs = c(0.025, 0.975))))
    t2 = Sys.time()
    return(list(beta_hat = beta_hat,
                CIs = CIs,
                time = as.numeric(difftime(t2, t1, units = 'secs'))))
  }
  
  t2 = Sys.time()
  return(list(beta_hat = beta_hat, time = as.numeric(difftime(t2, t1, units = 'secs'))))
}
####

## Doubly Debiased Procedures ##
DD.bd.fit = function(X, Y, transform=trim_transform, n_samples = 1000, lambda = NA, k = 1){
  # Fit the I-SVB method on the first k coordinates.
  t_1 = Sys.time()
  n = dim(X)[1]
  p = dim(X)[2]

  # Execute first debiasing step based on transform supplied - independent of k
  if (!missing(transform)){
    F_mat = transform(X)
    X = F_mat %*% X
    Y = F_mat %*% Y
  }else{
    # redundant lines but placed here for clarity
    X = X
    Y = Y
  }

  # Execute second debiasing procedure and return quantities to perform UQ on beta_{1:k}.
  # Note - here this is done based on the I-SVB method of our other paper.
  if(k == 1){
    # Treat k == 1 case separately because here we will use quantiles for cred interval
    X1 = X[,1]
    X1_norm_sq = sum(X1 * X1)
    H = X1 %*% t(X1) / X1_norm_sq
    #and the projection matrix onto span(X1)^perp
    I_minus_H = diag(rep(1,n)) - H

    #make the matrix P, consisting of basis vectors of span(X1)^perp
    svd_temp = svd(I_minus_H)
    U = svd_temp$u[,1:(n-1)]
    P = t(U)

    #make W_check and Y_check
    W_check = P %*% I_minus_H %*% X[,2:p]
    Y_check = P %*% I_minus_H %*% Y
    #apply svb package to the check model. Note can specify lambda via prior_scale arg
    if(is.na(lambda)){
      vbL_W<-svb.fit(W_check,Y_check,tol = 10e-5, max_iter = 10^4, family ="linear")  
    }else{
      vbL_W<-svb.fit(W_check,Y_check,tol = 10e-5, max_iter = 10^4, family ="linear",
                       lambda = lambda)
    }
    mu_hat = vbL_W$mu
      sigmas_hat = abs(vbL_W$sigma)
    gammas_hat = vbL_W$gamma
    #sample from the variational posterior of beta_{-1}
    beta_minus_1_samples = sample_from_VB_posterior(n_samples,
                                                    mu = mu_hat,
                                                    sigma = sigmas_hat,
                                                    gamma = gammas_hat)
    
    
    #compute gammas (from DY paper, so call them gamma_prime)
    gamma_prime = apply(X[,2:p], 2, function(x) sum(X1*x)/X1_norm_sq)
    #diffs is a temporary variable which we subtract from beta_1^* to get beta_1
    diffs = beta_minus_1_samples %*% gamma_prime
    
    # improper specific part
    posterior_means_beta_1 = rep(1/X1_norm_sq * t(X1) %*% Y, n_samples) - diffs
      posterior_V_beta_1 = 1/X1_norm_sq
      posterior_samples_beta_1 = rnorm(n_samples, posterior_means_beta_1, 
                      sd=sqrt(posterior_V_beta_1))
      # end improper specific part

      beta_hat = mean(posterior_samples_beta_1)
      credible_interval = quantile(posterior_samples_beta_1, probs = c(0.025, 0.975))
      t_2 = Sys.time()
      return(list(beta_hat=beta_hat,
            CI=as.numeric(credible_interval),
            fit_time=as.numeric(difftime(t_2, t_1, units = 'secs')))) 
  }else{
    A_k = X[,1:k]
    L = chol(t(A_k)%*% A_k)
    H = A_k %*% solve(t(A_k)%*%A_k, t(A_k))
    I_minus_H = diag(rep(1,n)) - H
    #make the matrix P, consisting of basis vectors of span(X1)^perp
    svd_temp = svd(I_minus_H)
    U = svd_temp$u[,1:(n-k)]
    P = t(U)
    W_check = P %*% I_minus_H %*% X[,(k+1):p]
    Y_check = P %*% I_minus_H %*% Y
    #apply svb package to the check model. Note can specify lambda via prior_scale arg
    if(is.na(lambda)){
      vbL_W<-svb.fit(W_check,Y_check,tol = 10e-5, max_iter = 10^4, family ="linear")  
    }else{
      vbL_W<-svb.fit(W_check,Y_check,tol = 10e-5, max_iter = 10^4, family ="linear",
                     lambda = lambda)
    }
    #extract relevant params from the fit
    mu_hat = vbL_W$mu
    sigmas_hat = abs(vbL_W$sigma)
    gammas_hat = vbL_W$gamma
    beta_minus_k_samples = sample_from_VB_posterior(n_samples,
                                                    mu = mu_hat,
                                                    sigma = sigmas_hat,
                                                    gamma = gammas_hat)
    
    
    Gamma_mat = solve(t(A_k)%*%A_k, t(A_k)) %*% X[,(k+1):p]
    
    diffs = beta_minus_k_samples %*% t(Gamma_mat)

    # Improper specific part
    Sigma_p = solve(t(A_k) %*% A_k, diag(k))  
    Mu_p = Sigma_p %*% t(A_k) %*% Y
    beta_star_k_samples = rmvnorm(n_samples, mean = Mu_p, sigma = Sigma_p)
    beta_k_samples = beta_star_k_samples - diffs
    # Improper specific part end

    beta_hat = apply(beta_k_samples, 2, mean)
    cov_hat = empirical_covariance(beta_k_samples)
    cov_hat_beta_star = Sigma_p 
    cov_hat_diffs = Gamma_mat %*% diag(sigmas_hat**2 * gammas_hat**2) %*% t(Gamma_mat)

    t_2 = Sys.time()
    return(list(beta_hat=beta_hat,
          cov_hat=cov_hat,
          cov_hat_beta_star = cov_hat_beta_star,
          cov_hat_nuisance = cov_hat_diffs,
          fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))
  }
}

DD.gcb.fit = function(X, Y, transform=trim_transform, oracle_lambda=FALSE, k = 1){
  # Fit the method suggested by GCB (Doubly Debiased Lasso)
  # Fit the ZZ method
  if(k != 1){
    stop('Only implemented for k = 1')
  }
  n = dim(X)[1]
  p = dim(X)[2]
  t_1 = Sys.time()

  # First debiasing step based on the transform
  F_mat = transform(X)
  X_tilde = F_mat %*% X
  Y_tilde = F_mat %*% Y

  # Zhang and Zhang debiasing procedure
  gamma = 0.95
  cv_fit = cv.glmnet(X, Y)
  beta_lasso = glmnet(X, Y,lambda =cv_fit$lambda.min)$beta
  X_minus_1 = X[,-1]
  cv_fit_gamma = cv.glmnet(x = X_minus_1, y = X[,1])
  gamma_hat = glmnet(x = X_minus_1, y = X[,1], lambda = cv_fit_gamma$lambda.min)$beta
  
  z_1 = as.matrix(X[,1] - X_minus_1 %*% gamma_hat)
  beta_hat = as.numeric(beta_lasso[1] + t(z_1)%*%(Y - X%*%beta_lasso)/(t(z_1)%*%X[,1]))
  CI = c(beta_hat - qnorm((1+gamma)/2)/sqrt(t(z_1) %*% z_1),
       beta_hat + qnorm((1+gamma)/2)/sqrt(t(z_1) %*% z_1))
  t_2 = Sys.time()
  return(list(beta_hat=beta_hat,
          CI=CI,
          fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))
}

DDL.fit = function(X, Y, k = 1){
  if(k != 1){
    stop('Only implemented for k = 1')
  }
  n = dim(X)[1]
  p = dim(X)[2]
  t_1 = Sys.time()
  ddl_fit = DDL(X, Y, index=1)
  est_ddl = ddl_fit$est_ddl
  se = ddl_fit$se

  beta_hat = est_ddl
  CI = c(est_ddl - qnorm(1 - 0.05 / 2)*se, est_ddl + qnorm(1 - 0.05 / 2)*se)

  t_2 = Sys.time()
  return(list(beta_hat=beta_hat,
          CI=CI,
          fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))
}



# for refererence
isvb.fit = function(X, Y, n_samples = 1000, lambda = NA, k = 1){
  # Fit the I-SVB method on the first k coordinates.
  t_1 = Sys.time()
  n = dim(X)[1]
  p = dim(X)[2]
  if(k == 1){
    # Treat k == 1 case separately because here we will use quantiles for cred interval
    X1 = X[,1]
    X1_norm_sq = sum(X1 * X1)
    H = X1 %*% t(X1) / X1_norm_sq
    #and the projection matrix onto span(X1)^perp
    I_minus_H = diag(rep(1,n)) - H

    #make the matrix P, consisting of basis vectors of span(X1)^perp
    svd_temp = svd(I_minus_H)
    U = svd_temp$u[,1:(n-1)]
    P = t(U)

    #make W_check and Y_check
    W_check = P %*% I_minus_H %*% X[,2:p]
    Y_check = P %*% I_minus_H %*% Y
    #apply svb package to the check model. Note can specify lambda via prior_scale arg
    if(is.na(lambda)){
      vbL_W<-svb.fit(W_check,Y_check,tol = 10e-5, max_iter = 10^4, family ="linear")  
    }else{
      vbL_W<-svb.fit(W_check,Y_check,tol = 10e-5, max_iter = 10^4, family ="linear",
                       lambda = lambda)
    }
    mu_hat = vbL_W$mu
      sigmas_hat = abs(vbL_W$sigma)
    gammas_hat = vbL_W$gamma
    #sample from the variational posterior of beta_{-1}
    beta_minus_1_samples = sample_from_VB_posterior(n_samples,
                                                    mu = mu_hat,
                                                    sigma = sigmas_hat,
                                                    gamma = gammas_hat)
    
    
    #compute gammas (from DY paper, so call them gamma_prime)
    gamma_prime = apply(X[,2:p], 2, function(x) sum(X1*x)/X1_norm_sq)
    #diffs is a temporary variable which we subtract from beta_1^* to get beta_1
    diffs = beta_minus_1_samples %*% gamma_prime
    
    # improper specific part
    posterior_means_beta_1 = rep(1/X1_norm_sq * t(X1) %*% Y, n_samples) - diffs
      posterior_V_beta_1 = 1/X1_norm_sq
      posterior_samples_beta_1 = rnorm(n_samples, posterior_means_beta_1, 
                      sd=sqrt(posterior_V_beta_1))
      # end improper specific part

      beta_hat = mean(posterior_samples_beta_1)
      credible_interval = quantile(posterior_samples_beta_1, probs = c(0.025, 0.975))
      t_2 = Sys.time()
      return(list(beta_hat=beta_hat,
            CI=as.numeric(credible_interval),
            fit_time=as.numeric(difftime(t_2, t_1, units = 'secs')))) 
  }else{
    A_k = X[,1:k]
    L = chol(t(A_k)%*% A_k)
    H = A_k %*% solve(t(A_k)%*%A_k, t(A_k))
    I_minus_H = diag(rep(1,n)) - H
    #make the matrix P, consisting of basis vectors of span(X1)^perp
    svd_temp = svd(I_minus_H)
    U = svd_temp$u[,1:(n-k)]
    P = t(U)
    W_check = P %*% I_minus_H %*% X[,(k+1):p]
    Y_check = P %*% I_minus_H %*% Y
    #apply svb package to the check model. Note can specify lambda via prior_scale arg
    if(is.na(lambda)){
      vbL_W<-svb.fit(W_check,Y_check,tol = 10e-5, max_iter = 10^4, family ="linear")  
    }else{
      vbL_W<-svb.fit(W_check,Y_check,tol = 10e-5, max_iter = 10^4, family ="linear",
                     lambda = lambda)
    }
    #extract relevant params from the fit
    mu_hat = vbL_W$mu
    sigmas_hat = abs(vbL_W$sigma)
    gammas_hat = vbL_W$gamma
    beta_minus_k_samples = sample_from_VB_posterior(n_samples,
                                                    mu = mu_hat,
                                                    sigma = sigmas_hat,
                                                    gamma = gammas_hat)
    
    
    Gamma_mat = solve(t(A_k)%*%A_k, t(A_k)) %*% X[,(k+1):p]
    
    diffs = beta_minus_k_samples %*% t(Gamma_mat)

    # Improper specific part
    Sigma_p = solve(t(A_k) %*% A_k, diag(k))  
    Mu_p = Sigma_p %*% t(A_k) %*% Y
    beta_star_k_samples = rmvnorm(n_samples, mean = Mu_p, sigma = Sigma_p)
    beta_k_samples = beta_star_k_samples - diffs
    # Improper specific part end

    beta_hat = apply(beta_k_samples, 2, mean)
    cov_hat = empirical_covariance(beta_k_samples)
    cov_hat_beta_star = Sigma_p 
    cov_hat_diffs = Gamma_mat %*% diag(sigmas_hat**2 * gammas_hat**2) %*% t(Gamma_mat)

    t_2 = Sys.time()
    return(list(beta_hat=beta_hat,
          cov_hat=cov_hat,
          cov_hat_beta_star = cov_hat_beta_star,
          cov_hat_nuisance = cov_hat_diffs,
          fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))
  }
}
####

#### Single Debiased Estimation Loops ####
sample_loss = function(n, p, s0, q, signal_size=NA, fits='all',
                       mcmc_burn=1000, mcmc_samples=5000, rho=0, dataset='generated', AR=FALSE, confounding_layers=1){
    
  if(dataset=='cbm'){
    dat = make_data_CBM(n, p, s0, q)  
    oracle_lambda=FALSE
    cbm_lambda=TRUE
  }else if(dataset == 'generated'){
    dat = make_data(n, p, s0, q, signal_size=signal_size, AR=AR, confounding_layers=confounding_layers)  
    oracle_lambda=FALSE
    cbm_lambda=FALSE
  }else if(dataset == 'riboflavin'){
    dat = make_data(n, p, s0, q, signal_size=signal_size, AR=AR)
    oracle_lambda=FALSE
    cbm_lambda=FALSE
  }else{
    stop('Dataset must be one of {cbm, generated, riboflavin}')
  }
  X = dat$X
  Y = dat$Y
  beta_0 = dat$beta_0
  if(fits == 'all'){
    fits = list(CBM_trim = fit_cbm(X, Y, transform = trim_transform, oracle_lambda),
                CBM_lava = fit_cbm(X, Y, transform = function(x) lava_transform(x, cbm_lambda=cbm_lambda), oracle_lambda),
                BD_trim_vb = fit_vb_posterior(X, Y, transform = trim_transform),
                BD_lava_vb = fit_vb_posterior(X, Y, transform = lava_transform),
                BD_trim_mcmc = fit_mcmc_posterior(X, Y, transform = trim_transform,
                                                  burnin = mcmc_burn, n_samples = mcmc_samples),
                BD_lava_mcmc = fit_mcmc_posterior(X, Y, transform = lava_transform,
                                                  burnin = mcmc_burn, n_samples = mcmc_samples),
                SAS_vb = fit_vb_posterior(X, Y),
                SAS_mcmc = fit_mcmc_posterior(X, Y,
                                              burnin = mcmc_burn, n_samples = mcmc_samples)
    )
  }else if(fits == 'mcmc'){
    fits = list(BD_trim_mcmc = fit_mcmc_posterior(X, Y, transform = trim_transform,
                                                  burnin = mcmc_burn, n_samples = mcmc_samples),
                BD_lava_mcmc = fit_mcmc_posterior(X, Y, transform = lava_transform,
                                                  burnin = mcmc_burn, n_samples = mcmc_samples),
                SAS_mcmc = fit_mcmc_posterior(X, Y,
                                              burnin = mcmc_burn, n_samples = mcmc_samples)
    )
  }else if(fits == 'no_mcmc'){
    fits = list(CBM_trim = fit_cbm(X, Y, transform = trim_transform, oracle_lambda),
                CBM_lava = fit_cbm(X, Y, transform = function(x) lava_transform(x, cbm_lambda=cbm_lambda), oracle_lambda),
                # CBM_ORT = fit_cbm(X, Y, transform = function(x) oracle_removal_transform(x, q=q)),
                BD_trim_vb = fit_vb_posterior(X, Y, transform = trim_transform),
                BD_lava_vb = fit_vb_posterior(X, Y, transform = lava_transform),
                # BD_trim_ORT = fit_vb_posterior(X, Y, function(x) oracle_removal_transform(x, q=q)),
                SAS_vb = fit_vb_posterior(X, Y)) 
  }else if(fits == 'cbm_comp'){
    fits = list(CBM_trim = fit_cbm(X, Y, transform = trim_transform, oracle_lambda),
                CBM_lava = fit_cbm(X, Y, transform = function(x) lava_transform(x, cbm_lambda=cbm_lambda), oracle_lambda),
                # CBM_puffer = fit_cbm(X, Y, transform = puffer_transform, oracle_lambda),
                LASSO = fit_cbm(X, Y, transform = identity_transform, oracle_lambda),
                BD_trim_vb = fit_vb_posterior(X, Y, transform = trim_transform),
                BD_lava_vb = fit_vb_posterior(X, Y, transform = lava_transform),
                SAS_vb = fit_vb_posterior(X, Y)) 
  }
  
  # Compute Metrics
  l2_loss = data.frame(lapply(fits, function(fit) l2_dist(fit$beta_hat, beta_0)))
  l1_loss = data.frame(lapply(fits, function(fit) l1_dist(fit$beta_hat, beta_0)))
  precision_score = data.frame(lapply(fits, function(fit) precision(fit$beta_hat, beta_0)))
  recall_score = data.frame(lapply(fits, function(fit) recall(fit$beta_hat, beta_0)))
  f1_score = data.frame(lapply(fits, function(fit) f1(fit$beta_hat, beta_0)))
  times = data.frame(lapply(fits, function(fit) fit$time))
  
  return(list(l2_loss = l2_loss, l1_loss = l1_loss,
              precision = precision_score, recall = recall_score, f1 = f1_score,
              time = times))
}

estimate_losses = function(n_replicates, n, p, s0, q, signal_size=NA, fits='all', rho=0, dataset='generated', AR=FALSE, confounding_layers=1,
                           mc.cores=1){
  cat('\n---------------------------',
    '\nRunning experiment with:',
    '\nn: ',n,
    '\np: ', p,'\ns0: ',s0,'\nq: ',q,
    '\nsignal_size: ',signal_size,
    '\nfeature_correlation: ', rho,'\nn_replicates: ',n_replicates,
    '\ncores: ',mc.cores,
    '\n----------------------------\n')
  # replicated = pbreplicate(n_replicates, sample_loss(n, p, s0, q, signal_size=signal_size, 
                                                     # fits = fits, rho = rho, scheme=scheme))
  if(mc.cores > 6){
    stop('More than 6 cores, commment this check to proceed.')
  }                                                     
  replicated = mc_replicate(n_replicates,
   sample_loss(n, p, s0, q, signal_size=signal_size, fits = fits, rho = rho, dataset=dataset, AR=AR, confounding_layers=confounding_layers),
   mc.cores=mc.cores
  )                                                     
  bound = apply(replicated, 1, function(x) do.call(rbind, x))
  means = lapply(bound, function(x) apply(x, 2, mean))
  sds = lapply(bound, function(x) apply(x, 2, sd))
  return(list(mean = means, sd = sds))
}
sample_coverage_old = function(n, p, s0, q, signal_size=NA, AR=FALSE,
                       mcmc_burn=1000, mcmc_samples=5000, fits='all', rho=0.0, confounding_layers=1){
  dat = make_data(n, p, s0, q, signal_size=signal_size, rho=rho, AR=AR, confounding_layers=confounding_layers)
  X = dat$X
  Y = dat$Y
  beta_0 = dat$beta_0
  if(fits == 'all'){
    fits = list(BD_trim_vb = fit_vb_posterior(X, Y,
                                              transform = trim_transform, credible_intervals = TRUE),
                BD_lava_vb = fit_vb_posterior(X, Y,
                                              transform = lava_transform, credible_intervals = TRUE),
                BD_trim_mcmc = fit_mcmc_posterior(X, Y,
                                                  transform = trim_transform,
                                                  burnin = mcmc_burn, n_samples = mcmc_samples,
                                                  credible_intervals = TRUE),
                BD_lava_mcmc = fit_mcmc_posterior(X, Y,
                                                  transform = lava_transform,
                                                  burnin = mcmc_burn, n_samples = mcmc_samples,
                                                  credible_intervals = TRUE),
                SAS_vb = fit_vb_posterior(X, Y, credible_intervals = TRUE),
                SAS_mcmc = fit_mcmc_posterior(X, Y,
                                              burnin = mcmc_burn, n_samples = mcmc_samples,
                                              credible_intervals = TRUE)
    )
  }else if(fits == 'with_ort'){
    fits = list(BD_trim_vb = fit_vb_posterior(X, Y,
                                              transform = trim_transform, credible_intervals = TRUE),
                BD_lava_vb = fit_vb_posterior(X, Y,
                                              transform = lava_transform, credible_intervals = TRUE),
                BD_ort_vb = fit_vb_posterior(X, Y,
                                              transform = function(x) oracle_removal_transform(x, q=q),
                                             credible_intervals = TRUE),
                BD_trim_mcmc = fit_mcmc_posterior(X, Y,
                                                  transform = trim_transform,
                                                  burnin = mcmc_burn, n_samples = mcmc_samples,
                                                  credible_intervals = TRUE),
                BD_lava_mcmc = fit_mcmc_posterior(X, Y,
                                                  transform = lava_transform,
                                                  burnin = mcmc_burn, n_samples = mcmc_samples,
                                                  credible_intervals = TRUE),
                BD_ort_mcmc = fit_mcmc_posterior(X, Y,
                                             transform = function(x) oracle_removal_transform(x, q),
                                             burnin = mcmc_burn, n_samples = mcmc_samples,
                                             credible_intervals = TRUE)
    )
  }else{
    fits = list(BD_trim_vb = fit_vb_posterior(X, Y,
                                              transform = trim_transform, credible_intervals = TRUE),
                BD_lava_vb = fit_vb_posterior(X, Y,
                                              transform = lava_transform, credible_intervals = TRUE),
                SAS_vb = fit_vb_posterior(X, Y, credible_intervals = TRUE)
    )
  }
  # Compute Metrics
  l2_loss = data.frame(lapply(fits, function(fit) l2_dist(fit$beta_hat, beta_0)))
  l1_loss = data.frame(lapply(fits, function(fit) l1_dist(fit$beta_hat, beta_0)))
  precision_score = data.frame(lapply(fits, function(fit) precision(fit$beta_hat, beta_0)))
  recall_score = data.frame(lapply(fits, function(fit) recall(fit$beta_hat, beta_0)))
  f1_score = data.frame(lapply(fits, function(fit) f1(fit$beta_hat, beta_0)))
  times = data.frame(lapply(fits, function(fit) fit$time))
  
  active_hit = data.frame(lapply(fits, function(fit) active_hit(fit, beta_0)))
  inactive_hit = data.frame(lapply(fits, function(fit) inactive_hit(fit, beta_0)))
  active_length = data.frame(lapply(fits, function(fit) active_length(fit, beta_0)))
  inactive_length = data.frame(lapply(fits, function(fit) inactive_length(fit, beta_0)))
  
  return(list(l2_loss = l2_loss, l1_loss = l1_loss,
              precision = precision_score, recall = recall_score, f1 = f1_score,
              time = times, active_hit = active_hit, inactive_hit = inactive_hit,
              active_length = active_length, inactive_length = inactive_length))
}

sample_coverage = function(n, p, s0, q, fits, signal_size=NA, AR=FALSE,
                       mcmc_burn=1000, mcmc_samples=5000, rho=0.0, confounding_layers=1, dataset='generated'){
  dat = make_data(n, p, s0, q, signal_size=signal_size, rho=rho, AR=AR, confounding_layers=confounding_layers, dataset=dataset)
  X = dat$X
  Y = dat$Y
  beta_0 = dat$beta_0
  fits = lapply(fits, function(fit) fit(X, Y))
  
  # Compute Metrics
  l2_loss = data.frame(lapply(fits, function(fit) l2_dist(fit$beta_hat, beta_0)))
  l1_loss = data.frame(lapply(fits, function(fit) l1_dist(fit$beta_hat, beta_0)))
  precision_score = data.frame(lapply(fits, function(fit) precision(fit$beta_hat, beta_0)))
  recall_score = data.frame(lapply(fits, function(fit) recall(fit$beta_hat, beta_0)))
  f1_score = data.frame(lapply(fits, function(fit) f1(fit$beta_hat, beta_0)))
  times = data.frame(lapply(fits, function(fit) fit$time))
  
  active_hit = data.frame(lapply(fits, function(fit) active_hit(fit, beta_0)))
  inactive_hit = data.frame(lapply(fits, function(fit) inactive_hit(fit, beta_0)))
  active_length = data.frame(lapply(fits, function(fit) active_length(fit, beta_0)))
  inactive_length = data.frame(lapply(fits, function(fit) inactive_length(fit, beta_0)))
  
  return(list(l2_loss = l2_loss, l1_loss = l1_loss,
              precision = precision_score, recall = recall_score, f1 = f1_score,
              time = times, active_hit = active_hit, inactive_hit = inactive_hit,
              active_length = active_length, inactive_length = inactive_length))
}

estimate_coverage = function(n_replicates, n, p, s0, q, signal_size=NA, fits = 'all', rho=0, AR=FALSE, confounding_layers=1, mc.cores=1, dataset='generated'){
  cat('\n---------------------------',
    '\nRunning experiment with:',
    '\nn: ',n,
    '\np: ', p,'\ns0: ',s0,'\nq: ',q,
    '\nsignal_size: ',signal_size,
    '\nfeature_correlation: ', rho,'\nn_replicates: ',n_replicates,
    '\ncores: ',mc.cores,
    '\n----------------------------\n')
  if(mc.cores > 6){
    stop('More than 6 cores, commment this check to proceed.')
  }
  replicated = mc_replicate(n_replicates,
   sample_coverage(n, p, s0, q, signal_size = signal_size, fits = fits, rho = rho, AR=AR, confounding_layers=confounding_layers, dataset=dataset),
   mc.cores=mc.cores
   )
  # replicated = pbreplicate(n_replicates, sample_coverage(n, p, s0, q, signal_size = signal_size, fits = fits, rho = rho))
  bound = apply(replicated, 1, function(x) do.call(rbind, x))
  means = lapply(bound, function(x) apply(x, 2, mean))
  sds = lapply(bound, function(x) apply(x, 2, sd))
  return(list(mean = means, sd = sds))
}
####


row_join = function(t_row){
  '
  Format output of estimate_losses for use in a table
  '
  out = do.call(c, lapply(t_row,
                          function(x) do.call(c, x)))
  # if want to reorder columns use [c(rbind(c(1:6), c(7:12)))]
  return(out)
}
sample_from_VB_posterior = function(n, mu, sigma, gamma){
  p = length(mu)
  samples = sapply(c(1:p), function(i) (runif(n) <= gamma[i])*rnorm(n,
                                                                    mean = mu[i],
                                                                    sd = sigma[i]))
  return(samples)
}

#### Double Debiased Estimation Loops ####
sample_coverage_DD = function(n, p, s0, q, signal_size=NA,
                              fits=list(DD.BD = DD.bd.fit, DD.GCB = DD.gcb.fit),
                              rho=0.0, AR=FALSE, dataset='generated', index=1, noise_var=1, confounding_layers=1, k=1){
  dat = make_data(n, p, s0, q, signal_size=signal_size, rho=rho, noise_var=noise_var, AR=AR, dataset=dataset, index=index, confounding_layers=confounding_layers)
  X = dat$X
  Y = dat$Y
  beta_0 = dat$beta_0

  if(noise_var != 1){
    mod = cv.glmnet(X, Y)
    fit = glmnet(X, Y, lambda = mod$lambda.min)
    yhat = X %*% fit$beta
    sigma_hat = sqrt((1/(n - s0)) * sum((yhat - Y)**2) )
    X_scaled = X/sigma_hat
    Y_scaled = Y/sigma_hat  
  }else{
    X_scaled = X
    Y_scaled = Y
  }

  fitted = lapply(fits, function(fn) fn(X_scaled, Y_scaled, k=k))
  
  if(k==1){
    fit_maes = data.frame(lapply(fitted, function(fit) MAE(fit, beta_0[1])))
    fit_hits = data.frame(lapply(fitted, function(fit) CI_hit(fit, beta_0[1])))
    fit_lengths = data.frame(lapply(fitted, function(fit) int_length(fit, beta_0[1])))
    fit_times = data.frame(lapply(fitted, function(fit) fit$fit_time))  
    return(list('hit'=fit_hits, 'mae'=fit_maes, 'length'=fit_lengths, 'time' = fit_times))
  }else{
    fit_l2s = data.frame(lapply(fitted, function(fit) L2(fit, beta_0)))
    fit_hits = data.frame(lapply(fitted, function(fit) kHit(fit, beta_0)))
    fit_vols = data.frame(lapply(fitted, function(fit) volume(fit)))
    fit_times = data.frame(lapply(fitted, function(fit) fit$fit_time))  
    return(list('hit'=fit_hits, 'l2'=fit_l2s, 'volume'=fit_vols, 'time' = fit_times))
  }
  
}

estimate_coverage_dd = function(n_replicates, n, p, s0, q, signal_size=NA,
                                fits=list(DD.BD = DD.bd.fit, DD.GCB = DD.gcb.fit),
                                rho=0, AR=FALSE, confounding_layers=1, dataset='generated', index=1, k=1, mc.cores=1){
  cat('\n---------------------------',  
    '\nRunning experiment with:',
    '\nn: ',n,
    '\np: ', p,'\ns0: ',s0,'\nq: ',q,
    '\nsignal_size: ',signal_size,
    '\nfeature_correlation: ', rho,'\nconfounding_layers: ', confounding_layers,'\nn_replicates: ',n_replicates,
    '\nindex: ', index,
    '\ncores: ',mc.cores,
    '\n----------------------------\n')
  if(mc.cores > 6){
    stop('More than 6 cores, commment this check to proceed.')
  }
  reps = mc_replicate(n_replicates,
   sample_coverage_DD(n, p, s0, q, signal_size=signal_size, fits=fits, rho=rho, AR=AR, confounding_layers=confounding_layers, dataset=dataset, index=index, k=k),
   mc.cores=mc.cores
   )
  print('Finished replicates.')
  if(k == 1){
    metrics =     c('hit', 'mae', 'length','time')
    compute_sd =  c(FALSE,  TRUE, TRUE,  TRUE)
  }else{
    metrics =     c('hit', 'l2', 'volume','time')
    compute_sd =  c(FALSE,  TRUE, TRUE,  TRUE)
  }
  results=list()
  results_sds=list()
  for(i in c(1:length(metrics))){
    metric=metrics[i]
    met_mean = apply(do.call(rbind, reps[metric,]), 2, mean)
    results = append(results, list(met_mean))
    
    if(compute_sd[i]){
      met_sd = apply(do.call(rbind, reps[metric,]), 2, sd)
      results_sds = append(results_sds, list(met_sd))
    }
  }
  names(results)=metrics
  names(results_sds) = metrics[compute_sd]
  return(list('mean'=results, 'sd'=results_sds))
}

format_results = function(res, round=3, methods=names(fits), k=1){
  mean_df = res$mean
  sd_df = res$sd
  
  if(k == 1){
    metrics =     c('hit', 'mae', 'length','time')
    metrics_sd =  c(FALSE,  TRUE, TRUE,  TRUE)
  }else{
    metrics =     c('hit', 'l2', 'volume','time')
    metrics_sd =  c(FALSE,  TRUE, TRUE,  TRUE)
  }
  formatted_results = methods
  for(i in c(1:length(metrics))){
    metric = metrics[i]
    metric_sd = metrics_sd[i]
    if(metric == 'volume'){
      # Then normalise
      if('oracle' %in% methods){
          normalisation = mean_df[[metric]]['oracle']
          if(normalisation == 0){normalisation = mean_df[[metric]]['isvb']}
          met_col = paste(format(round(mean_df[[metric]]/normalisation,round), round),
                      format(round(sd_df[[metric]]/normalisation,round), round), sep = ' ± ')   
        }else{
          met_col = paste(format(round(mean_df[[metric]]/mean_df[[metric]]['isvb'],round), round),
                      format(round(sd_df[[metric]]/mean_df[[metric]]['isvb'],round), round), sep = ' ± ')   
        }
    }else if(metric_sd){
      met_col = paste(format(round(mean_df[[metric]],round), round),
                      format(round(sd_df[[metric]],round), round), sep = ' ± ')
    }else{
      met_col = format(round(mean_df[[metric]],round), round)
    }
    formatted_results = cbind(formatted_results, met_col)
  }
  formatted_results = as.data.frame(formatted_results)
  colnames(formatted_results) = c('method', metrics)
  rownames(formatted_results) = c()
  return(formatted_results)
}

uq_format_results = function(res, round=3, include_sd=TRUE){
  mean_df = res$mean
  sd_df = res$sd
  metrics = names(mean_df)
  methods = names(res[[1]][[1]])
  formatted_results = methods

  for(i in c(1:length(metrics))){
    metric = metrics[i]
    if(include_sd){
      # met_col = paste(format(round(mean_df[[metric]],round), round),
      #                 format(round(sd_df[[metric]],round), round), sep = ' ± ')
      met_col = paste(format(round(mean_df[[metric]], round), digits = round, nsmall = round),
                      format(round(sd_df[[metric]], round), digits = round, nsmall = round), sep = ' ± ')
    }else{
      met_col = format(round(mean_df[[metric]], round), digits = round, nsmall = round)
    }
    formatted_results = cbind(formatted_results, met_col)
  }
  formatted_results = as.data.frame(t(formatted_results))
  rownames(formatted_results) = c('method', metrics)
  colnames(formatted_results) = c()

  return(t(formatted_results))
}

format_riboflavin_results = function(results, round=3, methods=names(fits), k=1){
  if(k != 1){
    stop('format_riboflavin_results() only implemented with k = 1.')
  }
  mean_data <- lapply(results, function(x) x$mean)
  sd_data <- lapply(results, function(x) x$sd)
  metrics = c('hit', 'mae', 'length','time')
  metrics_sd =  c(FALSE,  TRUE, TRUE,  TRUE)
  mean_results = list()
  sd_results = list()
  for(i in c(1:length(metrics))){
    metric = metrics[i]
    metric_sd  = metrics_sd[i]
    metric_list = lapply(mean_data, function(x) x[[metric]])
    metric_mean = colMeans(do.call(rbind, metric_list))  
    mean_results[[metric]] = metric_mean
    if(metric_sd){
      sd_metric_list = lapply(sd_data, function(x) x[[metric]])
      sd_metric_mean = colMeans(do.call(rbind, sd_metric_list))  
      sd_results[[metric]] = sd_metric_mean
    }
  }

  return(format_results(list('mean'=mean_results, 'sd'=sd_results)))
}