library(mcreplicate)
library(tidyverse)

#### Data Generation Functions ####

make_data_LinConfounding = function(n, p, xs, q, fn_0){

  E = matrix(rnorm(n*p), nrow = n, ncol = p)
  e = rnorm(n, sd=1)

  H = matrix(rnorm(n*q), nrow = n, ncol = q)
  Gamma = matrix(rnorm(q*p), nrow = q, ncol = p)
  delta = rnorm(q)  
  # b = 1/(t(Gamma) %*% Gamma + 1) * t(Gamma) %*% delta

  X = H %*% Gamma + E
  Y = apply(X, 1, fn_0) + H %*% delta + e

  fn_xs = apply(xs, 1, fn_0)

  return(list(X = X, Y = Y, fn_xs = fn_xs))
} 

make_data_NonLinPerturbation = function(n, xs, fn_0, omega){

  gn_0 = function(x){
    return(sin(2*pi*omega*x))
  }  

  e = rnorm(n, sd=1)
  X = matrix(runif(n, min = -2, max = 2), nrow=n)
  Y = apply(X, 1, fn_0) + gn_0(X) + e

  fn_xs = apply(xs, 1, fn_0)

  return(list(X = X, Y = Y, fn_xs = fn_xs))
} 

#### GP Fitting Functions ####

fit_gp = function(X, Y, xs, kernel){
  n = dim(X)[1]
  p = dim(X)[2]
  t1 = Sys.time()
  K_X_X = apply(X, 1, function(x) apply(X, 1, function(y) kernel(x, y)))
  K_xs_X = apply(X, 1, function(x) apply(xs, 1, function(y) kernel(x, y)))
  K_xs_xs = apply(xs, 1, function(x) apply(xs, 1, function(y) kernel(x, y)))
  
  mat_inv = solve(K_X_X + diag(rep(1, n)))
  
  post_mean = K_xs_X %*% mat_inv %*% Y
  post_cov = K_xs_xs - K_xs_X %*% mat_inv %*% t(K_xs_X)
  t2 = Sys.time()
  elapsed = difftime(t2, t1, units='secs')
  return(list(est = post_mean, pointwise_var = diag(post_cov), cov = post_cov, fit_time = elapsed))
}

fit_dgp = function(X, Y, xs, kernel, dSigma=NA, dkernel=NA, type='linear', sing_offset=0){
  n = dim(X)[1]
  p = dim(X)[2]
  t1 = Sys.time()
  K_X_X = apply(X, 1, function(x) apply(X, 1, function(y) kernel(x, y)))
  K_xs_X = apply(X, 1, function(x) apply(xs, 1, function(y) kernel(x, y)))
  K_xs_xs = apply(xs, 1, function(x) apply(xs, 1, function(y) kernel(x, y)))
  if(type=='linear'){
    if(any(is.na(dSigma))){
      dSigma = diag(rep(1, p))
    }
    mat_inv = solve(K_X_X + diag(rep(1, n)) + X %*% dSigma %*% t(X) )  
  }else if(type == 'non-linear'){
    H_X_X = apply(X, 1, function(x) apply(X, 1, function(y) dkernel(x, y))) + diag(sing_offset, n)
    mat_inv = solve(K_X_X + diag(rep(1, n)) + solve(H_X_X) )  
  }else{
    stop('Type must be one of {linear, non-linear}')
  }
  
  post_mean = K_xs_X %*% mat_inv %*% Y
  post_cov = K_xs_xs - K_xs_X %*% mat_inv %*% t(K_xs_X)
  t2 = Sys.time()
  elapsed = difftime(t2, t1, units='secs')
  return(list(est = post_mean, pointwise_var = diag(post_cov), cov = post_cov, fit_time=elapsed))
}

#### GP Kernels ####

mat_kern = function(x, y, l = 1){
  return( exp (-sqrt(sum( (x - y)**2 ))))
}

bm_kern = function(x, y, cn = 1){
  return( cn * min(x, y))
}


#### Losses ####
l1_l = function(fit, f0_xs){
  return( mean(abs(fit$est - f0_xs)) )
}

l2_l = function(fit, f0_xs){
  return( mean((fit$est - f0_xs)**2) )
}

hit_l = function(fit, f0_xs){
  # Assume 95% CIs
  scale = qnorm((1+0.95)/2)
  
  cis = cbind(fit$est - scale*sqrt(fit$pointwise_var), 
              fit$est + scale*sqrt(fit$pointwise_var))
  hits = (cis[,1] <= f0_xs) & (cis[,2] >= f0_xs)

  return( mean(hits) )
}

length_l = function(fit, f0_xs){
  scale = qnorm((1+0.95)/2)
  return(2*scale*sqrt(fit$pointwise_var))
}

time_l = function(fit, f0_xs=NA){
  return(as.numeric(fit$fit_time))
}

#### Experiment Functions ####

sample_loss = function(n, p, xs, make_data, fits){
  dat = make_data(n, p, xs)
  X = dat$X
  Y = dat$Y
  f0_xs = dat$fn_xs

  fitted = lapply(fits, function(fn) fn(X, Y, xs))
  loss_fns = list(l1 = l1_l, l2 = l2_l, hit = hit_l, length = length_l, time=time_l)
  # losses = lapply(fitted, function(fit) lapply(loss_fns, function(loss_fn) loss_fn(fit, f0_xs)))
  losses = lapply(loss_fns, function(loss_fn) lapply(fitted, function(fit) loss_fn(fit, f0_xs)))
  return(losses)
}

estimate_losses = function(n_replicates, n, p, xs, data_fn, fits,
                           mc.cores=1){
  replicated = mc_replicate(n_replicates,
                            sample_loss(n, p, xs, data_fn, fits),
                            mc.cores=mc.cores
  )                
  
  bound = apply(replicated, 1, function(x) do.call(rbind, x)) %>% 
    lapply(function(arr) apply(arr, 2, function(v) do.call(c, v)))
  means = bound %>% lapply(function(arr) apply(arr, 2, mean))
  sds = bound %>% lapply(function(arr) apply(arr, 2, sd))
  return(list(mean = means, sd = sds))
}


format_results = function(results, fit_names, round=3){
  means = lapply(results, function(res) res$mean)
  sds = lapply(results, function(res) res$sd)
  losses = names(means$p1)
  formatted = list()
  for(loss in losses){
    mean_losses = do.call(cbind, lapply(means, function(x) x[[loss]]))
    sd_losses = do.call(cbind, lapply(sds, function(x) x[[loss]]))
    joined = paste(format(round(mean_losses, 3), 3),
                 format(round(sd_losses, 3), 3),
                 sep=' ± ')
    loss_df = as.data.frame(matrix(joined, nrow = length(fit_names)))
    rownames(loss_df) = fit_names
    colnames(loss_df) = names(results)
    formatted[[loss]] = loss_df
  }
  return(formatted)
}