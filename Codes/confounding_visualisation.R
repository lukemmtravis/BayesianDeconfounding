rm(list=ls())
setwd('/Users/lmt15/Documents/phd/Causal Inference/Codes')
source('Deconfounding.R')
source('NonLinearDeconfounding.R')
library(latex2exp)
library(ggpubr)

xs = matrix(seq(-1, 1, length.out = 5), nrow = 5)
data_fn = function(n, p, xs) make_data_LinConfounding(n, p, xs, q=2, function(x) abs(x)^{0.5})

k = function(x, y){return(mat_kern(x, y, l = 1))}
h = function(x, y){return( exp(-(x-y)**2))}

fits = list(
  gp = function(X, Y, xs) fit_gp(X, Y, xs, kernel=k),
  dgp = function(X, Y, xs) fit_dgp(X, Y, xs, kernel=k, dkernel=h, type='non-linear', sing_offset = 1e-6)
)
fn_0 = function(x) sqrt(sum(x**2))^{0.5}
matern32_kernel <- function(x1, x2, length_scale = 1.0, variance = 1.0) {
  # Compute the Euclidean distance between the two points
  distance <- sqrt(sum((x1 - x2)^2))
  
  # Compute the Matern 3/2 kernel
  factor <- sqrt(3) * distance / length_scale
  kernel_value <- variance * (1 + factor) * exp(-factor)
  
  return(kernel_value)
}
bm_kern = function(x1, x2, cn = 1){
  return(cn * min(x1, x2))
}

#### Linear Deconfounding
{
# Define kernel for GP prior on f  
k = function(x, y){return(mat_kern(x, y, l = 1))}
# k = function(x, y){return(matern32_kernel(x, y))}

# Make data with fixed seed
set.seed(1)
n = 500
p = 1
q = 2
rho = 0.0

# Make unconfounded components
E = make_design(n, p, rho=rho)
e = rnorm(n, sd=1)

# Make confounding components
H = matrix(rnorm(n*q), nrow = n, ncol = q)
Gamma = matrix(rnorm(q*p), nrow = q, ncol = p)
delta = rnorm(q)

# Compute b for information, not actually used in the method 
b = 1/(t(Gamma) %*% Gamma + 1) * t(Gamma) %*% delta
cat("n: ", n, "||b^0||_2: ", sum(b*b), '\n')

# Compute observed data
X = H %*% Gamma + E
Y = apply(X, 1, fn_0) + H %*% delta + e

# Define prior covariance on b
Sigma = diag(rep(1, p))

# Define test points
xs = matrix(seq(min(X), max(X), length.out = 200), nrow = 200)

# Define GP fits
gp = fit_gp(X, Y, xs, kernel=k)
dgp = fit_dgp(X, Y, xs, kernel=k)
fits = list(GP = gp, 'Deconfounded GP' = dgp)

# Format data
gp_dat = data.frame()
for(i in c(1:length(fits))){
  fit = fits[[i]]
  df_temp = data.frame(mean = fit$est, var = fit$pointwise_var, x = xs)
  df_temp['Method'] = names(fits)[i]
  gp_dat = rbind(gp_dat, df_temp)
}
train_dat = data.frame(x = X, y = Y, f0 = apply(X, 1, fn_0), Xb = X %*% b)
train_dat_1 = data.frame(x = X, y = Y, f0 = apply(X, 1, fn_0), 'f0Xb' = X %*% b + apply(X, 1, fn_0)) %>% 
  gather(key='Function', value='fX', -c(x, y))

# Define lims to make plot clearer (plot all GP information but not all
# training data information).
lims = c(min(gp_dat$mean - 1.96*sqrt(gp_dat$var)), max(gp_dat$mean + 1.96*sqrt(gp_dat$var))) 

ggplot(train_dat_1)  + 
  geom_point(aes(x = x, y = y), alpha=0.25) + 
  geom_line(aes(x = x, y = fX, linetype=Function)) +
  ggtitle('Data and fits with one feature and two linear confounding variables') +
  theme(legend.position='bottom') +
  xlab(TeX(r'($x$)')) +
  ylab(TeX(r'($y$)')) +
  scale_linetype_discrete(name='Function',
                       labels=c(TeX(r'($f_0(x)$)'), TeX(r'($f_0(x) + Xb$)'))) +
  ylim(lims)
# ggsave('../Figures/example_confounded_data.pdf', width=10, height=6, units='in')

ggplot(train_dat)  + 
  geom_point(aes(x = x, y = y), alpha=0.1) + 
  geom_line(aes(x=x, y = f0)) +
  geom_line(aes(x=x, y = f0)) +
  geom_line(aes(x=x, y = Xb + f0), linetype='dashed') +
  # geom_line(aes(x=x, y = Xb), linetype='dashed') +
  geom_line(data=gp_dat, aes(x = x, y = mean, color=Method)) +
  geom_ribbon(data=gp_dat, aes(x = x, ymin = mean - 1.96*sqrt(var), ymax  = mean + 1.96*sqrt(var), fill=Method), alpha=0.25) +
  ggtitle('Data and fits with one feature and two linear confounding variables') +
  theme(legend.position='bottom') +
  xlab(TeX(r'($x$)')) +
  ylab(TeX(r'($y$)')) +
  ylim(lims)


}
ggsave('../Figures/example_confounded_gp_regression_with_dash.pdf', width=10, height=6, units='in')

# Non Linear Deconfounding
{
fn_0 = function(x){
  return(abs(x)^{0.5})
}
gn_0 = function(x, omega=0.5){
  return(sin(2*pi*omega*x))
}  
k = function(x, y){return(mat_kern(x, y, l = 1))}
# h = function(x, y){return(matern32_kernel(x, y))}
h = function(x, y){return( exp(-(x-y)**2))}
# h = function(x, y){return( min(x, y) )}
set.seed(1)
n = 500
p = 1
q = 0
rho = 0.0
E = make_design(n, p, rho=rho)
e = rnorm(n, sd=1)
H = matrix(rnorm(n*q), nrow = n, ncol = q)

Gamma = matrix(rnorm(q*p), nrow = q, ncol = p)
delta = rnorm(q)  

# X = H %*% Gamma + E
X = matrix(runif(n, min = -2, max = 2), nrow=n)
Y = fn_0(X) + gn_0(X, omega=8) + e

# temp_dat = data_fn_4(n, p, xs)
# X = temp_dat$X
# Y = temp_dat$Y

xs = matrix(seq(min(X), max(X), length.out = 200), nrow = 200)
gp = fit_gp(X, Y, xs, kernel=k)
dgp = fit_dgp(X, Y, xs, kernel=k, dkernel=h, type='non-linear', sing_offset = 1e-6)
# fits = list(gp = gp, dgp_01 = dgp_01, dgp_100 = dgp_100)
fits = list(GP = gp, 'Deconfounded GP' = dgp)
gp_dat = data.frame()
for(i in c(1:length(fits))){
  fit = fits[[i]]
  df_temp = data.frame(mean = fit$est, var = fit$pointwise_var, x = xs)
  df_temp['Method'] = names(fits)[i]
  gp_dat = rbind(gp_dat, df_temp)
}
train_dat = data.frame(x = X, y = Y, f0 = fn_0(X))
ggplot(train_dat)  + 
  geom_point(aes(x = x, y = y), alpha=0.1) + 
  geom_line(aes(x=x, y = f0)) +
  geom_line(data=gp_dat, aes(x = x, y = mean, color=Method)) +
  geom_ribbon(data=gp_dat, aes(x = x, ymin = mean - 1.96*sqrt(var), ymax  = mean + 1.96*sqrt(var), fill=Method), alpha=0.1) +
  ggtitle('Data and fits with one feature and a non-linear perturbation') +
  theme(legend.position='bottom') +
  xlab(TeX(r'($x$)')) +
  ylab(TeX(r'($y$)'))
  # ylim(-1.5, 3.0)
}
# ggsave('../Writing/Figures/example_non_linear_perturbation.pdf', width=10, height=6, units='in')



sample_omega_data = function(omega){
  fn_0 = function(x){
    return(abs(x)^{0.5})
  }
  gn_0 = function(x){
    return(sin(2*pi*omega*x))
  }  
  k = function(x, y){return(mat_kern(x, y, l = 1))}
  # h = function(x, y){return(matern32_kernel(x, y))}
  h = function(x, y){return( exp(-(x-y)**2))}
  # h = function(x, y){return(mat_kern(x, y))}
  set.seed(1)
  n = 1000
  p = 1
  q = 0
  rho = 0.0
  E = make_design(n, p, rho=rho)
  e = rnorm(n, sd=1)
  H = matrix(rnorm(n*q), nrow = n, ncol = q)
  
  Gamma = matrix(rnorm(q*p), nrow = q, ncol = p)
  delta = rnorm(q)  
  
  # X = H %*% Gamma + E
  X = matrix(runif(n, min = -2, max = 2), nrow=n)
  Y = fn_0(X) + gn_0(X) + e
  
  #### Computing Lg_0
  H_X_X = apply(X, 1, function(x) apply(X, 1, function(y) h(x, y))) + diag(1e-6, n)
  L_squared = diag(rep(1, n)) - solve( diag(rep(1, n)) + H_X_X )
  svd_ = svd(L_squared)
  L = svd_$v %*% diag(sqrt(svd_$d)) %*% t(svd_$v)
  Lg_0 = sqrt(mean( (L %*% gn_0(X))**2 ))
  ####
  
  xs = matrix(seq(min(X), max(X), length.out = 200), nrow = 200)
  gp = fit_gp(X, Y, xs, kernel=k)
  dgp = fit_dgp(X, Y, xs, kernel=k, dkernel=h, type='non-linear', sing_offset=1e-6)
  # fits = list(gp = gp, dgp_01 = dgp_01, dgp_100 = dgp_100)
  fits = list(GP = gp, 'Deconfounded GP' = dgp)
  gp_dat = data.frame()
  for(i in c(1:length(fits))){
    fit = fits[[i]]
    df_temp = data.frame(mean = fit$est, var = fit$pointwise_var, x = xs)
    df_temp['Method'] = names(fits)[i]
    df_temp['Omega'] = omega
    df_temp['Lg_norm'] = Lg_0
    gp_dat = rbind(gp_dat, df_temp)
  }
  train_dat = data.frame(x = X, y = Y, f0 = fn_0(X), Omega=omega)
  return(list(gp_dat = gp_dat, train_dat = train_dat))
}

datasets = lapply(c(0.5,2), sample_omega_data)

train_dat = do.call(rbind, lapply(datasets, function(out) out$train_dat))
gp_dat = do.call(rbind, lapply(datasets, function(out) out$gp_dat))

gp_dat['Label'] = paste0('omega: ', gp_dat$Omega, ',   ||Lg_0||: ', round(gp_dat$Lg_norm, 3))

ggplot(train_dat)  + 
  # geom_point(aes(x = x, y = y), alpha=0.1) +
  geom_line(aes(x=x, y = f0)) +
  geom_line(data=gp_dat, aes(x = x, y = mean, color=Method)) +
  geom_ribbon(data=gp_dat, aes(x = x, ymin = mean - 1.96*sqrt(var), ymax  = mean + 1.96*sqrt(var), fill=Method), alpha=0.15) +
  facet_wrap(~Label) +
  # ggtitle('Data and fits with one feature and a non-linear perturbation') +
  theme(legend.position='bottom') +
  xlab(TeX(r'($x$)')) +
  ylab(TeX(r'($y$)'))
ggsave('../Figures/example_omegas_talk.pdf', width=12, height=8, units='in')
# ggsave('../Writing/Figures/example_omegas_no_title.pdf', width=10, height=6, units='in')

# sample_omega_data_bm = function(omega){
#   fn_0 = function(x){
#     return(abs(x - 0.5)^{0.5})
#   }
#   gn_0 = function(x){
#     return(sin(pi*omega*x))
#   }  
#   k = function(x, y){return(mat_kern(x, y, l = 1))}
#   h = function(x, y){return(bm_kern(x, y))}
#   
#   set.seed(1)
#   n = 1000
#   p = 1
# 
#   X = matrix(c(1:n)/(n+0.5), nrow=n)
#   Y = fn_0(X) + gn_0(X) + rnorm(n)
#   
#   #### Computing Lg_0
#   H_X_X = apply(X, 1, function(x) apply(X, 1, function(y) h(x, y)))
#   L_squared = diag(rep(1, n)) - solve( diag(rep(1, n)) + H_X_X )
#   svd_ = svd(L_squared)
#   L = svd_$v %*% diag(sqrt(svd_$d)) %*% t(svd_$v)
#   Lg_0 = sqrt(mean( (L %*% gn_0(X))**2 ))
#   ####
#   
#   xs = matrix(seq(min(X), max(X), length.out = 200), nrow = 200)
#   gp = fit_gp(X, Y, xs, kernel=k)
#   dgp = fit_dgp(X, Y, xs, kernel=k, dkernel=h, type='non-linear', sing_offset=1e-6)
#   # fits = list(gp = gp, dgp_01 = dgp_01, dgp_100 = dgp_100)
#   fits = list(GP = gp, 'Deconfounded GP' = dgp)
#   gp_dat = data.frame()
#   for(i in c(1:length(fits))){
#     fit = fits[[i]]
#     df_temp = data.frame(mean = fit$est, var = fit$pointwise_var, x = xs)
#     df_temp['Method'] = names(fits)[i]
#     df_temp['Omega'] = omega
#     df_temp['Lg_norm'] = Lg_0
#     gp_dat = rbind(gp_dat, df_temp)
#   }
#   train_dat = data.frame(x = X, y = Y, f0 = fn_0(X), Omega=omega)
#   return(list(gp_dat = gp_dat, train_dat = train_dat))
# }
# 
# datasets = lapply(c(0.5, 2, 250, 251), sample_omega_data_bm)
# 
# train_dat = do.call(rbind, lapply(datasets, function(out) out$train_dat))
# gp_dat = do.call(rbind, lapply(datasets, function(out) out$gp_dat))
# 
# gp_dat['Label'] = paste0('omega: ', gp_dat$Omega, ',   ||Lg_0||: ', round(gp_dat$Lg_norm, 3))
# 
# ggplot(train_dat)  + 
#   # geom_point(aes(x = x, y = y), alpha=0.1) + 
#   geom_line(aes(x=x, y = f0)) +
#   geom_line(data=gp_dat, aes(x = x, y = mean, color=Method)) +
#   geom_ribbon(data=gp_dat, aes(x = x, ymin = mean - 1.96*sqrt(var), ymax  = mean + 1.96*sqrt(var), fill=Method), alpha=0.15) +
#   facet_wrap(~Label) +
#   # ggtitle('Data and fits with one feature and a non-linear perturbation') +
#   theme(legend.position='bottom') +
#   xlab(TeX(r'($x$)')) +
#   ylab(TeX(r'($y$)'))



#### Two Dimensional Linear Deconfounding
{
# Define kernel for GP prior on f  
k = function(x, y){return(mat_kern(x, y, l = 1))}
# k = function(x, y){return(matern32_kernel(x, y))}
fn_0 = function(x) sqrt(sum(x**2))^{0.5}
  
# Make data with fixed seed
set.seed(1)
n = 1000
p = 2
q = 2
rho = 0.0

# Make unconfounded components
E = make_design(n, p, rho=rho)
# E = matrix(runif(n, min = -2, max = 2), nrow = n)
e = rnorm(n, sd=1)

# Make confounding components
H = matrix(rnorm(n*q), nrow = n, ncol = q)
Gamma = matrix(rnorm(q*p), nrow = q, ncol = p)
delta = rnorm(q)  

# Compute b for information, not actually used in the method 
b = solve(t(Gamma) %*% Gamma + diag(p)) %*% t(Gamma) %*% delta

# Compute observed data
X = H %*% Gamma + E
Y = apply(X, 1, fn_0) + H %*% delta + e

# Define prior covariance on b
Sigma = diag(rep(1, p))

r = 2
thetas = seq(0, 2*pi, length.out = 100)
xs = cbind(r*cos(thetas), r*sin(thetas))

# Define GP fits
gp = fit_gp(X, Y, xs, kernel=k)

K_X_X = apply(X, 1, function(x) apply(X, 1, function(y) k(x, y)))

# dgp_01 = fit_dgp(X, Y, xs, kernel=k, dSigma=0.01*diag(rep(1, p)))
# dgp_100 = fit_dgp(X, Y, xs, kernel=k, dSigma=100*diag(rep(1, p)))
dgp = fit_dgp(X, Y, xs, kernel=k, dSigma=diag(rep(1, p)))
# fits = list(gp = gp, dgp_01 = dgp_01, dgp_1 = dgp, dgp_100 = dgp_100)
fits = list(GP = gp, 'Deconfounded GP' = dgp)

# Format data
gp_dat = data.frame()
for(i in c(1:length(fits))){
  fit = fits[[i]]
  df_temp = data.frame(mean = fit$est, var = fit$pointwise_var, theta = thetas)
  df_temp['Method'] = names(fits)[i]
  gp_dat = rbind(gp_dat, df_temp)
}
train_dat = data.frame(theta = thetas, f0 = apply(xs, 1, fn_0), 'Plus Perturbation' = apply(xs, 1, fn_0) + xs %*% b)
train_dat = train_dat %>% gather(key='True Stat.', value='Val', -theta)

# Define lims to make plot clearer (plot all GP information but not all
# training data information).
lims = c(min(gp_dat$mean - 1.96*sqrt(gp_dat$var)), max(gp_dat$mean + 1.96*sqrt(gp_dat$var))) 

ggplot(train_dat)  + geom_line(aes(x=theta, y=Val, linetype=`True Stat.`)) +
  # geom_line(aes(x=theta, y = f0)) + x
  # geom_line(aes(x=theta))
  # geom_line(aes(x=x, y = Xb), linetype='dashed') +
  geom_line(data=gp_dat, aes(x = theta, y = mean, color=Method)) +
  geom_ribbon(data=gp_dat, aes(x = theta, ymin = mean - 1.96*sqrt(var), ymax  = mean + 1.96*sqrt(var), fill=Method), alpha=0.1) +
  ggtitle('f_0(X(theta)), f_0(X(theta)) + X(theta) b and fits with two feature and two linear confounding variables') +
  theme(legend.position='bottom') +
  xlab(TeX(r'($ \theta$)')) +
  ylab(TeX(r'($y$)')) +
  ylim(lims)
}
ggsave('../Figures/2d_embedding_example.pdf', width=10, height=6, units='in')

# Non Linear Confounded Model
{
  fn_0 = function(x){
    return(abs(x)^{0.5})
  }
  phi = function(h, omega=0.5){
    r = sqrt(sum(h**2))
    return(2*sin(2*pi*omega*r))
  }  
  k = function(x, y){return(mat_kern(x, y, l = 1))}
  h = function(x, y){return( exp(-(x-y)**2))}
  set.seed(1)
  n = 2000
  p = 1
  q = 2
  rho = 0.0
  omega = 8
  
  Gamma = matrix(rnorm(q*p), nrow = q, ncol = p)
  e = rnorm(n, sd=1)
  
  # To randomly sample H then X | H
  
  # GAUSSIAN
  # H = matrix(rnorm(n*q), nrow = n, ncol = q)
  # E = make_design(n, p, rho=rho)
  
  # UNIFORM
  H = matrix(runif(n*q, -2, 2), nrow = n, ncol = q)
  # E = matrix(runif(n*p, -2, 2), nrow = n, ncol = p)
  E = matrix(rnorm(n*p), nrow = n, ncol = p)
  X = H %*% Gamma + E
  
  
  # To sample X uniformly, then generate L_2-approximating H.
  # E = make_design(n, p, rho=rho)
  # X = matrix(runif(n, min=-2, max = 2), nrow=n)
  # H = X %*% t(Gamma) %*% solve(Gamma %*% t(Gamma))
  
  Y = fn_0(X) + apply(H, 1, function(h) phi(h, omega=omega)) + e
  
  xs = matrix(seq(min(X), max(X), length.out = 200), nrow = 200)
  gp = fit_gp(X, Y, xs, kernel=k)
  dgp = fit_dgp(X, Y, xs, kernel=k, dkernel=h, type='non-linear', sing_offset = 1e-6)
  ldgp = fit_dgp(X, Y, xs, kernel=k, dSigma='trim', type='linear', sing_offset = 1e-6)
  # fits = list(gp = gp, dgp_01 = dgp_01, dgp_100 = dgp_100)
  fits = list(GP = gp, 'Non-linearly Deconfounded GP' = dgp, 'Linearly Deconfounded GP' = ldgp)
  gp_dat = data.frame()
  for(i in c(1:length(fits))){
    fit = fits[[i]]
    df_temp = data.frame(mean = fit$est, var = fit$pointwise_var, x = xs)
    df_temp['Method'] = names(fits)[i]
    gp_dat = rbind(gp_dat, df_temp)
  }
  train_dat = data.frame(x = X, y = Y, f0 = fn_0(X))
  
  H_dat = data.frame(X = X, HGamma = H %*% Gamma, g = apply(H, 1, function(h) phi(h, omega=omega)))
  
  plot1 = ggplot(train_dat)  + 
    # geom_point(aes(x = x, y = y), alpha=0.1) +
    geom_line(aes(x=x, y = f0)) +
    geom_line(data=gp_dat, aes(x = x, y = mean, color=Method)) +
    facet_wrap(~Method) + 
    geom_ribbon(data=gp_dat, aes(x = x, ymin = mean - 1.96*sqrt(var), ymax  = mean + 1.96*sqrt(var), fill=Method), alpha=0.1) +
    # geom_point(data = H_dat, aes(x = X, y = g), color='orange', alpha=0.5) +
    ggtitle('Non-linear perturbation (omega=8) of ||H|| (q = 2) with non-linear regression in X (p = 1)') +
    theme(legend.position='bottom') +
    xlab(TeX(r'($x$)')) +
    ylab(TeX(r'($y$)'))
  # ylim(-1.5, 3.0)
  
  plot2 = ggplot(train_dat)  + 
    # geom_point(aes(x = x, y = y), alpha=0.1) + 
    geom_line(aes(x=x, y = f0)) +
    geom_line(data=gp_dat, aes(x = x, y = mean, color=Method)) +
    geom_ribbon(data=gp_dat, aes(x = x, ymin = mean - 1.96*sqrt(var), ymax  = mean + 1.96*sqrt(var), fill=Method), alpha=0.1) +
    # geom_point(data = H_dat, aes(x = X, y = g), color='orange', alpha=0.5) +
    # ggtitle('Non-linear perturbation (omega=8) of ||H|| (q = 2) with non-linear regression in X (p = 1)') +
    theme(legend.position='bottom') +
    xlab(TeX(r'($x$)')) +
    ylab(TeX(r'($y$)'))
  # ylim(-1.5, 3.0)
  
  ggarrange(plot1, plot2, common.legend = TRUE, legend='bottom', nrow=2)
}
# ggsave('../Figures/example_non_linear_perturbation.pdf', width=10, height=6, units='in')


# Looking at different
plot_fits = function(X, Y, fn_0, title){
  n = dim(X)[1]
  p = dim(X)[2]
  xs = matrix(seq(min(X), max(X), length.out = 200), nrow = 200)
  # k = function(x, y) se_kern(x, y, l = 1)
  k_r = function(x, y) se_kern(x, y, l = 2^{1/4} * n^{-1/4})
  h = function(x, y){return( exp(-(x-y)**2))}
  
  gp = fit_gp(X, Y, xs, kernel=k_r)
  dgp = fit_dgp(X, Y, xs, kernel=k_r, dkernel=h, type='non-linear', sing_offset = 1e-6)
  ldgp = fit_dgp(X, Y, xs, kernel=k_r, dSigma='trim', type='linear', sing_offset = 1e-6)
  dgpDouble = fit_dgpDouble(X, Y, xs, kernel=k_r, dSigma ='trim', dkernel=h, sing_offset = 1e-6)
  fits = list(A = gp, B = ldgp, 
              C = dgp, D = dgpDouble)
  gp_dat = data.frame()
  for(i in c(1:length(fits))){
    fit = fits[[i]]
    df_temp = data.frame(mean = fit$est, var = fit$pointwise_var, x = xs)
    df_temp['Method'] = names(fits)[i]
    gp_dat = rbind(gp_dat, df_temp)
  }
  method.labs <- c("GP", "Deconfounded GP (Linear)", "Deconfounded GP (Non-linear)", "Deconfounded GP (Double)")
  names(method.labs) <- c("A", "B", "C", "D")
  train_dat = data.frame(x = X, y = Y, f0 = fn_0(X))
  plot1 = ggplot(train_dat)  + 
    geom_line(aes(x=x, y = f0)) +
    geom_line(data=gp_dat, aes(x = x, y = mean, color=Method)) +
    facet_wrap(~Method, nrow = 1, labeller = labeller(Method = method.labs)) + 
    geom_ribbon(data=gp_dat, aes(x = x, ymin = mean - 1.96*sqrt(var), ymax  = mean + 1.96*sqrt(var), fill=Method), alpha=0.1) +
    ggtitle(title) +
    theme(legend.position='bottom') +
    scale_color_discrete(name='Method', 
                         labels=c('GP', 'Deconfounded GP (Linear)', 
                                  'Deconfounded GP (Non-linear)', 'Deconfounded GP (Double)')) +
    scale_fill_discrete(name='Method', 
                         labels=c('GP', 'Deconfounded GP (Linear)', 
                                  'Deconfounded GP (Non-linear)', 'Deconfounded GP (Double)')) +
    xlab(TeX(r'($x$)')) +
    ylab(TeX(r'($f(x)$)'))
  return(plot1)
}
fn_0 = function(x){
  # return(abs(x)^{0.5} + abs(x + 1)^{0.75} - abs(x-1)^{0.75})
  return(abs(x-0.5)^{0.5})
}
phi = function(h, omega=0.5){
  r = sqrt(sum(h**2))
  return(2*sin(2*pi*omega*r))
}
set.seed(1)
n = 2000
p = 1
q = 2
rho = 0.0
omega = 4

Gamma = matrix(rnorm(q*p), nrow = q, ncol = p)
e = rnorm(n, sd=1)
H = matrix(runif(n*q, -2, 2), nrow = n, ncol = q)
E = matrix(rnorm(n*p), nrow = n, ncol = p)
X = H %*% Gamma + E

# No confounding
Y = fn_0(X) + e
no_confounding = plot_fits(X, Y, fn_0, title=TeX(r'($ Y = f(X) + \epsilon $)'))

# Linear confounding
delta = rnorm(q)
Y = fn_0(X) + H %*% delta + e
lin_confounding = plot_fits(X, Y, fn_0, title=TeX(r'($ Y = f(X) + H \delta + \epsilon $)'))

# Non-linear confounding
Y = fn_0(X) + apply(H, 1, function(h) phi(h, omega=omega)) + e
nonlin_confounding = plot_fits(X, Y, fn_0, title= TeX(r'($ Y = f(X) + \phi(H) + \epsilon $)'))

# Both forms of confounding
Y = fn_0(X) + H %*% delta + apply(H, 1, function(h) phi(h, omega=omega)) + e
double_confounding = plot_fits(X, Y, fn_0, title= TeX(r'($ Y = f(X) + H \delta + \phi(H) + \epsilon $)'))

ggarrange(no_confounding,
          lin_confounding,
          nonlin_confounding, 
          double_confounding,
          common.legend = TRUE, legend='bottom', nrow=4)
# ggsave('../Figures/ConfoundingSettingsSE.pdf', units='in', width=8, height=10)


### NOTE: IT SEEMS TO BE THE CASE THAT THE SE KERN IS BETTER FOR `RIGHT SMOOTHNESS`
### BUT MATERN KERNEL IS BETTER FOR MODELLING THE FUNCTION BECAUSE OF THE KINKS IN THE CURVE


compute_coords = function(n, p, q, s0){
  n = 100
  p = 200
  q = 5
  s0 = 5
  signal_size = log(n)
  noise_var = 1
  beta_0 = make_true_param(p, s0, signal_size=signal_size, include_first_param=TRUE)
  E = make_design(n, p, rho=0, AR=FALSE)
  e = rnorm(n, sd=sqrt(noise_var))
  H = matrix(rnorm(n*q), nrow = n, ncol = q)
  Gamma = matrix(rnorm(q*p), nrow = q, ncol = p)
  delta = rep(log(n), q)
  b = solve(t(Gamma)%*%Gamma + diag(1, p))%*%t(Gamma) %*% delta
  X = H %*% Gamma + E
  Y = X %*% beta_0 + H %*% delta + e
  svd_X = svd(X)
  V = svd_X$v
  d = svd_X$d
  
  b_coeffs = (t(V) %*% b)
  beta_coeffs = (t(V) %*% beta_0)
  d = svd_X$d
  d_trim = pmin(d, quantile(d, 0.5))
  mu = quantile(d, 0.5)/2
  d_lava = sqrt( n*mu*d**2/(n*mu + d**2) )

  return(list(b = abs(b_coeffs), beta = abs(beta_coeffs), 
              d=d, d_trim = d_trim, d_lava=d_lava
              ))
}

num_sims = 5000
results = replicate(num_sims, compute_coords(100, 200, 5, 5))
b = apply(do.call(cbind, results[1,]), 1, mean)
beta = apply(do.call(cbind, results[2,]), 1, mean)

d = apply(do.call(cbind, results[3,]), 1, mean)
d_trim = apply(do.call(cbind, results[4,]), 1, mean)
d_lava = apply(do.call(cbind, results[5,]), 1, mean)


coeffs_plot_dat = data.frame(index = c(1:100), b = b, beta = beta) %>% 
  gather(key = 'a', value = 'CoordinateSize', -index)
coeffs_plot = ggplot(coeffs_plot_dat, aes(x = index, y = CoordinateSize, color = a)) + geom_point() +
  geom_line() +
  scale_color_discrete(name=TeX(r'($\bf{a}$)'), 
                       labels=c(TeX(r'($b$)'), TeX(r'($\beta$)'))) +
  xlab(TeX(r'($k$)')) +
  ylab(TeX(r'($ \bf{v}_k^T \bf{a} $)')) +
  theme(legend.position = 'bottom')
# ggsave('../Figures/Coefficients.pdf', units='in', width=8, height=5)

plot_dat = data.frame(index = c(1:100), d = d, d_trim = d_trim, d_lava = d_lava) %>% 
  gather(key = 'Transform', value = 'SingValue', -index)
singvals_plot = ggplot(plot_dat, aes(x = index, y = SingValue, color = Transform)) + 
  geom_point() +
  geom_line() +
  scale_color_discrete(name='Transform',
                       labels=c('None', 'Lava', 'Trim')) +
  xlab(TeX(r'($k$)')) +
  ylab(TeX(r'($ d_k $)')) +
  theme(legend.position = 'bottom')

ggarrange(singvals_plot, coeffs_plot, legend='bottom', nrow=1)

ggsave('../Figures/SingValues.pdf', units='in', width=12, height=6)





bound = function(n, j, alpha, beta, gamma){
  if(abs(j) < n^{1/(1+2*alpha)}){
    return( (n^{2/(1+2*alpha)})^{-1/4 - gamma} )
  }else{
    return( (n^{2/(1+2*alpha)})^{-1/4 - gamma} * abs(j)^{beta} )
  }
}


bound_spectrum = function(n, alpha, beta, gamma){
  js = seq(1, n, by = 1)
  bounds = sapply(js, function(j) bound(n, j, alpha, beta , gamma))
  return(data.frame(j = js, bound = bounds, 
                    n = rep(n, length(js)), alpha = rep(alpha, length(js)), gamma = rep(gamma, length(js))))
}


n = 1000
alpha = 1
gamma = 1
beta = 2
test = bound_spectrum(n, alpha, beta, gamma)

ggplot(test, aes(x = j, y = bound)) + geom_line()

