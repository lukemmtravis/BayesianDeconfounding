rm(list=ls())
setwd('/Users/lmt15/Documents/phd/Causal Inference/Codes')
library(latex2exp)
source('Deconfounding.R')

n = 100
p = 200
s0 = 5
rho = 0.5

sample_sing_vals = function(n, p, s0, rho){
  X = make_design(n, p, rho=rho)
  svd_with_corr = svd(X)
  sing_vals_corr = svd_with_corr$d
  
  H = matrix(rnorm(n), nrow = n)
  Gamma = matrix(sqrt(1 - rho), nrow=1, ncol=p)
  E = matrix(rnorm(n*p, sd = sqrt(1-rho)), nrow=n, ncol=p)
  X = H %*% Gamma + E
  svd_proc = svd(X)
  sing_vals_proc = svd_proc$d
  return(list(corr = sing_vals_corr, proc = sing_vals_proc))
}



n_rep = 50
temp = pbreplicate(n_rep, sample_sing_vals(n, p, s0, rho))
mean_corr = apply(do.call(rbind, temp[1,]), 2, mean)
mean_proc = apply(do.call(rbind, temp[2,]), 2, mean)

df = data.frame(Sigma_rho=mean_corr, Confounding_Model=mean_proc, index = c(1:n)) %>% 
  gather(key='Process', value = 'SingVal', -index)
ggplot(df, aes(x = index, y = SingVal, color = Process)) + geom_point(alpha=0.75) +
  scale_color_discrete(labels=c('Confounding Model', 'Correlation Structure')) +
  theme(legend.position = 'top')
# ggsave('../Figures/corr_vs_confounding.pdf', units='in', width=6, height=6)


sample_correlation_effect = function(n, p, s0, rho){
  X = make_design(n, p, rho=rho)
  F_trim = trim_transform(X)
  F_lava = lava_transform(X)
  
  Fs = list(none = diag(n), trim = trim_transform(X), lava = lava_transform(X))
  X_tildes = lapply(Fs, function(Ftemp) Ftemp %*% X)
  corrs = lapply(X_tildes, function(Xt) (sum(cor(Xt)) - p)/(p**2 - p))
  return(corrs)
}

temp = t(pbreplicate(100, sample_correlation_effect(n, p, s0, rho)))
df = as.data.frame(t(apply(temp, 1, function(v) do.call(rbind, v))))
colnames(df) = c('Anone', 'Btrim', 'Clava')
df = df %>% gather(key='transform', value = 'pairwise_corr')

ggplot(df, aes(x = pairwise_corr, fill=transform, color=NULL)) +
  geom_histogram(position='identity', alpha = 0.5, bins=80) +
  geom_vline(xintercept=0.5, alpha = 0.2) + 
  geom_vline(xintercept=0, alpha=0.2) + 
  xlab(TeX(r'($ \hat{\rho} $)')) +
  scale_fill_discrete(labels=c('None', 'Trim', 'Lava')) +
  theme(legend.position='top')
# ggsave('../Figures/X_tilde_corr.pdf', units='in', width=6, height=6)









