#gen data

require(mvtnorm)
require(truncnorm)
require(tidyverse)


Nrespondents=5000
Qs=map(Nq_list,~profiles[sample(2^Nskill,.,replace=T),])
Xs=map(1:Nskill,~cbind(rep(1,Nrespondents),c(rep(0,Nrespondents/2),rep(1,Nrespondents/2))))

true_params=generate_params(Nrespondents,Qs,Xs)
true_alpha=sample_alpha_from_gamma(Nrespondents,true_params$gamma_list,Xs)
Ys=generate_data(Nrespondents,Nq_list,true_params$beta_mat,true_alpha)

sampler_out=sample_longitudinal(Ys,Xs,Qs,20,priors=list(beta_prior=.3,
                                                        gamma_prior=.1))

delta_cat=do.call(rbind,map(Qs,Q_to_delta))
plotdf=data.frame(postmean=sampler_out$beta_post$mean[delta_cat==1],
           postsd=sampler_out$beta_post$sd[delta_cat==1],
           true=true_params$beta_mat[delta_cat==1])%>%
  arrange(true)%>%
  mutate(i=1:n())
