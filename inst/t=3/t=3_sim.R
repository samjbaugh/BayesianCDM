#gen data
rm(list = ls())
set.seed(300)
setwd("~/Desktop/CDM/")
devtools::load_all()

require(mvtnorm)
require(truncnorm)
require(tidyverse)


M            = 3000
Nrespondents = 800
Nskill       = 3
Ntime        = 3
k            = 3
Nprofile     = 2^Nskill


profiles=map_dfr(1:Nprofile-1,
                 ~data.frame(t(rev(as.integer(intToBits(.))[1:Nskill]))))
profile_list = as.list(data.frame(t(profiles)))
# Alternative for "profiles"
# profile_list=list()
# for(ii in 1:(2^Nskill)){
#   tmp=rep(NA,Nskill)
#   for(jj in 1:Nskill){
#     tmp[Nskill-jj+1]=((ii-1)%/%(2^(jj-1)))%%2
#   }
#   profile_list[[ii]]=rev(tmp)
# }



# Q Matrix
Q=do.call(rbind,lapply(profile_list[-1],function(x)
  do.call(rbind,lapply(1:ceiling(16/(Nprofile-1)), function(i) x))))

# Random Q
# Q=map(Nq_list,~profiles[sample(2^Nskill,Nq,replace=T),]) 

# Manual Q
# q_rep = 3
# Q = matrix(c(rep(c(0,0,0),q_rep),
#              rep(c(0,0,1),q_rep),
#              rep(c(0,1,0),q_rep),
#              rep(c(1,0,0),q_rep),
#              rep(c(1,1,0),q_rep),
#              rep(c(1,0,1),q_rep),
#              rep(c(0,1,1),q_rep),
#              rep(c(1,1,1),q_rep)), q_rep*8, 3, byrow=T)
Qs      = list(Q,Q,Q)
Nq_list = map(Qs,~dim(.)[1])
Nq      = dim(Q)[1]
Nq_total=sum(unlist(Nq_list))



# Standard one-intervention covariate
Xbase=cbind(rep(1,Nrespondents),c(rep(0,Nrespondents/2),rep(1,Nrespondents/2)))
Xs=map(1:Nskill,list)
for(i in 1:Nskill){
  Xs[[i]][[2]]=Xbase
  Xs[[i]][[3]]=Xbase
  Xs[[i]][[4]]=Xbase
  Xs[[i]][[5]]=Xbase
  for(j in (1:(2^Ntime-1))[-c(2:5)]){
    Xs[[i]][[j]]=matrix(1,dim(Xbase)[1],1)
  }
}
# > T=3 transitions
# X1 X2 X3
# 1  0  0  0
# 2  0  0  1
# 3  0  1  0
# 4  0  1  1
# 5  1  0  0
# 6  1  0  1
# 7  1  1  0
# 8  1  1  1

# Complex covariates 
# Intervention = c(rep(0,Nrespondents/2),rep(1,Nrespondents/2))
# Education = sample(0:3, Nrespondents, replace=T)
# Income = rnorm(Nrespondents, 0, 5)
# Ethnicity = sample(0:5, Nrespondents, replace=T)
# Xbase=cbind(rep(1,Nrespondents),
#             Intervention,
#             Education,
#             Income,
#             Ethnicity)
# Xs=build_Xlist_01(Xbase,Nskill,Ntime)


delta=Q_to_delta(Q)
true_params=list()
true_params$beta_mat=
  cbind(rtruncnorm(Nq,mean=-1.5,sd=1,b=0),rtruncnorm(Nq,mean=3,sd=1,a=0),
        rtruncnorm(Nq,mean=2,sd=1,a=0),rtruncnorm(Nq,mean=2,sd=1,a=0), # intercept + base
        rtruncnorm(Nq,mean=0.5,sd=0.5,a=0),rtruncnorm(Nq,mean=0.5,sd=0.5,a=0),
        rtruncnorm(Nq,mean=0.5,sd=0.5,a=0),rtruncnorm(Nq,mean=0.5,sd=0.5,a=0))*   # interactions
  delta
true_params$gamma_list=map(Xs,~map(.,function(x) rnorm(dim(x)[2],sd=0.3)))

gamma_to_transprobs(true_params$gamma_list, Xs) 
plot(true_params$beta_mat[delta==1])



# a=c(.7,.05,.05,.2)
# b=c(.2,.55,.05,.2)
# g1=log(a[2:4]/(1-sum(a[2:4])))
# g2=log(b[2:4]/(1-sum(b[2:4])))-g1
# gamma_list_true=map(1:Nskill,~
#   rbind(g1,g2),
# )
# Xtmp=build_Xlist_01(rbind(c(1,0),c(1,1)),Nskill,Ntime)
# gamma_list_true=
#   list(list(c(-2,3),c(-1.5),c(-1.5)),
#        list(c(-2,3),c(-1.5),c(-1.5)),
#        list(c(-2,3),c(-1.5),c(-1.5)))
# gamma_list=gamma_list_true
# true_params$gamma_list=gamma_list
# 
# gamma_to_transprobs(gamma_list,Xtmp)

set.seed(1)
registerDoParallel(10)
Nsims=300

sim_results = foreach(sim = 1:Nsims) %dorng% {
  

group_is=Xbase[,2]+1
alpha_sampdf=sample_alpha_from_gamma(Nrespondents,Ntime,true_params$gamma_list,Xs)%>%
  data.frame()%>%
  set_names(paste0('a',1:Ntime))%>%
  mutate(group=group_is)%>%
  relocate(group,.before=a1)%>%
  group_by(group)%>%
  arrange(group,a1,a2)
true_alpha=alpha_sampdf%>%
  as.matrix()%>%
  .[,-1]

# simulate data
profiles=map_dfr(1:Nprofile-1,
                 ~data.frame(t(rev(as.integer(intToBits(.))[1:Nskill]))))
A=model.matrix(data=profiles,as.formula(paste0('~', paste0(names(profiles),collapse='*'))))

ltheta=true_params$beta_mat%*%t(A)
theta=logistic(ltheta)
Nq_total=sum(unlist(Nq_list))
qt_map=unlist(map(1:length(Nq_list),~rep(.,Nq_list[[.]])))

Ys=map(Nq_list,~matrix(NA,Nrespondents,.))
for (t in 1:Ntime){
  for(i in 1:Nrespondents){
    for(j in 1:Nq){
      p=theta[j,true_alpha[i,t]]
      Ys[[t]][i,j]=sample(c(0,1),1,prob=c(1-p,p))
    }
  }
}


priors=list(beta_prior=2.5,gamma_prior=1)
runtime = system.time({
sampler_out=fit_longitudinal_cdm_full(Ys,Xs,Qs,M,
                                      initparams = NULL,
                                      priors=priors,
                                      fixed_beta = TRUE)})

# rm burn-in
sampler_out$samples_burned_in = sampler_out$samples[500:M,]

list(sim = sim, true_alpha = true_alpha, Ys = Ys, 
     sampler_out = sampler_out)

}

#################################################################################################################################################################

true_alphas = lapply(1:Nsims, function(x) sim_results[[x]]$true_alpha)
Ys = lapply(1:Nsims, function(x) sim_results[[x]]$Ys)

# Beta cleaning
beta_mat=t(sapply(1:Nsims, function(x) sim_results[[x]]$sampler_out$samples_burned_in%>%
                    dplyr::select(contains('beta'))%>%
                    apply(2,mean)))
beta_sd=t(sapply(1:Nsims, function(x) sim_results[[x]]$sampler_out$samples_burned_in%>%
                   dplyr::select(contains('beta'))%>%
                   apply(2,sd)))
# Beta credible interval coverage
beta_upper = beta_mat + 2*beta_sd
beta_lower = beta_mat - 2*beta_sd
true_beta = true_params$beta_mat[delta==1]
beta_credible_post = matrix(0,Nsims,length(true_beta))
for(i in 1:length(true_beta))
{
  beta_credible_post[,i] = (beta_lower[,i] < true_beta[i]) & (true_beta[i] < beta_upper[,i])
}
barplot(round(apply(beta_credible_post,2,mean),3), main = "Beta Coverage",
        names = 1:length(true_beta), col = c(rep(7, 21), rep(2,36), rep(3,21)))

# Gamma cleaning
gamma_post=t(sapply(1:Nsims, function(x) sim_results[[x]]$sampler_out$samples_burned_in%>%
                      dplyr::select(contains('gamma'))%>%
                      apply(2,mean)))
gamma_sd=t(sapply(1:Nsims, function(x) sim_results[[x]]$sampler_out$samples_burned_in%>%
                    dplyr::select(contains('gamma'))%>%
                    apply(2,sd)))
# Gamma credible interval coverage
gamma_upper = gamma_post + 2*gamma_sd
gamma_lower = gamma_post - 2*gamma_sd
true_gamma = unlist(true_params$gamma_list)
gamma_credible_post = matrix(0,Nsims,length(true_gamma))
for(i in 1:length(true_gamma))
{
  gamma_credible_post[,i] = (gamma_lower[,i] < true_gamma[i]) & (true_gamma[i] < gamma_upper[,i])
}
barplot(round(apply(gamma_credible_post,2,mean),3), ylim = c(0,1), main = "Gamma Coverage",
        names = 1:length(true_gamma), col = c(rep(11, length(true_gamma)/3),
                                              rep(12, length(true_gamma)/3),
                                              rep(13, length(true_gamma)/3)))

### Beta posteriors
ggplot(data.frame(beta.post=c(beta_mat), q=rep(1:78, each = 300), 
                  beta.true=rep(true_params$beta_mat[delta==1], each=300))) + 
  geom_point(aes(x=q ,y=beta.true, col='betatrue'))+  
  geom_boxplot(aes(y=beta.post, x=q, group = q, col='beta'), outlier.size=0,
               fill = "white", position="identity", alpha=.1)+
  theme(legend.position = c(.9, .15),
        legend.background = element_rect(fill="white",
                                         linewidth=0.5, linetype="solid", 
                                         colour ="black"))+
  labs(color=NULL)

### Gamma posteriors
plotdf=data.frame(gamma_post=c(gamma_post),
                  gamma_true=rep(unlist(true_params$gamma_list),each=300),
                  i=1:length(c(gamma_post)))%>%
  mutate(skill=c(rep(1,11*300),rep(2,11*300),rep(3,11*300)),i=rep(c(1:11,1:11,1:11),each=300),
         q=rep(1:33,each=300))
pgamma=ggplot(plotdf)+
  # geom_errorbar(aes(x=i-.25,ymin=postmean-postsd*2,ymax=postmean+postsd*2,col='posterior mean'))+
  geom_boxplot(aes(y=gamma_post, x=i, group=q, col='posterior mean'), outlier.size=0,
               fill = "white", position="identity", alpha=.5)+  
  geom_point(aes(x=i,y=gamma_true,col='true'))+
  # geom_point(aes(x=i,y=gamma_true,col='true'))+
  # geom_point(aes(x=i+.25,y=pmean,col='posterior mean'))+
  # geom_errorbar(aes(x=i+.25,ymin=pmean-psd*2,ymax=pmean+psd*2,col='ss'))+
  facet_grid(~skill) +
  theme(legend.position = c(.9, .92),
        legend.background = element_rect(fill="white",
                                         linewidth=0.5, linetype="solid", 
                                         colour ="black"))+
  labs(color=NULL)
pgamma

{
  # delta=do.call(rbind,map(Qs,Q_to_delta))
  beta_mat=sampler_out$samples_burned_in%>%
    dplyr::select(contains('beta'))%>%
    apply(2,mean)
  beta_sd=sampler_out$samples_burned_in%>%
    dplyr::select(contains('beta'))%>%
    apply(2,sd)
  pbeta=data.frame(postmean=beta_mat,
                   betatrue=true_params$beta_mat[delta==1],
                   betasd=beta_sd)%>%
    # filter(postmean<8)%>%
    mutate(i=1:n())%>%
    ggplot()+
    geom_point(aes(x=i,y=postmean,col='postmean'))+
    geom_errorbar(aes(x=i,ymin=postmean-2*betasd,
                      ymax=postmean+2*betasd,col='postmean'))+
    geom_point(aes(x=i,y=betatrue,col='betatrue'))  
  # coord_cartesian(ylim = c(-5, 10))
  pbeta
}


plot(sampler_out$samples_burned_in[,32], type = "l")


{
  sampler_out = sim_results[[1]]$sampler_out
  true_alpha = true_alphas[[1]]
  alpha_post=sampler_out$samples_burned_in%>%
    dplyr::select(starts_with('alpha'))
  # M=20
  # data.frame(a=as.numeric(alpha_post[M,]),ii=1:(Nrespondents*Ntime))%>%
  #   mutate(i=ti_map[ii],t=it_map[ii])%>%
  #   mutate(true_alpha=c(true_alpha))%>%
  #   ggplot()+geom_point(aes(x=i,y=a,col=factor(true_alpha)))+
  #   facet_grid(~t)
  
  profiles=map_dfr(1:Nprofile-1,
                   ~data.frame(t(rev(as.integer(intToBits(.))[1:Nskill]))))
  M=dim(alpha_post)[1]
  it_map=unlist(map(1:Ntime,~rep(.,Nrespondents)))
  ti_map=unlist(map(1:Ntime,~1:Nrespondents))
  post_probs=alpha_post%>%
    pivot_longer(everything())%>%
    mutate(name=factor(name,levels=paste0('alpha_vec[',1:dim(alpha_post)[2],']')),
           value=factor(value,levels=1:Nprofile,labels=paste0('prof',1:Nprofile)))%>%
    group_by(name,value)%>%
    summarise(prop=n()/M)%>%
    complete(value,fill=list(prop=0))%>%mutate(ii=as.numeric(map(strsplit(as.character(name),'\\[|\\]'),2)),
                                               i=ti_map[ii],t=it_map[ii])%>%
    rowwise()%>%
    mutate(true_prof=true_alpha[i,t])
  true_probs=post_probs%>%
    group_by(name,true_prof,i,t)%>%
    summarise(prop_correct=prop[true_prof[1]])
  prof_labels=
    apply(profiles,1,function(x) paste(x,collapse=''))
  palpha=true_probs%>%
    mutate(t=factor(t,levels=1:Ntime,labels=paste0('time',1:Ntime)))%>%
    ggplot()+
    geom_point(aes(x=i,y=prop_correct,color=factor(true_prof)))+
    facet_grid(~t)+
    xlab('Person ID')+
    scale_color_discrete(name='profile',labels=prof_labels)
  palpha
}
# scale_y_log10()


{
  gamma_postmean=sampler_out$samples_burned_in%>%
    dplyr::select(contains('gamma'))%>%
    apply(2,mean)
  gamma_postsd=sampler_out$samples_burned_in%>%
    dplyr::select(contains('gamma'))%>%
    apply(2,sd)
  
  trans_mat=alpha_to_transitions(true_alpha,profiles)
  priorsd_gamma=map(Xs,~map(.,function(x) priors$gamma_prior*diag(dim(x)[2])))
  tmp=sample_gamma(true_params$gamma_list,trans_mat,Xs,priorsd_gamma,retmean=T)
  plotdf=data.frame(gamma_true=unlist(true_params$gamma_list),
                    postmean=gamma_postmean,
                    postsd=gamma_postsd,
                    pmean=unlist(tmp$gamma_mean),
                    psd=sqrt(unlist(map(tmp$gamma_var,~map(.,diag)))),
                    i=1:length(gamma_postmean))%>%
    mutate(skill=c(rep(1,11),rep(2,11),rep(3,11)),i=c(1:11,1:11,1:11))
  pgamma=ggplot(plotdf)+
    geom_point(aes(x=i-.25,y=postmean,col='posterior mean'))+
    geom_errorbar(aes(x=i-.25,ymin=postmean-postsd*2,ymax=postmean+postsd*2,col='posterior mean'))+
    geom_point(aes(x=i,y=gamma_true,col='true'))+
    # geom_point(aes(x=i+.25,y=pmean,col='posterior mean'))+
    # geom_errorbar(aes(x=i+.25,ymin=pmean-psd*2,ymax=pmean+psd*2,col='ss'))+
    facet_grid(~skill)
  pgamma
}

pbeta
palpha
pgamma
ggsave(pbeta,file='../Figures/beta_recover.png',height=5,width=10)
ggsave(palpha,file='../Figures/alpha_recover.png',height=5,width=10)
ggsave(pgamma,file='../Figures/gamma_recover.png',height=5,width=10)




if(F){
  A=model.matrix(data=profiles,
                 as.formula(paste0('~', paste0(names(profiles),collapse='*'))))
  
  gamma_list=list(gamma_list_true[[1]][1,,drop=F],
                  gamma_list_true[[2]][1,,drop=F])
  Xs=list(matrix(rep(1,Nrespondents),Nrespondents,1),
          matrix(rep(1,Nrespondents),Nrespondents,1))
  gamma_to_probs(gamma_list,Xtmp)
  true_alpha=sample_alpha_from_gamma(Nrespondents,Ntime,gamma_list,Xs)
  true_trans=alpha_to_transitions(true_alpha,profiles)
  trans_mat=alpha_to_transitions(true_alpha,profiles)
  
  eta=map2(gamma_list,Xs,~(.y%*%.x)%>%
             {cbind(rep(0,Nrespondents),.)}%>%
             apply(1,function(x) x-log(sum(exp(x))-exp(x)))%>%
             t())
  tstar=map(1:Nskill,~matrix(NA,Nrespondents,2^Ntime))
  for(i in 1:Nskill){
    for(j in 1:Ntransiton){
      tstar[[i]][,j]=map_dbl(eta[[i]][,j],~rpg(num=1,h=1,z=.))
    }
  }
  
  kappa=map(1:Nskill,~matrix(NA,Nrespondents,Ntransition))
  for(i in 1:Nskill){
    for(j in 1:Ntransition){
      kappa[[i]][,j]=(trans_mat[,i]==j)-1/2
    }
  }
  
  C=map2(gamma_list,Xs,~(.y%*%.x)%>%
           {cbind(rep(0,Nrespondents),.)}%>%
           apply(1,function(x) log(sum(exp(x))-exp(x)))%>%
           t()
  )
  gamma_var=map(1:Nskill,list)
  priorsd_gamma=map(Xs,~map(.,function(x) priors$gamma_prior*diag(dim(x)[2])))
  
  gamma_out=gamma_list
  for(i in 1:Nskill){
    for(j in 2:Ntransition){
      priorsd_gamma[[i]]=10000
      gibbsvar_gamma=solve(crossprod(Xs[[i]]*tstar[[i]][,j,drop=F],Xs[[i]])+
                             solve(priorsd_gamma[[i]]))
      gibbsmean_gamma=gibbsvar_gamma%*%(t(Xs[[i]])%*%(kappa[[i]][,j]+tstar[[i]][,j,drop=F]*C[[i]][,j]))
      
      (t(Xs[[i]])%*%(kappa[[i]][,j]+tstar[[i]][,j,drop=F]*C[[i]][,j]))
      sum(kappa[[i]][,j]+tstar[[i]][,j,drop=F]*C[[i]][,j])
      
      sum(tstar[[i]][,j,drop=F]*C[[i]][,j])/sum(tstar[[i]][,j])
      
      gamma_out[[i]][,j-1]=sum(kappa[[i]][,j])/sum(tstar[[i]][,j])+
        mean(C[[i]][,j])
    }
  }
  
  gamma_list[[1]]
  gamma_out[[1]]
  Xtmp=list(matrix(c(1,1),2,1),
            matrix(c(1,1),2,1))
  gamma_to_probs(gamma_list,Xtmp)[[1]][1,]
  gamma_to_probs(gamma_out,Xtmp)[[1]][1,]
  rbind(prop.table(table(true_trans[,1])))
  
}

