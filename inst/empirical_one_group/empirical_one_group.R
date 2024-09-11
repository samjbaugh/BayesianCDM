
#
#  Single-Group Empirical Dataset from Madison et al., 2018 (a)
#


  require(mvtnorm)
  require(truncnorm)
  require(tidyverse)
  
  setwd("~/Desktop/CDM/")
  devtools::load_all()
  setwd("~/Desktop/CDM/inst/empirical_one_group/")
  dat = read.table("prepost1.txt")[,-1]
  colnames(dat) = c(paste("Pre",1:21),paste("Post",1:21)) #gives items names for mirt
  dat[dat == 99] = NA                                     #replace 99 with NA's
  complete.dat     = dat[complete.cases(dat),]            #remove NA's
  complete.dat.mat = as.matrix(complete.dat)
  Ys            = list(pre = complete.dat.mat[,1:21], post = complete.dat.mat[,22:42])
  
  Ntime        = 2
  Nrespondents = nrow(complete.dat.mat)
  Nskill = k   = 4
  
  Xbase=cbind(rep(1,Nrespondents))
  # Xs=build_Xlist_01(Xbase,Nskill,Ntime)
  Xs=build_Xlist(Xbase,Nskill,Ntime)
  
  
  Q = matrix(0, nrow = ncol(complete.dat.mat)/2, ncol = 4)
  Q[c(1,13:14,17), 1] = Q[c(2:3,9:12), 2] =
    Q[c(4:8), 3] = Q[c(15:16, 18:21), 4] = 1
  Qs = list(Q, Q)
  
  Nq_list=map(Qs,~dim(.)[1])
  Nq_total=sum(unlist(Nq_list))
  Nprofile=2^Nskill
  
  runtime = system.time({
  sampler_out=fit_longitudinal_cdm_full(Ys,Xs,Qs,M=3000,
                                  initparams = NULL,
                                  priors=list(beta_prior=1,gamma_prior=0.5),
                                  fixed_beta = T)})



#################################################################################################################################################################
  
# rm burn-in
sampler_out$samples_burned_in = sampler_out$samples[500:3000,]

# Trace Plots
beta_samples = sampler_out$samples_burned_in%>%dplyr::select(contains('beta'))
gamma_samples = sampler_out$samples_burned_in%>%dplyr::select(contains('gamma'))

ntrace = 9
par(mfrow = c(ntrace,1))
set.seed(1)
beta_trace_i = sample(dim(beta_samples)[2], ntrace, replace=F)
beta_trace_samples = data.frame(Index = rep(paste0("Beta ",beta_trace_i), each=2501),
                                Value = unlist(beta_samples[,beta_trace_i]),
                                Iterations = rep(1:2501, ntrace))
ggplot(beta_trace_samples, aes(x = Iterations, y = Value)) + 
  theme_bw() + geom_line() + facet_wrap(~ Index)


set.seed(10)
gamma_trace_i = sample(dim(gamma_samples)[2], ntrace, replace=F)
gamma_trace_samples = data.frame(Index = rep(paste0("Gamma ", gamma_trace_i), each=2501),
                                 Value = unlist(gamma_samples[,gamma_trace_i]),
                                 Iterations = rep(1:2501, ntrace))
ggplot(gamma_trace_samples, aes(x = Iterations, y = Value)) + 
  theme_bw() + geom_line() + facet_wrap(~ Index)


par(mfrow = c(1,1))

  


#####################
#  Beta Posteriors  #
#####################
  {
    delta_cat=do.call(rbind,map(Qs,Q_to_delta))
    beta_samples = sampler_out$samples_burned_in%>%dplyr::select(contains('beta'))
    beta_mat=sampler_out$samples_burned_in%>%
      dplyr::select(contains('beta'))%>%
      apply(2,mean)
    beta_sd=sampler_out$samples_burned_in%>%
      dplyr::select(contains('beta'))%>%
      apply(2,sd)
    pbeta=
      data.frame(int=beta_mat[1:21],me=beta_mat[22:42],
                 int_sd=beta_sd[1:21],me_sd=beta_sd[22:42],i=1:21)%>%
      ggplot()+
      geom_point(aes(x=i,y=int,col='Intercepts'), shape = 1, size = 2)+
      geom_point(aes(x=i,y=me,col='Main Effects'), shape = 15, size = 2)+
      # geom_point(aes(x=i,y=postmean,col='postmean'))+
      geom_errorbar(aes(x=i,ymin=int-2*int_sd, ymax=int+2*int_sd,col='Intercepts')) +
      geom_errorbar(aes(x=i,ymin=me-2*me_sd, ymax=me+2*me_sd,col='Main Effects')) +
      # geom_hline(yintercept=0, col = "dodgerblue", size = 1.5) +
      coord_cartesian(xlim = c(0,21), ylim = c(-6.5, 7)) +
      theme(legend.position = "none",
            text = element_text(size = 17)) +
      xlab("Item Number") +
      ylab(expression(Beta ~ "Estimated Value"))
    # +geom_point(aes(x=i,y=betatrue,col='betatrue'))
    pbeta
  }
  beta_mat

##### beta point estimates
beta_order = c(1,13,14,17,2,3,9,10,11,12,4,5,6,7,8,15,16,18,19,20,21)
# plot(beta_mat[1:21], ylim = c(-6,6), xlab = "Item Number", ylab = "Estimated Value")
# points(beta_mat[22:42][order(beta_order)], pch = 20)
# plot(sampler_out$samples_burned_in[,39], type = "l")

betaplotdf = data.frame(int=beta_mat[1:21],
                        me=beta_mat[22:42][order(beta_order)],
                        i=1:21)
ggplot(betaplotdf)+
  geom_point(aes(x=i,y=int,col='Intercepts'), shape = 1, size = 3)+
  geom_point(aes(x=i,y=me,col='Main Effects'), shape = 15, size = 3)+
  coord_cartesian(ylim = c(-5, 5.5))+
  theme(text = element_text(size = 15))+
  labs(y = "Estimated Values", x = "Item Number", color = "")


################################
#  Transition Prob Posteriors  #
################################

##### conditional transition probabilities
gamma_samples = sampler_out$samples_burned_in%>%dplyr::select(contains('gamma'))
gamma_means = colMeans(gamma_samples)
gamma_mean = list(list(array(gamma_means[1], dim = c(1,1)),
                       array(gamma_means[2], dim = c(1,1)),
                       array(gamma_means[3], dim = c(1,1))),
                  list(array(gamma_means[4], dim = c(1,1)),
                       array(gamma_means[5], dim = c(1,1)),
                       array(gamma_means[6], dim = c(1,1))),
                  list(array(gamma_means[7], dim = c(1,1)),
                       array(gamma_means[8], dim = c(1,1)),
                       array(gamma_means[9], dim = c(1,1))),
                  list(array(gamma_means[10], dim = c(1,1)),
                       array(gamma_means[11], dim = c(1,1)),
                       array(gamma_means[12], dim = c(1,1))))

post_probs = lapply(1:Nskill, function(s) gamma_to_transprobs(gamma_mean, Xs)[[s]][1,])
cond_post_probs = list()
for (i in 1:Nskill)
{
  post_probs0 = sum(post_probs[[i]][1:2])
  post_probs1 = sum(post_probs[[i]][3:4])
  cond_post_probs[[i]] = post_probs[[i]] / c(post_probs0, post_probs0, post_probs1, post_probs1)
}
cond_post_probs


# prob dist
##### conditional transition probabilities
gamma_samples = sampler_out$samples_burned_in%>%dplyr::select(contains('gamma'))
prob_samples = matrix(nrow = nrow(sampler_out$samples_burned_in), ncol = 12)

for(i in 1:nrow(gamma_samples))
{
  gamma_i = unlist(gamma_samples[i,])
  gamma_i = list(list(array(c(gamma_i[1], gamma_i[2]), dim = c(2,1)),
                      array(gamma_i[3], dim = c(1,1)),
                      array(gamma_i[4], dim = c(1,1))),
                 list(array(c(gamma_i[5], gamma_i[6]), dim = c(2,1)),
                      array(gamma_i[7], dim = c(1,1)),
                      array(gamma_i[8], dim = c(1,1))),
                 list(array(c(gamma_i[9], gamma_i[10]), dim = c(2,1)),
                      array(gamma_i[11], dim = c(1,1)),
                      array(gamma_i[12], dim = c(1,1))),
                 list(array(c(gamma_i[13], gamma_i[14]), dim = c(2,1)),
                      array(gamma_i[15], dim = c(1,1)),
                      array(gamma_i[16], dim = c(1,1))))
  
  post_probs = lapply(1:Nskill, function(s) gamma_to_transprobs(gamma_i, Xs)[[s]][1,])
  cond_post_probs = list()
  for (j in 1:Nskill)
  {
    post_probs0 = sum(post_probs[[j]][1:2])
    post_probs1 = sum(post_probs[[j]][3:4])
    cond_post_probs[[j]] = post_probs[[j]] / c(post_probs0, post_probs0, post_probs1, post_probs1)
  }
  prob_samples[i,] = unlist(cond_post_probs)
}
apply(prob_samples, 2, sd)
colMeans(prob_samples)




######################
#  Gamma Posteriors  #
######################
{
  gamma_postmean=sampler_out$samples_burned_in%>%
    dplyr::select(contains('gamma'))%>%
    apply(2,mean)
  gamma_postsd=sampler_out$samples_burned_in%>%
    dplyr::select(contains('gamma'))%>%
    apply(2,sd)
  
  # trans_mat=alpha_to_transitions(true_alpha,profiles)
  # priorsd_gamma=map(Xs,~map(.,function(x) priors$gamma_prior*diag(dim(x)[2])))
  # tmp=sample_gamma(gamma_list,trans_mat,Xs,priorsd_gamma,retmean=T)
  plotdf=data.frame(
    # gamma_true=unlist(gamma_list_true),
    postmean=gamma_postmean,
    postsd=gamma_postsd,
    # pmean=unlist(tmp$gamma_mean),
    # psd=sqrt(unlist(map(tmp$gamma_var,~map(.,diag)))),
    i=1:length(gamma_postmean))%>%
    mutate(skill=c(rep("RPR",3),rep("MD",3),rep("NF",3),rep("GG",3)),
           i=factor(c(1:3,1:3,1:3,1:3)))
  pgamma=ggplot(plotdf)+
    geom_point(aes(x=i,y=postmean,col='posterior mean'))+
    geom_errorbar(aes(x=i,ymin=postmean-postsd*2,ymax=postmean+postsd*2,col='posterior mean'))+
    # geom_point(aes(x=i,y=gamma_true,col='true'))+
    # geom_point(aes(x=i+.25,y=pmean,col='posterior mean'))+
    # geom_errorbar(aes(x=i+.25,ymin=pmean-psd*2,ymax=pmean+psd*2,col='ss'))+
    scale_x_discrete(guide=guide_axis(angle = 90),
                     labels=c(expression(atop(NA, atop(0%->%1,"Intercept"))),
                              expression(atop(NA, atop(1%->%0,"Intercept"))),
                              expression(atop(NA, atop(1%->%1,"Intercept")))))+
    facet_grid(~factor(skill, levels=c("RPR","MD","NF","GG"))) +
    theme(legend.position = "none",
          text = element_text(size = 17)) +
    xlab("") +
    ylab(expression(Gamma ~ "Estimated Value"))
  pgamma
}




####################### 
#  Predicative Check  #
#######################

# respondent profile categorization at each time point
alpha_post=sampler_out$samples[2001:3000,]%>%
  dplyr::select(starts_with('alpha'))
alpha_post = list(alpha_post[,1:Nrespondents], alpha_post[,(Nrespondents+1):(2*Nrespondents)])

set.seed(1)
pred_data = 
  lapply(1:1000, function(iter) {
    lapply(1:Ntime, function(time) {
      t(sapply(1:Nrespondents, function(i) {
        rbinom(Nq_list[[1]], 1, prob = sampler_out$theta_list[[iter]][ ,alpha_post[[time]][iter,i]])} )
      )})})

pred_correct = 
  lapply(1:Ntime, function(time) {
    lapply(1:1000, function(iter) pred_data[[iter]][[time]] == Ys[[time]])
  })

# By time:
pred_qprob = map(pred_correct, ~colMeans(t(simplify2array(lapply(1:1000, function(m) colMeans(.[[m]]))))))
byq = map(pred_qprob, ~round(.,2))
byq
mean(pred_qprob[[2]])


# By respondents:
pred_iprob = map(pred_correct, ~colMeans(t(simplify2array(lapply(1:1000, function(m) rowMeans(.[[m]]))))))
resp_correct = data.frame(Prop = c(pred_iprob[[1]], pred_iprob[[2]]),
                          Time = c(rep("Time1", Nrespondents), rep("Time2", Nrespondents)))
ggplot(resp_correct, aes(Prop))+
  geom_histogram(fill = "deepskyblue", col = "black",alpha=0.3,binwidth=.03,position="identity")+
  facet_grid(Time ~ .) +
  xlab("Proportion Matching") +
  ylab("Respondent Count")

# By questions (prob correct instead of prob matching):
true_byq = map(Ys, ~colMeans(.))
pred_qprob = map(pred_data, ~map(.,~colMeans(.)))
pred_byq = list(t(sapply(1:1000, function(x) pred_qprob[[x]][[1]])),
                t(sapply(1:1000, function(x) pred_qprob[[x]][[2]])))
#time1
ggplot(data.frame(Question = rep(1:21,each=1000), Probability = c(pred_byq[[1]]), 
                  True = rep(true_byq[[1]],each=1000)), 
       aes(x=Question, y=Probability, group=Question)) +
  geom_boxplot(outlier.size=.05) +
  geom_point(aes(x=Question, y=True),col="red") +
  ylim(0,0.8) +
  ggtitle("Time 1") +
  theme(plot.title = element_text(colour = "steelblue", face = "bold", family = "Helvetica"),
        text = element_text(size = 17))
#time2
ggplot(data.frame(Question = rep(1:21,each=1000), Probability = c(pred_byq[[2]]), 
                  True = rep(true_byq[[2]],each=1000)), 
       aes(x=Question, y=Probability, group=Question)) +
  geom_boxplot(outlier.size=.05) +
  geom_point(aes(x=Question, y=True),col="red") +
  ylim(0,0.8) +
  ggtitle("Time 2") +
  theme(plot.title = element_text(colour = "steelblue", face = "bold", family = "Helvetica"),
        text = element_text(size = 17))

# By respondents (prob correct instead of prob matching):
# true_byr = map(Ys, ~rowMeans(.))
# pred_rprob = map(pred_data, ~map(.,~rowMeans(.)))
# pred_byr = list(t(sapply(1:1000, function(x) pred_rprob[[x]][[1]])),
#                 t(sapply(1:1000, function(x) pred_rprob[[x]][[2]])))
# #time1
# ggplot(data.frame(Respondents = rep(1:Nrespondents,each=1000), Probability = c(pred_byr[[1]]), 
#                   True = rep(true_byr[[1]],each=1000)), 
#        aes(x=Respondents, y=Probability, group=Respondents)) +
#   geom_boxplot(outlier.size=.05) +
#   geom_point(aes(x=Respondents, y=True),col="red")
# #time2
# ggplot(data.frame(Respondents = rep(1:Nrespondents,each=1000), Probability = c(pred_byr[[2]]), 
#                   True = rep(true_byr[[2]],each=1000)), 
#        aes(x=Respondents, y=Probability, group=Respondents)) +
#   geom_boxplot(outlier.size=.05) +
#   geom_point(aes(x=Respondents, y=True),col="red")


### CURRENTLY IN PAPER
# # thetas (from final 1000 iterations)
theta = Reduce("+", sampler_out$theta_list) / 1000

# sampler_out$theta_list

alpha_post=sampler_out$samples[2001:3000,]%>%
  dplyr::select(starts_with('alpha'))
profiles=map_dfr(1:Nprofile-1,
                 ~data.frame(t(rev(as.integer(intToBits(.))[1:Nskill]))))
M=dim(alpha_post)[1]
it_map=unlist(map(1:Ntime,~rep(.,Nrespondents)))
ti_map=unlist(map(1:Ntime,~1:Nrespondents))
post_probs=alpha_post%>%
  pivot_longer(everything()) %>%
  mutate(name=factor(name,levels=paste0('alpha_vec[',1:dim(alpha_post)[2],']')),
         value=factor(value,levels=1:Nprofile)) %>%
  group_by(name,value) %>%
  summarise(prop=n()/M) %>%
  complete(value,fill=list(prop=0))%>%mutate(ii=as.numeric(map(strsplit(as.character(name),'\\[|\\]'),2)),
                                             i=ti_map[ii],t=it_map[ii])%>%
  rowwise()
max_prob_index = sapply(1:(2*Nrespondents), function(x)
                        which.max(post_probs$prop[post_probs$ii==x])) +
                 seq(0,nrow(post_probs)-1, by = 16)
profile_id_mat = post_probs[max_prob_index,] %>%
                 mutate(value = as.numeric(value))
profile_id = list(profile_id_mat[1:Nrespondents,],
                  profile_id_mat[(Nrespondents+1):(2*Nrespondents),])
# 
# sampler_out$theta_list[[iter]][ ,alpha_post[[time]][iter,i]]
# sampler_out$theta_list[[1]][ ,alpha_post[[1]][1,1]]

# pred_data[[respondent]][[time]][iteration, question]
set.seed(1)
fitted_data = lapply(1:Ntime, function(time)
  t(sapply(1:Nrespondents, function (i) rbinom(Nq_list[[1]],1,theta[ ,profile_id[[time]][i,]$value])    )))

# # Merge
pred_correct = lapply(1:Ntime, function(time) fitted_data[[time]] == Ys[[time]])
item_correct = t(sapply(1:Ntime, function(t) colMeans(pred_correct[[t]])))
item_correct
round(item_correct, 2)
# 
resp_correct = data.frame(Prop = c(rowMeans(pred_correct[[1]]), rowMeans(pred_correct[[2]])),
                          Time = c(rep("Time1", Nrespondents), rep("Time2", Nrespondents)))
ggplot(resp_correct, aes(Prop))+
  geom_histogram(fill = "deepskyblue", col = "black",alpha=0.3,binwidth=.05,position="identity")+
  facet_grid(Time ~ .) +
  theme(text = element_text(size = 17))+
  xlab("Proportion Matching") +
  ylab("Respondent Count")

  


######################
#  Alpha Posteriors  #
######################  
{
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
    pivot_longer(everything()) %>%
    mutate(name=factor(name,levels=paste0('alpha_vec[',1:dim(alpha_post)[2],']')),
           value=factor(value,levels=1:Nprofile,labels=paste0('prof',1:Nprofile))) %>%
    group_by(name,value) %>%
    summarise(prop=n()/M) %>%
    complete(value,fill=list(prop=0))%>%mutate(ii=as.numeric(map(strsplit(as.character(name),'\\[|\\]'),2)),
                                               i=ti_map[ii],t=it_map[ii])%>%
    rowwise()
    # mutate(true_prof=true_alpha[i,t])
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