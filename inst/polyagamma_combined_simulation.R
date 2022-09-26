# simulation for the polyagamma combined model
require(tidyverse)
require(truncnorm)
require(dplyr)

setwd("~/Documents/GitHub/BayesianCDM/")
devtools::load_all()

set.seed(1000)
#########################
#                       #
#     Simulate Data     #
#                       #
#########################

# Data parameters
Nrespondents=200
Nquestions=21
Ngroup=1
Ntime=2
myseed=2171506
Nskill=3
Nprofile=2^Nskill
Nrespcov=0
Ngroupcov=0
Q=gen_profile_list(Nprofile)%>%.[-1]%>%
  {do.call(rbind,lapply(.,function(x)
    do.call(rbind,lapply(1:ceiling(Nquestions/(Nprofile-1)),function(i) x))))%>%
      .[1:Nquestions,]}
desmat=list(model.matrix((1:Nrespondents)~1),
            model.matrix((1:Nrespondents)~1),
            model.matrix((1:Nrespondents)~1),
            model.matrix((1:Nrespondents)~1))
desmat=lapply(1:Nskill, function(s) desmat)
tunits=c("4", "3", "2", "1")
J=length(tunits)

# How true_gamma is generated in gen_data_wrapper
# true_gamma=map(desmat[[1]][1:(J-1)],~rnorm(dim(.)[2]))%>%
#   set_names(paste0('gamma',tunits[1:(J-1)]))
# true_logits=bind_cols(map2(desmat[[1]][1:(J-1)],true_gamma,~as.numeric(.x%*%.y)))%>%
#   cbind(rep(0,Nrespondents))
# true_probs=t(apply(true_logits,1,function(x) exp(x-max(x))/sum(exp(x-max(x)))))

# Simulate data using parameters
simdata=gen_data_wrapper(Nrespondents=Nrespondents,
                         Nquestions=Nquestions,
                         Nskill=Nskill,
                         Ntime=Ntime,
                         Ngroup=Ngroup,
                         Nrespcov=Nrespcov,
                         Ngroupcov=Ngroupcov,
                         multinomial = T)

Xdata      = simdata$Xdata
sum(Xdata$Xs[[1]])
sum(Xdata$Xs[[2]])


initparams = params = simdata$true_params
q_info     = gen_q_info(simdata$Xdata$Q)
X          = list(pre = as.matrix(Xdata$Xs[[1]]), post = as.matrix(Xdata$Xs[[2]]))
Ns         = list(Nrespondents=Nrespondents,Ntime=Ntime,
                  Nquestions=Nquestions,Nrespcov=Nrespcov,
                  Ngroupcov=Ngroupcov,Ngroup=Ngroup,
                  Nskill=Nskill,Nprofile=Nprofile)

# gamma columns matching transition type
# [,1] 0->1
# [,2] 1->0
# [,3] 1->1
# [,4] 0->0

alpha=gen_alpha(Nprofile,Nskill)
A=gen_A(alpha)

Delta=Q_to_delta(Q)
beta_mat=abs(matrix(rnorm(Nquestions*Nprofile,sd=.5),Nquestions,Nprofile))*Delta
beta_mat[,1]=rnorm(Nquestions,sd=.25)
beta_vec=beta_mat[Delta==1]

psi=beta_mat%*%t(A)
theta=logistic(psi)
sigma_jp=matrix(1,Nquestions,Nprofile)
#Enforce positivity in all but intercept, as is done in the paper
Ljp=matrix(-Inf,Nquestions,Nprofile)
Ljp[,1]=-Inf



##########################################
#                                        #
#     Run Polya-Gamma Integrated DCM     #
#                                        #
##########################################

myt=system.time({
  M=1000                                 # iterations for item-response model
  Nsim=10                                # iterations for transition model
  beta_samples=matrix(NA,length(beta_vec),M)
  beta_samples[,1]=c(beta_vec)

  prof_samples=list(matrix(NA,Nrespondents,M),
                    matrix(NA,Nrespondents,M))
  prof_samples[[1]][,1]=random_profiles()
  prof_samples[[2]][,1]=random_profiles()
  trans_m = list(matrix(NA,nrow=M,ncol=J-1),
                 matrix(NA,nrow=M,ncol=J-1),
                 matrix(NA,nrow=M,ncol=J-1))

  for(m in 2:M){

    #compute transition probabilities
    if(m==2){
      gamma_init=lapply(1:Nskill, function(s) map(desmat[[s]][1:(J-1)],~rnorm(dim(.)[2])))
    } else {
      gamma_init=lapply(1:Nskill, function(s) as.list(trans_params[[s]][Nsim,]))
    }
    trans_params=lapply(1:Nskill, function(s)
      multinom_sim(Nsim,m,Nrespondents,s=s,
      tunits=tunits,desmat[[s]],
      gamma_init=gamma_init[[s]]))
    prof_probs=multinom_transition_wrapper(trans_params)

    for (s in 1:Nskill) {
      trans_m[[s]][m,] = colMeans(trans_params[[s]][(Nsim/2):Nsim,])
    }

    for(t in 1:Ntime){

      #sample profile
      nci=gen_nci(prof_samples[[t]][,m-1])
      prof_sample=sapply(1:Nrespondents,function(r) sample_profile2(r,t,theta,prof_probs))
      prof_samples[[t]][,m]=prof_sample

      #sample augmented data
      nc = as.numeric(sapply(1:Nprofile, function(p) sum(prof_sample==p)))
      ystar=sample_ystar(psi,nc)

      #sample beta
      njc=gen_njc(prof_samples[[t]][,m],X[[t]])
      kappa=sweep(njc,2,nc/2,'-') #t(apply(nmat,1,function(x) x-sum(x)/2))
      z=kappa/ystar

      vjp=get_vjp(ystar,A,sigma_jp)
      mjp=get_mjp(vjp,z,beta_mat)

      beta_mat=matrix(rtruncnorm(length(mjp),mean=c(mjp),sd=c(vjp),a=c(Ljp)),
                      dim(beta_mat)[1],dim(beta_mat)[2])*Delta


      psi=beta_mat%*%t(A)
      theta=logistic(psi)

      beta_vec=beta_mat[Delta==1]
      beta_samples[,m]=c(beta_vec)
    }
    print(m)
  }
})

save.image(file='M1000_Nsim10.RData')

# True gammas converted to probabilies
lapply(1:Nskill, function(s) simdata$true_params$true_gamma_probs[[s]][1,])
# Result gammas converted to probabilities
result_gamma = lapply(1:Nskill, function(s) as.list(colMeans(trans_params[[s]])))
result_gamma_logits = lapply(1:Nskill, function(s) bind_cols(map2(desmat[[1]][1:(J-1)],
                               result_gamma[[s]],~as.numeric(.x%*%.y)))%>%cbind(rep(0,Nrespondents)))
result_gamma_probs = lapply(1:Nskill, function(s) t(apply(result_gamma_logits[[s]],1,function(x)
  exp(x-max(x))/sum(exp(x-max(x))))))
lapply(1:Nskill, function(s) result_gamma_probs[[s]][1,])

# Using trans_m to see how gamma changes with m
select_skill = 2
select_gamma = 3
plot(trans_m[[select_skill]][,select_gamma])
plot(beta_samples[6,])
plot(trans_params[[1]][,3])

# Compare points estimate beta results with true betas
plot(colMeans(t(beta_samples)), ylim = c(-3,5))
points(simdata$true_params$beta, pch = 20)




# baseline_category=paste0('beta',tunits[J])
# beta_df=trans_params[[1]]%>%
#   data.frame()%>%
#   set_names(paste0('beta',tunits[1:(J-1)]))%>%
#   cbind(set_names(data.frame(tmp=0),baseline_category))%>%
#   mutate(iter=1:Nsim)%>%
#   na.omit()%>%
#   pivot_longer(!iter,values_to='beta')
# ggplot(data=beta_df%>%filter(name!=baseline_category))+
#   aes(x=iter,y=beta)+geom_point()+
#   facet_grid(~name)
# Tmat = trans_mat(m)[[3]][,1]
# emp_prob=prop.table(table(Tmat))%>%
#   as_tibble()%>%
#   mutate(name=paste0('beta',Tmat))%>%
#   rename(emp_prob='n')%>%
#   select(name,emp_prob)
# postsumm=beta_df%>%
#   group_by(iter)%>%
#   mutate(prob=exp(beta)/sum(exp(beta)))%>%
#   group_by(name)%>%
#   summarise(mbeta=mean(beta),
#             varbeta=var(beta),
#             mprob=mean(prob),
#             q05prob=quantile(prob,.05),
#             q95prob=quantile(prob,.95))%>%
#   right_join(emp_prob,by='name')%>%
#   mutate(id="trial1")
