require(tidyverse)

Nrespondents=1000

sample_ystar_multinom=function(xiTbeta,ncs){
  ystar=list()
  for(j in 1:J){
    ystar[[j]]=map_dbl(xiTbeta[[j]],
                       ~BayesLogit::rpg(num=1,h=ncs[[j]],z=.))

  }
  ystar
}

#functions
gen_C=function(xiTbeta){
  C=matrix(NA,Nrespondents,J)
  for(i in 1:Nrespondents){
    for(j in 1:J){
      C[i,j]=log_sum_exp(map_dbl(xiTbeta,i)[-j])
    }
  }
  C
}

gen_mVj=function(ystar,kappa){
  Vj=list()
  mj=list()
  for(j in 1:(J-1)){
    Omegaj=diag(ystar[[j]])
    Vj[[j]]=solve(t(desmat[[j]])%*%Omegaj%*%desmat[[j]]+solve(V0j[[j]]))
    mj[[j]]=Vj[[j]]%*%(t(desmat[[j]])%*%(kappa[[j]]-Omegaj%*%C[,j]))
  }
  return(list(Vj=Vj,mj=mj))
}



#multinomial example (one time point four categories for simplicity)
respondent_covariates=rnorm(Nrespondents)
tunits=c('01','10','11','00')
J=length(tunits)

desmat=list(model.matrix((1:Nrespondents)~1),
            model.matrix((1:Nrespondents)~1),
            model.matrix((1:Nrespondents)~1))

true_beta=map(desmat,~rnorm(dim(.)[2]))
true_logits=bind_cols(map2(desmat,true_beta,~as.numeric(.x%*%.y)))%>%
  cbind(rep(0,Nrespondents))
true_probs=t(apply(true_logits,1,function(x) exp(x-max(x))/sum(exp(x-max(x)))))

respondent_covariates=rnorm(Nrespondents)

beta_init=map(desmat,~rnorm(dim(.)[2]))
beta=beta_init



M=2000

V0j=map(desmat,~diag(dim(.)[2]))
ncs=as.numeric(table(Tmat))

all_beta=matrix(NA,sum(map_dbl(beta,length)),M)
for(m in 1:M){
  print(m)
  suppressMessages({mylogits=map2(desmat,beta,~as.numeric(.x%*%.y))%>%
    bind_cols()%>%
    cbind(0)})
  C=gen_C(mylogits)
  #ni is 1/2 because only one transition is observed
  kappa=map(1:J,~((Tmat==tunits[[.]])-1/2))
  ystar=sample_ystar_multinom(mylogits,ncs)
  mVj=gen_mVj(ystar,kappa)

  beta=map(1:(J-1),~c(mvtnorm::rmvnorm(1,mVj$mj[[.]],mVj$Vj[[.]])))
  all_beta[,m]=unlist(beta)
}
plotdf=all_beta%>%t()%>%
  data.frame()%>%
  set_names(c('beta01','beta10','beta11-base'))%>%
  mutate(t=1:M)%>%
  na.omit()%>%
  pivot_longer(!t)
true_vals=true_beta%>%unlist%>%t%>%
  data.frame()%>%
  set_names(c('beta01','beta10','beta11-base'))%>%
  pivot_longer(everything())

ggplot(data=plotdf)+
  geom_histogram(aes(x=value,col=name,fill=name))+
  geom_vline(data=true_vals,aes(xintercept=value,col=name))+
  facet_grid(~name,scales='free')

