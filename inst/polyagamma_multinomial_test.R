require(tidyverse)

Nrespondents=100

sample_ystar_multinom=function(logits){
  ystar=list()
  for(j in 1:J){
    ystar[[j]]=map_dbl(logits[[j]],
                       ~BayesLogit::rpg(num=1,h=1,z=.))

  }
  ystar
}

require(BayesLogit)

multinom_sim=function(M,Nrespondents,tunits,desmat,id=''){
  J=length(tunits)
  desmat=desmat[1:(J-1)]

  true_beta=map(desmat,~rnorm(dim(.)[2]))%>%
    set_names(paste0('beta',tunits[1:(J-1)]))
  true_logits=bind_cols(map2(desmat,true_beta,~as.numeric(.x%*%.y)))%>%
    cbind(rep(0,Nrespondents))
  true_probs=t(apply(true_logits,1,function(x) exp(x-max(x))/sum(exp(x-max(x)))))

  Tmat=t(t(apply(true_probs,1,function(x) sample(tunits,1,prob=x))))

  beta_init=map(desmat,~rnorm(dim(.)[2]))
  beta=true_beta

  V0j=map(desmat,~diag(dim(.)[2]))
  ncs=as.numeric(table(Tmat))

  all_beta=matrix(NA,M,sum(map_dbl(beta,length)))
  for(m in 1:M){
    print(m)
    suppressMessages({
      mylogits=map2(desmat,beta,~as.numeric(.x%*%.y))%>%
        bind_cols()%>%
        cbind(0)
    })
    C=t(apply(mylogits,1,function(x) log(sum(exp(x))-exp(x))))
    kappa=map(1:J,~((Tmat==tunits[[.]])-1/2))

    eta=mylogits-C
    ystar=sample_ystar_multinom(eta)

    Vj=map(1:(J-1),
           ~solve(crossprod(desmat[[.]]*ystar[[.]],desmat[[.]])+solve(V0j[[.]])))
    mj=map(1:(J-1),
           ~Vj[[.]]%*%(t(desmat[[.]])%*%(kappa[[.]]+ystar[[.]]*C[,.,drop=F])))

    beta=map2(mj,Vj,~mvtnorm::rmvnorm(1,.x,.y))
    all_beta[m,]=unlist(mj)
  }
  baseline_category=paste0('beta',tunits[J])
  beta_df=all_beta%>%
    data.frame()%>%
    set_names(paste0('beta',tunits[1:(J-1)]))%>%
    cbind(set_names(data.frame(tmp=0),baseline_category))%>%
    mutate(iter=1:M)%>%
    na.omit()%>%
    pivot_longer(!iter,values_to='beta')
  ggplot(data=beta_df%>%filter(name!=baseline_category))+
    aes(x=iter,y=beta)+geom_point()+
    facet_grid(~name)
  true_beta_df=true_beta%>%
    bind_rows()%>%
    pivot_longer(everything(),values_to='trueval')%>%
    rbind(data.frame(name=baseline_category,trueval=0))%>%
    mutate(trueprob=exp(trueval)/sum(exp(trueval)))
  emp_prob=prop.table(table(Tmat))%>%
    as_tibble()%>%
    mutate(name=paste0('beta',Tmat))%>%
    rename(emp_prob='n')%>%
    select(name,emp_prob)
  postsumm=beta_df%>%
    group_by(iter)%>%
    mutate(prob=exp(beta)/sum(exp(beta)))%>%
    group_by(name)%>%
    summarise(mbeta=mean(beta),
              varbeta=var(beta),
              mprob=mean(prob),
              q05prob=quantile(prob,.05),
              q95prob=quantile(prob,.95))%>%
    right_join(true_beta_df,by='name')%>%
    right_join(emp_prob,by='name')%>%
    mutate(id="trial1")
  return(postsumm)
}

desmat=list(model.matrix((1:Nrespondents)~1),
            model.matrix((1:Nrespondents)~1),
            model.matrix((1:Nrespondents)~1))
J=3
simresults=
  map_dfr(paste0('trial',1:5),~multinom_sim(M=10,Nrespondents=100,
                                            tunits=c('1','2','3','4')[1:J], #'11','00'),
                                            desmat[1:(J-1)],id=.))

simresults%>%
  ggplot()+
  geom_errorbar(aes(x=id,ymin=q05prob,ymax=q95prob,col='postsamp'))+
  geom_point(aes(x=id,y=mprob,col='postsamp'))+
  geom_point(aes(x=id,y=emp_prob,col='emp_prob'))+
  geom_point(aes(x=id,y=trueprob,col='trueprob'))+
  facet_grid(~name)
