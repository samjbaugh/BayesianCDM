require(tidyverse)

sample_ystar_multinom=function(logits){
  ystar=list()
  for(j in 1:J){
    ystar[[j]]=map_dbl(logits[[j]],
                       ~BayesLogit::rpg(num=1,h=1,z=.))

  }
  ystar
}

require(BayesLogit)

multinom_sim=function(Nsim,m,Nrespondents,s,tunits,desmat,beta_init){
  desmat=desmat[1:(J-1)]

  Tmat = trans_mat(m)[[3]][,s]

  beta = beta_init
  V0j  = map(desmat,~diag(dim(.)[2]))
  ncs  = as.numeric(table(Tmat))

  all_beta=matrix(NA,Nsim,sum(map_dbl(beta,length)))
  all_beta[1,] = as.numeric(beta_init)
  for(nsim in 2:Nsim){
    #print(nsim)
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
    all_beta[nsim,]=unlist(mj)
  }
  return(all_beta)
}


multinom_transition_wrapper = function(trans_params)
{
  current_beta = lapply(1:Nskill, function(s) colMeans(trans_params[[s]][(Nsim/2):Nsim,]))
  suppressMessages({
  current_logits = lapply(1:Nskill, function(s)
                 bind_cols(map2(desmat[[s]][1:(J-1)],current_beta[[s]],~as.numeric(.x%*%.y)))%>%
                            cbind(rep(0,Nrespondents)))
  })
  current_prob = lapply(1:Nskill, function(s)
    t(apply(current_logits[[s]],1,function(x) exp(x-max(x))/sum(exp(x-max(x))))))

  no_skill  = list(lapply(1:Nskill, function(s) current_prob[[s]][,1] + current_prob[[s]][,4]),
                   lapply(1:Nskill, function(s) current_prob[[s]][,1] + current_prob[[s]][,3]))
  has_skill = list(lapply(1:Nskill, function(s) current_prob[[s]][,2] + current_prob[[s]][,3]),
                   lapply(1:Nskill, function(s) current_prob[[s]][,2] + current_prob[[s]][,4]))

  prof_probs = list(time1 = matrix(nrow = Nrespondents, ncol = Nprofile),
                    time2 = matrix(nrow = Nrespondents, ncol = Nprofile))
  for (t in 1:Ntime){
    for (p in 1:Nprofile){
      profilep_probs = c()
      for (s in 1:Nskill){
        profilep_probs = c(profilep_probs,
                           max((q_info$skill_in_profile[s,p]==0)*no_skill[[t]][[s]],
                                (q_info$skill_in_profile[s,p]==1)*has_skill[[t]][[s]]))
      }
      prof_probs[[t]][,p] = prod(profilep_probs)
    }
  }

  return(prof_probs)
}



# desmat=list(model.matrix((1:Nrespondents)~1),
#             model.matrix((1:Nrespondents)~1),
#             model.matrix((1:Nrespondents)~1),
#             model.matrix((1:Nrespondents)~1))
# desmat=lapply(1:Nskill, function(s) desmat)
# tunits=c("4", "3", "2", "1")
# J=length(tunits)
# trans_probs=lapply(1:Nskill, function(s) multinom_sim(M=100,Nrespondents,s=s,
#                                             tunits=tunits,desmat[[s]]))







