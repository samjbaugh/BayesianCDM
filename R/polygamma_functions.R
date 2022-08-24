
gen_njc=function(profiles,y){
  njc=matrix(NA,Nquestions,Nprofile)
  for(j in 1:Nquestions){
    for(c in 1:Nprofile){
      njc[j,c]=sum(y[profiles==c,j])
    }
  }
  njc
}

gen_nci=function(profiles){
  nci=matrix(NA,Nprofile,Nrespondents)
  for(c in 1:Nprofile){
    for(i in 1:Nrespondents){
      nci[c,i]=sum(profiles[-i]==c)
    }
  }
  nci
}

sample_ystar=function(psi,nc){
  ystar=matrix(NA,Nquestions,Nprofile)
  for(j in 1:Nquestions){
    for(c in 1:Nprofile){
      ystar[j,c]=BayesLogit::rpg(h=nc[c],z=psi[j,c])
    }
  }
  ystar
}

sample_profile=function(r,t,theta,nmat){
  logp=sapply(1:Nprofile,function(p)
    sum(dbinom(X[[t]][r,],1,theta[,p],log=T)))
  probs=(nmat[,r]+1)*exp(logp-max(logp))/
    sum((nmat[,r]+1)*exp(logp-max(logp)))
  return(sample(Nprofile,1,prob=probs))
}

random_profiles=function(){
  return(sample(Nprofile,Nrespondents,rep=T))
}

beta_to_theta=function(x){
  logit(generate_logits_discrete(list(Ns=initparams$Ns,
                                      theta=x,
                                      value_key=initparams$value_key),
                                 q_info))
}

