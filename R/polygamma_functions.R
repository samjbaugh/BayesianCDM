

# %%  modulo, remainder
# %/% division rounded down
gen_alpha=function(Nprofile,Nskill){
  alpha=matrix(NA,Nprofile,Nskill)
  for(ii in 1:Nprofile){
    tmp=rep(NA,Nskill)
    for(jj in 1:Nskill){
      tmp[Nskill-jj+1]=((ii-1)%/%(2^(jj-1)))%%2
    }
    alpha[ii,]=rev(tmp)
  }
  alpha
}

Q_to_delta=function(Q,maxiter='not yet'){
  Q=data.frame(Q)
  formula=paste0('~', paste0(names(Q),collapse='*'))
  as.matrix(model.matrix(data=Q,as.formula(formula)))
}

gen_A=function(alpha){
  alpha=data.frame(alpha)
  formula=paste0('~', paste0(names(alpha),collapse='*'))
  model.matrix(data=alpha,as.formula(formula))
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

sample_profile1=function(r,t,theta,nmat){
  logp=sapply(1:Nprofile,function(p)
    sum(dbinom(X[[t]][r,],1,theta[,p],log=T)))
  probs=(nmat[,r]+1)*exp(logp-max(logp))/
    sum((nmat[,r]+1)*exp(logp-max(logp)))
  return(sample(Nprofile,1,prob=probs))
}

sample_profile2=function(r,t,theta,prior){
  logp=sapply(1:Nprofile,function(p)
    sum(dbinom(X[[t]][r,],1,theta[,p],log=F)))
  # prior = nmat[,r]+1 from original paper,      nmat = nci
  # change to prof_probs here using transition model
  # probs=(nmat[,r]+1)*exp(logp-max(logp))/
  #   sum(nmat[,r]+1)*exp(logp-max(logp)))
  probs=(exp(logp-max(logp))/exp(logp-max(logp))) * prior[[t]][r,]
  return(sample(Nprofile,1,prob=probs))
}

sample_profile=function(r,t,theta,prior){
  logp=sapply(1:Nprofile,function(p)
    sum(dbinom(X[[t]][r,],1,theta[,p],log=F)))
  # prior = nmat[,r]+1 from original paper,      nmat = nci
  # change to prof_probs here using transition model
  # probs=(nmat[,r]+1)*exp(logp-max(logp))/
  #   sum(nmat[,r]+1)*exp(logp-max(logp)))
  probs=(exp(logp-max(logp))/exp(logp-max(logp))) * prior[[t]][r,]
  return(sample(Nprofile,1,prob=probs))
}

sample_ystar=function(psi,nc){
  ystar=matrix(NA,Nquestions,Nprofile)
  for(j in 1:Nquestions){
    for(c in 1:Nprofile){
      ystar[j,c]=BayesLogit::rpg(h=max(nc[c],1),z=psi[j,c])
    }
  }
  ystar
}

gen_njc=function(profiles,y){
  njc=matrix(NA,Nquestions,Nprofile)
  for(j in 1:Nquestions){
    for(c in 1:Nprofile){
      njc[j,c]=sum(y[profiles==c,j])
    }
  }
  njc
}

random_profiles=function(){
  return(sample(Nprofile,Nrespondents,rep=T))
}


get_vjp=function(ystar,A,sigma_jp){
  vjp=matrix(NA,Nquestions,Nprofile)
  for(j in 1:Nquestions){
    for(p in 1:Nprofile){
      omegaj=diag(ystar[j,])
      vjp[j,p]=sqrt(1/(t(A[,p])%*%(ystar[j,]*A[,p])+1/(sigma_jp[j,p]^2)))
    }
  }
  vjp
}

get_mjp=function(vjp,z,beta_mat){
  mjp=matrix(NA,Nquestions,Nprofile)
  for(j in 1:Nquestions){
    for(p in 1:Nprofile){
      omegaj=diag(ystar[j,])
      ztildej=z[j,]-A[,-p]%*%beta_mat[j,-p]
      mjp[j,p]=vjp[j,p]^2*t(A[,p])%*%(ystar[j,]*ztildej)
    }
  }
  mjp
}

logistic=function(x) {
  1/(1+exp(-x))
}

get_V_multinom=function(ystar,A,sigma_jp){
  vjp=matrix(NA,Nquestions,Nprofile)
  for(j in 1:Nquestions){
    for(p in 1:Nprofile){
      omegaj=diag(ystar[j,])
      vjp[j,p]=sqrt(1/(t(A[,p])%*%(ystar[j,]*A[,p])+1/(sigma_jp[j,p]^2)))
    }
  }
  vjp
}


trans_mat = function(m) {
  time1 = sapply(1:Nskill, function(s) prof_samples[[1]][,m] %in%
                   q_info$which_skill_profile[,s])
  time2 = sapply(1:Nskill, function(s) prof_samples[[2]][,m] %in%
                   q_info$which_skill_profile[,s])
  forward_mat = (time1 == 0) & (time2 == 1)
  backward_mat = (time1 == 1) & (time2 == 0)
  stay_same = forward_mat == 0 & backward_mat == 0
  forward_mat[stay_same] = backward_mat[stay_same] = NA
  mult_trans = sapply(1:Nskill, function(s)
    rowSums(matrix(c((time1[,s]==0 & time2[,s]==0)*1,
                     (time1[,s]==1 & time2[,s]==1)*2,
                     (time1[,s]==1 & time2[,s]==0)*3,
                     (time1[,s]==0 & time2[,s]==1)*4), ncol = Nskill)))
  return(list(forward_mat, backward_mat, mult_trans))
}


prof_transition = function(t = 2){

}



