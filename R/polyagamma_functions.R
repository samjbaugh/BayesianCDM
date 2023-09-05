#'  Build design matrix
#'  with no covariates
#'
#' @description Add detail
#' @param Xbase Design matrix for (0,1) for each skill
#' @param Nskill Nskill
#' @param Ntime Ntime
build_Xlist=function(Xbase,Nskill,Ntime){
  X1=matrix(1,dim(Xbase)[1],1)
  retval=map(1:Nskill,list)
  for(i in 1:Nskill){
    for(j in 1:(2^Ntime-1)){
      retval[[i]][[j]]=X1
    }
  }
  retval
}

#'  Build design matrix
#'  with only 0->1 covariates
#'
#' @description Add detail
#' @param Xbase Design matrix for (0,1) for each skill
#' @param Nskill Nskill
#' @param Ntime Ntime
build_Xlist_01=function(Xbase,Nskill,Ntime){
  X1=matrix(1,dim(Xbase)[1],1)
  retval=map(1:Nskill,list)
  for(i in 1:Nskill){
    retval[[i]][[1]]=Xbase
    for(j in 2:(2^Ntime-1)){
      retval[[i]][[j]]=X1
    }
  }
  retval
}

#'  Build design matrix
#'  with 0->1, 1->0 covariates
#'
#' @description Add detail
#' @param Xbase Design matrix for (0,1) for each skill
#' @param Nskill Nskill
#' @param Ntime Ntime
build_Xlist_01_10=function(Xbase,Nskill,Ntime){
  X1=matrix(1,dim(Xbase)[1],1)
  retval=map(1:Nskill,list)
  for(i in 1:Nskill){
    retval[[i]][[1]]=Xbase
    retval[[i]][[2]]=Xbase
    for(j in 3:(2^Ntime-1)){
      retval[[i]][[j]]=X1
    }
  }
  retval
}

#'  Convert transition matrix to alpha matrix
#'
#' @description Add detail
#' @param trans_mat transition matrix
#' @param transitions list of transitions
transitions_to_alpha=function(trans_mat,transitions){
  Nrespondents=dim(trans_mat)[1]
  Nskill=dim(trans_mat)[2]
  Ntime=dim(transitions)[2]
  retval=matrix(NA,Nrespondents,Ntime)
  v=2^(Nskill:1-1)
  for(i in 1:Nrespondents){
    tmpmat=matrix(NA,Nskill,Ntime)
    for(k in 1:Nskill){
      tmpmat[k,]=as.numeric(transitions[trans_mat[i,k],])
    }
    retval[i,]=t(tmpmat)%*%v+1
  }
  return(retval)
}

#'  Build design matrix
#'  with full transition covariates
#'
#' @description Add detail
#' @param Xbase Design matrix for (0,1) for each skill
#' @param Nskill Nskill
#' @param Ntime Ntime
build_Xlist_full=function(Xbase,Nskill,Ntime){
  X1=matrix(1,dim(Xbase)[1],1)
  retval=map(1:Nskill,list)
  for(i in 1:Nskill){
    retval[[i]][[1]]=Xbase
    retval[[i]][[2]]=Xbase
    retval[[i]][[3]]=Xbase
  }
  retval
}

#'  Convert alpha matrix to transition matrix
#'
#' @description Add detail
#' @param alpha_mat alpha_mat
#' @param profiles Profile list
alpha_to_transitions=function(alpha_mat,profiles){
  Nrespondents=dim(alpha_mat)[1]
  Ntime=dim(alpha_mat)[2]
  Nskill=dim(profiles)[2]
  retval=matrix(NA,Nrespondents,Nskill)
  nu=2^(Ntime:1-1)
  for(i in 1:Nrespondents){
    tmpmat=matrix(NA,Ntime,Nskill)
    for(k in 1:Ntime){
      tmpmat[k,]=as.numeric(profiles[alpha_mat[i,k],])
    }
    retval[i,]=t(tmpmat)%*%nu+1
  }
  return(retval)
}

#'  Convert Q to delta
#'
#' @description Add detail
#' @param Q Q matrix
#' @importFrom stats model.matrix as.formula
Q_to_delta=function(Q){
  Q=data.frame(Q)
  formula=paste0('~', paste0(names(Q),collapse='*'))
  as.matrix(model.matrix(data=Q,as.formula(formula)))
}

#'  Sample polya-gamma auxillary variables
#'
#' @description Add detail
#' @param psi See Balamuta paper
#' @param nct See Balamuta paper
#' @param qt_map Map from q to t
#' @importFrom BayesLogit rpg
sample_ystar=function(psi,nct,qt_map){
  Nq_total=dim(psi)[1]
  Nprofile=dim(psi)[2]
  ystar=matrix(NA,Nq_total,Nprofile)
  for(j in 1:Nq_total){
    for(c in 1:Nprofile){
      ystar[j,c]=rpg(h=max(nct[c,qt_map[j]],1),z=psi[j,c])
    }
  }
  ystar
}

#'  Convert gamma to transition probabilities
#'
#' @description Add detail
#' @param gamma_list List of gamma parameters
#' @param Xs list of design matrices
gamma_to_transprobs=function(gamma_list,Xs) apply_gamma_X(gamma_list,Xs,function(x) exp(x)/sum(exp(x)))

apply_gamma_X=function(gamma_list,Xs,f){
  Ntransition=length(gamma_list[[1]])+1
  Nskill=length(Xs)
  Nr=dim(Xs[[1]][[1]])[1]
  retval=list()
  for(k in 1:Nskill){
    a=do.call(cbind,map(1:(Ntransition-1),~Xs[[k]][[.]]%*%gamma_list[[k]][[.]]))
    tmp=cbind(rep(0,Nr),a)
    retval[[k]]=t(apply(tmp,1,f))
  }
  retval
}
# gamma_to_transprobs=function(gamma_list,Xs){
#   Nrespondents=dim(Xs[[1]])[1]
#   map2(gamma_list,Xs,~(.y%*%.x)%>%
#          {cbind(rep(0,Nrespondents),.)}%>%
#          apply(1,function(x) exp(x)/sum(exp(x)))%>%
#          t())
# }



#'  Sample alpha matrix
#'
#' @description Add detail
#' @param alpha_mat Alpha matrix
#' @param theta Response probabilities
#' @param trans_probs Transition probabilities
#' @param qt_map Map between questions and time
#' @param profiles List of profiles
#' @param Ys Data
sample_alpha_mat=function(alpha_mat,theta,trans_probs,qt_map,profiles,Ys){
  Nskill=dim(profiles)[2]
  Ntime=length(Ys)
  nu=2^(Ntime:1-1)
  Nprofile=dim(theta)[2]
  for(i in 1:dim(alpha_mat)[1]){
    for(t in 1:dim(alpha_mat)[2]){
      # theta_time=theta[qt_map==t,]
      theta_time = theta
      joint_prob_alpha=rep(NA,Nprofile)
      logcond_prob_Y=rep(NA,Nprofile)
      for(c in 1:Nprofile){
        profile_sequence=alpha_mat[i,]
        profile_sequence[t]=c
        prof_seq=profiles[profile_sequence,]
        trans_ints=t(prof_seq)%*%nu+1
        joint_prob_alpha[c]=
          prod(map_dbl(1:Nskill,~trans_probs[[.]][i,trans_ints][.]))
        logcond_prob_Y[c]=sum((1-Ys[[t]][i,])*log(1-theta_time[,c])+
                                Ys[[t]][i,]*log(theta_time[,c]))
      }
      cond_prob_alpha=joint_prob_alpha/sum(joint_prob_alpha)
      punorm=log(cond_prob_alpha)+logcond_prob_Y
      samp_prob_alpha=exp(punorm-max(punorm))/sum(exp(punorm-max(punorm)))
      alpha_mat[i,t]=sample(1:Nprofile,1,prob=samp_prob_alpha)
    }
  }
  alpha_mat
}

#'  Return mean and variance of beta's Gibbs distribution
#'
#' @description Add detail
#' @param ystar sampled auxiliary values
#' @param A Profile design matrix
#' @param z z values
#' @param beta_mat Old beta_mat
#' @param priorsd_beta Variance priors
beta_gibbs_dist=function(ystar,A,z,beta_mat,priorsd_beta){
  Nq_total=dim(ystar)[1]
  Nprofile=dim(ystar)[2]
  condsd=matrix(NA,Nq_total,Nprofile)
  for(j in 1:Nq_total){
    for(p in 1:Nprofile){
      omegaj=diag(ystar[j,])
      condsd[j,p]=sqrt(1/(t(A[,p])%*%(ystar[j,]*A[,p])+1/(priorsd_beta[j,p]^2)))
    }
  }

  condmean=matrix(NA,Nq_total,Nprofile)
  for(j in 1:Nq_total){
    for(p in 1:Nprofile){
      omegaj=diag(ystar[j,])
      ztildej=z[j,]-A[,-p]%*%beta_mat[j,-p]
      condmean[j,p]=condsd[j,p]^2*t(A[,p])%*%(ystar[j,]*ztildej)
    }
  }
  return(list(mean=condmean,sd=condsd))
}


#'  Sample gamma
#'
#' @description Add detail
#' @param gamma_list gamma_list
#' @param trans_mat trans_mat
#' @param Xs Xs
#' @param priorsd_gamma priorsd_gamma
#' @param retmean retmean
#' @importFrom BayesLogit rpg
sample_gamma=function(gamma_list,trans_mat,Xs,priorsd_gamma,retmean=F){
  Nskill=length(gamma_list)
  Nrespondents=dim(Xs[[1]][[1]])[1]
  Ntransition=length(gamma_list[[1]])+1
  Ntime=log2(Ntransition)
  
  eta=apply_gamma_X(gamma_list,Xs,function(x) x-log(sum(exp(x))-exp(x)))
  tstar=map(1:Nskill,~matrix(NA,Nrespondents,2^Ntime))
  for(i in 1:Nskill){
    for(j in 1:Ntransition){
      tstar[[i]][,j]=map_dbl(eta[[i]][,j],~rpg(num=1,h=1,z=.))
    }
  }
  
  kappa=map(1:Nskill,~matrix(NA,Nrespondents,Ntransition))
  for(i in 1:Nskill){
    for(j in 1:Ntransition){
      kappa[[i]][,j]=(trans_mat[,i]==j)-1/2
    }
  }
  
  C=apply_gamma_X(gamma_list,Xs,function(x) log(sum(exp(x))-exp(x)))
  gamma_var=map(1:Nskill,list)
  for(i in 1:Nskill){
    for(j in 2:Ntransition){
      gibbsvar_gamma=solve(crossprod(Xs[[i]][[j-1]]*tstar[[i]][,j],Xs[[i]][[j-1]])+
                             solve(priorsd_gamma[[i]][[j-1]]))
      gibbsmean_gamma=gibbsvar_gamma%*%(t(Xs[[i]][[j-1]])%*%(kappa[[i]][,j]+tstar[[i]][,j]*C[[i]][,j]))
      if(retmean){
        gamma_var[[i]][[j-1]]=gibbsvar_gamma
        gamma_list[[i]][[j-1]]=gibbsmean_gamma
      }else{
        gamma_list[[i]][[j-1]]=c(rmvnorm(1,gibbsmean_gamma,gibbsvar_gamma))
      }
    }
  }
  if(retmean){
    return(list(gamma_mean=gamma_list,gamma_var=gamma_var))
  }else{
    return(gamma_list)
  }
}


# tt=function(ystar,A,z,beta_mat,priorsd_beta){
#   Nq_total=dim(ystar)[1]
#   Nprofile=dim(ystar)[2]
#   condsd=matrix(NA,Nq_total,Nprofile)
#   for(j in 1:Nq_total){
#     for(p in 1:Nprofile){
#       omegaj=diag(ystar[j,])
#       condsd[j,p]=sqrt(1/(t(A[,p])%*%(ystar[j,]*A[,p])+1/(priorsd_beta[j,p]^2)))
#     }
#   }
# 
#   condmean=matrix(NA,Nq_total,Nprofile)
#   for(j in 1:Nq_total){
#     for(p in 1:Nprofile){
#       omegaj=diag(ystar[j,])
#       ztildej=z[j,]-A[,-p]%*%beta_mat[j,-p]
#       condmean[j,p]=condsd[j,p]^2*t(A[,p])%*%(ystar[j,]*ztildej)
#     }
#   }
#   return(list(mean=condmean,sd=condsd))
# }


#'  Logistic function
#'
#' @description Add detail
#' @param x input
logistic=function(x) {
  1/(1+exp(-x))
}




################## fixed_beta functions ################## 

#'  Sample polya-gamma auxillary variables (fixed_beta)
#'
#' @description Add detail
#' @param psi See Balamuta paper
#' @param nct See Balamuta paper
#' @importFrom BayesLogit rpg
sample_ystar_fb=function(psi,nct){
  Nq=dim(psi)[1]
  Nprofile=dim(psi)[2]
  ystar=map(1:Ntime,~matrix(NA,Nq,Nprofile))
  for(t in 1:Ntime){
    for(j in 1:Nq){
      for(c in 1:Nprofile){
        ystar[[t]][j,c]=rpg(h=max(nct[c,t],1),z=psi[j,c])
      }
    }
  }
  ystar
}

#'  Return mean and variance of beta's Gibbs distribution (fixed_beta)
#'
#' @description Add detail
#' @param ystar sampled auxiliary values
#' @param A Profile design matrix
#' @param z z values
#' @param beta_mat Old beta_mat
#' @param priorsd_beta Variance priors
beta_gibbs_dist_fb=function(ystar,A,z,beta_mat,priorsd_beta){
  Nq=dim(ystar[[1]])[1]
  Nprofile=dim(ystar[[1]])[2]
  condsd=matrix(NA,Nq,Nprofile)
  for(j in 1:Nq){
    for(p in 1:Nprofile){
      omegaj=unlist(map(1:Ntime,~ystar[[.]][j,]))
      Aconcat=do.call(rbind, replicate(Ntime, A, simplify=FALSE)) 
      condsd[j,p]=sqrt(1/(t(Aconcat[,p])%*%(omegaj*Aconcat[,p])+1/(priorsd_beta[j,p]^2)))
    }
  }
  condmean=matrix(NA,Nq,Nprofile)
  for(j in 1:Nq){
    for(p in 1:Nprofile){
      omegaj=unlist(map(1:Ntime,~ystar[[.]][j,]))
      zj=c(unlist(map(z,~.[j,])))
      ztildej=zj-Aconcat[,-p]%*%beta_mat[j,-p]
      condmean[j,p]=condsd[j,p]^2*t(Aconcat[,p])%*%(omegaj*ztildej)
    }
  }
  return(list(mean=condmean,sd=condsd))
}















