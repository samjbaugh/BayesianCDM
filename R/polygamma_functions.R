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

#'  Convert alpha matrix to alpha matrix
#'
#' @description Add detail
#' @param alpha_mat alpha_mat
#' @param profiles Profile list
alpha_to_transitions=function(alpha_mat,profiles){
  Nrespondents=dim(alpha_mat)[1]
  Ntime=dim(alpha_mat)[2]
  Nskill=dim(profiles)[1]
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

#'  Sample tstar for transition model
#'
#' @description Add detail
#' @param eta See Polson paper
#' @importFrom BayesLogit rpg
sample_tstar=function(eta){
  Nskill=length(eta)
  Ntransiton=dim(eta[[1]])[2]
  tstar=map(1:Nskill,~matrix(NA,Nrespondents,2^Ntime))
  for(i in 1:Nskill){
    for(j in 1:Ntransiton){
      tstar[[i]][,j]=map_dbl(eta[[i]][,j],~rpg(num=1,h=1,z=.))
    }
  }
  tstar
}

#'  Convert gamma to transition probabilities
#'
#' @description Add detail
#' @param gamma_list List of gamma parameters
#' @param Xs list of design matrices
gamma_to_transprobs=function(gamma_list,Xs){
  Nrespondents=dim(Xs[[1]])[1]
  map2(gamma_list,Xs,~(.y%*%.x)%>%
         {cbind(rep(0,Nrespondents),.)}%>%
         apply(1,function(x) exp(x)/sum(exp(x)))%>%
         t())
}


#'  Sample gamma
#'
#' @description Add detail
#' @param tstar tstar
#' @param gamma_list gamma_list
#' @param trans_mat trans_mat
#' @param Xs Xs
#' @param priorsd_gamma priorsd_gamma
#' @param retmean retmean
#' @importFrom BayesLogit rpg
sample_gamma=function(tstar,gamma_list,trans_mat,Xs,priorsd_gamma,retmean=F){
  Nskill=length(tstar)
  Nrespondents=dim(tstar[[1]])[1]
  Ntransition=dim(tstar[[1]])[2]
  kappa=map(1:Nskill,~matrix(NA,Nrespondents,Ntransition))
  for(i in 1:Nskill){
    for(j in 1:Ntransition){
      kappa[[i]][,j]=(trans_mat[,i]==j)-1/2
    }
  }

  C=map2(gamma_list,Xs,~(.y%*%.x)%>%
             {cbind(rep(0,Nrespondents),.)}%>%
             apply(1,function(x) log(sum(exp(x))-exp(x)))%>%
             t())
  gamma_var=map(1:Nskill,list)
  for(i in 1:Nskill){
    for(j in 1:(Ntransition-1)){
      gibbsvar_gamma=solve(crossprod(Xs[[i]]*tstar[[i]][,j],Xs[[i]])+solve(priorsd_gamma[[i]]))
      gibbsmean_gamma=gibbsvar_gamma%*%(t(Xs[[i]])%*%(kappa[[i]][,j]+tstar[[i]][,j]*C[[i]][,j]))
      if(retmean){
        gamma_var[[i]][[j]]=gibbsvar_gamma
        gamma_list[[i]][,j]=gibbsmean_gamma
      }else{
        gamma_list[[i]][,j]=rmvnorm(1,gibbsmean_gamma,gibbsvar_gamma)
      }
    }
  }
  if(retmean){
    return(list(gamma_mean=gamma_list,gamma_var=gamma_var))
  }else{
    return(gamma_list)
  }
}

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
      iprof=alpha_mat[i,]
      theta_time=theta[qt_map==t,]
      joint_prob_alpha=rep(NA,Nprofile)
      logcond_prob_Y=rep(NA,Nprofile)
      for(c in 1:Nprofile){
        profile_sequence=alpha_mat[i,]
        profile_sequence[t]=c
        prof_seq=profiles[profile_sequence,]
        trans_ints=t(prof_seq)%*%nu+1
        joint_prob_alpha[c]=prod(map_dbl(1:Nskill,~trans_probs[[.]][i,trans_ints][.]))
        logcond_prob_Y[c]=sum((1-Ys[[t]][i,])*log(theta_time[,c])+
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


#'  Logistic function
#'
#' @description Add detail
#' @param x input
logistic=function(x) {
  1/(1+exp(-x))
}

