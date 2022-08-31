
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




