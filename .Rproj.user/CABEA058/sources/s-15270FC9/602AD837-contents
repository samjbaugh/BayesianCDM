#' Generate level 2 residual values from beta and gamma parameters
#'
#' @param params Input values
#' @param Xdata Data list
#' @return List containing "uforward" and "ubackward" matrices
#' @export
gen_u=function(params,Xdata){
  gammaw=matrix(NA,params$Ns$Ngroup,params$Ns$Nrespcov+1)
  for(v in 1:(params$Ns$Nrespcov+1)){
    for(g in 1:params$Ns$Ngroup){
      gammaw[g,v]=params$gamma[v,]%*%Xdata$group_designmat[g,]
    }
  }
  uforward=array(apply(params$forward_betas,c(3,4),function(x) x-gammaw),
                 c(Ngroup,Nrespcov+1,Nskill,Ntime-1))
  ubackward=array(apply(params$forward_betas,c(3,4),function(x) x-gammaw),
                  c(Ngroup,Nrespcov+1,Nskill,Ntime-1))
  return(list(uforward=uforward,ubackward=ubackward))
}

#' Convert parameter list to vector
#'
#' @param params params
#' @export
convert_params_to_vector=function(params){
  fbeta_names=paste0('fbeta',c(params$beta_names))
  bbeta_names=paste0('bbeta',c(params$beta_names))

  vector= c(params$theta, c(params$forward_betas),c(params$backward_betas),c(params$gamma))
  varnames= c(params$theta_names,c(fbeta_names),c(bbeta_names),c(params$gamma_names))

  names(vector)=varnames
  return(vector)
}

#' Convert vector to parameter list
#'
#' @param vector params
#' @param oldparams Only needs to contain "Ns" entry
#' @param varnames Optional varnames
#' @export
convert_vector_to_params=function(vector,oldparams,varnames=NULL){
  if(is.null(varnames)){
    varnames=names(vector)
  }
  Nskill=oldparams$Ns$Nskill
  Ntime=oldparams$Ns$Ntime
  Nrespcov=oldparams$Ns$Nrespcov
  Ngroup=oldparams$Ns$Ngroup

  newparams=oldparams
  newparams[['theta']]=vector[startsWith(varnames,'theta')]
  newparams[['forward_betas']]=
    array(vector[startsWith(varnames,'fbeta')],c(Ngroup,Nrespcov+1,Nskill,Ntime-1))
  newparams[['backward_betas']]=
    array(vector[startsWith(varnames,'bbeta')],c(Ngroup,Nrespcov+1,Nskill,Ntime-1))
  newparams[['gamma']]=
    array(vector[startsWith(varnames,'gamma')],c(Nrespcov+1,Ngroupcov+1))

  return(newparams)
}


#' Generate list of profiles
#'
#' @param Nprofile Nprofile
gen_profile_list=function(Nprofile){
  profile_list=list()
  for(ii in 1:(2^Nskill)){
    tmp=rep(NA,Nskill)
    for(jj in 1:Nskill){
      tmp[Nskill-jj+1]=((ii-1)%/%(2^(jj-1)))%%2
    }
    profile_list[[ii]]=rev(tmp)
  }
  return(profile_list)
}

#' Logit function
#'
#' @param x x
logit=function(x) {
  1/(1+exp(-x))
}

#' Log logit
#'
#' @param x x
llogit=function(x) {
  log(1/(1+exp(-x)))
}

#' Log-sum-exp function
#'
#' @param x x
log_sum_exp<-function(x){
  c=max(x)
  return(c+log(sum(exp(x-c))))
}
