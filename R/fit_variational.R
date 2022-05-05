#' Generate initial values from longitudinal cdm model
#'
#' @param variational_dist Variational distributions, including means and vars
#' @param Xdata Data list
#' @param Q Q-matrix
#' @param init_vals Initial values
#' @return Returns value of ELBO function at variational dists
#' @export
cdm_ELBO_logistic=function(vb_vec,Xdata,paraminfo,Q,nsamples=20){
  nparams=length(vb_vec)/2
  prior_means=prior_info$prior_mean
  prior_vars=prior_info$prior_var
  vb_mean=vb_vec[1:nparams]
  vb_var=vb_vec[(nparams+1):(2*nparams)]
  q_info=gen_q_info(Q)

  samp_ELBO=function(vec){
    myvec=rnorm(nparams,mean=vb_mean,sd=sqrt(vb_var))%>%
      set_names(names(vb_mean))
    likelihood_master(convert_vector_to_params(myvec,paraminfo),Xdata,q_info)
  }
  logliks=map_dbl(1:nsamples,samp_ELBO)

  lik_term=sum(logliks)
  variational_term=-sum((vb_mean^2+vb_var)/(vb_var)+log(vb_var))/2
  prior_term=-sum((vb_mean*(vb_mean-2*prior_means)+vb_var)/(prior_vars))/2

  return(-(lik_term+variational_term+prior_term))
}

cdm_ELBO_probit=function(vb_vec,Xdata,paraminfo,Q,prior_info=NULL,nsamples=20){
  nparams=length(vb_vec)/2
  prior_means=prior_info$prior_mean
  prior_vars=prior_info$prior_var
  vb_mean=vb_vec[1:nparams]
  vb_var=vb_vec[(nparams+1):(2*nparams)]
  q_info=gen_q_info(Q)

  Xs=Xdata$Xs
  person_lik=function(r){
    lvals=array(NA,c(Nquestions,Nprofile,Nprofile,Ntime-1))
    for(t in 2:Ntime){
      for(p1 in 1:Nprofile){
        for(p2 in 1:Nprofile){
          lvals[,p1,p2,t-1]=
            log_trans_probs[[t-1]][[r]][p1,p2]+
            dbinom(Xs[[t-1]][r,],1,response_probs[,p1],log=T)+
            dbinom(Xs[[t]][r,],1,response_probs[,p2],log=T)
        }
      }
    }
    retval=log_sum_exp(apply(lvals,c(2,3),sum))
    return(sum(retval))
  }
  retval=sum(map_dbl(1:Nrespondents,person_lik))

  lik_term=sum(logliks)
  variational_term=-sum((vb_mean^2+vb_var)/(vb_var)+log(vb_var))/2
  prior_term=-sum((vb_mean*(vb_mean-2*prior_means)+vb_var)/(prior_var))/2

  return(-(lik_term+variational_term+prior_term))
}


variational_fit=function(initial_params,Xdata,maxiter=100){
  paramvec=convert_params_to_vector(initial_params)
  init_vb_vec=c(convert_params_to_vector(initial_params),rep(.1,length(paramvec)))
  variational_vec=cbind()
  Nparam=length(paramvec)
  priors=list(prior_means=rep(0,Nparam),
              prior_vars=rep(1,Nparam))
  ELBO_wrapper=function(vec){
    cdm_ELBO_logistic(vec,Xdata,initial_params,Q,priors,nsamples=5)
  }
  ELBO_wrapper(init_vb_vec)
  optout=optim(vb_vec,ELBO_wrapper,control=list(maxit=10000,trace=6))
  npar=length(init_vb_vec)/2
  return(list(vb_mean=optout$par[1:npar],vb_var=optout$par[(npar+1):(2*npar)]))
}

variational_fit_probit=function(initial_params,Xdata,maxiter=100){
  paramvec=convert_params_to_vector(initial_params)
  init_vb_vec=c(convert_params_to_vector(initial_params),rep(.1,length(paramvec)))
  #initialize

  variational_vec=cbind()
  Nparam=length(paramvec)
  priors=list(prior_means=rep(0,Nparam),
              prior_vars=rep(1,Nparam))
  ELBO_wrapper=function(vec){
    cdm_ELBO_logistic(vec,Xdata,initial_params,Q,priors,nsamples=5)
  }
  ELBO_wrapper(init_vb_vec)
  optout=optim(vb_vec,ELBO_wrapper,control=list(maxit=10000,trace=6))
  npar=length(init_vb_vec)/2
  return(list(vb_mean=optout$par[1:npar],vb_var=optout$par[(npar+1):(2*npar)]))
}
