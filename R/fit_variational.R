#' Probit ELBO for cdm model
#'
#' @export
cdm_ELBO=function(vb_vec,Xdata,param_info,Q,prior_info,g,nsamples=20){
  nparams=length(vb_vec)/2
  prior_means=prior_info$prior_mean
  prior_vars=prior_info$prior_var
  vb_mean=vb_vec[1:nparams]
  vb_var=vb_vec[(nparams+1):(2*nparams)]
  q_info=gen_q_info(Q)

  samp_ELBO=function(vec){
    myvec=rnorm(nparams,mean=vb_mean,sd=sqrt(vb_var))%>%
      set_names(names(vb_mean))
    likelihood_master(convert_vector_to_params(myvec,param_info),Xdata,q_info,g=g)
  }
  logliks=map_dbl(1:nsamples,samp_ELBO)

  lik_term=sum(logliks)
  variational_term=-sum((vb_mean^2+vb_var)/(vb_var)+log(vb_var))/2
  prior_term=-sum((vb_mean*(vb_mean-2*prior_means)+vb_var)/(prior_vars))/2

  return(-(lik_term+variational_term+prior_term))
}

#' Probit variational fit for cdm model
#'
#' @export
variational_fit=function(Xdata,initial_params,g,maxiter=10000){
  #obtain metadata from init vector
  param_info=initial_params
  paramvec=convert_params_to_vector(initial_params)
  init_vb_vec=c(convert_params_to_vector(initial_params),rep(.1,length(paramvec)))
  variational_vec=cbind()
  Nparam=length(paramvec)
  prior_info=list(prior_means=rep(0,Nparam),
                  prior_vars=rep(1,Nparam))
  Q=Xdata$Q
  ELBO_wrapper=function(vec){
    cdm_ELBO(vec,Xdata,param_info,Q,prior_info,g,nsamples=5)
  }
  # ELBO_wrapper(init_vb_vec)
  optout=optim(init_vb_vec,ELBO_wrapper,control=list(maxit=maxiter,trace=6))
  npar=length(init_vb_vec)/2
  return(list(vb_mean=optout$par[1:npar],vb_var=optout$par[(npar+1):(2*npar)]))
}

