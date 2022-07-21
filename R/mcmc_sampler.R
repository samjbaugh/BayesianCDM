#' Run MCMC sampler for longitudinal cdm model
#'
#' @param Xdata Data list
#' @param Q Q-matrix
#' @param initparams Initial Parameters
#' @return Returns dataframe of MCMC samples
#' @export
mcmc_sampler_main=function(Xdata,initparams,M){
  init_paramvec=convert_params_to_vector(initparams)

  Nparams=length(init_paramvec)
  mcmc_samples=matrix(NA,M,Nparams)%>%
    data.frame()%>%
    set_names(names(init_paramvec))
  mcmc_samples[1,]=init_paramvec
  q_info=gen_q_info(Xdata$Q)

  #prior variance for u is a 3x3 matrix since we have two
  #individual level covariates+intercept
  Nu=initparams$Ns$Nrespcov+1
  init_ucovmat=matrix(rWishart(1,Nu,diag(Nu))[,,1],Nu,Nu)
  prior_vars=list('theta'=1,
                    'gamma'=1,
                    'u'=init_ucovmat)

  paraminfo=initparams[c('Ns','log_lambda','value_key','beta_names','gamma_names')]
  lik_wrapper=function(vec){
    likelihood_master(convert_vector_to_params(vec,paraminfo),Xdata,q_info)
  }

  prior_wrapper=function(vec){
    theta_prior=dnorm(vec[startsWith(names(vec),'theta')],mean=0,sd=prior_vars$theta,log=T)
    gamma_prior=dnorm(vec[startsWith(names(vec),'gamma')],mean=0,sd=prior_vars$gamma,log=T)
    #for beta, need to calculate u
    params=convert_vector_to_params(vec,initparams)
    u=gen_u(params,Xdata)
    uf_prior=apply(u$uforward,c(1,3,4),function(x) dmvnorm(x,mean=rep(0,Nu),sigma=prior_vars$u))
    ub_prior=apply(u$ubackward,c(1,3,4),function(x) dmvnorm(x,mean=rep(0,Nu),sigma=prior_vars$u))
    return(sum(theta_prior)+sum(gamma_prior)+sum(uf_prior)+sum(ub_prior))
  }

  proposal_sds=rep(.2,Nparams)
  current_paramvec=init_paramvec
  current_likelihood=lik_wrapper(init_paramvec)
  current_prior=prior_wrapper(init_paramvec)

  likelihoods=rep(NA,M)
  likelihoods[1]=current_likelihood
  accepts=rep(NA,M)
  accepts[1]=T

  for(i in 2:M){
     print(i)
     proposed_paramvec=current_paramvec+rnorm(Nparams,0,sd=proposal_sds)
     proposed_likelihood=lik_wrapper(proposed_paramvec)
     proposed_prior=prior_wrapper(proposed_paramvec)

     mhval=proposed_likelihood-current_likelihood+
       proposed_prior-current_prior
     if(exp(mhval)>runif(1)){
       current_paramvec=proposed_paramvec
       current_likelihood=proposed_likelihood
       current_prior=proposed_prior
       accepts[i]=T
     }else{
       accepts[i]=F
     }
     mcmc_samples[i,]=current_paramvec
  }
  return(mcmc_samples)
}

