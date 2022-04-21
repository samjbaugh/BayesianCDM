#' Generate initial values from longitudinal cdm model
#'
#' @param variational_dist Variational distributions, including means and vars
#' @param Xdata Data list
#' @param Q Q-matrix
#' @param init_vals Initial values
#' @return Returns value of ELBO function at variational dists
#' @export
cdm_ELBO=function(variational_dist,Xdata,Q,init_vals){
  nsamples=5
  value_vec=rep(NA,nsamples)
  q_info=generate_q_info(Q)
  nlambda=length(variational_dist$lambda_mean)
  samp_ELBO=function(x){
    sampled_lambda=rnorm(nlambda,
                         mean=variational_dist$lambda_mean,
                         sd=variational_dist$lambda_sd)
    sample_vals=list('lambda'=sampled_lambda,
                     'log_theta'=init_vals$log_theta,
                     'value_key'=init_vals$value_key)
    return(likelihood_discrete(sample_vals,X,q_info))
  }
  value_vec=map_dbl(1:nsamples,samp_ELBO)
  #compute expectation
  term1=mean(value_vec)
  term2=sum(log(variational_dist$lambda_sd))
  prior_term=-sum(variational_dist$lambda_mean^2)-
    sum(variational_dist$lambda_sd^2)
  return(-(term1+term2+prior_term))
}
