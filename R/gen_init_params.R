#' Generate initial values from longitudinal cdm model
#'
#' @param Nrespondents Input values
#' @param Q Q-matrix
#' @param Ntime Number of time points
#' @param Nrespcov Number of respondent-level covariates
#' @param Ngroupcov Number of group-level covariates
#' @param Ngroup Number of groups
#' @param seed Optional seed
#' @param onelvl Whether model only consider one level
#' @return Returns randomly-generated parameter values
#' @export
gen_initial_values_longitudinal<-function(Nrespondents,Q,Ntime=2,
                                          Nrespcov=0,Ngroupcov=0,
                                          Ngroup=2,seed=NULL){
  if(!is.null(seed)){
    set.seed(seed)
  }
  Nskill=dim(Q)[2]
  Nquestions=dim(Q)[1]
  Nprofile=2^Nskill

  q_info=gen_q_info(Q)

  alpha=gen_alpha(Nprofile,Nskill)
  A=gen_A(alpha)

  Delta=Q_to_delta(Q)
  beta_mat=matrix(NA,Nquestions,Nprofile)*Delta
  # set intercepts
  beta_mat[,1]=rnorm(Nquestions,mean=-2,sd=.5)
  # set base effects
  beta_mat[,2:(Nskill+1)]=rnorm(Nquestions*Nskill,mean=2,sd=.5)
  # set interactions
  beta_mat[,(Nskill+2):ncol(beta_mat)]=rnorm(Nquestions*length((Nskill+2):ncol(beta_mat)),mean=0.5,sd=.5)

  beta_mat=beta_mat*Delta
  beta_vec=beta_mat[Delta==1]

  intercepts   = beta_mat[,1]
  base_effects = beta_mat[,2:(Nskill+1)]
  interactions = beta_mat[,(Nskill+2):ncol(beta_mat)]

  psi=beta_mat%*%t(A)
  theta=logistic(psi)

  value_key=c(rep('intercepts',length(beta_mat[,1])),
              rep('base_effects',length(c(beta_mat[,2:(Nskill+1)]))),
              rep('interactions',length(c(beta_mat[,(Nskill+2):ncol(beta_mat)]))))

  base_effect_tags=c(outer(1:Nquestions,1:Nskill,function(j,s) paste0('_j',j,'_s',s)))
  int_tags = c(outer(1:Nquestions,1:ncol(interactions),function(j,int) paste0('_j',j,'_int',colnames(interactions[,int]))))
  beta_names=c(paste0('beta0_j',1:length(intercepts)),
               paste0('beta',base_effect_tags),
               paste0('beta_',int_tags))
  theta_names=c(outer(1:Nquestions,1:Nprofile,function(j,c) paste0('_j',j,'_c',c)))

  params_keep = c(Delta==1) ###### nonzero params kept


  init_vals=list(Ns=list(Nrespondents=Nrespondents,Ntime=Ntime,
                           Nquestions=Nquestions,Nrespcov=Nrespcov,
                           Ngroupcov=Ngroupcov,Ngroup=Ngroup,
                           Nskill=Nskill,Nprofile=Nprofile),
                 beta=unname(c(intercepts,base_effects,interactions))[params_keep],
                 log_lambda=log(rep(1,Nprofile)/Nprofile),
                 value_key=value_key[params_keep],
                 beta_names=beta_names,
                 theta=theta,
                 theta_names=theta_names)
  return(init_vals)
}
