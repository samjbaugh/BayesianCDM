#' Generate initial values from longitudinal cdm model
#'
#' @param Nrespondents Input values
#' @param Q Q-matrix
#' @param Ntime Number of time points
#' @param Nrespcov Number of respondent-level covariates
#' @param Ngroupcov Number of group-level covariates
#' @param Ngroup Number of groups
#' @param seed Optional seed
#' @return Returns randomly-generated parameter values
#' @export
gen_initial_values_longitudinal<-function(Nrespondents,Q,Ntime=2,
                                          Nrespcov=0,Ngroupcov=0,
                                          Ngroup=2,seed=NULL){
  if(!is.null(seed)){
    set.seed(seed)
  }
  # if(Ngroup>2){
  #   print("For Ngroup>2 not implemented")
  #   return(NA)
  # }
  Nskill=dim(Q)[2]
  Nquestions=dim(Q)[1]
  Nprofile=2^Nskill

  q_info=gen_q_info(Q)

  intercepts=rnorm(Nquestions,mean=-3,sd=.1)
  Ninteraction=length(q_info$interaction_qids)*
    (q_info$interaction_qids[1]!=0)
  Ninteraction_unique=length(unique(q_info$interaction_qids))*
    (q_info$interaction_qids[1]!=0)
  base_effects=matrix(NA,Nquestions,Nskill)
  base_effects[1:Nquestions,]=rnorm(Nskill*Nquestions,mean=6,sd=.1)
  base_effects[unique(q_info$interaction_qids),]=
    rnorm(Nskill*Ninteraction_unique,mean=3,sd=.1)

  interactions=rnorm(Ninteraction,mean=0,sd=.1)

  gen_forward=function(s) {
    set.seed(123)
    forward_beta=array(rnorm(Ngroup*(Nrespcov+1),sd=.5),c(Ngroup,Nrespcov+1))
    #simulate different group base effects
    forward_beta[,1]=seq(-2,2,length=Ngroup)
    return(forward_beta)
  }
  gen_backward=function(s){
    return(matrix(-2,Ngroup,Nrespcov+1))
  }
  gen_names=function(t,s){
    outer(1:Ngroup,1:(Nrespcov+1),function(g,v) paste0('_g',g,'_v',v-1,'_s',s,'_t',t))
  }
  #first index is over time, second is over skills

  #can recover vector from array(vec,c(Ngroup,Nrespcov+1,Nskill,Ntime-1))
  forward_betas=simplify2array(lapply(1:(Ntime-1),function(t) simplify2array(lapply(1:Nskill,gen_forward))))
  backward_betas=simplify2array(lapply(1:(Ntime-1),function(t) simplify2array(lapply(1:Nskill,gen_backward))))
  beta_names=simplify2array(lapply(1:(Ntime-1),function(t) simplify2array(lapply(1:Nskill,function(s) gen_names(t,s)))))

  #can recover vector from array(vec,c(Nrespcov+1,Ngroupcov+1))
  gamma_mat=matrix(rnorm((Nrespcov+1)*(Ngroupcov+1)),Nrespcov+1,Ngroupcov+1)
  gamma_names=outer(1:(Nrespcov+1),1:(Ngroupcov+1),function(v,w) paste0('gamma_v',v-1,'_w',w-1))

  value_key=c(rep('intercepts',length(intercepts)),
              rep('base_effects',length(c(base_effects))),
              rep('interactions',Ninteraction))

  base_effect_tags=c(outer(1:Nskill,1:Nquestions,function(s,i) paste0('_s',s,'_i',i)))
  int_tags=q_info$interaction_list%>%map('interaction')%>%
    map_chr(function(x) paste0('(',x[1],',',x[2],')'))

  theta_names=c(paste0('theta0_i',1:length(intercepts)),
                paste0('theta',base_effect_tags))
  if(Ninteraction>0){
    theta_names=c(theta_names,paste0('theta_',int_tags))
  }


  init_vals=list(Ns=list(Nrespondents=Nrespondents,Ntime=Ntime,
                           Nquestions=Nquestions,Nrespcov=Nrespcov,
                           Ngroupcov=Ngroupcov,Ngroup=Ngroup,
                           Ninteraction=Ninteraction,
                           Nskill=Nskill,Nprofile=Nprofile),
                 theta=c(intercepts,base_effects,interactions),
                 log_lambda=log(rep(1,Nprofile)/Nprofile),
                 value_key=value_key,
                 beta_names=beta_names,
                 gamma_names=gamma_names,
                 gamma=gamma_mat,
                 forward_betas=forward_betas,
                 backward_betas=backward_betas,
                 theta_names=theta_names)
  return(init_vals)
}
