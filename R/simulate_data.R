#' Simulate CDM Xdata
#'
#'@param Nskill Number of skills to simulate
#'@param seed Fix seed for reproducibility
#'@return Returns simulated X and Q-matrix
#'@export
simulate_cdm_data<-function(Nrespondents=20,
                            Nquestions=22,
                            Nskill=3,
                            Ntime=1,
                            Ngroup=1,
                            respondent_designmat=Matrix(rep(1,Nrespondents),Nrespondents,1),
                            group_designmat=Matrix(rep(1,Ngroup),Ngroup,1),
                            group_assignments=NULL,
                            seed=NULL,
                            true_params=NULL,
                            true_skill_transprob=NULL,
                            Q=NULL,
                            true_profiles=NULL){

  if(!is.null(seed)){
    set.seed(seed)
  }

  Nprofile=2^Nskill

  if(is.null(group_assignments)){
    #if not specified assume only one group
    group_assignments=c(1,Nrespondents)
  }
  Ngroup=length(unique(group_assignments))
  if(is.null(Q)){
    profile_list=gen_profile_list(Nprofile)
    #Some profiles may be excluded if Nquestions/Nprofile is not an integer
    Q=do.call(rbind,lapply(profile_list[-1],function(x)
      do.call(rbind,lapply(1:ceiling(Nquestions/(Nprofile-1)),function(i) x))))%>%
      .[1:Nquestions,]
  }

  tmpdata=list(Ns=list(Nrespondents=Nrespondents,Ntime=Ntime,
                     Nquestions=Nquestions,Nrespcov=Nrespcov,
                     Ngroupcov=Ngroupcov,Ngroup=Ngroup,
                     Nskill=Nskill,Nprofile=Nprofile),
             Q=Q,
             group_designmat=group_designmat,
             respondent_designmat=respondent_designmat,
             group_assignments=group_assignments)

  if(is.null(true_params)){
    true_params=gen_initial_values_longitudinal(Nrespondents,Q,
                                                Nrespcov=Nrespcov,Ngroupcov=Ngroupcov,
                                                Ngroup=Ngroup,seed=myseed)
  }

  #generate transition matrices
  transition_probabilities=gen_trans_probs_mult(true_params,tmpdata,
                                                ret_prof_trans = T)
  if(is.null(true_profiles)){
    t1_profiles = c()
    t2_profiles = c()
    for (r in 1:Nrespondents){
      prof_sample = sample(1:Nprofile^2, 1, prob = c(transition_probabilities$profile[[r]]))
      t1_profiles = c(t1_profiles, ((prof_sample-1) %% 8) + 1)  # row
      t2_profiles = c(t2_profiles, ((prof_sample-1) %/% 8) + 1) # col
    }
    true_profiles = list(t1_profiles, t2_profiles)
  }

  q_info=gen_q_info(Q)
  interaction_qids=q_info$interaction_qids

  true_probs=logit(generate_logits_discrete(true_params,q_info))
  Xs=map(1:Ntime,
         function(t) t(sapply(true_profiles[[t]],
         function(y) sapply(true_probs[,y],
         function(p) sample(c(0,1),size=1,prob=c(1-p,p))))))
  # Xdata[['Xs']]=Xs
  Xdata=list(Ns=list(Nrespondents=Nrespondents,Ntime=Ntime,
                     Nquestions=Nquestions,Nrespcov=Nrespcov,
                     Ngroupcov=Ngroupcov,Ngroup=Ngroup,
                     Ninteraction=true_params$Ns$Ninteraction,
                     Nskill=Nskill,Nprofile=Nprofile),
             Q=Q,
             group_designmat=group_designmat,
             group_assignments=group_assignments,
             respondent_designmat=respondent_designmat,
             Xs=Xs)
  return(list(Xdata=Xdata,true_profiles=true_profiles,true_params=true_params))
}

#' Wrapper for quick data generation
#'
#'@param Nskill Number of skills to simulate
#'@param seed Fix seed for reproducibility
#'@return Returns simulated X and Q-matrix
#'@export
gen_data_wrapper=function(Nrespondents=20,
                          Nquestions=22,
                          Nskill=3,
                          Ntime=1,
                          Ngroup=1,
                          Nrespcov=2,
                          Ngroupcov=2,
                          myseed=NULL,
                          multinomial = F){
  Nprofile=2^Nskill
  group_assignments=map(1:Ngroup,function(i)
      rep(i,ceiling(Nrespondents/Ngroup)))%>%{do.call(c,.)}%>%.[1:Nrespondents]
  Q=gen_profile_list(Nprofile)%>%.[-1]%>%
    {do.call(rbind,lapply(.,function(x)
      do.call(rbind,lapply(1:ceiling(Nquestions/(Nprofile-1)),function(i) x))))%>%
        .[1:Nquestions,]}
  # respondent_designmat=cbind(rep(1,Nrespondents),
  #                            matrix(rnorm(Nrespondents*Ngroup,sd=.2),
  #                                   Nrespondents,Nrespcov))
  respondent_designmat = cbind(rep(1,Nrespondents), rep(1,Nrespondents))              ######### no cov
  # group_designmat=cbind(rep(1,Ngroup),
  #                       matrix(rnorm(Ngroup*Ngroupcov,sd=.2),
  #                              Ngroup,Ngroupcov))
  group_designmat = cbind(rep(1,Ngroup))              ######### no cov

  true_params=gen_initial_values_longitudinal(Nrespondents,Q,
                                              Ngroup=Ngroup,seed=myseed,
                                              Nrespcov=Nrespcov,
                                              Ngroupcov=Ngroupcov)

  if (multinomial == TRUE)
  {
    true_beta = lapply(1:Nskill, function(s) map(desmat[[1]][1:(J-1)],
                       ~rnorm(dim(.)[2]))%>%set_names(paste0('beta',tunits[1:(J-1)])))
    true_logits = lapply(1:Nskill, function(s) bind_cols(map2(desmat[[1]][1:(J-1)],
                         true_beta[[s]],~as.numeric(.x%*%.y)))%>%cbind(rep(0,Nrespondents)))
    true_probs = lapply(1:Nskill, function(s) t(apply(true_logits[[s]],1,function(x)
                        exp(x-max(x))/sum(exp(x-max(x))))))
    true_params$true_beta = true_beta
    true_params$true_logits = true_logits
    true_params$true_probs = true_probs
  }

  simulate_cdm_data(Nrespondents=Nrespondents,
                    Nquestions=Nquestions,
                    Nskill=Nskill,
                    Ntime=Ntime,
                    Ngroup=Ngroup,
                    group_assignments=group_assignments,
                    seed=myseed,
                    true_params=true_params,
                    Q=Q,
                    respondent_designmat=respondent_designmat,
                    group_designmat=group_designmat)
}
