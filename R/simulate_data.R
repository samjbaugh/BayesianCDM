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
    Q=do.call(rbind,lapply(profile_list,function(x)
      do.call(rbind,lapply(1:ceiling(Nquestions/(Nprofile)),function(i) x))))%>%
      .[1:Nquestions,]
  }

  tmpdata=list(Ns=list(Nrespondents=Nrespondents,Ntime=Ntime,
                     Nquestions=Nquestions,Nrespcov=Nrespcov,
                     Ngroupcov=Ngroupcov,Ngroup=Ngroup,
                     Nskill=Nskill,Nprofile=Nprofile),
             Q=Q,
             group_designmat=group_designmat,
             respondent_designmat=respondent_designmat)

  if(is.null(true_params)){
    true_params=gen_initial_values_longitudinal(Nrespondents,Q,
                                                Nrespcov=Nrespcov,Ngroupcov=Ngroupcov,
                                                Ngroup=Ngroup,seed=myseed)
  }

  #generate transition matrices
  transition_probabilities=gen_trans_probs(true_params,tmpdata,ret_prof_trans = T)
  if(is.null(true_profiles)){
    #give an even distribution for profiles at time 1
    true_profiles=list(
      c(sapply(1:Nprofile,function(x) rep(x,ceiling(Nrespondents/Nprofile))))[1:Nrespondents]
    )
    #Use transition probability matrix at future times:
    trans_probs=gen_trans_probs(true_params,tmpdata,ret_prof_trans=T)$profile
    if(Ntime>1){
      for(t in 2:Ntime){
        true_profiles[[t]]=map_int(1:Nrespondents,function(r) sample(Nprofile,1,prob=trans_probs[[t-1]][[r]][true_profiles[[t-1]][r],]))
      }
    }
  }

  q_info=gen_q_info(Q)
  interaction_qids=q_info$interaction_qids

  true_probs=logit(generate_logits_discrete(true_params,q_info))
  Xs=map(1:Ntime,function(t) t(sapply(true_profiles[[t]],function(y) sapply(true_probs[,y],function(p) sample(c(0,1),size=1,prob=c(1-p,p))))))
  # Xdata[['Xs']]=Xs
  Xdata=list(Ns=list(Nrespondents=Nrespondents,Ntime=Ntime,
                     Nquestions=Nquestions,Nrespcov=Nrespcov,
                     Ngroupcov=Ngroupcov,Ngroup=Ngroup,
                     Nskill=Nskill,Nprofile=Nprofile),
             Q=Q,
             group_designmat=group_designmat,
             respondent_designmat=respondent_designmat,
             Xs=Xs)
  return(list(Xdata=Xdata,true_profiles=true_profiles,true_params=true_params))
}
