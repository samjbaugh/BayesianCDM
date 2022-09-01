#' Extracts information from Q-matrix to be used in fitting the model
#'
#' @param Q Q-matrix
#' @return Returns a list containing interaction_in_profile,
#' interaction_qids, interaction_list, skill_in_profile, profile_list,
#' and q_profiles
#' @export
gen_q_info<-function(Q){
  Nq=dim(Q)[1]
  Nskill=dim(Q)[2]
  Nprofile=2^Nskill

  profile_list=list()
  for(ii in 1:(2^Nskill)){
    tmp=rep(NA,Nskill)
    for(jj in 1:Nskill){
      tmp[Nskill-jj+1]=((ii-1)%/%(2^(jj-1)))%%2
    }
    profile_list[[ii]]=rev(tmp)
  }

  skill_in_profile=sapply(1:Nprofile, function(i)
    sapply(1:Nskill, function(j) profile_list[[i]][j]==1))

  which_skill_profile=sapply(1:Nskill, function(s)
    which(skill_in_profile[s,]))


  #list of interaction names for each question
  interaction_list=list()
  counter=0
  for(iquestion in 1:Nq){
    for(iprofile in 1:Nprofile){
      prof=profile_list[[iprofile]]
      if(sum(prof)>1 &
         all(Q[iquestion,][prof==1]==1)){
        counter=counter+1
        interaction_list[[counter]]=list('q'=iquestion,'profn'=iprofile,
                                         interaction=which(prof==1))
      }
    }
  }
  if(length(interaction_list)>0){
    interaction_qids=sapply(interaction_list,function(x) x$q)
  }else{
    interaction_qids=c(0)
  }


  Ninteraction=counter#length(interaction_list)
  interaction_in_profile=array(NA,c(Ninteraction,Nprofile))
  for(iprofile in 1:Nprofile){
    if(Ninteraction>0){
      for(i_interaction in 1:Ninteraction){
        interaction_in_profile[i_interaction,iprofile]=
          all(profile_list[[iprofile]][interaction_list[[i_interaction]]$interaction]==1)
      }
    }
  }
  q_profiles=apply(Q,1,function(y)which(sapply(profile_list,function(x) all(x==y))))

  return(list(Q=Q,
              interaction_in_profile=interaction_in_profile,
              interaction_qids=interaction_qids,
              interaction_list=interaction_list,
              skill_in_profile=skill_in_profile,
              which_skill_profile=which_skill_profile,
              profile_list=profile_list,
              q_profiles=q_profiles))
}
