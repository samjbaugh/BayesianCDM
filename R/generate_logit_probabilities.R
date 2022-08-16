generate_logits_discrete<-function(params,q_info){
  Nquestions=params$Ns$Nquestions
  Nskill=params$Ns$Nskill
  Nprofile=params$Ns$Nprofile
  Nquestions=params$Ns$Nquestions
  value_key=params$value_key
  Ninteraction=params$Ns$Ninteraction

  intercepts=params$theta[value_key=='intercepts']
  base_effects=matrix(0,Nquestions,Nskill)
  base_effects[Q==1]=params$theta[params$value_key=='base_effects']
  interactions=params$theta[params$value_key=='interactions']


  ret_logits=matrix(0,Nquestions,Nprofile)
  for(iquestion in 1:Nquestions){
    for(iprofile in 1:Nprofile){
      myprob =
        intercepts[iquestion];
      for(iskill in 1:Nskill){
        if(q_info$Q[iquestion,iskill]==1 &
           q_info$skill_in_profile[iskill,iprofile]==1){
          myprob=myprob +
            base_effects[iquestion,iskill];
        }
      }
      if(Ninteraction>0){
        for(i_interaction in 1:Ninteraction){
          if(q_info$interaction_qids[i_interaction]==iquestion &
             q_info$interaction_in_profile[i_interaction,iprofile]==1){
            myprob=myprob + interactions[i_interaction];
          }
        }
      }
      ret_logits[iquestion,iprofile]=myprob
    }
  }

  return(ret_logits)
}
