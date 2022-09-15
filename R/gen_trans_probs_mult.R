#' Logit function
#' @param params parameter object
#' @param Xdata data (list form)
#' @param g link function
#' @param ret_prof_trans should profiles be returned
#' @export
gen_trans_probs_mult=function(params,Xdata,ret_prof_trans=F){
  transitions01 = sapply(1:Nskill, function(s) params$true_probs[[s]][,1])
  transitions10 = sapply(1:Nskill, function(s) params$true_probs[[s]][,2])
  transitions11 = sapply(1:Nskill, function(s) params$true_probs[[s]][,3])
  transitions00 = sapply(1:Nskill, function(s) params$true_probs[[s]][,4])

  if(ret_prof_trans){
    get_transmat=function(r){
      retval=matrix(NA,Xdata$Ns$Nprofile,Xdata$Ns$Nprofile)
      for(p1 in 1:Xdata$Ns$Nprofile){
        for(p2 in 1:Xdata$Ns$Nprofile){
          profile1=profile_list[[p1]]
          profile2=profile_list[[p2]]
          retval[p1,p2]=prod(colSums(matrix(c(
            (transitions11[r,])*(profile1*profile2),
            (transitions01[r,])*(1-profile1)*profile2,
            (transitions10[r,])*(1-profile2)*profile1,
            (transitions00[r,])*(1-profile2)*(1-profile1)), byrow = T, ncol = Nskill)))
        }
      }
      return(retval)
    }
    profile_list=gen_q_info(Xdata$Q)$profile_list
    profile_transitions=lapply(1:Xdata$Ns$Nrespondents,function(r) get_transmat(r))
  }else{
    profile_transitions=NA
  }
  return(list(profile=profile_transitions))
}
