#' Logit function
#' @param params parameter object
#' @param Xdata data (list form)
#' @param g link function
#' @param ret_prof_trans should profiles be returned
#' @export
gen_trans_probs=function(params,Xdata,g=logit,ret_prof_trans=F){
  forward_transitions=lapply(2:params$Ns$Ntime,function(t)
    g(sapply(1:params$Ns$Nskill,function(s)
      apply(params$forward_betas[Xdata$group_assignments,,s,t-1]*
              Xdata$respondent_designmat,1,sum))))
  backward_transitions=lapply(2:params$Ns$Ntime,function(t)
    g(sapply(1:params$Ns$Nskill,function(s) apply(params$backward_betas[Xdata$group_assignments,,s,t-1]*Xdata$respondent_designmat,1,sum))))

  if(ret_prof_trans){
    get_transmat=function(t,r){
      retval=matrix(NA,Xdata$Ns$Nprofile,Xdata$Ns$Nprofile)
      for(p1 in 1:Xdata$Ns$Nprofile){
        for(p2 in 1:Xdata$Ns$Nprofile){
          profile1=profile_list[[p1]]
          profile2=profile_list[[p2]]
          retval[p1,p2]=exp(sum(
            log(1-backward_transitions[[t-1]][r,])*(profile1*profile2) +
              log(forward_transitions[[t-1]][r,])*(1-profile1)*profile2 +
              log(backward_transitions[[t-1]][r,])*(1-profile2)*profile1 +
              log(1-forward_transitions[[t-1]][r,])*(1-profile2)*(1-profile1)))
        }
      }
      return(retval)
    }
    profile_list=gen_q_info(Xdata$Q)$profile_list
    profile_transitions=lapply(2:Xdata$Ns$Ntime,function(t) lapply(1:Xdata$Ns$Nrespondents,function(r) get_transmat(t,r)))
  }else{
    profile_transitions=NA
  }
  return(list(forward=forward_transitions,
              backward=backward_transitions,
              profile=profile_transitions))
}
