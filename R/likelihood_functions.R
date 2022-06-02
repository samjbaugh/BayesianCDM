#' Likelihood function for generic longitudinal multi-level cdm model
#'
#' @param params parameter list
#' @param Xdata data list
#' @param q_info Result from applying gen_qinfo with Q matrix
#' @return Returns likelihood value
#' @export
likelihood_master<-function(params,Xdata,q_info){
  response_probs=logit(generate_logits_discrete(params,q_info))
  log_trans_probs=gen_trans_probs(params,Xdata,ret_prof_trans=T)$profile
  Nprofile=params$Ns$Nprofile
  Ntime=params$Ns$Ntime
  Nrespondents=params$Ns$Nrespondents
  Nquestions=params$Ns$Nquestions
  #only for two time points for now

  #calculate data likelihood
  Xs=Xdata$Xs
  person_lik=function(r){
    lvals=array(0,c(Nquestions,Nprofile,Nprofile,Ntime-1))
    #first time point (no transitions)
    for(p1 in 1:Nprofile){
      lvals[,p1,p1,t-1]=
        dbinom(Xs[[1]][r,],1,response_probs[,p1],log=T)
    }
    #further time points (include transition probabilities)
    for(t in 2:Ntime){
      for(p1 in 1:Nprofile){
        for(p2 in 1:Nprofile){
          lvals[,p1,p2,t-1]=
            log_trans_probs[[t-1]][[r]][p1,p2]+
            dbinom(Xs[[t-1]][r,],1,response_probs[,p1],log=T)+
            dbinom(Xs[[t]][r,],1,response_probs[,p2],log=T)
        }
      }
    }
    retval=log_sum_exp(apply(lvals,c(2,3),sum))
    return(sum(retval))
  }
  retval=sum(map_dbl(1:Nrespondents,person_lik))

  return(sum(retval))
}
