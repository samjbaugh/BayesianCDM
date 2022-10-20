#' MCMC Sampler
#'
#' @description Add detail
#' @param Ys List of response data
#' @param Xs Design matrices for gamma
#' @param Qs List of Q matrices
#' @param M Number of MCMC iterations
#' @param priors priors
#' @param initparams Initial parameters (optional)
#' @importFrom magrittr %>%
#' @importFrom stats model.matrix as.formula rnorm
#' @import dplyr purrr truncnorm
#' @importFrom mvtnorm rmvnorm
#' @export

sample_longitudinal=function(Ys,Xs,Qs,M,
                             priors=list(beta_prior=1,gamma_prior=1),
                             initparams=NULL,
                             fixed_beta=F){
    gen_init=is.null(initparams)

    # print(paste('gen_init',gen_init))
    # print('asfdasdf')
    # priors=list(beta_prior=500,gamma_prior=1)
    priorsd_beta=matrix(priors$beta_prior,Nq_total,Nprofile)
    priorsd_gamma=map(Xs,~map(.,function(x) priors$gamma_prior*diag(dim(x)[2])))

    #setup
    {
      Yconcat=do.call(cbind,Ys)
      Nrespondents=dim(Yconcat)[1]
      Nq_list=map(Qs,~dim(.)[1])
      Nq_total=sum(unlist(Nq_list))
      Nskill=dim(Qs[[1]])[2]
      Ntime=length(Ys)
      Nprofile=2^Nskill
      Ntransition=2^Ntime
      Nq=Nq_list[[1]]                    # new
      single_delta=Q_to_delta(Qs[[1]])   # new

      profiles=map_dfr(1:Nprofile-1,
                       ~data.frame(t(rev(as.integer(intToBits(.))[1:Nskill]))))
      A=model.matrix(data=profiles,
                     as.formula(paste0('~', paste0(names(profiles),collapse='*'))))
      transitions=map_dfr(1:(2^Ntime)-1,
                          ~data.frame(t(rev(as.integer(intToBits(.))[1:Ntime]))))

      delta_cat=do.call(rbind,map(Qs,Q_to_delta))

      qt_map=unlist(map(1:length(Nq_list),~rep(.,Nq_list[[.]])))
    }


    gen_init=T
    if(gen_init){
      beta_mat=matrix(rnorm(Nprofile*Nq_total),Nq_total,Nprofile)*
        delta_cat
      beta_mat[,1]=-abs(beta_mat[,1])
      beta_mat[,2]=abs(beta_mat[,2])
      beta_mat[,3]=abs(beta_mat[,3])
      gamma_list=map(Xs,~map(.,function(x) rnorm(dim(x)[2])))
    }else{
      initparams=true_params
      beta_mat=initparams$beta_mat
      gamma_list=initparams$gamma_list
    }
    # gamma_list=gamma_list_true
    beta_vec=c(beta_mat[delta_cat==1])
    # beta_mat=true_params$beta_mat
    # gamma_list=gamma_list_true
    # alpha_mat=true_alpha
    gamma_vec=unlist(gamma_list)

    trans_probs=gamma_to_transprobs(gamma_list,Xs)          # logistic(psi) 
    trans_mat=matrix(NA,Nrespondents,Nskill)                # rho
    for(i in 1:Nrespondents){
      for(j in 1:Nskill){
        trans_mat[i,j]=sample(1:(2^Ntime),1,prob=trans_probs[[j]][i,])
      }
    }
    alpha_mat=transitions_to_alpha(trans_mat,transitions)  
    # print(mean(alpha_mat==true_alpha))
    alpha_vec=c(alpha_mat)

    #MCMC setup
    {
      varnames=c('beta','gamma','alpha')
      sampnames=varnames%>%paste0('_vec')%>%
        sapply(function(x) paste(x,'[',1:length(get(x)),']',sep=''))%>%
        unlist()
      samples=matrix(NA,M,length(sampnames))%>%
        as.data.frame()%>%
        set_names(sampnames)
      samples[1,]=c(beta_vec,unlist(gamma_list),alpha_vec)

      Ljp=-Inf
    }
    
    #################### Same beta over time ####################
    if (fixed_beta)
    {
      for(m in 2:M){
        print(m)
        ltheta=beta_mat%*%t(A)
        theta=logistic(ltheta)
      
        #sample alpha
        trans_probs=gamma_to_transprobs(gamma_list,Xs)
        #Do this to avoid NA's:
        theta[theta==1]=.99
        alpha_mat=sample_alpha_mat(alpha_mat,theta,trans_probs,qt_map,profiles,Ys)
        # print(table(alpha_mat))
        # print(mean(alpha_mat==true_alpha))
      
        #calculate new transitions
        trans_mat=alpha_to_transitions(alpha_mat,profiles)
      
        #sample ystar
        nct = matrix(NA,Nprofile,Ntime)
        for(c in 1:Nprofile){
          for(t in 1:Ntime){
            nct[c,t]=sum(alpha_mat[,t]==c)
          }
        }
        ystar=map(1:Ntime,~matrix(NA,Nq,Nprofile))
        for(t in 1:Ntime){
          for(j in 1:Nq){
            for(c in 1:Nprofile){
              ystar[[t]][j,c]=rpg(h=max(nct[c,t],1),z=ltheta[j,c])
            }
          }
        }
      
        #sample beta
        kappa = map(1:Ntime,~matrix(NA,Nq,Nprofile))
        for(t in 1:Ntime){
          for(j in 1:Nq){
            for(c in 1:Nprofile){
              kappa[[t]][j,c]=sum(Ys[[t]][alpha_mat[,t]==c,j])-nct[c,t]/2
            }
          }
        }
        z=map2(kappa,ystar,~.x/.y)
        condsd=matrix(NA,Nq,Nprofile)
        for(j in 1:Nq){
          for(p in 1:Nprofile){
            omegaj=unlist(map(1:Ntime,~ystar[[.]][j,]))
            Aconcat=do.call(rbind, replicate(Ntime, A, simplify=FALSE)) 
            condsd[j,p]=sqrt(1/(t(Aconcat[,p])%*%(omegaj*Aconcat[,p])+1/(priorsd_beta[j,p]^2)))
          }
        }
        condmean=matrix(NA,Nq,Nprofile)
        for(j in 1:Nq){
          for(p in 1:Nprofile){
            omegaj=unlist(map(1:Ntime,~ystar[[.]][j,]))
            zj=c(unlist(map(z,~.[j,])))
            ztildej=zj-Aconcat[,-p]%*%beta_mat[j,-p]
            condmean[j,p]=condsd[j,p]^2*t(Aconcat[,p])%*%(omegaj*ztildej)
          }
        }
        beta_post = list(mean=condmean,sd=condsd)
        beta_mat=matrix(rtruncnorm(length(beta_post$mean),mean=c(beta_post$mean),
                                   sd=c(beta_post$sd),a=c(Ljp)),
                        dim(beta_post$mean)[1],dim(beta_post$mean)[2])*single_delta
        beta_vec=c(beta_mat[single_delta==1])
        beta_mat=do.call(rbind, replicate(Ntime, beta_mat, simplify=FALSE)) 
        beta_vec = rep(beta_vec, Ntime)
        
        #sample gamma
        gamma_list=sample_gamma(gamma_list,trans_mat,Xs,priorsd_gamma)
        
        samples[m,]=c(beta_vec,unlist(gamma_list),alpha_mat)
      }
      

    } else {
      
      
    #################### Varying beta over time ####################
    for(m in 2:M){
      print(m)
      ltheta=beta_mat%*%t(A)
      theta=logistic(ltheta)

      #sample alpha
      trans_probs=gamma_to_transprobs(gamma_list,Xs)
      #Do this to avoid NA's:
      theta[theta==1]=.99
      alpha_mat=sample_alpha_mat(alpha_mat,theta,trans_probs,qt_map,profiles,Ys)
      # print(table(alpha_mat))
      # print(mean(alpha_mat==true_alpha))

      #calculate new transitions
      trans_mat=alpha_to_transitions(alpha_mat,profiles)

      #sample ystar
      nct = matrix(NA,Nprofile,Ntime)
      for(c in 1:Nprofile){
        for(t in 1:Ntime){
          nct[c,t]=sum(alpha_mat[,t]==c)
        }
      }
      ystar=sample_ystar(ltheta,nct,qt_map)#*delta_cat

      #sample beta
      kappa = matrix(NA,Nq_total,Nprofile)
      for(j in 1:Nq_total){
        for(c in 1:Nprofile){
          t=qt_map[j]
          kappa[j,c]=sum(Yconcat[alpha_mat[,t]==c,j])-nct[c,t]/2
        }
      }
      # kappa=sweep(njc,2,nc/2,'-') #t(apply(nmat,1,function(x) x-sum(x)/2))
      z=kappa/ystar

      # if(F){
      beta_post=beta_gibbs_dist(ystar,A,z,beta_mat,priorsd_beta)
      beta_mat=matrix(rtruncnorm(length(beta_mat),mean=c(beta_post$mean),sd=c(beta_post$sd),a=c(Ljp)),
                      dim(beta_mat)[1],dim(beta_mat)[2])*delta_cat
      beta_vec=c(beta_mat[delta_cat==1])


      #sample gamma
      gamma_list=sample_gamma(gamma_list,trans_mat,Xs,priorsd_gamma)

      samples[m,]=c(beta_vec,unlist(gamma_list),alpha_mat)
    }
    }

    gamma_post=sample_gamma(gamma_list,trans_mat,Xs,priorsd_gamma,retmean=T)

    retval=list(samples=samples,
                beta_post=beta_vec,
                gamma_post=gamma_post,
                alpha_mat=alpha_mat)

  return(retval)
}
