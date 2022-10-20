#gen data
rm(list = ls())
setwd("~/Desktop/BayesianCDM-package_branch/")
devtools::load_all()

require(mvtnorm)
require(truncnorm)
require(tidyverse)
devtools::load_all()


M = 1000
Nrespondents=200
Nskill=3
Ntime=2
k = 3
# Qs=map(Nq_list,~profiles[sample(2^Nskill,.,replace=T),])
Xbase=cbind(rep(1,Nrespondents),c(rep(0,Nrespondents/2),rep(1,Nrespondents/2)))
Xs=build_Xlist_01(Xbase)

Qs=list(matrix(c(rep(c(0,0,0),7),
                   rep(c(0,0,1),7),
                   rep(c(0,1,0),7),
                   rep(c(1,0,0),7),
                   rep(c(1,1,0),7),
                   rep(c(1,0,1),7),
                   rep(c(0,1,1),7),
                   rep(c(1,1,1),7)), 56, 3, byrow=T),
        matrix(c(rep(c(0,0,0),7),
                 rep(c(0,0,1),7),
                 rep(c(0,1,0),7),
                 rep(c(1,0,0),7),
                 rep(c(1,1,0),7),
                 rep(c(1,0,1),7),
                 rep(c(0,1,1),7),
                 rep(c(1,1,1),7)), 56, 3, byrow=T)
)


Nq_list=map(Qs,~dim(.)[1])

# gamma_list=list(
#
# )

delta_cat=do.call(rbind,map(Qs,Q_to_delta))
true_params=generate_params(Nrespondents,Qs,Xs)
Nq_total=sum(unlist(Nq_list))
true_params$beta_mat=
 cbind(rnorm(Nq_total,-1.7,0.4),rnorm(Nq_total,2,0.4),rnorm(Nq_total,2,0.4),rnorm(Nq_total,2,0.4),
       rnorm(Nq_total,-.7,0.4),rnorm(Nq_total,-.7,0.4),rnorm(Nq_total,-.7,0.4),rnorm(Nq_total,-.7,0.4))*
 delta_cat

# a=c(.7,.05,.05,.2)
# b=c(.2,.55,.05,.2)
# g1=log(a[2:4]/(1-sum(a[2:4])))
# g2=log(b[2:4]/(1-sum(b[2:4])))-g1
# gamma_list_true=map(1:Nskill,~
#   rbind(g1,g2),
# )
Xtmp=build_Xlist_01(rbind(c(1,0),c(1,1)))
gamma_list_true=
  list(list(c(-2,3),c(-1.5),c(-1.5)),
       list(c(-2,3),c(-1.5),c(-1.5)),
       list(c(-2,3),c(-1.5),c(-1.5)))
gamma_list=gamma_list_true
true_params$gamma_list=gamma_list

gamma_to_transprobs(gamma_list,Xtmp)

Ntime=2
Nprofile=2^Nskill
group_is=Xs[[1]][[1]][,2]+1
alpha_sampdf=sample_alpha_from_gamma(Nrespondents,Ntime,gamma_list,Xs)%>%
  data.frame()%>%
  set_names(paste0('a',1:Ntime))%>%
  mutate(group=group_is)%>%
  relocate(group,.before=a1)%>%
  group_by(group)%>%
  arrange(group,a1,a2)
true_alpha=alpha_sampdf%>%
  as.matrix()%>%
  .[,-1]

if(F){
  trans_mat=alpha_to_transitions(true_alpha,profiles)
  priorsd_gamma=map(Xs,~map(.,function(x) priors$gamma_prior*diag(dim(x)[2])))
  gamma_out=sample_gamma(gamma_list,trans_mat,Xs,priorsd_gamma)
  
  gamma_to_transprobs(gamma_list,Xtmp)
  gamma_to_transprobs(gamma_out,Xtmp)
}

priors=list(beta_prior=500,gamma_prior=.1)

Ys=generate_data(Nrespondents,Nq_list,true_params$beta_mat,true_alpha)
sampler_out=sample_longitudinal(Ys,Xs,Qs,M,
                                initparams = NULL,
                                priors=list(beta_prior=500,gamma_prior=1))


################################################################################



{
  delta_cat=do.call(rbind,map(Qs,Q_to_delta))
  beta_mat=sampler_out$samples%>%
    dplyr::select(contains('beta'))%>%
    apply(2,mean)
  beta_sd=sampler_out$samples%>%
    dplyr::select(contains('beta'))%>%
    apply(2,sd)
  pbeta=data.frame(postmean=beta_mat,
                   betatrue=true_params$beta_mat[delta_cat==1],
                   betasd=beta_sd)%>%
    filter(postmean<8)%>%
    mutate(i=1:n())%>%
    ggplot()+
    geom_point(aes(x=i,y=postmean,col='postmean'))+
    # geom_errorbar(aes(x=i,ymin=postmean-2*betasd,
    #                   ymax=postmean+2*betasd,col='postmean'))+
    geom_point(aes(x=i,y=betatrue,col='betatrue')) + 
    coord_cartesian(ylim = c(-5, 5), xlim = c(0,224)) 
  pbeta
}

{
  
  alpha_post=sampler_out$samples%>%
    dplyr::select(starts_with('alpha'))
  # M=20
  # data.frame(a=as.numeric(alpha_post[M,]),ii=1:(Nrespondents*Ntime))%>%
  #   mutate(i=ti_map[ii],t=it_map[ii])%>%
  #   mutate(true_alpha=c(true_alpha))%>%
  #   ggplot()+geom_point(aes(x=i,y=a,col=factor(true_alpha)))+
  #   facet_grid(~t)
  
  profiles=map_dfr(1:Nprofile-1,
                   ~data.frame(t(rev(as.integer(intToBits(.))[1:Nskill]))))
  M=dim(alpha_post)[1]
  it_map=unlist(map(1:Ntime,~rep(.,Nrespondents)))
  ti_map=unlist(map(1:Ntime,~1:Nrespondents))
  post_probs=alpha_post%>%
    pivot_longer(everything())%>%
    mutate(name=factor(name,levels=paste0('alpha_vec[',1:dim(alpha_post)[2],']')),
           value=factor(value,levels=1:Nprofile,labels=paste0('prof',1:Nprofile)))%>%
    group_by(name,value)%>%
    summarise(prop=n()/M)%>%
    complete(value,fill=list(prop=0))%>%mutate(ii=as.numeric(map(strsplit(as.character(name),'\\[|\\]'),2)),
                                               i=ti_map[ii],t=it_map[ii])%>%
    rowwise()%>%
    mutate(true_prof=true_alpha[i,t])
  true_probs=post_probs%>%
    group_by(name,true_prof,i,t)%>%
    summarise(prop_correct=prop[true_prof[1]])
  prof_labels=
    apply(profiles,1,function(x) paste(x,collapse=''))
  palpha=true_probs%>%
    mutate(t=factor(t,levels=1:Ntime,labels=paste0('time',1:Ntime)))%>%
    ggplot()+
    geom_point(aes(x=i,y=prop_correct,color=factor(true_prof)))+
    facet_grid(~t)+
    xlab('Person ID')+
    scale_color_discrete(name='profile',labels=prof_labels)
  palpha
}
# scale_y_log10()

{
  gamma_postmean=sampler_out$samples%>%
    dplyr::select(contains('gamma'))%>%
    apply(2,mean)
  gamma_postsd=sampler_out$samples%>%
    dplyr::select(contains('gamma'))%>%
    apply(2,sd)
  
  trans_mat=alpha_to_transitions(true_alpha,profiles)
  priorsd_gamma=map(Xs,~map(.,function(x) priors$gamma_prior*diag(dim(x)[2])))
  tmp=sample_gamma(gamma_list,trans_mat,Xs,priorsd_gamma,retmean=T)
  plotdf=data.frame(gamma_true=unlist(gamma_list_true),
                    postmean=gamma_postmean,
                    postsd=gamma_postsd,
                    pmean=unlist(tmp$gamma_mean),
                    psd=sqrt(unlist(map(tmp$gamma_var,~map(.,diag)))),
                    i=1:length(gamma_postmean))%>%
    mutate(skill=c(rep(1,4),rep(2,4),rep(3,4)),i=c(1:4,1:4,1:4))
  pgamma=ggplot(plotdf)+
    geom_point(aes(x=i-.25,y=postmean,col='posterior mean'))+
    geom_errorbar(aes(x=i-.25,ymin=postmean-postsd*2,ymax=postmean+postsd*2,col='posterior mean'))+
    geom_point(aes(x=i,y=gamma_true,col='true'))+
    # geom_point(aes(x=i+.25,y=pmean,col='posterior mean'))+
    # geom_errorbar(aes(x=i+.25,ymin=pmean-psd*2,ymax=pmean+psd*2,col='ss'))+
    facet_grid(~skill)
  pgamma
}

pbeta
palpha
pgamma
ggsave(pbeta,file='../Figures/beta_recover.png',height=5,width=10)
ggsave(palpha,file='../Figures/alpha_recover.png',height=5,width=10)
ggsave(pgamma,file='../Figures/gamma_recover.png',height=5,width=10)




if(F){
  A=model.matrix(data=profiles,
                 as.formula(paste0('~', paste0(names(profiles),collapse='*'))))
  
  gamma_list=list(gamma_list_true[[1]][1,,drop=F],
                  gamma_list_true[[2]][1,,drop=F])
  Xs=list(matrix(rep(1,Nrespondents),Nrespondents,1),
          matrix(rep(1,Nrespondents),Nrespondents,1))
  gamma_to_probs(gamma_list,Xtmp)
  true_alpha=sample_alpha_from_gamma(Nrespondents,Ntime,gamma_list,Xs)
  true_trans=alpha_to_transitions(true_alpha,profiles)
  trans_mat=alpha_to_transitions(true_alpha,profiles)
  
  eta=map2(gamma_list,Xs,~(.y%*%.x)%>%
             {cbind(rep(0,Nrespondents),.)}%>%
             apply(1,function(x) x-log(sum(exp(x))-exp(x)))%>%
             t())
  tstar=map(1:Nskill,~matrix(NA,Nrespondents,2^Ntime))
  for(i in 1:Nskill){
    for(j in 1:Ntransiton){
      tstar[[i]][,j]=map_dbl(eta[[i]][,j],~rpg(num=1,h=1,z=.))
    }
  }
  
  kappa=map(1:Nskill,~matrix(NA,Nrespondents,Ntransition))
  for(i in 1:Nskill){
    for(j in 1:Ntransition){
      kappa[[i]][,j]=(trans_mat[,i]==j)-1/2
    }
  }
  
  C=map2(gamma_list,Xs,~(.y%*%.x)%>%
           {cbind(rep(0,Nrespondents),.)}%>%
           apply(1,function(x) log(sum(exp(x))-exp(x)))%>%
           t()
  )
  gamma_var=map(1:Nskill,list)
  priorsd_gamma=map(Xs,~map(.,function(x) priors$gamma_prior*diag(dim(x)[2])))
  
  gamma_out=gamma_list
  for(i in 1:Nskill){
    for(j in 2:Ntransition){
      priorsd_gamma[[i]]=10000
      gibbsvar_gamma=solve(crossprod(Xs[[i]]*tstar[[i]][,j,drop=F],Xs[[i]])+
                             solve(priorsd_gamma[[i]]))
      gibbsmean_gamma=gibbsvar_gamma%*%(t(Xs[[i]])%*%(kappa[[i]][,j]+tstar[[i]][,j,drop=F]*C[[i]][,j]))
      
      (t(Xs[[i]])%*%(kappa[[i]][,j]+tstar[[i]][,j,drop=F]*C[[i]][,j]))
      sum(kappa[[i]][,j]+tstar[[i]][,j,drop=F]*C[[i]][,j])
      
      sum(tstar[[i]][,j,drop=F]*C[[i]][,j])/sum(tstar[[i]][,j])
      
      gamma_out[[i]][,j-1]=sum(kappa[[i]][,j])/sum(tstar[[i]][,j])+
        mean(C[[i]][,j])
    }
  }
  
  gamma_list[[1]]
  gamma_out[[1]]
  Xtmp=list(matrix(c(1,1),2,1),
            matrix(c(1,1),2,1))
  gamma_to_probs(gamma_list,Xtmp)[[1]][1,]
  gamma_to_probs(gamma_out,Xtmp)[[1]][1,]
  rbind(prop.table(table(true_trans[,1])))
  
}

