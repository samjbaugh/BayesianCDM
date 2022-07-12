require(tidyverse)
devtools::load_all()

###small data version
Ngroup=2
Ntime=2
myseed=2171506
Nskill=2
Nrespcov=2
Ngroupcov=2

Nrespondents=10
Nquestions=4



time_list=list(mcmc=rep(NA,ntrials),
               vs=rep(NA,ntrials))

M=10000
methods=c('mcmc','varprobit_simple')
fitfun=list(mcmc=function(x,y) mcmc_sampler_main(x,y,M),
            varprobit_simple=function(x,y) variational_fit(x,y,g=pnorm))
ntrials=4
nrs=round(seq(10,40,length=ntrials))
nqs=round(seq(4,20,length=ntrials))
ntotal=nrs*nqs*Ntime
ntrials=length(ntotal)
time_list[['n']]=ntotal
for(itrial in 1:ntrials){
  for(method in methods){
    Nrespondents=nrs[itrial]
    Nquestions=nqs[itrial]
    simdata=gen_data_wrapper(Nrespondents=Nrespondents,
                             Nquestions=Nquestions,
                             Nskill=Nskill,
                             Ntime=Ntime,
                             Ngroup=Ngroup,
                             Nrespcov,
                             Ngroupcov)
    t=system.time({
      results=fitfun[[method]](simdata$Xdata,simdata$true_params)
    })
    time_list[[method]][itrial]=t[1]
  }
}

bind_rows(time_list)%>%
  dplyr::select('mcmc_10000','varprobit_simple','n')%>%
  pivot_longer(!n,values_to='time',names_to='method')%>%
  rename(ntotal='n')%>%
  ggplot()+
  geom_point(aes(x=ntotal,y=time,col=method))
do.call(rbind,time_list)




Xdata=simdata$Xdata
Q=simdata$Xdata$Q
q_info=gen_q_info(Q)
likelihood_master(simdata$true_params,
                  simdata$Xdata,q_info)

plot_mcmc_posteriors=function(mcmc_samples,
                              plot_indices=1:20,
                              burnin=1000){
  plotdf=mcmc_samples[burnin:dim(mcmc_samples)[1],plot_indices]%>%
    pivot_longer(everything())%>%
    mutate(type=map_chr(name,~strsplit(.,'_')[[1]][1]))

  ggplot(plotdf) +
    geom_histogram(aes(x = value,fill=type), bins = 15) +
    facet_wrap(~ name)
}

plot_variational_posteriors=function(mcmc_samples){
  plotdf=mcmc_samples[1001:2000,1:20]%>%
    pivot_longer(everything())%>%
    mutate(type=map_chr(name,~strsplit(.,'_')[[1]][1]))

  ggplot(plotdf) +
    geom_histogram(aes(x = value,fill=type), bins = 15) +
    facet_wrap(~ name)
}

