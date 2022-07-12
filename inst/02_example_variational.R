require(tidyverse)
require(cdmfits)
Nrespondents=50
Nquestions=21
#two group assignments
Ngroup=2
Ntime=2
group_assignments=map(1:Ngroup,function(i) rep(i,ceiling(Nrespondents/Ngroup)))%>%
  {do.call(c,.)}%>%.[1:Nrespondents]

myseed=2171506
Nskill=3
Nprofile=2^Nskill
Q=gen_profile_list(Nprofile)%>%
  {do.call(rbind,lapply(.,function(x)
    do.call(rbind,lapply(1:ceiling(Nquestions/(Nprofile)),function(i) x))))%>%
      .[1:Nquestions,]}


Nrespcov=2
Ngroupcov=2
true_params=gen_initial_values_longitudinal(Nrespondents,Q,Ngroup=2,seed=myseed,
                                            Nrespcov=Nrespcov,Ngroupcov=Ngroupcov)

respondent_designmat=cbind(rep(1,Nrespondents),matrix(rnorm(Nrespondents*Ngroup,sd=.2),Nrespondents,Nrespcov))
group_designmat=cbind(rep(1,Ngroup),
                      matrix(rnorm(Ngroup*Ngroupcov,sd=.2),
                             Ngroup,Ngroupcov))

Xdata=simulate_cdm_data(Nrespondents=Nrespondents,
                        Nquestions=Nquestions,
                        Nskill=Nskill,
                        Ntime=Ntime,
                        Ngroup=Ngroup,
                        group_assignments=group_assignments,
                        seed=myseed,
                        true_params=true_params,
                        Q=Q,
                        respondent_designmat=respondent_designmat,
                        group_designmat=group_designmat)$Xdata
initial_params=gen_initial_values_longitudinal(Nrespondents,Q=Q,
                                               Nrespcov=Nrespcov,
                                               Ngroupcov=Ngroupcov,
                                               Ntime=Ntime)
tvar=system.time({variational_out=variational_fit(initial_params,maxiter=100)})

nresamp=1000
variational_resamp=
  data.frame(rbind(variational_out$vb_mean,variational_out$vb_var))%>%
  mutate(type=c('mean','var'))%>%
  pivot_longer(!'type',names_to='d')%>%
  group_by(d)%>%
  nest()%>%
  # rowwise()%>%
  mutate(data=map(data,~data.frame(i=1:nresamp,vals=rnorm(nresamp,.$value[1],.$value[2]))))%>%
  unnest(data)%>%
  pivot_wider(names_from=d,values_from=vals)
# p(variational_resamp,names_from=d)

names(variational_resamp)[1+c(1:4,22:23,30:32+60,94:97,107:110+5,130:133)]
plotdf=variational_resamp[,1+c(1:4,22:23,30:31+60,94:97,107:110+5,130:133)]%>%
  pivot_longer(everything())%>%
  mutate(type=map_chr(name,~strsplit(.,'_')[[1]][1]))

ggplot(plotdf) +
  geom_density(aes(x = value,fill=type)) +
  facet_wrap(~ name)

# M=2000
# initparams=true_params
# require(mvtnorm)
# mcmc_samples=mcmc_sampler_main(Xdata,initparams,M)
#
# names(mcmc_samples)[c(1:4,22:23,30:32+60,94:97,107:110+5,130:133)]
# plotdf=mcmc_samples[1001:2000,c(1:4,22:23,30:31+60,94:97,107:110+5,130:133)]%>%
#   pivot_longer(everything())%>%
#   mutate(type=map_chr(name,~strsplit(.,'_')[[1]][1]))
#
# ggplot(plotdf) +
#   geom_histogram(aes(x = value,fill=type), bins = 15) +
#   facet_wrap(~ name)
#
