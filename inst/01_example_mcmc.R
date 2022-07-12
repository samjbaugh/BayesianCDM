require(tidyverse)
Nrespondents=50
Nquestions=21
#two group assignments
Ngroup=2
Ntime=2

myseed=2171506
Nskill=3
Nprofile=2^Nskill
Nrespcov=2
Ngroupcov=2

simdata=gen_data_wrapper(Nrespondents=Nrespondents,
                    Nquestions=Nquestions,
                    Nskill=Nskill,
                    Ntime=Ntime,
                    Ngroup=Ngroup,
                    Nrespcov,
                    Ngroupcov)
M=2000
require(mvtnorm)
system.time({
  mcmc_samples=mcmc_sampler_main(simdata$Xdata,simdata$true_params,M)
})


#sample analysis of parameters
viz_indices=c(1:4,22:23,30:32+60,94:97,107:110+5,130:133)
names(mcmc_samples)[viz_indices]
plotdf=mcmc_samples[1001:2000,c(1:4,22:23,30:31+60,94:97,107:110+5,130:133)]%>%
  pivot_longer(everything())%>%
  mutate(type=map_chr(name,~strsplit(.,'_')[[1]][1]))

ggplot(plotdf) +
  geom_histogram(aes(x = value,fill=type), bins = 15) +
  facet_wrap(~ name)

