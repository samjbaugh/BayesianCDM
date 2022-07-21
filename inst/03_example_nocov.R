require(tidyverse)
require(mvtnorm)
devtools::load_all()

Xs=readRDS('inst/Xs.rds')
Q=readRDS('inst/Q.rds')
group_assignments=readRDS('inst/group_assignments.rds')

Nrespcov=0
Ngroupcov=0
Nrespondents=dim(Xs[[1]])[1]
Ngroup=length(unique(group_assignments))

group_designmat=matrix(1,Ngroup,1)
respondent_designmat=matrix(1,Nrespondents,1)

input_data=list(Xs=Xs,
                Q=Q,
                group_assignments=group_assignments,
                group_designmat=group_designmat,
                respondent_designmat=respondent_designmat)
initparams=gen_initial_values_longitudinal(Nrespondents,
                                           Q,
                                           Ngroup=Ngroup,
                                           Nrespcov=Nrespcov,
                                           Ngroupcov=Ngroupcov)
input_data$Ns=initparams$Ns
mcmc_sampler_main(input_data,initparams,100)
