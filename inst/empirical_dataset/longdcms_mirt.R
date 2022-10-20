## Author: Matthew Madison
## Date: 4/22/2022
## Purpose: long dcms using mirt package

library(dplyr)
library(purrr)
setwd("~/Desktop/CDM/BayesianCDM")
devtools::load_all()
setwd("~/Desktop/CDM/BayesianCDM/single_group")

#Load packages
library(mirt) #for DCM estimation with constraints


##Step 1: read in pre - post data, 21 items, 4 attributes
dat = read.table("prepost1.txt")[,-1]
colnames(dat) = c(paste("Pre",1:21),paste("Post",1:21)) #gives items names for mirt
dat[dat == 99] = NA #replace 99 with na

###### Remove missing
complete.dat = dat[complete.cases(dat),]


#########################################################################################################################################################


##Step 2: specify the Q-matrix
#a1 is the "intercept"
#a2-a5 are attributes 1-4 at pre
#a6-a9 are attributes 1-4 at post
#constrain portion specifies measurement invariance
model = mirt.model('a1 = 1-42
a2 = 1,13-14,17
a3 = 2-3,9-12
a4 = 4-8
a5 = 15-16,18-21
a6 = 22,34-35,38
a7 = 23-24,30-33
a8 = 25-29
a9 = 36-37,39-42
CONSTRAIN = (1,22,a2,a6),(1,22,a1),
(2,23,a3,a7),(2,23,a1),
(3,24,a3,a7),(3,24,a1),
(4,25,a4,a8),(4,25,a1),
(5,26,a4,a8),(5,26,a1),
(6,27,a4,a8),(6,27,a1),
(7,28,a4,a8),(7,28,a1),
(8,29,a4,a8),(8,29,a1),
(9,30,a3,a7),(9,30,a1),
(10,31,a3,a7),(10,31,a1),
(11,32,a3,a7),(11,32,a1),
(12,33,a3,a7),(12,33,a1),
(13,34,a2,a6),(13,34,a1),
(14,35,a2,a6),(14,35,a1),
(15,36,a5,a9),(15,36,a1),
(16,37,a5,a9),(16,37,a1),
(17,38,a2,a6),(17,38,a1),
(18,39,a5,a9),(18,39,a1),
(19,40,a5,a9),(19,40,a1),
(20,41,a5,a9),(20,41,a1),
(21,42,a5,a9),(21,42,a1)') #constrain pre-post

#Specify the profiles (intercept in col1, proficiency status in col2)
theta <- cbind(1, thetaComb(0:1, 8))

#calibrate the model
m1 <- mdirt(dat, model, customTheta = theta)
coef(m1) #examine item parameter estimates

#assess model fit, relative and absolute
itemfit(m1, na.rm = T)
extract.mirt(m1, "AIC")
extract.mirt(m1, "BIC")
extract.mirt(m1, "SABIC")
extract.mirt(m1, "logLik")
extract.mirt(m1, "df")
extract.mirt(m1, "nest")

#summarize classification results
posteriors1 = as.matrix(fscores(m1)[,2:9]) #obtain posterior probs of proficiency
#growth in attribute mastery
colMeans(posteriors1)[5]-colMeans(posteriors1)[1]
colMeans(posteriors1)[6]-colMeans(posteriors1)[2]
colMeans(posteriors1)[7]-colMeans(posteriors1)[3]
colMeans(posteriors1)[8]-colMeans(posteriors1)[4]

###########################################################
###########################################################

##Assess measurement invariance
##Fit freely estimated model

set.seed(100)

#Step 2: specify the Q-matrix
#measurement invariance not assumed
model2 = mirt.model('a1 = 1-42
a2 = 1,13-14,17
a3 = 2-3,9-12
a4 = 4-8
a5 = 15-16,18-21
a6 = 22,34-35,38
a7 = 23-24,30-33
a8 = 25-29
a9 = 36-37,39-42')

#Specify the profiles (intercept in col1, proficiency status in col2)
theta <- cbind(1, thetaComb(0:1, 8))

#calibrate the model
system.time({
  mirt_results <- mdirt(dat, model2, customTheta = theta)
})
mirt_results = readRDS("~/Desktop/CDM/BayesianCDM/inst/mirt_results.Rds")

# user  system elapsed
# 91.348   6.219  98.739

mirtcoef = coef(mirt_results)
mirtcoef[[43]] = NULL
mirtcoef = matrix(unlist(mirtcoef), nrow = 42, byrow = T)[22:42,c(1,6:9)]
colnames(mirtcoef) = c("intercept", "s1", "s2", "s3", "s4")
mirtcoef_vec = c(mirtcoef)[c(mirtcoef)!=0]

posteriors = as.matrix(fscores(mirt_results)[,2:9])

plot(mirtcoef_vec[1:21], ylim = c(-9,8))
points(mirtcoef_vec[22:42], pch = 19)




#compare fit to MI model
anova(m1,m2)

#see Madison and Bradshaw for identical results from Mplus



#######################################################################################################################




################# MCMC and Variation CDM #################
set.seed(100)

complete.dat.mat = as.matrix(complete.dat)
Ntime        = 2
Ngroup       = 1
Nrespondents = nrow(complete.dat.mat)
Nrespcov     = 0
Ngroupcov    = 0
Nskill       = 4

Q = matrix(0, nrow = ncol(complete.dat.mat)/2, ncol = 4)
Q[c(1,13:14,17), 1] = Q[c(2:3,9:12), 2] =
Q[c(4:8), 3] = Q[c(15:16, 18:21), 4] = 1
qinfo = q_info = gen_q_info(Q)

group_assignments     = rep(1, Nrespondents)
group_designmat       = matrix(1,Ngroup,1)
respondent_designmat  = matrix(1,Nrespondents,1)
Xs                    = list(pre = complete.dat.mat[,1:21], post = complete.dat.mat[,22:42])


Xdata = list(Xs=Xs,
             Q=Q,
             group_assignments=group_assignments,
             group_designmat=group_designmat,
             respondent_designmat=respondent_designmat)
initparams = gen_initial_values_longitudinal(Nrespondents = Nrespondents,
                                             Q            = Q,
                                             Ntime        = 2,
                                             Nrespcov     = Nrespcov,
                                             Ngroupcov    = Ngroupcov,
                                             Ngroup       = Ngroup,
                                             seed         = NULL,
                                             onelvl       = T)
Xdata$Ns = initparams$Ns
params   = initparams
onelvl   = T
g        = logit

### MCMC
system.time({
  mcmc_results=mcmc_sampler_main(Xdata      = Xdata,
                                 initparams = initparams,
                                 M          = 10000,
                                 onelvl     = T)
})
# user   system  elapsed
# 3664.009  133.855 4010.336
# saveRDS(mcmc_results, "~/Desktop/CDM/BayesianCDM/single_group/mcmc_results.Rds")
# mcmc_results = readRDS("~/Desktop/CDM/BayesianCDM/inst/mcmc_results.Rds")
plot(colMeans(mcmc_results[5001:10000,])[1:21], ylim = c(-9,8))
points(colMeans(mcmc_results[5001:10000,])[22:42], pch = 19)

logit(colMeans(mcmc_results[5001:10000,])[43:50])


### Variational Logit
set.seed(100)
system.time({var_results=variational_fit(Xdata          = Xdata,
                                         initial_params = initparams,
                                         g              = logit,
                                         maxiter        = 5000)
})
# user    system   elapsed
# 45444.764  1784.031 92946.523
saveRDS(var_results, "~/Desktop/CDM/BayesianCDM/single_group/var_results.Rds")
# var_results = readRDS("~/Desktop/CDM/BayesianCDM/single_group/var_results.Rds")
plot(var_results$vb_mean[1:21], ylim = c(-3.1,6.2))
points(var_results$vb_mean[22:42], pch = 19)

nresamp=1000
variational_resamp=
  data.frame(rbind(var_results$vb_mean,var_results$vb_var))%>%
  mutate(type=c('mean','var'))%>%
  pivot_longer(!'type',names_to='d')%>%
  group_by(d)%>%
  nest()%>%
  # rowwise()%>%
  mutate(data=map(data,~data.frame(i=1:nresamp,vals=rnorm(nresamp,.$value[1],.$value[2]))))%>%
  unnest(data)%>%
  pivot_wider(names_from=d,values_from=vals)
# p(variational_resamp,names_from=d)

names(variational_resamp)[1+viz_indices]
plotdf=variational_resamp[,1+viz_indices]%>%
  pivot_longer(everything())%>%
  mutate(type=map_chr(name,~strsplit(.,'_')[[1]][1]))

ggplot(plotdf) +
  geom_histogram(aes(x = value,fill=type), bins = 15) +
  facet_wrap(~ name)




# MLE
likelihood_wrapper = function(vec){
  likelihood_master(convert_vector_to_params(vec, initparams), Xdata,qinfo,g=logit)
}
system.time({mle_results = optim(convert_params_to_vector(initparams),
                                 likelihood_wrapper,
                                 control=list(maxit=10000,trace=6))
})
#saveRDS(mle_results, "~/Desktop/CDM/BayesianCDM/inst/mle_results.Rds")
mle_results = readRDS("~/Desktop/CDM/BayesianCDM/inst/mle_results.Rds")
# user   system  elapsed
# 1798.748   66.782 1868.938
plot(mle_results$par[1:21], ylim = c(-9, 8))
points(mle_results$par[22:42], pch = 19)

plot(convert_params_to_vector(initparams)[1:21], ylim = c(-9, 8))
points(convert_params_to_vector(initparams)[22:42], pch = 19)

logit(mle_results$par[43:50])





plot(colMeans(mcmc_results))
