#Input:
#y: binary vector of transitions
#
sample_ml_pg=function(y,X,Z,beta,gamma,Sigma,group_assignments,priors){
  invSigma=solve(Sigma)
  kappa=y-0.5

  nG=dim(Z)[1]
  nX=dim(X)[2]
  nZ=dim(Z)[2]

  gamma_mean=list()

  #sample beta, for each group
  beta_samp_list=list()
  eta=apply(X*beta[group_assignments,],1,sum)
  for(j in 1:nG){
    Xj=X[group_assignments==j,]
    omega_j=map_dbl(eta[group_assignments==j],~rpg(1,h=1,z=.))
    beta_post_prec=t(Xj)%*%(c(omega_j)*Xj)+invSigma
    beta_post_mu=
      solve(beta_post_prec,t(X[group_assignments==j,])%*%
              (omega_j*kappa[group_assignments==j])+
              invSigma%*%t(gamma)%*%Z[j,])
    beta_samp_list[[j]]=spam::rmvnorm(1,beta_post_mu,solve(beta_post_prec))
  }
  new_beta=do.call(rbind,beta_samp_list)


  beta_stack=c(beta)
  gamma_stack=c(gamma)
  Z_stack=kronecker(diag(3),Z)
  gammaZ_stack=Z_stack%*%gamma_stack
  Sigma_stack=kronecker(Sigma,diag(5))
  prec_stack=solve(Sigma_stack)
  gamma_post_prec=t(Z_stack)%*%prec_stack%*%Z_stack+priors$sigma_gamma^2*diag(nX*nZ)
  gamma_post_mu=solve(gamma_post_prec,c(beta_stack%*%Z_stack))
  gamma_samp=spam::rmvnorm(1,gamma_post_mu,solve(gamma_post_prec))
  #convert back to matrix
  new_gamma=matrix(gamma_samp,nZ,nX)

  residuals=matrix(beta_stack-Z_stack%*%c(gamma_samp),nG,nX)
  new_Sigma=rWishart(1,priors$V_prior+t(residuals)%*%residuals,df=priors$v_prior+nG)[,,1]

  return(list(beta=new_beta,gamma=new_gamma,Sigma=new_Sigma))


}

#For testing the above function:
# {
#   nI=99
#   nG=5
#   group_assignments=sort(1:nI%%nG)+1
#   nX=3
#   nZ=4
#   X=cbind(rep(1,nI),matrix(rnorm(nI*(nX-1)),nI,nX-1))
#   Z=cbind(rep(1,nG),matrix(rnorm(nG*(nZ-1)),nG,nZ-1))
#   sigma_gamma=2.5
#   Sigma=rWishart(1,diag(nX),df=3)[,,1]
#   g=function(x) 1/(1+exp(-x))
#   ginv=function(x) -log(1/x-1)
#   gamma=matrix(rnorm(nX*nZ,sd=1),nZ,nX)
#   beta2=map(1:nG,~spam::rmvnorm(1,(Z%*%gamma)[.,],Sigma=1e-16*diag(nX)))%>%
#     {do.call(rbind,.)}
#
#   eta=apply(X*beta[group_assignments,],1,sum)
#   yprob=g(eta)
#   y=map_lgl(yprob,~sample(c(0,1),1,prob=c(1-.,.)))+0
#
#   priors=list(
#     V_prior=diag(nX),
#     v_prior=1,
#     sigma_gamma=sigma_gamma
#   )
# }
# rm(list=setdiff(ls(),c("y","X","Z",'beta','gamma','Sigma','priors','group_assignments','sample_ml_pg')))
# sample_ml_pg(y,X,Z,beta,gamma,Sigma,group_assignments,priors)
# #
