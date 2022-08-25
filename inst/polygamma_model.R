require(tidyverse)

{
  dat = read.table("single_group/prepost1.txt")[,-1]
  colnames(dat) = c(paste("Pre",1:21),paste("Post",1:21)) #gives items names for mirt
  dat[dat == 99] = NA #replace 99 with na

  ###### Remove missing
  complete.dat = dat[complete.cases(dat),]

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

  Nquestions=dim(Q)[1]
  Nprofile=2^Nskill

  Y=list(as.matrix(complete.dat[,1:21]),
         as.matrix(complete.dat[,22:42]))

  Ns=list(Nrespondents=Nrespondents,Ntime=Ntime,
          Nquestions=Nquestions,Nrespcov=Nrespcov,
          Ngroupcov=Ngroupcov,Ngroup=Ngroup,
          Nskill=Nskill,Nprofile=Nprofile)
}

#initialize
alpha=gen_alpha(Nprofile,Nskill)
A=gen_A(alpha)

Delta=Q_to_delta(Q)
beta_mat=abs(matrix(rnorm(Nquestions*Nprofile,sd=.5),Nquestions,Nprofile))*
  Delta
beta_mat[,1]=rnorm(Nquestions,sd=.25)
beta_vec=beta_mat[Delta==1]

psi=beta_mat%*%t(A)
theta=logit(psi)
sigma_jp=matrix(1,Nquestions,Nprofile)
#Enforce positivity in all but intercept, as is done in the paper
Ljp=matrix(0,Nquestions,Nprofile)
Ljp[,1]=-Inf
require(truncnorm)

myt=system.time({
  M=1000

  beta_samples=matrix(NA,length(beta_vec),M)
  beta_samples[,1]=c(beta_vec)

  prof_samples=list(matrix(NA,Nrespondents,M),
                    matrix(NA,Nrespondents,M))
  prof_samples[[1]][,1]=random_profiles()
  prof_samples[[2]][,1]=random_profiles()

  for(m in 2:M){
    for(t in 1:Ntime){
      #sample profile
      nci=gen_nci(prof_samples[[t]][,m-1])
      prof_sample=sapply(1:Nrespondents,function(r) sample_profile(r,t,theta,nci))
      prof_samples[[t]][,m]=prof_sample

      #sample augmented data
      nc=as.numeric(table(prof_sample))
      ystar=sample_ystar(psi,nc)

      #sample beta
      njc=gen_njc(prof_samples[[t]][,m],X[[t]])
      kappa=sweep(njc,2,nc/2,'-') #t(apply(nmat,1,function(x) x-sum(x)/2))
      z=kappa/ystar

      vjp=get_vjp(ystar,A,sigma_jp)
      mjp=get_mjp(vjp,z,beta_mat)

      beta_mat=matrix(rtruncnorm(length(mjp),mean=c(mjp),sd=c(vjp),a=c(Ljp)),
                      dim(beta_mat)[1],dim(beta_mat)[2])*Delta


      psi=beta_mat%*%t(A)
      theta=logit(psi)

      beta_vec=beta_mat[Delta==1]
      beta_samples[,1]=c(beta_vec)
    }
    print(m)
  }
})

profmat1=prof_samples[[1]][,1:(m-1)]
profmat2=prof_samples[[2]][,1:(m-1)]

forward_transitions=rep(NA,Nskill)
for(s in 1:Nskill){
  profiles=which(q_info$skill_in_profile[s,])
  t1_havent=which(!(profmat1 %in% profiles))
  t1_have=which((profmat1 %in% profiles))
  forward_transitions[s]=mean(profmat2[t1_havent]%in%profiles)
  backward_transitions[s]=mean(!(profmat2[t1_have]%in%profiles))
}

forward_transitions
backward_transitions
