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
}


beta=rnorm(Nquestions*2,sd=.1)
Ns=list(Nrespondents=Nrespondents,Ntime=Ntime,
        Nquestions=Nquestions,Nrespcov=Nrespcov,
        Ngroupcov=Ngroupcov,Ngroup=Ngroup,
        Nskill=Nskill,Nprofile=Nprofile)
value_key=c(rep("intercepts",Nquestions),
            rep("base_effects",Nquestions))
myparams=list(Ns=Ns,
              theta=beta,
              value_key=value_key)

theta=beta_to_theta(beta)
X=list(as.matrix(complete.dat[,1:21]),
       as.matrix(complete.dat[,22:42]))

Nskill=dim(Q)[2]
Nprofile=2^Nskill
ilogit=function(x) -log(1/x-1)

myt=system.time({
  M=1000

  beta_samples=matrix(NA,length(init_beta),M)
  beta_samples[,1]=c(beta)

  prof_samples=list(matrix(NA,Nrespondents,M),
                    matrix(NA,Nrespondents,M))
  prof_samples[[1]][,1]=random_profiles()
  prof_samples[[2]][,1]=random_profiles()

  for(m in 2:M){
    for(t in 1:Ntime){
      #sample profile
      nci=gen_nci(prof_samples[[t]][,m-1])
      profs_t=sapply(1:Nrespondents,function(r) sample_profile(r,t,theta,nci))
      prof_samples[[t]][,m]=profs_t

      #sample augmented data
      psi=ilogit(theta)
      nc=as.numeric(table(profs_t))
      ystar=sample_ystar(psi,nc)
      matrix(NA,Nquestions,Nprofile)

      #sample beta
      njc=gen_njc(prof_samples[[t]][,m],X[[t]])
      kappa=sweep(njc,2,nc/2,'-') #t(apply(nmat,1,function(x) x-sum(x)/2))
      z=kappa/ystar

      A=model.matrix(~factor(profs_t)-1)

      beta0=beta[value_key=='intercepts']
      beta1=beta[value_key=='base_effects']
      sigmaj0=1
      sigmaj1=1

      #p=0
      ztilde0=z[,1]-beta1
      ApT_omegaj0=ystar[,1] #sapply(1:Nquestions,function(i) ystar[i,1])
      vj0=1/(ApT_omegaj0+1/sigmaj0^2)
      mj0=vj0^2 * c(ApT_omegaj0*ztilde0)

      #p=1
      tmp=apply(Q==1,1,which)+1
      ztilde1=sapply(1:Nquestions,function(i) z[i,tmp[i]])-beta0
      ApT_omegaj1=sapply(1:Nquestions,function(i) ystar[i,tmp[i]])
      vj1=1/(ApT_omegaj1+1/sigmaj1^2)
      mj1=vj1^2 * c(ApT_omegaj1*ztilde1)

      beta=c(rnorm(Nquestions,mj0,vj0),
                rnorm(Nquestions,mj1,vj1))
      theta=beta_to_theta(beta)
      beta_samples[,m]=c(beta)
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
