
#just take the first "Q":
Nq=Nq_list[[1]]
singe_delta=Q_to_delta(Qs[[1]])

beta_mat=beta_mat*delta
ltheta=beta_mat%*%A

theta=logistic(ltheta)
Nprofile=2^Nskill

#sample ystar
nct = matrix(NA,Nprofile,Ntime)
for(c in 1:Nprofile){
  for(t in 1:Ntime){
    nct[c,t]=sum(alpha_mat[,t]==c)
  }
}
ystar=map(1:Ntime,~matrix(NA,Nq,Nprofile))
# matrix(NA,Nq_total,Nprofile)
for(t in 1:Ntime){
  for(j in 1:Nq){
    for(c in 1:Nprofile){
      ystar[[t]][j,c]=rpg(h=max(nct[c,t],1),z=ltheta[j,c])
    }
  }
}


#sample beta
kappa = map(1:Ntime,~matrix(NA,Nq_total,Nprofile))
for(t in 1:Ntime){
  for(j in 1:Nq_total){
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
    Aconcat=rbind(A,A)
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

