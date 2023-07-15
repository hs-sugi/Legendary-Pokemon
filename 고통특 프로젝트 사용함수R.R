# boxcox
#######################################################################
#                  function 1: evaluate likelihood
#######################################################################
likelihood = function(x,lam,eps){
  xlam=bocotranmat(x,lam,eps);
  xmar=standmar(x,eps)
  mu=apply(xlam,2,mean)
  sig=var(xlam)
  n=dim(x)[1]
  loglik=0
  for(i in 1:n){
    loglik = loglik + (-1/2) * (xlam[i,]-mu) %*% solve(sig) %*% (xlam[i,]-mu) +
      (-1/2) * log( det(sig)) +
      sum((lam - 1) * log(xmar[i,]))}
  return(c(loglik) )
}
#######################################################################
#                  function 2: standardize a matrix
#                              treating each row as a random vector
#                              in an iid sample
#######################################################################
standmat = function(x){
  mu = apply(x,2,mean)
  sig = var(x)
  signrt = matpower(sig,-1/2)
  return(t(t(x) - mu)%*%signrt)
}
#######################################################################
#                  function 3: marginally standardize a matrix
#                              treating each row as a random vector
#                              in an iid sample, and then transform 
#                              it to be positive
#######################################################################
standmar = function(x,eps){
  mu = apply(x,2,mean)
  sig = diag(diag(var(x)))
  signrt = matpower(sig,-1/2)
  x1 = t(t(x) - mu)%*%signrt
  x2 = t(t(x1) - apply(x1,2,min))+eps
  return(x2)
}
#######################################################################
#                  function 4: standardize a vector
#######################################################################
standvec = function(x){
  return((x - mean(x))/sd(x))}
#######################################################################
#                  function 5: boxcox transform a vector
#                              x is n dim vector
#                              lam is box cox parameter
#                              eps is distance above 0
#######################################################################
bocotranvec = function(x,lam,eps){
  n = length(x)
  x1 = standvec(x)
  x2 = x1 - min(x1) + eps
  if(abs(lam)< 10^(-10)) 
    x3 = log(x2) else
      x3 = (x2^lam - 1)/lam
  return(x3)
}
#######################################################################
#                  function 6: boxcox transform a matrix
#                              x is n x p matrix
#                              lam is p vector
#                              eps is 0.5 say
#######################################################################
bocotranmat = function(x,lam,eps){
  n = dim(x)[1]
  p = dim(x)[2]
  xlam = numeric()
  for(i in 1:p){
    xlam = cbind(xlam,bocotranvec(x[,i],lam[i],eps))}
  return(xlam)
}
#######################################################################
#                  function 7: find maximizer
#######################################################################
argmax = function(x,y){
  n = length(x)
  return(x[order(y)[n]])}
#######################################################################
#       function 8: maximizing over one lambda in gauss seidel
#######################################################################
gaussonestep=function(x,lam,eps,ilam,mlam){
  onelam=seq(from=-2,to=2,length=mlam)
  loglik = numeric() 
  for(k in 1:mlam){
    lam[ilam]=onelam[k]
    loglik=c(loglik,likelihood(x,lam,eps))}
  lamopt=argmax(onelam,loglik)
  lamout=lam 
  lamout[ilam]=lamopt
  return(lamout)
}
#######################################################################
#                function 9:  gauss seidel iterations
#######################################################################
gauss=function(x,lam,eps,mlam){
  matpower=function(a,alpha){
    a=(a+t(a))/2;tmp=eigen(a)
    return(tmp$vectors%*%diag((tmp$values)^alpha)%*%t(tmp$vectors))}
  standvec=function(x) return((x - mean(x))/sd(x)) 
  standmar=function(x,eps){
    mu=apply(x,2,mean);sig=diag(diag(var(x)));signrt=matpower(sig,-1/2)
    x1=t(t(x)-mu)%*%signrt;x2=t(t(x1)-apply(x1,2,min))+eps
    return(x2)}
  bocotranvec=function(x,lam,eps){
    n=length(x);x1=standvec(x);x2=x1-min(x1)+eps
    if(abs(lam)<10^(-10)) x3=log(x2)  
    if(abs(lam)>=10^(-10)) x3=(x2^lam - 1)/lam
    return(x3)}
  bocotranmat=function(x,lam,eps){
    n=dim(x)[1];p=dim(x)[2];xlam=numeric()
    for(i in 1:p) xlam=cbind(xlam,bocotranvec(x[,i],lam[i],eps)) 
    return(xlam)}
  likelihood=function(x,lam,eps){
    xlam=bocotranmat(x,lam,eps);xmar=standmar(x,eps)
    mu=apply(xlam,2,mean);sig=var(xlam);n=dim(x)[1]
    loglik=0
    for(i in 1:n){
      loglik = loglik+(-1/2)*(xlam[i,]-mu)%*%solve(sig)%*%(xlam[i,]-mu)+
        (-1/2)*log(det(sig))+sum((lam-1)*log(xmar[i,]))}
    return(c(loglik))}
  argmax = function(x,y)  return(x[order(y)[length(x)]])
  
  gaussonestep=function(x,lam,eps,ilam,mlam){
    onelam=seq(from=-2,to=2,length=mlam)
    loglik=numeric() 
    for(k in 1:mlam){lam[ilam]=onelam[k];loglik=c(loglik,likelihood(x,lam,eps))}
    lamopt=argmax(onelam,loglik);lamout=lam;lamout[ilam]=lamopt
    return(lamout)}
  
  p = length(lam)
  for(igauss in 1:5) for(i in 1:p){
    xlam=bocotranmat(x,lam,eps)
    #pairs(xlam) #첫뻔쟨 일리티컬 만족 못함
    lam1=gaussonestep(x,lam,eps,i,mlam) #이니셜밸류
    lam=lam1
  }
  return(lam)
}

# base
################################################################
#          Power of matrix
################################################################
matpower=function(a,alpha){
  a = (a + t(a))/2
  tmp = eigen(a)
  return(tmp$vectors%*%diag((tmp$values)^alpha)%*%
           t(tmp$vectors))}

##############################################################
#      Center X (n by p matrix)        
##############################################################
center = function(x){
  return(t(t(x)-apply(x,2,mean)))}


################################################################
#           discretize
################################################################
discretize=function(y,h){
  n=length(y);m=round(n/h)
  y=y+.00001*mean(y)*rnorm(n)
  yord = y[order(y)]
  divpt=numeric();for(i in 1:(h-1)) divpt = c(divpt,yord[i*m+1])
  y1=rep(0,n);y1[y<divpt[1]]=1;y1[y>=divpt[h-1]]=h
  for(i in 2:(h-1)) y1[(y>=divpt[i-1])&(y<divpt[i])]=i
  return(y1)
}

##############################################################
#                   symmtrize a matrix
############################################################## 
symmetry = function(a){
  return((a + t(a))/2)}



# ladle
###############################################################
#          ladle estimator 
###############################################################
ladle=function(x,y,h,nboot,method,ytype){
  r=2;n=dim(x)[1];p=dim(x)[2] 
  if(p<=10) kmax=p-2;if(p>10) kmax=floor(p/log(p)) 
  candmat=function(x,y,h,r,ytype,method){
    if(method=="sir") mat=sir_up(x,y,h,r,ytype)$sirmat
    if(method=="save") mat=save_up(x,y,h,r,ytype)$savemat
    if(method=="dr") mat=dr_up(x,y,h,r,ytype)$drmat
    return(mat)}
  phin=function(kmax,eval){den=1+sum(eval[1:(kmax+1)]);return(eval[1:(kmax+1)]/den)}
  out=candmat(x,y,h,r,ytype,method)
  eval.full=eigen(out)$values;evec.full=eigen(out)$vectors
  pn=phin(kmax,eval.full)
  
  prefn0=function(kmax,evec1,evec2){ 
    out=numeric();for(k in 0:kmax){
      if(k==0) out=c(out,0)
      if(k==1) out=c(out,1-abs(t(evec1[,1])%*%evec2[,1]))
      if(k!=0&k!=1) out=c(out,1-abs(det(t(evec1[,1:k])%*%evec2[,1:k])))}
    return(out)}
  
  fn0=0
  for(iboot in 1:nboot){ 
    bootindex=round(runif(n,min=-0.5,max=n+0.5))
    xs=x[bootindex,];ys=y[bootindex]
    mat=candmat(xs,ys,h,r,ytype,method);eval=eigen(mat)$values;evec=eigen(mat)$vectors
    fn0=fn0+prefn0(kmax,evec.full,evec)/nboot}
  
  minimizer=function(a,b) return(a[order(b)][1])
  
  fn=fn0/(1+sum(fn0));gn=pn+fn;rhat=minimizer(0:kmax,gn) 
  return(list(kset=(0:kmax),gn=gn,rhat=rhat))
}

############# update sir, save, dr are appended below

sir_up=function(x,y,h,r,ytype){
  p=ncol(x);n=nrow(x)
  signrt=matpower(var(x),-1/2)
  xc=t(t(x)-apply(x,2,mean))
  xst=xc%*%signrt
  if(ytype=="continuous") ydis=discretize(y,h)
  if(ytype=="categorical") ydis=y
  yless=ydis;ylabel=numeric()
  for(i in 1:n) {if(var(yless)!=0) {ylabel=c(ylabel,yless[1]);yless=yless[yless!=yless[1]]}}
  ylabel=c(ylabel,yless[1])
  prob=numeric();exy=numeric()
  for(i in 1:h) prob=c(prob,length(ydis[ydis==ylabel[i]])/n) 
  for(i in 1:h) exy=rbind(exy,apply(xst[ydis==ylabel[i],],2,mean))
  sirmat=t(exy)%*%diag(prob)%*%exy
  return(list(beta=signrt%*%eigen(sirmat)$vectors[,1:r],sirmat=sirmat))}


save_up=function(x,y,h,r,ytype){
  p=ncol(x);n=nrow(x)
  signrt=matpower(var(x),-1/2)
  xc=t(t(x)-apply(x,2,mean))
  xst=xc%*%signrt
  if(ytype=="continuous") ydis=discretize(y,h)
  if(ytype=="categorical") ydis=y
  ylabel=unique(ydis)
  prob=numeric() 
  for(i in 1:h) prob=c(prob,length(ydis[ydis==ylabel[i]])/n)
  vxy = array(0,c(p,p,h))
  for(i in 1:h) vxy[,,i] = var(xst[ydis==ylabel[i],]) 
  savemat=0
  for(i in 1:h){
    savemat=savemat+prob[i]*(vxy[,,i]-diag(p))%*%(vxy[,,i]-diag(p))}
  return(list(beta=signrt%*%eigen(savemat)$vectors[,1:r],savemat=savemat))}


dr_up=function(x,y,h,r,ytype){
  p=ncol(x);n=nrow(x)
  signrt=matpower(var(x),-1/2)
  xc=t(t(x)-apply(x,2,mean))
  xst=xc%*%signrt
  if(ytype=="continuous") ydis=discretize(y,h)
  if(ytype=="categorical") ydis=y
  ylabel=unique(ydis)
  prob=numeric() 
  for(i in 1:h) prob=c(prob,length(ydis[ydis==ylabel[i]])/n)
  vxy = array(0,c(p,p,h));exy=numeric()
  for(i in 1:h) {
    vxy[,,i]=var(xst[ydis==ylabel[i],])
    exy=rbind(exy,apply(xst[ydis==ylabel[i],],2,mean))}
  mat1 = matrix(0,p,p);mat2 = matrix(0,p,p)
  for(i in 1:h){
    mat1 = mat1+prob[i]*(vxy[,,i]+exy[i,]%*%t(exy[i,]))%*%
      (vxy[,,i]+exy[i,]%*%t(exy[i,]))
    mat2 = mat2+prob[i]*exy[i,]%*%t(exy[i,])}
  out = 2*mat1+2*mat2%*%mat2+2*sum(diag(mat2))*mat2-2*diag(p)
  return(list(beta=signrt%*%eigen(out)$vectors[,1:r],drmat=out))
}



################################################################
#                           sir
################################################################
sir=function(x,y,h,r,ytype){
  p=ncol(x);n=nrow(x)
  signrt=matpower(var(x),-1/2)
  xc=t(t(x)-apply(x,2,mean))
  xst=xc%*%signrt
  if(ytype=="continuous") ydis=discretize(y,h)
  if(ytype=="categorical") ydis=y
  ylabel=unique(ydis)
  prob=numeric();exy=numeric()
  for(i in 1:h) prob=c(prob,length(ydis[ydis==ylabel[i]])/n) 
  for(i in 1:h) exy=rbind(exy,apply(xst[ydis==ylabel[i],],2,mean))
  sirmat=t(exy)%*%diag(prob)%*%exy
  return(signrt%*%eigen(sirmat)$vectors[,1:r])}


################################################################
#                          save
################################################################
save=function(x,y,h,r,ytype){
  p=ncol(x);n=nrow(x)
  signrt=matpower(var(x),-1/2)
  xc=t(t(x)-apply(x,2,mean))
  xst=xc%*%signrt
  if(ytype=="continuous") ydis=discretize(y,h)
  if(ytype=="categorical") ydis=y
  ylabel=unique(ydis)
  prob=numeric() 
  for(i in 1:h) prob=c(prob,length(ydis[ydis==ylabel[i]])/n)
  vxy = array(0,c(p,p,h))
  for(i in 1:h) vxy[,,i] = var(xst[ydis==ylabel[i],]) 
  savemat=0
  for(i in 1:h){
    savemat=savemat+prob[i]*(vxy[,,i]-diag(p))%*%(vxy[,,i]-diag(p))}
  return(signrt%*%eigen(savemat)$vectors[,1:r])}

################################################################
#                                 dr 
################################################################
dr=function(x,y,h,r,ytype){
  p=ncol(x);n=nrow(x)
  signrt=matpower(var(x),-1/2)
  xc=t(t(x)-apply(x,2,mean))
  xst=xc%*%signrt
  if(ytype=="continuous") ydis=discretize(y,h)
  if(ytype=="categorical") ydis=y
  ylabel=unique(ydis)
  prob=numeric() 
  for(i in 1:h) prob=c(prob,length(ydis[ydis==ylabel[i]])/n) 
  
  vxy = array(0,c(p,p,h));exy=numeric()
  for(i in 1:h) {
    vxy[,,i]=var(xst[ydis==ylabel[i],]) 
    exy=rbind(exy,apply(xst[ydis==ylabel[i],],2,mean))} 
  
  mat1 = matrix(0,p,p);mat2 = matrix(0,p,p)
  for(i in 1:h){
    mat1 = mat1+prob[i]*(vxy[,,i]+exy[i,]%*%t(exy[i,]))%*% 
      (vxy[,,i]+exy[i,]%*%t(exy[i,]))
    mat2 = mat2+prob[i]*exy[i,]%*%t(exy[i,])}
  out = 2*mat1+2*mat2%*%mat2+2*sum(diag(mat2))*mat2-2*diag(p) 
  return(signrt%*%eigen(out)$vectors[,1:r])
}

