#####

## measurement:
##    z_t     =  mH_t * h_t
## transition:
##    h_{t+1} =  mF * h_t + mB * x_t + mQ * e_t
## note e_t follows (0, I)

##### to sample from the smoothed distribution
smooth_samp <- function(mZ,mX,lH,lH0=NULL,mF,mB,mQ,iT,ip,iq,is,h0,P0=NULL,X0)
## Inputs
## mZ : T   by p   matrix representing observations
## mX : T   by s   matrix representing regressors
## lH : T list, p by q
## mF : q by q
## mB : q by s
## mQ : q by q for w_t
## h0 : q   by 1   matrix ... accurate initial state
## P0 : q   by q   matrix ... accurate initial state's cov
## X0 : s   by 1   matrix ... x at time step 0
## outputs
## KF part
## h1 : T   by q   expectation of h a priori
## h2 : T   by q   expectation of h a posteriori
## P1 : Tqq array  variance of h a priori
## P2 : Tqq array  variance of h a posteriori
## mv  : T   by p   residuals in measurement eq.
## aS  : Tpp array  covariance of z
## aK  : Tqp array  Kalman gain
## SM part
## mu : T   by p+q E[eps |Z X]
## u0 : p+q by 1
## mr : T   by q
## r0 : q   by 1
## me : T   by p
{


  # simulation



  mE = matrix(rnorm(iT*iq),iT,iq)
  mZZ = matrix(0,iT,ip)
  mhh = matrix(0,iT,iq)
  if(is.null(P0)){
    hh0 = h0; P0=matrix(0,iq,iq)
  } else{
      hh0 = h0 + t(chol(P0))%*%rnorm(iq)
  }

  mH = lH[[1]]
  mhh[1,] = mF %*% hh0 + mB %*% X0 + mQ%*%mE[1,]
  mZZ[1,] = mH %*% mhh[1,]
  for(iter in 2:iT){
    mH = lH[[iter]]
    mhh[iter,] = mF %*% mhh[iter-1,] + mB %*% mX[iter-1,] + mQ%*%mE[iter,]
    mZZ[iter,] = mH %*% mhh[iter,]
  }

  # mZZ, mhh and mE from 0 to iT-1
  mE = rbind(mE,0)

  if(is.null(lH0)) lH0=lH
  u1 = smoothing(mZ=mZ,mX=mX,lH=lH0,mF=mF,mB=mB,mQ=mQ,iT=iT,ip=ip,iq=iq,is=is,h0=h0,P0=P0,X0=X0)
  u2 = smoothing(mZ=mZZ,mX=mX,lH=lH0,mF=mF,mB=mB,mQ=mQ,iT=iT,ip=ip,iq=iq,is=is,h0=h0,P0=P0,X0=X0)
  mU = mE - u2 + u1
  # from 0 to iT

  mH = lH[[1]]
  mhh[1,] = mF %*% hh0 + mB %*% X0 + mQ%*%mU[1,]
  mZZ[1,] = mH %*% mhh[1,]
  for(iter in 2:iT){
    mH = lH[[iter]]
    mhh[iter,] = mF %*% mhh[iter-1,] + mB %*% mX[iter-1,] + mQ%*%mU[iter,]
    mZZ[iter,] = mH %*% mhh[iter,]
  }

  return(list(mh=mhh,mZ=mZZ))
}

smoothing <- function(mZ,mX,lH,mF,mB,mQ,iT,ip,iq,is,h0,P0,X0)
{
  QQ = tcrossprod(mQ)

  # filtering

  #mv=array(0,dim=c(iT,ip))# residual in z eq.
  mv=list(); length(mv)=iT
  #IS=array(0,dim=c(iT,ip,ip))
  IS=list(); length(IS)=iT
  aK=list(); length(aK)=iT
  #aK=array(0,dim=c(iT,iq,ip))# Kalman gain

  #Predict
  h1 = mF%*%h0 + mB%*%X0
  P1 = mF%*%P0%*%t(mF) + QQ
  #Update
  mH = matrix(c(lH[[1]]),ncol=iq); vz=mZ[1,!is.na(mZ[1,])]
  mv[[1]] = vz-mH%*%h1 ## smoothing
  aS = mH%*%P1%*%t(mH); IS[[1]] = chol2inv(chol(aS)) ## smoothing
  aK[[1]] = P1%*%t(mH)%*%IS[[1]] ## smoothing
  h2 = h1+c(aK[[1]]%*%mv[[1]])
  P2 = (diag(1,iq)-aK[[1]]%*%mH) %*% P1

  for(iter in 2:iT){
    #Predict
    h1 = mF%*%h2 + mB%*%mX[(iter-1),]
    P1 = mF%*%P2%*%t(mF) + QQ
    #Update
    mH = matrix(c(lH[[iter]]),ncol=iq); vz=mZ[iter,!is.na(mZ[iter,])]
    mv[[iter]] = vz-mH%*%h1 ## smoothing
    aS = mH%*%P1%*%t(mH); IS[[iter]] = chol2inv(chol(aS)) ## smoothing
    aK[[iter]] = P1%*%t(mH)%*%IS[[iter]] ## smoothing
    h2 = h1+c(aK[[iter]]%*%mv[[iter]])
    P2 = (diag(1,iq)-aK[[iter]]%*%mH) %*% P1
  }

  # smoothing

  me = list(); length(me)=iT
  mr = matrix(0,iT,iq)
  mu = matrix(0,iT,iq)

  for(iter in iT:2){
    mH = matrix(c(lH[[iter]]),ncol=iq)
    fk = mF %*% aK[[iter]]
    me[[iter]] = IS[[iter]]%*%mv[[iter]] - t(fk)%*%mr[iter,]
    mr[iter-1,] = crossprod(mH,me[[iter]]) + crossprod(mF,mr[iter,])
    mu[iter,] = crossprod(mQ,mr[iter,])
  }

  mH = matrix(c(lH[[1]]),ncol=iq)
  fk = mF %*% aK[[1]]
  me[[1]] = IS[[1]]%*%mv[[1]] - t(fk)%*%mr[1,]
  r0 = crossprod(mH,me[[1]]) + crossprod(mF,mr[1,])
  mu[1,] = crossprod(mQ,mr[1,])

  #u0=t(mQ)%*%r0
  mu = rbind(c(t(mQ)%*%r0), mu)
  ## from 0 to iT

  return(mu)
}


SimKF<-function(mX,lH,mF,mB,mQ,iT,ip,iq,is,h0,P0=NULL,X0)
## inputs
## mX  : T   by s   matrix representing regressors
## lH : T list, p by q
## mF : q by q
## mB : q by s
## mQ : q by q for w_t
## mG  : p+q by p+q matrix ... parameters of errors (R' Q')'
## h0 : q   by 1   matrix ... initial state
## P0 : q   by q   matrix ... accurate initial state's cov
## X0 : s   by 1   matrix ... x at time step 0
## output
## mZ  : T   by p   observations
## mh  : T   by q   underlying process
{
  mh=matrix(0,iT,iq)# true states
  me=matrix(rnorm(iT*iq),iT,iq)# innovations
  mZ=matrix(0,iT,ip) # data
  if(is.null(P0)){h0 = h0; P0=matrix(0,iq,iq)}
  else{h0 = h0 + t(chol(P0))%*%rnorm(iq)}

  mH = lH[[1]]; dim(mH) = c(ip,iq)
  mh[1,] = mF %*% h0 + mB %*% X0 + mQ%*%me[1,]
  mZ[1,] = mH %*% mh[1,]
  for(iter in 2:iT){
    mH = lH[[iter]]; dim(mH) = c(ip,iq)
    mh[iter,] = mF %*% mh[iter-1,] + mB %*% mX[iter-1,] + mQ%*%me[iter,]
    mZ[iter,] = mH %*% mh[iter,]
  }
  me = rbind(me,0) # from 0 to iT

  return(list(mZ=mZ,mh=mh,me=me))
}
