smooth_samp_xx <- function(mZ,mX,lH,lH0=NULL,mF,mB,mQ,iT,ip,iq,is,h0,P0=NULL,X0,Lambda)
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
  Y <- mZ-mZZ
  # from 0 to iT

  mH = lH[[1]]
  mhhE <- mhh*0
  mhhE[1, ] = mF %*% hh0 + mB %*% X0 + mQ%*%mE[1,]
  for(iter in 2:iT){
    mH = lH[[iter]]
    mhhE[iter,] = mF %*% mhhE[iter-1,] + mB %*% mX[iter-1,] + mQ%*%mE[iter,]
  }

  Pi_comp <- mF
  Q_comp <- mQ
  n_T <- iT
  n_vars <- ip
  n_comp <- iq
  h00 <- h0-h0

  u4 <- mfbvar::smoother(Y, Lambda, Pi_comp, Q_comp, n_T, n_vars, n_comp, h00, P0)

  return(list(mh=u4+mhhE))
}

