par_fun_top <- function(mvn_fun) {
  function(j, XX, startlatent, D, latent_nofac) {
    mvn_fun(XX/exp(startlatent[,j]*0.5), D[,j], latent_nofac[,j,drop=FALSE]/exp(startlatent[,j]*0.5))
  }
}

par_fun_AR1 <- function(j, XX, startlatent, D, latent_nofac, prior_Pi_AR1) {
    rmvn_ccm(XX/exp(startlatent[,j]*0.5), D[,j], latent_nofac[,j,drop=FALSE]/exp(startlatent[,j]*0.5), prior_Pi_AR1[j], j)
}
