// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#define _USE_MATH_DEFINES // for C++
#include <cmath>


class KF {
protected:
  arma::mat r, P_TT;
public:
  arma::mat y, G, H, a, a_tt, a_tT, Tt, c, intercept, d, a1, P1, Z, a_t1, P_t1, v;
  arma::cube P, FF_inv, L, N;
  unsigned int n_T, n_vars, n_state;
  void set_pars(arma::mat y_, arma::mat Z_, arma::mat c_, arma::mat G_, arma::mat Tt_, arma::mat d_, arma::mat H_, arma::mat a1_, arma::mat P1_, arma::mat intercept_);
  void filter();
  void smoother(arma::mat r_T);
  void simulator();
};

void KF::set_pars(arma::mat y_, arma::mat Z_, arma::mat c_, arma::mat G_, arma::mat Tt_, arma::mat d_, arma::mat H_, arma::mat a1_, arma::mat P1_, arma::mat intercept_) {
  y = y_;
  Z = Z_;
  c = c_;
  G = G_;
  Tt = Tt_;
  d = d_;
  H = H_;
  a1 = a1_;
  P1 = P1_;
  n_T = y_.n_rows;
  n_vars = y_.n_cols;
  n_state = Tt_.n_rows;
  intercept = intercept_;

  v = arma::mat(n_T, n_vars).fill(NA_REAL);
  a = arma::mat(n_T, n_state).fill(NA_REAL);
  a_tt = arma::mat(n_T, n_state).fill(NA_REAL);
  r = arma::mat(n_T, n_state).fill(0);
  a_tT = arma::mat(n_T, n_state).fill(NA_REAL);

  P = arma::cube(n_state, n_state, n_T).fill(NA_REAL);
  FF_inv = arma::cube(n_vars, n_vars, n_T).fill(NA_REAL);
  L = arma::cube(n_state, n_state, n_T).fill(NA_REAL);
  N = arma::cube(n_state, n_state, n_T).fill(NA_REAL);
}

void KF::filter() {

  arma::uvec obs_vars(n_vars);
  arma::mat Zt;
  arma::mat Gt;
  arma::uvec t_vec(1);
  arma::mat M_t;
  arma::mat FF_inv_t;
  arma::mat K_t;
  arma::mat v_t;
  arma::mat a_t = a1.t();
  arma::mat P_t = P1;

  a.row(0) = a_t;
  P.slice(0) = P_t;


  for (uword t = 0; t < n_T; t++) {
    obs_vars = find_finite(y.row(t));
    t_vec(0) = t;

    Zt = Z.rows(obs_vars);
    Gt = G.rows(obs_vars);

    v_t = y.submat(t_vec, obs_vars) - a_t * Zt.t() - c.submat(t_vec, obs_vars) - intercept.cols(obs_vars);
    v.submat(t_vec, obs_vars) = v_t;

    M_t = P_t * Zt.t() + H * Gt.t();

    FF_inv_t = inv_sympd(symmatu(Zt * M_t + Gt * trans(Gt + Zt * H)));
    FF_inv.slice(t).submat(obs_vars, obs_vars) = FF_inv_t;

    K_t = Tt * M_t * FF_inv_t;

    L.slice(t) = Tt - K_t * Zt;
    N.slice(t) = P_t * L.slice(t).t() - H * Gt.t() * K_t.t();
    a_tt.row(t) = a_t + v_t * FF_inv_t * M_t.t();

    if (t < n_T - 1) {
      a_t = a_t * Tt.t() + d.row(t) + v_t * K_t.t();
      a.row(t+1) = a_t;
      P_t = Tt * N.slice(t) + H * H.t();
      P_t = symmatu(P_t);
      P.slice(t+1) = P_t;
    } else {
      a_t1 = a_t * Tt.t() + d.row(t) + v_t * K_t.t();
      P_t1 = Tt * N.slice(t) + H * H.t();
      P_t1 = symmatu(P_t1);
      P_TT = P_t - M_t * FF_inv_t * M_t.t();
    }
  }


}

void KF::smoother(arma::mat r_T) {
  arma::uvec obs_vars(n_vars);
  arma::mat Zt;
  arma::uvec t_vec(1);
  arma::mat FF_inv_t;
  arma::mat r_t = r_T;
  r.row(y.n_rows-1) = r_t;
  for (uword t = y.n_rows - 1; t >= 1; t--) {
    obs_vars = find_finite(y.row(t));
    t_vec(0) = t;
    Zt = Z.rows(obs_vars);
    a_tT.row(t) = a_tt.row(t) + r_t * trans(N.slice(t));
    r_t = v.submat(t_vec, obs_vars) * trans(FF_inv.slice(t).submat(obs_vars, obs_vars)) * Zt + r_t * L.slice(t);
    r.row(t-1) = r_t;
  }

  obs_vars = find_finite(y.row(0));
  Zt = Z.rows(obs_vars);
  a_tT.row(0) = a_tt.row(0) + r_t * trans(N.slice(0));
}

void KF::simulator() {
  arma::uvec obs_vars;
  arma::mat y_sim = arma::mat(size(y)).fill(NA_REAL);
  arma::mat Zt;
  arma::mat Gt;
  arma::uvec t_vec(1);
  arma::mat M_t = arma::mat(n_state, n_vars).fill(NA_REAL);
  arma::mat FF_inv_t = arma::mat(n_vars, n_vars).fill(NA_REAL);
  arma::mat K_t = arma::mat(n_state, n_vars).fill(NA_REAL);

  c = arma::mat(size(y)).fill(NA_REAL);
  d = arma::mat(n_T, n_state);

  a.row(0) = a1.t();
  P.slice(0) = P1;

  for (uword t = 0; t < n_T; t++) {
    obs_vars = find_finite(y.row(t));
    t_vec(0) = t;

    Zt = Z.rows(obs_vars);
    Gt = G.rows(obs_vars);
    // obtain c_t and d
    // generate alpha_t
    // generate y_t
  }
}


class KF_ragged: public KF {

public:
  arma::mat a_TbTb, P_TbTb, Phi, F_Phi, Sigma, Lambda, y_Tb, Omega, W, F_Phi_c, a_Tb1, P_Tb1, Lambda_companion, Sigma_chol, Z1;
  unsigned int n_q, T_b, n_lags, n_m;
  void set_ragged_pars(arma::mat Phi_, arma::mat Sigma_, arma::mat Lambda_, int n_q_, int T_b_, arma::mat Z1_);
  void compact_to_companion(arma::mat Lambda_);
  void original_to_compact(arma::mat y_Tb_);
  arma::mat create_d(int T_end_);
};

void KF_ragged::set_ragged_pars(arma::mat Phi_, arma::mat Sigma_, arma::mat Lambda_, int n_q_, int T_b_, arma::mat Z1_) {
  Phi = Phi_;
  Sigma = Sigma_;
  Sigma_chol = trans(chol(Sigma));
  Lambda = Lambda_;
  n_q = n_q_;
  T_b = T_b_;
  Z1 = Z1_;
  n_vars = Phi.n_rows;
  n_m = n_vars - n_q;
  n_lags = (Phi.n_cols - 1)/n_vars;
  n_state = n_q*(n_lags + 1);
  a_TbTb = arma::mat(1, (n_lags+1)*n_vars);
  P_TbTb = arma::mat((n_lags+1)*n_vars, n_vars*(n_lags+1));
  F_Phi = arma::mat((n_lags+1)*n_vars, n_vars*(n_lags+1));
  a_TbTb.fill(0);
  P_TbTb.fill(0);
  F_Phi.fill(0);
  F_Phi_c.fill(0);

  v = arma::mat(T_b, n_vars).fill(NA_REAL);
  a = arma::mat(T_b, n_state).fill(NA_REAL);
  a_tt = arma::mat(T_b, n_state).fill(NA_REAL);
  r = arma::mat(T_b, n_state).fill(0);
  a_tT = arma::mat(T_b, n_state).fill(NA_REAL);

  P = arma::cube(n_state, n_state, T_b).fill(NA_REAL);
  FF_inv = arma::cube(n_vars, n_vars, T_b).fill(NA_REAL);
  L = arma::cube(n_state, n_state, T_b).fill(NA_REAL);
  N = arma::cube(n_state, n_state, T_b).fill(NA_REAL);
}

void KF_ragged::compact_to_companion(arma::mat Lambda_) {
  F_Phi = arma::mat(n_vars*(n_lags+1), n_vars*(n_lags+1), fill::zeros);
  F_Phi.submat(0, 0, n_vars - 1, n_vars*n_lags-1) = Phi.cols(0, n_vars*n_lags-1);
  F_Phi.submat(n_vars, 0, n_vars*(n_lags+1)-1, n_vars*n_lags - 1) = arma::eye(n_vars*n_lags, n_vars*n_lags);
  F_Phi_c = arma::mat(n_vars*(n_lags+1), 1, fill::zeros);
  F_Phi_c.rows(0, n_vars - 1) = Phi.col(n_vars*n_lags);
  Omega = arma::mat(n_vars*(n_lags+1), n_vars*(n_lags+1), fill::zeros);
  Omega.submat(0, 0, n_vars-1, n_vars-1) = Sigma_chol;
  arma::mat X_mat = y.submat(T_b - n_lags - 1, 0, T_b - 1, n_vars - n_q - 1);
  X_mat = arma::trans(arma::flipud(X_mat));
  arma::mat a_mat = arma::reshape(a_tt.row(T_b-1), n_q, n_lags + 1);
  arma::mat Xa_mat = arma::join_cols(X_mat, a_mat);
  a_TbTb = arma::reshape(Xa_mat, n_vars*(n_lags+1), 1);
  P_TbTb = arma::mat(n_vars*(n_lags+1), n_vars*(n_lags+1));
  P_TbTb.fill(0);
  int n_m = n_vars - n_q;
  for (unsigned int i = 0; i < n_lags + 1; i++) {
    for (unsigned int j = i; j < n_lags + 1; j++) {
      P_TbTb.submat(n_m+j*n_vars, n_m+i*n_vars, n_vars-1+j*n_vars, n_vars-1+i*n_vars) = P_TT.submat(j*n_q, i*n_q, n_q-1+j*n_q, n_q-1+i*n_q);
    }
  }
  P_TbTb = symmatl(P_TbTb);

  a_Tb1 = F_Phi * a_TbTb + F_Phi_c;
  P_Tb1 = F_Phi * P_TbTb * F_Phi.t() + Omega*Omega.t();

  Lambda_companion = arma::mat(n_vars, (n_vars*(n_lags+1)), fill::zeros);
  Lambda_companion.submat(0, 0, n_m-1, n_m-1) = arma::eye(n_m, n_m);

  for (uword i = 0; i < Lambda_.n_cols/n_q; i++) {
    Lambda_companion(span(n_m, n_vars - 1), span(n_m + i*n_vars, (i+1)*n_vars-1)) = Lambda_.cols(i*n_q, (i+1)*(n_q)-1);
  }
}

void KF_ragged::original_to_compact(arma::mat y_Tb_) {
  y_Tb = y_Tb_;
  y = y_Tb_;
  n_T = T_b;
  arma::mat Phi_mm(n_m, n_m*n_lags);
  arma::mat Phi_mq(n_m, n_q*n_lags);
  arma::mat Phi_qm(n_q, n_m*n_lags);
  arma::mat Phi_qq(n_q, n_q*n_lags);

  for (unsigned int i = 0; i < n_lags; i++) {
    Phi_mm.cols(i*n_m, (i+1)*n_m - 1) = Phi.submat(0, i*n_vars, n_m - 1, i*n_vars + n_m - 1);
    Phi_mq.cols(i*n_q, (i+1)*n_q - 1) = Phi.submat(0, i*n_vars+n_m, n_m - 1, (i+1)*n_vars - 1);
    Phi_qm.cols(i*n_m, (i+1)*n_m - 1) = Phi.submat(n_m, i*n_vars, n_vars - 1, i*n_vars + n_m - 1);
    Phi_qq.cols(i*n_q, (i+1)*n_q - 1) = Phi.submat(n_m, i*n_vars+n_m, n_vars - 1, (i+1)*n_vars - 1);
  }

  Z = arma::mat(n_vars, n_q*(n_lags + 1));
  Z.fill(0);
  Z.submat(0, n_q, n_m - 1, n_q*(n_lags + 1) - 1) = Phi_mq;
  Z.submat(n_m, 0 ,n_vars - 1, Lambda.n_cols - 1) = Lambda;

  Tt = arma::mat(n_q*(n_lags + 1), n_q*(n_lags + 1));
  Tt.fill(0);
  Tt.submat(0, 0, n_q - 1, n_q*n_lags - 1) = Phi_qq;
  Tt.submat(n_q, 0, n_q*(n_lags+1)-1, n_q*n_lags - 1) = eye(n_lags*n_q, n_lags*n_q);

  G = arma::mat(n_vars, n_vars);
  G.fill(0);
  G.rows(0, n_m - 1) = Sigma_chol.rows(0, n_m - 1);

  H = arma::mat(n_q*(n_lags + 1), n_vars);
  H.fill(0);
  H.rows(0, n_q - 1) = Sigma_chol.rows(n_m, n_vars - 1);

  arma::mat X(T_b, n_m*n_lags);
  X.fill(0);
  X.row(0) = reshape(trans(flipud(Z1.cols(0, n_m-1))), 1, n_lags*n_m);
  for (unsigned int i = 1; i < n_lags; i++) {
    X.row(i).cols(0, i*n_m - 1) = reshape(trans(flipud(y_Tb(span(0, i-1), span(0, n_m - 1)))), 1, i*n_m);
    X.row(i).cols(i*n_m, n_lags*n_m - 1) = reshape(trans(flipud(Z1(span(i, n_lags - 1), span(0, n_m - 1)))), 1, (n_lags-i)*n_m);
  }
  for (unsigned int i = n_lags; i < T_b - 1; i++) {
    X.row(i) = reshape(trans(flipud(y_Tb(span(i-n_lags, i-1), span(0, n_m - 1)))), 1, n_lags*n_m);
  }

  c = arma::mat(T_b, n_vars, fill::zeros);
  c.cols(0, n_m - 1) = X * trans(Phi_mm);
  intercept = arma::mat(1, n_vars, fill::zeros);
  intercept.cols(0, n_m - 1) = trans(Phi.submat(0, n_vars*n_lags, n_m - 1, n_vars*n_lags));

  W = arma::mat(T_b, n_m*n_lags + 1, fill::ones);
  W(span(0, T_b - 2), span(0, n_m*n_lags - 1)) = X.rows(1, T_b - 1);
  W.row(T_b-1).cols(0, n_m*n_lags - 1) = reshape(trans(flipud(y_Tb.submat(T_b-n_lags, 0, T_b-1, n_m - 1))), 1, n_lags*n_m);
  d = arma::mat(T_b, n_q*(n_lags + 1), fill::zeros);
  arma::mat Beta_W = join_rows(Phi_qm, Phi.submat(n_m, n_vars*n_lags, n_vars-1, n_vars*n_lags));
  d.cols(0, n_q - 1) = W * trans(Beta_W);

  arma::mat means = solve(eye(n_vars, n_vars) - Phi.cols(0, n_vars*n_lags - 1) * repmat(eye(n_vars, n_vars), n_lags, 1), Phi.col(n_vars*n_lags));
  a1 = arma::mat(n_q*(n_lags+1), 1, fill::zeros);
  a1.rows(0, n_q*n_lags - 1) = reshape(trans(flipud(Z1.cols(n_m, n_vars - 1))), 1, n_lags*n_q).t();
  arma::mat W0 = arma::mat(1, n_m*n_lags + 1, fill::ones);
  W0.cols(0, n_m*n_lags - 1) = X.row(0);
  arma::mat d0 = arma::mat(1, n_q*(n_lags + 1), fill::zeros);
  d0.cols(0, n_q - 1) = W0 * trans(Beta_W);
  a1 = Tt * a1 + d0.t();
  P1 = H * H.t();
}

arma::mat KF_ragged::create_d(int T_end_) {
  arma::mat d = arma::join_rows(trans(Phi.submat(0, n_vars*n_lags, n_vars - 1, n_vars*n_lags)), arma::mat(1, n_vars*n_lags, fill::zeros));
  d = repmat(d, T_end_, 1);
  return d;
}

arma::mat companion_reshaper(arma::mat obj_, unsigned int n_m_, unsigned int n_q_, unsigned int n_T_, unsigned int n_lags_) {
  unsigned int n_vars_ = n_m_ + n_q_;
  arma::mat temp = arma::reshape(obj_, n_vars_*n_T_, n_lags_ + 1);
  return arma::reshape(temp.rows(n_m_*n_T_, n_vars_*n_T_-1), n_T_, n_q_ * (n_lags_ + 1));
}

//' @title Kalman filter and smoother
//'
//' @description Kalman filter and smoother (\code{kf_ragged}) and simulation smoother (\code{kf_sim_smooth}) for mixed-frequency data with ragged edges. This function is more computationally efficient than using a companion form representation.
//' @param y_ matrix with the data
//' @param Phi_ matrix with the autoregressive parameters, where the last column is the intercept
//' @param Sigma_ error covariance matrix
//' @param Lambda_ aggregation matrix (for quarterly variables only)
//' @param n_q_ number of quarterly variables
//' @param T_b_ final time period where all monthly variables are observed
//' @keywords internal
//' @return For \code{kf_ragged}, a list with elements:
//' \item{a}{The one-step predictions (for the compact form)}
//' \item{a_tt}{The filtered estimates (for the compact form)}
//' \item{a_tT}{The smoothed estimates (for the compact form)}
//' \item{Z_tT}{The smoothed estimated (for the original form)}
//' @details The returned matrices have the same number of rows as \code{y_}, but the first \code{n_lags} rows are zero.
// [[Rcpp::export]]
Rcpp::List kf_ragged(arma::mat y_, arma::mat Phi_, arma::mat Sigma_, arma::mat Lambda_, arma::mat Z1_, int n_q_, unsigned int T_b_) {

  // Initialization of variables
  KF_ragged kf_obj;
  KF kf_end;

  // Initialization of filter
  arma::mat y_Tb = y_.rows(0, T_b_ - 1);
  kf_obj.set_ragged_pars(Phi_, Sigma_, Lambda_, n_q_, T_b_, Z1_);
  kf_obj.original_to_compact(y_Tb);

  unsigned int n_vars, n_lags, n_m, n_q, T_b, T_end, T_full;
  n_vars  = kf_obj.n_vars;
  n_lags  = kf_obj.n_lags;
  n_m     = kf_obj.n_m;
  n_q     = kf_obj.n_q;
  T_b     = kf_obj.T_b;
  T_full  = y_.n_rows;
  T_end   = T_full - T_b;

  arma::mat alpha_t1, r_T;
  arma::mat a = arma::mat(T_full, n_q*(n_lags+1)).fill(NA_REAL);
  arma::mat a_tt = arma::mat(T_full, n_q*(n_lags+1)).fill(NA_REAL);
  arma::mat a_tT = arma::mat(T_full, n_q*(n_lags+1)).fill(NA_REAL);
  arma::mat Z_tT = arma::mat(T_full, n_vars).fill(NA_REAL);

  kf_obj.filter();
  a.rows(0, T_b-1) = kf_obj.a;
  a_tt.rows(0, T_b-1) = kf_obj.a_tt;

  if (T_b_ < T_full) {

    kf_obj.compact_to_companion(Lambda_);
    arma::mat y = y_.rows(T_b_, T_full-1);

    kf_end.set_pars(y,                                                   // y
                    kf_obj.Lambda_companion,                             // Z
                    arma::mat(size(y), fill::zeros),                     // c
                    arma::mat(n_vars, n_vars*(n_lags+1), fill::zeros),   // G
                    kf_obj.F_Phi,                                        // T
                    kf_obj.create_d(T_end),                              // d
                    kf_obj.Omega,                                        // H
                    kf_obj.a_Tb1,                                        // a1
                    kf_obj.P_Tb1,                                        // P1
                    arma::mat(1, n_vars, fill::zeros));                  // intercept

    kf_end.filter();

    // Fill in one-step predictions and filtered estimates
    a.rows(T_b, T_full - 1)    = companion_reshaper(kf_end.a,    n_m, n_q, T_end, n_lags);
    a_tt.rows(T_b, T_full - 1) = companion_reshaper(kf_end.a_tt, n_m, n_q, T_end, n_lags);

    kf_end.smoother(arma::mat(1, n_vars*(n_lags + 1), fill::zeros));
    alpha_t1 = companion_reshaper(kf_end.a_tT.row(0), n_m, n_q, 1, n_lags);
    r_T = (alpha_t1 - kf_obj.a_t1) * arma::pinv(kf_obj.P_t1);

    kf_obj.smoother(r_T);
    a_tT.rows(T_b, T_full - 1) = companion_reshaper(kf_end.a_tT, n_m, n_q, T_end, n_lags);
    Z_tT.rows(T_b, T_full - 1) = kf_end.a_tT.cols(0, n_vars - 1);


  } else {
    kf_obj.smoother(arma::mat(1, n_q*(n_lags + 1), fill::zeros));
  }

  a_tT.rows(0, T_b - 1) = kf_obj.a_tT;
  Z_tT(span(0, T_b - 1), span(0, n_m-1)) = y_Tb.cols(0, n_m - 1);
  Z_tT(span(0, T_b - 1), span(n_m, n_vars - 1)) = kf_obj.a_tT.cols(0, n_q - 1);

  return Rcpp::List::create(Rcpp::Named("a") = a,
                            Rcpp::Named("a_tt") = a_tt,
                            Rcpp::Named("a_tT") = a_tT,
                            Rcpp::Named("Z_tT") = Z_tT);
}

//' @describeIn kf_ragged Simulation smoother
//' @param Z1 initial values, with \code{n_lags} rows and same number of columns as \code{y_}
//' @return For \code{kf_sim_smooth}, a matrix with the draw from the posterior distribution.
// [[Rcpp::export]]
arma::mat kf_sim_smooth(arma::mat y_, arma::mat Phi_, arma::mat Sigma_, arma::mat Lambda_, arma::mat Z1_, int n_q_, unsigned int T_b_) {

  unsigned int n_vars, n_lags, n_m, n_q, T_full;
  n_vars  = y_.n_cols;
  n_lags  = (Phi_.n_cols-1)/n_vars;
  n_q     = n_q_;
  n_m     = n_vars - n_q;
  T_full  = y_.n_rows;

  arma::mat Sigma_chol = chol(Sigma_).t();

  // Initialize Z
  // Instead of usin n_vars*n_lags columns, use n_vars and then extract multiple rows
  arma::mat Z = arma::mat(T_full, n_vars).fill(NA_REAL);
  Z.rows(0, n_lags - 1) = Z1_;

  // Create y_sim
  arma::mat y_sim = arma::mat(T_full, n_vars).fill(NA_REAL);

  // Draw errors
  arma::mat epsilon = arma::mat(T_full - n_lags, n_vars);
  for (uword i = 0; i < n_vars; i++) {
    epsilon.col(i) = as<arma::vec>(rnorm(T_full - n_lags));
  }

  epsilon = epsilon * Sigma_chol.t();

  arma::rowvec Z_t1(n_vars*n_lags);
  arma::mat Phi_no_c = Phi_.cols(0, n_vars*n_lags-1).t();
  arma::rowvec Phi_c = arma::vectorise(Phi_.col(n_vars*n_lags), 1);
  arma::rowvec Z_t(n_vars);
  arma::rowvec y_all(n_vars);
  arma::uvec obs_vars;
  arma::uvec t_vec(1);
  arma::mat Lambda_t = Lambda_.t();
  arma::uword agg_length = Lambda_t.n_rows/Lambda_t.n_cols;
  arma::mat Z_mean = arma::mat(T_full, n_vars, fill::zeros);
  y_sim.rows(0, n_lags - 1) = y_.rows(0, n_lags - 1);

  for (uword t = n_lags; t < T_full; t++) {
    obs_vars = find_finite(y_.row(t));
    t_vec(0) = t;
    Z_t1 = arma::vectorise(arma::fliplr(Z.rows(t-n_lags, t-1).t())).t();
    Z_t = Z_t1 * Phi_no_c + Phi_c + epsilon.row(t - n_lags);
    Z.row(t) = Z_t;
    y_all.cols(0, n_m - 1) = Z_t.cols(0, n_m - 1);
    y_all.cols(n_m, n_vars - 1) = arma::join_rows(Z_t.cols(n_m, n_vars - 1), arma::vectorise(arma::fliplr(Z.submat(t-agg_length+1, n_m, t-1, n_vars - 1).t())).t()) * Lambda_t;
    y_sim(t_vec, obs_vars) = y_all.cols(obs_vars);
  }

  arma::mat Phi_diff = Phi_;
  Phi_diff.col(n_vars*n_lags) = arma::mat(n_vars, 1, fill::zeros);
  arma::mat y_diff = y_.rows(n_lags, T_full - 1) - y_sim.rows(n_lags, T_full - 1);
  arma::mat Z1_diff = arma::mat(arma::size(Z1_), fill::zeros);
  Rcpp::List smooth_diff = kf_ragged(y_diff, Phi_diff, Sigma_, Lambda_, Z1_diff, n_q_, T_b_ - n_lags);

  arma::mat Z_draw = Z.rows(n_lags, T_full - 1) + as<arma::mat>(smooth_diff["Z_tT"]);

  return Z_draw;
}


