#include "mfbvar.h"

arma::mat companion_reshaper(arma::mat obj_, unsigned int n_m_, unsigned int n_q_, unsigned int n_T_, unsigned int n_lags_) {
  unsigned int n_vars_ = n_m_ + n_q_;
  arma::mat temp = arma::reshape(obj_, n_vars_*n_T_, n_lags_ + 1);
  return arma::reshape(temp.rows(n_m_*n_T_, n_vars_*n_T_-1), n_T_, n_q_ * (n_lags_ + 1));
}

arma::mat create_Phi_uu(const arma::mat &Phi, arma::uword n_vars, arma::uword n_q, arma::uword n_m, arma::uword n_om,
                        arma::uword n_om2, arma::uword n_lags, arma::uvec non_obs_m, arma::uvec non_obs_m2) {
  arma::mat Phi_uu = arma::mat(n_q+n_m-n_om, (n_q+n_m-n_om2)*n_lags, arma::fill::zeros);
  // Monthly part
  for (arma::uword i = 0; i < (n_m-n_om)+n_q; i++) {
    for (arma::uword j = 0; j < n_lags; j++) {
      for (arma::uword k = 0; k < (n_m-n_om2); k++) {
        if (i < (n_m-n_om)) {
          Phi_uu(i, k + (n_m-n_om2+n_q)*j) = Phi(non_obs_m(i), non_obs_m2(k) + j*n_vars + 1);
        } else {
          Phi_uu(i, k + (n_m-n_om2+n_q)*j) = Phi(i+n_om, non_obs_m2(k) + j*n_vars + 1);
        }

      }
    }
  }
  // Quarterly part
  for (arma::uword i = 0; i < (n_m-n_om)+n_q; i++) {
    for (arma::uword j = 0; j < n_lags; j++) {
      for (arma::uword k = 0; k < n_q; k++) {
        if (i < (n_m-n_om)) {
          Phi_uu(i, k + (n_m-n_om2) + (n_m-n_om2+n_q)*j) = Phi(non_obs_m(i), n_m + k + j*n_vars + 1);
        } else {
          Phi_uu(i, k + (n_m-n_om2) + (n_m-n_om2+n_q)*j) = Phi(i+n_om, n_m + k + j*n_vars + 1);
        }
      }
    }
  }
  return Phi_uu;
}

arma::mat create_Phi_uom(const arma::mat &Phi, arma::uword n_vars, arma::uword n_q, arma::uword n_m, arma::uword n_om,
                         arma::uword n_om2, arma::uword n_lags, arma::uvec non_obs_m, arma::uvec obs_m2) {
  arma::mat Phi_uom = arma::mat(n_m-n_om+n_q, n_om2*n_lags, arma::fill::zeros);
  // Monthly part
  for (arma::uword i = 0; i < (n_m-n_om); i++) {
    for (arma::uword j = 0; j < n_lags; j++) {
      for (arma::uword k = 0; k < n_om2; k++) {
        Phi_uom(i, k+j*n_om2) = Phi(non_obs_m(i), obs_m2(k)+j*n_vars+1);
      }
    }
  }
  if (n_q > 0) {
    // Quarterly part
    for (arma::uword i = 0; i < n_q; i++) {
      for (arma::uword j = 0; j < n_lags; j++) {
        for (arma::uword k = 0; k < n_om2; k++) {
          Phi_uom(i+n_m-n_om, k+j*n_om2) = Phi(n_m+i, obs_m2(k)+j*n_vars+1);
        }
      }
    }
  }
  return Phi_uom;
}

arma::mat create_Phi_omu(const arma::mat &Phi, arma::uword n_vars, arma::uword n_q, arma::uword n_m, arma::uword n_om,
                         arma::uword n_om2, arma::uword n_lags, arma::uvec non_obs_m, arma::uvec obs_m, arma::uvec obs_vars) {
  arma::mat Phi_omu = arma::mat(obs_vars.n_elem, (n_q+n_m-n_om)*n_lags, arma::fill::zeros);
  // Monthly part
  for (arma::uword i = 0; i < n_om; i++) {
    for (arma::uword j = 0; j < n_lags; j++) {
      for (arma::uword k = 0; k < (n_m-n_om); k++) {
        Phi_omu(i, k+j*(n_m-n_om+n_q)) = Phi(obs_m(i), non_obs_m(k)+j*n_vars+1);
      }
    }
  }
  // Quarterly part
  for (arma::uword i = 0; i < n_om; i++) {
    for (arma::uword j = 0; j < n_lags; j++) {
      for (arma::uword k = 0; k < n_q; k++) {
        Phi_omu(i, (n_m-n_om)+k+j*(n_m-n_om+n_q)) = Phi(obs_m(i), n_m+k+j*n_vars+1);
      }
    }
  }
  return Phi_omu;
}

arma::mat create_Phi_omom(const arma::mat &Phi, arma::uword n_vars, arma::uword n_om, arma::uword n_om2,
                          arma::uword n_lags, arma::uvec obs_m, arma::uvec obs_m2) {
  // arma::mat Phi_omom = arma::mat(n_om, n_om2*n_lags);
  // for (arma::uword i = 0; i < n_om; i++) {
    //   for (arma::uword j = 0; j < n_lags; j++) {
      //     for (arma::uword k = 0; k < n_om2; k++) {
        //       Phi_omom(i, k+j*n_om2) = Phi(obs_m(i), obs_m2(k)+j*n_vars+1);
        //     }
      //   }
    // }
  arma::mat Phi_omom = arma::mat(n_om, n_om*n_lags, arma::fill::zeros);
  for (arma::uword i = 0; i < n_om; i++) {
    for (arma::uword j = 0; j < n_lags; j++) {
      for (arma::uword k = 0; k < n_om; k++) {
        Phi_omom(i, k+j*n_om) = Phi(obs_m(i), obs_m(k)+j*n_vars+1);
      }
    }
  }
  return Phi_omom;
}

void create_Tt_d(arma::mat & Tt, arma::mat & d, const arma::mat & Phi_uu, arma::uword t,
                 const arma::mat & y_, arma::uword n_m, arma::uword n_q,
                 arma::uword n_om, arma::uword n_om2, arma::uword n_lags,
                 arma::uvec obs_m2, arma::uvec non_obs_m, const arma::mat & y_tpt2,
                 const arma::mat & Phi_uom) {
  Tt.fill(0);
  if (!Phi_uu.is_empty()) {
    Tt.submat(0, 0, n_q+n_m-n_om-1, (n_q+n_m-n_om2)*n_lags-1) = Phi_uu;
  }

  int counter;
  for (arma::uword j = 0; j < n_lags; j++) {
    counter = 0;
    for (arma::uword i = 0; i < (n_m - n_om); i++) {
      if (std::find(obs_m2.begin(), obs_m2.end(), non_obs_m(i)) == obs_m2.end()) {
        Tt(i+(n_q+n_m-n_om)+j*(n_q+n_m-n_om), counter+j*(n_q+n_m-n_om2)) = 1;
        counter += 1;
      } else {
        d(i+(n_q+n_m-n_om)+j*(n_q+n_m-n_om)) = y_(t-j, non_obs_m(i));
      }
    }
    for (arma::uword i = 0; i < n_q; i++) {
      Tt(i+(n_m-n_om) + (n_q+n_m-n_om)+j*(n_q+n_m-n_om), i+(n_m-n_om2)+j*(n_q+n_m-n_om2)) = 1;
    }
  }

  arma::mat W = arma::mat(1, n_om2*n_lags+1, arma::fill::ones);
  if (n_om2*n_lags >= 1) {
    W.cols(0, n_om2*n_lags - 1) = reshape(trans(flipud(y_tpt2.rows(1, n_lags))), 1, n_lags*n_om2);
  }
  if (n_m-n_om+n_q >= 1) {
    d.cols(0, n_m-n_om+n_q-1) = W.cols(0, n_om2*n_lags - 1) * arma::trans(Phi_uom);
  }
}

void create_Zt(arma::mat & Zt, const arma::mat & Phi_omu, const arma::mat & Lambda,
               arma::uword n_ovars, arma::uword n_m, arma::uword n_om, arma::uword n_om2,
               arma::uword n_q, arma::uword n_oq, arma::uword n_lags, const arma::uvec & obs_q) {
  if (!Zt.is_empty()) {
    Zt(0, n_q+n_m-n_om, arma::size(n_ovars, (n_q+n_m-n_om)*n_lags)) = Phi_omu;
    if (n_oq != 0) {
      for (arma::uword i = 0; i < n_oq; i++) {
        for (arma::uword j = 0; j < Lambda.n_cols/n_q; j++) {
          Zt(i+n_om, obs_q(i)+n_m-n_om+(n_m-n_om+n_q)*j) = Lambda(obs_q(i), obs_q(i)+j*n_q);
        }
      }
    }
  }
}

void update_missing(arma::mat & y_t, arma::uvec & obs_vars, arma::uvec & obs_q,
                    arma::uword & n_ovars, arma::uword & n_oq,
                    arma::uvec & obs_m, arma::uword & n_om, arma::uvec & non_obs_m,
                    arma::uvec & obs_m2, arma::uword & n_om2, arma::uvec & non_obs_m2,
                    arma::mat & y_tpt, arma::mat & y_tpt2,
                    arma::uword t, const arma::mat & y_, arma::uword n_vars, arma::uword n_m, arma::uword n_lags) {
  y_t = y_.row(t);
  obs_vars = find_finite(y_t);
  n_ovars = obs_vars.n_elem;
  if (n_m <= n_vars - 1) {
    obs_q = find_finite(y_t.cols(n_m, n_vars - 1));
    n_oq = obs_q.n_elem;
  }

  obs_m = find_finite(y_t.cols(0, n_m - 1));
  n_om = obs_m.n_elem;
  non_obs_m = find_nonfinite(y_t.cols(0, n_m - 1));

  obs_m2 = find_finite(y_.row(t-1).cols(0, n_m - 1));
  n_om2 = obs_m2.n_elem;
  non_obs_m2 = find_nonfinite(y_.row(t-1).cols(0, n_m - 1));

  y_tpt = y_.rows(t-n_lags, t);
  y_tpt2 = y_.rows(t-n_lags-1, t-1);
  y_tpt2 = y_tpt2.cols(obs_m2);
  y_tpt = y_tpt.cols(obs_m);
}

void store_a(arma::mat & a, arma::mat a_t, arma::uword t, arma::uvec non_obs_m,
             arma::uword n_m, arma::uword n_om, arma::uword n_q) {
  for (arma::uword i = 0; i < (n_m - n_om); i++) {
    a(t, non_obs_m(i)) = a_t(i);
  }
  for (arma::uword i = 0; i < n_q; i++) {
    a(t, n_m + i) = a_t(i+(n_m-n_om));
  }
}

void compact_smoothing(arma::mat & a_tT_y, arma::mat & r,
                       const arma::mat & a_tt_compact, const arma::mat & Z,
                       const arma::mat & y, const arma::cube & N_compact,
                       const arma::cube & L_compact, const arma::mat & v_FF_inv,
                       arma::uword T_b, arma::uword n_vars, arma::uword n_m,
                       arma::uword n_q) {
  arma::mat a_tT, Zt;
  arma::uvec obs_vars;
  arma::uvec t_vec(1);
  for (arma::uword t = T_b - 1; t >= 1; t--) {
    obs_vars = arma::find_finite(y.row(t));
    t_vec(0) = t;
    Zt = Z.rows(obs_vars);
    a_tT = a_tt_compact.row(t) + r * N_compact.slice(t).t();
    a_tT_y.row(t).cols(n_m, n_vars - 1) = a_tT.cols(0, n_q-1);

    r = v_FF_inv.submat(t_vec, obs_vars) * Zt + r * L_compact.slice(t);
  }

  obs_vars = find_finite(y.row(0));
  Zt = Z.rows(obs_vars);
  a_tT = a_tt_compact.row(0) + r * N_compact.slice(0).t();
  a_tT_y.row(0).cols(n_m, n_vars - 1) = a_tT.cols(0, n_q-1);
}

arma::mat adaptive_to_compact_smoothing(arma::mat & a_tT_y,
                                        const arma::mat &a_tt_y,
                                        const arma::mat & a_tt,
                                        const arma::mat & a_tt_compact,
                                        const arma::field<arma::mat> & a_tt_out,
                                        const arma::field<arma::mat> & N_store,
                                        const arma::field<arma::mat> & L_store,
                                        const arma::field<arma::mat> & ZFv,
                                        const arma::mat y_, const arma::mat & a_t1,
                                        const arma::mat & P_t1, arma::uword n_vars,
                                        arma::uword n_m, arma::uword n_q,
                                        arma::uword n_T, arma::uword T_b, arma::uword n_lags,
                                        arma::uword n_om) {
  arma::mat r, a_tT, y_t;
  arma::uvec obs_m, non_obs_m;
  if (n_T - 1 > T_b) {
    r = ZFv(n_T-T_b-1, 0);
    for (arma::uword t = n_T - T_b - 1; t >= 1; t--) {
      a_tT = a_tt_out(t-1, 0) + r * N_store(t-1, 0).t();
      r = ZFv(t-1, 0) + r * L_store(t-1, 0);

      y_t = y_.row(t+T_b-1);
      obs_m = find_finite(y_t.cols(0, n_m - 1));
      n_om = obs_m.n_elem;
      non_obs_m = find_nonfinite(y_t.cols(0, n_m - 1));

      store_a(a_tT_y, a_tT, t+T_b-1, non_obs_m, n_m, n_om, n_q);
    }
  } else {
    a_tT = a_tt;
  }
  arma::mat alpha_t1 = arma::mat(1, n_q*(n_lags+1));
  if (T_b < n_T) {
    a_tT_y.row(n_T-1) = a_tt_y.row(n_T-1);
    for (arma::uword i = 0; i < n_lags + 1; i++) {
      for (arma::uword j = 0; j < n_q; j++) {
        alpha_t1(j + i*n_q) = a_tT(j + n_m - n_om + i*(n_q+n_m-n_om));
      }
    }
    if (alpha_t1.is_empty() || a_t1.is_empty()) {
      r = (alpha_t1 - a_t1) * arma::pinv(P_t1);
    }
  } else {
    r = arma::mat(1, n_q*(n_lags+1), arma::fill::zeros);
    a_tT_y.row(n_T-1).cols(n_m, n_vars - 1) = a_tt_compact.row(n_T-1).cols(0, n_q-1);
  }
  return r;
}


arma::uvec create_quarterly_indexes(const arma::mat & Lambda, arma::uword n_m,
                                    arma::uword n_q, arma::uword n_vars){
  arma::uvec quarterly_indexes = arma::uvec(Lambda.n_cols);
  arma::uword count = 0;
  for (arma::uword i = 0; i < Lambda.n_cols/n_q; i++) {
    for (arma::uword j = 0; j < n_q; j++) {
      quarterly_indexes(count) = n_vars*i + n_m + j;
      count += 1;
    }
  }
  return quarterly_indexes;
}

void create_sim(arma::mat & Z_gen, arma::mat & y_sim,
                const arma::mat & Phi_no_c, const arma::mat & Phi_c,
                const arma::mat & epsilon, const arma::mat & y_,
                const arma::mat & Z1, const arma::mat & Lambda,
                const arma::uvec quarterly_indexes,
                arma::uword n_vars, arma::uword n_m, arma::uword n_q,
                arma::uword n_T, arma::uword n_lags) {
  Z_gen.row(0) = arma::reshape(arma::trans(arma::flipud(Z1)), 1, n_lags*n_vars); // Initial value
  arma::mat Z_gen_t = arma::mat(1, n_vars * n_lags, arma::fill::zeros);
  arma::uvec obs_vars;
  arma::uvec t_vec(1);
  for (arma::uword i = 1; i <= n_T; i++) {
    Z_gen_t.cols(n_vars, n_vars*n_lags - 1) = Z_gen.submat(i-1, 0, i-1, n_vars*(n_lags-1) - 1);
    Z_gen_t.cols(0, n_vars - 1) = Z_gen.row(i-1) * Phi_no_c.t() + Phi_c.t() + epsilon.row(i-1);
    Z_gen.row(i) = Z_gen_t;

    obs_vars = find_nonfinite(y_.row(i-1));
    t_vec(0) = i-1;
    if (n_m > 0) {
      y_sim.submat(i-1, 0, i-1, n_m - 1) = Z_gen_t.cols(0, n_m - 1);
    }
    if (n_q > 0) {
      y_sim.submat(i-1, n_m, i-1, n_vars - 1) = Z_gen_t(quarterly_indexes).t() * Lambda.t();
    }
    y_sim.submat(t_vec, obs_vars) = arma::mat(1, obs_vars.n_elem).fill(NA_REAL);
  }
}

void create_matrices(arma::mat & Phi_mm, arma::mat & Phi_mq, arma::mat & Phi_qm,
                     arma::mat & Phi_qq, arma::mat & Z, arma::mat & Tt, arma::mat & intercept,
                     arma::mat & a1, arma::mat & W, arma::mat & d, arma::mat & Beta_W,
                     arma::uword n_vars, arma::uword n_m, arma::uword n_q,
                     arma::uword n_lags, const arma::mat & Lambda,
                     const arma::mat & Z1, const arma::mat & Phi) {
  for (unsigned int i = 0; i < n_lags; i++) {
    Rcpp::Rcout << "create mat mm" << std::endl;
    if (n_m > 0) {
      Phi_mm.cols(i*n_m, (i+1)*n_m - 1) = Phi.submat(0, i*n_vars+1, n_m - 1, i*n_vars + n_m);
    }
    if (n_m > 0 && n_q > 0) {
      Rcpp::Rcout << "create mat mq" << std::endl;
      Phi_mq.cols(i*n_q, (i+1)*n_q - 1) = Phi.submat(0, i*n_vars+n_m+1, n_m - 1, (i+1)*n_vars);
      Rcpp::Rcout << "create mat qm" << std::endl;
      Phi_qm.cols(i*n_m, (i+1)*n_m - 1) = Phi.submat(n_m, i*n_vars+1, n_vars - 1, i*n_vars + n_m);
    }
    if (n_q > 0) {
      Rcpp::Rcout << "create mat qq" << std::endl;
      Phi_qq.cols(i*n_q, (i+1)*n_q - 1) = Phi.submat(n_m, i*n_vars+n_m+1, n_vars - 1, (i+1)*n_vars);
    }
  }

  if (n_m > 0 && n_q > 0) {
    Z.submat(0, n_q, n_m - 1, n_q*(n_lags + 1) - 1) = Phi_mq;
  }
  if (n_m > 0) {
    intercept.cols(0, n_m - 1) = arma::trans(Phi.submat(0, 0, n_m - 1, 0));
    W.cols(0, n_m*n_lags - 1) = reshape(trans(flipud(Z1.cols(0, n_m-1))), 1, n_lags*n_m); //X.row(0)
    if (n_q > 0) {
      Beta_W = join_rows(Phi_qm, Phi.submat(n_m, 0, n_vars-1, 0));
      d.cols(0, n_q - 1) = W * trans(Beta_W);
      Z.submat(n_m, 0, n_vars - 1, Lambda.n_cols - 1) = Lambda;
    }
  }

  if (n_q > 0) {
    Tt.submat(0, 0, n_q - 1, n_q*n_lags - 1) = Phi_qq;
    Tt.submat(n_q, 0, n_q*(n_lags+1)-1, n_q*n_lags - 1) = arma::eye(n_lags*n_q, n_lags*n_q);
    a1.rows(0, n_q*n_lags - 1) = reshape(trans(flipud(Z1.cols(n_m, n_vars - 1))), 1, n_lags*n_q).t();
    a1 = Tt * a1;
  }

}

void prepare_filtering_t(arma::uvec & obs_vars, arma::uvec & t_vec, arma::mat & X,
                         arma::mat & W, arma::mat & c, arma::mat & d, arma::uword t,
                         const arma::mat & y, const arma::mat & Z1,
                         const arma::mat & Phi_mm, const arma::mat & Beta_W,
                         arma::uword n_m, arma::uword n_q, arma::uword n_lags) {
  obs_vars = find_finite(y.row(t));
  t_vec(0) = t;

  // Exogenous part
  X = W.cols(0, n_m*n_lags - 1);
  if (t >= n_lags - 1) {
    W.cols(0, n_m*n_lags - 1) = arma::reshape(arma::trans(arma::flipud(y(arma::span(t-n_lags+1, t), arma::span(0, n_m - 1)))), 1, n_lags*n_m);
  } else {
    W.cols(0, (t+1)*n_m - 1) = arma::reshape(arma::trans(arma::flipud(y(arma::span(0, t), arma::span(0, n_m - 1)))), 1, (t+1)*n_m);
    W.cols((t+1)*n_m, n_lags*n_m - 1) = reshape(trans(flipud(Z1(arma::span(t+1, n_lags - 1), arma::span(0, n_m - 1)))), 1, (n_lags-t-1)*n_m);
  }

  c.cols(0, n_m - 1) = X * trans(Phi_mm);
  d.cols(0, n_q - 1) = W * trans(Beta_W);
}
