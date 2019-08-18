#include "simsm_utils.h"

#ifndef MFBVAR_SIMSM_ADAPTIVE_SV_H
#define MFBVAR_SIMSM_ADAPTIVE_SV_H

inline arma::mat simsm_adaptive_sv(arma::mat y_, arma::mat Phi, arma::cube Sigma_chol,
                         arma::mat Lambda, arma::mat Z1, arma::uword n_q_, arma::uword T_b) {
  // intercept is first column
  arma::uword n_vars = y_.n_cols;
  arma::uword n_lags = Z1.n_rows;
  arma::uword n_T = y_.n_rows;
  arma::uword n_q = n_q_;
  arma::uword n_m = n_vars - n_q;

  ///////////////////////////////////////////////
  //               SIMULATING                  //
  ///////////////////////////////////////////////
  arma::mat Phi_no_c = Phi.cols(1, n_vars * n_lags);
  arma::mat Phi_c = Phi.col(0);
  arma::uvec quarterly_indexes = create_quarterly_indexes(Lambda, n_m, n_q, n_vars);

  // Generating errors
  arma::mat epsilon = arma::mat(n_T, n_vars, arma::fill::zeros);
  std::generate(epsilon.begin(), epsilon.end(), ::norm_rand);
  arma::mat Z_gen = arma::mat(n_T+1, n_vars*n_lags, arma::fill::zeros);
  arma::mat y_sim = arma::mat(n_T, n_vars).fill(NA_REAL);

  arma::uvec obs_vars;
  arma::uvec t_vec(1);
  arma::mat Zt;

  for (arma::uword i = 1; i <= n_T; i++) {
    epsilon.row(i-1) = epsilon.row(i-1) * Sigma_chol.slice(i-1).t();
  }

  create_sim(Z_gen, y_sim, Phi_no_c, Phi_c, epsilon, y_, Z1, Lambda,
             quarterly_indexes, n_vars, n_m, n_q, n_T, n_lags);

  ///////////////////////////////////////////////
  //                 COMPACT                   //
  ///////////////////////////////////////////////
  arma::mat y_orig = y_;
  y_ = y_ - y_sim;
  arma::mat Z1_orig = Z1;
  Z1 = arma::mat(arma::size(Z1), arma::fill::zeros);
  Phi.col(0) = arma::mat(n_vars, 1, arma::fill::zeros);
  Phi_c = Phi_c.fill(0.0);


  ///////////////////////////////////////////////
  //                 COMPACT                   //
  ///////////////////////////////////////////////
  arma::mat y_Tb = y_.rows(0, T_b - 1);
  arma::mat y = y_Tb;

  arma::mat Phi_mm(n_m, n_m*n_lags);
  arma::mat Phi_mq(n_m, n_q*n_lags);
  arma::mat Phi_qm(n_q, n_m*n_lags);
  arma::mat Phi_qq(n_q, n_q*n_lags);

  arma::mat Z = arma::mat(n_vars, n_q*(n_lags + 1), arma::fill::zeros);
  arma::mat Tt = arma::mat(n_q*(n_lags + 1), n_q*(n_lags + 1), arma::fill::zeros);

  arma::cube G = arma::cube(n_vars, n_vars, Sigma_chol.n_slices, arma::fill::zeros);
  arma::cube H = arma::cube(n_q*(n_lags + 1), n_vars, Sigma_chol.n_slices, arma::fill::zeros);
  arma::mat HT = arma::mat(n_q*(n_lags + 1), n_vars, arma::fill::zeros);

  arma::mat X = arma::mat(1, n_m*n_lags, arma::fill::zeros);
  arma::mat c = arma::mat(1, n_vars, arma::fill::zeros);
  arma::mat intercept = arma::mat(1, n_vars, arma::fill::zeros);

  arma::mat W = arma::mat(1, n_m*n_lags + 1, arma::fill::ones);
  arma::mat d = arma::mat(1, n_q*(n_lags + 1), arma::fill::zeros);
  arma::mat Beta_W = arma::mat(n_q, n_m*n_lags+1, arma::fill::zeros);

  arma::mat a1 = arma::mat(n_q*(n_lags+1), 1, arma::fill::zeros);

  create_matrices(Phi_mm, Phi_mq, Phi_qm, Phi_qq, Z, Tt, intercept, a1, W, d,
                  Beta_W, n_vars, n_m, n_q, n_lags, Lambda, Z1, Phi);

  G.tube(arma::span(0, n_m - 1), arma::span::all) = Sigma_chol(arma::span(0, n_m - 1), arma::span::all, arma::span(0, Sigma_chol.n_slices - 1));
  H.tube(arma::span(0, n_q - 1), arma::span::all) = Sigma_chol(arma::span(n_m, n_vars - 1), arma::span::all, arma::span(0, Sigma_chol.n_slices - 1));
  // Only if there are ragged edges do we need HT
  if (T_b < n_T) {
    HT(arma::span(0, n_q - 1), arma::span::all) = Sigma_chol.slice(T_b)(arma::span(n_m, n_vars - 1), arma::span::all);
  }

  a1 += d.t();
  arma::mat P1 = H.slice(0) * H.slice(0).t();
  ///////////////////////////////////////////////
  //             COMPACT FILTERING             //
  ///////////////////////////////////////////////
  arma::mat r_T;
  arma::mat a_tt_compact = arma::mat(n_T, n_q*(n_lags+1)).fill(NA_REAL);
  arma::mat a = arma::mat(n_T, n_q*(n_lags+1)).fill(NA_REAL);
  arma::mat a_tT = arma::mat(n_T, n_q*(n_lags+1)).fill(NA_REAL);
  arma::mat Z_tT = arma::mat(n_T, n_vars).fill(NA_REAL);

  arma::mat a_tt_y = y_;
  arma::mat a_tT_y = y_;
  arma::mat a_y = arma::mat(arma::size(y_)).fill(NA_REAL);

  arma::mat Gt, Ht, M_t, FF_inv_t, K_t, v_t, a_t1, P_t1, P_TT;
  arma::mat a_t = a1.t();
  arma::mat P_t = P1;
  Ht = H.slice(0);

  arma::mat v_FF_inv = arma::mat(n_T, n_vars, arma::fill::zeros);
  arma::mat v = arma::mat(n_T, n_vars, arma::fill::zeros);
  arma::mat L, N;
  arma::mat v_FF_inv_t, F_t, GZH;

  arma::cube L_compact = arma::cube(n_q*(n_lags+1), n_q*(n_lags+1), T_b).fill(NA_REAL);
  arma::cube N_compact = arma::cube(n_q*(n_lags+1), n_q*(n_lags+1), T_b+1).fill(NA_REAL);

  for (arma::uword t = 0; t < T_b; t++) {
    prepare_filtering_t(obs_vars, t_vec, X, W, c, d, t, y, Z1, Phi_mm, Beta_W,
                        n_m, n_q, n_lags);

    Zt = Z.rows(obs_vars);
    Gt = G.slice(t).rows(obs_vars);

    v_t = y.submat(t_vec, obs_vars) - a_t * Zt.t() - c.cols(obs_vars) - intercept.cols(obs_vars);
    M_t = P_t * Zt.t() + Ht * Gt.t();
    GZH = Gt + (Zt * Ht);
    F_t = Zt * M_t + Gt * GZH.t();
    FF_inv_t = arma::inv_sympd(arma::symmatu(F_t));
    v_FF_inv_t = v_t * FF_inv_t;
    v_FF_inv.submat(t_vec, obs_vars) = v_FF_inv_t;
    v.submat(t_vec, obs_vars) = v_t;

    K_t = Tt * M_t * FF_inv_t;
    L = Tt - K_t * Zt;
    N = P_t * L.t() - Ht * Gt.t() * K_t.t();
    a.row(t) = a_t;
    a_tt_compact.row(t) = a_t + v_FF_inv_t * M_t.t();
    a_y.row(t).cols(n_m, n_vars - 1) = a.row(t).cols(0, n_q - 1);
    a_tt_y.row(t).cols(n_m, n_vars - 1) = a_tt_compact.row(t).cols(0, n_q - 1);

    L_compact.slice(t) = L;
    N_compact.slice(t) = N;

    if (t < T_b - 1) {
      a_t = a_tt_compact.row(t) * Tt.t() + d;
      Ht = H.slice(t+1);
      P_t = Tt * N;
      P_t += Ht * Ht.t();
      P_t = arma::symmatu(P_t);
    } else {
      a_t1 = a_t * Tt.t() + d + v_t * K_t.t();
      if (T_b < n_T) {
        Ht = HT;
        P_t1 = Tt * N + Ht * Ht.t();
        P_t1 = arma::symmatu(P_t1);
        P_TT = P_t - M_t * FF_inv_t * M_t.t();
      }
    }

  }

  ///////////////////////////////////////////////
  //                ADAPTIVE                   //
  ///////////////////////////////////////////////
  arma::field<arma::mat> a_tt_out(n_T-T_b, 1);
  arma::field<arma::mat> L_store(n_T-T_b, 1);
  arma::field<arma::mat> N_store(n_T-T_b, 1);
  arma::field<arma::mat> ZFv(n_T-T_b, 1);

  arma::mat Phi_umum, Phi_omom, Phi_omu, W_intercept, Phi_uom, y_t, y_tpt, Phi_uu, y_tpt2, a_tt;
  arma::uvec obs_m, obs_q, non_obs_m, obs_m2, non_obs_m2;
  arma::uword n_ovars, n_oq, n_om2;
  arma::uword n_om = 0;

  if (T_b < n_T) {
    update_missing(y_t, obs_vars, obs_q, n_ovars, n_oq, obs_m, n_om, non_obs_m, obs_m2,
                   n_om2, non_obs_m2, y_tpt, y_tpt2, T_b, y_, n_vars, n_m, n_lags);

    Ht = H.slice(T_b);
    Ht = arma::mat((n_m-n_om+n_q)*(n_lags+1), n_vars, arma::fill::zeros);
    if (n_m > n_om) {
      Ht.rows(0, n_m - n_om - 1) = Sigma_chol.slice(T_b).rows(non_obs_m);
    }
    Ht.rows(n_m - n_om, n_m - n_om + n_q - 1) = Sigma_chol.slice(T_b).rows(n_m, n_m+n_q-1);

    a_tt = a_tt_compact.row(T_b-1);

    Phi_uom = create_Phi_uom(Phi, n_vars, n_q, n_m, n_om, n_om2, n_lags, non_obs_m, obs_m2);
    Phi_uu = create_Phi_uu(Phi, n_vars, n_q, n_m, n_om, n_om2, n_lags, non_obs_m, non_obs_m2);

    // Update Tt
    Tt = arma::mat((n_q+n_m-n_om)*(n_lags + 1), (n_q+n_m-n_om2)*(n_lags + 1));
    d = arma::mat(1, (n_m-n_om+n_q)*(n_lags + 1), arma::fill::zeros);
    create_Tt_d(Tt, d, Phi_uu, T_b-1, y_, n_m, n_q, n_om,
                n_om2, n_lags, obs_m2, non_obs_m, y_tpt2, Phi_uom);

    a_t = a_tt * Tt.t() + d;
    P_t = (Tt * P_TT) * Tt.t() + Ht * Ht.t();
    store_a(a_y, a_t, T_b, non_obs_m, n_m, n_om, n_q);

    for (arma::uword t = T_b; t < n_T; t++) {
      t_vec(0) = t;

      X = arma::mat(1, n_om*n_lags, arma::fill::ones);
      X.cols(0, n_om*n_lags - 1) = reshape(trans(flipud(y_tpt.rows(0, n_lags-1))), 1, n_lags*n_om);

      W_intercept = arma::mat(n_m-n_om+n_q, 1);
      W_intercept.rows(0, n_m-n_om-1) = intercept.cols(non_obs_m).t();
      W_intercept.rows(n_m-n_om, n_m-n_om+n_q-1) = Phi.col(0).rows(n_m, n_vars - 1);

      Phi_omom = create_Phi_omom(Phi, n_vars, n_om, n_om2, n_lags, obs_m, obs_m2);
      Phi_omu = create_Phi_omu(Phi, n_vars, n_q, n_m, n_om, n_om2, n_lags, non_obs_m, obs_m, obs_vars);

      c = arma::mat(1, n_ovars, arma::fill::zeros);
      c.cols(0, n_om - 1) = X * trans(Phi_omom);

      Zt = arma::mat(obs_vars.n_elem, (n_q+n_m-n_om)*(n_lags+1), arma::fill::zeros);
      create_Zt(Zt, Phi_omu, Lambda, n_ovars, n_m, n_om, n_om2, n_q, n_oq, n_lags, obs_q);

      Gt = G.slice(t).rows(obs_vars);
      v_t = y_.submat(t_vec, obs_vars) - a_t * Zt.t() - c - intercept.cols(obs_vars);
      M_t = P_t * Zt.t() + Ht * Gt.t();
      FF_inv_t = inv_sympd(symmatu(Zt * M_t + Gt * trans(Gt + Zt * Ht)));
      v_FF_inv_t = v_t * FF_inv_t;
      v_FF_inv.submat(t_vec, obs_vars) = v_FF_inv_t;
      v.submat(t_vec, obs_vars) = v_t;
      ZFv(t-T_b, 0) = v_FF_inv_t * Zt;

      if (t < n_T - 1) {
        update_missing(y_t, obs_vars, obs_q, n_ovars, n_oq, obs_m, n_om, non_obs_m, obs_m2,
                       n_om2, non_obs_m2, y_tpt, y_tpt2, t+1, y_, n_vars, n_m, n_lags);
      } else {
        obs_m2 = obs_m;
        n_om2 = n_om;
        non_obs_m2 = non_obs_m;
      }

      Phi_uom = create_Phi_uom(Phi, n_vars, n_q, n_m, n_om, n_om2, n_lags, non_obs_m, obs_m2);
      Phi_uu = create_Phi_uu(Phi, n_vars, n_q, n_m, n_om, n_om2, n_lags, non_obs_m, non_obs_m2);

      if (t < n_T - 1) {
        Tt = arma::mat((n_q+n_m-n_om)*(n_lags + 1), (n_q+n_m-n_om2)*(n_lags + 1));
        d = arma::mat(1, (n_m-n_om+n_q)*(n_lags + 1), arma::fill::zeros);
        create_Tt_d(Tt, d, Phi_uu, t, y_, n_m, n_q, n_om,
                    n_om2, n_lags, obs_m2, non_obs_m, y_tpt2, Phi_uom);
        K_t = Tt * M_t * FF_inv_t;
        L = Tt - K_t * Zt;
        N = P_t * L.t() - Ht * Gt.t() * K_t.t();

        L_store(t-T_b, 0) = L;
        N_store(t-T_b, 0) = N;
      }

      a_tt = a_t + v_FF_inv_t * M_t.t();
      a_tt_out(t-T_b, 0) = a_tt;
      store_a(a_tt_y, a_tt, t, non_obs_m2, n_m, n_om2, n_q);
      if (t < n_T - 1) {
        Ht = arma::mat((n_m-n_om+n_q)*(n_lags+1), n_vars, arma::fill::zeros);
        Ht.rows(0, n_m - n_om - 1) = Sigma_chol.slice(t+1).rows(non_obs_m);
        Ht.rows(n_m - n_om, n_m - n_om + n_q - 1) = Sigma_chol.slice(t+1).rows(n_m, n_m+n_q-1);

        a_t = a_tt * Tt.t() + d;
        P_t = Tt * N + Ht * Ht.t();
        P_t = arma::symmatu(P_t);

        store_a(a_y, a_t, t+1, non_obs_m, n_m, n_om, n_q);
      }
    }
  }
  ///////////////////////////////////////////////
  //               SMOOTHING                   //
  ///////////////////////////////////////////////
  arma::mat r = adaptive_to_compact_smoothing(a_tT_y, a_tt_y, a_tt, a_tt_compact, a_tt_out,
                                              N_store, L_store, ZFv, y_, a_t1, P_t1, n_vars,
                                              n_m, n_q, n_T, T_b, n_lags, n_om);
  compact_smoothing(a_tT_y, r, a_tt_compact, Z, y, N_compact, L_compact, v_FF_inv,
                    T_b, n_vars, n_m, n_q);
  arma::mat Z_rand = Z_gen.rows(1, n_T).cols(0, n_vars - 1);
  arma::mat Z_draw = a_tT_y + Z_rand;
  return Z_draw;
}

#endif
