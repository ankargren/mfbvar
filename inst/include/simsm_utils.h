#ifndef MFBVAR_SIMSM_UTILS_H
#define MFBVAR_SIMSM_UTILS_H

arma::mat create_Phi_uu(const arma::mat &Phi, arma::uword n_vars, arma::uword n_q, arma::uword n_m, arma::uword n_om,
                        arma::uword n_om2, arma::uword n_lags, arma::uvec non_obs_m, arma::uvec non_obs_m2);
arma::mat create_Phi_uom(const arma::mat &Phi, arma::uword n_vars, arma::uword n_q, arma::uword n_m, arma::uword n_om,
                         arma::uword n_om2, arma::uword n_lags, arma::uvec non_obs_m, arma::uvec obs_m2);
arma::mat create_Phi_omu(const arma::mat &Phi, arma::uword n_vars, arma::uword n_q, arma::uword n_m, arma::uword n_om,
                         arma::uword n_om2, arma::uword n_lags, arma::uvec non_obs_m, arma::uvec obs_m, arma::uvec obs_vars);
arma::mat create_Phi_omom(const arma::mat &Phi, arma::uword n_vars, arma::uword n_om, arma::uword n_om2,
                          arma::uword n_lags, arma::uvec obs_m, arma::uvec obs_m2);
void create_Tt_d(arma::mat & Tt, arma::mat & d, const arma::mat & Phi_uu, arma::uword t,
                 const arma::mat & y_, arma::uword n_m, arma::uword n_q,
                 arma::uword n_om, arma::uword n_om2, arma::uword n_lags,
                 arma::uvec obs_m2, arma::uvec non_obs_m, const arma::mat & y_tpt2,
                 const arma::mat & Phi_uom);
void create_Zt(arma::mat & Zt, const arma::mat & Phi_omu, const arma::mat & Lambda,
               arma::uword n_ovars, arma::uword n_m, arma::uword n_om, arma::uword n_om2,
               arma::uword n_q, arma::uword n_oq, arma::uword n_lags, const arma::uvec & obs_q);
void update_missing(arma::mat & y_t, arma::uvec & obs_vars, arma::uvec & obs_q,
                    arma::uword & n_ovars, arma::uword & n_oq,
                    arma::uvec & obs_m, arma::uword & n_om, arma::uvec & non_obs_m,
                    arma::uvec & obs_m2, arma::uword & n_om2, arma::uvec & non_obs_m2,
                    arma::mat & y_tpt, arma::mat & y_tpt2,
                    arma::uword t, const arma::mat & y_, arma::uword n_vars, arma::uword n_m, arma::uword n_lags);
void store_a(arma::mat & a, arma::mat a_t, arma::uword t, arma::uvec non_obs_m,
             arma::uword n_m, arma::uword n_om, arma::uword n_q);
arma::mat companion_reshaper(arma::mat obj_, unsigned int n_m_, unsigned int n_q_, unsigned int n_T_, unsigned int n_lags_);
void compact_smoothing(arma::mat & a_tT_y, arma::mat & r, const arma::mat & a_tt_compact, const arma::mat & Z,
                       const arma::mat & y, const arma::cube & N_compact,
                       const arma::cube & L_compact, const arma::mat & v_FF_inv,
                       arma::uword T_b, arma::uword n_vars, arma::uword n_m,
                       arma::uword n_q);
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
                                        arma::uword n_om);
arma::uvec create_quarterly_indexes(const arma::mat & Lambda, arma::uword n_m,
                                    arma::uword n_q, arma::uword n_vars);
void create_sim(arma::mat & Z_gen, arma::mat & y_sim,
                const arma::mat & Phi_no_c, const arma::mat & Phi_c,
                const arma::mat & epsilon, const arma::mat & y_,
                const arma::mat & Z1, const arma::mat & Lambda,
                const arma::uvec quarterly_indexes,
                arma::uword n_vars, arma::uword n_m, arma::uword n_q,
                arma::uword n_T, arma::uword n_lags);
void create_matrices(arma::mat & Phi_mm, arma::mat & Phi_mq, arma::mat & Phi_qm,
                     arma::mat & Phi_qq, arma::mat & Z, arma::mat & Tt, arma::mat & intercept,
                     arma::mat & a1, arma::mat & W, arma::mat & d, arma::mat & Beta_W,
                     arma::uword n_vars, arma::uword n_m, arma::uword n_q,
                     arma::uword n_lags, const arma::mat & Lambda,
                     const arma::mat & Z1, const arma::mat & Phi);
void prepare_filtering_t(arma::uvec & obs_vars, arma::uvec & t_vec, arma::mat & X,
                         arma::mat & W, arma::mat & c, arma::mat & d, arma::uword t,
                         const arma::mat & y, const arma::mat & Z1,
                         const arma::mat & Phi_mm, const arma::mat & Beta_W,
                         arma::uword n_m, arma::uword n_q, arma::uword n_lags);
#endif
