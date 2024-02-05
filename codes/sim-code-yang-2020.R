## simulation with analytical standard errors

yang_sim <- function(k) {
  
  sample_prob <- pop_data[sample(1:N, n),]
  sample_prob$w <- N/n
  ## calculate quantiles
  sample_p1 <- joint_calib_create_matrix(sample_prob[, c("x1", "x2")], 
                                         N, 
                                         pop_quantiles = list(x1=quantile(sample_prob$x1, p_quantiles1), 
                                                              x2=quantile(sample_prob$x2, p_quantiles1)))
  
  sample_p2 <- joint_calib_create_matrix(sample_prob[, c("x1", "x2")], 
                                         N, 
                                         pop_quantiles = list(x1=quantile(sample_prob$x1, p_quantiles2), 
                                                              x2=quantile(sample_prob$x2, p_quantiles2)))
   
  colnames(sample_p1) <- c(paste0("x1_", p_quantiles1), paste0("x2_", p_quantiles1))
  colnames(sample_p2) <- c(paste0("x1_", p_quantiles2), paste0("x2_", p_quantiles2))
  
  sample_prob <- cbind(sample_prob, sample_p1, sample_p2)
  sample_prob_svy <- svydesign(ids=~1, weights = ~w, data = sample_prob)
  
  q1_est <- svyquantile( ~ x1 + x2, sample_prob_svy, p_quantiles1)
  q2_est <- svyquantile( ~ x1 + x2, sample_prob_svy, p_quantiles2)
  
  x_totals <- svytotal( ~ x1 + x2, sample_prob_svy)
  
  sample_bd1 <- pop_data[rbinom(N,1,pop_data$p1)==1, ]
  sample_bd1$w_naive <- N/nrow(sample_bd1)
  sample_bd2 <- pop_data[rbinom(N,1,pop_data$p2)==1, ]
  sample_bd2$w_naive <- N/nrow(sample_bd2)
  
  
  ##########################################
  # linear inclusion --------------------------------------------------------
  ##########################################
  ## calibration  --------------------------------------------------------
  ## calibration with totals only (CAL) 
  bd1_w_cal_0 <- calib(Xs = with(sample_bd1, cbind(1, x1,x2)), 
                       d = sample_bd1$w_naive, 
                       total = c(N, x_totals), method = "raking")
  
  bd1_w_el_0 <- calib_el(X = with(sample_bd1, cbind(1, x1,x2)), 
                         d = sample_bd1$w_naive, 
                         totals = c(N, x_totals))
  
  ## calibration with quantiles only (QCAL1A)
  bd1_w_cal_1a <- joint_calib(formula_quantiles = ~ x1 + x2, data = sample_bd1, dweights = sample_bd1$w_naive,
                              N = N,  pop_quantiles = q1_est, method = "raking", backend = "sampling")
  
  ## calibration with quantiles and totals (QCAL1B)
  bd1_w_cal_1b <- joint_calib(formula_quantiles = ~ x1 + x2, formula_totals = ~ x1 + x2, data = sample_bd1,
                              dweights = sample_bd1$w_naive,  N = N,  pop_totals = x_totals,
                              pop_quantiles = q1_est, method = "raking", backend = "sampling")
  
  ## calibration with quantiles only (QCAL2A)
  bd1_w_cal_2a <- joint_calib(formula_quantiles = ~ x1 + x2, data = sample_bd1, dweights = sample_bd1$w_naive,
                              N = N,  pop_quantiles = q2_est, method = "raking", backend = "sampling")
  
  ## calibration with quantiles and totals (QCAL2B)
  bd1_w_cal_2b <- joint_calib(formula_quantiles = ~ x1 + x2, formula_totals = ~ x1 + x2, data = sample_bd1,
                              dweights = sample_bd1$w_naive,  N = N,  pop_totals = x_totals,
                              pop_quantiles = q2_est, method = "raking", backend = "sampling")
  
  ## add weights
  sample_bd1[ , ":="(w_cal_0  = w_naive*bd1_w_cal_0,    w_el_0   = w_naive*bd1_w_el_0,
                     w_cal_1a = w_naive*bd1_w_cal_1a$g, w_cal_1b = w_naive*bd1_w_cal_1b$g, 
                     w_cal_2a = w_naive*bd1_w_cal_2a$g, w_cal_2b = w_naive*bd1_w_cal_2b$g)]
  
  
  ## add XQ variables to sample_bd1
  bd1_p1x <- bd1_w_cal_1a$Xs[, -1]
  colnames(bd1_p1x) <- colnames(sample_p1)
  bd1_p2x <- bd1_w_cal_2a$Xs[, -1]
  colnames(bd1_p2x) <- colnames(sample_p2)
  sample_bd1 <- cbind(sample_bd1, bd1_p1x, bd1_p2x[, setdiff(colnames(bd1_p2x), colnames(bd1_p1x))])
  
  ## IPW  --------------------------------------------------------
  
  bd1_ipw_0 <- nonprob(selection = ~ x1 + x2, 
                       target = ~ y11 + y12 + y21 + y22, svydesign = sample_prob_svy,
                       data = sample_bd1,  
                       se = TRUE,
                       control_selection = controlSel(est_method_sel = "gee", h = 1))

  bd1_ipw_1a <- nonprob(selection = formula_xs_p1a,
                        target = ~ y11 + y12 + y21 + y22,
                        svydesign = sample_prob_svy,
                        data = sample_bd1,
                        se = TRUE,
                        start_selection = rep(0, ncol(bd1_w_cal_1a$Xs)),
                        control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd1_ipw_1b <- nonprob(selection = formula_xs_p1b,
                        target = ~ y11 + y12 + y21 + y22, svydesign = sample_prob_svy,
                        data = sample_bd1,
                        se = TRUE,
                        start_selection = rep(0, ncol(bd1_w_cal_1b$Xs)),
                        control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  
  bd1_ipw_2a <- nonprob(selection = formula_xs_p2a,
                        target = ~ y11 + y12 + y21 + y22, svydesign = sample_prob_svy,
                        data = sample_bd1,
                        se = TRUE,
                        start_selection = rep(0, ncol(bd1_w_cal_2a$Xs)),
                        control_selection = controlSel(est_method_sel = "gee", h = 1))


  bd1_ipw_2b <- nonprobsvy:::nonprob(selection = formula_xs_p2b,
                                     target = ~ y11 + y12 + y21 + y22, svydesign = sample_prob_svy,
                                     data = sample_bd1,
                                     se = TRUE,
                                     start_selection = rep(0, ncol(bd1_w_cal_2b$Xs)),
                                     control_selection = controlSel(est_method_sel = "gee",h = 1))

  
  ## mass imputation ---------------------------------------------------------
  ## mi for y11 (nn)
  bd1_mi_y11 <- nonprob(outcome = y11 ~ x1 + x2, svydesign = sample_prob_svy, data = sample_bd1, 
                        method_outcome = "nn", se = TRUE)
  bd1_mi_y12 <- nonprob(outcome = y12 ~ x1 + x2, svydesign = sample_prob_svy, data = sample_bd1, 
                        method_outcome = "nn", se = TRUE)
  bd1_mi_y21 <- nonprob(outcome = y21 ~ x1 + x2, svydesign = sample_prob_svy, data = sample_bd1, 
                        method_outcome = "nn", se = TRUE)
  bd1_mi_y22 <- nonprob(outcome = y22 ~ x1 + x2, svydesign = sample_prob_svy, data = sample_bd1, 
                        method_outcome = "nn", se = TRUE)
  
  ## mi for y11 (glm)
  bd1_mi_y11_glm <- nonprob(outcome = y11 ~ x1 + x2, svydesign = sample_prob_svy, data = sample_bd1, 
                            method_outcome = "glm", se = TRUE)
  bd1_mi_y12_glm <- nonprob(outcome = y12 ~ x1 + x2, svydesign = sample_prob_svy, data = sample_bd1, 
                            method_outcome = "glm", se = TRUE)
  bd1_mi_y21_glm <- nonprob(outcome = y21 ~ x1 + x2, svydesign = sample_prob_svy, data = sample_bd1, 
                            method_outcome = "glm", family_outcome = "binomial", se = TRUE)
  bd1_mi_y22_glm <- nonprob(outcome = y22 ~ x1 + x2, svydesign = sample_prob_svy, data = sample_bd1, 
                            method_outcome = "glm", family_outcome = "binomial", se = TRUE)
  
  ## with quantiles
  ## mi for y11 (glm)
  bd1_mi_y11_glm_p1a <- nonprob(outcome = formula_y11_p1a, svydesign = sample_prob_svy, data = sample_bd1,
                                method_outcome = "glm", se = TRUE)
  bd1_mi_y12_glm_p1a <- nonprob(outcome = formula_y12_p1a, svydesign = sample_prob_svy, data = sample_bd1, 
                                method_outcome = "glm", se = TRUE)
  bd1_mi_y21_glm_p1a <- nonprob(outcome = formula_y21_p1a, svydesign = sample_prob_svy, data = sample_bd1, 
                                method_outcome = "glm", family_outcome = "binomial", se = TRUE)
  bd1_mi_y22_glm_p1a <- nonprob(outcome = formula_y22_p1a, svydesign = sample_prob_svy, data = sample_bd1, 
                                method_outcome = "glm", family_outcome = "binomial", se = TRUE)
  
  bd1_mi_y11_glm_p1b <- nonprob(outcome = formula_y11_p1b, svydesign = sample_prob_svy, data = sample_bd1, 
                                method_outcome = "glm", se = TRUE)
  bd1_mi_y12_glm_p1b <- nonprob(outcome = formula_y12_p1b, svydesign = sample_prob_svy, data = sample_bd1, 
                                method_outcome = "glm", se = TRUE)
  bd1_mi_y21_glm_p1b <- nonprob(outcome = formula_y21_p1b, svydesign = sample_prob_svy, data = sample_bd1, 
                                method_outcome = "glm", family_outcome = "binomial", se = TRUE)
  bd1_mi_y22_glm_p1b <- nonprob(outcome = formula_y22_p1b, svydesign = sample_prob_svy, data = sample_bd1, 
                                method_outcome = "glm", family_outcome = "binomial", se = TRUE)
  ## with percentiles
  
  ## mi for y11 (glm)
  bd1_mi_y11_glm_p2a <- nonprob(outcome = formula_y11_p2a, svydesign = sample_prob_svy, data = sample_bd1, 
                                method_outcome = "glm", se = TRUE)
  bd1_mi_y12_glm_p2a <- nonprob(outcome = formula_y12_p2a, svydesign = sample_prob_svy, data = sample_bd1, 
                                method_outcome = "glm", se = TRUE)
  bd1_mi_y21_glm_p2a <- nonprob(outcome = formula_y21_p2a, svydesign = sample_prob_svy, data = sample_bd1, 
                                method_outcome = "glm", family_outcome = "binomial", se = TRUE)
  bd1_mi_y22_glm_p2a <- nonprob(outcome = formula_y22_p2a, svydesign = sample_prob_svy, data = sample_bd1, 
                                method_outcome = "glm", family_outcome = "binomial", se = TRUE)
  
  bd1_mi_y11_glm_p2b <- nonprob(outcome = formula_y11_p2b, svydesign = sample_prob_svy, data = sample_bd1, 
                                method_outcome = "glm", se = TRUE)
  bd1_mi_y12_glm_p2b <- nonprob(outcome = formula_y12_p2b, svydesign = sample_prob_svy, data = sample_bd1, 
                                method_outcome = "glm", se = TRUE)
  bd1_mi_y21_glm_p2b <- nonprob(outcome = formula_y21_p2b, svydesign = sample_prob_svy, data = sample_bd1, 
                                method_outcome = "glm", family_outcome = "binomial", se = TRUE)
  bd1_mi_y22_glm_p2b <- nonprob(outcome = formula_y22_p2b, svydesign = sample_prob_svy, data = sample_bd1, 
                                method_outcome = "glm", family_outcome = "binomial", se = TRUE)
  
  ## doubly robust -----------------------------------------------------------
  
  
  ## standard
  bd1_dr_y11_0 <- nonprob(selection = ~ x1 + x2, outcome = y11 ~ x1 + x2,  svydesign = sample_prob_svy, 
                          data = sample_bd1,
                          se = TRUE,
                          control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd1_dr_y12_0 <- nonprob(selection = ~ x1 + x2, outcome = y12 ~ x1 + x2,   svydesign = sample_prob_svy,
                          data = sample_bd1,
                          se = TRUE,
                          control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd1_dr_y21_0 <- nonprob(selection = ~ x1 + x2, outcome = y21 ~ x1 + x2, svydesign = sample_prob_svy,
                          data = sample_bd1,
                          family_outcome = "binomial",se = TRUE,
                          control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd1_dr_y22_0 <- nonprob(selection = ~ x1 + x2, outcome = y22 ~ x1 + x2, svydesign = sample_prob_svy,
                          data = sample_bd1, 
                          family_outcome = "binomial",se = TRUE,
                          control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  ## with quartiles
  bd1_dr_y11_1a <- nonprob(selection = formula_xs_p1a, outcome = formula_y11_p1a,  svydesign = sample_prob_svy, 
                           data = sample_bd1, 
                           se = TRUE,
                           start_selection = rep(0, ncol(bd1_w_cal_1a$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd1_dr_y12_1a <- nonprob(selection = formula_xs_p1a, outcome = formula_y12_p1a, svydesign = sample_prob_svy,
                           data = sample_bd1, 
                           se = TRUE,
                           start_selection = rep(0, ncol(bd1_w_cal_1a$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd1_dr_y21_1a <- nonprob(selection = formula_xs_p1a, outcome = formula_y21_p1a, svydesign = sample_prob_svy,
                           data = sample_bd1, se = TRUE,
                           family_outcome = "binomial",
                           start_selection = rep(0, ncol(bd1_w_cal_1a$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd1_dr_y22_1a <- nonprob(selection = formula_xs_p1a, outcome = formula_y22_p1a, svydesign = sample_prob_svy,
                           data = sample_bd1, se = TRUE,
                           family_outcome = "binomial",
                           start_selection = rep(0, ncol(bd1_w_cal_1a$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd1_dr_y11_1b <- nonprob(selection = formula_xs_p1b, outcome = formula_y11_p1b,  svydesign = sample_prob_svy, 
                           se = TRUE,
                           data = sample_bd1,
                           start_selection = rep(0, ncol(bd1_w_cal_1b$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  
  bd1_dr_y12_1b <- nonprob(selection = formula_xs_p1b, outcome = formula_y12_p1b,   svydesign = sample_prob_svy,
                           se = TRUE,
                           data = sample_bd1,
                           start_selection = rep(0, ncol(bd1_w_cal_1b$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd1_dr_y21_1b <- nonprob(selection = formula_xs_p1b, outcome = formula_y21_p1b, svydesign = sample_prob_svy,
                           data = sample_bd1, se = TRUE,
                           family_outcome = "binomial",
                           start_selection = rep(0, ncol(bd1_w_cal_1b$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd1_dr_y22_1b <- nonprob(selection = formula_xs_p1b, outcome = formula_y22_p1b, svydesign = sample_prob_svy,
                           data = sample_bd1, se = TRUE,
                           family_outcome = "binomial",
                           start_selection = rep(0, ncol(bd1_w_cal_1b$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  
  ## with deciles
  bd1_dr_y11_2a <- nonprob(selection = formula_xs_p2a, outcome = formula_y11_p2a,  svydesign = sample_prob_svy, 
                           data = sample_bd1, se = TRUE,
                           start_selection = rep(0, ncol(bd1_w_cal_2a$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd1_dr_y12_2a <- nonprob(selection = formula_xs_p2a, outcome = formula_y12_p2a,   svydesign = sample_prob_svy,
                           se = TRUE, data = sample_bd1,
                           start_selection = rep(0, ncol(bd1_w_cal_2a$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd1_dr_y21_2a <- nonprob(selection = formula_xs_p2a, outcome = formula_y21_p2a, svydesign = sample_prob_svy,
                           data = sample_bd1, se = TRUE, family_outcome = "binomial",
                           start_selection = rep(0, ncol(bd1_w_cal_2a$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd1_dr_y22_2a <- nonprob(selection = formula_xs_p2a, outcome = formula_y22_p2a, svydesign = sample_prob_svy,
                           data = sample_bd1, se = TRUE, family_outcome = "binomial",
                           start_selection = rep(0, ncol(bd1_w_cal_2a$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd1_dr_y11_2b <- nonprob(selection = formula_xs_p2b, outcome = formula_y11_p2b,  svydesign = sample_prob_svy, 
                           se = TRUE, data = sample_bd1,
                           start_selection = rep(0, ncol(bd1_w_cal_2b$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd1_dr_y12_2b <- nonprob(selection = formula_xs_p2b, outcome = formula_y12_p2b,   svydesign = sample_prob_svy,
                           se = TRUE, data = sample_bd1,
                           start_selection = rep(0, ncol(bd1_w_cal_2b$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd1_dr_y21_2b <- nonprob(selection = formula_xs_p2b, outcome = formula_y21_p2b, svydesign = sample_prob_svy,
                           data = sample_bd1, se = TRUE, family_outcome = "binomial",
                           start_selection = rep(0, ncol(bd1_w_cal_2b$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd1_dr_y22_2b <- nonprob(selection = formula_xs_p2b, outcome = formula_y22_p2b, svydesign = sample_prob_svy,
                           data = sample_bd1, se = TRUE, family_outcome = "binomial",
                           start_selection = rep(0, ncol(bd1_w_cal_2b$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  
  ############################################
  # non-linear inclusion ----------------------------------------------------
  ############################################
  
  ## calibration  --------------------------------------------------------
  ## calibration with totals only (CAL)
  bd2_w_cal_0 <- calib(Xs = with(sample_bd2, cbind(1, x1,x2)), 
                       d = sample_bd2$w_naive, total = c(N,x_totals), method = "raking")
  bd2_w_el_0 <- calib_el(X = with(sample_bd2, cbind(1, x1,x2)), 
                         d = sample_bd2$w_naive, totals = c(N, x_totals))
  
  
  ## calibration with quantiles only (QCAL1A)
  bd2_w_cal_1a <- joint_calib(formula_quantiles = ~ x1 + x2, data = sample_bd2, dweights = sample_bd2$w_naive,
                              N = N,  pop_quantiles = q1_est, method = "raking", backend = "sampling")
  
  ## calibration with quantiles and totals (QCAL1B)
  bd2_w_cal_1b <- joint_calib(formula_quantiles = ~ x1 + x2, formula_totals = ~ x1 + x2, data = sample_bd2,
                              dweights = sample_bd2$w_naive,  N = N,  pop_totals = x_totals,
                              pop_quantiles = q1_est, method = "raking", backend = "sampling")
  
  ## calibration with quantiles only (QCAL2A)
  bd2_w_cal_2a <- joint_calib(formula_quantiles = ~ x1 + x2, data = sample_bd2, dweights = sample_bd2$w_naive,
                              N = N,  pop_quantiles = q2_est, method = "raking", backend = "sampling")
  
  ## calibration with quantiles and totals (QCAL2B)
  bd2_w_cal_2b <- joint_calib(formula_quantiles = ~ x1 + x2, formula_totals = ~ x1 + x2, data = sample_bd2,
                              dweights = sample_bd2$w_naive,  N = N,  pop_totals = x_totals,
                              pop_quantiles = q2_est, method = "raking", backend = "sampling")
  
  ## add weights
  sample_bd2[ , ":="(w_cal_0  = w_naive*bd2_w_cal_0,    w_el_0   = w_naive*bd2_w_el_0,
                     w_cal_1a = w_naive*bd2_w_cal_1a$g, w_cal_1b = w_naive*bd2_w_cal_1b$g, 
                     w_cal_2a = w_naive*bd2_w_cal_2a$g, w_cal_2b = w_naive*bd2_w_cal_2b$g)]
  
  
  ## add XQ variables to sample_bd2
  bd2_p1x <- bd2_w_cal_1a$Xs[, -1]
  colnames(bd2_p1x) <- colnames(sample_p1)
  bd2_p2x <- bd2_w_cal_2a$Xs[, -1]
  colnames(bd2_p2x) <- colnames(sample_p2)
  sample_bd2 <- cbind(sample_bd2, bd2_p1x, bd2_p2x[, setdiff(colnames(bd2_p2x), colnames(bd2_p1x))])
  
  ## IPW  --------------------------------------------------------
  
  bd2_ipw_0 <- nonprob(selection = ~ x1 + x2, target = ~ y11 + y12 + y21 + y22, svydesign = sample_prob_svy,
                       data = sample_bd2,  
                       se = TRUE,
                       control_selection = controlSel(est_method_sel = "gee", h = 1, start_type = "naive"))
 
  
  bd2_ipw_1a <- nonprob(selection = formula_xs_p1a,
                        target = ~ y11 + y12 + y21 + y22,
                        svydesign = sample_prob_svy,
                        data = sample_bd2,
                        start_selection = rep(0, ncol(bd2_w_cal_1a$Xs)),
                        se = TRUE,
                        control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd2_ipw_1b <- nonprob(selection = formula_xs_p1b,
                        target = ~ y11 + y12 + y21 + y22, svydesign = sample_prob_svy,
                        data = sample_bd2,
                        start_selection = rep(0, ncol(bd2_w_cal_1b$Xs)),
                        se = TRUE,
                        control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd2_ipw_2a <- nonprob(selection = formula_xs_p2a,
                        target = ~ y11 + y12 + y21 + y22, svydesign = sample_prob_svy,
                        data = sample_bd2,
                        start_selection = rep(0, ncol(bd2_w_cal_2a$Xs)),
                        se = TRUE,
                        control_selection = controlSel(est_method_sel = "gee", h = 1))

  bd2_ipw_2b <- nonprob(selection = formula_xs_p2b,
                        target = ~ y11 + y12 + y21 + y22, svydesign = sample_prob_svy,
                        data = sample_bd2,
                        start_selection = rep(0, ncol(bd2_w_cal_2b$Xs)),
                        se = TRUE,
                        control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  ## mass imputation ---------------------------------------------------------
  ## mi for y11 (nn)
  bd2_mi_y11 <- nonprob(outcome = y11 ~ x1 + x2, svydesign = sample_prob_svy, data = sample_bd2, 
                        method_outcome = "nn", se = TRUE)
  bd2_mi_y12 <- nonprob(outcome = y12 ~ x1 + x2, svydesign = sample_prob_svy, data = sample_bd2, 
                        method_outcome = "nn", se = TRUE)
  bd2_mi_y21 <- nonprob(outcome = y21 ~ x1 + x2, svydesign = sample_prob_svy, data = sample_bd2, 
                        method_outcome = "nn", se = TRUE)
  bd2_mi_y22 <- nonprob(outcome = y22 ~ x1 + x2, svydesign = sample_prob_svy, data = sample_bd2, 
                        method_outcome = "nn", se = TRUE)
  
  ## mi for y11 (glm)
  bd2_mi_y11_glm <- nonprob(outcome = y11 ~ x1 + x2, svydesign = sample_prob_svy, data = sample_bd2, 
                            method_outcome = "glm", se = TRUE)
  bd2_mi_y12_glm <- nonprob(outcome = y12 ~ x1 + x2, svydesign = sample_prob_svy, data = sample_bd2, 
                            method_outcome = "glm", se = TRUE)
  bd2_mi_y21_glm <- nonprob(outcome = y21 ~ x1 + x2, svydesign = sample_prob_svy, data = sample_bd2, 
                            method_outcome = "glm", family_outcome = "binomial", se = TRUE)
  bd2_mi_y22_glm <- nonprob(outcome = y22 ~ x1 + x2, svydesign = sample_prob_svy, data = sample_bd2, 
                            method_outcome = "glm", family_outcome = "binomial", se = TRUE)
  
  ## with quantiles
  bd2_mi_y11_glm_p1a <- nonprob(outcome = formula_y11_p1a, svydesign = sample_prob_svy, data = sample_bd2,
                                method_outcome = "glm", se = TRUE)
  bd2_mi_y12_glm_p1a <- nonprob(outcome = formula_y12_p1a, svydesign = sample_prob_svy, data = sample_bd2, 
                                method_outcome = "glm", se = TRUE)
  bd2_mi_y21_glm_p1a <- nonprob(outcome = formula_y21_p1a, svydesign = sample_prob_svy, data = sample_bd2, 
                                method_outcome = "glm", family_outcome = "binomial", se = TRUE)
  bd2_mi_y22_glm_p1a <- nonprob(outcome = formula_y22_p1a, svydesign = sample_prob_svy, data = sample_bd2, 
                                method_outcome = "glm", family_outcome = "binomial", se = TRUE)
  
  bd2_mi_y11_glm_p1b <- nonprob(outcome = formula_y11_p1b, svydesign = sample_prob_svy, data = sample_bd2, 
                                method_outcome = "glm", se = TRUE)
  bd2_mi_y12_glm_p1b <- nonprob(outcome = formula_y12_p1b, svydesign = sample_prob_svy, data = sample_bd2, 
                                method_outcome = "glm", se = TRUE)
  bd2_mi_y21_glm_p1b <- nonprob(outcome = formula_y21_p1b, svydesign = sample_prob_svy, data = sample_bd2, 
                                method_outcome = "glm", family_outcome = "binomial", se = TRUE)
  bd2_mi_y22_glm_p1b <- nonprob(outcome = formula_y22_p1b, svydesign = sample_prob_svy, data = sample_bd2, 
                                method_outcome = "glm", family_outcome = "binomial", se = TRUE)
  ## with percentiles
  bd2_mi_y11_glm_p2a <- nonprob(outcome = formula_y11_p2a, svydesign = sample_prob_svy, data = sample_bd2, 
                                method_outcome = "glm", se = TRUE)
  bd2_mi_y12_glm_p2a <- nonprob(outcome = formula_y12_p2a, svydesign = sample_prob_svy, data = sample_bd2, 
                                method_outcome = "glm", se = TRUE)
  bd2_mi_y21_glm_p2a <- nonprob(outcome = formula_y21_p2a, svydesign = sample_prob_svy, data = sample_bd2, 
                                method_outcome = "glm", family_outcome = "binomial", se = TRUE)
  bd2_mi_y22_glm_p2a <- nonprob(outcome = formula_y22_p2a, svydesign = sample_prob_svy, data = sample_bd2, 
                                method_outcome = "glm", family_outcome = "binomial", se = TRUE)
  
  bd2_mi_y11_glm_p2b <- nonprob(outcome = formula_y11_p2b, svydesign = sample_prob_svy, data = sample_bd2, 
                                method_outcome = "glm", se = TRUE)
  bd2_mi_y12_glm_p2b <- nonprob(outcome = formula_y12_p2b, svydesign = sample_prob_svy, data = sample_bd2, 
                                method_outcome = "glm", se = TRUE)
  bd2_mi_y21_glm_p2b <- nonprob(outcome = formula_y21_p2b, svydesign = sample_prob_svy, data = sample_bd2, 
                                method_outcome = "glm", family_outcome = "binomial", se = TRUE)
  bd2_mi_y22_glm_p2b <- nonprob(outcome = formula_y22_p2b, svydesign = sample_prob_svy, data = sample_bd2, 
                                method_outcome = "glm", family_outcome = "binomial", se = TRUE)
  
  ## doubly robust -----------------------------------------------------------
  
  ## standard
  bd2_dr_y11_0 <- nonprob(selection = ~ x1 + x2, outcome = y11 ~ x1 + x2,  svydesign = sample_prob_svy, 
                          data = sample_bd2,
                          se = TRUE,
                          control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd2_dr_y12_0 <- nonprob(selection = ~ x1 + x2, outcome = y12 ~ x1 + x2,   
                          svydesign = sample_prob_svy,
                          data = sample_bd2,
                          se = TRUE,
                          control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd2_dr_y21_0 <- nonprob(selection = ~ x1 + x2, outcome = y21 ~ x1 + x2, svydesign = sample_prob_svy,
                          data = sample_bd2,
                          family_outcome = "binomial",se = TRUE,
                          control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd2_dr_y22_0 <- nonprob(selection = ~ x1 + x2, outcome = y22 ~ x1 + x2, svydesign = sample_prob_svy,
                          data = sample_bd2, 
                          family_outcome = "binomial",se = TRUE,
                          control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  ## with quartiles
  bd2_dr_y11_1a <- nonprob(selection = formula_xs_p1a, 
                           outcome = formula_y11_p1a,  
                           svydesign = sample_prob_svy, 
                           data = sample_bd2,  se = TRUE,
                           start_selection = rep(0, ncol(bd2_w_cal_1a$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd2_dr_y12_1a <- nonprob(selection = formula_xs_p1a, 
                           outcome = formula_y12_p1a, svydesign = sample_prob_svy,
                           data = sample_bd2,  se = TRUE,
                           start_selection = rep(0, ncol(bd2_w_cal_1a$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd2_dr_y21_1a <- nonprob(selection = formula_xs_p1a, outcome = formula_y21_p1a, svydesign = sample_prob_svy,
                           data = sample_bd2, se = TRUE, family_outcome = "binomial",
                           start_selection = rep(0, ncol(bd2_w_cal_1a$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd2_dr_y22_1a <- nonprob(selection = formula_xs_p1a, outcome = formula_y22_p1a, svydesign = sample_prob_svy,
                           data = sample_bd2, se = TRUE, family_outcome = "binomial",
                           start_selection = rep(0, ncol(bd2_w_cal_1a$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd2_dr_y11_1b <- nonprob(selection = formula_xs_p1b, outcome = formula_y11_p1b,  svydesign = sample_prob_svy, 
                           se = TRUE, data = sample_bd2,
                           start_selection = rep(0, ncol(bd2_w_cal_1b$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd2_dr_y12_1b <- nonprob(selection = formula_xs_p1b, outcome = formula_y12_p1b,   svydesign = sample_prob_svy,
                           se = TRUE, data = sample_bd2,
                           start_selection = rep(0, ncol(bd2_w_cal_1b$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd2_dr_y21_1b <- nonprob(selection = formula_xs_p1b, outcome = formula_y21_p1b, svydesign = sample_prob_svy,
                           data = sample_bd2, se = TRUE, family_outcome = "binomial",
                           start_selection = rep(0, ncol(bd2_w_cal_1b$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd2_dr_y22_1b <- nonprob(selection = formula_xs_p1b, outcome = formula_y22_p1b, svydesign = sample_prob_svy,
                           data = sample_bd2, se = TRUE, family_outcome = "binomial",
                           start_selection = rep(0, ncol(bd2_w_cal_1b$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  
  ## with deciles
  bd2_dr_y11_2a <- nonprob(selection = formula_xs_p2a, 
                           outcome = formula_y11_p2a,  
                           svydesign = sample_prob_svy, 
                           data = sample_bd2, se = TRUE,
                           start_selection = rep(0, ncol(bd2_w_cal_2a$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd2_dr_y12_2a <- nonprob(selection = formula_xs_p2a, outcome = formula_y12_p2a,   svydesign = sample_prob_svy,
                           se = TRUE, data = sample_bd2,
                           start_selection = rep(0, ncol(bd2_w_cal_2a$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd2_dr_y21_2a <- nonprob(selection = formula_xs_p2a, outcome = formula_y21_p2a, svydesign = sample_prob_svy,
                           data = sample_bd2, se = TRUE, family_outcome = "binomial",
                           start_selection = rep(0, ncol(bd2_w_cal_2a$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd2_dr_y22_2a <- nonprob(selection = formula_xs_p2a, outcome = formula_y22_p2a, svydesign = sample_prob_svy,
                           data = sample_bd2, se = TRUE, family_outcome = "binomial",
                           start_selection = rep(0, ncol(bd2_w_cal_2a$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd2_dr_y11_2b <- nonprob(selection = formula_xs_p2b, outcome = formula_y11_p2b,  svydesign = sample_prob_svy, 
                           se = TRUE, data = sample_bd2,
                           start_selection = rep(0, ncol(bd2_w_cal_2b$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd2_dr_y12_2b <- nonprob(selection = formula_xs_p2b, outcome = formula_y12_p2b,   svydesign = sample_prob_svy,
                           se = TRUE, data = sample_bd2,
                           start_selection = rep(0, ncol(bd2_w_cal_2b$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd2_dr_y21_2b <- nonprob(selection = formula_xs_p2b, outcome = formula_y21_p2b, svydesign = sample_prob_svy,
                           data = sample_bd2, se = TRUE, family_outcome = "binomial",
                           start_selection = rep(0, ncol(bd2_w_cal_2b$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd2_dr_y22_2b <- nonprob(selection = formula_xs_p2b, outcome = formula_y22_p2b, svydesign = sample_prob_svy,
                           data = sample_bd2, se = TRUE, family_outcome = "binomial",
                           start_selection = rep(0, ncol(bd2_w_cal_2b$Xs)),
                           control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  
  # save results ------------------------------------------------------------
  
  
  ## save results for ipw, mi, dr
  
  df1 <- data.frame(est = c(
    rep("naive", 4),
    rep("cal0", 4), rep("cal1a", 4), rep("cal1b", 4), rep("cal2a", 4), rep("cal2b", 4), 
    rep("ipw0", 4), 
    rep("ipw1a", 4), rep("ipw1b", 4), rep("ipw2a", 4), rep("ipw2b", 4),
    rep("mi_nn", 4), 
    rep("mi_glm0", 4),
    rep("mi_glm1a", 4), rep("mi_glm1b", 4), rep("mi_glm2a", 4), rep("mi_glm2b", 4),
    rep("dr0", 4), 
    rep("dr1a", 4), rep("dr1b", 4),  rep("dr2a", 4), rep("dr2b", 4)),
    y = rep(c("y11", "y12", "y21", "y22"), times = 22),
    rbind(
      data.frame(mean = unlist(sample_bd1[, lapply(.SD, mean), .SDcols = patterns("y")]), 
                 SE=NA, lower_bound=NA, upper_bound = NA),
      data.frame(mean = unlist(sample_bd1[, lapply(.SD, weighted.mean, w = w_cal_0), .SDcols = patterns("y")]), 
                 SE=NA, lower_bound=NA, upper_bound = NA),
      data.frame(mean = unlist(sample_bd1[, lapply(.SD, weighted.mean, w = w_cal_1a), .SDcols = patterns("y")]), 
                 SE=NA, lower_bound=NA, upper_bound = NA),
      data.frame(mean = unlist(sample_bd1[, lapply(.SD, weighted.mean, w = w_cal_1b), .SDcols = patterns("y")]), 
                 SE=NA, lower_bound=NA, upper_bound = NA),
      data.frame(mean = unlist(sample_bd1[, lapply(.SD, weighted.mean, w = w_cal_2a), .SDcols = patterns("y")]), 
                 SE=NA, lower_bound=NA, upper_bound = NA),
      data.frame(mean = unlist(sample_bd1[, lapply(.SD, weighted.mean, w = w_cal_2b), .SDcols = patterns("y")]), 
                 SE=NA, lower_bound=NA, upper_bound = NA),                     
      cbind(bd1_ipw_0$output,      bd1_ipw_0$confidence_interval),
      cbind(bd1_ipw_1a$output,      bd1_ipw_1a$confidence_interval),
      cbind(bd1_ipw_1b$output,      bd1_ipw_1b$confidence_interval),
      cbind(bd1_ipw_2a$output,      bd1_ipw_2a$confidence_interval),
      cbind(bd1_ipw_2b$output,      bd1_ipw_2b$confidence_interval),
      cbind(bd1_mi_y11$output,     bd1_mi_y11$confidence_interval),
      cbind(bd1_mi_y12$output,     bd1_mi_y12$confidence_interval),
      cbind(bd1_mi_y21$output,     bd1_mi_y21$confidence_interval),
      cbind(bd1_mi_y22$output,     bd1_mi_y22$confidence_interval),
      cbind(bd1_mi_y11_glm$output, bd1_mi_y11_glm$confidence_interval),
      cbind(bd1_mi_y12_glm$output, bd1_mi_y12_glm$confidence_interval),
      cbind(bd1_mi_y21_glm$output, bd1_mi_y21_glm$confidence_interval),
      cbind(bd1_mi_y22_glm$output, bd1_mi_y22_glm$confidence_interval),
      
      cbind(bd1_mi_y11_glm_p1a$output, bd1_mi_y11_glm_p1a$confidence_interval),
      cbind(bd1_mi_y12_glm_p1a$output, bd1_mi_y12_glm_p1a$confidence_interval),
      cbind(bd1_mi_y21_glm_p1a$output, bd1_mi_y21_glm_p1a$confidence_interval),
      cbind(bd1_mi_y22_glm_p1a$output, bd1_mi_y22_glm_p1a$confidence_interval),
      cbind(bd1_mi_y11_glm_p1b$output, bd1_mi_y11_glm_p1b$confidence_interval),
      cbind(bd1_mi_y12_glm_p1b$output, bd1_mi_y12_glm_p1b$confidence_interval),
      cbind(bd1_mi_y21_glm_p1b$output, bd1_mi_y21_glm_p1b$confidence_interval),
      cbind(bd1_mi_y22_glm_p1b$output, bd1_mi_y22_glm_p1b$confidence_interval),
      
      cbind(bd1_mi_y11_glm_p2a$output, bd1_mi_y11_glm_p2a$confidence_interval),
      cbind(bd1_mi_y12_glm_p2a$output, bd1_mi_y12_glm_p2a$confidence_interval),
      cbind(bd1_mi_y21_glm_p2a$output, bd1_mi_y21_glm_p2a$confidence_interval),
      cbind(bd1_mi_y22_glm_p2a$output, bd1_mi_y22_glm_p2a$confidence_interval),
      cbind(bd1_mi_y11_glm_p2b$output, bd1_mi_y11_glm_p2b$confidence_interval),
      cbind(bd1_mi_y12_glm_p2b$output, bd1_mi_y12_glm_p2b$confidence_interval),
      cbind(bd1_mi_y21_glm_p2b$output, bd1_mi_y21_glm_p2b$confidence_interval),
      cbind(bd1_mi_y22_glm_p2b$output, bd1_mi_y22_glm_p2b$confidence_interval),
      
      cbind(bd1_dr_y11_0$output, bd1_dr_y11_0$confidence_interval),
      cbind(bd1_dr_y12_0$output, bd1_dr_y12_0$confidence_interval),
      cbind(bd1_dr_y21_0$output, bd1_dr_y21_0$confidence_interval),
      cbind(bd1_dr_y22_0$output, bd1_dr_y22_0$confidence_interval),
      
      cbind(bd1_dr_y11_1a$output, bd1_dr_y11_1a$confidence_interval),
      cbind(bd1_dr_y12_1a$output, bd1_dr_y12_1a$confidence_interval),
      cbind(bd1_dr_y21_1a$output, bd1_dr_y21_1a$confidence_interval),
      cbind(bd1_dr_y22_1a$output, bd1_dr_y22_1a$confidence_interval),
      
      cbind(bd1_dr_y11_1b$output, bd1_dr_y11_1b$confidence_interval),
      cbind(bd1_dr_y12_1b$output, bd1_dr_y12_1b$confidence_interval),
      cbind(bd1_dr_y21_1b$output, bd1_dr_y21_1b$confidence_interval),
      cbind(bd1_dr_y22_1b$output, bd1_dr_y22_1b$confidence_interval),
      
      cbind(bd1_dr_y11_2a$output, bd1_dr_y11_2a$confidence_interval),
      cbind(bd1_dr_y12_2a$output, bd1_dr_y12_2a$confidence_interval),
      cbind(bd1_dr_y21_2a$output, bd1_dr_y21_2a$confidence_interval),
      cbind(bd1_dr_y22_2a$output, bd1_dr_y22_2a$confidence_interval),
      
      cbind(bd1_dr_y11_2b$output, bd1_dr_y11_2b$confidence_interval),
      cbind(bd1_dr_y12_2b$output, bd1_dr_y12_2b$confidence_interval),
      cbind(bd1_dr_y21_2b$output, bd1_dr_y21_2b$confidence_interval),
      cbind(bd1_dr_y22_2b$output, bd1_dr_y22_2b$confidence_interval)
    ))
  
  df1$bd <- 1
  
  df2 <- data.frame(est = c(
    rep("naive", 4),
    rep("cal0", 4), rep("cal1a", 4), rep("cal1b", 4), rep("cal2a", 4), rep("cal2b", 4), 
    rep("ipw0", 4), 
    rep("ipw1a", 4), rep("ipw1b", 4), rep("ipw2a", 4), rep("ipw2b", 4),
    rep("mi_nn", 4), 
    rep("mi_glm0", 4),
    rep("mi_glm1a", 4), rep("mi_glm1b", 4), rep("mi_glm2a", 4), rep("mi_glm2b", 4),
    rep("dr0", 4), 
    rep("dr1a", 4), rep("dr1b", 4),  rep("dr2a", 4), rep("dr2b", 4)),
    y = rep(c("y11", "y12", "y21", "y22"), times = 22),
    rbind(
      data.frame(mean = unlist(sample_bd2[, lapply(.SD, mean), .SDcols = patterns("y")]), 
                 SE=NA, lower_bound=NA, upper_bound = NA),
      data.frame(mean = unlist(sample_bd2[, lapply(.SD, weighted.mean, w = w_cal_0), .SDcols = patterns("y")]), 
                 SE=NA, lower_bound=NA, upper_bound = NA),
      data.frame(mean = unlist(sample_bd2[, lapply(.SD, weighted.mean, w = w_cal_1a), .SDcols = patterns("y")]), 
                 SE=NA, lower_bound=NA, upper_bound = NA),
      data.frame(mean = unlist(sample_bd2[, lapply(.SD, weighted.mean, w = w_cal_1b), .SDcols = patterns("y")]), 
                 SE=NA, lower_bound=NA, upper_bound = NA),
      data.frame(mean = unlist(sample_bd2[, lapply(.SD, weighted.mean, w = w_cal_2a), .SDcols = patterns("y")]), 
                 SE=NA, lower_bound=NA, upper_bound = NA),
      data.frame(mean = unlist(sample_bd2[, lapply(.SD, weighted.mean, w = w_cal_2b), .SDcols = patterns("y")]), 
                 SE=NA, lower_bound=NA, upper_bound = NA),                    
      cbind(bd2_ipw_0$output,  bd2_ipw_0$confidence_interval),
      cbind(bd2_ipw_1a$output, bd2_ipw_1a$confidence_interval),
      cbind(bd2_ipw_1b$output, bd2_ipw_1b$confidence_interval),
      cbind(bd2_ipw_2a$output, bd2_ipw_2a$confidence_interval),
      cbind(bd2_ipw_2b$output, bd2_ipw_2b$confidence_interval),
      cbind(bd2_mi_y11$output, bd2_mi_y11$confidence_interval),
      cbind(bd2_mi_y12$output, bd2_mi_y12$confidence_interval),
      cbind(bd2_mi_y21$output, bd2_mi_y21$confidence_interval),
      cbind(bd2_mi_y22$output, bd2_mi_y22$confidence_interval),
      
      cbind(bd2_mi_y11_glm$output, bd2_mi_y11_glm$confidence_interval),
      cbind(bd2_mi_y12_glm$output, bd2_mi_y12_glm$confidence_interval),
      cbind(bd2_mi_y21_glm$output, bd2_mi_y21_glm$confidence_interval),
      cbind(bd2_mi_y22_glm$output, bd2_mi_y22_glm$confidence_interval),
      
      cbind(bd2_mi_y11_glm_p1a$output, bd2_mi_y11_glm_p1a$confidence_interval),
      cbind(bd2_mi_y12_glm_p1a$output, bd2_mi_y12_glm_p1a$confidence_interval),
      cbind(bd2_mi_y21_glm_p1a$output, bd2_mi_y21_glm_p1a$confidence_interval),
      cbind(bd2_mi_y22_glm_p1a$output, bd2_mi_y22_glm_p1a$confidence_interval),
      
      cbind(bd2_mi_y11_glm_p1b$output, bd2_mi_y11_glm_p1b$confidence_interval),
      cbind(bd2_mi_y12_glm_p1b$output, bd2_mi_y12_glm_p1b$confidence_interval),
      cbind(bd2_mi_y21_glm_p1b$output, bd2_mi_y21_glm_p1b$confidence_interval),
      cbind(bd2_mi_y22_glm_p1b$output, bd2_mi_y22_glm_p1b$confidence_interval),
      
      cbind(bd2_mi_y11_glm_p2a$output, bd2_mi_y11_glm_p2a$confidence_interval),
      cbind(bd2_mi_y12_glm_p2a$output, bd2_mi_y12_glm_p2a$confidence_interval),
      cbind(bd2_mi_y21_glm_p2a$output, bd2_mi_y21_glm_p2a$confidence_interval),
      cbind(bd2_mi_y22_glm_p2a$output, bd2_mi_y22_glm_p2a$confidence_interval),
      
      cbind(bd2_mi_y11_glm_p2b$output, bd2_mi_y11_glm_p2b$confidence_interval),
      cbind(bd2_mi_y12_glm_p2b$output, bd2_mi_y12_glm_p2b$confidence_interval),
      cbind(bd2_mi_y21_glm_p2b$output, bd2_mi_y21_glm_p2b$confidence_interval),
      cbind(bd2_mi_y22_glm_p2b$output, bd2_mi_y22_glm_p2b$confidence_interval),
      
      cbind(bd2_dr_y11_0$output, bd2_dr_y11_0$confidence_interval),
      cbind(bd2_dr_y12_0$output, bd2_dr_y12_0$confidence_interval),
      cbind(bd2_dr_y21_0$output, bd2_dr_y21_0$confidence_interval),
      cbind(bd2_dr_y22_0$output, bd2_dr_y22_0$confidence_interval),
      
      cbind(bd2_dr_y11_1a$output, bd2_dr_y11_1a$confidence_interval),
      cbind(bd2_dr_y12_1a$output, bd2_dr_y12_1a$confidence_interval),
      cbind(bd2_dr_y21_1a$output, bd2_dr_y21_1a$confidence_interval),
      cbind(bd2_dr_y22_1a$output, bd2_dr_y22_1a$confidence_interval),
      
      cbind(bd2_dr_y11_1b$output, bd2_dr_y11_1b$confidence_interval),
      cbind(bd2_dr_y12_1b$output, bd2_dr_y12_1b$confidence_interval),
      cbind(bd2_dr_y21_1b$output, bd2_dr_y21_1b$confidence_interval),
      cbind(bd2_dr_y22_1b$output, bd2_dr_y22_1b$confidence_interval),
      
      cbind(bd2_dr_y11_2a$output, bd2_dr_y11_2a$confidence_interval),
      cbind(bd2_dr_y12_2a$output, bd2_dr_y12_2a$confidence_interval),
      cbind(bd2_dr_y21_2a$output, bd2_dr_y21_2a$confidence_interval),
      cbind(bd2_dr_y22_2a$output, bd2_dr_y22_2a$confidence_interval),
      
      cbind(bd2_dr_y11_2b$output, bd2_dr_y11_2b$confidence_interval),
      cbind(bd2_dr_y12_2b$output, bd2_dr_y12_2b$confidence_interval),
      cbind(bd2_dr_y21_2b$output, bd2_dr_y21_2b$confidence_interval),
      cbind(bd2_dr_y22_2b$output, bd2_dr_y22_2b$confidence_interval)
    ))
  df2$bd <- 2
  result <- rbind(df1, df2)
  result$k <- k
  result
}

