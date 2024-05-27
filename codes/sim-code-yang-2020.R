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
  sample_bd2 <- pop_data[rbinom(N,1,pop_data$p2)==1, ]
  
  ##########################################
  # linear inclusion --------------------------------------------------------
  ##########################################
  sample_np1 <- joint_calib_create_matrix(sample_bd1[, c("x1", "x2")], 
                                         sum(sample_prob$w), 
                                         pop_quantiles = list(x1=quantile(sample_prob$x1, p_quantiles1), 
                                                              x2=quantile(sample_prob$x2, p_quantiles1)))
  
  sample_np2 <- joint_calib_create_matrix(sample_bd1[, c("x1", "x2")], 
                                          sum(sample_prob$w),
                                         pop_quantiles = list(x1=quantile(sample_prob$x1, p_quantiles2), 
                                                              x2=quantile(sample_prob$x2, p_quantiles2)))
  
  ## add XQ variables to sample_bd1
  colnames(sample_np1) <- colnames(sample_p1)
  colnames(sample_np2) <- colnames(sample_p2)
  sample_bd1 <- cbind(sample_bd1, sample_np1, sample_np2[, setdiff(colnames(sample_np2), colnames(sample_np1))])
  
  ## IPW  not_calibrated --------------------------------------------------------
  bd1_ipw_standard_0 <- nonprob(selection = ~ x1 + x2, target = ~ y11 + y12 + y21 + y22, 
                                svydesign = sample_prob_svy, data = sample_bd1, se = TRUE)
  
  bd1_ipw_standard_1a <- nonprob(selection = formula_xs_p1a, target = ~ y11 + y12 + y21 + y22, 
                                 svydesign = sample_prob_svy, data = sample_bd1, se = TRUE, 
                                 start_selection = rep(0, length(attr(terms(formula_xs_p1a), "term.labels")) +1))
  
  bd1_ipw_standard_1b <- nonprob(selection = formula_xs_p1b, target = ~ y11 + y12 + y21 + y22, 
                                 svydesign = sample_prob_svy, data = sample_bd1, se = TRUE, 
                                 start_selection = rep(0, length(attr(terms(formula_xs_p1b), "term.labels")) + 1))
  
  bd1_ipw_standard_2a <- nonprob(selection = formula_xs_p2a, target = ~ y11 + y12 + y21 + y22, 
                                 svydesign = sample_prob_svy, data = sample_bd1, se = TRUE, 
                                 start_selection = rep(0, length(attr(terms(formula_xs_p2a), "term.labels")) + 1))
  
  bd1_ipw_standard_2b <- nonprob(selection = formula_xs_p2b, target = ~ y11 + y12 + y21 + y22, 
                                 svydesign = sample_prob_svy, data = sample_bd1, se = TRUE, 
                                 start_selection = rep(0, length(attr(terms(formula_xs_p2b), "term.labels")) + 1))
  
  ## IPW  calibrated --------------------------------------------------------
  
  bd1_ipw_0 <- nonprob(selection = ~ x1 + x2, target = ~ y11 + y12 + y21 + y22, 
                       svydesign = sample_prob_svy, data = sample_bd1, se = TRUE, 
                       control_selection = controlSel(est_method_sel = "gee", h = 1))

  bd1_ipw_1a <- nonprob(selection = formula_xs_p1a, target = ~ y11 + y12 + y21 + y22, 
                        svydesign = sample_prob_svy, data = sample_bd1, se = TRUE, 
                        start_selection = rep(0, length(attr(terms(formula_xs_p1a), "term.labels")) +1),
                        control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd1_ipw_1b <- nonprob(selection = formula_xs_p1b, target = ~ y11 + y12 + y21 + y22, 
                        svydesign = sample_prob_svy, data = sample_bd1, se = TRUE, 
                        start_selection = rep(0, length(attr(terms(formula_xs_p1b), "term.labels")) +1),
                        control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd1_ipw_2a <- nonprob(selection = formula_xs_p2a, target = ~ y11 + y12 + y21 + y22, 
                        svydesign = sample_prob_svy, data = sample_bd1, se = TRUE, 
                        start_selection = rep(0, length(attr(terms(formula_xs_p2a), "term.labels")) +1),
                        control_selection = controlSel(est_method_sel = "gee", h = 1))

  bd1_ipw_2b <- nonprob(selection = formula_xs_p2b, target = ~ y11 + y12 + y21 + y22, 
                        svydesign = sample_prob_svy, data = sample_bd1, se = TRUE, 
                        start_selection = rep(0, length(attr(terms(formula_xs_p2b), "term.labels")) +1),
                        control_selection = controlSel(est_method_sel = "gee", h = 1))

  
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
  
  ## doubly robust with IPW not calibrated -----------------------------------------------------------
  
  bd1_dr_y11_standard_0 <- nonprob(selection = ~ x1 + x2, outcome = y11 ~ x1 + x2,  svydesign = sample_prob_svy, 
                          data = sample_bd1,
                          se = TRUE)
  
  bd1_dr_y12_standard_0 <- nonprob(selection = ~ x1 + x2, outcome = y12 ~ x1 + x2,   svydesign = sample_prob_svy,
                          data = sample_bd1,
                          se = TRUE)
  
  bd1_dr_y21_standard_0 <- nonprob(selection = ~ x1 + x2, outcome = y21 ~ x1 + x2, svydesign = sample_prob_svy,
                          data = sample_bd1,
                          family_outcome = "binomial",se = TRUE)
  
  bd1_dr_y22_standard_0 <- nonprob(selection = ~ x1 + x2, outcome = y22 ~ x1 + x2, svydesign = sample_prob_svy,
                          data = sample_bd1, 
                          family_outcome = "binomial",se = TRUE)
  
  ## doubly robust with IPW with calibration -----------------------------------

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
  
  ############################################
  # non-linear inclusion ----------------------------------------------------
  ############################################
  
  sample_np1 <- joint_calib_create_matrix(sample_bd2[, c("x1", "x2")], 
                                          sum(sample_prob$w), 
                                          pop_quantiles = list(x1=quantile(sample_prob$x1, p_quantiles1), 
                                                               x2=quantile(sample_prob$x2, p_quantiles1)))
  
  sample_np2 <- joint_calib_create_matrix(sample_bd2[, c("x1", "x2")], 
                                          sum(sample_prob$w),
                                          pop_quantiles = list(x1=quantile(sample_prob$x1, p_quantiles2), 
                                                               x2=quantile(sample_prob$x2, p_quantiles2)))
  
  ## add XQ variables to sample_bd1
  colnames(sample_np1) <- colnames(sample_p1)
  colnames(sample_np2) <- colnames(sample_p2)
  sample_bd2 <- cbind(sample_bd2, sample_np1, sample_np2[, setdiff(colnames(sample_np2), colnames(sample_np1))])
  
  ## IPW  not_calibrated --------------------------------------------------------
  bd2_ipw_standard_0 <- nonprob(selection = ~ x1 + x2, target = ~ y11 + y12 + y21 + y22, 
                                svydesign = sample_prob_svy, data = sample_bd2, se = TRUE)
  
  bd2_ipw_standard_1a <- nonprob(selection = formula_xs_p1a, target = ~ y11 + y12 + y21 + y22, 
                                 svydesign = sample_prob_svy, data = sample_bd2, se = TRUE, 
                                 start_selection = rep(0, length(attr(terms(formula_xs_p1a), "term.labels")) +1))
  
  bd2_ipw_standard_1b <- nonprob(selection = formula_xs_p1b, target = ~ y11 + y12 + y21 + y22, 
                                 svydesign = sample_prob_svy, data = sample_bd2, se = TRUE, 
                                 start_selection = rep(0, length(attr(terms(formula_xs_p1b), "term.labels")) + 1))
  
  bd2_ipw_standard_2a <- nonprob(selection = formula_xs_p2a, target = ~ y11 + y12 + y21 + y22, 
                                 svydesign = sample_prob_svy, data = sample_bd2, se = TRUE, 
                                 start_selection = rep(0, length(attr(terms(formula_xs_p2a), "term.labels")) + 1))
  
  bd2_ipw_standard_2b <- nonprob(selection = formula_xs_p2b, target = ~ y11 + y12 + y21 + y22, 
                                 svydesign = sample_prob_svy, data = sample_bd2, se = TRUE, 
                                 start_selection = rep(0, length(attr(terms(formula_xs_p2b), "term.labels")) + 1))
  
  ## IPW  calibrated --------------------------------------------------------
  
  bd2_ipw_0 <- nonprob(selection = ~ x1 + x2, target = ~ y11 + y12 + y21 + y22, 
                       svydesign = sample_prob_svy, data = sample_bd2, se = TRUE, 
                       control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd2_ipw_1a <- nonprob(selection = formula_xs_p1a, target = ~ y11 + y12 + y21 + y22, 
                        svydesign = sample_prob_svy, data = sample_bd2, se = TRUE, 
                        start_selection = rep(0, length(attr(terms(formula_xs_p1a), "term.labels")) +1),
                        control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd2_ipw_1b <- nonprob(selection = formula_xs_p1b, target = ~ y11 + y12 + y21 + y22, 
                        svydesign = sample_prob_svy, data = sample_bd2, se = TRUE, 
                        start_selection = rep(0, length(attr(terms(formula_xs_p1b), "term.labels")) +1),
                        control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd2_ipw_2a <- nonprob(selection = formula_xs_p2a, target = ~ y11 + y12 + y21 + y22, 
                        svydesign = sample_prob_svy, data = sample_bd2, se = TRUE, 
                        start_selection = rep(0, length(attr(terms(formula_xs_p2a), "term.labels")) +1),
                        control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd2_ipw_2b <- nonprob(selection = formula_xs_p2b, target = ~ y11 + y12 + y21 + y22, 
                        svydesign = sample_prob_svy, data = sample_bd2, se = TRUE, 
                        start_selection = rep(0, length(attr(terms(formula_xs_p2b), "term.labels")) +1),
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
  
  ## doubly robust with IPW not calibrated -----------------------------------------------------------
  
  bd2_dr_y11_standard_0 <- nonprob(selection = ~ x1 + x2, outcome = y11 ~ x1 + x2,  svydesign = sample_prob_svy, 
                                   data = sample_bd2, se = TRUE)
  
  bd2_dr_y12_standard_0 <- nonprob(selection = ~ x1 + x2, outcome = y12 ~ x1 + x2,   svydesign = sample_prob_svy,
                                   data = sample_bd2, se = TRUE)
  
  bd2_dr_y21_standard_0 <- nonprob(selection = ~ x1 + x2, outcome = y21 ~ x1 + x2, svydesign = sample_prob_svy,
                                   data = sample_bd2,
                                   family_outcome = "binomial",se = TRUE)
  
  bd2_dr_y22_standard_0 <- nonprob(selection = ~ x1 + x2, outcome = y22 ~ x1 + x2, svydesign = sample_prob_svy,
                                   data = sample_bd2, 
                                   family_outcome = "binomial",se = TRUE)
  
  ## doubly robust with IPW with calibration -----------------------------------
  
  bd2_dr_y11_0 <- nonprob(selection = ~ x1 + x2, outcome = y11 ~ x1 + x2,  svydesign = sample_prob_svy, 
                          data = sample_bd2,
                          se = TRUE,
                          control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd2_dr_y12_0 <- nonprob(selection = ~ x1 + x2, outcome = y12 ~ x1 + x2,   svydesign = sample_prob_svy,
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
  
  # save results ------------------------------------------------------------
  
  
  ## save results for ipw, mi, dr
  
  df1 <- data.frame(est = c(
    rep("naive", 4),
    
    rep("ipwst0", 4), 
    rep("ipwst1a", 4), rep("ipwst1b", 4), rep("ipwst2a", 4), rep("ipwst2b", 4),
    rep("ipw0", 4), 
    rep("ipw1a", 4), rep("ipw1b", 4), rep("ipw2a", 4), rep("ipw2b", 4),
    rep("mi_nn", 4), 
    rep("mi_glm0", 4),
    rep("drst0", 4),
    rep("dr0", 4)),
    y = rep(c("y11", "y12", "y21", "y22"), times = 15),
    rbind(
      data.frame(mean = unlist(sample_bd1[, lapply(.SD, mean), .SDcols = patterns("y")]), 
                 SE=NA, lower_bound=NA, upper_bound = NA),
      
      cbind(bd1_ipw_standard_0$output,       bd1_ipw_standard_0$confidence_interval),
      cbind(bd1_ipw_standard_1a$output,      bd1_ipw_standard_1a$confidence_interval),
      cbind(bd1_ipw_standard_1b$output,      bd1_ipw_standard_1b$confidence_interval),
      cbind(bd1_ipw_standard_2a$output,      bd1_ipw_standard_2a$confidence_interval),
      cbind(bd1_ipw_standard_2b$output,      bd1_ipw_standard_2b$confidence_interval),
      
      cbind(bd1_ipw_0$output,       bd1_ipw_0$confidence_interval),
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
      
      cbind(bd1_dr_y11_standard_0$output, bd1_dr_y11_standard_0$confidence_interval),
      cbind(bd1_dr_y12_standard_0$output, bd1_dr_y12_standard_0$confidence_interval),
      cbind(bd1_dr_y21_standard_0$output, bd1_dr_y21_standard_0$confidence_interval),
      cbind(bd1_dr_y22_standard_0$output, bd1_dr_y22_standard_0$confidence_interval),
      
      cbind(bd1_dr_y11_0$output, bd1_dr_y11_0$confidence_interval),
      cbind(bd1_dr_y12_0$output, bd1_dr_y12_0$confidence_interval),
      cbind(bd1_dr_y21_0$output, bd1_dr_y21_0$confidence_interval),
      cbind(bd1_dr_y22_0$output, bd1_dr_y22_0$confidence_interval)
      ))
  
  df1$bd <- 1
  
  df2 <- data.frame(est = c(
    rep("naive", 4),
    
    rep("ipwst0", 4), 
    rep("ipwst1a", 4), rep("ipwst1b", 4), rep("ipwst2a", 4), rep("ipwst2b", 4),
    rep("ipw0", 4), 
    rep("ipw1a", 4), rep("ipw1b", 4), rep("ipw2a", 4), rep("ipw2b", 4),
    rep("mi_nn", 4), 
    rep("mi_glm0", 4),
    rep("drst0", 4),
    rep("dr0", 4)),
    y = rep(c("y11", "y12", "y21", "y22"), times = 15),
    rbind(
      data.frame(mean = unlist(sample_bd1[, lapply(.SD, mean), .SDcols = patterns("y")]), 
                 SE=NA, lower_bound=NA, upper_bound = NA),
      
      cbind(bd2_ipw_standard_0$output,       bd2_ipw_standard_0$confidence_interval),
      cbind(bd2_ipw_standard_1a$output,      bd2_ipw_standard_1a$confidence_interval),
      cbind(bd2_ipw_standard_1b$output,      bd2_ipw_standard_1b$confidence_interval),
      cbind(bd2_ipw_standard_2a$output,      bd2_ipw_standard_2a$confidence_interval),
      cbind(bd2_ipw_standard_2b$output,      bd2_ipw_standard_2b$confidence_interval),
      
      cbind(bd2_ipw_0$output,       bd2_ipw_0$confidence_interval),
      cbind(bd2_ipw_1a$output,      bd2_ipw_1a$confidence_interval),
      cbind(bd2_ipw_1b$output,      bd2_ipw_1b$confidence_interval),
      cbind(bd2_ipw_2a$output,      bd2_ipw_2a$confidence_interval),
      cbind(bd2_ipw_2b$output,      bd2_ipw_2b$confidence_interval),
      
      cbind(bd2_mi_y11$output,     bd2_mi_y11$confidence_interval),
      cbind(bd2_mi_y12$output,     bd2_mi_y12$confidence_interval),
      cbind(bd2_mi_y21$output,     bd2_mi_y21$confidence_interval),
      cbind(bd2_mi_y22$output,     bd2_mi_y22$confidence_interval),
      
      cbind(bd2_mi_y11_glm$output, bd2_mi_y11_glm$confidence_interval),
      cbind(bd2_mi_y12_glm$output, bd2_mi_y12_glm$confidence_interval),
      cbind(bd2_mi_y21_glm$output, bd2_mi_y21_glm$confidence_interval),
      cbind(bd2_mi_y22_glm$output, bd2_mi_y22_glm$confidence_interval),
      
      cbind(bd2_dr_y11_standard_0$output, bd2_dr_y11_standard_0$confidence_interval),
      cbind(bd2_dr_y12_standard_0$output, bd2_dr_y12_standard_0$confidence_interval),
      cbind(bd2_dr_y21_standard_0$output, bd2_dr_y21_standard_0$confidence_interval),
      cbind(bd2_dr_y22_standard_0$output, bd2_dr_y22_standard_0$confidence_interval),
      
      cbind(bd2_dr_y11_0$output, bd2_dr_y11_0$confidence_interval),
      cbind(bd2_dr_y12_0$output, bd2_dr_y12_0$confidence_interval),
      cbind(bd2_dr_y21_0$output, bd2_dr_y21_0$confidence_interval),
      cbind(bd2_dr_y22_0$output, bd2_dr_y22_0$confidence_interval)
    ))
  
  df2$bd <- 2
  
  result <- rbind(df1, df2)
  result$k <- k
  result
}

