## simulation with analytical standard errors

yang_sim <- function(k) {
  
  ### probability sample
  sample_prob <- pop_data[sample(1:N, n),]
  ### d-weights 
  sample_prob$w <- N/n
  
  ## create A matrix for quartiles
  sample_p1 <- joint_calib_create_matrix(sample_prob[, c("x1", "x2")], 
                                         N, 
                                         pop_quantiles = list(x1=quantile(sample_prob$x1, p_quantiles1), 
                                                              x2=quantile(sample_prob$x2, p_quantiles1)))
  ### create A matrix for deciles
  sample_p2 <- joint_calib_create_matrix(sample_prob[, c("x1", "x2")], 
                                         N, 
                                         pop_quantiles = list(x1=quantile(sample_prob$x1, p_quantiles2), 
                                                              x2=quantile(sample_prob$x2, p_quantiles2)))
  ### column names for A matrices
  colnames(sample_p1) <- c(paste0("x1_", p_quantiles1), paste0("x2_", p_quantiles1))
  colnames(sample_p2) <- c(paste0("x1_", p_quantiles2), paste0("x2_", p_quantiles2))
  
  ### add A matrix to probability sample
  sample_prob <- cbind(sample_prob, sample_p1, sample_p2)
  
  ### create survey design object
  sample_prob_svy <- svydesign(ids=~1, weights = ~w, data = sample_prob)
  
  ### estimate population quantiles and totals
  q1_est <- svyquantile( ~ x1 + x2, sample_prob_svy, p_quantiles1)
  q2_est <- svyquantile( ~ x1 + x2, sample_prob_svy, p_quantiles2)
  
  ### for comparison (quartiles and deciles x1, quartiles x2, deciles x2)
  p_comp <- c(q1_est$x1[,1],q2_est$x1[,1][-5], q1_est$x2[,1],q2_est$x2[,1][-5])
  p_comp_quants <- as.numeric(names(p_comp)[1:11])
  
  x_totals <- svytotal( ~ x1 + x2, sample_prob_svy)
  
  ### select nonprobabilistic samples
  sample_bd1 <- pop_data[rbinom(N,1,pop_data$p1)==1, ]
  sample_bd2 <- pop_data[rbinom(N,1,pop_data$p2)==1, ]
  
  ##########################################
  # linear inclusion --------------------------------------------------------
  ##########################################
  
  ### create A matrix for nonprobabilistic samples
  sample_np1 <- joint_calib_create_matrix(sample_bd1[, c("x1", "x2")], 
                                         sum(sample_prob$w), 
                                         pop_quantiles = list(x1=quantile(sample_prob$x1, p_quantiles1), 
                                                              x2=quantile(sample_prob$x2, p_quantiles1)))
  
  sample_np2 <- joint_calib_create_matrix(sample_bd1[, c("x1", "x2")], 
                                          sum(sample_prob$w),
                                         pop_quantiles = list(x1=quantile(sample_prob$x1, p_quantiles2), 
                                                              x2=quantile(sample_prob$x2, p_quantiles2)))
  
  ### add column names
  colnames(sample_np1) <- colnames(sample_p1)
  colnames(sample_np2) <- colnames(sample_p2)
  ### add A matrices to probability samples
  sample_bd1 <- cbind(sample_bd1, sample_np1, sample_np2[, setdiff(colnames(sample_np2), colnames(sample_np1))])
  
  ### IPW not_calibrated 
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
  
  ### quality metrics for IPW not calibrated
  bd1_ipw_standard_n <- c(sum(bd1_ipw_standard_0$weights), 
                          sum(bd1_ipw_standard_1a$weights), 
                          sum(bd1_ipw_standard_1b$weights), 
                          sum(bd1_ipw_standard_2a$weights), 
                          sum(bd1_ipw_standard_2b$weights))
  
  bd1_ipw_standard_tau <- c(
    sqrt(sum((colSums(sample_bd1*bd1_ipw_standard_0$weights)[1:2] - x_totals[1:2])^2)),
    sqrt(sum((colSums(sample_bd1*bd1_ipw_standard_1a$weights)[1:2] - x_totals[1:2])^2)),
    sqrt(sum((colSums(sample_bd1*bd1_ipw_standard_1b$weights)[1:2] - x_totals[1:2])^2)),
    sqrt(sum((colSums(sample_bd1*bd1_ipw_standard_2a$weights)[1:2] - x_totals[1:2])^2)),
    sqrt(sum((colSums(sample_bd1*bd1_ipw_standard_2b$weights)[1:2] - x_totals[1:2])^2))
  )
  
  bd1_ipw_standard_alphas_0 <-  c(weightedQuantile(x = sample_bd1$x1, probs = p_comp_quants, weights = bd1_ipw_standard_0$weights),
                                  weightedQuantile(x = sample_bd1$x2, probs = p_comp_quants, weights = bd1_ipw_standard_0$weights)) 
  bd1_ipw_standard_alphas_1a <- c(weightedQuantile(x = sample_bd1$x1, probs = p_comp_quants, weights = bd1_ipw_standard_1a$weights),
                                  weightedQuantile(x = sample_bd1$x2, probs = p_comp_quants, weights = bd1_ipw_standard_1a$weights))
  bd1_ipw_standard_alphas_1b <- c(weightedQuantile(x = sample_bd1$x1, probs = p_comp_quants, weights = bd1_ipw_standard_1b$weights),
                                  weightedQuantile(x = sample_bd1$x2, probs = p_comp_quants, weights = bd1_ipw_standard_1b$weights))
  bd1_ipw_standard_alphas_2a <- c(weightedQuantile(x = sample_bd1$x1, probs = p_comp_quants, weights = bd1_ipw_standard_2a$weights),
                                  weightedQuantile(x = sample_bd1$x2, probs = p_comp_quants, weights = bd1_ipw_standard_2a$weights))
  bd1_ipw_standard_alphas_2b <- c(weightedQuantile(x = sample_bd1$x1, probs = p_comp_quants, weights = bd1_ipw_standard_2b$weights),
                                  weightedQuantile(x = sample_bd1$x2, probs = p_comp_quants, weights = bd1_ipw_standard_2b$weights))
  
  bd1_ipw_standard_alphas <- c(
    sqrt(sum((bd1_ipw_standard_alphas_0  - p_comp)^2)),
    sqrt(sum((bd1_ipw_standard_alphas_1a - p_comp)^2)),
    sqrt(sum((bd1_ipw_standard_alphas_1b - p_comp)^2)),
    sqrt(sum((bd1_ipw_standard_alphas_2a - p_comp)^2)),
    sqrt(sum((bd1_ipw_standard_alphas_2b - p_comp)^2))
  )
  
  ## IPW calibrated 
  
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

  ### quality metrics for IPW calibrated
  bd1_ipw_n <- c(sum(bd1_ipw_0$weights), 
                 sum(bd1_ipw_1a$weights), 
                 sum(bd1_ipw_1b$weights), 
                 sum(bd1_ipw_2a$weights), 
                 sum(bd1_ipw_2b$weights))
  

  bd1_ipw_tau <- c(
    sqrt(sum((colSums(sample_bd1*bd1_ipw_0$weights)[1:2] - x_totals[1:2])^2)),
    sqrt(sum((colSums(sample_bd1*bd1_ipw_1a$weights)[1:2] - x_totals[1:2])^2)),
    sqrt(sum((colSums(sample_bd1*bd1_ipw_1b$weights)[1:2] - x_totals[1:2])^2)),
    sqrt(sum((colSums(sample_bd1*bd1_ipw_2a$weights)[1:2] - x_totals[1:2])^2)),
    sqrt(sum((colSums(sample_bd1*bd1_ipw_2b$weights)[1:2] - x_totals[1:2])^2))
  )
  
  bd1_ipw_alphas_0 <-  c(weightedQuantile(x = sample_bd1$x1, probs = p_comp_quants, weights = bd1_ipw_0$weights),
                                  weightedQuantile(x = sample_bd1$x2, probs = p_comp_quants, weights = bd1_ipw_0$weights)) 
  bd1_ipw_alphas_1a <- c(weightedQuantile(x = sample_bd1$x1, probs = p_comp_quants, weights = bd1_ipw_1a$weights),
                                  weightedQuantile(x = sample_bd1$x2, probs = p_comp_quants, weights = bd1_ipw_1a$weights))
  bd1_ipw_alphas_1b <- c(weightedQuantile(x = sample_bd1$x1, probs = p_comp_quants, weights = bd1_ipw_1b$weights),
                                  weightedQuantile(x = sample_bd1$x2, probs = p_comp_quants, weights = bd1_ipw_1b$weights))
  bd1_ipw_alphas_2a <- c(weightedQuantile(x = sample_bd1$x1, probs = p_comp_quants, weights = bd1_ipw_2a$weights),
                                  weightedQuantile(x = sample_bd1$x2, probs = p_comp_quants, weights = bd1_ipw_2a$weights))
  bd1_ipw_alphas_2b <- c(weightedQuantile(x = sample_bd1$x1, probs = p_comp_quants, weights = bd1_ipw_2b$weights),
                                  weightedQuantile(x = sample_bd1$x2, probs = p_comp_quants, weights = bd1_ipw_2b$weights))
  
  bd1_ipw_alphas <- c(
    sqrt(sum((bd1_ipw_alphas_0  - p_comp)^2)),
    sqrt(sum((bd1_ipw_alphas_1a - p_comp)^2)),
    sqrt(sum((bd1_ipw_alphas_1b - p_comp)^2)),
    sqrt(sum((bd1_ipw_alphas_2a - p_comp)^2)),
    sqrt(sum((bd1_ipw_alphas_2b - p_comp)^2))
  )
  
  
  ## mass imputation (nearest neighbours)
  bd1_mi <- nonprob(outcome = y11 + y12 + y21 + y22 ~ x1 + x2, 
                    svydesign = sample_prob_svy, data = sample_bd1, 
                    method_outcome = "nn", se = TRUE)
  
  ## mass imputation (glm)
  bd1_mi_c_glm <- nonprob(outcome = y11 + y12 ~ x1 + x2, svydesign = sample_prob_svy, data = sample_bd1, 
                            method_outcome = "glm", se = TRUE)
  bd1_mi_b_glm <- nonprob(outcome = y21 + y22 ~ x1 + x2, svydesign = sample_prob_svy, data = sample_bd1, 
                            method_outcome = "glm", family_outcome = "binomial", se = TRUE)
  
  ## doubly robust with IPW MLE 
  bd1_dr_c_standard_0 <- nonprob(selection = ~ x1 + x2, outcome = y11 + y12 ~ x1 + x2,  svydesign = sample_prob_svy, 
                                   data = sample_bd1,
                                   se = TRUE)
  
  bd1_dr_b_standard_0 <- nonprob(selection = ~ x1 + x2, outcome = y21 + y22 ~ x1 + x2, svydesign = sample_prob_svy,
                          data = sample_bd1,
                          family_outcome = "binomial",se = TRUE)

  
  ## doubly robust with IPW GEE 
  bd1_dr_c_0 <- nonprob(selection = ~ x1 + x2, outcome = y11 + y12 ~ x1 + x2,  svydesign = sample_prob_svy, 
                          data = sample_bd1,
                          se = TRUE,
                          control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd1_dr_b_0 <- nonprob(selection = ~ x1 + x2, outcome = y21 + y22 ~ x1 + x2, svydesign = sample_prob_svy,
                          data = sample_bd1,
                          family_outcome = "binomial",se = TRUE,
                          control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  
  ############################################
  # non-linear inclusion (the same comments applies as above)
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
  
  ### IPW not_calibrated 
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
  
  ### quality metrics for IPW not calibrated
  bd2_ipw_standard_n <- c(sum(bd2_ipw_standard_0$weights), 
                          sum(bd2_ipw_standard_1a$weights), 
                          sum(bd2_ipw_standard_1b$weights), 
                          sum(bd2_ipw_standard_2a$weights), 
                          sum(bd2_ipw_standard_2b$weights))
  
  bd2_ipw_standard_tau <- c(
    sqrt(sum((colSums(sample_bd2*bd2_ipw_standard_0$weights)[1:2] - x_totals[1:2])^2)),
    sqrt(sum((colSums(sample_bd2*bd2_ipw_standard_1a$weights)[1:2] - x_totals[1:2])^2)),
    sqrt(sum((colSums(sample_bd2*bd2_ipw_standard_1b$weights)[1:2] - x_totals[1:2])^2)),
    sqrt(sum((colSums(sample_bd2*bd2_ipw_standard_2a$weights)[1:2] - x_totals[1:2])^2)),
    sqrt(sum((colSums(sample_bd2*bd2_ipw_standard_2b$weights)[1:2] - x_totals[1:2])^2))
  )
  
  bd2_ipw_standard_alphas_0 <-  c(weightedQuantile(x = sample_bd2$x1, probs = p_comp_quants, weights = bd2_ipw_standard_0$weights),
                                  weightedQuantile(x = sample_bd2$x2, probs = p_comp_quants, weights = bd2_ipw_standard_0$weights)) 
  bd2_ipw_standard_alphas_1a <- c(weightedQuantile(x = sample_bd2$x1, probs = p_comp_quants, weights = bd2_ipw_standard_1a$weights),
                                  weightedQuantile(x = sample_bd2$x2, probs = p_comp_quants, weights = bd2_ipw_standard_1a$weights))
  bd2_ipw_standard_alphas_1b <- c(weightedQuantile(x = sample_bd2$x1, probs = p_comp_quants, weights = bd2_ipw_standard_1b$weights),
                                  weightedQuantile(x = sample_bd2$x2, probs = p_comp_quants, weights = bd2_ipw_standard_1b$weights))
  bd2_ipw_standard_alphas_2a <- c(weightedQuantile(x = sample_bd2$x1, probs = p_comp_quants, weights = bd2_ipw_standard_2a$weights),
                                  weightedQuantile(x = sample_bd2$x2, probs = p_comp_quants, weights = bd2_ipw_standard_2a$weights))
  bd2_ipw_standard_alphas_2b <- c(weightedQuantile(x = sample_bd2$x1, probs = p_comp_quants, weights = bd2_ipw_standard_2b$weights),
                                  weightedQuantile(x = sample_bd2$x2, probs = p_comp_quants, weights = bd2_ipw_standard_2b$weights))
  
  bd2_ipw_standard_alphas <- c(
    sqrt(sum((bd2_ipw_standard_alphas_0  - p_comp)^2)),
    sqrt(sum((bd2_ipw_standard_alphas_1a - p_comp)^2)),
    sqrt(sum((bd2_ipw_standard_alphas_1b - p_comp)^2)),
    sqrt(sum((bd2_ipw_standard_alphas_2a - p_comp)^2)),
    sqrt(sum((bd2_ipw_standard_alphas_2b - p_comp)^2))
  )
  
  ## IPW calibrated 
  
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
  
  ### quality metrics for IPW calibrated
  bd2_ipw_n <- c(sum(bd2_ipw_0$weights), 
                 sum(bd2_ipw_1a$weights), 
                 sum(bd2_ipw_1b$weights), 
                 sum(bd2_ipw_2a$weights), 
                 sum(bd2_ipw_2b$weights))
  
  
  bd2_ipw_tau <- c(
    sqrt(sum((colSums(sample_bd2*bd2_ipw_0$weights)[1:2] - x_totals[1:2])^2)),
    sqrt(sum((colSums(sample_bd2*bd2_ipw_1a$weights)[1:2] - x_totals[1:2])^2)),
    sqrt(sum((colSums(sample_bd2*bd2_ipw_1b$weights)[1:2] - x_totals[1:2])^2)),
    sqrt(sum((colSums(sample_bd2*bd2_ipw_2a$weights)[1:2] - x_totals[1:2])^2)),
    sqrt(sum((colSums(sample_bd2*bd2_ipw_2b$weights)[1:2] - x_totals[1:2])^2))
  )
  
  bd2_ipw_alphas_0 <-  c(weightedQuantile(x = sample_bd2$x1, probs = p_comp_quants, weights = bd2_ipw_0$weights),
                         weightedQuantile(x = sample_bd2$x2, probs = p_comp_quants, weights = bd2_ipw_0$weights)) 
  bd2_ipw_alphas_1a <- c(weightedQuantile(x = sample_bd2$x1, probs = p_comp_quants, weights = bd2_ipw_1a$weights),
                         weightedQuantile(x = sample_bd2$x2, probs = p_comp_quants, weights = bd2_ipw_1a$weights))
  bd2_ipw_alphas_1b <- c(weightedQuantile(x = sample_bd2$x1, probs = p_comp_quants, weights = bd2_ipw_1b$weights),
                         weightedQuantile(x = sample_bd2$x2, probs = p_comp_quants, weights = bd2_ipw_1b$weights))
  bd2_ipw_alphas_2a <- c(weightedQuantile(x = sample_bd2$x1, probs = p_comp_quants, weights = bd2_ipw_2a$weights),
                         weightedQuantile(x = sample_bd2$x2, probs = p_comp_quants, weights = bd2_ipw_2a$weights))
  bd2_ipw_alphas_2b <- c(weightedQuantile(x = sample_bd2$x1, probs = p_comp_quants, weights = bd2_ipw_2b$weights),
                         weightedQuantile(x = sample_bd2$x2, probs = p_comp_quants, weights = bd2_ipw_2b$weights))
  
  bd2_ipw_alphas <- c(
    sqrt(sum((bd2_ipw_alphas_0  - p_comp)^2)),
    sqrt(sum((bd2_ipw_alphas_1a - p_comp)^2)),
    sqrt(sum((bd2_ipw_alphas_1b - p_comp)^2)),
    sqrt(sum((bd2_ipw_alphas_2a - p_comp)^2)),
    sqrt(sum((bd2_ipw_alphas_2b - p_comp)^2))
  )
  
  ## mass imputation (nearest neighbours)
  bd2_mi <- nonprob(outcome = y11 + y12 + y21 + y22 ~ x1 + x2, svydesign = sample_prob_svy, data = sample_bd2, 
                    method_outcome = "nn", se = TRUE)
  
  ## mass imputation (glm)
  bd2_mi_c_glm <- nonprob(outcome = y11 + y12 ~ x1 + x2, svydesign = sample_prob_svy, data = sample_bd2, 
                          method_outcome = "glm", se = TRUE)
  bd2_mi_b_glm <- nonprob(outcome = y21 + y22 ~ x1 + x2, svydesign = sample_prob_svy, data = sample_bd2, 
                          method_outcome = "glm", family_outcome = "binomial", se = TRUE)
  
  ## doubly robust with IPW MLE 
  bd2_dr_c_standard_0 <- nonprob(selection = ~ x1 + x2, outcome = y11 + y12 ~ x1 + x2,  svydesign = sample_prob_svy, 
                                 data = sample_bd2,
                                 se = TRUE)
  
  bd2_dr_b_standard_0 <- nonprob(selection = ~ x1 + x2, outcome = y21 + y22 ~ x1 + x2, svydesign = sample_prob_svy,
                                 data = sample_bd2,
                                 family_outcome = "binomial",se = TRUE)
  
  
  ## doubly robust with IPW GEE 
  bd2_dr_c_0 <- nonprob(selection = ~ x1 + x2, outcome = y11 + y12 ~ x1 + x2,  svydesign = sample_prob_svy, 
                        data = sample_bd2,
                        se = TRUE,
                        control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  bd2_dr_b_0 <- nonprob(selection = ~ x1 + x2, outcome = y21 + y22 ~ x1 + x2, svydesign = sample_prob_svy,
                        data = sample_bd2,
                        family_outcome = "binomial",se = TRUE,
                        control_selection = controlSel(est_method_sel = "gee", h = 1))
  
  # save results ------------------------------------------------------------
  
  df1 <- data.frame(est = c(
    rep("naive", 4),
    
    rep("ipwst0", 4), 
    rep("ipwst1a", 4), rep("ipwst1b", 4), rep("ipwst2a", 4), rep("ipwst2b", 4),
    rep("ipw0", 4), 
    rep("ipw1a", 4), rep("ipw1b", 4), rep("ipw2a", 4), rep("ipw2b", 4),
    rep("mi_nn", 4), 
    rep("mi_glm", 4),
    rep("drst", 4),
    rep("dr", 4)),
    y = rep(c("y11", "y12", "y21", "y22"), times = 15),
    rbind(
      data.frame(mean = unlist(sample_bd1[, lapply(.SD, mean), .SDcols = patterns("y")]), 
                 SE=NA, lower_bound=NA, upper_bound = NA, 
                 N = NA, taus = NA, alphas = NA),
      
      cbind(bd1_ipw_standard_0$output,       bd1_ipw_standard_0$confidence_interval,  N=bd1_ipw_standard_n[1], taus = bd1_ipw_standard_tau[1], alphas = bd1_ipw_standard_alphas[1]),
      cbind(bd1_ipw_standard_1a$output,      bd1_ipw_standard_1a$confidence_interval, N=bd1_ipw_standard_n[2], taus = bd1_ipw_standard_tau[2], alphas = bd1_ipw_standard_alphas[2]),
      cbind(bd1_ipw_standard_1b$output,      bd1_ipw_standard_1b$confidence_interval, N=bd1_ipw_standard_n[3], taus = bd1_ipw_standard_tau[3], alphas = bd1_ipw_standard_alphas[3]),
      cbind(bd1_ipw_standard_2a$output,      bd1_ipw_standard_2a$confidence_interval, N=bd1_ipw_standard_n[4], taus = bd1_ipw_standard_tau[4], alphas = bd1_ipw_standard_alphas[4]),
      cbind(bd1_ipw_standard_2b$output,      bd1_ipw_standard_2b$confidence_interval, N=bd1_ipw_standard_n[5], taus = bd1_ipw_standard_tau[5], alphas = bd1_ipw_standard_alphas[5]),
      
      cbind(bd1_ipw_0$output,       bd1_ipw_0$confidence_interval, N=bd1_ipw_n[1], taus = bd1_ipw_tau[1], alphas = bd1_ipw_alphas[1]),
      cbind(bd1_ipw_1a$output,      bd1_ipw_1a$confidence_interval,N=bd1_ipw_n[2], taus = bd1_ipw_tau[2], alphas = bd1_ipw_alphas[2]),
      cbind(bd1_ipw_1b$output,      bd1_ipw_1b$confidence_interval,N=bd1_ipw_n[3], taus = bd1_ipw_tau[3], alphas = bd1_ipw_alphas[3]),
      cbind(bd1_ipw_2a$output,      bd1_ipw_2a$confidence_interval,N=bd1_ipw_n[4], taus = bd1_ipw_tau[4], alphas = bd1_ipw_alphas[4]),
      cbind(bd1_ipw_2b$output,      bd1_ipw_2b$confidence_interval,N=bd1_ipw_n[5], taus = bd1_ipw_tau[5], alphas = bd1_ipw_alphas[5]),
      
      cbind(bd1_mi$output,     bd1_mi$confidence_interval,N = NA, taus = NA, alphas = NA),
      cbind(bd1_mi_c_glm$output, bd1_mi_c_glm$confidence_interval,N = NA, taus = NA, alphas = NA),
      cbind(bd1_mi_b_glm$output, bd1_mi_b_glm$confidence_interval,N = NA, taus = NA, alphas = NA),
      
      cbind(bd1_dr_c_standard_0$output, bd1_dr_c_standard_0$confidence_interval,N = NA, taus = NA, alphas = NA),
      cbind(bd1_dr_b_standard_0$output, bd1_dr_b_standard_0$confidence_interval,N = NA, taus = NA, alphas = NA),
      
      cbind(bd1_dr_c_0$output, bd1_dr_c_0$confidence_interval,N = NA, taus = NA, alphas = NA),
      cbind(bd1_dr_b_0$output, bd1_dr_b_0$confidence_interval,N = NA, taus = NA, alphas = NA)
      ))
  
  df1$bd <- 1
  
  df2 <- data.frame(est = c(
    rep("naive", 4),
    
    rep("ipwst0", 4), 
    rep("ipwst1a", 4), rep("ipwst1b", 4), rep("ipwst2a", 4), rep("ipwst2b", 4),
    rep("ipw0", 4), 
    rep("ipw1a", 4), rep("ipw1b", 4), rep("ipw2a", 4), rep("ipw2b", 4),
    rep("mi_nn", 4), 
    rep("mi_glm", 4),
    rep("drst", 4),
    rep("dr", 4)),
    y = rep(c("y11", "y12", "y21", "y22"), times = 15),
    rbind(
      data.frame(mean = unlist(sample_bd1[, lapply(.SD, mean), .SDcols = patterns("y")]), 
                 SE=NA, lower_bound=NA, upper_bound = NA,
                 N = NA, taus = NA, alphas = NA),
      
      cbind(bd2_ipw_standard_0$output,       bd2_ipw_standard_0$confidence_interval, N=bd2_ipw_standard_n[1], taus = bd2_ipw_standard_tau[1], alphas = bd2_ipw_standard_alphas[1]),
      cbind(bd2_ipw_standard_1a$output,      bd2_ipw_standard_1a$confidence_interval,N=bd2_ipw_standard_n[2], taus = bd2_ipw_standard_tau[2], alphas = bd2_ipw_standard_alphas[2]),
      cbind(bd2_ipw_standard_1b$output,      bd2_ipw_standard_1b$confidence_interval,N=bd2_ipw_standard_n[3], taus = bd2_ipw_standard_tau[3], alphas = bd2_ipw_standard_alphas[3]),
      cbind(bd2_ipw_standard_2a$output,      bd2_ipw_standard_2a$confidence_interval,N=bd2_ipw_standard_n[4], taus = bd2_ipw_standard_tau[4], alphas = bd2_ipw_standard_alphas[4]),
      cbind(bd2_ipw_standard_2b$output,      bd2_ipw_standard_2b$confidence_interval,N=bd2_ipw_standard_n[5], taus = bd2_ipw_standard_tau[5], alphas = bd2_ipw_standard_alphas[5]),
      
      cbind(bd2_ipw_0$output,       bd2_ipw_0$confidence_interval,  N=bd2_ipw_n[1], taus = bd2_ipw_tau[1], alphas = bd2_ipw_alphas[1]),
      cbind(bd2_ipw_1a$output,      bd2_ipw_1a$confidence_interval, N=bd2_ipw_n[2], taus = bd2_ipw_tau[2], alphas = bd2_ipw_alphas[2]),
      cbind(bd2_ipw_1b$output,      bd2_ipw_1b$confidence_interval, N=bd2_ipw_n[3], taus = bd2_ipw_tau[3], alphas = bd2_ipw_alphas[3]),
      cbind(bd2_ipw_2a$output,      bd2_ipw_2a$confidence_interval, N=bd2_ipw_n[4], taus = bd2_ipw_tau[4], alphas = bd2_ipw_alphas[4]),
      cbind(bd2_ipw_2b$output,      bd2_ipw_2b$confidence_interval, N=bd2_ipw_n[5], taus = bd2_ipw_tau[5], alphas = bd2_ipw_alphas[5]),
      
      cbind(bd2_mi$output,     bd2_mi$confidence_interval,N = NA, taus = NA, alphas = NA),
      
      cbind(bd2_mi_c_glm$output, bd2_mi_c_glm$confidence_interval,N = NA, taus = NA, alphas = NA),
      cbind(bd2_mi_b_glm$output, bd2_mi_b_glm$confidence_interval,N = NA, taus = NA, alphas = NA),
    
      cbind(bd2_dr_c_standard_0$output, bd2_dr_c_standard_0$confidence_interval,N = NA, taus = NA, alphas = NA),
      cbind(bd2_dr_b_standard_0$output, bd2_dr_b_standard_0$confidence_interval,N = NA, taus = NA, alphas = NA),
      
      cbind(bd2_dr_c_0$output, bd2_dr_c_0$confidence_interval,N = NA, taus = NA, alphas = NA),
      cbind(bd2_dr_b_0$output, bd2_dr_b_0$confidence_interval,N = NA, taus = NA, alphas = NA)
    ))
  
  df2$bd <- 2
  
  result <- rbind(df1, df2)
  result$k <- k
  result
}

