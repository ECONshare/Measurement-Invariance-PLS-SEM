# The following code replicates the empirical example in the paper
# "Measurement invariance testin in partial least squares structural equation modeling"
# For custom applications contact the Author on email: benlien@econ.au.dk

# To run this example, install the following packages
# install.packages(c("usethis", "testthat", "roxygen2", "pkgbuild"))
# install.packages("devtools")
library(devtools)
# install_bitbucket("richierocks/assertive")
# install_github("mronkko/matrixpls")
# install.packages("tidyverse")
# install.packages("lubridate")
# install.packages("Matrix")
# install.packages("car")

source(file = "functions.R")
# Load data ---------------------------------------------------------------
# Read cleaned data
treat_group <- read_csv(file = "treatment.csv")
control_group <- read_csv(file = "control.csv")
# MI - full -----------------------------------------------------------
# t = 3 is pre-intervention control group
dat_long <- transform_to_long(treat_group, control_group, type = "full")
# Specify null hypothesis 
load_null_ref1 <-  matrix(0, nrow = 12, 
                          dimnames = list(c("LV_in_k1_t1", "LV_in_k1_t2", "LV_in_k1_t4",
                                            "LV_in_k2_t1", "LV_in_k2_t2", "LV_in_k2_t4",
                                            "LV_in_k3_t1", "LV_in_k3_t2", "LV_in_k3_t4",
                                            "LV_in_k4_t1", "LV_in_k4_t2", "LV_in_k4_t4")))

# Intercept null hyp: time 3 is reference time and k3 is marker variable
intercept_null_ref1 <-  matrix(0, nrow = 9, 
                               dimnames = list(c("LV_in_k1_t1", "LV_in_k1_t2", "LV_in_k1_t4",
                                                 "LV_in_k2_t1", "LV_in_k2_t2", "LV_in_k2_t4",
                                                 "LV_in_k4_t1", "LV_in_k4_t2", "LV_in_k4_t4")))

MI_res_full <- MI_block(sepDataOrg = dat_long$dat_long[,c("id", "time","POSR1", "POSR2" ,"POSR3", "POSR4", "POSR")],
                        SepIndNamesOrg = c("POSR1", "POSR2", "POSR3", "POSR4"),
                        SepLVNameOrg = c("POSR"),
                        ManRobCov = F,
                        load_null = load_null_ref1,
                        intercept_null = intercept_null_ref1,
                        marker_nam = "POSR3",
                        reference_time = 3)

# New loading and intercept null hypothesis
load_null_boot <- as.matrix(MI_res_full$Loadings)
intercept_null_boot <- as.matrix(MI_res_full$Intercepts)

set.seed(42)
B <- 5000
boot_res <- matrix(0, ncol = 2, nrow = B)
colnames(boot_res) <- c("chisq_load", "chisq_intercept")
for (b in 1:B) {
  boot_treat <- treat_group[sample(1:nrow(treat_group),replace = TRUE),]
  boot_control <- control_group[sample(1:nrow(control_group),replace = TRUE),]
  boot_dat_long <- transform_to_long(treat = boot_treat, control = boot_control, type = "full")
  tmp_res <- MI_block(sepDataOrg = boot_dat_long$dat_long[,c("id", "time","POSR1", "POSR2" ,"POSR3", "POSR4", "POSR")],
                      SepIndNamesOrg = c("POSR1", "POSR2", "POSR3", "POSR4"),
                      SepLVNameOrg = c("POSR"),
                      ManRobCov = F,
                      load_null = load_null_boot,
                      intercept_null = intercept_null_boot,
                      marker_nam = "POSR3",
                      reference_time = 3)
  boot_res[b,c("chisq_load", "chisq_intercept")] <- c(tmp_res$TestLoadings_Chisq, tmp_res$TestIntercepts_Chisq)
}

boot_full <- boot_p_val(boot_Chisq = boot_res, B = B, Org_res = MI_res_full)
boot_full

# Test if K4 has significantly different intercept across groups ------------------------------
# t = 3 is pre-intervention control group
dat_long <- transform_to_long(treat_group, control_group, type = "full")

# Specify null hypothesis 
load_null_ref1 <-  matrix(0, nrow = 12, 
                          dimnames = list(c("LV_in_k1_t1", "LV_in_k1_t2", "LV_in_k1_t4",
                                            "LV_in_k2_t1", "LV_in_k2_t2", "LV_in_k2_t4",
                                            "LV_in_k3_t1", "LV_in_k3_t2", "LV_in_k3_t4",
                                            "LV_in_k4_t1", "LV_in_k4_t2", "LV_in_k4_t4")))

# Intercept null hyp: time 3 is reference time and k3 is marker variable. Test if K4 deviates across groups
intercept_null_ref1 <-  matrix(0, nrow = 9, 
                               dimnames = list(c("LV_in_k1_t1", "LV_in_k1_t2", "LV_in_k1_t4",
                                                 "LV_in_k2_t1", "LV_in_k2_t2", "LV_in_k2_t4",
                                                 "LV_in_k4_t1", "LV_in_k4_t2", "LV_in_k4_t4")))

intercept_null_for_test <-  matrix(0, nrow = 3, 
                                   dimnames = list(c("LV_in_k4_t1", "LV_in_k4_t2", "LV_in_k4_t4")))

MI_res_k4 <- MI_block_one_indicator_inter(sepDataOrg = dat_long$dat_long[,c("id", "time","POSR1", "POSR2" ,"POSR3", "POSR4", "POSR")],
                                          SepIndNamesOrg = c("POSR1", "POSR2", "POSR3", "POSR4"),
                                          SepLVNameOrg = c("POSR"),
                                          ManRobCov = F,
                                          load_null = load_null_ref1,
                                          intercept_null = intercept_null_ref1,
                                          intercept_null_for_test = intercept_null_for_test,
                                          marker_nam = "POSR3",
                                          reference_time = 3)

# New loading and intercept null hypothesis
load_null_boot <- as.matrix(MI_res_k4$Loadings)
intercept_null_boot <- as.matrix(MI_res_k4$Intercepts)
intercept_null_for_test <- as.matrix(MI_res_k4$Intercepts[7:9])

set.seed(42)
B <- 5000
boot_res <- matrix(0, ncol = 2, nrow = B)
colnames(boot_res) <- c("chisq_load", "chisq_intercept")
for (b in 1:B) {
  boot_treat <- treat_group[sample(1:nrow(treat_group),replace = TRUE),]
  boot_control <- control_group[sample(1:nrow(control_group),replace = TRUE),]
  boot_dat_long <- transform_to_long(treat = boot_treat, control = boot_control, type = "full")
  tmp_res <- MI_block_one_indicator_inter(sepDataOrg = boot_dat_long$dat_long[,c("id", "time","POSR1", "POSR2" ,"POSR3", "POSR4", "POSR")],
                                          SepIndNamesOrg = c("POSR1", "POSR2", "POSR3", "POSR4"),
                                          SepLVNameOrg = c("POSR"),
                                          ManRobCov = F,
                                          load_null = load_null_boot,
                                          intercept_null = intercept_null_boot,
                                          intercept_null_for_test = intercept_null_for_test,
                                          marker_nam = "POSR3",
                                          reference_time = 3)
  boot_res[b,c("chisq_load", "chisq_intercept")] <- c(tmp_res$TestLoadings_Chisq, tmp_res$TestIntercepts_Chisq)
}

boot_k4 <- boot_p_val(boot_Chisq = boot_res, B = B, Org_res = MI_res_full)
boot_k4


# MI exclude POSR4 -----------------------------------------------------------
dat_long <- transform_to_long(treat_group, control_group, type = "exclude_POSR4")
# Specify null hypothesis 
load_null_ref1 <-  matrix(0, nrow = 9, 
                          dimnames = list(c("LV_in_k1_t1", "LV_in_k1_t2", "LV_in_k1_t4",
                                            "LV_in_k2_t1", "LV_in_k2_t2", "LV_in_k2_t4",
                                            "LV_in_k3_t1", "LV_in_k3_t2", "LV_in_k3_t4")))

# Intercept null hyp: time 3 is reference time and k3 is marker variable
intercept_null_ref1 <-  matrix(0, nrow = 6, 
                               dimnames = list(c("LV_in_k1_t1", "LV_in_k1_t2", "LV_in_k1_t4",
                                                 "LV_in_k2_t1", "LV_in_k2_t2", "LV_in_k2_t4")))

dat_long <- transform_to_long(treat_group, control_group, type = "exclude_POSR4")

MI_res_full <- MI_block(sepDataOrg = dat_long$dat_long[,c("id", "time","POSR1", "POSR2", "POSR3", "POSR")],
                        SepIndNamesOrg = c("POSR1", "POSR2", "POSR3"),
                        SepLVNameOrg = c("POSR"),
                        ManRobCov = F,
                        load_null = load_null_ref1,
                        intercept_null = intercept_null_ref1,
                        marker_nam = "POSR3",
                        reference_time = 3)

# New loading and intercept null hypothesis
load_null_boot <- as.matrix(MI_res_full$Loadings)
intercept_null_boot <- as.matrix(MI_res_full$Intercepts)

set.seed(42)
B <- 5000
boot_res <- matrix(0, ncol = 2, nrow = B)
colnames(boot_res) <- c("chisq_load", "chisq_intercept")
for (b in 1:B) {
  boot_treat <- treat_group[sample(1:nrow(treat_group),replace = TRUE),]
  boot_control <- control_group[sample(1:nrow(control_group),replace = TRUE),]
  boot_dat_long <- transform_to_long(treat = boot_treat, control = boot_control, type = "exclude_POSR4")
  tmp_res <- MI_block(sepDataOrg = boot_dat_long$dat_long[,c("id", "time","POSR1", "POSR2", "POSR3", "POSR")],
                      SepIndNamesOrg = c("POSR1", "POSR2", "POSR3"),
                      SepLVNameOrg = c("POSR"),
                      ManRobCov = F,
                      load_null = load_null_boot,
                      intercept_null = intercept_null_boot,
                      marker_nam = "POSR3",
                      reference_time = 3)
  boot_res[b,c("chisq_load", "chisq_intercept")] <- c(tmp_res$TestLoadings_Chisq, tmp_res$TestIntercepts_Chisq)
}

boot_without_POSR4 <- boot_p_val(boot_Chisq = boot_res, B = B, Org_res = MI_res_full)
boot_without_POSR4
# Latent Analysis of Covariance -------------------------------------------
# Combine data for ANCOVA 
dat_ancova <- rbind(cbind(treat = 1, treat_group), 
                    cbind(treat = 0, control_group)) %>% 
  as_tibble()
# Remove POSR4
dat_ancova <- dat_ancova %>%
  dplyr::select(-POSR1_4, -POSR2_4)
# Estimate PLS model
matrixplsmodel <- pls_specification(model_type = "ancova", 
                                    LVnames = c("POSR_t1","POSR_t2"), 
                                    indicator_nam = c("POSR1_1", "POSR1_2", "POSR1_3","POSR2_1", "POSR2_2", "POSR2_3"))

pls_data <- matrixpls(S = cov(dat_ancova[,c("POSR1_1", "POSR1_2", "POSR1_3","POSR2_1", "POSR2_2", "POSR2_3")]), 
                      model = matrixplsmodel, 
                      parametersInner = estimator.ols, 
                      standardize = F)
score_data <- as.matrix(dat_ancova[,c("POSR1_1", "POSR1_2", "POSR1_3","POSR2_1", "POSR2_2", "POSR2_3")])%*%t(attr(pls_data, "W"))
# Relate scores to treatment indicator
dat_ancova_org <- cbind(dat_ancova, score_data)
# Pre-treatment deviation from mean
dat_ancova_org <- dat_ancova_org %>%
  dplyr::mutate(pre_dev_mean = POSR_t1 - mean(POSR_t1))
# ANCOVA regression
ancova_reg <- lm(POSR_t2 ~ treat*pre_dev_mean, data = dat_ancova_org) 
sum_ancova <- summary(ancova_reg)
# Testing hypothesis using Chisq test
org_ancova_chisq <- c(intercept = linearHypothesis(model = ancova_reg,
                                                   hypothesis.matrix = "(Intercept)=0",
                                                   test = c("Chisq"))$Chisq[2],
                      treat = linearHypothesis(model = ancova_reg,
                                               hypothesis.matrix = "treat=0",
                                               test = c("Chisq"))$Chisq[2],
                      pre_dev_mean = linearHypothesis(model = ancova_reg,
                                                      hypothesis.matrix = "pre_dev_mean=0",
                                                      test = c("Chisq"))$Chisq[2],
                      `treat:pre_dev_mean` = linearHypothesis(model = ancova_reg,
                                                              hypothesis.matrix = "treat:pre_dev_mean=0",
                                                              test = c("Chisq"))$Chisq[2])



# Bootstrapped ANCOVA
boot_null_hyp_intercept <- paste("(Intercept)",sum_ancova$coefficients[1], sep="=")
boot_null_hyp_treat <- paste("treat",sum_ancova$coefficients[2], sep="=")
boot_null_hyp_pre_dev_mean <- paste("pre_dev_mean",sum_ancova$coefficients[3], sep = "=")
boot_null_hyp_interact <- paste("treat:pre_dev_mean",sum_ancova$coefficients[4], sep = "=")

## Bootstrapping the ANCOVA using Chisq test statistic ----
set.seed(42)
B <-5000
boot_res <- matrix(0, nrow = B, ncol = 8)
colnames(boot_res) <- c("boot_intercept_chisq", "boot_intercept_coef",
                        "boot_treat_chisq", "boot_treat_coef",
                        "boot_pre_dev_mean_chisq", "boot_pre_dev_mean_coef",
                        "boot_interact_chisq", "boot_interact_coef")
for (b in 1:B) {
  # Bootstrap groups same size as original data size
  dat_ancova_boot_treat <- dat_ancova[dat_ancova$treat == 1,]
  dat_ancova_boot_control <- dat_ancova[dat_ancova$treat == 0,]
  dat_ancova_boot <- rbind(dat_ancova_boot_treat[sample(1:nrow(dat_ancova_boot_treat), size = nrow(dat_ancova_boot_treat), replace = T),],
                           dat_ancova_boot_control[sample(1:nrow(dat_ancova_boot_control), size = nrow(dat_ancova_boot_control), replace = T),])
  
  pls_data <- matrixpls(S = cov(dat_ancova_boot[,c("POSR1_1", "POSR1_2", "POSR1_3","POSR2_1", "POSR2_2", "POSR2_3")]), 
                        model = matrixplsmodel, 
                        parametersInner = estimator.ols, 
                        standardize = F)
  score_data <- as.matrix(dat_ancova_boot[,c("POSR1_1", "POSR1_2", "POSR1_3","POSR2_1", "POSR2_2", "POSR2_3")])%*%t(attr(pls_data, "W"))
  # Relate scores to treatment indicator
  dat_ancova_boot <- cbind(dat_ancova_boot, score_data)
  # Pre-treatment deviation from mean
  dat_ancova_boot <- dat_ancova_boot %>%
    dplyr::mutate(pre_dev_mean = POSR_t1 - mean(POSR_t1))
  # ANCOVA regression
  ancova_reg_boot <- lm(POSR_t2 ~ treat*pre_dev_mean, data = dat_ancova_boot) 
  sum_ancova_boot <- summary(ancova_reg_boot)
  # Testing hypothesis
  boot_res[b,"boot_intercept_chisq"] <- linearHypothesis(model = ancova_reg_boot,hypothesis.matrix = boot_null_hyp_intercept,
                                                         test = c("Chisq"))$Chisq[2]
  boot_res[b,"boot_intercept_coef"] <- sum_ancova_boot$coefficients[1,1]
  
  boot_res[b,"boot_treat_chisq"] <- linearHypothesis(model = ancova_reg_boot,hypothesis.matrix = boot_null_hyp_treat,
                                                     test = c("Chisq"))$Chisq[2]
  boot_res[b,"boot_treat_coef"] <- sum_ancova_boot$coefficients[2,1]
  
  boot_res[b,"boot_pre_dev_mean_chisq"] <- linearHypothesis(model = ancova_reg_boot,hypothesis.matrix = boot_null_hyp_pre_dev_mean,
                                                            test = c("Chisq"))$Chisq[2]
  boot_res[b,"boot_pre_dev_mean_coef"] <- sum_ancova_boot$coefficients[3,1]
  
  boot_res[b,"boot_interact_chisq"] <- linearHypothesis(model = ancova_reg_boot,hypothesis.matrix = boot_null_hyp_interact,
                                                        test = c("Chisq"))$Chisq[2]
  boot_res[b,"boot_interact_coef"] <- sum_ancova_boot$coefficients[4,1]
}
## P-values
mean(org_ancova_chisq["intercept"] <= boot_res[,"boot_intercept_chisq"])
mean(org_ancova_chisq["treat"] <= boot_res[,"boot_treat_chisq"])
mean(org_ancova_chisq["pre_dev_mean"] <= boot_res[,"boot_pre_dev_mean_chisq"])
mean(org_ancova_chisq["treat:pre_dev_mean"] <= boot_res[,"boot_interact_chisq"])
