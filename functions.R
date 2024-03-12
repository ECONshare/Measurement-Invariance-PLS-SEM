library(tidyverse)
library(lubridate)
require(devtools)
library(matrixpls)
library(Matrix)
library(dplyr)
library(car)


# Used in example ---------------------------------------------------------
# Calculate PLS-SEM and transform to long format
transform_to_long <- function(treat5, control5, type){
  if (type == "full") {
    # Estimate PLS models
    matrixpls_model <- pls_specification(model_type = "simple", LVnames = c("POSR_t1", "POSR_t2"),
                                         indicator_nam = c("POSR1_1", "POSR1_2", "POSR1_3", "POSR1_4",
                                                           "POSR2_1", "POSR2_2", "POSR2_3", "POSR2_4"))
    pls_treat <-matrixpls(cov(treat5),matrixpls_model,
                          parametersInner = estimator.ols, standardize = F)
    pls_control <-matrixpls(cov(control5),matrixpls_model,
                            parametersInner = estimator.ols, standardize = F)
    # Get LV scores
    lv_score_treat <- as.matrix(treat5) %*% t(attr(pls_treat, "W"))
    lv_score_control <- as.matrix(control5) %*% t(attr(pls_control, "W"))
    # Arrange scores in long format
    long_treat_LV1 <- cbind(id = 1:nrow(lv_score_treat), 
                            time = 1,
                            lv_score = lv_score_treat[, "POSR_t1" ],
                            treat5[,c("POSR1_1", "POSR1_2", "POSR1_3", "POSR1_4")])
    colnames(long_treat_LV1) <- c("id", "time", "POSR", "POSR1", "POSR2", "POSR3", "POSR4")
    long_treat_LV2 <- cbind(id = 1:nrow(lv_score_treat), 
                            time = 2,
                            lv_score = lv_score_treat[, "POSR_t2" ],
                            treat5[,c("POSR2_1", "POSR2_2", "POSR2_3", "POSR2_4")])
    colnames(long_treat_LV2) <- c("id", "time", "POSR", "POSR1", "POSR2", "POSR3", "POSR4")
    long_control_LV1 <- cbind(id = (nrow(lv_score_treat)+1):(nrow(lv_score_control)+nrow(lv_score_treat)), 
                              time = 3,
                              lv_score = lv_score_control[, "POSR_t1" ],
                              control5[,c("POSR1_1", "POSR1_2", "POSR1_3", "POSR1_4")])
    colnames(long_control_LV1) <- c("id", "time", "POSR", "POSR1", "POSR2", "POSR3", "POSR4")
    long_control_LV2 <- cbind(id = (nrow(lv_score_treat)+1):(nrow(lv_score_control)+nrow(lv_score_treat)), 
                              time = 4,
                              lv_score = lv_score_control[, "POSR_t2" ],
                              control5[,c("POSR2_1", "POSR2_2", "POSR2_3", "POSR2_4")])
    colnames(long_control_LV2) <- c("id", "time", "POSR", "POSR1", "POSR2", "POSR3", "POSR4")
    
    dat_long = rbind(long_treat_LV1, long_treat_LV2, long_control_LV1, long_control_LV2)
  }
  if (type == "exclude_POSR4") {
    treat5_reduced <- treat5[, c("POSR1_1", "POSR1_2", "POSR1_3",
                                 "POSR2_1", "POSR2_2", "POSR2_3")]
    control5_reduced <- control5[, c("POSR1_1", "POSR1_2", "POSR1_3",
                                     "POSR2_1", "POSR2_2", "POSR2_3")]
    # Estimate PLS models
    matrixpls_model <- pls_specification(model_type = "reduced", LVnames = c("POSR_t1", "POSR_t2"),
                                         indicator_nam = c("POSR1_1", "POSR1_2", "POSR1_3",
                                                           "POSR2_1", "POSR2_2", "POSR2_3"))
    pls_treat <-matrixpls(cov(treat5_reduced),matrixpls_model,
                          parametersInner = estimator.ols, standardize = F)
    pls_control <-matrixpls(cov(control5_reduced),matrixpls_model,
                            parametersInner = estimator.ols, standardize = F)
    # Get LV scores
    lv_score_treat <- as.matrix(treat5_reduced) %*% t(attr(pls_treat, "W"))
    lv_score_control <- as.matrix(control5_reduced) %*% t(attr(pls_control, "W"))
    # Arrange scores in long format
    long_treat_LV1 <- cbind(id = 1:nrow(lv_score_treat), 
                            time = 1,
                            lv_score = lv_score_treat[, "POSR_t1" ],
                            treat5[,c("POSR1_1", "POSR1_2", "POSR1_3")])
    colnames(long_treat_LV1) <- c("id", "time", "POSR", "POSR1", "POSR2", "POSR3")
    long_treat_LV2 <- cbind(id = 1:nrow(lv_score_treat), 
                            time = 2,
                            lv_score = lv_score_treat[, "POSR_t2" ],
                            treat5[,c("POSR2_1", "POSR2_2", "POSR2_3")])
    colnames(long_treat_LV2) <- c("id", "time", "POSR", "POSR1", "POSR2", "POSR3")
    long_control_LV1 <- cbind(id = (nrow(lv_score_treat)+1):(nrow(lv_score_control)+nrow(lv_score_treat)), 
                              time = 3,
                              lv_score = lv_score_control[, "POSR_t1" ],
                              control5[,c("POSR1_1", "POSR1_2", "POSR1_3")])
    colnames(long_control_LV1) <- c("id", "time", "POSR", "POSR1", "POSR2", "POSR3")
    long_control_LV2 <- cbind(id = (nrow(lv_score_treat)+1):(nrow(lv_score_control)+nrow(lv_score_treat)), 
                              time = 4,
                              lv_score = lv_score_control[, "POSR_t2" ],
                              control5[,c("POSR2_1", "POSR2_2", "POSR2_3")])
    colnames(long_control_LV2) <- c("id", "time", "POSR", "POSR1", "POSR2", "POSR3")
    
    dat_long = rbind(long_treat_LV1, long_treat_LV2, long_control_LV1, long_control_LV2)
  }
  return(list(dat_long = dat_long, matrixpls_model = matrixpls_model))
}
## PLS-SEM specifications seperate models for each group
pls_specification <- function(model_type, LVnames, indicator_nam){
  if (model_type == "simple") {
    # Inner model
    LV1 = c(0,0)
    LV2 = c(1,0)
    inner = rbind(LV1,LV2)
    colnames(inner) <- rownames(inner) <- LVnames
    
    # Reflective
    reflective <- matrix(0,nrow = 8, ncol =2)
    rownames(reflective)<-indicator_nam
    colnames(reflective)<-colnames(inner)
    reflective[1:4,1]<-1
    reflective[5:8,2]<-1
    # Formative
    formative <- t(reflective)*0
    
    matrixplsmodel <- list(inner = inner,reflective = reflective,
                           formative = formative)
    return(matrixplsmodel)
  }
  if (model_type == "reduced") {
    # Inner model
    LV1 = c(0,0)
    LV2 = c(1,0)
    inner = rbind(LV1,LV2)
    colnames(inner) <- rownames(inner) <- LVnames
    
    # Reflective
    reflective <- matrix(0,nrow = 6, ncol =2)
    rownames(reflective)<-indicator_nam
    colnames(reflective)<-colnames(inner)
    reflective[1:3,1]<-1
    reflective[4:6,2]<-1
    # Formative
    formative <- t(reflective)*0
    
    matrixplsmodel <- list(inner = inner,reflective = reflective,
                           formative = formative)
    return(matrixplsmodel)
  }
  if (model_type == "reduced1") {
    # Inner model
    LV1 = c(0,0)
    LV2 = c(1,0)
    inner = rbind(LV1,LV2)
    colnames(inner) <- rownames(inner) <- LVnames
    
    # Reflective
    reflective <- matrix(0,nrow = 4, ncol =2)
    rownames(reflective)<-indicator_nam
    colnames(reflective)<-colnames(inner)
    reflective[1:2,1]<-1
    reflective[3:4,2]<-1
    # Formative
    formative <- t(reflective)*0
    
    matrixplsmodel <- list(inner = inner,reflective = reflective,
                           formative = formative)
    return(matrixplsmodel)
  }
  if (model_type == "ancova") {
    # Inner model
    LV1 = c(0,0)
    LV2 = c(1,0)
    inner = rbind(LV1,LV2)
    colnames(inner) <- rownames(inner) <- LVnames
    # Reflective
    reflective <- matrix(0,nrow = 6, ncol =2)
    rownames(reflective)<-indicator_nam
    colnames(reflective)<-colnames(inner)
    reflective[1:3,1]<-1
    reflective[4:6,2]<-1
    # Formative
    formative <- t(reflective)*0
    # Combined model
    matrixplsmodel <- list(inner = inner,reflective = reflective,
                           formative = formative)
    return(matrixplsmodel)
  }
  
}
# MI-block's
MI_block<-function(sepDataOrg,SepIndNamesOrg,SepLVNameOrg,ManRobCov,load_null,intercept_null, marker_nam, reference_time){
  # To meet assumption of equal LV variance across groups: scale latent variable to have unit variance
  sepDataOrg[,SepLVNameOrg] <- scale(sepDataOrg[,SepLVNameOrg], center = F, scale = sqrt(var(sepDataOrg[,SepLVNameOrg])))
  # Make variable names with time indices
  all_names_list<-TimeIndicesVarNames(SepLVNameOrg=SepLVNameOrg,SepIndNamesOrg=SepIndNamesOrg,sepDataOrg=sepDataOrg)
  # Original indicator and LV names
  OrgNames<-colnames(sepDataOrg)[-c(1,2)]
  # Number of observations and time periods
  N_each <- split(sepDataOrg[,"id"],sepDataOrg[,"time"])
  N<-sum(unlist(lapply(N_each, function(x) {length(x)})))
  time<-length(unique(sepDataOrg[,"time"]))
  # Mean center for loadings null hypothesis
  sepDataOrg1 <- sepDataOrg
  sepDataOrg1_tmp <- vector(mode = "list", length = time)
  
  for (t in 1:time){
    sepDataOrg1_tmp[[t]] <- as.data.frame(lapply(sepDataOrg1[sepDataOrg1$time == t,-c(1,2)], function(x) {scale(x, center = T, scale = FALSE)} ))
  }
  sepDataOrg1[-c(1,2)] <- do.call(rbind,sepDataOrg1_tmp)
  ################## Test Loadings
  # Create matrix with Y_star in diagonal
  # Split LV by group/time
  split_lv_by_time <- lapply(split(sepDataOrg1[,SepLVNameOrg], sepDataOrg1[,"time"]), function(x){x <- as.matrix(x)})
  # Create block diagonal matrix with LV scores on diagonal
  Y_star <- bdiag(split_lv_by_time)
  # Put all LV scores in first column
  Y_star[,reference_time] <- sepDataOrg1[,SepLVNameOrg]
  # Make block diagonal matrix with as many blocks as indicators
  var_dum_list <- vector(mode = "list", length = length(SepIndNamesOrg))
  Y_full_matrix <- bdiag(lapply(var_dum_list, function(x) {x <- Y_star}))
  # Naming the block diagonal matrix columns using the specified loading null hypothesis matrix
  split_null_hyp_names <- split(rownames(load_null), ceiling(seq_along(rownames(load_null))/(time-1)))
  tmp_test_names <- Y_names_from_hypothesis(split_null_hyp_names, reference_time, time)
  colnames(Y_full_matrix) <- tmp_test_names
  # Stack indicators
  stacked_indicators <- stack(sepDataOrg1[,SepIndNamesOrg])[,1]
  # Regression
  reg_load <- lm(stacked_indicators ~ 0 + . ,data = as.data.frame(as.matrix(Y_full_matrix)))
  # Testing loadings (metric invariance)
  load_null_1 <- paste(rownames(load_null), load_null[,1], sep = "=")
  Hyp_load<-linearHypothesis(reg_load, hypothesis.matrix = load_null_1,test=c("Chisq"))
  p_val_load<-Hyp_load$`Pr(>Chisq)`[2]
  Chisq_load<-Hyp_load$Chisq[2]
  
  ################## Test intercepts
  # Getting marker variable
  marker_nam_load <- lm(sepDataOrg[,marker_nam] ~ sepDataOrg[,SepLVNameOrg])$coefficients[2]
  marker_var <- sepDataOrg[,c("time",marker_nam)]
  # fill data matrix for lambda
  for_anova <- data.frame(dependent = rep(0,N*(length(SepIndNamesOrg)-1)), ind_name = "empty", stringsAsFactors = F)
  index_begin <- 1
  index_end <- N
  non_mark_nam <- SepIndNamesOrg[SepIndNamesOrg!=marker_nam]
  for (ind_nam in non_mark_nam) {
    # Loading restricted to be equal over groups
    tmp_load <- lm(sepDataOrg[,ind_nam] ~ sepDataOrg[,SepLVNameOrg])$coefficients[2]
    # Calculating lambda
    for_anova[index_begin:index_end,"dependent"] <- sepDataOrg[,ind_nam] - (tmp_load/marker_nam_load) * sepDataOrg[,marker_nam]
    for_anova[index_begin:index_end,"ind_name"] <- rep(ind_nam,N)
    index_begin <- index_begin + N
    index_end <- index_end + N
  }
  # Dummy variable matrix
  time_dum_list <- vector(mode = "list", length =time)
  samp_each_group <- N
  D <- bdiag(lapply(N_each, function(x) as.matrix(rep(1,length(x)))))
  # Reference group
  D[,reference_time] <- reference_time
  var_dum_list <- vector(mode = "list", length = length(SepIndNamesOrg)-1)
  D_full_mat <- bdiag(lapply(var_dum_list, function(x) {x <- D}))
  
  split_null_hyp_names <- split(rownames(intercept_null), ceiling(seq_along(rownames(intercept_null))/(time-1)))
  # aux variables is variables for which the coefficients are not reported (only the differences is reported)
  tmp_test_names <- Y_names_from_hypothesis(split_null_hyp_names, reference_time, time)
  colnames(D_full_mat) <- tmp_test_names
  # Intercept testing regression
  intercept_reg <- lm(for_anova$dependent ~ 0 + ., data= as.data.frame(as.matrix(D_full_mat)))
  # Test hypothesis
  intercept_null_1 <- paste(rownames(intercept_null), intercept_null[,1], sep = "=")
  Hyp_inter<-linearHypothesis(intercept_reg, hypothesis.matrix = intercept_null_1,test=c("Chisq"))
  p_val_inter<-Hyp_inter$`Pr(>Chisq)`[2]
  Chisq_inter<-Hyp_inter$Chisq[2]
  
  # Output
  out<-list(Intercepts=intercept_reg$coefficients[rownames(intercept_null)],
            TestIntercepts_Chisq=Chisq_inter,
            TestIntercepts_p=p_val_inter,
            Loadings=reg_load$coefficients[rownames(load_null)],
            TestLoadings_Chisq=Chisq_load,
            TestLoadings_p=p_val_load)
  return(out)
}
MI_block_one_indicator_inter<-function(sepDataOrg,SepIndNamesOrg,SepLVNameOrg,ManRobCov,load_null,intercept_null, intercept_null_for_test, marker_nam, reference_time){
  # To meet assumption of equal LV variance across groups: scale latent variable to have unit variance
  sepDataOrg[,SepLVNameOrg] <- scale(sepDataOrg[,SepLVNameOrg], center = F, scale = sqrt(var(sepDataOrg[,SepLVNameOrg])))
  # Make variable names with time indices
  all_names_list<-TimeIndicesVarNames(SepLVNameOrg=SepLVNameOrg,SepIndNamesOrg=SepIndNamesOrg,sepDataOrg=sepDataOrg)
  # Original indicator and LV names
  OrgNames<-colnames(sepDataOrg)[-c(1,2)]
  # Number of observations and time periods
  N_each <- split(sepDataOrg[,"id"],sepDataOrg[,"time"])
  N<-sum(unlist(lapply(N_each, function(x) {length(x)})))
  time<-length(unique(sepDataOrg[,"time"]))
  # Mean center for loadings null hypothesis
  sepDataOrg1 <- sepDataOrg
  sepDataOrg1_tmp <- vector(mode = "list", length = time)
  
  for (t in 1:time){
    sepDataOrg1_tmp[[t]] <- as.data.frame(lapply(sepDataOrg1[sepDataOrg1$time == t,-c(1,2)], function(x) {scale(x, center = T, scale = FALSE)} ))
  }
  sepDataOrg1[-c(1,2)] <- do.call(rbind,sepDataOrg1_tmp)
  ################## Test Loadings
  # Create matrix with Y_star in diagonal
  # Split LV by group/time
  split_lv_by_time <- lapply(split(sepDataOrg1[,SepLVNameOrg], sepDataOrg1[,"time"]), function(x){x <- as.matrix(x)})
  # Create block diagonal matrix with LV scores on diagonal
  Y_star <- bdiag(split_lv_by_time)
  # Put all LV scores in first column
  Y_star[,reference_time] <- sepDataOrg1[,SepLVNameOrg]
  # Make block diagonal matrix with as many blocks as indicators
  var_dum_list <- vector(mode = "list", length = length(SepIndNamesOrg))
  Y_full_matrix <- bdiag(lapply(var_dum_list, function(x) {x <- Y_star}))
  # Naming the block diagonal matrix columns using the specified loading null hypothesis matrix
  split_null_hyp_names <- split(rownames(load_null), ceiling(seq_along(rownames(load_null))/(time-1)))
  tmp_test_names <- Y_names_from_hypothesis(split_null_hyp_names, reference_time, time)
  colnames(Y_full_matrix) <- tmp_test_names
  # Stack indicators
  stacked_indicators <- stack(sepDataOrg1[,SepIndNamesOrg])[,1]
  # Regression
  reg_load <- lm(stacked_indicators ~ 0 + . ,data = as.data.frame(as.matrix(Y_full_matrix)))
  # Testing loadings (metric invariance)
  load_null_1 <- paste(rownames(load_null), load_null[,1], sep = "=")
  Hyp_load<-linearHypothesis(reg_load, hypothesis.matrix = load_null_1,test=c("Chisq"))
  p_val_load<-Hyp_load$`Pr(>Chisq)`[2]
  Chisq_load<-Hyp_load$Chisq[2]
  
  ################## Test intercepts
  # Getting marker variable
  marker_nam_load <- lm(sepDataOrg[,marker_nam] ~ sepDataOrg[,SepLVNameOrg])$coefficients[2]
  marker_var <- sepDataOrg[,c("time",marker_nam)]
  # fill data matrix for lambda
  for_anova <- data.frame(dependent = rep(0,N*(length(SepIndNamesOrg)-1)), ind_name = "empty", stringsAsFactors = F)
  index_begin <- 1
  index_end <- N
  non_mark_nam <- SepIndNamesOrg[SepIndNamesOrg!=marker_nam]
  for (ind_nam in non_mark_nam) {
    # Loading restricted to be equal over groups
    tmp_load <- lm(sepDataOrg[,ind_nam] ~ sepDataOrg[,SepLVNameOrg])$coefficients[2]
    # Calculating lambda
    for_anova[index_begin:index_end,"dependent"] <- sepDataOrg[,ind_nam] - (tmp_load/marker_nam_load) * sepDataOrg[,marker_nam]
    for_anova[index_begin:index_end,"ind_name"] <- rep(ind_nam,N)
    index_begin <- index_begin + N
    index_end <- index_end + N
  }
  # Dummy variable matrix
  time_dum_list <- vector(mode = "list", length =time)
  D <- bdiag(lapply(N_each, function(x) as.matrix(rep(1,length(x)))))
  # Reference group
  D[,reference_time] <- reference_time
  var_dum_list <- vector(mode = "list", length = length(SepIndNamesOrg)-1)
  D_full_mat <- bdiag(lapply(var_dum_list, function(x) {x <- D}))
  
  split_null_hyp_names <- split(rownames(intercept_null), ceiling(seq_along(rownames(intercept_null))/(time-1)))
  # aux variables is variables for which the coefficients are not reported (only the differences is reported)
  tmp_test_names <- Y_names_from_hypothesis(split_null_hyp_names, reference_time, time)
  colnames(D_full_mat) <- tmp_test_names
  # Intercept testing regression
  intercept_reg <- lm(for_anova$dependent ~ 0 + ., data= as.data.frame(as.matrix(D_full_mat)))
  # Test hypothesis
  intercept_null_linear_hyp <- paste(rownames(intercept_null_for_test), intercept_null_for_test[,1], sep = "=")
  Hyp_inter<-linearHypothesis(intercept_reg, hypothesis.matrix = intercept_null_linear_hyp,test=c("Chisq"))
  p_val_inter<-Hyp_inter$`Pr(>Chisq)`[2]
  Chisq_inter<-Hyp_inter$Chisq[2]
  
  # Output
  out<-list(Intercepts=intercept_reg$coefficients[rownames(intercept_null)],
            TestIntercepts_Chisq=Chisq_inter,
            TestIntercepts_p=p_val_inter,
            Loadings=reg_load$coefficients[rownames(load_null)],
            TestLoadings_Chisq=Chisq_load,
            TestLoadings_p=p_val_load)
  return(out)
}
TimeIndicesVarNames <- function(SepLVNameOrg=SepLVNameOrg,SepIndNamesOrg=SepIndNamesOrg,sepDataOrg=sepDataOrg){
  LVName_t<-paste(SepLVNameOrg,unique(sepDataOrg$time),sep="_")
  IndNames_t<-lapply(SepIndNamesOrg,function(x) {paste(x,unique(sepDataOrg$time),sep="_")})
  all_names<-sort(c(unlist(IndNames_t),LVName_t))
  return(list(all_names=all_names,IndNames_t=IndNames_t,LVName_t=LVName_t))
}
Y_names_from_hypothesis <- function(split_null_hyp_names, reference_time, time){
  if (reference_time == 1) {
    out_nam <- unlist(mapply(split_null_hyp_names, names(split_null_hyp_names), 
                             FUN = function(x,y){x <- c(paste("aux_k",y,sep=""),x)},
                             SIMPLIFY = FALSE))
  }
  if (reference_time == time) {
    out_nam <- unlist(mapply(split_null_hyp_names, names(split_null_hyp_names), 
                             FUN = function(x,y){x <- c(x, paste("aux_k",y,sep=""))},
                             SIMPLIFY = FALSE))
  }
  if ((reference_time != 1) & (reference_time != time)) {
    out_nam <- unlist(mapply(split_null_hyp_names, names(split_null_hyp_names), 
                             FUN = function(x,y){x <- c(x[1:(reference_time-1)], 
                                                        paste("aux_k",y,sep=""),
                                                        x[(reference_time):(length(x))])
                             },
                             SIMPLIFY = FALSE))
  }
  return(out_nam)
  
}
# Get p-values from bootstrapping
boot_p_val<-function(boot_Chisq, B, Org_res){
  # Matrix to hold results
  p_val <- matrix(0, ncol = 2)
  colnames(p_val) <-  c("p_val_loadings", "p_val_intercept")
  # Sort bootstrapped chisq values
  sort_inter <- sort(boot_Chisq[,"chisq_intercept"])
  names(sort_inter)<-1:B
  # Number of bootstrapped chi_sq values above original chi_square
  sort_inter<-Org_res$TestIntercepts_Chisq<=sort_inter
  true_sort_inter <- !(sort_inter == TRUE)
  # If no bootstrapped chi_sq above --> p-val = 0
  if (sum(true_sort_inter)==B) {
    p_val[,"p_val_intercept"]<- 0
  } else { # p-val calculated as wooldridge (panel data and cross-section)
    min_num_inter<-min(which(sort_inter == TRUE))
    p_val[,"p_val_intercept"]<- 1-min_num_inter/(B+1)
  }
  # Sort bootstrapped chisq values
  sort_load <- sort(boot_Chisq[,"chisq_load"])
  names(sort_load)<-1:B
  # Number of bootstrapped chi_sq values above original chi_square
  sort_load<-Org_res$TestLoadings_Chisq<=sort_load
  true_sort_load <- !(sort_load == TRUE)
  # If no bootstrapped chi_sq above --> p-val = 0
  if (sum(true_sort_load)==B) {
    p_val[,"p_val_loadings"] <- 0
  } else {# p-val calculated as wooldridge (panel data and cross-section)
    min_num_load<-min(which(sort_load == TRUE))
    p_val[,"p_val_loadings"]<- 1-min_num_load/(B+1)
  }
  return(p_val)
  
}



