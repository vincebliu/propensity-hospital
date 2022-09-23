library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(forcats)
library(MatchIt)
library(nnet)
library(rpart)
library(glmnet)
library(caret)
library(purrr)
library(rsample)
library(recipes)
library(pROC)


theme_set(theme_bw(base_size = 20))

dat<-read_csv("sample50pct.csv") # read the Data

########################################################################################
# make the PAC column more in line with the paper
########################################################################################
clean_dat=rename(dat, PAC=dischargegrp)
# remove AMA and "Others", which are both infrequent categories
clean_dat=clean_dat %>%
  mutate(PAC=case_when(PAC=="5.Discharged to home/ self-care"~"Home",
                       PAC=="2.Discharge to SNF" ~ "SNF",
                       PAC=="4.Discharged to Rehab" ~ "IRF",
                       PAC=="1.DISCH/TRANS TO HOME CARE OR HOME HEALTH" ~ "HHA",
                       PAC=="3. Discharged to LTC, federal hospital, p" ~ "LTCH",
                       TRUE ~ "Remove")
  )
clean_dat=clean_dat %>% filter(PAC!="Remove")
clean_dat=clean_dat %>% relocate(PAC)

#Relocate the response variable to the second column
clean_dat=clean_dat %>% relocate(readmit30, .after=PAC)

clean_dat=clean_dat[complete.cases(clean_dat),]

# split into train/test

set.seed(2012021)
train_test_split <- rsample::initial_split(clean_dat, prop = 3/4, strata = PAC)
train_dat <- rsample::training(train_test_split)
test_dat <- rsample::testing(train_test_split)

########################################################################################
# Helper Functions
########################################################################################

preprocess=function(pac, df) {
  df=filter(df, PAC %in% c("Home", pac))
  trained_rec=recipe(readmit30 ~ ., data=train_dat%>%select(-PAC)) %>%
    step_nzv(all_predictors(), freq_cut=50/5, unique_cut=50) %>%
    step_dummy(all_predictors(), -all_numeric()) %>%
    step_corr(all_predictors(), threshold=.8) %>%
    step_nzv(all_predictors(), freq_cut=50/5, unique_cut=50) %>%
    step_center(all_predictors()) %>%
    step_scale(all_predictors()) %>%
    prep(training=train_dat)
  
  dat=bake(trained_rec, new_data=df) %>% mutate(PAC=df$PAC)
  return(dat)
}



#' Do propensity score matching for a particular PAC.
#' Estimate propensity scores
#' Binary model comparing "Home care" against different PACs.
#' @param df is a dataframe
#' @param pac is a string, the PAC category (e.g. HHA)
#' @param dist is a string, the method used to estimate propensity scores (comes from matchit, e.g. nnet, glm, randomforest)
#' @return a list of 1. matched df, 2. propensity score vector, 3. fitted model
est_psm <- function(pac, df, dist){
  message(str_glue('PAC is {pac}')) 
  df_pac <- df %>%
    filter(PAC %in% c("Home", pac))
  df_pac=preprocess(pac, df)
  df_pac=df_pac %>%
    mutate(is_pac = ifelse(PAC=="Home", 0, 1)) %>%
    relocate(is_pac, .after=PAC)
  
  if (dist=="nnet") {
    # rename the levels to be "yes/no" instead of 1/0
    df_pac$is_pac <- as.factor(df_pac$is_pac) %>%
      fct_recode(yes = '1', no = '0')
    
    
    p_mod <- caret::train(is_pac ~ . - PAC - readmit30, df_pac, 
                          method='nnet', metric = 'ROC',
                          #MaxNWts=3700,
                          trControl = trainControl(method = 'cv', number = 5, classProbs = T,
                                                   summaryFunction = twoClassSummary),
                          tuneGrid=expand.grid(.size=c(1,5),
                                               .decay=c(0.001,0.1))
    ) 
    prop_scores <- as.numeric(fitted(p_mod))
    
    matches_pac <- matchit(is_pac ~ . - PAC - readmit30, 
                           method="nearest", data=df_pac, 
                           distance = prop_scores)
    
  } else if(dist == 'lasso'){
    # preprocessing x for form of glmnet
    lasso_x <- df_pac %>%
      select(-c('PAC', 'readmit30', 'is_pac')) %>% # get rid of response cols
      #mutate(across(where(is.numeric), ~ scale(.x, center = T, scale = T))) %>% # scale
      as.matrix.data.frame() %>%
      scale(center = T, scale = T)
    
    
    
    # fit lasso model
    p_mod <- glmnet::cv.glmnet(x = lasso_x, y = df_pac$is_pac, 
                               nfolds = 5, family = 'binomial', alpha = 1)
    
    # get estiamted propensity scores
    prop_scores <- predict(p_mod, newx = lasso_x, 
                           s = 'lambda.min', type = 'response')
    
    # use propensity scores to get a match
    matches_pac <- matchit(is_pac ~ . - PAC - readmit30, 
                           method="nearest", data=df_pac, 
                           distance=as.numeric(prop_scores))
    
  } else if(dist == 'rf'){
    # rename the levels to be "yes/no" instead of 1/0
    df_pac$is_pac <- as.factor(df_pac$is_pac) %>%
      fct_recode(yes = '1', no = '0')
    
    p_mod <- caret::train(is_pac ~ . - PAC - readmit30, df_pac, 
                          method='ranger', metric = 'ROC',
                          num.trees = 100,
                          trControl = trainControl(method = 'cv', number = 5, classProbs = T,
                                                   summaryFunction = twoClassSummary),
                          tuneGrid= expand.grid(.mtry = c(3,4),
                                                .splitrule = 'gini',
                                                .min.node.size = 10)
    ) 
    prop_scores <- predict(p_mod, type = 'prob')$yes
    
    matches_pac <- matchit(is_pac ~ . - PAC - readmit30, 
                           method="nearest", data=df_pac, 
                           distance = prop_scores)
  } else if(dist == 'glm'){
    matches_pac <- matchit(is_pac ~ . - PAC - readmit30, 
                           method="nearest", data=df_pac, 
                           distance = 'glm', distance.options = list(maxit = 100))
    p_mod <- matches_pac$model
  } else{
    stop('invalid distance. should be one of "rf", "nn", "glm", "lasso"')
  }
  df_pac_matched <- match.data(matches_pac)
  
  # return matched data, as well as the scores themselves
  list(
    match_df = df_pac_matched,
    pscore = matches_pac$distance,
    pmod = p_mod
  )
}


#' scale column (with name) of df_to_scale by subtracting the mean and dividng sd of same col in df_for_ref 
scale_data <- function(col_to_scale, df_for_ref, name){
  
  metadata_cols <- c("PAC", 'readmit30', 'is_pac')
  # if it's metadata, return the original column without any changes.
  if(name %in% metadata_cols){
    return(col_to_scale)
  }
  
  col_for_ref <- df_for_ref %>% pull(name)
  
  # we need at least 2 observations in the ref col to get a standard deviation
  # so if the whole column has < 2 observations in df_for_ref, then return the original column without any changes
  if(sum(!is.na(col_for_ref)) < 2){
    return(col_to_scale)
  }
  
  (col_to_scale - mean(col_for_ref)) / sd(col_for_ref)
}


#function to return fpr, tpr given prediction and true label
#label ordering goes c(negative class, positive class)
#' @param pred are predictions
#' @param label is the true values
#' @param y is the metric on y axis of desired plot (fpr for roc)
#' @param x is the metric on x axis of desired plot (fpr for roc)
#' Can also make precision recall curve by letting x = rec (or tpr), y = prec (or ppv)
fpr_tpr <- function(pred, label, y = 'tpr', x = 'fpr', label_ordering = NULL){
  rocpred <- ROCR::prediction(pred, label, label.ordering = label_ordering)
  rocfpr_tpr <- ROCR::performance(rocpred, measure = y, x.measure = x)
  
  # if plotting roc, get auc. if plotting precision-recall, get roprc
  if(y == 'tpr'){
    rocauc <- ROCR::performance(rocpred, measure = 'auc')
  } else{
    rocauc <- ROCR::performance(rocpred, measure = 'aucpr')
  }
  
  return(tibble(x = deframe(rocfpr_tpr@x.values), 
                y = deframe(rocfpr_tpr@y.values),
                auc = as.numeric(rocauc@y.values)))
}

#' get predictions on test set
#' @param outcome_mod separates whether we're doing PS models or outcomes models
#' @pac is a string, one of the four PAC abbreivations
#' @dist is the method used to calculate PS.
#' @param psm_list_item is one element of the output of est_psm() if outcome_mod = F, 
#'                      and one element of the output of outcome_mod() if outcome_mod = T
#' @output a list of 2 tables: one with the predictions, and one with fpr/tpr values for plotting ROC
i=0
test_set_pred <- function(psm_list_item, pac, dist, outcome_mod = F){
  message(pac)
  test_pac <- test_dat %>%
    filter(PAC %in% c("Home", pac))
  test_pac=preprocess(pac, test_pac)
  test_pac=test_pac %>%
    mutate(is_pac = ifelse(PAC=="Home", 0, 1)) %>%
    relocate(is_pac, .after=PAC)
  test_pac=test_pac[complete.cases(test_pac),]
  
  psm_truth <- test_pac$is_pac
  
  # if this is for the outcome models, only 1 method (logit glm), 
  #so no need to go through all these shenenanigans
  if(outcome_mod){
    truth=test_pac$readmit30
    psm_pred <- predict(psm_list_item, test_pac, type = 'response')
    
    i=which(is.na(psm_pred))
    if (length(i)>0)
      return(data.frame(pred=psm_pred[-i], truth=truth[-i]))
    return(data.frame(pred=psm_pred, truth=truth))
  }
  
  # requires scaled data, special preprocessing
  if(dist == 'lasso'){
    # get training data for reference in scaling
    train_pac <- train_dat %>%
      filter(PAC %in% c("Home", pac))
    train_pac=preprocess(pac, train_pac)
    train_pac=train_pac%>%
      mutate(is_pac = ifelse(PAC=="Home", 0, 1)) %>%
      relocate(is_pac, .after=PAC)
    
    # preprocessing x for form of glmnet
    # scale using train column mean/max
    lasso_test <- test_pac %>%
      select(-c('PAC', 'readmit30', 'is_pac')) %>% # get rid of response cols
      purrr::map2_dfc(names(.), ~scale_data(.x, train_pac, .y)) %>%
      as.matrix.data.frame() 
    
    psm_pred <- predict(psm_list_item$pmod, lasso_test, s = 'lambda.min', type= 'response')
  } else if(dist == 'glm'){
    psm_pred <- predict(psm_list_item$pmod, test_pac, type = 'response')
  } else{
    psm_pred <- predict(psm_list_item$pmod, test_pac, type = 'prob')$yes
  }
  
  #pred_table <- tibble(pac = pac,
  #                     method = dist,
  #                     truth = psm_truth,
  #                     pred = psm_pred
  #)
  i=which(is.na(psm_pred))
  if (length(i)>0)
    return(data.frame(pred=psm_pred[-i], truth=psm_truth[-i]))
  return(data.frame(pred=psm_pred, truth=psm_truth))
  
  #roc_table <- fpr_tpr(pred_table$pred[-i], pred_table$truth[-i]) %>%
  #  mutate(pac = pac, method = dist)
  
  #list(pred_df = pred_table,
  #     roc_df = roc_table)
}

# preproc <- recipe(readmit30 ~ . - PAC - distance - weights - subclass, data=df_pac_matched) %>%
#   step_normalize(all_numeric(), -all_outcomes())

########################################################################################
# Part 1: Trying different methods of estimating propensity scores
########################################################################################

pac_names <- c('HHA', 'SNF', 'IRF', 'LTCH')
method_names <- c('glm', 'rf', 'nnet', 'lasso', "no_psm")

# LOGIT PSM (done in the original paper)----- 
time=Sys.time(); logit_psm_list <- purrr::map(pac_names, ~est_psm(.x, train_dat, 'glm')); print(Sys.time()-time)
#Time difference of 5.057268 mins

# random forest PSM
time=Sys.time(); rf_psm_list <- purrr::map(pac_names, ~est_psm(.x, train_dat, 'rf')); print(Sys.time()-time)


# one layer nn PSM
time=Sys.time(); nn_psm_list <- purrr::map(pac_names, ~est_psm(.x, train_dat, 'nnet')); print(Sys.time()-time)
#Time difference of 2.741376 mins

# LASSO psm
time=Sys.time(); lasso_psm_list <- purrr::map(pac_names, ~est_psm(.x, train_dat, 'lasso')); print(Sys.time()-time)
#Time difference of  mins

unmatch_list=purrr::map(pac_names, function(x) {
  df=filter(train_dat, PAC %in% c("Home", x))
  df=preprocess(x, train_dat)
  return (  list(
    match_df = df,
    pscore = NULL,
    pmod = NULL
  ))
})

## the above code takes a long time to run. so save output
saveRDS(logit_psm_list, file = here::here('STAT201', '201B_project', 'logit_psm_list.Rds'))
saveRDS(rf_psm_list, file = here::here('STAT201', '201B_project', 'rf_psm_list.Rds'))
saveRDS(nn_psm_list, file = here::here('STAT201', '201B_project', 'nn_psm_list.Rds'))
saveRDS(lasso_psm_list, file = here::here('STAT201', '201B_project', 'lasso_psm_list.Rds'))

## (how to read in the saved output (see email for the files))
logit_psm_list <- readRDS(file = here::here('STAT201', '201B_project', 'logit_psm_list.Rds'))
rf_psm_list <- readRDS(file = here::here('STAT201', '201B_project', 'rf_psm_list.Rds'))
nn_psm_list <- readRDS(file = here::here('STAT201', '201B_project', 'nn_psm_list.Rds'))
lasso_psm_list <- readRDS(file = here::here('STAT201', '201B_project', 'lasso_psm_list.Rds'))


# get predictions + roc tables
logit_test_out <- imap(logit_psm_list, ~test_set_pred(.x, pac_names[.y], 'glm'))
rf_test_out <- imap(rf_psm_list, ~test_set_pred(.x, pac_names[.y], 'rf'))
nn_test_out <- imap(nn_psm_list, ~test_set_pred(.x, pac_names[.y], 'nnet'))
lasso_test_out <- imap(lasso_psm_list, ~test_set_pred(.x, pac_names[.y], 'lasso'))

test_roc_dfs=data.frame()
for (i in 1:4) {
  lasso_test_out[[i]]=rename(lasso_test_out[[i]], pred=X1)
  a=roc(predictor=logit_test_out[[i]]$pred, response=logit_test_out[[i]]$truth, positive=1)
  b=roc(predictor=rf_test_out[[i]]$pred, response=rf_test_out[[i]]$truth, positive=1)
  c=roc(predictor=nn_test_out[[i]]$pred, response=nn_test_out[[i]]$truth, positive=1)
  d=roc(predictor=lasso_test_out[[i]]$pred, response=lasso_test_out[[i]]$truth, positive=1)
  test_roc_dfs=bind_rows(test_roc_dfs,
    data.frame(sensitivity=a$sensitivities, specificity=a$specificities, pac=pac_names[i], method="glm"),
    data.frame(sensitivity=b$sensitivities, specificity=b$specificities, pac=pac_names[i], method="rf"),
    data.frame(sensitivity=c$sensitivities, specificity=c$specificities, pac=pac_names[i], method="nnet"),
    data.frame(sensitivity=d$sensitivities, specificity=d$specificities, pac=pac_names[i], method="lasso")
  )
}



# pull out just the roc tables
test_roc_dfs <- list(logit_test_out, rf_test_out, 
                     nn_test_out, lasso_test_out) %>%
  map_depth(2, ~.x[[2]]) %>%
  bind_rows()

# get aucs
oos_auc_df <- test_roc_dfs %>%
  group_by(pac, method) %>%
  summarize(auc = auc[1], .groups = 'drop')



# plot roc curves
ggplot(test_roc_dfs %>% mutate(y=specificity, x=1-sensitivity)) +
  geom_line(aes(x,y, color = method), size = 1) + 
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~pac) +
  labs(
    title = 'ROC curves for predicting treatment against being sent home',
    x = 'False Positive Rate',
    y = 'True Positive Rate'
  ) #+
  #theme(legend.text=element_text(size=rel(2)))

ggsave(here::here('STAT201', '201B_project', 'ps_roc.jpg'), height = 10, width = 14)

# # 
# orhha=exp(modelhha$coef["hha"])
# orsnf=exp(modelsnf$coef["snf"])
# orirf=exp(modelirf$coef["irf"])
# orltch=exp(modelltch$coef["ltch"])
# 
# cat("      Replicated FromPaper\n",
#     "HHA:  ", orhha, " 1.26\n",
#     "SNF:  ", orsnf, " 1.25\n",
#     "IRF:  ", orirf, " .77\n",
#     "LTCH: ", orltch, " .76\n") 
# 


########################################################################################
# Part 2: Using propensity score MATCHING in an outcome model
########################################################################################

#' function to fit a logistic glm after computing propensity scores
#' @param psm_list_item is a list, one element of the output of est_psm()
#' @param pac a string, a PAC to run the model on (compared to home)
#' @param use_weights a binary, 
#'                    if TRUE, will use FULL training data to fit model, WITH inv prop weights
#'                    if FALSE, will use MATCHED training data to fit model, WITHOUT weights
outcome_mod <- function(psm_list_item, pac, use_weights = F){
  message(str_glue('PAC is {pac}')) 
  
  if(use_weights){
    # use full training data
    df_pac <- train_dat %>%
      filter(PAC %in% c("Home", pac)) %>%
      mutate(PAC = as.factor(PAC) %>%
               fct_relevel('Home', after = 0))
    
    pscores <- psm_list_item$pscore
    pweights <- with(df_pac, ifelse(PAC == 'Home', 1/(1-pscores), 1/pscores))
    
    
    mod_pac <- glm(readmit30 ~ ., data=df_pac, family=quasibinomial(link="logit"), 
                   weights = pweights)
  } else{
    # use matched data
    df_pac <- psm_list_item$match_df %>%
      select(-c('distance', 'weights', 'subclass')) %>%
      mutate(is_pac=ifelse(PAC==pac, 1, 0)) %>%
      select(-PAC)
    
    mod_pac <- glm(readmit30 ~ ., data=df_pac, family=quasibinomial(link="logit"))
  }
  
  mod_pac
}


##################################
# Outcome models using matched data
### NATHAN NOTE: GO BACK AND FIX CONVERGENCE ISSUES AFTER PRESENTATION
##################################

# logit is original paper
logit_match_mods <- imap(logit_psm_list, ~outcome_mod(.x, pac_names[.y], use_weights = F))
rf_match_mods <- imap(rf_psm_list, ~outcome_mod(., pac_names[.y], use_weights = F))
nn_match_mods <- imap(nn_psm_list, ~outcome_mod(.x, pac_names[.y], use_weights = F))
lasso_match_mods <- imap(lasso_psm_list, ~outcome_mod(., pac_names[.y], use_weights =F))
unmatch_mods <- purrr::map(unmatch_list, function(x) {
  x=mutate(x[[1]], is_pac=ifelse(PAC=="Home", 0, 1)) %>%
    select(-PAC)
  glm(readmit30 ~ ., data=x, family=quasibinomial(link="logit"))
})

## the above code takes a long time to run. so save output
saveRDS(logit_match_mods, file = here::here('STAT201', '201B_project', 'logit_match_mods.Rds'))
saveRDS(rf_match_mods, file = here::here('STAT201', '201B_project', 'rf_match_mods.Rds'))
saveRDS(nn_match_mods, file = here::here('STAT201', '201B_project', 'nn_match_mods.Rds'))
saveRDS(lasso_match_mods, file = here::here('STAT201', '201B_project', 'lasso_match_mods.Rds'))


# get tables of odds ratios
match_or_df <- list(logit_match_mods, rf_match_mods,
                       nn_match_mods, lasso_match_mods, unmatch_mods) %>%
  map_dfr(function(y) 
    exp(map_dbl(y, ~coef(.x)["is_pac"])) %>%
      set_names(nm = pac_names)
    ) %>%
  bind_cols(method = method_names) %>%
  relocate(method, .before = 1)

saveRDS(match_or_df, here::here('STAT201', '201B_project', 'match_or_df.Rds'))

# get predictions + roc tables
logit_test_match <- imap(logit_match_mods, ~test_set_pred(.x, pac_names[.y], 'glm', outcome_mod = T))
rf_test_match <- imap(rf_match_mods, ~test_set_pred(.x, pac_names[.y], 'rf', outcome_mod = T))
nn_test_match <- imap(nn_match_mods, ~test_set_pred(.x, pac_names[.y], 'nnet', outcome_mod = T))
lasso_test_match <- imap(lasso_match_mods, ~test_set_pred(.x, pac_names[.y], 'lasso', outcome_mod = T))
unmatch_test_match <- imap(unmatch_match_mods, ~test_set_pred(.x, pac_names[.y], 'no_psm', outcome_mod = T))

test_match_roc_dfs=data.frame()
for (i in 1:4) {
  a=roc(predictor=logit_test_match[[i]]$pred, response=logit_test_match[[i]]$truth, positive=1)
  b=roc(predictor=rf_test_match[[i]]$pred, response=rf_test_match[[i]]$truth, positive=1)
  c=roc(predictor=nn_test_match[[i]]$pred, response=nn_test_match[[i]]$truth, positive=1)
  d=roc(predictor=lasso_test_match[[i]]$pred, response=lasso_test_match[[i]]$truth, positive=1)
  e=roc(predictor=unmatch_test_match[[i]]$pred, response=unmatch_test_match[[i]]$truth, positive=1)
  test_match_roc_dfs=bind_rows(test_match_roc_dfs,
                         data.frame(sensitivity=a$sensitivities, specificity=a$specificities, pac=pac_names[i], method="glm"),
                         data.frame(sensitivity=b$sensitivities, specificity=b$specificities, pac=pac_names[i], method="rf"),
                         data.frame(sensitivity=c$sensitivities, specificity=c$specificities, pac=pac_names[i], method="nnet"),
                         data.frame(sensitivity=d$sensitivities, specificity=d$specificities, pac=pac_names[i], method="lasso"),
                         data.frame(sensitivity=e$sensitivities, specificity=e$specificities, pac=pac_names[i], method="no_psm")
  )
}

# pull out just the roc tables
# test_match_roc_dfs <- list(logit_test_match, rf_test_match,
#                      nn_test_match, lasso_test_match) %>%
#   map_depth(2, ~.x[[2]]) %>%
#   bind_rows() #%>%
#   # ## NATHAN FIX THIS LATER
#   # mutate(y = ifelse(pac == 'IRF' & method == 'rf', 1-y, y))
# 
# # get aucs
# oos_auc_match_df <- test_match_roc_dfs %>%
#   group_by(pac, method) %>%
#   summarize(auc = auc[1], .groups = 'drop')
# 
# # plot roc curves
 ggplot(test_match_roc_dfs %>% mutate(y=specificity, x=1-sensitivity)) +
   geom_line(aes(x,y, color = method), size = 1) +
   geom_abline(slope = 1, intercept = 0) +
   facet_wrap(~pac) +
   labs(
     title = 'ROC curves for predicting hospital readmission',
     x = 'False Positive Rate',
     y = 'True Positive Rate'
   ) +
   scale_color_manual(values=c("black",
                              "red",
                              "green",
                              "purple",
                              "yellow"))
 #+
#theme(legend.text=element_text(size=rel(2)))
# 
# ggsave(here::here('STAT201', '201B_project', 'outcome_match_roc.jpg'), height = 10, width = 14)



##################################
# Outcome models using PS weights
##################################

logit_weight_mods <- imap(logit_psm_list, ~outcome_mod(.x, pac_names[.y], use_weights = T))
nn_weight_mods <- imap(nn_psm_list, ~outcome_mod(.x, pac_names[.y], use_weights = T))
rf_weight_mods <- imap(rf_psm_list, ~outcome_mod(., pac_names[.y], use_weights = T))
lasso_weight_mods <- imap(lasso_psm_list, ~outcome_mod(., pac_names[.y], use_weights =T))


## the above code takes a long time to run. so save output
saveRDS(logit_weight_mods, file = here::here('STAT201', '201B_project', 'logit_weight_mods.Rds'))
saveRDS(rf_weight_mods, file = here::here('STAT201', '201B_project', 'rf_weight_mods.Rds'))
saveRDS(nn_weight_mods, file = here::here('STAT201', '201B_project', 'nn_weight_mods.Rds'))
saveRDS(lasso_weight_mods, file = here::here('STAT201', '201B_project', 'lasso_weight_mods.Rds'))


# get tables of odds ratios
weight_or_df <- list(logit_weight_mods, rf_weight_mods,
                       nn_weight_mods, lasso_weight_mods) %>%
  map_dfr(function(y) 
    exp(map_dbl(y, ~coef(.x)[2])) %>%
      set_names(nm = pac_names)
  ) %>%
  bind_cols(method = method_names) %>%
  relocate(method, .before = 1)

saveRDS(weight_or_df, here::here('STAT201', '201B_project', 'weight_or_df.Rds'))


# get predictions + roc tables
logit_test_weight <- imap(logit_weight_mods, ~test_set_pred(.x, pac_names[.y], 'glm', outcome_mod = T))
rf_test_weight <- imap(rf_weight_mods, ~test_set_pred(.x, pac_names[.y], 'rf', outcome_mod = T))
nn_test_weight <- imap(nn_weight_mods, ~test_set_pred(.x, pac_names[.y], 'nnet', outcome_mod = T))
lasso_test_weight <- imap(lasso_weight_mods, ~test_set_pred(.x, pac_names[.y], 'lasso', outcome_mod = T))


# pull out just the roc tables
test_weight_roc_dfs <- list(logit_test_weight, rf_test_weight,
                     nn_test_weight, lasso_test_weight) %>%
  map_depth(2, ~.x[[2]]) %>%
  bind_rows() #%>%
  # ## NATHAN FIX THIS LATER
  # mutate(y = ifelse(pac == 'IRF' & method == 'rf', 1-y, y))

# get aucs
oos_auc_weight_df <- rocs %>%
  group_by(pac, method) %>%
  summarize(auc = auc[1], .groups = 'drop')

rocs=rbind(a1$roc_df, b1$roc_df, c1$roc_df, d1$roc_df, e1$roc_df, f1$roc_df, g1$roc_df)

# plot roc curves
ggplot(rocs) +
  geom_line(aes(x,y, color = method), size = 1) +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~pac) +
  labs(
    title = 'ROC curves for predicting hospital readmission',
    subtitle = 'weighted regression',
    x = 'False Positive Rate',
    y = 'True Positive Rate'
  ) #+
theme(legend.text=element_text(size=rel(2)))
# 
# ggsave(here::here('STAT201', '201B_project', 'outcome_weight_roc.jpg'), height = 10, width = 14)
