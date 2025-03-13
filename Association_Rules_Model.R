####LOADING THE 'arules' PACKAGE - VERSION 1.5-2####
# Define the arules package ver 1.5-2 url.
url <- "https://cran.r-project.org/src/contrib/Archive/arules/arules_1.5-2.tar.gz"
# Install the arules package ver 1.5-2.
install.packages(url, repos=NULL, type="source")
# load into environment
remotes::install_version("arules", version = "1.5-2", repos = "http://cran.us.r-project.org")

####LOAD THE PACKAGES####
library(readr)
library(dplyr)
library(tidyverse)
library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(arules)
library(reshape2)
library(epitools)
library(coxme)
library(forcats)
library(ggtext)

####SET SEED####
set.seed(7895)

####LOAD THE DATA####

# Ensure the analysed data is tidy, with categorical factor variables. 
# Presence of a comorbidity should be encoded as '1', and absence as '0'.
# Split the data into 'case' and 'control' cohorts, based on the outcome variable.

####DEFINE ASSOCIATION RULES FUNCTIONS####
get_most_frequent_patterns <- function(data_in) {
  cnames <- names(data_in)
  col_inds <- -1 * c(1:6)
  
  data_matrix_1 <- as.matrix(data_in[,col_inds])
  data_transactions_1 <- as(data_matrix_1, "transactions")
  
  most_frequent_patterns <-  eclat(data_transactions_1, 
                                   parameter=new("ECparameter", 
                                                 tidLists=TRUE, 
                                                 support = 0.005, 
                                                 ext = TRUE, 
                                                 minlen=1, 
                                                 maxlen=ncol(data_in)), 
                                   control= new("ECcontrol", verbose = FALSE)
  )
  return(most_frequent_patterns)
}

select_and_multiply <- function(conditional_frequent_sets, TF_interestingness, values) {
  interesting_sets <- subset(conditional_frequent_sets, TF_interestingness)
  out <- matrix(nrow=0, ncol=0)
  if (length(interesting_sets) > 0){
    items_in_interesting_sets <- as(items(x=interesting_sets), "matrix")
    out <- diag(values[TF_interestingness]) %*% items_in_interesting_sets
  }
  return(out)
}

combination_matrix <- function(most_frequent_patterns) {	
  matrix_N_Cases <- select_and_multiply(conditional_frequent_sets = most_frequent_patterns, 
                                        TF_interestingness = rep(TRUE, length(most_frequent_patterns)),
                                        values = quality(most_frequent_patterns)$transIdenticalToItemsets)
  
  combn_mtrx <- as.data.frame(matrix_N_Cases)
  combn_mtrx[combn_mtrx > 0] <- 1
  
  return(combn_mtrx)
}

get_combo_names <- function(most_frequent_patterns)  {
  combn_mtrx <- combination_matrix(most_frequent_patterns)
  
  str_combn <- function(x){
    acr_ids <- colnames(combn_mtrx)[x == 1]
    str <- paste(acr_ids, collapse = "+")
    return(str)
  }
  
  combs <- apply(combn_mtrx, 1, FUN = str_combn)
  return(combs)
}

get_comb_freq <- function(most_frequent_patterns) {	
  matrix_N_Cases <- select_and_multiply(conditional_frequent_sets = most_frequent_patterns, 
                                        TF_interestingness = rep(TRUE, length(most_frequent_patterns)),
                                        values = quality(most_frequent_patterns)$transIdenticalToItemsets)
  freqs <- apply(matrix_N_Cases, 1, max)
  return(freqs)
}

interestingness_scores <- function(most_frequent_patterns, data_in) {
  combn_mtrx <- combination_matrix(most_frequent_patterns)
  
  cnames <- names(data_in)
  col_inds <- -1 * c(1:6)
  n1 <- dim(data_in)[1]
  prevalence <- colSums(data_in[,col_inds] / n1)
  
  prob_col <- function(i) {
    v <- as.numeric(combn_mtrx[,i])
    pv <- rep(0, length(v))
    pv[v == 0] <- 1
    pv[v > 0] <- as.numeric(prevalence[i])
    return(pv)
  }
  
  supp <- quality(most_frequent_patterns)$support
  prob_mtrx <- do.call(cbind, lapply(1:length(prevalence), prob_col))
  null_prob <- apply(prob_mtrx, 1, FUN = function(x){exp(sum(log(x)))})
  lift <- supp / null_prob
  
  return(lift)
}

####ASSOCIATION RULES ANALYSIS####

# Define the interestingness score cut off.
iscore_cutoff <- 2.0

# Case.
most_frequent_patterns_case <- get_most_frequent_patterns(df_case)
combn_mtrx_case <- combination_matrix(most_frequent_patterns_case)
freqs_case <- get_comb_freq(most_frequent_patterns_case)
combn_name_case <- get_combo_names(most_frequent_patterns_case)
iscore_case <- interestingness_scores(most_frequent_patterns_case, df_case)
df_result_case <- data.frame(ID=1:(length(freqs_case)),
                             COMBN=as.character(combn_name_case),
                             ISCORE=iscore_case,
                             LOG_ISCORE=log(iscore_case),
                             FREQ=freqs_case)
df_result_case <- df_result_case[order(df_result_case$ISCORE),]
df_result_case_select <- df_result_case[df_result_case$ISCORE >= iscore_cutoff,]
nrow_case <- dim(df_result_case_select)[1]
df_result_case_select_long <- data.frame(COMBN=rep(df_result_case_select$COMBN, 2), DATA=c(df_result_case_select$LOG_ISCORE, df_result_case_select$FREQ), TYPE=rep(c("Log(O/E)", "Number of occurrence"), each=nrow_case))
df_result_case_select_long$COMBN <- factor(df_result_case_select_long$COMBN, levels=df_result_case_select$COMBN)

# Control.
most_frequent_patterns_control <- get_most_frequent_patterns(df_control)
combn_mtrx_control <- combination_matrix(most_frequent_patterns_control)
freqs_control <- get_comb_freq(most_frequent_patterns_control)
combn_name_control <- get_combo_names(most_frequent_patterns_control)
iscore_control <- interestingness_scores(most_frequent_patterns_control, df_control)
df_result_control <- data.frame(ID=1:(length(freqs_control)),
                                COMBN=as.character(combn_name_control),
                                ISCORE=iscore_control,
                                LOG_ISCORE=log(iscore_control),
                                FREQ=freqs_control)
df_result_control <- df_result_control[order(df_result_control$ISCORE),]
df_result_control_select <- df_result_control[df_result_control$ISCORE >= iscore_cutoff,]
nrow_control <- dim(df_result_control_select)[1]
df_result_control_select_long <- data.frame(COMBN=rep(df_result_control_select$COMBN, 2), DATA=c(df_result_control_select$LOG_ISCORE, df_result_control_select$FREQ), TYPE=rep(c("Log(O/E)", "Number of occurrence"), each=nrow_control))
df_result_control_select_long$COMBN <- factor(df_result_control_select_long$COMBN, levels=df_result_control_select$COMBN)