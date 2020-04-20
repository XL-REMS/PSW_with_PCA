##########################################
###  Project: Data Analysis for PSM    ###
###  Author: Seungman Kim, Yifang Zeng ###
###  Last Modified:                    ###
##########################################

# 03-26-2020

######################################################################################
##Load pkg
######################################################################################
library(data.table)
library(MatchIt)
library(foreach)
library(dplyr)
######################################################################################
##Set up
######################################################################################
rm(list = ls())

names_XX <- c("XX20_", "XX50_", "XX70_")
names_TS <- c("TS10", "TS25", "TS40")
names_TY <- c("TY00_", "TY20_", "TY70_")
folderName <- c()

for (j in 1:3){ # Start loop #1
  for (k in 1:3){ # Start loop #2
    for (l in 1:3){ # Start loop #3
      
      # Create folder name
      folderName <- c(folderName, paste0(names_XX[j],names_TY[l],names_TS[k]))
      
    } # End loop #3
  } # End loop #2
} # End loop #1

# Find the directory in where this file is located
# cur_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

cur_dir <- getwd()
cur_dir

#####################################################################################
# Global Vars
#####################################################################################
varcut <- 0.4 # Threshold for cumulative variance 
act_list <- list() # Active list
col_names <- c(paste0('X', c(2:3, 5:6, 8:9)), "Y", "T", 
               paste0('X', c(2:3, 5:6, 8:9),'SQ'), 
               paste0('X2X', c(3, 5:6, 8:9)), 
               paste0('X3X', c(5:6, 8:9)), 
               paste0('X5X', c(6, 8:9)), 
               paste0('X6X', 8:9),
               'X8X9')

formula <- paste("T", paste(col_names[-c(7,8)], collapse=" + "), sep = " ~ ")

######################################################################################
# Function to fill active_list
######################################################################################

fill_active_list <- function(){
  # fread in each csv file
  act_list[[i]] <<- fread(file_list[i], header=TRUE)
  
  # combind 3 parts
  act_list[[i]] <<- cbind(act_list[[i]][,c(2:3, 5:6, 8:11)], # 1st part: x2, X3, X5, X6, X8, x9, Y, T
                          act_list[[i]][,c(2:3, 5:6, 8:9)]^2,  # 2nd part: all the square term of x2sq, X3sq, X5sq, X6sq, X8sq, x9sq
                          model.matrix( ~.^2, data=act_list[[i]][,c(2:3, 5:6, 8:9)])[,c(8:22)] # 3rd part: all the possible two-way interaction term among x2, X3, X5, X6, X8, x9
  )
  # assigning proper column names for combined matrix
  colnames(act_list[[i]]) <<- col_names
}


#####################################################################################
# Main for loop
#####################################################################################
# ptm <- proc.time()
# wrapitr <- 1 
for (wrapitr in 1:length(folderName)){
  
  # Set the working directory
  
  # Test set
  setwd(paste0("C:/Users/Johnny/Dropbox/Research Projects/Current/PSM Group/P2_PSW_PCA/",folderName[wrapitr]))
  
  # Local
  # setwd(paste0("D:/PSM_PC/DataSet/",folderName[wrapitr]))
  
  # Quanah
  # setwd(paste0("/home/seungmki/pcm_psm/DataSetUnix/",folderName[wrapitr]))
  
  
  ######################################################################################
  ##Read in list of csv
  ######################################################################################
  # create an empty list of files
  file_list_all <- list.files()
  
  # get every file with "XX"
  file_list <- file_list_all[grepl("XX", file_list_all)]
  
  #############################################
  # loop over i files for calculations
  # This loop is used to fill act_list
  #############################################
  foreach (i = 1:length(file_list)) %dopar% { # start of i loop
    fill_active_list()
    print(i)  
  } # end of i loop
  
  ######################################################################################
  ##Analysis for Outcome CSV File
  ######################################################################################
  # create an empty matrix for outcomes
  dat_ana <- matrix(NA, nrow = length(act_list), ncol = 32) # B+ R+ CI
  nums <- matrix(NA, nrow = length(act_list), ncol = 1) # Number of PCs
  
  for (i in 1:length(act_list)){
    print(i)
    # 1. Genral info. for dataset
    dat_ana[i,1] <- as.numeric(substr(file_list[i], start = 3, stop = 4))
    dat_ana[i,2] <- as.numeric(substr(file_list[i], start = 8, stop = 9))
    dat_ana[i,3] <- as.numeric(substr(file_list[i], start = 13, stop = 14))
    dat_ana[i,4] <- as.numeric(substr(file_list[i], start = 16, stop = 19))
    
    # No ps case ----------------------------------------------------------------------------
    # Linear regressoin 
    if (dat_ana[i, 2] != 0){
      out_nops <- lm(Y~T, data=act_list[[i]])
      no_ps_b_hat <- out_nops$coefficients[2] # b_hat
      dat_ana[i, 5] <- no_ps_b_hat # b_hat
      dat_ana[i, 6] <- summary(out_nops)$r.squared # r squared
      #dat_ana[i, 6] <- summary(out_nops)$adj.r.squared # adjusted r squared
      dat_ana[i, 7] <- no_ps_b_hat - 1.96*summary(out_nops)$coefficients[4] # Lower CI
      dat_ana[i, 8] <- no_ps_b_hat + 1.96*summary(out_nops)$coefficients[4] # Upper CI
      
      # traditional ps----------------------------------------------------------------------------
      param_ps <- glm(as.formula(formula), data=act_list[[i]],family="binomial")
      act_list[[i]] <- act_list[[i]] %>%
        mutate(ps = predict(param_ps, type = "response"),
               odds = 1*T+ ps*(1-T)/(1- ps),
               iptw = ((mean(ps))/ ps)*T+((mean(1- ps))/(1- ps))*(1-T))
      
      # Linear regression that is weighted by odds
      out_odds <- lm(Y ~ T, data = act_list[[i]], weights = odds)
      ps_odd_b_hat <- out_odds$coefficients[2] # b_hat
      dat_ana[i, 9] <- ps_odd_b_hat
      dat_ana[i, 10] <- summary(out_odds)$r.squared # r squared
      #dat_ana[i, 10] <- summary(out_odds)$adj.r.squared # adjusted r squared
      dat_ana[i, 11] <- ps_odd_b_hat - 1.96*summary(out_odds)$coefficients[4] # Lower CI
      dat_ana[i, 12] <- ps_odd_b_hat + 1.96*summary(out_odds)$coefficients[4] # Upper CI  
      
      # Linear regression that is weighted bt iptw
      out_iptw <- lm(Y~ T, data = act_list[[i]], weights = iptw)
      ps_iptw_b_hat <- out_iptw$coefficients[2] # b_hat
      dat_ana[i, 13] <- ps_iptw_b_hat
      dat_ana[i, 14] <- summary(out_iptw)$r.squared # r squared
      #dat_ana[i, 14] <- summary(out_iptw)$adj.r.squared # adjusted r squared
      dat_ana[i, 15] <- ps_iptw_b_hat - 1.96*summary(out_iptw)$coefficients[4] # Lower CI
      dat_ana[i, 16] <- ps_iptw_b_hat + 1.96*summary(out_iptw)$coefficients[4] # Upper CI
      
      # PC1------------------------------------------------------------------------------------------
      # Calculate principal component score
      pca <- princomp(act_list[[i]][, -c(7, 8, 30:32)]) # princomp to get pc; 2:9 col been X2-X9, 12:56 col been all the square and interaction terms
      act_list[[i]] <- cbind(act_list[[i]], pca$scores) # cbind act_list with the Principal Component Scores
      
      # propensity score for Odds & iptw
      pc1 <- glm(T ~ Comp.1, data = act_list[[i]],family="binomial")
      act_list[[i]] <- act_list[[i]] %>%
        mutate(ps_pc1 = predict(pc1, type = "response"),
               odds_pc1 = 1*T+ ps_pc1*(1-T)/(1- ps_pc1),
               iptw_pc1 = ((mean(ps_pc1))/ ps_pc1)*T+((mean(1- ps_pc1))/(1- ps_pc1))*(1-T))
      
      # linear regression that is weighted by odds for pc1
      out_odds_pc1 <- lm(Y ~ T, data = act_list[[i]], weights = odds_pc1)
      pc1_odd_b_hat <- out_odds_pc1$coefficients[2] # b_hat
      dat_ana[i, 17] <- pc1_odd_b_hat
      dat_ana[i, 18] <- summary(out_odds_pc1)$r.squared # r squared
      #dat_ana[i, 18] <- summary(out_odds_pc1)$adj.r.squared # adjusted r squared
      dat_ana[i, 19] <- pc1_odd_b_hat - 1.96*summary(out_odds_pc1)$coefficients[4] # Lower CI
      dat_ana[i, 20] <- pc1_odd_b_hat + 1.96*summary(out_odds_pc1)$coefficients[4] # Upper CI  
      
      # linear regression that is weighted bt iptw
      out_iptw_pc1 <- lm(Y~ T, data = act_list[[i]], weights = iptw_pc1)
      pc1_iptw_b_hat <- out_iptw_pc1$coefficients[2] # b_hat
      dat_ana[i, 21] <- pc1_iptw_b_hat
      dat_ana[i, 22] <- summary(out_iptw_pc1)$r.squared # r squared
      #dat_ana[i, 22] <- summary(out_iptw_pc1)$adj.r.squared # adjusted r squared
      dat_ana[i, 23] <- pc1_iptw_b_hat - 1.96*summary(out_iptw_pc1)$coefficients[4] # Lower CI
      dat_ana[i, 24] <- pc1_iptw_b_hat + 1.96*summary(out_iptw_pc1)$coefficients[4] # Upper CI
      
      # PC2------------------------------------------------------------------------------------------
      # var over 40%
      var <- apply(pca$scores, 2, var)              # get the var
      props <- var/sum(var)                         # calculate proportion of var
      varpct <- cumsum(props)                       # calculate cumulative propotion of var
      num <- sum(varpct<varcut)                     # number of principal component socres will be used beyond Comp.1
      varname <- names(act_list[[i]])[33:(33+num)]  # extract names for Propensity Score Matching with PC2 (40%)
      
      # f<-T ~ Comp.1
      f <- as.formula(paste("T",paste(varname, collapse = " + "), sep =" ~ ")) # make the extraction as formula
      nums[i,] <- num + 1
      
      # odds & iptw using pc account for more than 40% variance
      pc2 <- glm(f, data = act_list[[i]],family="binomial")
      act_list[[i]] <- act_list[[i]] %>%
        mutate(ps_pc2 = predict(pc2, type = "response"),
               odds_pc2 = 1*T+ (ps_pc2/(1- ps_pc2))*(1-T),
               iptw_pc2 = ((mean(ps_pc2))/ ps_pc2)*T+((mean(1- ps_pc2))/(1- ps_pc2))*(1-T))
      
      # linear regression that is weighted by odds for pc2
      out_odds_pc2 <- lm(Y ~ T, data = act_list[[i]])#, weights = odds_pc2)
      pc2_odd_b_hat <- out_odds_pc2$coefficients[2] # b_hat
      dat_ana[i, 25] <- pc2_odd_b_hat
      dat_ana[i, 26] <- summary(out_odds_pc2)$r.squared # r squared
      #dat_ana[i, 26] <- summary(out_odds_pc2)$adj.r.squared # adjusted r squared
      dat_ana[i, 27] <- pc2_odd_b_hat - 1.96*summary(out_odds_pc2)$coefficients[4] # Lower CI
      dat_ana[i, 28] <- pc2_odd_b_hat + 1.96*summary(out_odds_pc2)$coefficients[4] # Upper CI  
      
      # linear regression that is weighted bt iptw
      out_iptw_pc2 <- lm(Y~ T, data = act_list[[i]], weights = iptw_pc2)
      pc2_iptw_b_hat <- out_iptw_pc2$coefficients[2] # b_hat
      dat_ana[i, 29] <- pc2_iptw_b_hat
      dat_ana[i, 30] <- summary(out_iptw_pc2)$r.squared # r squared
      #dat_ana[i, 30] <- summary(out_iptw_pc2)$adj.r.squared # adjusted r squared
      dat_ana[i, 31] <- pc2_iptw_b_hat - 1.96*summary(out_iptw_pc2)$coefficients[4] # Lower CI
      dat_ana[i, 32] <- pc2_iptw_b_hat + 1.96*summary(out_iptw_pc2)$coefficients[4] # Upper CI
      
    } else {
      dat_ana[i,9:32] <- 0
    } # end of loop j
  }
  
  # assigning proper column names for dat_ana
  colnames(dat_ana) <- c("XX", "TY", "TS", "IT", "Nops_b","Nops_r2",
                         "Nops_low","Nops_up","Odds_b","Odds_r2",
                         "Odds_low","Odds_up","Iptw","Iptw_r2",
                         "Iptw_low","Iptw_up", "Odds_b_PC1", "Odds_r2_PC1",
                         "Odds_low_PC1", "Odds_up_PC1", "Iptw_b_PC1",
                         "Iptw_r2_PC1","Iptw_low_PC1","Iptw_up_PC1",
                         "Odds_b_PC2", "Odds_r2_PC2",
                         "Odds_loW_PC2", "Odds_up_PC2", "Iptw_b_PC1",
                         "Iptw_r2_PC2","Iptw_low_PC2",	"Iptw_up_PC2"
  )
  
  # change the working directory as cur_dir to save result files in this directory
  setwd(cur_dir)

    
  write.csv(dat_ana, paste0("dat_ana_weight_",folderName[wrapitr],".csv"))
  write.csv(nums, paste0("numOfPCs_weight_",folderName[wrapitr],".csv")) 
}


# Initializing
act_list <<- list()
rm(list = ls(all.names = TRUE)) # clear all objects includes hidden objects.
gc() # free up memrory and report the memory usage.
