


# Initial setup/ REMOVE VARIABLES ----


rm(list = ls())    # remove any variables in R's memory 

# First, make sure you have the VGAM package installed and loaded
#install.packages("VGAM")
#install.packages("DirichletReg")

# Load packages

library(dplyr)
library(tidyr)
library(reshape2)   # For melting data
library(ggplot2)    # For plotting
library(ggrepel)    # For plotting
library(gridExtra)  # For plotting
library(ellipse)    # For plotting
library(scales)     # For dollar signs and commas
library(patchwork)  # For combining ggplot2 figures
library(dampack)  # Uncomment to use CEA and PSA visualization functionality from dampack instead of the functions included in this repository
library(darthtools) # Uncomment to use WCC, parameter transformation, and matrix checks from darthtools instead of the functions included in this repository
library(doParallel) # For running PSA in parallel
library(MASS) 
library(MCMCpack)
# Load supplementary functions

#source ("/Users/tessavernooij/Documents/EUHEM/Thesis/Malaria")


# Model input ----

## General setup ----


cycle_length <- 5/365 # cycle length 5 Days  
n_age_init <- 1 # age at baseline --> maybe start at age where children are completely vaccinated ? which is at 1 year and 7 months 
n_age_max  <- 5 # maximum age of follow up
n_cycles <- (n_age_max - n_age_init)/cycle_length # time horizon, number of cycles
## Age labels 
v_age_names <- paste(rep(n_age_init:(n_age_max-1), each = 1/cycle_length), 
                     1:(1/cycle_length), 
                     sep = ".")

## Health states of the model:
v_names_states <- c("W",  #well 
                    "MM", # Mild Malaria (MM)
                    "CM", # Cerebral Malaria (CM)
                    "SMA", # Severe Malaria Anemia (SMA)  
                    "D") # Death other causes  


n_states <- length(v_names_states)   # number of health states 

## Tunnel inputs ----

#Number of cycles for each tunnel state
n_tunnel_size_CM <- 2
n_tunnel_size_SMA <- 2

#### Vector with cycles for cerebral malaria tunnel
v_cycles_tunnel_CM <- 1:n_tunnel_size_CM

#### Vector with names for tunnel states of cerebral malaria
v_CM_tunnel <- paste("CM", seq(1, n_tunnel_size_CM), "", sep = "")

#### Vector with cycles for severe anemic malaria tunnel
v_cycles_tunnel_SMA <- 1:n_tunnel_size_SMA

#### Vector with names for tunnel states of severe anemic malaria
v_SMA_tunnel <- paste("SMA", seq(1, n_tunnel_size_SMA), "", sep = "")

#### Create variables for model with tunnels
v_names_states_tunnels <- c("W", "MM", v_CM_tunnel, v_SMA_tunnel, "D") # health state names
n_states_tunnels <- length(v_names_states_tunnels)         # number of health states

## Discounting factors ----

d_c <- 0.03# annual discount rate for costs 
d_e <- 0.03 # annual discount rate for DALYs



## Strategies ----

v_names_str <- c("LLIN", 
                 "R21") 
n_str       <- length(v_names_str)        # number of strategieS


gen_wcc <- function(n_cycles, method = "half-cycle") {
  if (method == "half-cycle") {
    v_wcc <- rep(1, n_cycles + 1)
    v_wcc[1] <- v_wcc[n_cycles + 1] <- 0.5
    return(v_wcc)  # Return v_wcc
  } else {
    stop("Invalid method. Please choose 'half-cycle'.")
  }
}

### detailed Costs ----
c_consult_HF <- 3.33 
c_consult_CF <- 0.80
diagnost_HF_RDT <- 2.54
Diagnostic_HF_M <- 0.92
Diagnostic_CHW_RDT <- 2.54
Diagnostic_CHW_M <- 0
HF_treatconsult <- 3.33
CHW_treatconsult <- 0.80 
Follow_up_con <- 8.27
Anti_mal_treat_MM <- 0.188
Anti_mal_treat_SM <- 8.0251
CM_treat <- 1.49
SMA_treat <- 8.74
Trans_pat <- 3.59
Trans_acomp <- 0.83
wage_loss <-11.57
time_loss_work_MM <- 6.9
time_loss_work_SM <- 10
LLIN_unit <- 0.68
R21_MatrixM <- 12

LLIN_duration <- 5 
cov_consult_HF <- 0.83
cov_consult_CF <- 0.16
cov_diagnostic_HF_RDT <-0.17
cov_diagnostic_HF_M <- 0.66
cov_diagnostic_CHW_RDT <- 0.14
cov_HFtreatconsult <- 0.83
cov_CHW_treatconsult <- 0.16
cov_LLIN <- 0.37
cov_R21 <-0.83


c_W <- ((LLIN_unit)*cov_LLIN)

c_MM <- (c_consult_HF*cov_consult_HF) + (c_consult_CF*cov_consult_CF) + (diagnost_HF_RDT*cov_diagnostic_HF_RDT) +
  (Diagnostic_HF_M*cov_diagnostic_HF_M) + (Diagnostic_CHW_RDT*cov_diagnostic_CHW_RDT) +Diagnostic_CHW_M + 
  (HF_treatconsult*cov_HFtreatconsult)+ (CHW_treatconsult*cov_CHW_treatconsult )+ Follow_up_con+ 
  Anti_mal_treat_MM+ Trans_pat+ Trans_acomp+(wage_loss*time_loss_work_MM)       

c_SMA <- (c_consult_HF*cov_consult_HF) + (c_consult_CF*cov_consult_CF) + (diagnost_HF_RDT*cov_diagnostic_HF_RDT) + 
  (Diagnostic_HF_M*cov_diagnostic_HF_M) + (Diagnostic_CHW_RDT*cov_diagnostic_CHW_RDT)+Diagnostic_CHW_M + 
  (HF_treatconsult*cov_HFtreatconsult)+ (CHW_treatconsult*cov_CHW_treatconsult )+ Follow_up_con+ Anti_mal_treat_SM+ 
  SMA_treat+ Trans_pat+ Trans_acomp 


c_CM <- (c_consult_HF*cov_consult_HF) + (c_consult_CF*cov_consult_CF) + (diagnost_HF_RDT*cov_diagnostic_HF_RDT) + 
  (Diagnostic_HF_M*cov_diagnostic_HF_M) + (Diagnostic_CHW_RDT*cov_diagnostic_CHW_RDT)+Diagnostic_CHW_M +
  (HF_treatconsult*cov_HFtreatconsult)+ (CHW_treatconsult*cov_CHW_treatconsult )+ Follow_up_con+ Anti_mal_treat_SM+
  CM_treat+ Trans_pat+ Trans_acomp 

c_W_R21 <-  ((LLIN_unit)*cov_LLIN) + (R21_MatrixM*cov_R21)

c_MM_R21 <- (c_consult_HF*cov_consult_HF) + (c_consult_CF*cov_consult_CF) + (diagnost_HF_RDT*cov_diagnostic_HF_RDT) +
  (Diagnostic_HF_M*cov_diagnostic_HF_M) + (Diagnostic_CHW_RDT*cov_diagnostic_CHW_RDT)+Diagnostic_CHW_M + 
  (HF_treatconsult*cov_HFtreatconsult)+ (CHW_treatconsult*cov_CHW_treatconsult)+ Follow_up_con+ Anti_mal_treat_MM+ 
  Trans_pat+ Trans_acomp+(wage_loss*time_loss_work_MM)

c_SMA_R21 <- (c_consult_HF*cov_consult_HF) + (c_consult_CF*cov_consult_CF) + (diagnost_HF_RDT*cov_diagnostic_HF_RDT) + 
  (Diagnostic_HF_M*cov_diagnostic_HF_M) + (Diagnostic_CHW_RDT*cov_diagnostic_CHW_RDT)+Diagnostic_CHW_M + 
  (HF_treatconsult*cov_HFtreatconsult)+ (CHW_treatconsult*cov_CHW_treatconsult )+ Follow_up_con+ Anti_mal_treat_SM+ 
  SMA_treat+ Trans_pat+ Trans_acomp

c_CM_R21 <- (c_consult_HF*cov_consult_HF) + (c_consult_CF*cov_consult_CF) + (diagnost_HF_RDT*cov_diagnostic_HF_RDT) + 
  (Diagnostic_HF_M*cov_diagnostic_HF_M) + (Diagnostic_CHW_RDT*cov_diagnostic_CHW_RDT)+Diagnostic_CHW_M +
  (HF_treatconsult*cov_HFtreatconsult)+ (CHW_treatconsult*cov_CHW_treatconsult )+ Follow_up_con+ Anti_mal_treat_SM+ 
  CM_treat+ Trans_pat+ Trans_acomp


## from the malawi paper, most of the time you leave its suppose to last one cycle but of course poeple mihght stay a few more days in that health state 
efficacy_R21_MM <- 0.73 
efficacy_R21_SM <- 0.67 
efficacy_LLIN <- 0.55

### DALYs ---- overall disease burden 
W_MM <- 0.19
W_CM <-0.47
W_SMA <-0.47
Duration_MM <-0.014
Duration_CM <-0.041
Duration_SMA <-0.041
N_cases_MM <-801
N_cases_CM <-88.11
N_cases_SMA <-104.13

W_MM_R21 <-0.19
W_CM_R21 <-0.47
W_SMA_R21 <-0.47
Duration_MM_R21 <-0.014
Duration_CM_R21 <-0.041
Duration_SMA_R21 <-0.041
N_cases_MM_R21 <- 519
N_cases_CM_R21 <-57
N_cases_SMA_R21 <-67

DALY_W   <- 0  # you are healthy so 0 represents perfect health 
DALY_MM   <- W_MM * Duration_MM * N_cases_MM  # annual DALY of having moderate malaria
DALY_CM   <- W_CM * Duration_CM * N_cases_CM   # annual DALY of having cerebral malaria
DALY_SMA   <- W_SMA *Duration_SMA * N_cases_SMA    # annual DALY of being severe malaria anemia   
DALY_W_R21  <- 0  # you are healthy so 0 represents perfect health 
DALY_MM_R21   <- W_MM_R21 * Duration_MM_R21 * N_cases_MM_R21# annual DALY of having moderate malaria
DALY_CM_R21   <- W_CM_R21 * Duration_CM_R21 * N_cases_CM_R21# annual DALY of having cerebral malaria
DALY_SMA_R21   <- W_SMA_R21 * Duration_SMA_R21*  N_cases_SMA_R21   # annual DALY of being severe malaria anemia   


##Discount weight for costs and effects ---
v_dwc  <- 1 / ((1 + (d_e * cycle_length)) ^ (0:n_cycles))
v_dwe  <- 1 / ((1 + (d_c * cycle_length)) ^ (0:n_cycles))


p_WMM  <- 0.0001169687# constant annual probability of getting mild malaria when Well conditional on surviving 
p_MMCM  <- 0.00160
p_MMSMA  <- 0.00853 
#p_MMSMA  <- 0.00853*2
p_CMW <- 0.85322 ##also use malawi text 
p_SMAW <- 0.89122 ##also use malawi text 

p_WMM_R21 <- p_WMM *(1-(cov_R21*efficacy_R21_MM))
p_MMCM_R21 <- p_MMCM *(1-(cov_R21*efficacy_R21_SM))
p_MMSMA_R21 <-p_MMSMA *(1-(cov_R21*efficacy_R21_SM))

p_CMW_R21 <- 0.85322 ##
p_SMAW_R21 <- 0.89122 ## 
p_D_natural <- (0.08 /365)*5 ## 

# Construct state-transition models ---- see table already
## Initial state vector ----
v_m_init_tunnels <- c(1000, 0, rep(0, n_tunnel_size_CM + n_tunnel_size_SMA), 0) 

## Initialize cohort traces ----
### Initialize cohort trace for state-residece dependent cSTM under SoC----
m_M_tunnels_LLIN <- matrix(0, 
                           nrow   = (n_cycles + 1), ncol = n_states_tunnels, 
                           dimnames = list(0:n_cycles, v_names_states_tunnels))

#* Store the initial state vector in the first row of the cohort trace
m_M_tunnels_LLIN[1, ] <- v_m_init_tunnels

### Initialize cohort trace for strategies R21
#* Structure and initial states are the same as for SoC
m_M_tunnels_R21  <- m_M_tunnels_LLIN # Strategy R21

## Create transition probability arrays for strategy no intervention ----
### Initialize transition probability array for strategy SoC ----
##'s a three-dimensional array with dimensions matching the number of health states (n_states_tunnels) and cycles (n_cycles). The dimnames argument provides labels for the dimensions, which is good practice for clarity and readability.
#* All transitions to a non-death state are assumed to be conditional on survival
a_P_tunnels_LLIN <- array(0, 
                          dim = c(n_states_tunnels, n_states_tunnels, n_cycles),
                          dimnames = list(v_names_states_tunnels, 
                                          v_names_states_tunnels, 
                                          0:(n_cycles-1)))



###################################################################
###  Initialize transition probability array for strategy in array LLIN

## From W 
a_P_tunnels_LLIN ["W", "W", ]               <-((1-p_D_natural) *(1- p_WMM))
a_P_tunnels_LLIN ["W", "MM", ]              <-(1-p_D_natural)* p_WMM 
a_P_tunnels_LLIN ["W", v_CM_tunnel, ]       <- 0
a_P_tunnels_LLIN ["W", v_SMA_tunnel, ]      <- 0
a_P_tunnels_LLIN ["W", "D", ]               <- p_D_natural

## From MM
a_P_tunnels_LLIN["MM", "W", ]               <- ((1-p_D_natural)*(1-(p_MMCM+p_MMSMA)))
a_P_tunnels_LLIN["MM", "MM", ]              <- 0 
a_P_tunnels_LLIN["MM", v_CM_tunnel[1], ]    <- (1-p_D_natural)*p_MMCM
a_P_tunnels_LLIN["MM", v_SMA_tunnel[1], ]   <- (1-p_D_natural) *p_MMSMA
a_P_tunnels_LLIN["MM", v_CM_tunnel[2], ]    <- 0  # Transition probability to the second tunnel set to 0
a_P_tunnels_LLIN["MM", v_SMA_tunnel[2], ]   <- 0  # Transition probability to the second tunnel set to 0
a_P_tunnels_LLIN ["MM", "D", ]              <- p_D_natural


## From CM
# Transition probability within the first CM tunnel NOT NEEDED BECAUSE IT'S A TUNNEL; YOU ALWAYS MOVE.
# Transition probability from the first CM tunnel to the "well" state
a_P_tunnels_LLIN[v_CM_tunnel[1], "W", ]                   <- (1-p_D_natural)* p_CMW
a_P_tunnels_LLIN[v_CM_tunnel[1], "MM", ]                   <- 0
a_P_tunnels_LLIN[v_CM_tunnel[1], v_CM_tunnel[2], ]        <-  (1-p_D_natural)*(1-p_CMW) 
a_P_tunnels_LLIN[v_CM_tunnel[1], v_CM_tunnel[1], ]        <- 0
a_P_tunnels_LLIN[v_CM_tunnel[1], v_SMA_tunnel[1], ]        <- 0
a_P_tunnels_LLIN[v_CM_tunnel[1], v_SMA_tunnel[2], ]        <- 0
a_P_tunnels_LLIN[v_CM_tunnel[1], "D", ]        <- p_D_natural

# Transition probability from the second CM tunnel to the "well" state
a_P_tunnels_LLIN[v_CM_tunnel[2], "W", ]                   <- (1-p_D_natural)
a_P_tunnels_LLIN[v_CM_tunnel[2], "MM", ]                   <- 0
a_P_tunnels_LLIN[v_CM_tunnel[2], v_CM_tunnel[2], ]        <- 0
a_P_tunnels_LLIN[v_CM_tunnel[2], v_CM_tunnel[1], ]          <- 0
a_P_tunnels_LLIN[v_CM_tunnel[2], v_SMA_tunnel[1], ]          <- 0
a_P_tunnels_LLIN[v_CM_tunnel[2], v_SMA_tunnel[2], ]          <- 0
a_P_tunnels_LLIN[v_CM_tunnel[2], "D", ]        <- p_D_natural


## From SMA
# Transition probability within the first CM tunnel NOT NEEDED BECAUSE IT'S A TUNNEL; YOU ALWAYS MOVE.
# Transition probability from the first CM tunnel to the "well" state
a_P_tunnels_LLIN[v_SMA_tunnel[1], "W", ]                <- (1-p_D_natural)*p_SMAW
a_P_tunnels_LLIN[v_SMA_tunnel[1], "MM", ]               <- 0
a_P_tunnels_LLIN[v_SMA_tunnel[1], v_SMA_tunnel[2], ]    <- (1-p_D_natural)* (1-p_SMAW ) 
a_P_tunnels_LLIN[v_SMA_tunnel[1], v_SMA_tunnel[1], ]    <- 0 
a_P_tunnels_LLIN[v_SMA_tunnel[1], v_CM_tunnel[1], ]     <- 0
a_P_tunnels_LLIN[v_SMA_tunnel[1], v_CM_tunnel[2], ]     <- 0
a_P_tunnels_LLIN[v_SMA_tunnel[1], "D", ]        <- p_D_natural

# Transition probability from the second CM tunnel to the "well" state
a_P_tunnels_LLIN[v_SMA_tunnel[2], "W", ]         <-  (1-p_D_natural) * 1
a_P_tunnels_LLIN[v_SMA_tunnel[2], "MM", ]        <- 0
a_P_tunnels_LLIN[v_SMA_tunnel[2], v_SMA_tunnel[2], ]   <- 0
a_P_tunnels_LLIN[v_SMA_tunnel[2], v_SMA_tunnel[1], ] <- 0 
a_P_tunnels_LLIN[v_SMA_tunnel[2], v_CM_tunnel[1], ]    <- 0
a_P_tunnels_LLIN[v_SMA_tunnel[2], v_CM_tunnel[2], ]    <- 0
a_P_tunnels_LLIN[v_SMA_tunnel[2], "D", ]        <- p_D_natural

## From D
a_P_tunnels_LLIN["D", "W", ]              <- 0
a_P_tunnels_LLIN["D", "MM", ]             <- 0
a_P_tunnels_LLIN["D", "D", ]              <- 1
a_P_tunnels_LLIN["D", v_CM_tunnel[2], ]   <- 0
a_P_tunnels_LLIN["D", v_CM_tunnel[1], ]   <- 0
a_P_tunnels_LLIN["D", v_SMA_tunnel[1], ]  <- 0
a_P_tunnels_LLIN["D", v_SMA_tunnel[2], ]  <- 0


print("Transition Probabilities:")
print(a_P_tunnels_LLIN)

## Check that transition probabilities are [0, 1]
check_transition_probability(a_P_tunnels_LLIN, verbose = TRUE)
# Check the sum of transition probabilities
check_sum_of_transition_array(a_P_tunnels_LLIN, n_states = n_states_tunnels, n_cycles = n_cycles, verbose = TRUE)


###################################################################

# Initialize transition probability array for strategy R21  ----
a_P_tunnels_R21 <- a_P_tunnels_LLIN
###################################################################
###  Initialize transition probability array for strategy in array LLIN
## From W 
a_P_tunnels_R21 ["W", "W", ]               <- (1-p_D_natural)*(1 - p_WMM_R21) 
a_P_tunnels_R21 ["W", "MM", ]              <- (1-p_D_natural)*p_WMM_R21 
a_P_tunnels_R21 ["W", v_CM_tunnel, ]       <- 0
a_P_tunnels_R21 ["W", v_SMA_tunnel, ]      <- 0
a_P_tunnels_R21 ["W", "D", ]               <- p_D_natural


## From MM
a_P_tunnels_R21["MM", "W", ]               <- (1-p_D_natural)*(1-(p_MMCM_R21+p_MMSMA_R21))
a_P_tunnels_R21["MM", "MM", ]              <- 0
a_P_tunnels_R21["MM", v_CM_tunnel[1], ]    <-  (1-p_D_natural)*p_MMCM_R21
a_P_tunnels_R21["MM", v_SMA_tunnel[1], ]   <-  (1-p_D_natural)*p_MMSMA_R21
a_P_tunnels_R21["MM", v_CM_tunnel[2], ]    <- 0  # Transition probability to the second tunnel set to 0
a_P_tunnels_R21["MM", v_SMA_tunnel[2], ]   <- 0  # Transition probability to the second tunnel set to 0
a_P_tunnels_R21 ["MM", "D", ]               <- p_D_natural

## From CM
# Transition probability within the first CM tunnel NOT NEEDED BECAUSE IT'S A TUNNEL; YOU ALWAYS MOVE.
# Transition probability from the first CM tunnel to the "well" state
a_P_tunnels_R21[v_CM_tunnel[1], "W", ]                   <- (1-p_D_natural)*  p_CMW_R21
a_P_tunnels_R21[v_CM_tunnel[1], "MM", ]                   <- 0
a_P_tunnels_R21[v_CM_tunnel[1], v_CM_tunnel[2], ]        <-  (1-p_D_natural)* (1-p_CMW_R21) 
a_P_tunnels_R21[v_CM_tunnel[1], v_CM_tunnel[1], ]        <- 0
a_P_tunnels_R21[v_CM_tunnel[1], v_SMA_tunnel[1], ]        <- 0
a_P_tunnels_R21[v_CM_tunnel[1], v_SMA_tunnel[2], ]        <- 0
a_P_tunnels_R21[v_CM_tunnel[1], "D", ]        <- p_D_natural

# Transition probability from the second CM tunnel to the "well" state
a_P_tunnels_R21[v_CM_tunnel[2], "W", ]                   <- (1-p_D_natural)* 1
a_P_tunnels_R21[v_CM_tunnel[2], "MM", ]                   <- 0
a_P_tunnels_R21[v_CM_tunnel[2], v_CM_tunnel[2], ]        <- 0
a_P_tunnels_R21[v_CM_tunnel[2], v_CM_tunnel[1], ]          <- 0
a_P_tunnels_R21[v_CM_tunnel[2], v_SMA_tunnel[1], ]          <- 0
a_P_tunnels_R21[v_CM_tunnel[2], v_SMA_tunnel[2], ]          <- 0
a_P_tunnels_R21[v_CM_tunnel[2], "D", ]        <- p_D_natural

## From SMA
# Transition probability within the first CM tunnel NOT NEEDED BECAUSE IT'S A TUNNEL; YOU ALWAYS MOVE.
# Transition probability from the first CM tunnel to the "well" state
a_P_tunnels_R21[v_SMA_tunnel[1], "W", ]                <-  (1-p_D_natural)* p_SMAW_R21
a_P_tunnels_R21[v_SMA_tunnel[1], "MM", ]               <- 0
a_P_tunnels_R21[v_SMA_tunnel[1], v_SMA_tunnel[2], ]    <- (1-p_D_natural)* (1-p_SMAW_R21 ) 
a_P_tunnels_R21[v_SMA_tunnel[1], v_SMA_tunnel[1], ]    <- 0 
a_P_tunnels_R21[v_SMA_tunnel[1], v_CM_tunnel[1], ]     <- 0
a_P_tunnels_R21[v_SMA_tunnel[1], v_CM_tunnel[2], ]     <- 0
a_P_tunnels_R21[v_SMA_tunnel[1], "D", ]        <- p_D_natural

# Transition probability from the second CM tunnel to the "well" state
a_P_tunnels_R21[v_SMA_tunnel[2], "W", ]         <- (1-p_D_natural)* 1
a_P_tunnels_R21[v_SMA_tunnel[2], "MM", ]        <- 0
a_P_tunnels_R21[v_SMA_tunnel[2], v_SMA_tunnel[2], ]   <-  0
a_P_tunnels_R21[v_SMA_tunnel[2], v_SMA_tunnel[1], ] <- 0 
a_P_tunnels_R21[v_SMA_tunnel[2], v_CM_tunnel[1], ]    <- 0
a_P_tunnels_R21[v_SMA_tunnel[2], v_CM_tunnel[2], ]    <- 0
a_P_tunnels_R21[v_SMA_tunnel[2], "D", ]        <- p_D_natural

## From D
a_P_tunnels_R21["D", "W", ]              <- 0
a_P_tunnels_R21["D", "MM", ]             <- 0
a_P_tunnels_R21["D", "D", ]              <- 1
a_P_tunnels_R21["D", v_CM_tunnel[2], ]   <- 0
a_P_tunnels_R21["D", v_CM_tunnel[1], ]   <- 0
a_P_tunnels_R21["D", v_SMA_tunnel[1], ]  <- 0
a_P_tunnels_R21["D", v_SMA_tunnel[2], ]  <- 0


print("Transition Probabilities:")
print(a_P_tunnels_R21)

## Check that transition probabilities are [0, 1]
check_transition_probability(a_P_tunnels_R21, verbose = TRUE)
# Check the sum of transition probabilities
check_sum_of_transition_array(a_P_tunnels_R21, n_states = n_states_tunnels, n_cycles = n_cycles, verbose = TRUE)


# Create transition dynamics arrays ----

#* These arrays will capture transitions from each state to another over time 
### Initialize transition dynamics array for strategy SoC ----
a_A_tunnels_LLIN <- array(0,
                          dim = c(n_states_tunnels, n_states_tunnels, n_cycles + 1),
                          dimnames = list(v_names_states_tunnels, v_names_states_tunnels, 0:n_cycles))
#* Set first slice of Array with the initial state vector in its diagonal
diag(a_A_tunnels_LLIN[, , 1]) <- v_m_init_tunnels
### Initialize transition-dynamics array  ----
#* Structure and initial states are the same as for SoC
a_A_tunnels_R21  <- a_A_tunnels_LLIN


#  Run Markov model ----updates cohort traces and transition-dynamics arrays for multiple interventions in 
# a Markov model, calculating the next state probabilities based on transition probabilities and the current
# state of the cohort
#* Iterative solution of state-residence dependent 
#* cycle length , states dimes.
for (t in 1:n_cycles) {
  ## Fill in cohort trace
  # For no intervention 
  m_M_tunnels_LLIN[t + 1, ]   <- m_M_tunnels_LLIN[t, ] %*% a_P_tunnels_LLIN [, , t]
  # For R21
  m_M_tunnels_R21[t + 1, ]      <- m_M_tunnels_R21[t, ] %*% a_P_tunnels_R21 [, , t]
  
  ## Fill in transition-dynamics array
  # For no intervention
  a_A_tunnels_LLIN[, , t + 1]   <- diag(m_M_tunnels_LLIN[t, ]) %*% a_P_tunnels_LLIN [, , t]
  # For R21
  a_A_tunnels_R21[, , t + 1]   <- diag(m_M_tunnels_R21[t, ]) %*% a_P_tunnels_R21[, , t]
  
}

m_M_tunnels_LLIN_sum <- cbind( W   = m_M_tunnels_LLIN[, "W"], 
                               MM  = m_M_tunnels_LLIN[, "MM"], 
                               CM  = rowSums ( m_M_tunnels_LLIN [, 3:4]), 
                               SMA = rowSums ( m_M_tunnels_LLIN [, 5:6]),
                               D= m_M_tunnels_LLIN[, "D"])

m_M_tunnels_R21_sum <- cbind( W   = m_M_tunnels_R21[, "W"], 
                              MM  = m_M_tunnels_R21[, "MM"], 
                              CM  = rowSums(m_M_tunnels_R21[, 3:4]), 
                              SMA = rowSums(m_M_tunnels_R21[, 5:6]),
                              D= m_M_tunnels_R21[, "D"])

## Store the cohort traces in a list ----
l_m_M <- list(m_M_tunnels_LLIN_sum,      
              m_M_tunnels_R21_sum)
names(l_m_M) <- v_names_str

## Store the transition dynamics array for each strategy in a list ----
l_a_A <- list(a_A_tunnels_LLIN,
              a_A_tunnels_R21)
names(l_a_A) <- v_names_str


## Plot the cohort trace for strategy SoC ----
plot_trace(m_M_tunnels_LLIN_sum)
plot_trace(m_M_tunnels_R21_sum)
## Plot the cohort trace for all strategies ----
plot_trace_strategy(l_m_M)



# Vector of DALY per cycle for SMA under strategy LLIN
v_DALY_SMA_LLIN <- rep(DALY_SMA, n_tunnel_size_SMA)
names(v_DALY_SMA_LLIN) <- v_SMA_tunnel[-1]  # Exclude the first tunnel state/ exclusion of the first tunnel state from the names() assignment is to ensure that the names assigned to the elements of the vector correspond correctly to the tunnel states.

# Vector of DALY per cycle for CM under strategy LLIN
v_DALY_CM_LLIN <- rep(DALY_CM, n_tunnel_size_CM)
names(v_DALY_CM_LLIN) <- v_CM_tunnel[-1]  # Exclude the first tunnel state

# Vector of DALY under strategy LLIN
v_DALY_LLIN <- c(W = DALY_W, 
                 MM= DALY_MM,
                 SMA = v_DALY_CM_LLIN,
                 CM = v_DALY_SMA_LLIN,
                 D=0) 



# Initialize vector to store wage loss adjustments for each cycle in SMA tunnel
wage_loss_per_tunnel_SM <- rep((wage_loss * time_loss_work_SM) / 2, n_tunnel_size_SMA)

# Adjust wage loss for the second tunnel
wage_loss_per_tunnel_SM[2] <- (wage_loss * time_loss_work_SM) / 2

# Calculate cost vectors for SMA without including wage loss
v_c_SMA_LLIN_without_wageloss <- rep(c_SMA, n_tunnel_size_SMA)

# Calculate total costs for SMA, including wage loss adjustments for each cycle
v_c_SMA_LLIN <- ifelse(1:n_tunnel_size_SMA == 1, wage_loss_per_tunnel_SM[1] + c_SMA, wage_loss_per_tunnel_SM) 


# Initialize vector to store wage loss adjustments for each cycle in SMA tunnel
wage_loss_per_tunnel_SM <- rep((wage_loss * time_loss_work_SM) / 2, n_tunnel_size_CM)

# Adjust wage loss for the second tunnel
wage_loss_per_tunnel_SM[2] <- (wage_loss * time_loss_work_SM) / 2

# Calculate cost vectors for SMA without including wage loss
v_c_CM_LLIN_without_wageloss <- rep(c_CM, n_tunnel_size_CM)

# Calculate total costs for SMA, including wage loss adjustments for each cycle
v_c_CM_LLIN <- ifelse(1:n_tunnel_size_CM == 1, wage_loss_per_tunnel_SM[1] + c_SMA, wage_loss_per_tunnel_SM) 


# Define cost vector for individuals who are well and using LLIN
v_c_W <- c_W

# Vector of costs under strategy No Intervention
v_c_LLIN <- c(W = v_c_W, 
              MM = c_MM,
              SMA = v_c_SMA_LLIN,
              CM = v_c_CM_LLIN, 
              D = 0)



# Now, let's create vectors for "R21" strategy
# (You can adapt this part according to your specific DALY and cost values for "R21")

# Vector of DALYS per cycle for SMA under strategy R21
v_DALY_SMA_R21 <- rep(DALY_SMA_R21, n_tunnel_size_SMA)
names(v_DALY_SMA_R21) <- v_SMA_tunnel[-1]  # Exclude the first tunnel state

# Vector of DALYS per cycle for CM under strategy R21
v_DALY_CM_R21 <- rep(DALY_CM_R21, n_tunnel_size_CM)
names(v_DALY_CM_R21) <- v_CM_tunnel[-1]  # Exclude the first tunnel state

# Vector of DALYS under strategy R21
v_DALY_R21 <- c(W = DALY_W_R21, 
                MM= DALY_MM_R21,
                SMA = v_DALY_SMA_R21,
                CM = v_DALY_CM_R21, 
                D= 0 ) 



# Initialize vector to store wage loss adjustments for each cycle in SMA tunnel
wage_loss_per_tunnel_SM_R21 <- rep((wage_loss * time_loss_work_SM) / 2, n_tunnel_size_SMA)

# Adjust wage loss for the second tunnel
wage_loss_per_tunnel_SM_R21[2] <- (wage_loss * time_loss_work_SM) / 2

# Calculate cost vectors for SMA without including wage loss
v_c_SMA_without_wageloss_R21 <- rep(c_SMA_R21, n_tunnel_size_SMA)

# Calculate total costs for SMA, including wage loss adjustments for each cycle
v_c_SMA_R21 <- ifelse(1:n_tunnel_size_SMA == 1, wage_loss_per_tunnel_SM[1] + c_SMA_R21, wage_loss_per_tunnel_SM) 

# Initialize vector to store wage loss adjustments for each cycle in SMA tunnel
wage_loss_per_tunnel_SM <- rep((wage_loss * time_loss_work_SM) / 2, n_tunnel_size_CM)

# Adjust wage loss for the second tunnel
wage_loss_per_tunnel_SM[2] <- (wage_loss * time_loss_work_SM) / 2

# Calculate cost vectors for SMA without including wage loss
v_c_CM_without_wageloss_R21 <- rep(c_CM_R21, n_tunnel_size_CM)

# Calculate total costs for SMA, including wage loss adjustments for each cycle
v_c_CM_R21 <- ifelse(1:n_tunnel_size_CM == 1, wage_loss_per_tunnel_SM[1] + c_SMA_R21, wage_loss_per_tunnel_SM) 

v_c_W_R21 <- c_W_R21

# Vector of costs under strategy R21
v_c_R21 <- c(W =v_c_W_R21, 
             MM= c_MM_R21,
             SMA = v_c_SMA_R21,
             CM = v_c_CM_R21,
             D= 0)  


## Store state rewards ----
#* Store the vectors of state DALY for each strategy in a list 
l_DALY   <- list(v_DALY_LLIN,
                 v_DALY_R21)
#* Store the vectors of state cost for each strategy in a list 
l_c   <- list(v_c_LLIN,
              v_c_R21)
#* assign strategy names to matching items in the lists
names(l_DALY) <- names(l_c) <- v_names_str


# Compute expected outcomes ----
#* Create empty vectors to store total DALYs and costs 
v_tot_DALY <- v_tot_cost <- vector(mode = "numeric", length = n_str)
names(v_tot_DALY) <- names(v_tot_cost) <- v_names_str


## Loop through each strategy and calculate total DALY and costs ----
for (i in 1:n_str) { # i <- 1
  v_DALY_str         <- l_DALY[[i]]   # select the vector of state utilities for the i-th strategy
  v_c_str         <- l_c[[i]]   # select the vector of state costs for the i-th strategy
  a_A_tunnels_str <- l_a_A[[i]] # select the transition array for the i-th strategy
  
  #* Create transition matrices of state DALY and state costs for the i-th strategy
  m_DALY_str  <- matrix(v_DALY_str, nrow = n_states_tunnels, ncol = n_states_tunnels, byrow = T)
  m_c_str  <- matrix(v_c_str, nrow = n_states_tunnels, ncol = n_states_tunnels, byrow = T)
  
  #* Expand the transition matrix of state utilities across cycles to form a transition array of state DALYs
  a_R_DALY_str <- array(m_DALY_str, 
                        dim      = c(n_states_tunnels, n_states_tunnels, n_cycles + 1),
                        dimnames = list(v_names_states_tunnels, v_names_states_tunnels, 0:n_cycles))
  # Expand the transition matrix of state costs across cycles to form a transition array of state costs
  a_R_c_str <- array(m_c_str, 
                     dim      = c(n_states_tunnels, n_states_tunnels, n_cycles + 1),
                     dimnames = list(v_names_states_tunnels, v_names_states_tunnels, 0:n_cycles))
  
  
  ###* Expected DALYs and costs per cycle
  ##* Vector of DALYs and costs
  v_DALY_str <- apply(a_A_tunnels_str* a_R_DALY_str , 3, sum) # sum the proportion of the cohort across transitions 
  v_c_str <- apply(a_A_tunnels_str*a_R_c_str, 3, sum) # sum the proportion of the cohort across transitions
  
  ##* Discounted total expected DALYs and Costs per strategy & half cycle correction 
  #* DALYs
  v_tot_DALY[i] <- t(v_DALY_str) %*% (v_dwe* gen_wcc(n_cycles, method = "half-cycle"))
  #* Costs
  v_tot_cost[i] <- t(v_c_str) %*% (v_dwc* gen_wcc(n_cycles, method = "half-cycle"))
  
}

# Cost-effectiveness analysis (CEA) ----
## Incremental cost-effectiveness ratios (ICERs) ----
#* Function included in "R/Functions.R"; depends on the `dplyr` package
#* The latest version can be found in `dampack` package
df_cea <- calculate_icers(cost = v_tot_cost, 
                          effect = v_tot_DALY,
                          strategies = v_names_str)

## CEA table in proper format ----
table_cea <- format_table_cea(df_cea)
table_cea



###################################################################################
############################# Insidence of malaria ################################
###################################################################################

# Define parameters
efficacy <- 0.75
coverage <- 0.82

# Original probabilities
p_WMM  <- 0.0001169687 
p_MMCM  <- 0.00160
p_MMSMA  <- 0.00853 

# Adjusted probabilities considering vaccination
p_WMM_adjusted <- p_WMM * (1 - (efficacy * coverage))
p_MMCM_adjusted <- p_MMCM * (1 - (efficacy * coverage))
p_MMSMA_adjusted <- p_MMSMA * (1 - (efficacy * coverage))

# Population and infection data
total_population <- 153810
total_infected_year <- 801
total_infected_week <- 15.4

# Weekly infections #extracted from data set
infections_per_week <- c(12, 10, 11, 7, 11, 9, 20, 21, 25, 27, 18, 15, 16, 11, 20, 16, 15, 15, 6, 11, 17, 15, 20, 18, 12,
                         4, 9, 18, 23, 15, 12, 20, 12, 22, 14, 21, 29, 7, 10, 17, 20, 17, 15, 12, 20, 14, 12, 20, 15, 20, 14, 11)

# Weeks
weeks <- 1:length(infections_per_week)

# Calculate the number of cases for each type of malaria before vaccination
total_uncomplicated <- total_infected_year
total_cerebral <- total_infected_year * p_MMCM
total_anemic <- total_infected_year * p_MMSMA

# Calculate the number of cases prevented by vaccination
cases_uncomplicated_prevented <- total_uncomplicated * (efficacy * coverage)
cases_CM_prevented <- total_cerebral * (efficacy * coverage)
cases_SMA_prevented <- total_anemic * (efficacy * coverage)

# Update the total number of infected individuals after vaccination
total_uncomplicated_after_vaccination <- total_uncomplicated - cases_uncomplicated_prevented
total_cerebral_after_vaccination <- total_cerebral - cases_CM_prevented
total_anemic_after_vaccination <- total_anemic - cases_SMA_prevented

# Calculate weekly infections before and after vaccination
weekly_uncomplicated_before <- infections_per_week
weekly_uncomplicated_after <- infections_per_week * (1 - (efficacy * coverage))

weekly_cerebral_before <- infections_per_week * p_MMCM
weekly_cerebral_after <- infections_per_week * p_MMCM_adjusted

weekly_anemic_before <- infections_per_week * p_MMSMA
weekly_anemic_after <- infections_per_week * p_MMSMA_adjusted

#######################################################################################################
############################## Create a table of cases before and after vaccination
#######################################################################################################
cases_table <- data.frame(
  Week = weeks,
  Uncomplicated_Before = weekly_uncomplicated_before,
  Uncomplicated_After = weekly_uncomplicated_after,
  Cerebral_Before = weekly_cerebral_before,
  Cerebral_After = weekly_cerebral_after,
  Anemic_Before = weekly_anemic_before,
  Anemic_After = weekly_anemic_after
)
print(cases_table)
# Calculate the average cases before and after vaccination
average_uncomplicated_before <- mean(weekly_uncomplicated_before)
average_uncomplicated_after <- mean(weekly_uncomplicated_after)
average_cerebral_before <- mean(weekly_cerebral_before)
average_cerebral_after <- mean(weekly_cerebral_after)
average_anemic_before <- mean(weekly_anemic_before)
average_anemic_after <- mean(weekly_anemic_after)


#########################################################################################
########################### Create a summary table of average cases
#########################################################################################
average_cases_table <- data.frame(
  Type = c("Uncomplicated", "Cerebral", "Anemic"),
  Before_Vaccination = c(average_uncomplicated_before, average_cerebral_before, average_anemic_before),
  After_Vaccination = c(average_uncomplicated_after, average_cerebral_after, average_anemic_after)
)

# Print the average cases table
print(average_cases_table)


##########################################################################################
########################## Plot the trend of infections per week over time
##########################################################################################
plot(weeks, infections_per_week, type = "l", 
     xlab = "Week", ylab = "Number of infections", 
     main = "Trend of Infections per Week Over Time")

# Add lines for different types of malaria cases before vaccination
lines(weeks, weekly_uncomplicated_before, col = "blue")
#lines(weeks, weekly_cerebral_before, col = "pink")
#lines(weeks, weekly_anemic_before, col = "green")

# Add lines for different types of malaria cases after vaccination
lines(weeks, weekly_uncomplicated_after, col = "black", lty = 2)
#lines(weeks, weekly_cerebral_after, col = "red", lty = 2)
#lines(weeks, weekly_anemic_after, col = "darkgreen", lty = 2)

# Add a legend
legend("topright", legend = c("Uncomplicated Before", 
                              "Uncomplicated After" ),
       col = c("blue", "black"), lty = c(1, 2), cex = 0.8)



# Determine the range for the y-axis to include all lines
y_max <- max(c(weekly_anemic_before, weekly_cerebral_before, weekly_anemic_after, weekly_cerebral_after))
y_min <- min(c(weekly_anemic_before, weekly_cerebral_before, weekly_anemic_after, weekly_cerebral_after))


plot(weeks, weekly_anemic_before, type = "l", 
     xlab = "Week", ylab = "Number of Infections", 
     main = "Trend of Severe Infections per Week Over Time",
     col = "green", ylim = c(y_min, y_max))


#lines(weeks, weekly_uncomplicated_before, col = "blue")
lines(weeks, weekly_cerebral_before, col = "pink")
lines(weeks, weekly_anemic_before, col = "green")

# Add lines for different types of malaria cases after vaccination
#lines(weeks, weekly_uncomplicated_after, col = "black", lty = 2)
lines(weeks, weekly_cerebral_after, col = "red", lty = 2)
lines(weeks, weekly_anemic_after, col = "darkgreen", lty = 2)

# Add a legend
legend("topright", legend = c("Cerebral Before", "Anemic Before", 
                              "Cerebral After", "Anemic After"), 
       col = c( "pink", "green",  "red", "darkgreen"), lty = c(1, 1, 2, 2), cex = 0.8)

##########################################################################################
# Create a table of probabilities before and after vaccination############################
##########################################################################################
probabilities_table <- data.frame(
  Type = c("Mild Malaria", "Cerebral Malaria", "Anemic Malaria"),
  Probability_Before_Vaccination = c(p_WMM, p_MMCM, p_MMSMA),
  Probability_After_Vaccination = c(p_WMM_adjusted, p_MMCM_adjusted, p_MMSMA_adjusted)
)

# Print the probabilities table
print(probabilities_table)

##########################################################################################



######################################################################################
###########EVPI######################################################################
######################################################################################

## ---------------------------------------------
generate_psa_params <- function(n_sim= 1000, seed = 071818){
  set.seed(seed) # set a seed to be able to reproduce the same results
  
  mean_p_WMM <- 0.0001169687
  mean_p_MMCM <- 0.00160
  mean_p_MMSMA <- 0.00853
  mean_p_CMW <- 0.85322
  mean_p_SMAW <- 0.89122
  mean_p_D_natural <- (0.08 / 365) * 5
  
  mean_p_WMM_R21 <- mean_p_WMM * (1 - (cov_R21 * efficacy_R21_MM))
  mean_p_MMCM_R21 <- mean_p_MMCM * (1 - (cov_R21 * efficacy_R21_SM))
  mean_p_MMSMA_R21 <- mean_p_MMSMA * (1 - (cov_R21 * efficacy_R21_SM))
  mean_p_CMW_R21 <- mean_p_CMW
  mean_p_SMAW_R21 <- mean_p_SMAW
  
  mean_Duration_MM <- 0.014
  mean_Duration_CM <- 0.041
  mean_Duration_SMA <- 0.041
  
  mean_Duration_MM_R21 <- mean_Duration_MM
  mean_Duration_CM_R21 <- mean_Duration_CM
  mean_Duration_SMA_R21 <- mean_Duration_SMA
  
  # Calculate min and max values based on 10% increase/decrease
  min_p_WMM <- mean_p_WMM * 0.9
  max_p_WMM <- mean_p_WMM * 1.1
  min_p_MMCM <- mean_p_MMCM * 0.9
  max_p_MMCM <- mean_p_MMCM * 1.1
  min_p_MMSMA <- mean_p_MMSMA * 0.9
  max_p_MMSMA <- mean_p_MMSMA * 1.1
  min_p_CMW <- mean_p_CMW * 0.9
  max_p_CMW <- mean_p_CMW * 1.1
  min_p_SMAW <- mean_p_SMAW * 0.9
  max_p_SMAW <- mean_p_SMAW * 1.1
  min_p_D_natural <- mean_p_D_natural * 0.9
  max_p_D_natural <- mean_p_D_natural * 1.1
  
  min_p_WMM_R21 <- mean_p_WMM_R21 * 0.9
  max_p_WMM_R21 <- mean_p_WMM_R21 * 1.1
  min_p_MMCM_R21 <- mean_p_MMCM_R21 * 0.9
  max_p_MMCM_R21 <- mean_p_MMCM_R21 * 1.1
  min_p_MMSMA_R21 <- mean_p_MMSMA_R21 * 0.9
  max_p_MMSMA_R21 <- mean_p_MMSMA_R21 * 1.1
  min_p_CMW_R21 <- mean_p_CMW_R21 * 0.9
  max_p_CMW_R21 <- mean_p_CMW_R21 * 1.1
  min_p_SMAW_R21 <- mean_p_SMAW_R21 * 0.9
  max_p_SMAW_R21 <- mean_p_SMAW_R21 * 1.1
  
  # Calculate min and max values based on 10% increase/decrease for durations
  min_Duration_MM <- mean_Duration_MM * 0.9
  max_Duration_MM <- mean_Duration_MM * 1.1
  min_Duration_CM <- mean_Duration_CM * 0.9
  max_Duration_CM <- mean_Duration_CM * 1.1
  min_Duration_SMA <- mean_Duration_SMA * 0.9
  max_Duration_SMA <- mean_Duration_SMA * 1.1
  
  min_Duration_MM_R21 <- mean_Duration_MM_R21 * 0.9
  max_Duration_MM_R21 <- mean_Duration_MM_R21 * 1.1
  min_Duration_CM_R21 <- mean_Duration_CM_R21 * 0.9
  max_Duration_CM_R21 <- mean_Duration_CM_R21 * 1.1
  min_Duration_SMA_R21 <- mean_Duration_SMA_R21 * 0.9
  max_Duration_SMA_R21 <- mean_Duration_SMA_R21 * 1.1
  
  
  df_psa <- data.frame(
    
    p_WMM = runif(n_sim, min = min_p_WMM, max = max_p_WMM),
    p_MMCM = runif(n_sim, min = min_p_MMCM, max = max_p_MMCM),
    p_MMSMA = runif(n_sim, min = min_p_MMSMA, max = max_p_MMSMA),
    p_CMW = runif(n_sim, min = min_p_CMW, max = max_p_CMW),
    p_SMAW = runif(n_sim, min = min_p_SMAW, max = max_p_SMAW),
    p_D_natural = runif(n_sim, min = min_p_D_natural, max = max_p_D_natural),
    
    p_WMM_R21 = runif(n_sim, min = min_p_WMM_R21, max = max_p_WMM_R21),
    p_MMCM_R21 = runif(n_sim, min = min_p_MMCM_R21, max = max_p_MMCM_R21),
    p_MMSMA_R21 = runif(n_sim, min = min_p_MMSMA_R21, max = max_p_MMSMA_R21),
    p_CMW_R21 = runif(n_sim, min = min_p_CMW_R21, max = max_p_CMW_R21),
    p_SMAW_R21 = runif(n_sim, min = min_p_SMAW_R21, max = max_p_SMAW_R21),
    
    Duration_MM = runif(n_sim, min = min_Duration_MM, max = max_Duration_MM),
    Duration_CM = runif(n_sim, min = min_Duration_CM, max = max_Duration_CM),
    Duration_SMA = runif(n_sim, min = min_Duration_SMA, max = max_Duration_SMA),
    
    Duration_MM_R21 = runif(n_sim, min = min_Duration_MM_R21, max = max_Duration_MM_R21),
    Duration_CM_R21 = runif(n_sim, min = min_Duration_CM_R21, max = max_Duration_CM_R21),
    Duration_SMA_R21 = runif(n_sim, min = min_Duration_SMA_R21, max = max_Duration_SMA_R21),
  
    
    # Costs
    c_consult_HF     = rgamma (n_sim, shape = 114.3138, scale = 0.0291),    
    c_consult_CF     = rgamma (n_sim, shape = 1.5864, scale = 0.5042), 
    diagnost_HF_RDT  = rgamma (n_sim, shape = 50.8907, scale = 0.0499), 
    Diagnostic_HF_M  = rgamma (n_sim, shape = 2.3824, scale = 0.3844), 
    Diagnostic_CHW_RDT = rgamma (n_sim, shape = 50.8907, scale = 0.0499), 
    HF_treatconsult    = rgamma (n_sim, shape = 114.3138, scale = 0.0291), 
    CHW_treatconsult   = rgamma (n_sim, shape = 1.5864, scale = 0.5042), 
    Follow_up_con      = rgamma (n_sim, shape = 2263.4629, scale = 0.0037), 
    Anti_mal_treat_MM  = rgamma (n_sim, shape = 3.5457, scale = 0.0531), 
    Anti_mal_treat_SM  = rgamma (n_sim, shape = 6440.1861, scale = 0.0012), 
    CM_treat    = rgamma (n_sim, shape = 222.0100, scale = 0.0067), 
    SMA_treat   = rgamma (n_sim, shape = 7643.1306, scale = 0.0011), 
    Trans_pat   = rgamma (n_sim, shape = 5.6457, scale = 0.4924), 
    Trans_acomp = rgamma (n_sim, shape = 12.9070, scale = 0.0640), 
    wage_loss   = rgamma (n_sim, shape = 13386.4900, scale = 0.0009), 
    time_loss_work_MM = rgamma (n_sim, shape = 4761.0000, scale = 0.0014), 
    time_loss_work_SM = rgamma (n_sim, shape = 10000.0000, scale = 0.0010), 
    LLIN_unit   = rgamma (n_sim, shape = 46.2400, scale = 0.0147), 
    R21_MatrixM = rgamma (n_sim, shape = 14400.0000, scale = 0.0008),
    cov_consult_HF  = rbeta (n_sim, shape1 =  11.7113, shape2 = 2.3987), 
    cov_consult_CF  = rbeta (n_sim, shape1 = 2.1504, shape2 = 11.2896), 
    cov_diagnostic_HF_RDT  = rbeta (n_sim, shape1 = 2.3987, shape2 = 11.7113), 
    cov_diagnostic_HF_M    = rbeta (n_sim, shape1 = 14.8104, shape2 = 7.6296), 
    cov_diagnostic_CHW_RDT = rbeta (n_sim, shape1 = 1.6856, shape2 = 10.3544), 
    cov_HFtreatconsult   = rbeta (n_sim, shape1 = 11.7113, shape2 = 2.3987), 
    cov_CHW_treatconsult = rbeta (n_sim, shape1 = 2.1504, shape2 = 11.2896), 
    cov_LLIN = rbeta (n_sim, shape1 = 8.6247, shape2 = 14.6853), 
    cov_R21  = rbeta (n_sim, shape1 = 3.9521, shape2 = 0.8675), 
    
    # DALYS
    W_MM =  rbeta (n_sim, shape1 = 77.6151, shape2 = 328.7465), 
    W_CM =  rbeta (n_sim, shape1 = 28.6509, shape2 = 32.1790), 
    W_SMA = rbeta (n_sim, shape1 = 28.609, shape2 =32.1790), 
    
    N_cases_MM =  rgamma (n_sim, shape = 64160100.000, scale = 0.000012),   
    N_cases_CM =  rgamma (n_sim, shape = 776337.210, scale = 0.000113),   
    N_cases_SMA = rgamma (n_sim, shape = 1084305.690, scale = 0.000096),   
    
    W_MM_R21  = rbeta (n_sim, shape1 = 77.6151, shape2 = 328.7465), 
    W_CM_R21  = rbeta (n_sim, shape1 = 28.6509, shape2 = 32.1790), 
    W_SMA_R21 = rbeta (n_sim, shape1 = 28.6509, shape2 = 32.1790), 
    
    N_cases_MM_R21  = rgamma (n_sim, shape = 26936100.000, scale = 0.000019),   
    N_cases_CM_R21  = rgamma (n_sim, shape = 325926.810, scale = 0.000175),   
    N_cases_SMA_R21 = rgamma (n_sim, shape = 455220.090, scale = 0.000148)   
    
  )
  return(df_psa)
}


l_params_all <- list( # Transition probabilities (per cycle), hazard ratios
  c_consult_HF <- 3.33 ,
  c_consult_CF <- 0.80,
  diagnost_HF_RDT <- 2.54,
  Diagnostic_HF_M <- 0.92,
  Diagnostic_CHW_RDT <- 2.54,
  Diagnostic_CHW_M <- 0,
  HF_treatconsult <- 3.33,
  CHW_treatconsult <- 0.80 ,
  Follow_up_con <- 8.27,
  Anti_mal_treat_MM <- 0.188,
  Anti_mal_treat_SM <- 8.0251,
  CM_treat <- 1.49,
  SMA_treat <- 8.74,
  Trans_pat <- 3.59,
  Trans_acomp <- 0.83,
  wage_loss <-11.57,
  time_loss_work_MM <- 6.9,
  time_loss_work_SM <- 10,
  LLIN_unit <- 0.68,
  R21_MatrixM <- 12,
  
  LLIN_duration <- 5,
  cov_consult_HF <- 0.83,
  cov_consult_CF <- 0.16,
  cov_diagnostic_HF_RDT <-0.17,
  cov_diagnostic_HF_M <- 0.66,
  cov_diagnostic_CHW_RDT <- 0.14,
  cov_HFtreatconsult <- 0.83,
  cov_CHW_treatconsult <- 0.16,
  cov_LLIN <- 0.37,
  cov_R21 <-  0.82,
  
  efficacy_R21_MM <- 0.73 ,
  efficacy_R21_SM <- 0.67 ,
  efficacy_LLIN <- 0.55,
  
  ### DALYs ---- overall disease burden 
  W_MM <- 0.19,
  W_CM <-0.47,
  W_SMA <-0.47,
  Duration_MM <-0.014,
  Duration_CM <-0.041,
  Duration_SMA <-0.041,
  
  N_cases_MM <-801,
  N_cases_CM <-88.11,
  N_cases_SMA <-104.13,
  
  W_MM_R21 <-0.19,
  W_CM_R21 <-0.47,
  W_SMA_R21 <-0.47,
  Duration_MM_R21 <-0.014,
  Duration_CM_R21 <-0.041,
  Duration_SMA_R21 <-0.041,
  N_cases_MM_R21 <- 519,
  N_cases_CM_R21 <-57,
  N_cases_SMA_R21 <-67,
  
  
  p_WMM  <- 0.0001169687, # constant annual probability of getting mild malaria when Well conditional on surviving 
  p_MMCM  <- 0.00160,
  p_MMSMA  <- 0.00853, 
  p_CMW <- 0.85322, ##also use malawi text 
  p_SMAW <- 0.89122, ##also use malawi text 
  
  p_WMM_R21 <- p_WMM *(1-(cov_R21*efficacy_R21_MM)),
  p_MMCM_R21 <- p_MMCM *(1-(cov_R21*efficacy_R21_SM)),
  p_MMSMA_R21 <-p_MMSMA *(1-(cov_R21*efficacy_R21_SM)),
  #p_WMM_R21 <- 0.0000469397 #Adjusted Probability=Original Probability×(1−(Efficacy×Coverage))
  #p_MMCM_R21 <- 0.00072096  #Adjusted Probability=Original Probability×(1−(Efficacy×Coverage))
  #p_MMSMA_R21 <- 0.00384362   #Adjusted Probability=Original Probability×(1−(Efficacy×Coverage))
  
  p_CMW_R21 <- 0.85322, ##also use malawi text 
  p_SMAW_R21 <- 0.89122, ##also use malawi text  the going back to health stays the same because the vaccine helps not to get sick but not to recover faster 
  p_D_natural <- (0.08 /365)*5, ## from the IHN page death fro other reasons 
  
  d_c <- 0.03, # annual discount rate for costs 
  d_e <- 0.03, # annual discount rate for DALYs
  
  cycle_length <- 5/365, # cycle length 5 Days  
  n_age_init <- 1, # age at baseline --> maybe start at age where children are completely vaccinated ? which is at 1 year and 7 months 
  n_age_max  <- 5, # maximum age of follow up
  
  #Number of cycles for each tunnel state
  n_tunnel_size_CM <- 2,
  n_tunnel_size_SMA <- 2,
  
  DALY_W   <- 0, 
  DALY_W_R21  <- 0  # you are healthy so 0 represents perfect health 
)


#* Store the parameter names into a vector
v_names_params <- names(l_params_all)

decision_model <- function(l_params_all, verbose = FALSE) {
  with(as.list(l_params_all), {
    ########################### Process model inputs ###########################
    ### Process model inputs
    ## Number of cycles
    
    n_cycles <- (n_age_max - n_age_init)/cycle_length # time horizon, number of cycles
    ## Age labels 
    v_age_names <- paste(rep(n_age_init:(n_age_max-1), each = 1/cycle_length), 
                         1:(1/cycle_length), 
                         sep = ".")
    
    #costs
    c_W <- ((LLIN_unit)*cov_LLIN)
    
    c_MM <- (c_consult_HF*cov_consult_HF) + (c_consult_CF*cov_consult_CF) + (diagnost_HF_RDT*cov_diagnostic_HF_RDT) +
      (Diagnostic_HF_M*cov_diagnostic_HF_M) + (Diagnostic_CHW_RDT*cov_diagnostic_CHW_RDT) +Diagnostic_CHW_M + 
      (HF_treatconsult*cov_HFtreatconsult)+ (CHW_treatconsult*cov_CHW_treatconsult )+ Follow_up_con+ 
      Anti_mal_treat_MM+ Trans_pat+ Trans_acomp+(wage_loss*time_loss_work_MM)       
    
    c_SMA <- (c_consult_HF*cov_consult_HF) + (c_consult_CF*cov_consult_CF) + (diagnost_HF_RDT*cov_diagnostic_HF_RDT) + 
      (Diagnostic_HF_M*cov_diagnostic_HF_M) + (Diagnostic_CHW_RDT*cov_diagnostic_CHW_RDT)+Diagnostic_CHW_M + 
      (HF_treatconsult*cov_HFtreatconsult)+ (CHW_treatconsult*cov_CHW_treatconsult )+ Follow_up_con+ Anti_mal_treat_SM+ 
      SMA_treat+ Trans_pat+ Trans_acomp 
    
    
    c_CM <- (c_consult_HF*cov_consult_HF) + (c_consult_CF*cov_consult_CF) + (diagnost_HF_RDT*cov_diagnostic_HF_RDT) + 
      (Diagnostic_HF_M*cov_diagnostic_HF_M) + (Diagnostic_CHW_RDT*cov_diagnostic_CHW_RDT)+Diagnostic_CHW_M +
      (HF_treatconsult*cov_HFtreatconsult)+ (CHW_treatconsult*cov_CHW_treatconsult )+ Follow_up_con+ Anti_mal_treat_SM+
      CM_treat+ Trans_pat+ Trans_acomp 
    
    c_W_R21 <-  ((LLIN_unit)*cov_LLIN) + (R21_MatrixM*cov_R21)
    
    c_MM_R21 <- (c_consult_HF*cov_consult_HF) + (c_consult_CF*cov_consult_CF) + (diagnost_HF_RDT*cov_diagnostic_HF_RDT) +
      (Diagnostic_HF_M*cov_diagnostic_HF_M) + (Diagnostic_CHW_RDT*cov_diagnostic_CHW_RDT)+Diagnostic_CHW_M + 
      (HF_treatconsult*cov_HFtreatconsult)+ (CHW_treatconsult*cov_CHW_treatconsult)+ Follow_up_con+ Anti_mal_treat_MM+ 
      Trans_pat+ Trans_acomp+(wage_loss*time_loss_work_MM)
    
    c_SMA_R21 <- (c_consult_HF*cov_consult_HF) + (c_consult_CF*cov_consult_CF) + (diagnost_HF_RDT*cov_diagnostic_HF_RDT) + 
      (Diagnostic_HF_M*cov_diagnostic_HF_M) + (Diagnostic_CHW_RDT*cov_diagnostic_CHW_RDT)+Diagnostic_CHW_M + 
      (HF_treatconsult*cov_HFtreatconsult)+ (CHW_treatconsult*cov_CHW_treatconsult )+ Follow_up_con+ Anti_mal_treat_SM+ 
      SMA_treat+ Trans_pat+ Trans_acomp
    
    c_CM_R21 <- (c_consult_HF*cov_consult_HF) + (c_consult_CF*cov_consult_CF) + (diagnost_HF_RDT*cov_diagnostic_HF_RDT) + 
      (Diagnostic_HF_M*cov_diagnostic_HF_M) + (Diagnostic_CHW_RDT*cov_diagnostic_CHW_RDT)+Diagnostic_CHW_M +
      (HF_treatconsult*cov_HFtreatconsult)+ (CHW_treatconsult*cov_CHW_treatconsult )+ Follow_up_con+ Anti_mal_treat_SM+ 
      CM_treat+ Trans_pat+ Trans_acomp
    
    
    DALY_MM   <- W_MM * Duration_MM * N_cases_MM  # annual DALY of having moderate malaria
    DALY_CM   <- W_CM * Duration_CM * N_cases_CM   # annual DALY of having cerebral malaria
    DALY_SMA   <- W_SMA *Duration_SMA * N_cases_SMA    # annual DALY of being severe malaria anemia   
    
    DALY_MM_R21   <- W_MM_R21 * Duration_MM_R21 * N_cases_MM_R21# annual DALY of having moderate malaria
    DALY_CM_R21   <- W_CM_R21 * Duration_CM_R21 * N_cases_CM_R21# annual DALY of having cerebral malaria
    DALY_SMA_R21   <- W_SMA_R21 * Duration_SMA_R21*  N_cases_SMA_R21   # annual DALY of being severe malaria anemia   
    
    
    
    ##################### Construct state-transition models ####################
    ## Initial state vector
    # All starting healthy
    ## Initial state vector ----
    
    ## Health states of the model:
    v_names_states <- c("W",   #well 
                        "MM",  # Mild Malaria (MM)
                        "CM",  # Cerebral Malaria (CM)
                        "SMA", # Severe Malaria Anemia (SMA)  
                        "D")   # Death other causes  
    n_states <- length(v_names_states)   # number of health states 
    
    
    #### Vector with cycles for cerebral malaria tunnel
    v_cycles_tunnel_CM <- 1:n_tunnel_size_CM
    #### Vector with names for tunnel states of cerebral malaria
    v_CM_tunnel <- paste("CM", seq(1, n_tunnel_size_CM), "", sep = "")
    
    #### Vector with cycles for severe anemic malaria tunnel
    v_cycles_tunnel_SMA <- 1:n_tunnel_size_SMA
    #### Vector with names for tunnel states of severe anemic malaria
    v_SMA_tunnel <- paste("SMA", seq(1, n_tunnel_size_SMA), "", sep = "")
    
    #### Create variables for model with tunnels
    v_names_states_tunnels <- c("W", "MM", v_CM_tunnel, v_SMA_tunnel, "D") # health state names
    n_states_tunnels <- length(v_names_states_tunnels)   
    
    
    ## Strategies ----
    
    v_names_str <- c("LLIN", 
                     "R21") 
    n_str       <- length(v_names_str)        # number of strategies
    
    
    gen_wcc <- function(n_cycles, method = "half-cycle") {
      if (method == "half-cycle") {
        v_wcc <- rep(1, n_cycles + 1)
        v_wcc[1] <- v_wcc[n_cycles + 1] <- 0.5
        return(v_wcc)  # Return v_wcc
      } else {
        stop("Invalid method. Please choose 'half-cycle'.")
      }
    }
    
    
    ##Discount weight for costs and effects ---
    v_dwc  <- 1 / ((1 + (d_e * cycle_length)) ^ (0:n_cycles))
    v_dwe  <- 1 / ((1 + (d_c * cycle_length)) ^ (0:n_cycles))
    
    # Construct state-transition models ---- see table already
    ## Initial state vector ----
    v_m_init_tunnels <- c(1000, 0, rep(0, n_tunnel_size_CM + n_tunnel_size_SMA), 0) 
    
    
    ## Initialize cohort traces ----
    ### Initialize cohort trace for state-residece dependent cSTM under SoC----
    m_M_tunnels_LLIN <- matrix(0, 
                               nrow   = (n_cycles + 1), ncol = n_states_tunnels, 
                               dimnames = list(0:n_cycles, v_names_states_tunnels))
    
    #* Store the initial state vector in the first row of the cohort trace
    m_M_tunnels_LLIN[1, ] <- v_m_init_tunnels
    
    ### Initialize cohort trace for strategies R21
    #* Structure and initial states are the same as for SoC
    m_M_tunnels_R21  <- m_M_tunnels_LLIN # Strategy R21
    
    
    ## Create transition probability arrays for strategy no intervention ----
    ### Initialize transition probability array for strategy SoC ----
    ##'s a three-dimensional array with dimensions matching the number of health states (n_states_tunnels) and cycles (n_cycles). The dimnames argument provides labels for the dimensions, which is good practice for clarity and readability.
    #* All transitions to a non-death state are assumed to be conditional on survival
    a_P_tunnels_LLIN <- array(0, 
                              dim = c(n_states_tunnels, n_states_tunnels, n_cycles),
                              dimnames = list(v_names_states_tunnels, 
                                              v_names_states_tunnels, 
                                              0:(n_cycles-1)))
    
    
    
    ###################################################################
    ###  Initialize transition probability array for strategy in array LLIN
    
    ## From W 
    a_P_tunnels_LLIN ["W", "W", ]               <-((1-p_D_natural) *(1- p_WMM))
    a_P_tunnels_LLIN ["W", "MM", ]              <-(1-p_D_natural)* p_WMM 
    a_P_tunnels_LLIN ["W", v_CM_tunnel, ]       <- 0
    a_P_tunnels_LLIN ["W", v_SMA_tunnel, ]      <- 0
    a_P_tunnels_LLIN ["W", "D", ]               <- p_D_natural
    
    ## From MM
    a_P_tunnels_LLIN["MM", "W", ]               <- ((1-p_D_natural)*(1-(p_MMCM+p_MMSMA)))
    a_P_tunnels_LLIN["MM", "MM", ]              <- 0
    a_P_tunnels_LLIN["MM", v_CM_tunnel[1], ]    <- (1-p_D_natural)*p_MMCM
    a_P_tunnels_LLIN["MM", v_SMA_tunnel[1], ]   <- (1-p_D_natural) *p_MMSMA
    a_P_tunnels_LLIN["MM", v_CM_tunnel[2], ]    <- 0  # Transition probability to the second tunnel set to 0
    a_P_tunnels_LLIN["MM", v_SMA_tunnel[2], ]   <- 0  # Transition probability to the second tunnel set to 0
    a_P_tunnels_LLIN ["MM", "D", ]               <- p_D_natural
    
    
    ## From CM
    # Transition probability within the first CM tunnel NOT NEEDED BECAUSE IT'S A TUNNEL; YOU ALWAYS MOVE.
    # Transition probability from the first CM tunnel to the "well" state
    a_P_tunnels_LLIN[v_CM_tunnel[1], "W", ]                   <- (1-p_D_natural)* p_CMW
    a_P_tunnels_LLIN[v_CM_tunnel[1], "MM", ]                   <- 0
    a_P_tunnels_LLIN[v_CM_tunnel[1], v_CM_tunnel[2], ]        <-  (1-p_D_natural)*(1-p_CMW) 
    a_P_tunnels_LLIN[v_CM_tunnel[1], v_CM_tunnel[1], ]        <- 0
    a_P_tunnels_LLIN[v_CM_tunnel[1], v_SMA_tunnel[1], ]        <- 0
    a_P_tunnels_LLIN[v_CM_tunnel[1], v_SMA_tunnel[2], ]        <- 0
    a_P_tunnels_LLIN[v_CM_tunnel[1], "D", ]        <- p_D_natural
    
    # Transition probability from the second CM tunnel to the "well" state
    a_P_tunnels_LLIN[v_CM_tunnel[2], "W", ]                   <- (1-p_D_natural)
    a_P_tunnels_LLIN[v_CM_tunnel[2], "MM", ]                   <- 0
    a_P_tunnels_LLIN[v_CM_tunnel[2], v_CM_tunnel[2], ]        <- 0
    a_P_tunnels_LLIN[v_CM_tunnel[2], v_CM_tunnel[1], ]          <- 0
    a_P_tunnels_LLIN[v_CM_tunnel[2], v_SMA_tunnel[1], ]          <- 0
    a_P_tunnels_LLIN[v_CM_tunnel[2], v_SMA_tunnel[2], ]          <- 0
    a_P_tunnels_LLIN[v_CM_tunnel[2], "D", ]        <- p_D_natural
    
    
    ## From SMA
    # Transition probability within the first CM tunnel NOT NEEDED BECAUSE IT'S A TUNNEL; YOU ALWAYS MOVE.
    # Transition probability from the first CM tunnel to the "well" state
    a_P_tunnels_LLIN[v_SMA_tunnel[1], "W", ]                <- (1-p_D_natural)*p_SMAW
    a_P_tunnels_LLIN[v_SMA_tunnel[1], "MM", ]               <- 0
    a_P_tunnels_LLIN[v_SMA_tunnel[1], v_SMA_tunnel[2], ]    <- (1-p_D_natural)* (1-p_SMAW ) 
    a_P_tunnels_LLIN[v_SMA_tunnel[1], v_SMA_tunnel[1], ]    <- 0 
    a_P_tunnels_LLIN[v_SMA_tunnel[1], v_CM_tunnel[1], ]     <- 0
    a_P_tunnels_LLIN[v_SMA_tunnel[1], v_CM_tunnel[2], ]     <- 0
    a_P_tunnels_LLIN[v_SMA_tunnel[1], "D", ]        <- p_D_natural
    
    # Transition probability from the second CM tunnel to the "well" state
    a_P_tunnels_LLIN[v_SMA_tunnel[2], "W", ]         <-  (1-p_D_natural) 
    a_P_tunnels_LLIN[v_SMA_tunnel[2], "MM", ]        <- 0
    a_P_tunnels_LLIN[v_SMA_tunnel[2], v_SMA_tunnel[2], ]   <- 0
    a_P_tunnels_LLIN[v_SMA_tunnel[2], v_SMA_tunnel[1], ] <- 0 
    a_P_tunnels_LLIN[v_SMA_tunnel[2], v_CM_tunnel[1], ]    <- 0
    a_P_tunnels_LLIN[v_SMA_tunnel[2], v_CM_tunnel[2], ]    <- 0
    a_P_tunnels_LLIN[v_SMA_tunnel[2], "D", ]        <- p_D_natural
    
    ## From D
    a_P_tunnels_LLIN["D", "W", ]              <- 0
    a_P_tunnels_LLIN["D", "MM", ]             <- 0
    a_P_tunnels_LLIN["D", "D", ]              <- 1
    a_P_tunnels_LLIN["D", v_CM_tunnel[2], ]   <- 0
    a_P_tunnels_LLIN["D", v_CM_tunnel[1], ]   <- 0
    a_P_tunnels_LLIN["D", v_SMA_tunnel[1], ]  <- 0
    a_P_tunnels_LLIN["D", v_SMA_tunnel[2], ]  <- 0
    
    
    print("Transition Probabilities:")
    print(a_P_tunnels_LLIN)
    
    ## Check that transition probabilities are [0, 1]
    check_transition_probability(a_P_tunnels_LLIN, verbose = TRUE)
    # Check the sum of transition probabilities
    check_sum_of_transition_array(a_P_tunnels_LLIN, n_states = n_states_tunnels, n_cycles = n_cycles, verbose = TRUE)
    
    
    ###################################################################
    
    # Initialize transition probability array for strategy R21  ----
    a_P_tunnels_R21 <- a_P_tunnels_LLIN
    ###################################################################
    ###  Initialize transition probability array for strategy in array LLIN
    ## From W 
    a_P_tunnels_R21 ["W", "W", ]               <- (1-p_D_natural)*(1 - p_WMM_R21) 
    a_P_tunnels_R21 ["W", "MM", ]              <- (1-p_D_natural)*p_WMM_R21 
    a_P_tunnels_R21 ["W", v_CM_tunnel, ]       <- 0
    a_P_tunnels_R21 ["W", v_SMA_tunnel, ]      <- 0
    a_P_tunnels_R21 ["W", "D", ]               <- p_D_natural
    
    
    ## From MM
    a_P_tunnels_R21["MM", "W", ]               <- (1-p_D_natural)*(1-(p_MMCM_R21+p_MMSMA_R21))
    a_P_tunnels_R21["MM", "MM", ]              <- 0
    a_P_tunnels_R21["MM", v_CM_tunnel[1], ]    <-  (1-p_D_natural)*p_MMCM_R21
    a_P_tunnels_R21["MM", v_SMA_tunnel[1], ]   <-  (1-p_D_natural)*p_MMSMA_R21
    a_P_tunnels_R21["MM", v_CM_tunnel[2], ]    <- 0  # Transition probability to the second tunnel set to 0
    a_P_tunnels_R21["MM", v_SMA_tunnel[2], ]   <- 0  # Transition probability to the second tunnel set to 0
    a_P_tunnels_R21 ["MM", "D", ]              <- p_D_natural
    
    ## From CM
    # Transition probability within the first CM tunnel NOT NEEDED BECAUSE IT'S A TUNNEL; YOU ALWAYS MOVE.
    # Transition probability from the first CM tunnel to the "well" state
    a_P_tunnels_R21[v_CM_tunnel[1], "W", ]                   <- (1-p_D_natural)*  p_CMW_R21
    a_P_tunnels_R21[v_CM_tunnel[1], "MM", ]                   <- 0
    a_P_tunnels_R21[v_CM_tunnel[1], v_CM_tunnel[2], ]        <-  (1-p_D_natural)* (1-p_CMW_R21) 
    a_P_tunnels_R21[v_CM_tunnel[1], v_CM_tunnel[1], ]        <- 0
    a_P_tunnels_R21[v_CM_tunnel[1], v_SMA_tunnel[1], ]        <- 0
    a_P_tunnels_R21[v_CM_tunnel[1], v_SMA_tunnel[2], ]        <- 0
    a_P_tunnels_R21[v_CM_tunnel[1], "D", ]        <- p_D_natural
    
    # Transition probability from the second CM tunnel to the "well" state
    a_P_tunnels_R21[v_CM_tunnel[2], "W", ]                   <- (1-p_D_natural)* 1
    a_P_tunnels_R21[v_CM_tunnel[2], "MM", ]                   <- 0
    a_P_tunnels_R21[v_CM_tunnel[2], v_CM_tunnel[2], ]        <- 0
    a_P_tunnels_R21[v_CM_tunnel[2], v_CM_tunnel[1], ]          <- 0
    a_P_tunnels_R21[v_CM_tunnel[2], v_SMA_tunnel[1], ]          <- 0
    a_P_tunnels_R21[v_CM_tunnel[2], v_SMA_tunnel[2], ]          <- 0
    a_P_tunnels_R21[v_CM_tunnel[2], "D", ]        <- p_D_natural
    
    ## From SMA
    # Transition probability within the first CM tunnel NOT NEEDED BECAUSE IT'S A TUNNEL; YOU ALWAYS MOVE.
    # Transition probability from the first CM tunnel to the "well" state
    a_P_tunnels_R21[v_SMA_tunnel[1], "W", ]                <-  (1-p_D_natural)* p_SMAW_R21
    a_P_tunnels_R21[v_SMA_tunnel[1], "MM", ]               <- 0
    a_P_tunnels_R21[v_SMA_tunnel[1], v_SMA_tunnel[2], ]    <- (1-p_D_natural)* (1-p_SMAW_R21 ) 
    a_P_tunnels_R21[v_SMA_tunnel[1], v_SMA_tunnel[1], ]    <- 0 
    a_P_tunnels_R21[v_SMA_tunnel[1], v_CM_tunnel[1], ]     <- 0
    a_P_tunnels_R21[v_SMA_tunnel[1], v_CM_tunnel[2], ]     <- 0
    a_P_tunnels_R21[v_SMA_tunnel[1], "D", ]        <- p_D_natural
    
    # Transition probability from the second CM tunnel to the "well" state
    a_P_tunnels_R21[v_SMA_tunnel[2], "W", ]         <- (1-p_D_natural)
    a_P_tunnels_R21[v_SMA_tunnel[2], "MM", ]        <- 0
    a_P_tunnels_R21[v_SMA_tunnel[2], v_SMA_tunnel[2], ]   <-  0
    a_P_tunnels_R21[v_SMA_tunnel[2], v_SMA_tunnel[1], ] <- 0 
    a_P_tunnels_R21[v_SMA_tunnel[2], v_CM_tunnel[1], ]    <- 0
    a_P_tunnels_R21[v_SMA_tunnel[2], v_CM_tunnel[2], ]    <- 0
    a_P_tunnels_R21[v_SMA_tunnel[2], "D", ]        <- p_D_natural
    
    ## From D
    a_P_tunnels_R21["D", "W", ]              <- 0
    a_P_tunnels_R21["D", "MM", ]             <- 0
    a_P_tunnels_R21["D", "D", ]              <- 1
    a_P_tunnels_R21["D", v_CM_tunnel[2], ]   <- 0
    a_P_tunnels_R21["D", v_CM_tunnel[1], ]   <- 0
    a_P_tunnels_R21["D", v_SMA_tunnel[1], ]  <- 0
    a_P_tunnels_R21["D", v_SMA_tunnel[2], ]  <- 0
    
    
    print("Transition Probabilities:")
    print(a_P_tunnels_R21)
    
    ## Check that transition probabilities are [0, 1]
    check_transition_probability(a_P_tunnels_R21, verbose = TRUE)
    # Check the sum of transition probabilities
    check_sum_of_transition_array(a_P_tunnels_R21, n_states = n_states_tunnels, n_cycles = n_cycles, verbose = TRUE)
    
    
    # Create transition dynamics arrays ----
    
    #* These arrays will capture transitions from each state to another over time 
    ### Initialize transition dynamics array for strategy SoC ----
    a_A_tunnels_LLIN <- array(0,
                              dim = c(n_states_tunnels, n_states_tunnels, n_cycles + 1),
                              dimnames = list(v_names_states_tunnels, v_names_states_tunnels, 0:n_cycles))
    #* Set first slice of Array with the initial state vector in its diagonal
    diag(a_A_tunnels_LLIN[, , 1]) <- v_m_init_tunnels
    ### Initialize transition-dynamics array  ----
    #* Structure and initial states are the same as for SoC
    a_A_tunnels_R21  <- a_A_tunnels_LLIN
    
    
    #  Run Markov model ----updates cohort traces and transition-dynamics arrays for multiple interventions in 
    # a Markov model, calculating the next state probabilities based on transition probabilities and the current
    # state of the cohort
    #* Iterative solution of state-residence dependent 
    #* cycle length , states dimes.
    for (t in 1:n_cycles) {
      ## Fill in cohort trace
      # For no intervention 
      m_M_tunnels_LLIN[t + 1, ]   <- m_M_tunnels_LLIN[t, ] %*% a_P_tunnels_LLIN [, , t]
      # For R21
      m_M_tunnels_R21[t + 1, ]      <- m_M_tunnels_R21[t, ] %*% a_P_tunnels_R21 [, , t]
      
      ## Fill in transition-dynamics array
      # For no intervention
      a_A_tunnels_LLIN[, , t + 1]   <- diag(m_M_tunnels_LLIN[t, ]) %*% a_P_tunnels_LLIN [, , t]
      # For R21
      a_A_tunnels_R21[, , t + 1]   <- diag(m_M_tunnels_R21[t, ]) %*% a_P_tunnels_R21[, , t]
    }
    
    m_M_tunnels_LLIN_sum <- cbind( W   = m_M_tunnels_LLIN[, "W"], 
                                   MM  = m_M_tunnels_LLIN[, "MM"], 
                                   CM  = rowSums ( m_M_tunnels_LLIN [, 3:4]), 
                                   SMA = rowSums ( m_M_tunnels_LLIN [, 5:6]),
                                   D= m_M_tunnels_LLIN[, "D"])
    
    m_M_tunnels_R21_sum <- cbind( W   = m_M_tunnels_R21[, "W"], 
                                  MM  = m_M_tunnels_R21[, "MM"], 
                                  CM  = rowSums(m_M_tunnels_R21[, 3:4]), 
                                  SMA = rowSums(m_M_tunnels_R21[, 5:6]),
                                  D= m_M_tunnels_R21[, "D"])
    
    ## Store the cohort traces in a list ----
    l_m_M <- list(m_M_tunnels_LLIN_sum,      
                  m_M_tunnels_R21_sum)
    names(l_m_M) <- v_names_str
    
    ## Store the transition dynamics array for each strategy in a list ----
    l_a_A <- list(a_A_tunnels_LLIN,
                  a_A_tunnels_R21)
    names(l_a_A) <- v_names_str
    
    
    ########################################## RETURN OUTPUT  ##########################################
    out <- list(n_cycles = n_cycles,
                n_tunnel_size = n_cycles,
                l_m_M = l_m_M,
                l_a_A = l_a_A)
    
    return(out)
  }
  )
} 

calculate_ce_out <- function(l_params_all, n_wtp = 1000){ # User defined
  with(as.list(l_params_all), {
    
    n_cycles <- (n_age_max - n_age_init)/cycle_length # time horizon, number of cycles
    ## Age labels 
    v_age_names <- paste(rep(n_age_init:(n_age_max-1), each = 1/cycle_length), 
                         1:(1/cycle_length), 
                         sep = ".")
    
    
    #### Vector with cycles for cerebral malaria tunnel
    v_cycles_tunnel_CM <- 1:n_tunnel_size_CM
    
    #### Vector with names for tunnel states of cerebral malaria
    v_CM_tunnel <- paste("CM", seq(1, n_tunnel_size_CM), "", sep = "")
    
    #### Vector with cycles for severe anemic malaria tunnel
    v_cycles_tunnel_SMA <- 1:n_tunnel_size_SMA
    
    #### Vector with names for tunnel states of severe anemic malaria
    v_SMA_tunnel <- paste("SMA", seq(1, n_tunnel_size_SMA), "", sep = "")
    
    #### Create variables for model with tunnels
    v_names_states_tunnels <- c("W", "MM", v_CM_tunnel, v_SMA_tunnel, "D") # health state names
    
    
    model <- decision_model(l_params_all = l_params_all)
    l_a_A <- model[["l_a_A"]]
    
    
    # Vector of DALY per cycle for SMA under strategy No Intervention
    v_DALY_SMA_LLIN <- rep(DALY_SMA, n_tunnel_size_SMA)
    names(v_DALY_SMA_LLIN) <- v_SMA_tunnel[-1]  # Exclude the first tunnel state/ exclusion of the first tunnel state from the names() assignment is to ensure that the names assigned to the elements of the vector correspond correctly to the tunnel states.
    
    # Vector of DALY per cycle for CM under strategy No Intervention
    v_DALY_CM_LLIN <- rep(DALY_CM, n_tunnel_size_CM)
    names(v_DALY_CM_LLIN) <- v_CM_tunnel[-1]  # Exclude the first tunnel state
    
    # Vector of DALY under strategy No Intervention
    v_DALY_LLIN <- c(W = DALY_W, 
                     MM= DALY_MM,
                     SMA = v_DALY_CM_LLIN,
                     CM = v_DALY_SMA_LLIN,
                     D=0) 
    
    
    
    # Initialize vector to store wage loss adjustments for each cycle in SMA tunnel
    wage_loss_per_tunnel_SM <- rep((wage_loss * time_loss_work_SM) / 2, n_tunnel_size_SMA)
    
    # Adjust wage loss for the second tunnel
    wage_loss_per_tunnel_SM[2] <- (wage_loss * time_loss_work_SM) / 2
    
    # Calculate cost vectors for SMA without including wage loss
    v_c_SMA_LLIN_without_wageloss <- rep(c_SMA, n_tunnel_size_SMA)
    
    # Calculate total costs for SMA, including wage loss adjustments for each cycle
    v_c_SMA_LLIN <- ifelse(1:n_tunnel_size_SMA == 1, wage_loss_per_tunnel_SM[1] + c_SMA, wage_loss_per_tunnel_SM) 
    
    
    # Initialize vector to store wage loss adjustments for each cycle in SMA tunnel
    wage_loss_per_tunnel_SM <- rep((wage_loss * time_loss_work_SM) / 2, n_tunnel_size_CM)
    
    # Adjust wage loss for the second tunnel
    wage_loss_per_tunnel_SM[2] <- (wage_loss * time_loss_work_SM) / 2
    
    # Calculate cost vectors for SMA without including wage loss
    v_c_CM_LLIN_without_wageloss <- rep(c_CM, n_tunnel_size_CM)
    
    # Calculate total costs for SMA, including wage loss adjustments for each cycle
    v_c_CM_LLIN <- ifelse(1:n_tunnel_size_CM == 1, wage_loss_per_tunnel_SM[1] + c_SMA, wage_loss_per_tunnel_SM) 
    
    
    # Define cost vector for individuals who are well and using LLIN
    v_c_W <- c_W
    
    # Vector of costs under strategy No Intervention
    v_c_LLIN <- c(W = v_c_W, 
                  MM = c_MM,
                  SMA = v_c_SMA_LLIN,
                  CM = v_c_CM_LLIN, 
                  D = 0)
    
    
    # Now, let's create vectors for "R21" strategy
    # (You can adapt this part according to your specific DALY and cost values for "R21")
    
    # Vector of DALYS per cycle for SMA under strategy R21
    v_DALY_SMA_R21 <- rep(DALY_SMA_R21, n_tunnel_size_SMA)
    names(v_DALY_SMA_R21) <- v_SMA_tunnel[-1]  # Exclude the first tunnel state
    
    # Vector of DALYS per cycle for CM under strategy R21
    v_DALY_CM_R21 <- rep(DALY_CM_R21, n_tunnel_size_CM)
    names(v_DALY_CM_R21) <- v_CM_tunnel[-1]  # Exclude the first tunnel state
    
    # Vector of DALYS under strategy R21
    v_DALY_R21 <- c(W = DALY_W_R21, 
                    MM= DALY_MM_R21,
                    SMA = v_DALY_SMA_R21,
                    CM = v_DALY_CM_R21, 
                    D= 0 ) 
    
    
    
    # Initialize vector to store wage loss adjustments for each cycle in SMA tunnel
    wage_loss_per_tunnel_SM_R21 <- rep((wage_loss * time_loss_work_SM) / 2, n_tunnel_size_SMA)
    
    # Adjust wage loss for the second tunnel
    wage_loss_per_tunnel_SM_R21[2] <- (wage_loss * time_loss_work_SM) / 2
    
    # Calculate cost vectors for SMA without including wage loss
    v_c_SMA_without_wageloss_R21 <- rep(c_SMA_R21, n_tunnel_size_SMA)
    
    # Calculate total costs for SMA, including wage loss adjustments for each cycle
    v_c_SMA_R21 <- ifelse(1:n_tunnel_size_SMA == 1, wage_loss_per_tunnel_SM[1] + c_SMA_R21, wage_loss_per_tunnel_SM) 
    
    # Initialize vector to store wage loss adjustments for each cycle in SMA tunnel
    wage_loss_per_tunnel_SM <- rep((wage_loss * time_loss_work_SM) / 2, n_tunnel_size_CM)
    
    # Adjust wage loss for the second tunnel
    wage_loss_per_tunnel_SM[2] <- (wage_loss * time_loss_work_SM) / 2
    
    # Calculate cost vectors for SMA without including wage loss
    v_c_CM_without_wageloss_R21 <- rep(c_CM_R21, n_tunnel_size_CM)
    
    # Calculate total costs for SMA, including wage loss adjustments for each cycle
    v_c_CM_R21 <- ifelse(1:n_tunnel_size_CM == 1, wage_loss_per_tunnel_SM[1] + c_SMA_R21, wage_loss_per_tunnel_SM) 
    
    v_c_W_R21 <- c_W_R21
    
    # Vector of costs under strategy R21
    v_c_R21 <- c(W =v_c_W_R21, 
                 MM= c_MM_R21,
                 SMA = v_c_SMA_R21,
                 CM = v_c_CM_R21,
                 D= 0)  
    
    
    ## Store state rewards ----
    #* Store the vectors of state DALY for each strategy in a list 
    l_DALY   <- list(v_DALY_LLIN,
                     v_DALY_R21)
    #* Store the vectors of state cost for each strategy in a list 
    l_c   <- list(v_c_LLIN,
                  v_c_R21)
    #* assign strategy names to matching items in the lists
    names(l_DALY) <- names(l_c) <- v_names_str
    
    # Compute expected outcomes ----
    #* Create empty vectors to store total DALYs and costs 
    v_tot_DALY <- v_tot_cost <- vector(mode = "numeric", length = n_str)
    names(v_tot_DALY) <- names(v_tot_cost) <- v_names_str
    
    
    ## Loop through each strategy and calculate total DALY and costs ----
    for (i in 1:n_str) { # i <- 1
      v_DALY_str         <- l_DALY[[i]]   # select the vector of state utilities for the i-th strategy
      v_c_str         <- l_c[[i]]   # select the vector of state costs for the i-th strategy
      a_A_tunnels_str <- l_a_A[[i]] # select the transition array for the i-th strategy
      
      #* Create transition matrices of state DALY and state costs for the i-th strategy
      m_DALY_str  <- matrix(v_DALY_str, nrow = n_states_tunnels, ncol = n_states_tunnels, byrow = T)
      m_c_str  <- matrix(v_c_str, nrow = n_states_tunnels, ncol = n_states_tunnels, byrow = T)
      
      #* Expand the transition matrix of state utilities across cycles to form a transition array of state DALYs
      a_R_DALY_str <- array(m_DALY_str, 
                            dim      = c(n_states_tunnels, n_states_tunnels, n_cycles + 1),
                            dimnames = list(v_names_states_tunnels, v_names_states_tunnels, 0:n_cycles))
      # Expand the transition matrix of state costs across cycles to form a transition array of state costs
      a_R_c_str <- array(m_c_str, 
                         dim      = c(n_states_tunnels, n_states_tunnels, n_cycles + 1),
                         dimnames = list(v_names_states_tunnels, v_names_states_tunnels, 0:n_cycles))
      
      
      ###* Expected DALYs and costs per cycle
      ##* Vector of DALYs and costs
      v_DALY_str <- apply(a_A_tunnels_str* a_R_DALY_str , 3, sum) # sum the proportion of the cohort across transitions 
      v_c_str <- apply(a_A_tunnels_str*a_R_c_str, 3, sum) # sum the proportion of the cohort across transitions
      
      ##* Discounted total expected DALYs and Costs per strategy & half cycle correction 
      #* DALYs
      v_tot_DALY[i] <- t(v_DALY_str) %*% (v_dwe* gen_wcc(n_cycles, method = "half-cycle"))
      #* Costs
      v_tot_cost[i] <- t(v_c_str) %*% (v_dwc* gen_wcc(n_cycles, method = "half-cycle"))
      
    }
    
    ## Vector with discounted net monetary benefits (NMB)
    v_nmb <- v_tot_DALY * n_wtp - v_tot_cost
    
    ## data.frame with discounted costs, effectiveness and NMB
    #df_ce <- data.frame(Strategy = v_names_str,
    # Cost     = v_tot_cost,
    #Effect   = v_tot_DALY,
    #NMB      = v_nmb)
    
    #v_averted_DALY <- v_tot_DALY[1] - v_tot_DALY
    ## data.frame with discounted costs, effectiveness and NMB
    df_ce <- data.frame(Strategy = v_names_str,
                        Cost     = v_tot_cost,
                        Effect   = v_tot_DALY,
                        NMB      = v_nmb)
    
    return(df_ce)
  }
  )
}

calculate_ce_out(l_params_all)
#* Test function to generate PSA input dataset
generate_psa_params(1000) # Function included in "R/Functions_cSTM_time_dep_state_residence.R"
n_sim <- 1000


#* Generate PSA input dataset
df_psa_input <- generate_psa_params(n_sim = n_sim)
#* First six observations
head(df_psa_input)


### Histogram of parameters ----

ggplot(melt(df_psa_input, variable.name = "Parameter"), aes(x = value)) +
  facet_wrap(~Parameter, scales = "free") +
  geom_histogram(aes(y = ..density..)) +
  ylab("") +
  theme_bw(base_size = 16) + 
  theme(axis.text = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank())



##run PSA-----
#* data.frame of costs
df_c <- as.data.frame(matrix(0, 
                             nrow = n_sim,
                             ncol = n_str))
colnames(df_c) <- v_names_str
#* data.frame of effectiveness
df_e <- as.data.frame(matrix(0, 
                             nrow = n_sim,
                             ncol = n_str))
colnames(df_e) <- v_names_str


#* data.frame of prevalence of LLIN
m_prev <- matrix(NA, nrow = n_sim, ncol = (n_cycles + 1), 
                 dimnames = list(1:n_sim, 0:n_cycles))
df_prevLLIN <- data.frame(States = "LLIN",
                          m_prev, check.names = "FALSE")
#* data.frame of prevalence of R21
df_prevR21 <- data.frame(States = "R21",
                         m_prev, check.names = "FALSE")


#* Conduct probabilistic sensitivity analysis
#* Run Markov model on each parameter set of PSA input dataset
n_time_init_psa_series <- Sys.time()
for (i in 1:n_sim) {# i <- 1
  l_psa_input <- update_param_list(l_params_all, df_psa_input[i,])
  
  # Economics Measures
  l_out_ce_temp  <- calculate_ce_out(l_psa_input)
  df_c[i, ]  <- l_out_ce_temp$Cost  
  df_e[i, ]  <- l_out_ce_temp$Effect
  
  # Display simulation progress
  if (i/(n_sim/100) == round(i/(n_sim/100), 0)) { # display progress every 5%
    cat('\r', paste(i/n_sim * 100, "% done", sep = " "))
  }
}



############################################################################
##############################PSA PLOT Fucntions ###############
############################################################################

plot.psa <- function(x,
                     center = TRUE, ellipse = TRUE,
                     alpha = 0.2, txtsize = 12, col = c("full", "bw"),
                     n_x_ticks = 6, n_y_ticks = 6,
                     xbreaks = NULL,
                     ybreaks = NULL,
                     xlim = NULL,
                     ylim = NULL,
                     ...) {
  
  effectiveness <- x$effectiveness
  cost <- x$cost
  strategies <- x$strategies
  currency <- x$currency
  
  # expect that effectiveness and costs have strategy column names
  # removes confusing 'No id variables; using all as measure variables'
  df_cost <- suppressMessages(
    pivot_longer(cost,
                 everything(),
                 names_to = "Strategy",
                 values_to = "Cost")
  )
  df_effect <- suppressMessages(
    pivot_longer(effectiveness,
                 cols = everything(),
                 names_to = "Strategy",
                 values_to = "Effectiveness")
  )
  ce_df <- data.frame("Strategy" = df_cost$Strategy,
                      "Cost" = df_cost$Cost,
                      "Effectiveness" = df_effect$Effectiveness)
  
  # make strategies in psa object into ordered factors
  ce_df$Strategy <- factor(ce_df$Strategy, levels = strategies, ordered = TRUE)
  
  psa_plot <- ggplot(ce_df, aes_string(x = "Effectiveness", y = "Cost", color = "Strategy")) +
    geom_point(size = 0.7, alpha = alpha, shape = 21) +
    ylab(paste("Cost (", currency, ")", sep = ""))
  
  # define strategy-specific means for the center of the ellipse
  if (center) {
    strat_means <- ce_df %>%
      group_by(.data$Strategy) %>%
      summarize(Cost.mean = mean(.data$Cost),
                Eff.mean = mean(.data$Effectiveness))
    # make strategies in psa object into ordered factors
    strat_means$Strategy <- factor(strat_means$Strategy, levels = strategies, ordered = TRUE)
    psa_plot <- psa_plot +
      geom_point(data = strat_means,
                 aes_string(x = "Eff.mean", y = "Cost.mean", fill = "Strategy"),
                 size = 8, shape = 21, color = "black")
  }
  
  if (ellipse) {
    # make points for ellipse plotting
    df_list_ell <- lapply(strategies, function(s) {
      strat_specific_df <- ce_df[ce_df$Strategy == s, ]
      els <-  with(strat_specific_df,
                   ellipse::ellipse(cor(Effectiveness, Cost),
                                    scale = c(sd(Effectiveness), sd(Cost)),
                                    centre = c(mean(Effectiveness), mean(Cost))))
      data.frame(els, group = s, stringsAsFactors = FALSE)
    })
    df_ell <- bind_rows(df_list_ell)
    # draw ellipse lines
    psa_plot <- psa_plot + geom_path(data = df_ell,
                                     aes_string(x = "x", y = "y", colour = "group"),
                                     size = 1, linetype = 2, alpha = 1)
  }
  
  # add common theme
  col <- match.arg(col)
  add_common_aes(psa_plot, txtsize, col = col, col_aes = c("color", "fill"),
                 continuous = c("x", "y"),
                 n_x_ticks = n_x_ticks, n_y_ticks = n_y_ticks,
                 xbreaks = xbreaks, ybreaks = ybreaks,
                 xlim = xlim, ylim = ylim)
}





create_sa <- function(parameters, parnames, effectiveness, strategies,
                      cost, currency, other_outcome) {
  # checks that each is a dataframe
  if (!is.null(cost)) {
    cost <- check_df_and_coerce(cost)
  }
  
  if (!is.null(other_outcome)) {
    other_outcome <- check_df_and_coerce(other_outcome)
  }
  
  if (!is.null(effectiveness)) {
    effectiveness <- check_df_and_coerce(effectiveness)
  }
  
  if (!is.null(parameters)) {
    parameters <- check_df_and_coerce(parameters)
  }
  
  ### argument checks and definitions of other variables ###
  
  # costs, effectiveness, and parameters have same number of rows
  n_sim_ls <- list(effectiveness, cost, parameters, other_outcome)
  if (length(unique(unlist(lapply(n_sim_ls[!unlist(lapply(n_sim_ls, is.null))], nrow)))) != 1) {
    stop("Among those provided, the cost, effectiveness, parameter,
         and other_outcome dataframes must all have the same number of rows.")
  }
  
  # define n_sim
  n_sim <- unique(unlist(lapply(n_sim_ls[!unlist(lapply(n_sim_ls, is.null))], nrow)))
  
  # costs and effectiveness have same number of columns (strategies)
  n_strategies_ls <- list(effectiveness, cost, other_outcome)
  if (length(unique(unlist(lapply(n_strategies_ls[!unlist(lapply(n_strategies_ls, is.null))], ncol)))) != 1) {
    stop("Among those provided, the cost, effectiveness,
         and other_outcome dataframes must all have the same number of columns.")
  }
  
  # define n_strategies
  n_strategies <- unique(unlist(lapply(n_strategies_ls[!unlist(lapply(n_strategies_ls, is.null))], ncol)))
  
  # If the strategy names are not provided, generate a generic vector
  # with strategy names
  if (is.null(strategies)) {
    strategies <- paste(rep("Strategy_", n_strategies), seq(1, n_strategies), sep = "")
  } else {
    # correct strategy names. they are used as data.frame column names and in lm()
    # so they need to be syntactically valid
    new_strategies <- make.names(strategies, unique = TRUE)
    
    # write warning to console, so user knows that strategy name was changed
    for (i in 1:n_strategies) {
      old_strat <- strategies[i]
      new_strat <- new_strategies[i]
      if (new_strat != old_strat) {
        warning(paste0("strategy name '", old_strat, "' was converted to '", new_strat,
                       "' for compatibility. See ?make.names"), call. = FALSE)
      }
    }
    # update strategies
    strategies <- new_strategies
    
    # make sure strategies is the same length as the number of columns
    if (n_strategies != length(strategies)) {
      stop(
        paste0("The number of columns in the cost and effectiveness",
               "matrices is different from the number of strategies provided"))
    }
  }
  
  # define cost and effectiveness column names using strategies
  if (!is.null(cost)) {
    names(cost) <- strategies
  }
  if (!is.null(effectiveness)) {
    names(effectiveness) <- strategies
  }
  
  # define sa as a named list
  sa <- list("n_strategies" = n_strategies,
             "strategies" = strategies,
             "n_sim" = n_sim,
             "cost" = cost,
             "effectiveness" = effectiveness,
             "other_outcome" = other_outcome,
             "parameters" = parameters,
             "parnames" = parnames,
             "currency" = currency)
  class(sa) <- "sa"
  return(sa)
}

make_psa_obj <- function(cost, effectiveness, parameters = NULL,
                         strategies = NULL, currency = "$", other_outcome = NULL) {
  
  # parameter names
  parnames <- names(parameters)
  
  # define psa as a named list
  psa_obj <- create_sa(parameters, parnames, effectiveness, strategies,
                       cost, currency, other_outcome)
  
  # give classes "psa" and "sa"
  class(psa_obj) <- c("psa", class(psa_obj))
  return(psa_obj)
}


check_psa_object <- function(psa) {
  if (!inherits(psa, "psa")) {
    stop(paste0("The psa results parameter must be an object of class `psa`.\n",
                "Please run the make_psa() function to create this object."))
  }
}

check_df_and_coerce <- function(obj) {
  obj_name <- deparse(substitute(obj))
  if (!inherits(obj, "data.frame")) {
    warning(paste0("\'", obj_name, "\'", " is not a data frame. coercing to data frame"))
    df <- as.data.frame(obj)
  } else {
    df <- as.data.frame(obj)
  }
  return(df)
}



check_psa_object <- function(psa) {
  if (!inherits(psa, "psa")) {
    stop(paste0("The psa results parameter must be an object of class `psa`.\n",
                "Please run the make_psa() function to create this object."))
  }
}
check_df_and_coerce <- function(obj) {
  obj_name <- deparse(substitute(obj))
  if (!inherits(obj, "data.frame")) {
    warning(paste0("\'", obj_name, "\'", " is not a data frame. coercing to data frame"))
    df <- as.data.frame(obj)
  } else {
    df <- as.data.frame(obj)
  }
  return(df)
}


calculate_icers <- function(cost, effect, strategy) {
  # Check input lengths
  if (length(cost) != length(effect) || length(cost) != length(strategy)) {
    stop("cost, effect, and strategy must have the same length", call. = FALSE)
  }
  
  # Calculate ICER
  inc_cost <- c(NA, diff(cost))
  inc_effect <- c(NA, abs(diff(effect)))
  icers <- inc_cost / inc_effect
  
  # Create data frame
  results <- data.frame("Strategy" = strategy,
                        "Cost" = cost,
                        "Effect" = effect,
                        "Inc_Cost" = inc_cost,
                        "Inc_Effect" = inc_effect,
                        "ICER" = icers)
  
  return(results)
}




plot.icers <- function(x,
                       txtsize = 12,
                       currency = "$",
                       effect_units = "DALYs",
                       label = c("frontier", "all", "none"),
                       label_max_char = NULL,
                       plot_frontier_only = FALSE,
                       alpha = 1,
                       n_x_ticks = 6,
                       n_y_ticks = 6,
                       xbreaks = NULL,
                       ybreaks = NULL,
                       xlim = NULL,
                       ylim = NULL,
                       xexpand = expansion(0.1),
                       yexpand = expansion(0.1),
                       max.iter = 20000,
                       ...) {
  if (ncol(x) > 7) {
    # reformat icers class object if uncertainty bounds are present
    x <- x %>%
      select(.data$Strategy, .data$Cost, .data$Effect,
             .data$Inc_Cost, .data$Inc_Effect,
             .data$ICER, .data$Status)
  }
  
  # type checking
  label <- match.arg(label)
  
  # this is so non-dominated strategies are plotted last (on top)
  x <- arrange(x, .data$Status)
  
  # change status text in data frame for plotting
  d_name <- "Dominated"
  ed_name <- "Weakly Dominated"
  nd_name <- "Efficient Frontier"
  
  status_expand <- c("D" = d_name, "ED" = ed_name,
                     "ND" = nd_name, "ref" = nd_name)
  x$Status <- factor(status_expand[x$Status], ordered = FALSE,
                     levels = c(d_name, ed_name, nd_name))
  
  # linetype
  plot_lines <- c("Dominated" = "blank",
                  "Weakly Dominated" = "blank",
                  "Efficient Frontier" = "solid")
  
  # names to refer to in aes_
  stat_name <- "Status"
  strat_name <- "Strategy"
  eff_name <- "Effect"
  cost_name <- "Cost"
  
  # frontier only
  if (plot_frontier_only) {
    plt_data <- x[x$Status == nd_name, ]
  } else {
    plt_data <- x
  }
  
  # make plot
  icer_plot <- ggplot(plt_data, aes_(x = as.name(eff_name), y = as.name(cost_name),
                                     shape = as.name(stat_name))) +
    geom_point(alpha = alpha, size = 2) +
    geom_line(aes_(linetype = as.name(stat_name), group = as.name(stat_name))) +
    scale_linetype_manual(name = NULL, values = plot_lines) +
    scale_shape_discrete(name = NULL) +
    labs(x = paste0("Effect (", effect_units, ")"),
         y = paste0("Cost (", currency, ")"))
  
  icer_plot <- add_common_aes(icer_plot, txtsize, col = "none",
                              continuous = c("x", "y"),
                              n_x_ticks = n_x_ticks, n_y_ticks = n_y_ticks,
                              xbreaks = xbreaks, ybreaks = ybreaks,
                              xlim = xlim, ylim = ylim,
                              xexpand = xexpand, yexpand = yexpand)
  
  # labeling
  if (label != "none") {
    if (!is.null(label_max_char)) {
      plt_data[, strat_name] <- str_sub(plt_data[, strat_name],
                                        start = 1L, end = label_max_char)
    }
    if (label == "all") {
      lab_data <- plt_data
    }
    if (label == "frontier") {
      lab_data <- plt_data[plt_data$Status == nd_name, ]
    }
    
    icer_plot <- icer_plot +
      ggrepel::geom_label_repel(data = lab_data,
                                aes_(label = as.name(strat_name)),
                                size = 3,
                                show.legend = FALSE,
                                max.iter = max.iter,
                                direction = "both")
  }
  return(icer_plot)
}
plot.ceac <- function(x,
                      frontier = TRUE,
                      points = TRUE,
                      currency = "$",
                      min_prob = 0,
                      txtsize = 12,
                      n_x_ticks = 10,
                      n_y_ticks = 8,
                      xbreaks = NULL,
                      ybreaks = NULL,
                      ylim = NULL,
                      xlim = c(0, NA),
                      col = c("full", "bw"),
                      ...) {
  wtp_name <- "WTP"
  prop_name <- "Proportion"
  strat_name <- "Strategy"
  x$WTP_thou <- x[, wtp_name] / 1000
  
  # removing strategies with probabilities always below `min_prob`
  # get group-wise max probability
  if (min_prob > 0) {
    max_prob <- x %>%
      group_by(.data$Strategy) %>%
      summarize(maxpr = max(.data$Proportion)) %>%
      filter(.data$maxpr >= min_prob)
    strat_to_keep <- max_prob$Strategy
    if (length(strat_to_keep) == 0) {
      stop(
        paste("no strategies remaining. you may want to lower your min_prob value (currently ",
              min_prob, ")", sep = "")
      )
    }
    # report filtered out strategies
    old_strat <- unique(x$Strategy)
    diff_strat <- setdiff(old_strat, strat_to_keep)
    n_diff_strat <- length(diff_strat)
    if (n_diff_strat > 0) {
      # report strategies filtered out
      cat("filtered out ", n_diff_strat, " strategies with max prob below ", min_prob, ":\n",
          paste(diff_strat, collapse = ","), "\n", sep = "")
      
      # report if any filtered strategies are on the frontier
      df_filt <- filter(x, .data$Strategy %in% diff_strat & .data$On_Frontier)
      if (nrow(df_filt) > 0) {
        cat(paste0("WARNING - some strategies that were filtered out are on the frontier:\n",
                   paste(unique(df_filt$Strategy), collapse = ","), "\n"))
      }
    }
    
    # filter dataframe
    x <- filter(x, .data$Strategy %in% strat_to_keep)
  }
  
  # Drop unused strategy names
  x$Strategy <- droplevels(x$Strategy)
  
  p <- ggplot(data = x, aes_(x = as.name("WTP_thou"),
                             y = as.name(prop_name),
                             color = as.name(strat_name))) +
    geom_line() +
    xlab(paste("Willingness to Pay (Thousand ", currency, " / DALY)", sep = "")) +
    ylab("Pr Cost-Effective")
  
  if (points) {
    p <- p + geom_point(aes_(color = as.name(strat_name)))
  }
  
  if (frontier) {
    front <- x[x$On_Frontier, ]
    p <- p + geom_point(data = front, aes_(x = as.name("WTP_thou"),
                                           y = as.name(prop_name),
                                           shape = as.name("On_Frontier")),
                        size = 3, stroke = 1, color = "black") +
      scale_shape_manual(name = NULL, values = 0, labels = "Frontier") +
      guides(color = guide_legend(order = 1),
             shape = guide_legend(order = 2))
  }
  col <- match.arg(col)
  add_common_aes(p, txtsize, col = col, col_aes = "color",
                 continuous = c("x", "y"), n_x_ticks = n_x_ticks, n_y_ticks = n_y_ticks,
                 xbreaks = xbreaks, ybreaks = ybreaks,
                 ylim = ylim, xlim = xlim)
}
calculate_outcome <- function(outcome = c("nhb", "nmb", "eff", "cost", "nhb_loss",
                                          "nmb_loss", "nhb_loss_voi", "nmb_loss_voi"),
                              cost, effect, wtp) {
  outcome <- match.arg(outcome)
  n_sim <- nrow(cost)
  if (outcome == "eff") {
    y <- effect
  } else if (outcome == "cost") {
    y <- cost
  } else {
    if (is.null(wtp)) {
      # the call. = FALSE makes the error message more clear
      stop("wtp must be provided for NHB and NMB",  call. = FALSE)
    }
    if (is.null(cost)) {
      stop("must provide cost for NHB and NMB.",  call. = FALSE)
    }
    if (outcome == "nhb") {
      y <- effect - cost / wtp
    }
    if (outcome == "nmb") {
      y <- effect * wtp - cost
    }
    if (outcome == "nhb_loss" | outcome == "nmb_loss") {
      if (outcome == "nhb_loss") {
        net_outcome <- "nhb"
      }
      if (outcome == "nmb_loss") {
        net_outcome <- "nmb"
      }
      netben <- calculate_outcome(net_outcome, cost, effect, wtp)
      max_str_rowwise <- max.col(netben)
      y <-  netben[cbind(1:n_sim, max_str_rowwise)] - netben
    }
    if (outcome == "nhb_loss_voi" | outcome == "nmb_loss_voi") {
      if (outcome == "nhb_loss_voi") {
        net_outcome <- "nhb"
      }
      if (outcome == "nmb_loss_voi") {
        net_outcome <- "nmb"
      }
      netben <- calculate_outcome(net_outcome, cost, effect, wtp)
      max_str <- which.max(colMeans(netben))
      y <- netben - netben[cbind(1:n_sim), max_str]
    }
  }
  return(y)
}

ceac <- function(wtp, psa) {
  # check that psa has class 'psa'
  check_psa_object(psa)
  
  # define needed variables
  strategies <- psa$strategies
  n_strategies <- psa$n_strategies
  effectiveness <- psa$effectiveness
  cost <- psa$cost
  n_sim <- psa$n_sim
  
  # number of willingness to pay thresholds
  n_wtps <- length(wtp)
  
  # matrix to store probability optimal for each strategy
  cea <- matrix(0, nrow = n_wtps, ncol = n_strategies)
  colnames(cea) <- strategies
  
  # vector to store strategy at the cost-effectiveness acceptability frontier
  frontv <- rep(0, n_wtps)
  
  for (l in 1:n_wtps) {
    # calculate net monetary benefit at wtp[l]
    lth_wtp <- wtp[l]
    nmb <-  calculate_outcome("nmb", cost, effectiveness, lth_wtp)
    
    # find the distribution of optimal strategies
    max.nmb <- max.col(nmb)
    opt <- table(max.nmb)
    cea[l, as.numeric(names(opt))] <- opt / n_sim
    
    # calculate point on CEAF
    # the strategy with the highest expected nmb
    frontv[l] <- which.max(colMeans(nmb))
  }
  
  # make cea df
  cea_df <- data.frame(wtp, cea, strategies[frontv],
                       stringsAsFactors = FALSE)
  colnames(cea_df) <- c("WTP", strategies, "fstrat")
  
  # Reformat df to long format
  ceac <- tidyr::pivot_longer(
    data = cea_df,
    cols = !c("WTP", "fstrat"),
    names_to = "Strategy",
    values_to = "Proportion"
  )
  
  # boolean for on frontier or not
  ceac$On_Frontier <- (ceac$fstrat == ceac$Strategy)
  
  # drop fstrat column
  ceac$fstrat <- NULL
  
  # order by WTP
  ceac <- ceac[order(ceac$WTP), ]
  
  # remove rownames
  rownames(ceac) <- NULL
  
  # make strategies in ceac object into ordered factors
  ceac$Strategy <- factor(ceac$Strategy, levels = strategies, ordered = TRUE)
  
  # define classes
  # defining data.frame as well allows the object to use print.data.frame, for example
  class(ceac) <- c("ceac", "data.frame")
  
  return(ceac)
}


############################################################################
################# PSA list scatterplot and CEA outcomes #################################
############################################################################

l_psa <- make_psa_obj(cost          = df_c, 
                      effectiveness = df_e, 
                      parameters    = df_psa_input, 
                      strategies    = v_names_str)
l_psa$strategies <- v_names_str
colnames(l_psa$effectiveness) <- v_names_str
colnames(l_psa$cost) <- v_names_str

#* Vector with willingness-to-pay (WTP) thresholds.
v_wtp <- seq(0, 200000, by = 5000)



### Cost-Effectiveness Scatter plot ----
txtsize <- 15

#* Function included in "R/Functions.R"; depends on `tidyr` and `ellipse` packages.
#* The latest version can be found in `dampack` package
gg_scattter <- plot.psa(l_psa, txtsize = txtsize) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  scale_y_continuous("Cost (Thousand $)", 
                     breaks = number_ticks(10),
                     labels = function(x) x/1000) +
  xlab("Effectiveness (QALYs)") +
  guides(col = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom")
gg_scattter


### Incremental cost-effectiveness ratios (ICERs) with probabilistic output ----
#* Compute expected costs and effects for each strategy from the PSA
#* Function included in "R/Functions.R". The latest version can be found in `dampack` package
df_out_ce_psa <- summary(l_psa)


#* Function included in "R/Functions.R"; depends on the `dplyr` package
#* The latest version can be found in `dampack` package
#* 

df_cea_psa <- calculate_icers(cost = df_out_ce_psa$meanCost, 
                              effect = df_out_ce_psa$meanEffect,
                              strategy = df_out_ce_psa$Strategy)
df_cea_psa



ggplot(df_cea_psa, aes(x = Inc_Effect, y = ICER)) +
  geom_point() +
  labs(x = "Incremental Effect (DALYs Averted)",
       y = "ICER",
       title = "Incremental Cost-Effectiveness Ratio (ICER) Scatter Plot")

### Plot cost-effectiveness frontier with probabilistic output ----
#* Function included in "R/Functions.R"; depends on the `ggplot2`  and `ggrepel` packages.
#* The latest version can be found in `dampack` package
plot.icers(df_cea_psa, label = "all", txtsize = txtsize) +
  expand_limits(x = max(table_cea$QALYs) + 0.1) +
  theme(legend.position = c(0.8, 0.2))


############################################################################
###############CEAC and CEAF calculations ###################################
############################################################################

## Cost-effectiveness acceptability curves (CEACs) and frontier (CEAF) ---
#* Functions included in "R/Functions.R". The latest versions can be found in `dampack` package
ceac_obj <- ceac(wtp = v_wtp, psa = l_psa)
#* Regions of highest probability of cost-effectiveness for each strategy
summary(ceac_obj)



#* CEAC & CEAF plot
gg_ceac <- plot.ceac(ceac_obj, txtsize = txtsize, xlim = c(0, 200), n_x_ticks = 14) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme(legend.position = c("bottom"))
gg_ceac


calculate_icers_psa <- function(psa, uncertainty = FALSE) {
  
  # check that psa has class 'psa'
  check_psa_object(psa)
  
  # Calculate mean outcome values
  psa_sum <- summary(psa)
  
  # Supply mean outcome values to calculate_icers
  icers <- calculate_icers(cost = psa_sum$meanCost,
                           effect = psa_sum$meanEffect,
                           strategies = psa_sum$Strategy)
  
  if (uncertainty == TRUE) {
    
    # extract cost and effect data.frames from psa object
    cost <- psa$cost
    effect <- psa$effectiveness
    
    # Calculate quantiles across costs and effects
    cost_bounds <- cost %>%
      pivot_longer(cols = everything(), names_to = "Strategy") %>%
      group_by(.data$Strategy) %>%
      summarize(Lower_95_Cost = quantile(.data$value, probs = 0.025, names = FALSE),
                Upper_95_Cost = quantile(.data$value, probs = 0.975, names = FALSE))
    
    effect_bounds <- effect %>%
      pivot_longer(cols = everything(), names_to = "Strategy") %>%
      group_by(.data$Strategy) %>%
      summarize(Lower_95_Effect = quantile(.data$value, probs = 0.025, names = FALSE),
                Upper_95_Effect = quantile(.data$value, probs = 0.975, names = FALSE))
    
    # merge bound data.frames into icers data.frame
    icers <- icers %>%
      left_join(cost_bounds, by = "Strategy") %>%
      left_join(effect_bounds, by = "Strategy") %>%
      select(.data$Strategy, .data$Cost, .data$Lower_95_Cost, .data$Upper_95_Cost,
             .data$Effect, .data$Lower_95_Effect, .data$Upper_95_Effect,
             .data$Inc_Cost, .data$Inc_Effect, .data$ICER, .data$Status)
  }
  
  return(icers)
}
plot.exp_loss <- function(x,
                          log_y = TRUE,
                          frontier = TRUE,
                          points = TRUE,
                          lsize = 1,
                          txtsize = 12,
                          currency = "$",
                          effect_units = "DALY",
                          n_y_ticks = 8,
                          n_x_ticks = 20,
                          xbreaks = NULL,
                          ybreaks = NULL,
                          xlim = c(0, NA),
                          ylim = NULL,
                          col = c("full", "bw"),
                          ...) {
  wtp_name <- "WTP_thou"
  loss_name <- "Expected_Loss"
  strat_name <- "Strategy"
  x[, wtp_name] <- x$WTP / 1000
  
  # split into on frontier and not on frontier
  nofront <- x
  front <- x[x$On_Frontier, ]
  
  # Drop unused levels from strategy names
  nofront$Strategy <- droplevels(nofront$Strategy)
  front$Strategy <- droplevels(front$Strategy)
  # formatting if logging the y axis
  if (log_y) {
    tr <- "log10"
  } else {
    tr <- "identity"
  }
  
  p <- ggplot(data = nofront, aes_(x = as.name(wtp_name),
                                   y = as.name(loss_name))) +
    xlab(paste0("Willingness to Pay (Thousand ", currency, "/", effect_units, ")")) +
    ylab(paste0("Expected Loss (", currency, ")"))
  
  # color
  col <- match.arg(col)
  ## change linetype too if color is black and white
  if (col == "full") {
    if (points) {
      p <- p + geom_point(aes_(color = as.name(strat_name)))
    }
    p <- p +
      geom_line(size = lsize, aes_(color = as.name(strat_name)))
    
  }
  if (col == "bw") {
    if (points) {
      p <- p + geom_point()
    }
    p <- p +
      geom_line(aes_(linetype = as.name(strat_name)))
  }
  
  p <- add_common_aes(p, txtsize, col = col, col_aes = c("color", "line"),
                      continuous = c("x", "y"),
                      n_x_ticks = n_x_ticks, n_y_ticks = n_y_ticks,
                      xbreaks = xbreaks, ybreaks = ybreaks,
                      xlim = xlim, ylim = ylim,
                      ytrans = tr)
  if (frontier) {
    p <- p + geom_point(data = front, aes_(x = as.name(wtp_name),
                                           y = as.name(loss_name),
                                           shape = as.name("On_Frontier")),
                        size = 3, stroke = 1, color = "black") +
      scale_shape_manual(name = NULL, values = 0, labels = "Frontier & EVPI") +
      guides(color = guide_legend(order = 1),
             linetype = guide_legend(order = 1),
             shape = guide_legend(order = 2))
  }
  return(p)
}
calc_evpi <- function(psa, wtp, pop = 1) {
  check_psa_object(psa)
  cost <- psa$cost
  effectiveness <- psa$effectiveness
  if (ncol(effectiveness) < 2) {
    stop("You need at least two different strategies to compute EVPI.")
  }
  # number of wtp thresholds
  n_wtps <- length(wtp)
  # vector to store evpi
  evpi <- rep(0, n_wtps)
  # Estimate the Loss matrix and EVPI at each WTP threshold
  for (l in 1:n_wtps) {
    ## Calculate the opportunity loss from choosing d.star for each strategy
    loss <- calculate_outcome("nmb_loss", cost, effectiveness, wtp[l])
    
    ## Compute EVPI
    evpi[l] <- min(apply(loss, 2, mean)) * pop
  }
  
  # Data frame to store EVPI for each WTP threshold
  df_evpi <- data.frame("WTP" = wtp, "EVPI" = evpi)
  
  # declare class as both evpi (plotting) and data.frame (printing)
  class(df_evpi) <- c("evpi", "data.frame")
  return(df_evpi)
}

plot.evpi <- function(x,
                      txtsize = 12,
                      currency = "$",
                      effect_units = "DALY",
                      n_y_ticks = 8,
                      n_x_ticks = 20,
                      xbreaks = NULL,
                      ybreaks = NULL,
                      xlim = c(0, NA),
                      ylim = NULL,
                      ...) {
  x$WTP_thou <- x$WTP / 1000
  g <- ggplot(data = x,
              aes_(x = as.name("WTP_thou"), y = as.name("EVPI"))) +
    geom_line() +
    xlab(paste("Willingness to Pay (Thousand ", currency, "/", effect_units, ")", sep = "")) +
    ylab(paste("EVPI (", currency, ")", sep = ""))
  add_common_aes(g, txtsize, continuous = c("x", "y"),
                 n_x_ticks = n_x_ticks, n_y_ticks = n_y_ticks,
                 xbreaks = xbreaks, ybreaks = ybreaks,
                 xlim = xlim, ylim = ylim)
}




############################################################################
##################################* CEAC & CEAF plot
#* ############################################################################
gg_ceac <- plot.ceac (ceac_obj, txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme(legend.position = c(0.8, 0.48))
gg_ceac



## Expected Loss Curves (ELCs) ----
#* Function included in "R/Functions.R".The latest version can be found in `dampack` package
elc_obj <- calc_exp_loss(wtp = v_wtp, psa = l_psa)
elc_obj

#* ELC plot
gg_elc <- plot.exp_loss(elc_obj, log_y = FALSE, 
                        txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14,
                        col = "full") +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  # geom_point(aes(shape = as.name("Strategy"))) +
  scale_y_continuous("Expected Loss (Thousand $)", 
                     breaks = number_ticks(10),
                     labels = function(x) x/1000) +
  theme(legend.position = c(0.4, 0.7),)
gg_elc

## Expected value of perfect information (EVPI) ----
#* Function included in "R/Functions.R". The latest version can be found in `dampack` package
evpi <- calc_evpi(wtp = v_wtp, psa = l_psa)
#* EVPI plot
gg_evpi <- plot.evpi(evpi, effect_units = "DALY", 
                     txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14) +
  scale_y_continuous("EVPI (Thousand $)", 
                     breaks = number_ticks(10),
                     labels = function(x) x/1000)
gg_evpi


### Combine all figures into one ----
patched_cea <- (gg_scattter +  gg_ceac + plot_layout(guides = "keep"))/(gg_elc + gg_evpi)
gg_psa_plots <- patched_cea + 
  plot_annotation(tag_levels = 'A')
gg_psa_plots


##############################################################################
###########one way sensitivity analysis ######################################
##############################################################################


##NORMAL PSA
sensitivity_results <- owsa(
  l_psa,
  params = NULL,
  ranges = NULL,
  nsamp = 100,
  outcome = c( "nmb"),
  wtp = 21377.49,
  strategies = NULL,
  poly.order = 2
)

tornado_plot <- owsa_tornado(
  owsa = sensitivity_results,
  return = "plot",  # Specify to return a plot
  txtsize = 8,
  min_rel_diff = 0,
  col = "full",
  n_y_ticks = 12,
  ylim = NULL,
  ybreaks = NULL
)
print (tornado_plot)
















