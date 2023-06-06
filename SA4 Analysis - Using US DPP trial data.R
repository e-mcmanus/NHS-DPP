########################################
# Markov model: Evaluating the NHS DPP #
# 4 modelling states:
# Normal glycaemic control
# NDH
# Type 2 Diabetes
# Death
#SA4: US DPP trial data used as parameter source in Herman (42)
######################################## 

getwd()
## SET WORKING DIRECTORY ##
setwd("")

library('ggplot2')
library('ggthemes')


#######
# AGE CAT 1
#######

set.seed(9009)
## model set-up ----

t_names <- c("without DPP", "with DPP")
n_treatments <- length(t_names)

s_names  <- c("Normal", "NDH", "T2D", "Dead")
n_states <- length(s_names)

n_cohort <- 58
n_cycles <- 36
Initial_age <- 35


#NUMBER OF PSA ITERATIONS
n_trials<-10000
costs <- matrix(NA,nrow=n_trials,ncol=n_treatments, dimnames=list(NULL,t_names))
qalys <- matrix(NA,nrow=n_trials,ncol=n_treatments, dimnames=list(NULL,t_names))


#Utility of states
# uNorm <- 0.8395923 / uNDH <- 0.8085401 / uT2D <- 0.708957 / uDeath <-0 / uDPP<-0 


#Cost of states
cNorm <- rgamma(1, 44.44444, scale=45.119475)
cNDH <- rgamma(1,  44.44444, scale=50.0571)
cT2D <- rgamma(1,  44.44444, scale=99.4628)
cDeath <- 0
cDPP <- 0

#Utility of states
uNorm <- rbeta(1,27.53076303,4.883160879)
uNDH <- rbeta(1,26.33007731, 6.334144666)
uT2D <- rbeta(1, 26.52841, 10.65937)
uDeath <-0
uDPP<-0 
#Discount rates
oDr <- 0.035
cDr <- 0.035

#utility gains from DPP
utility <- c(1-rbeta(n_cohort,801123.6265,3464.71833))
utility_sum <-sum(utility)

#costs of the DPP
dppcost <- c(rgamma(n_cohort,16.73964,scale=8.469118))
dppcost_sum <- sum(dppcost)

#Transition probabilities
#normal to NDH
tpProgNDH <- 0.074202731
#dying from NGT or NDH
tpDeath <- 0.05
#Progression to T2D from NDH
tpProgT2D <- 0.0249
#From NDH to NGT
tpNorm <- 0.0795
#From T2D to NDH
tpRegress <- 0.002796084
#Additional risk of death for T2D
tpExcessDeath <- 1.6

#Effectiveness of the programme
effect <- 0.8

# cost of staying in state

state_c_matrix <-
  matrix(c(cNorm, cNDH, cT2D, 0,
           cNorm, cNDH + cDPP, cT2D, 0),
         byrow = TRUE,
         nrow = n_treatments,
         dimnames = list(t_names,
                         s_names))

# qaly when staying in state
state_q_matrix <-
  matrix(c(uNorm, uNDH, uT2D, 0,
           uNorm, uNDH + uDPP, uT2D, 0),
         byrow = TRUE,
         nrow = n_treatments,
         dimnames = list(t_names,
                         s_names))

# cost of moving to a state
# same for both treatments
#only cost associated with death - no cost in this model
trans_c_matrix <-
  matrix(c(0, 0, 0, cDeath,
           0, 0, 0, cDeath,
           0, 0, 0, cDeath,
           0, 0, 0, 0),
         byrow = TRUE,
         nrow = n_states,
         dimnames = list(from = s_names,
                         to = s_names))

# Transition probabilities ---- 

# time-homogeneous
#not by row - by column
p_matrix <- array(data = c(1-tpProgNDH - tpDeath, tpNorm, 0, 0,
                           tpProgNDH, 1-tpProgT2D-tpDeath-tpNorm, tpRegress, 0,
                           0, tpProgT2D,  1-tpRegress - (tpDeath*tpExcessDeath), 0,
                           tpDeath, tpDeath, (tpDeath*tpExcessDeath), 1,
                           #below would be adjusted for effect of DPP
                           1-tpProgNDH - tpDeath, tpNorm, 0, 0,
                           tpProgNDH, 1-(tpProgT2D*effect)-tpDeath-tpNorm, tpRegress, 0,
                           0, (tpProgT2D*effect),  1-tpRegress - (tpDeath*tpExcessDeath), 0,
                           tpDeath, tpDeath, (tpDeath*tpExcessDeath), 1),
                  dim = c(n_states, n_states, n_treatments),
                  dimnames = list(from = s_names,
                                  to = s_names,
                                  t_names))

# Store population output for each cycle 

# state populations
pop <- array(data = NA,
             dim = c(n_states, n_cycles, n_treatments),
             dimnames = list(state = s_names,
                             cycle = NULL,
                             treatment = t_names))

pop["Normal", cycle = 1, ] <-0
pop["NDH", cycle = 1, ] <- n_cohort
pop["T2D", cycle = 1, ] <- 0
pop["Dead", cycle = 1, ] <- 0

halfpop<-array(data=NA, dim=c(n_states,n_cycles,n_treatments), 
               dimnames = list(state= s_names, cycle=NULL, treatment= t_names))

halfpop["Normal", cycle=1, treatment="without DPP"]<- 0
halfpop["NDH", cycle=1, treatment="without DPP"]<- 0
halfpop["T2D", cycle=1, treatment="without DPP"]<- 0
halfpop["Dead", cycle=1, treatment="without DPP"]<- 0
halfpop["Normal", cycle=1, treatment="with DPP"]<- 0
halfpop["NDH", cycle=1, treatment="with DPP"]<- 0
halfpop["T2D", cycle=1, treatment="with DPP"]<- 0
halfpop["Dead", cycle=1, treatment="with DPP"]<- 0

# _arrived_ state populations
trans <- array(data = NA,
               dim = c(n_states, n_cycles, n_treatments),
               dimnames = list(state = s_names,
                               cycle = NULL,
                               treatment = t_names))

trans[, cycle = 1, ] <- 0


# Sum costs and QALYs for each cycle at a time for each treatment arm

cycle_empty_array <-
  array(NA,
        dim = c(n_treatments, n_cycles),
        dimnames = list(treatment = t_names,
                        cycle = NULL))

cycle_state_costs <- cycle_trans_costs <- cycle_qalys <- cycle_empty_array 
cycle_costs <- cycle_QALYs <- cycle_empty_array
LE <- LYs <- cycle_empty_array    # life expectancy; life-years
cycle_QALE <- cycle_empty_array   # quality-adjusted life expectancy

total_costs <- setNames(c(NA, NA), t_names)
total_QALYs <- setNames(c(NA, NA), t_names)

# Transition probabilities
#normal to NDH
tpProgNDH <- 0.074202731
# dying from NGT or NDH
tpDeath <- 0.05
# Progression from NDH to T2D
tpProgT2D <- 0.0249
# From NDH to NGT
tpNorm <- 0.0795
# From T2D to NDH
tpRegress <- 0.002796084
# Additional risk of death for T2D
tpExcessDeath <- 1.6
#
#Effectiveness of the DPP
effect <- 0.8

tpNorm <- rbeta(1,335,3709)
tpProgT2D <- rbeta(1,25.27,974.73)
tpRegress <- rbeta(1,2.8,997.2)
tpProgNDH <- rbeta(1,77.1,922.9)
tpExcessDeath <- rlnorm(1, meanlog = 0.470004, sdlog = 0.064111)

# # Time-dependent probability matrix ----
# Time dependent transition - risk of death
p_matrix_cycle <- function(p_matrix, age, cycle,
                           tpProgNDH = rbeta(1,77.1,922.9),
                           #tpProgT2D = rbeta(1,25.27,974.73),
                           tpNorm = rbeta(1,335,3709),
                           tpRegress = rbeta(1,2.8,997.2),
                           tpExcessDeath = rlnorm(1, meanlog = 0.470004, sdlog = 0.064111)
) {
  tpDeath_lookup <-
    c("[30,34]" = 0.000674,
      "[35,39]" = 0.001007,
      "[40,44]" = 0.001499,
      "[45,49]" = 0.002298,
      "[50,54]" = 0.003359,
      "[55,59]" = 0.005056,
      "[60,64]" = 0.007940,
      "[65,69]" = 0.012389,
      "[70,74]" = 0.019640,
      "[75,79]" = 0.034483,
      "[80,84]" = 0.060938,
      "[85,89]" = 0.110627,
      "[90,94]" = 0.187364,
      "[95,120]" = 0.295984)
  
  age_grp <- cut(age, breaks = c(30,34,39,44,49,54,59,64,69,74,79, 84, 89, 94, 120))
  tpDeath <- tpDeath_lookup[age_grp]
  
  #Normal distribution rnorm(1,mean=0.108, sd=0.0095)
  tpProgT2D_lookup <-
    c("[30,39]" = rnorm(1,mean=0.108, sd=0.0095),
      "[40,49]" = rnorm(1,mean=0.108, sd=0.0095),
      "[50,59]" = rnorm(1,mean=0.108, sd=0.0095),
      "[60,69]" = rnorm(1,mean=0.108, sd=0.0095),
      "[70,79]" = rnorm(1,mean=0.108, sd=0.0095),
      "[80,120]" = rnorm(1,mean=0.108, sd=0.0095))
  
  age_grpT2D <- cut(age, breaks= c(30,39,49,59,69,79,120))
  tpProgT2D <-tpProgT2D_lookup[age_grpT2D]
  
  # [ includes, ( doesn't include
  effect_lookup <-
    c("(0,3]" = rlnorm(1, meanlog = -0.22314, sdlog = 0.044757),
      "(3,10]" = rnorm(1, mean = 1.000, sd = 0.00000001),
      "(10,15]"= rnorm(1, mean = 1.000, sd = 0.00000001),
      "(15,50]" = rnorm(1, mean = 1.000, sd = 0.00000001))
  
  cycle2 <- cut(cycle, breaks = c(0,3,10,15,50))
  
  effect <- effect_lookup[cycle2]
  
  #########
  
  p_matrix["Normal", "Normal", "without DPP"] <- 1 - tpProgNDH - tpDeath
  p_matrix["Normal", "NDH", "without DPP"] <- tpProgNDH 
  p_matrix["Normal", "Dead", "without DPP"] <- tpDeath
  p_matrix["NDH", "Normal", "without DPP"] <- tpNorm
  p_matrix["NDH", "NDH", "without DPP"] <-1-tpProgT2D-tpDeath-tpNorm
  p_matrix["NDH", "T2D", "without DPP"] <- tpProgT2D
  p_matrix["NDH", "Dead", "without DPP"] <- tpDeath
  p_matrix["T2D", "NDH", "without DPP"] <- tpRegress
  p_matrix["T2D", "T2D", "without DPP"] <- 1-tpRegress - (tpDeath*tpExcessDeath)
  p_matrix["T2D", "Dead", "without DPP"] <- (tpDeath)*tpExcessDeath
  p_matrix["Dead", "Dead", "without DPP"] <- 1
  
  p_matrix["Normal", "Normal", "with DPP"] <- 1 - tpProgNDH - tpDeath
  p_matrix["Normal", "NDH", "with DPP"] <- tpProgNDH
  p_matrix["Normal", "Dead", "with DPP"] <- tpDeath
  p_matrix["NDH", "Normal", "with DPP"] <- tpNorm
  p_matrix["NDH", "NDH", "with DPP"] <- 1 - (tpProgT2D*effect) - tpDeath - tpNorm
  p_matrix["NDH", "T2D", "with DPP"] <- (tpProgT2D*effect)
  p_matrix["NDH", "Dead", "with DPP"] <- tpDeath
  p_matrix["T2D", "NDH", "with DPP"] <- tpRegress
  p_matrix["T2D", "T2D", "with DPP"] <- 1 - tpRegress - (tpDeath*tpExcessDeath)
  p_matrix["T2D", "Dead", "with DPP"] <- tpDeath*tpExcessDeath
  p_matrix["Dead", "Dead", "with DPP"] <- 1
  
  return(p_matrix)
}


###########
#Making utility change with age
#state_q_matrix
state_q_matrix_cycle <- function(state_q_matrix, age)
{
  #uNorm 
  uNorm_lookup <-
    c("[30,39]" = rbeta(1,5243.210127, 647.8882745),
      "[40,49]" = rbeta(1,3749.390288, 552.6603545),
      "[50,59]" = rbeta(1,4317.392142, 733.4582261),
      "[60,69]" = rbeta(1,3275.399881, 605.1330164),
      "[70,79]" = rbeta(1,2793.5369, 592.8093771),
      "[80,120]" = rbeta(1,1048.602386, 338.6611273))
  
  age_grpnorm <- cut(age, breaks= c(30,39,49,59,69,79,120))
  uNorm <-uNorm_lookup[age_grpnorm]
  
  #uNDH
  uNDH_lookup <-
    c("[30,39]" = rbeta(1,18852.2607, 3585.449629),
      "[40,49]" = rbeta(1,40397.03381, 8321.148331),
      "[50,59]" = rbeta(1,79459.6746, 19951.65504),
      "[60,69]" = rbeta(1,137085.1892, 32219.07565),
      "[70,79]" = rbeta(1,164314.0338, 37615.0213),
      "[80,120]" = rbeta(1,53445.85673, 14946.63982))
  
  age_grpNDH <- cut(age, breaks= c(30,39,49,59,69,79,120))
  uNDH <-uNDH_lookup[age_grpNDH]
  
  #uT2D
  uT2D_lookup <-
    c("[30,39]" = rbeta(1,27.8648312, 8.05807557),
      "[40,49]" = rbeta(1,48.96668398, 23.6572448),
      "[50,59]" = rbeta(1,158.9763077, 60.70760667),
      "[60,69]" = rbeta(1,238.2666204, 102.8036175),
      "[70,79]" = rbeta(1,294.0078638, 100.8414828),
      "[80,120]" = rbeta(1,198.9254422, 97.37988722))
  
  age_grpT2D <- cut(age, breaks= c(30,39,49,59,69,79,120))
  uT2D <-uT2D_lookup[age_grpT2D]
  
  state_q_matrix[1,1]<- uNorm
  state_q_matrix["without DPP", "NDH"] <- uNDH
  state_q_matrix["without DPP", "T2D"] <- uT2D
  state_q_matrix["without DPP", "Dead"] <-0
  state_q_matrix["with DPP", "Normal"] <- uNorm
  state_q_matrix["with DPP", "NDH"] <- uNDH
  state_q_matrix["with DPP", "T2D"] <-uT2D
  state_q_matrix["with DPP", "Dead"]<-0
  
  
  return(state_q_matrix)
}




## Run model ----
for(n in 1:n_trials) {
  
  for (i in 1:n_treatments) {
    
    age <- Initial_age
    
    for (j in 2:n_cycles) {
      
      p_matrix <- p_matrix_cycle(p_matrix, age, j - 1)
      
      pop[, cycle = j, treatment = i] <-
        pop[, cycle = j - 1, treatment = i] %*% p_matrix[, , treatment = i]
      
      #Transition costs - none included in this model
      trans[, cycle = j, treatment = i] <-
        pop[, cycle = j - 1, treatment = i] %*% (trans_c_matrix * p_matrix[, , treatment = i])
      
      halfpop["Normal", cycle=j, treatment=i]<- 0.5*(pop["Normal", cycle=(j-1), treatment=i]+pop["Normal", cycle=j, treatment=i])
      halfpop["NDH", cycle=j, treatment=i]<- 0.5*(pop["NDH", cycle=(j-1), treatment=i]+pop["NDH", cycle=j, treatment=i])
      halfpop["T2D", cycle=j, treatment=i]<- 0.5*(pop["T2D", cycle=(j-1), treatment=i]+pop["T2D", cycle=j, treatment=i])
      halfpop["Dead", cycle=j, treatment=i]<- 0.5*(pop["Dead", cycle=(j-1), treatment=i]+pop["Dead", cycle=j, treatment=i])
      
      cycle_qalys[i, ] <-
        (state_q_matrix[treatment = i, ] %*% halfpop[, , treatment = i]) * (1/(1 + cDr))^(-1:(n_cycles-2))
      
      age <- age + 1
    }
    
    #Applying half cycle correction
    cycle_state_costs[i, ] <-
      (state_c_matrix[treatment = i, ] %*% halfpop[, , treatment = i]) * (1/(1 + cDr))^(-1:(n_cycles-2)) 
    
    cycle_costs[i, ] <- cycle_state_costs[i, ] #+ cycle_trans_costs[i, ]
    
    total_costs[i] <- sum(cycle_costs[treatment = i, -1])
    total_QALYs[i] <- sum(cycle_qalys[treatment = i, -1])
  }
  #adding in costs of DPP participation
  total_costs["with DPP"] <- total_costs["with DPP"] + dppcost_sum
  
  #adding in benefits from DPP participation
  #multipling by 0.5 as this is utilty gain - need to convert to QALY
  total_QALYs["with DPP"] <- total_QALYs["with DPP"] + 0.5*utility_sum
  
  
  costs[n, ] <- total_costs
  qalys[n,] <- total_QALYs
}

agecat1 <- data.frame(costs, qalys)
names(agecat1)[names(agecat1) == "without.DPP.1"] <- "Without.DPP.qalys.AGECAT1"
names(agecat1)[names(agecat1) == "with.DPP.1"] <- "With.DPP.qalys.AGECAT1"
names(agecat1)[names(agecat1) == "without.DPP"] <- "Without.DPP.costs.AGECAT1"
names(agecat1)[names(agecat1) == "with.DPP"] <- "With.DPP.costs.AGECAT1"
writexl::write_xlsx(agecat1, "age category1 - 58.xlsx")
save(agecat1,file="agecat1.Rda")

#example population trace
popagecat1 <- data.frame(pop)
writexl::write_xlsx(popagecat1, "age category1 trace.xlsx")
save(popagecat1,file="agecat1trace.Rda")

#######
# AGE CAT 2
#######
n_cohort <- 122
n_cycles <- 36
Initial_age <- 46


#NUMBER OF PSA ITERATIONS
n_trials<-10000
costs <- matrix(NA,nrow=n_trials,ncol=n_treatments, dimnames=list(NULL,t_names))
qalys <- matrix(NA,nrow=n_trials,ncol=n_treatments, dimnames=list(NULL,t_names))


#Utility of states
# uNorm <- 0.8395923 / uNDH <- 0.8085401 / uT2D <- 0.708957 / uDeath <-0 / uDPP<-0 

#Cost of states
cNorm <- rgamma(1, 44.44444, scale=45.119475)
cNDH <- rgamma(1,  44.44444, scale=50.0571)
cT2D <- rgamma(1,  44.44444, scale=99.4628)
cDeath <- 0
cDPP <- 0
#Utility of states
uNorm <- rbeta(1,27.53076303,4.883160879)
uNDH <- rbeta(1,26.33007731, 6.334144666)
uT2D <- rbeta(1, 26.52841, 10.65937)
uDeath <-0
uDPP<-0 
#Discount rates
oDr <- 0.035
cDr <- 0.035


#Transition probabilities
#normal to NDH
tpProgNDH <- 0.074202731
#dying from NGT or NDH
tpDeath <- 0.05
#Progression to T2D from NDH
tpProgT2D <- 0.0249
#From NDH to NGT
tpNorm <- 0.0795
#From T2D to NDH
tpRegress <- 0.002796084
#Additional risk of death for T2D
tpExcessDeath <- 1.6

#Effectiveness of the DPP
effect <- 0.8

#utility gains from DPP
utility <- c(1-rbeta(n_cohort,1653302.733,10122.60857))
utility_sum <-sum(utility)

#costs of the DPP
dppcost <- c(rgamma(n_cohort,16.73964,scale=8.469118))
dppcost_sum <- sum(dppcost)

# cost of staying in state

state_c_matrix <-
  matrix(c(cNorm, cNDH, cT2D, 0,
           cNorm, cNDH + cDPP, cT2D, 0),
         byrow = TRUE,
         nrow = n_treatments,
         dimnames = list(t_names,
                         s_names))

# qaly when staying in state
state_q_matrix <-
  matrix(c(uNorm, uNDH, uT2D, 0,
           uNorm, uNDH + uDPP, uT2D, 0),
         byrow = TRUE,
         nrow = n_treatments,
         dimnames = list(t_names,
                         s_names))

# cost of moving to a state
# same for both treatments
#only cost associated with death - no cost in this model
trans_c_matrix <-
  matrix(c(0, 0, 0, cDeath,
           0, 0, 0, cDeath,
           0, 0, 0, cDeath,
           0, 0, 0, 0),
         byrow = TRUE,
         nrow = n_states,
         dimnames = list(from = s_names,
                         to = s_names))

# Transition probabilities ---- 

# time-homogeneous
#not by row - by column
p_matrix <- array(data = c(1-tpProgNDH - tpDeath, tpNorm, 0, 0,
                           tpProgNDH, 1-tpProgT2D-tpDeath-tpNorm, tpRegress, 0,
                           0, tpProgT2D,  1-tpRegress - (tpDeath*tpExcessDeath), 0,
                           tpDeath, tpDeath, (tpDeath*tpExcessDeath), 1,
                           #below would be adjusted for effect of DPP
                           1-tpProgNDH - tpDeath, tpNorm, 0, 0,
                           tpProgNDH, 1-(tpProgT2D*effect)-tpDeath-tpNorm, tpRegress, 0,
                           0, (tpProgT2D*effect),  1-tpRegress - (tpDeath*tpExcessDeath), 0,
                           tpDeath, tpDeath, (tpDeath*tpExcessDeath), 1),
                  dim = c(n_states, n_states, n_treatments),
                  dimnames = list(from = s_names,
                                  to = s_names,
                                  t_names))

# Store population output for each cycle 

# state populations
pop <- array(data = NA,
             dim = c(n_states, n_cycles, n_treatments),
             dimnames = list(state = s_names,
                             cycle = NULL,
                             treatment = t_names))

pop["Normal", cycle = 1, ] <-0
pop["NDH", cycle = 1, ] <- n_cohort
pop["T2D", cycle = 1, ] <- 0
pop["Dead", cycle = 1, ] <- 0

halfpop<-array(data=NA, dim=c(n_states,n_cycles,n_treatments), 
               dimnames = list(state= s_names, cycle=NULL, treatment= t_names))

halfpop["Normal", cycle=1, treatment="without DPP"]<- 0
halfpop["NDH", cycle=1, treatment="without DPP"]<- 0
halfpop["T2D", cycle=1, treatment="without DPP"]<- 0
halfpop["Dead", cycle=1, treatment="without DPP"]<- 0
halfpop["Normal", cycle=1, treatment="with DPP"]<- 0
halfpop["NDH", cycle=1, treatment="with DPP"]<- 0
halfpop["T2D", cycle=1, treatment="with DPP"]<- 0
halfpop["Dead", cycle=1, treatment="with DPP"]<- 0

# _arrived_ state populations
trans <- array(data = NA,
               dim = c(n_states, n_cycles, n_treatments),
               dimnames = list(state = s_names,
                               cycle = NULL,
                               treatment = t_names))

trans[, cycle = 1, ] <- 0


# Sum costs and QALYs for each cycle at a time for each treatment arm

cycle_empty_array <-
  array(NA,
        dim = c(n_treatments, n_cycles),
        dimnames = list(treatment = t_names,
                        cycle = NULL))

cycle_state_costs <- cycle_trans_costs <- cycle_qalys <- cycle_empty_array 
cycle_costs <- cycle_QALYs <- cycle_empty_array
LE <- LYs <- cycle_empty_array    # life expectancy; life-years
cycle_QALE <- cycle_empty_array   # quality-adjusted life expectancy

total_costs <- setNames(c(NA, NA), t_names)
total_QALYs <- setNames(c(NA, NA), t_names)

# #Transition probabilities
#normal to NDH
tpProgNDH <- 0.074202731
# #dying from NGT or NDH
tpDeath <- 0.05
# #Progression from NDH to T2D
tpProgT2D <- 0.0249
# #From NDH to NGT
tpNorm <- 0.0795
# #From T2D to NDH
tpRegress <- 0.002796084
# #Additional risk of death for T2D
tpExcessDeath <- 1.6
#
# Effectiveness of the DPP
effect <- 0.8

tpNorm <- rbeta(1,335,3709)
tpProgT2D <- rbeta(1,25.27,974.73)
tpRegress <- rbeta(1,2.8,997.2)
tpProgNDH <- rbeta(1,77.1,922.9)
tpExcessDeath <- rlnorm(1, meanlog = 0.470004, sdlog = 0.064111)

# # Time-dependent probability matrix ----
# Time dependent transition - risk of death
p_matrix_cycle <- function(p_matrix, age, cycle,
                           tpProgNDH = rbeta(1,77.1,922.9),
                           #tpProgT2D = rbeta(1,25.27,974.73),
                           tpNorm = rbeta(1,335,3709),
                           tpRegress = rbeta(1,2.8,997.2),
                           tpExcessDeath = rlnorm(1, meanlog = 0.470004, sdlog = 0.064111)
) {
  tpDeath_lookup <-
    c("[30,34]" = 0.000674,
      "[35,39]" = 0.001007,
      "[40,44]" = 0.001499,
      "[45,49]" = 0.002298,
      "[50,54]" = 0.003359,
      "[55,59]" = 0.005056,
      "[60,64]" = 0.007940,
      "[65,69]" = 0.012389,
      "[70,74]" = 0.019640,
      "[75,79]" = 0.034483,
      "[80,84]" = 0.060938,
      "[85,89]" = 0.110627,
      "[90,94]" = 0.187364,
      "[95,120]" = 0.295984)
  
  age_grp <- cut(age, breaks = c(30,34,39,44,49,54,59,64,69,74,79, 84, 89, 94, 120))
  tpDeath <- tpDeath_lookup[age_grp]
  
  tpProgT2D_lookup <-
    c("[30,39]" = rnorm(1,mean=0.108, sd=0.0095),
      "[40,49]" = rnorm(1,mean=0.108, sd=0.0095),
      "[50,59]" = rnorm(1,mean=0.108, sd=0.0095),
      "[60,69]" = rnorm(1,mean=0.108, sd=0.0095),
      "[70,79]" = rnorm(1,mean=0.108, sd=0.0095),
      "[80,120]" = rnorm(1,mean=0.108, sd=0.0095))
  
  age_grpT2D <- cut(age, breaks= c(30,39,49,59,69,79,120))
  tpProgT2D <-tpProgT2D_lookup[age_grpT2D]
  
  # [ includes, ( doesn't include
  effect_lookup <-
    c("(0,3]" = rlnorm(1, meanlog = -0.22314, sdlog = 0.044757),
      "(3,10]" = rnorm(1, mean = 1.000, sd = 0.00000001),
      "(10,15]"= rnorm(1, mean = 1.000, sd = 0.00000001),
      "(15,50]" = rnorm(1, mean = 1.000, sd = 0.00000001))
  
  cycle2 <- cut(cycle, breaks = c(0,3,10,15,50))
  
  effect <- effect_lookup[cycle2]
  
  #########
  
  p_matrix["Normal", "Normal", "without DPP"] <- 1 - tpProgNDH - tpDeath
  p_matrix["Normal", "NDH", "without DPP"] <- tpProgNDH 
  p_matrix["Normal", "Dead", "without DPP"] <- tpDeath
  p_matrix["NDH", "Normal", "without DPP"] <- tpNorm
  p_matrix["NDH", "NDH", "without DPP"] <-1-tpProgT2D-tpDeath-tpNorm
  p_matrix["NDH", "T2D", "without DPP"] <- tpProgT2D
  p_matrix["NDH", "Dead", "without DPP"] <- tpDeath
  p_matrix["T2D", "NDH", "without DPP"] <- tpRegress
  p_matrix["T2D", "T2D", "without DPP"] <- 1-tpRegress - (tpDeath*tpExcessDeath)
  p_matrix["T2D", "Dead", "without DPP"] <- (tpDeath)*tpExcessDeath
  p_matrix["Dead", "Dead", "without DPP"] <- 1
  
  p_matrix["Normal", "Normal", "with DPP"] <- 1 - tpProgNDH - tpDeath
  p_matrix["Normal", "NDH", "with DPP"] <- tpProgNDH
  p_matrix["Normal", "Dead", "with DPP"] <- tpDeath
  p_matrix["NDH", "Normal", "with DPP"] <- tpNorm
  p_matrix["NDH", "NDH", "with DPP"] <- 1 - (tpProgT2D*effect) - tpDeath - tpNorm
  p_matrix["NDH", "T2D", "with DPP"] <- (tpProgT2D*effect)
  p_matrix["NDH", "Dead", "with DPP"] <- tpDeath
  p_matrix["T2D", "NDH", "with DPP"] <- tpRegress
  p_matrix["T2D", "T2D", "with DPP"] <- 1 - tpRegress - (tpDeath*tpExcessDeath)
  p_matrix["T2D", "Dead", "with DPP"] <- tpDeath*tpExcessDeath
  p_matrix["Dead", "Dead", "with DPP"] <- 1
  
  return(p_matrix)
}


###########
#Making utility change with age
#state_q_matrix
state_q_matrix_cycle <- function(state_q_matrix, age)
{
  #uNorm 
  uNorm_lookup <-
    c("[30,39]" = rbeta(1,5243.210127, 647.8882745),
      "[40,49]" = rbeta(1,3749.390288, 552.6603545),
      "[50,59]" = rbeta(1,4317.392142, 733.4582261),
      "[60,69]" = rbeta(1,3275.399881, 605.1330164),
      "[70,79]" = rbeta(1,2793.5369, 592.8093771),
      "[80,120]" = rbeta(1,1048.602386, 338.6611273))
  
  
  age_grpnorm <- cut(age, breaks= c(30,39,49,59,69,79,120))
  uNorm <-uNorm_lookup[age_grpnorm]
  
  #uNDH
  uNDH_lookup <-
    c("[30,39]" = rbeta(1,18852.2607, 3585.449629),
      "[40,49]" = rbeta(1,40397.03381, 8321.148331),
      "[50,59]" = rbeta(1,79459.6746, 19951.65504),
      "[60,69]" = rbeta(1,137085.1892, 32219.07565),
      "[70,79]" = rbeta(1,164314.0338, 37615.0213),
      "[80,120]" = rbeta(1,53445.85673, 14946.63982))
  
  
  age_grpNDH <- cut(age, breaks= c(30,39,49,59,69,79,120))
  uNDH <-uNDH_lookup[age_grpNDH]
  
  #uT2D
  uT2D_lookup <-
    c("[30,39]" = rbeta(1,27.8648312, 8.05807557),
      "[40,49]" = rbeta(1,48.96668398, 23.6572448),
      "[50,59]" = rbeta(1,158.9763077, 60.70760667),
      "[60,69]" = rbeta(1,238.2666204, 102.8036175),
      "[70,79]" = rbeta(1,294.0078638, 100.8414828),
      "[80,120]" = rbeta(1,198.9254422, 97.37988722))
  
  
  age_grpT2D <- cut(age, breaks= c(30,39,49,59,69,79,120))
  uT2D <-uT2D_lookup[age_grpT2D]
  
  state_q_matrix[1,1]<- uNorm
  state_q_matrix["without DPP", "NDH"] <- uNDH
  state_q_matrix["without DPP", "T2D"] <- uT2D
  state_q_matrix["without DPP", "Dead"] <-0
  state_q_matrix["with DPP", "Normal"] <- uNorm
  state_q_matrix["with DPP", "NDH"] <- uNDH
  state_q_matrix["with DPP", "T2D"] <-uT2D
  state_q_matrix["with DPP", "Dead"]<-0
  
  
  return(state_q_matrix)
}




## Run model ----
for(n in 1:n_trials) {
  
  for (i in 1:n_treatments) {
    
    age <- Initial_age
    
    for (j in 2:n_cycles) {
      
      p_matrix <- p_matrix_cycle(p_matrix, age, j - 1)
      
      pop[, cycle = j, treatment = i] <-
        pop[, cycle = j - 1, treatment = i] %*% p_matrix[, , treatment = i]
      
      #Transition costs - none included in this model
      trans[, cycle = j, treatment = i] <-
        pop[, cycle = j - 1, treatment = i] %*% (trans_c_matrix * p_matrix[, , treatment = i])
      
      halfpop["Normal", cycle=j, treatment=i]<- 0.5*(pop["Normal", cycle=(j-1), treatment=i]+pop["Normal", cycle=j, treatment=i])
      halfpop["NDH", cycle=j, treatment=i]<- 0.5*(pop["NDH", cycle=(j-1), treatment=i]+pop["NDH", cycle=j, treatment=i])
      halfpop["T2D", cycle=j, treatment=i]<- 0.5*(pop["T2D", cycle=(j-1), treatment=i]+pop["T2D", cycle=j, treatment=i])
      halfpop["Dead", cycle=j, treatment=i]<- 0.5*(pop["Dead", cycle=(j-1), treatment=i]+pop["Dead", cycle=j, treatment=i])
      
      cycle_qalys[i, ] <-
        (state_q_matrix[treatment = i, ] %*% halfpop[, , treatment = i]) * (1/(1 + cDr))^(-1:(n_cycles-2))
      
      age <- age + 1
    }
    
    #Applying half cycle correction
    cycle_state_costs[i, ] <-
      (state_c_matrix[treatment = i, ] %*% halfpop[, , treatment = i]) * (1/(1 + cDr))^(-1:(n_cycles-2)) 
    
    cycle_costs[i, ] <- cycle_state_costs[i, ] #+ cycle_trans_costs[i, ]
    
    total_costs[i] <- sum(cycle_costs[treatment = i, -1])
    total_QALYs[i] <- sum(cycle_qalys[treatment = i, -1])
  }
  #adding in costs of DPP participation
  total_costs["with DPP"] <- total_costs["with DPP"] + dppcost_sum
  
  #adding in benefits from DPP participation
  #multipling by 0.5 as this is utilty gain - need to convert to QALY
  #agecat1: rbeta(1,35.52768105,0.153650952)
  #agecat2: rbeta(1,34.56207, 0.211612)
  #agecat3: rbeta(1,34.18377, 0.307204)
  #agecat4: rbeta(1,36.65039, 0.515016)
  #agecat5: rbeta(1,37.76731, 0.58579)
  #agecat6: rbeta(1, 34.39323, 0.371137)
  
  total_QALYs["with DPP"] <- total_QALYs["with DPP"] + 0.5*utility_sum
  
  
  costs[n, ] <- total_costs
  qalys[n,] <- total_QALYs
}



agecat2 <- data.frame(costs, qalys)
names(agecat2)[names(agecat2) == "without.DPP.1"] <- "Without.DPP.qalys.AGECAT2"
names(agecat2)[names(agecat2) == "with.DPP.1"] <- "With.DPP.qalys.AGECAT2"
names(agecat2)[names(agecat2) == "without.DPP"] <- "Without.DPP.costs.AGECAT2"
names(agecat2)[names(agecat2) == "with.DPP"] <- "With.DPP.costs.AGECAT2"
writexl::write_xlsx(agecat2, "age category2 - 122.xlsx")
save(agecat2,file="agecat2.Rda")

popagecat2 <- data.frame(pop)
writexl::write_xlsx(popagecat2, "age category2 trace.xlsx")
save(popagecat2,file="agecat2trace.Rda")


#######
# AGE CAT 3
#######
n_cohort <- 219
n_cycles <- 36
Initial_age <- 56


#NUMBER OF PSA ITERATIONS
n_trials<-10000
costs <- matrix(NA,nrow=n_trials,ncol=n_treatments, dimnames=list(NULL,t_names))
qalys <- matrix(NA,nrow=n_trials,ncol=n_treatments, dimnames=list(NULL,t_names))


#Utility of states
# uNorm <- 0.8395923 / uNDH <- 0.8085401 / uT2D <- 0.708957 / uDeath <-0 / uDPP<-0 

#Cost of states
cNorm <- rgamma(1, 44.44444, scale=45.119475)
cNDH <- rgamma(1,  44.44444, scale=50.0571)
cT2D <- rgamma(1,  44.44444, scale=99.4628)
cDeath <- 0
cDPP <- 0
#Utility of states
uNorm <- rbeta(1,27.53076303,4.883160879)
uNDH <- rbeta(1,26.33007731, 6.334144666)
uT2D <- rbeta(1, 26.52841, 10.65937)
uDeath <-0
uDPP<-0 
#Discount rates
oDr <- 0.035
cDr <- 0.035


#Transition probabilities
#normal to NDH
tpProgNDH <- 0.074202731
#dying from NGT or NDH
tpDeath <- 0.05
#Progression to T2D from NDH
tpProgT2D <- 0.0249
#From NDH to NGT
tpNorm <- 0.0795
#From T2D to NDH
tpRegress <- 0.002796084
#Additional risk of death for T2D
tpExcessDeath <- 1.6

#Effectiveness of the DPP
effect <- 0.8


#utility gains from DPP
utility <- c(1-rbeta(n_cohort,2945489.923, 26470.65851))
utility_sum <-sum(utility)

#costs of the DPP
dppcost <- c(rgamma(n_cohort,16.73964,scale=8.469118))
dppcost_sum <- sum(dppcost)

# cost of staying in state

state_c_matrix <-
  matrix(c(cNorm, cNDH, cT2D, 0,
           cNorm, cNDH + cDPP, cT2D, 0),
         byrow = TRUE,
         nrow = n_treatments,
         dimnames = list(t_names,
                         s_names))

# qaly when staying in state
state_q_matrix <-
  matrix(c(uNorm, uNDH, uT2D, 0,
           uNorm, uNDH + uDPP, uT2D, 0),
         byrow = TRUE,
         nrow = n_treatments,
         dimnames = list(t_names,
                         s_names))

# cost of moving to a state
# same for both treatments
#only cost associated with death - no cost in this model
trans_c_matrix <-
  matrix(c(0, 0, 0, cDeath,
           0, 0, 0, cDeath,
           0, 0, 0, cDeath,
           0, 0, 0, 0),
         byrow = TRUE,
         nrow = n_states,
         dimnames = list(from = s_names,
                         to = s_names))

# Transition probabilities ---- 

# time-homogeneous
#not by row - by column
p_matrix <- array(data = c(1-tpProgNDH - tpDeath, tpNorm, 0, 0,
                           tpProgNDH, 1-tpProgT2D-tpDeath-tpNorm, tpRegress, 0,
                           0, tpProgT2D,  1-tpRegress - (tpDeath*tpExcessDeath), 0,
                           tpDeath, tpDeath, (tpDeath*tpExcessDeath), 1,
                           #below would be adjusted for effect of DPP
                           1-tpProgNDH - tpDeath, tpNorm, 0, 0,
                           tpProgNDH, 1-(tpProgT2D*effect)-tpDeath-tpNorm, tpRegress, 0,
                           0, (tpProgT2D*effect),  1-tpRegress - (tpDeath*tpExcessDeath), 0,
                           tpDeath, tpDeath, (tpDeath*tpExcessDeath), 1),
                  dim = c(n_states, n_states, n_treatments),
                  dimnames = list(from = s_names,
                                  to = s_names,
                                  t_names))

# Store population output for each cycle 

# state populations
pop <- array(data = NA,
             dim = c(n_states, n_cycles, n_treatments),
             dimnames = list(state = s_names,
                             cycle = NULL,
                             treatment = t_names))

pop["Normal", cycle = 1, ] <-0
pop["NDH", cycle = 1, ] <- n_cohort
pop["T2D", cycle = 1, ] <- 0
pop["Dead", cycle = 1, ] <- 0

halfpop<-array(data=NA, dim=c(n_states,n_cycles,n_treatments), 
               dimnames = list(state= s_names, cycle=NULL, treatment= t_names))

halfpop["Normal", cycle=1, treatment="without DPP"]<- 0
halfpop["NDH", cycle=1, treatment="without DPP"]<- 0
halfpop["T2D", cycle=1, treatment="without DPP"]<- 0
halfpop["Dead", cycle=1, treatment="without DPP"]<- 0
halfpop["Normal", cycle=1, treatment="with DPP"]<- 0
halfpop["NDH", cycle=1, treatment="with DPP"]<- 0
halfpop["T2D", cycle=1, treatment="with DPP"]<- 0
halfpop["Dead", cycle=1, treatment="with DPP"]<- 0

# _arrived_ state populations
trans <- array(data = NA,
               dim = c(n_states, n_cycles, n_treatments),
               dimnames = list(state = s_names,
                               cycle = NULL,
                               treatment = t_names))

trans[, cycle = 1, ] <- 0


# Sum costs and QALYs for each cycle at a time for each treatment arm

cycle_empty_array <-
  array(NA,
        dim = c(n_treatments, n_cycles),
        dimnames = list(treatment = t_names,
                        cycle = NULL))

cycle_state_costs <- cycle_trans_costs <- cycle_qalys <- cycle_empty_array 
cycle_costs <- cycle_QALYs <- cycle_empty_array
LE <- LYs <- cycle_empty_array    # life expectancy; life-years
cycle_QALE <- cycle_empty_array   # quality-adjusted life expectancy

total_costs <- setNames(c(NA, NA), t_names)
total_QALYs <- setNames(c(NA, NA), t_names)

# #Transition probabilities
#normal to NDH
tpProgNDH <- 0.074202731
# #dying from NGT or NDH
tpDeath <- 0.05
# #Progression from NDH to T2D
tpProgT2D <- 0.0249
# #From NDH to NGT
tpNorm <- 0.0795
# #From T2D to NDH
tpRegress <- 0.002796084
# #Additional risk of death for T2D
tpExcessDeath <- 1.6
#
#Effectiveness of the DPP
effect <- 0.8

tpNorm <- rbeta(1,335,3709)
tpProgT2D <- rbeta(1,25.27,974.73)
tpRegress <- rbeta(1,2.8,997.2)
tpProgNDH <- rbeta(1,77.1,922.9)
tpExcessDeath <- rlnorm(1, meanlog = 0.470004, sdlog = 0.064111)

# # Time-dependent probability matrix ----
# Time dependent transition - risk of death
p_matrix_cycle <- function(p_matrix, age, cycle,
                           tpProgNDH = rbeta(1,77.1,922.9),
                           #tpProgT2D = rbeta(1,25.27,974.73),
                           tpNorm = rbeta(1,335,3709),
                           tpRegress = rbeta(1,2.8,997.2),
                           tpExcessDeath = rlnorm(1, meanlog = 0.470004, sdlog = 0.064111)
) {
  tpDeath_lookup <-
    c("[30,34]" = 0.000674,
      "[35,39]" = 0.001007,
      "[40,44]" = 0.001499,
      "[45,49]" = 0.002298,
      "[50,54]" = 0.003359,
      "[55,59]" = 0.005056,
      "[60,64]" = 0.007940,
      "[65,69]" = 0.012389,
      "[70,74]" = 0.019640,
      "[75,79]" = 0.034483,
      "[80,84]" = 0.060938,
      "[85,89]" = 0.110627,
      "[90,94]" = 0.187364,
      "[95,120]" = 0.295984)
  
  age_grp <- cut(age, breaks = c(30,34,39,44,49,54,59,64,69,74,79, 84, 89, 94, 120))
  tpDeath <- tpDeath_lookup[age_grp]
  
  
  tpProgT2D_lookup <-
    c("[30,39]" = rnorm(1,mean=0.108, sd=0.0095),
      "[40,49]" = rnorm(1,mean=0.108, sd=0.0095),
      "[50,59]" = rnorm(1,mean=0.108, sd=0.0095),
      "[60,69]" = rnorm(1,mean=0.108, sd=0.0095),
      "[70,79]" = rnorm(1,mean=0.108, sd=0.0095),
      "[80,120]" = rnorm(1,mean=0.108, sd=0.0095))
  
  age_grpT2D <- cut(age, breaks= c(30,39,49,59,69,79,120))
  tpProgT2D <-tpProgT2D_lookup[age_grpT2D]
  
  # [ includes, ( doesn't include
  effect_lookup <-
    c("(0,3]" = rlnorm(1, meanlog = -0.22314, sdlog = 0.044757),
      "(3,10]" = rnorm(1, mean = 1.000, sd = 0.00000001),
      "(10,15]"= rnorm(1, mean = 1.000, sd = 0.00000001),
      "(15,50]" = rnorm(1, mean = 1.000, sd = 0.00000001))
  
  cycle2 <- cut(cycle, breaks = c(0,3,10,15,50))
  
  effect <- effect_lookup[cycle2]
  
  #########
  
  p_matrix["Normal", "Normal", "without DPP"] <- 1 - tpProgNDH - tpDeath
  p_matrix["Normal", "NDH", "without DPP"] <- tpProgNDH 
  p_matrix["Normal", "Dead", "without DPP"] <- tpDeath
  p_matrix["NDH", "Normal", "without DPP"] <- tpNorm
  p_matrix["NDH", "NDH", "without DPP"] <-1-tpProgT2D-tpDeath-tpNorm
  p_matrix["NDH", "T2D", "without DPP"] <- tpProgT2D
  p_matrix["NDH", "Dead", "without DPP"] <- tpDeath
  p_matrix["T2D", "NDH", "without DPP"] <- tpRegress
  p_matrix["T2D", "T2D", "without DPP"] <- 1-tpRegress - (tpDeath*tpExcessDeath)
  p_matrix["T2D", "Dead", "without DPP"] <- (tpDeath)*tpExcessDeath
  p_matrix["Dead", "Dead", "without DPP"] <- 1
  
  p_matrix["Normal", "Normal", "with DPP"] <- 1 - tpProgNDH - tpDeath
  p_matrix["Normal", "NDH", "with DPP"] <- tpProgNDH
  p_matrix["Normal", "Dead", "with DPP"] <- tpDeath
  p_matrix["NDH", "Normal", "with DPP"] <- tpNorm
  p_matrix["NDH", "NDH", "with DPP"] <- 1 - (tpProgT2D*effect) - tpDeath - tpNorm
  p_matrix["NDH", "T2D", "with DPP"] <- (tpProgT2D*effect)
  p_matrix["NDH", "Dead", "with DPP"] <- tpDeath
  p_matrix["T2D", "NDH", "with DPP"] <- tpRegress
  p_matrix["T2D", "T2D", "with DPP"] <- 1 - tpRegress - (tpDeath*tpExcessDeath)
  p_matrix["T2D", "Dead", "with DPP"] <- tpDeath*tpExcessDeath
  p_matrix["Dead", "Dead", "with DPP"] <- 1
  
  return(p_matrix)
}


###########
#Making utility change with age
#state_q_matrix
state_q_matrix_cycle <- function(state_q_matrix, age)
{
  #uNorm 
  uNorm_lookup <-
    c("[30,39]" = rbeta(1,5243.210127, 647.8882745),
      "[40,49]" = rbeta(1,3749.390288, 552.6603545),
      "[50,59]" = rbeta(1,4317.392142, 733.4582261),
      "[60,69]" = rbeta(1,3275.399881, 605.1330164),
      "[70,79]" = rbeta(1,2793.5369, 592.8093771),
      "[80,120]" = rbeta(1,1048.602386, 338.6611273))
  
  
  age_grpnorm <- cut(age, breaks= c(30,39,49,59,69,79,120))
  uNorm <-uNorm_lookup[age_grpnorm]
  
  #uNDH
  uNDH_lookup <-
    c("[30,39]" = rbeta(1,18852.2607, 3585.449629),
      "[40,49]" = rbeta(1,40397.03381, 8321.148331),
      "[50,59]" = rbeta(1,79459.6746, 19951.65504),
      "[60,69]" = rbeta(1,137085.1892, 32219.07565),
      "[70,79]" = rbeta(1,164314.0338, 37615.0213),
      "[80,120]" = rbeta(1,53445.85673, 14946.63982))
  
  
  age_grpNDH <- cut(age, breaks= c(30,39,49,59,69,79,120))
  uNDH <-uNDH_lookup[age_grpNDH]
  
  ##uT2D
  uT2D_lookup <-
    c("[30,39]" = rbeta(1,27.8648312, 8.05807557),
      "[40,49]" = rbeta(1,48.96668398, 23.6572448),
      "[50,59]" = rbeta(1,158.9763077, 60.70760667),
      "[60,69]" = rbeta(1,238.2666204, 102.8036175),
      "[70,79]" = rbeta(1,294.0078638, 100.8414828),
      "[80,120]" = rbeta(1,198.9254422, 97.37988722))
  
  
  age_grpT2D <- cut(age, breaks= c(30,39,49,59,69,79,120))
  uT2D <-uT2D_lookup[age_grpT2D]
  
  state_q_matrix[1,1]<- uNorm
  state_q_matrix["without DPP", "NDH"] <- uNDH
  state_q_matrix["without DPP", "T2D"] <- uT2D
  state_q_matrix["without DPP", "Dead"] <-0
  state_q_matrix["with DPP", "Normal"] <- uNorm
  state_q_matrix["with DPP", "NDH"] <- uNDH
  state_q_matrix["with DPP", "T2D"] <-uT2D
  state_q_matrix["with DPP", "Dead"]<-0
  
  
  return(state_q_matrix)
}




## Run model ----
for(n in 1:n_trials) {
  
  for (i in 1:n_treatments) {
    
    age <- Initial_age
    
    for (j in 2:n_cycles) {
      
      p_matrix <- p_matrix_cycle(p_matrix, age, j - 1)
      
      pop[, cycle = j, treatment = i] <-
        pop[, cycle = j - 1, treatment = i] %*% p_matrix[, , treatment = i]
      
      #Transition costs - none included in this model
      trans[, cycle = j, treatment = i] <-
        pop[, cycle = j - 1, treatment = i] %*% (trans_c_matrix * p_matrix[, , treatment = i])
      
      halfpop["Normal", cycle=j, treatment=i]<- 0.5*(pop["Normal", cycle=(j-1), treatment=i]+pop["Normal", cycle=j, treatment=i])
      halfpop["NDH", cycle=j, treatment=i]<- 0.5*(pop["NDH", cycle=(j-1), treatment=i]+pop["NDH", cycle=j, treatment=i])
      halfpop["T2D", cycle=j, treatment=i]<- 0.5*(pop["T2D", cycle=(j-1), treatment=i]+pop["T2D", cycle=j, treatment=i])
      halfpop["Dead", cycle=j, treatment=i]<- 0.5*(pop["Dead", cycle=(j-1), treatment=i]+pop["Dead", cycle=j, treatment=i])
      
      cycle_qalys[i, ] <-
        (state_q_matrix[treatment = i, ] %*% halfpop[, , treatment = i]) * (1/(1 + cDr))^(-1:(n_cycles-2))
      
      age <- age + 1
    }
    
    #Applying half cycle correction
    cycle_state_costs[i, ] <-
      (state_c_matrix[treatment = i, ] %*% halfpop[, , treatment = i]) * (1/(1 + cDr))^(-1:(n_cycles-2)) 
    
    cycle_costs[i, ] <- cycle_state_costs[i, ] #+ cycle_trans_costs[i, ]
    total_costs[i] <- sum(cycle_costs[treatment = i, -1])
    total_QALYs[i] <- sum(cycle_qalys[treatment = i, -1])
  }
  #adding in costs of DPP participation
  total_costs["with DPP"] <- total_costs["with DPP"] + dppcost_sum
  
  #adding in benefits from DPP participation
  #multipling by 0.5 as this is utilty gain - need to convert to QALY
  #agecat1: rbeta(1,35.52768105,0.153650952)
  #agecat2: rbeta(1,34.56207, 0.211612)
  #agecat3: rbeta(1,34.18377, 0.307204)
  #agecat4: rbeta(1,36.65039, 0.515016)
  #agecat5: rbeta(1,37.76731, 0.58579)
  #agecat6: rbeta(1, 34.39323, 0.371137)
  
  total_QALYs["with DPP"] <- total_QALYs["with DPP"] + 0.5*utility_sum
  
  
  costs[n, ] <- total_costs
  qalys[n,] <- total_QALYs
}



agecat3 <- data.frame(costs, qalys)
names(agecat3)[names(agecat3) == "without.DPP.1"] <- "Without.DPP.qalys.AGECAT3"
names(agecat3)[names(agecat3) == "with.DPP.1"] <- "With.DPP.qalys.AGECAT3"
names(agecat3)[names(agecat3) == "without.DPP"] <- "Without.DPP.costs.AGECAT3"
names(agecat3)[names(agecat3) == "with.DPP"] <- "With.DPP.costs.AGECAT3"
writexl::write_xlsx(agecat3, "age category3 - 219.xlsx")
save(agecat3,file="agecat3.Rda")

popagecat3 <- data.frame(pop)
writexl::write_xlsx(popagecat3, "age category3 trace.xlsx")
save(popagecat3,file="agecat3trace.Rda")


###############################
### AGE CAT 4
###############################
n_cohort <- 267
n_cycles <- 36
Initial_age <- 66


#NUMBER OF PSA ITERATIONS
n_trials<-10000
costs <- matrix(NA,nrow=n_trials,ncol=n_treatments, dimnames=list(NULL,t_names))
qalys <- matrix(NA,nrow=n_trials,ncol=n_treatments, dimnames=list(NULL,t_names))


#Utility of states
# uNorm <- 0.8395923 / uNDH <- 0.8085401 / uT2D <- 0.708957 / uDeath <-0 / uDPP<-0 

#Cost of states
cNorm <- rgamma(1, 44.44444, scale=45.119475)
cNDH <- rgamma(1,  44.44444, scale=50.0571)
cT2D <- rgamma(1,  44.44444, scale=99.4628)
cDeath <- 0
cDPP <- 0
#Utility of states
uNorm <- rbeta(1,27.53076303,4.883160879)
uNDH <- rbeta(1,26.33007731, 6.334144666)
uT2D <- rbeta(1, 26.52841, 10.65937)
uDeath <-0
uDPP<-0 
#Discount rates
oDr <- 0.035
cDr <- 0.035


#Transition probabilities
#normal to NDH
tpProgNDH <- 0.074202731
#dying from NGT or NDH
tpDeath <- 0.05
#Progression to T2D from NDH
tpProgT2D <- 0.0249
#From NDH to NGT
tpNorm <- 0.0795
#From T2D to NDH
tpRegress <- 0.002796084
#Additional risk of death for T2D
tpExcessDeath <- 1.6

#Effectiveness of the DPP
effect <- 0.8

#utility gains from DPP
utility <- c(1-rbeta(n_cohort,3884458.555, 54584.90078))
utility_sum <-sum(utility)

#costs of the DPP
dppcost <- c(rgamma(n_cohort,16.73964,scale=8.469118))
dppcost_sum <- sum(dppcost)

# cost of staying in state

state_c_matrix <-
  matrix(c(cNorm, cNDH, cT2D, 0,
           cNorm, cNDH + cDPP, cT2D, 0),
         byrow = TRUE,
         nrow = n_treatments,
         dimnames = list(t_names,
                         s_names))

# qaly when staying in state
state_q_matrix <-
  matrix(c(uNorm, uNDH, uT2D, 0,
           uNorm, uNDH + uDPP, uT2D, 0),
         byrow = TRUE,
         nrow = n_treatments,
         dimnames = list(t_names,
                         s_names))

# cost of moving to a state
# same for both treatments
#only cost associated with death - no cost in this model
trans_c_matrix <-
  matrix(c(0, 0, 0, cDeath,
           0, 0, 0, cDeath,
           0, 0, 0, cDeath,
           0, 0, 0, 0),
         byrow = TRUE,
         nrow = n_states,
         dimnames = list(from = s_names,
                         to = s_names))

# Transition probabilities ---- 

# time-homogeneous
#not by row - by column
p_matrix <- array(data = c(1-tpProgNDH - tpDeath, tpNorm, 0, 0,
                           tpProgNDH, 1-tpProgT2D-tpDeath-tpNorm, tpRegress, 0,
                           0, tpProgT2D,  1-tpRegress - (tpDeath*tpExcessDeath), 0,
                           tpDeath, tpDeath, (tpDeath*tpExcessDeath), 1,
                           #below would be adjusted for effect of DPP
                           1-tpProgNDH - tpDeath, tpNorm, 0, 0,
                           tpProgNDH, 1-(tpProgT2D*effect)-tpDeath-tpNorm, tpRegress, 0,
                           0, (tpProgT2D*effect),  1-tpRegress - (tpDeath*tpExcessDeath), 0,
                           tpDeath, tpDeath, (tpDeath*tpExcessDeath), 1),
                  dim = c(n_states, n_states, n_treatments),
                  dimnames = list(from = s_names,
                                  to = s_names,
                                  t_names))

# Store population output for each cycle 

# state populations
pop <- array(data = NA,
             dim = c(n_states, n_cycles, n_treatments),
             dimnames = list(state = s_names,
                             cycle = NULL,
                             treatment = t_names))

pop["Normal", cycle = 1, ] <-0
pop["NDH", cycle = 1, ] <- n_cohort
pop["T2D", cycle = 1, ] <- 0
pop["Dead", cycle = 1, ] <- 0

halfpop<-array(data=NA, dim=c(n_states,n_cycles,n_treatments), 
               dimnames = list(state= s_names, cycle=NULL, treatment= t_names))

halfpop["Normal", cycle=1, treatment="without DPP"]<- 0
halfpop["NDH", cycle=1, treatment="without DPP"]<- 0
halfpop["T2D", cycle=1, treatment="without DPP"]<- 0
halfpop["Dead", cycle=1, treatment="without DPP"]<- 0
halfpop["Normal", cycle=1, treatment="with DPP"]<- 0
halfpop["NDH", cycle=1, treatment="with DPP"]<- 0
halfpop["T2D", cycle=1, treatment="with DPP"]<- 0
halfpop["Dead", cycle=1, treatment="with DPP"]<- 0

# _arrived_ state populations
trans <- array(data = NA,
               dim = c(n_states, n_cycles, n_treatments),
               dimnames = list(state = s_names,
                               cycle = NULL,
                               treatment = t_names))

trans[, cycle = 1, ] <- 0


# Sum costs and QALYs for each cycle at a time for each treatment arm

cycle_empty_array <-
  array(NA,
        dim = c(n_treatments, n_cycles),
        dimnames = list(treatment = t_names,
                        cycle = NULL))

cycle_state_costs <- cycle_trans_costs <- cycle_qalys <- cycle_empty_array 
cycle_costs <- cycle_QALYs <- cycle_empty_array
LE <- LYs <- cycle_empty_array    # life expectancy; life-years
cycle_QALE <- cycle_empty_array   # quality-adjusted life expectancy

total_costs <- setNames(c(NA, NA), t_names)
total_QALYs <- setNames(c(NA, NA), t_names)

# #Transition probabilities
#normal to NDH
tpProgNDH <- 0.074202731
# #dying from NGT or NDH
tpDeath <- 0.05
# #Progression from NDH to T2D
tpProgT2D <- 0.0249
# #From NDH to NGT
tpNorm <- 0.0795
# #From T2D to NDH
tpRegress <- 0.002796084
# #Additional risk of death for T2D
tpExcessDeath <- 1.6
#
#Effectiveness of the DPP
effect <- 0.8

tpNorm <- rbeta(1,335,3709)
tpProgT2D <- rbeta(1,25.27,974.73)
tpRegress <- rbeta(1,2.8,997.2)
tpProgNDH <- rbeta(1,77.1,922.9)
tpExcessDeath <- rlnorm(1, meanlog = 0.470004, sdlog = 0.064111)

# # Time-dependent probability matrix ----
# Time dependent transition - risk of death
p_matrix_cycle <- function(p_matrix, age, cycle,
                           tpProgNDH = rbeta(1,77.1,922.9),
                           #tpProgT2D = rbeta(1,25.27,974.73),
                           tpNorm = rbeta(1,335,3709),
                           tpRegress = rbeta(1,2.8,997.2),
                           tpExcessDeath = rlnorm(1, meanlog = 0.470004, sdlog = 0.064111)
) {
  tpDeath_lookup <-
    c("[30,34]" = 0.000674,
      "[35,39]" = 0.001007,
      "[40,44]" = 0.001499,
      "[45,49]" = 0.002298,
      "[50,54]" = 0.003359,
      "[55,59]" = 0.005056,
      "[60,64]" = 0.007940,
      "[65,69]" = 0.012389,
      "[70,74]" = 0.019640,
      "[75,79]" = 0.034483,
      "[80,84]" = 0.060938,
      "[85,89]" = 0.110627,
      "[90,94]" = 0.187364,
      "[95,120]" = 0.295984)
  
  age_grp <- cut(age, breaks = c(30,34,39,44,49,54,59,64,69,74,79, 84, 89, 94, 120))
  tpDeath <- tpDeath_lookup[age_grp]
  
  
  tpProgT2D_lookup <-
    c("[30,39]" = rnorm(1,mean=0.108, sd=0.0095),
      "[40,49]" = rnorm(1,mean=0.108, sd=0.0095),
      "[50,59]" = rnorm(1,mean=0.108, sd=0.0095),
      "[60,69]" = rnorm(1,mean=0.108, sd=0.0095),
      "[70,79]" = rnorm(1,mean=0.108, sd=0.0095),
      "[80,120]" = rnorm(1,mean=0.108, sd=0.0095))
  
  age_grpT2D <- cut(age, breaks= c(30,39,49,59,69,79,120))
  tpProgT2D <-tpProgT2D_lookup[age_grpT2D]
  
  # [ includes, ( doesn't include
  effect_lookup <-
    c("(0,3]" = rlnorm(1, meanlog = -0.22314, sdlog = 0.044757),
      "(3,10]" = rnorm(1, mean = 1.000, sd = 0.00000001),
      "(10,15]"= rnorm(1, mean = 1.000, sd = 0.00000001),
      "(15,50]" = rnorm(1, mean = 1.000, sd = 0.00000001))
  
  cycle2 <- cut(cycle, breaks = c(0,3,10,15,50))
  
  effect <- effect_lookup[cycle2]
  
  #########
  
  p_matrix["Normal", "Normal", "without DPP"] <- 1 - tpProgNDH - tpDeath
  p_matrix["Normal", "NDH", "without DPP"] <- tpProgNDH 
  p_matrix["Normal", "Dead", "without DPP"] <- tpDeath
  p_matrix["NDH", "Normal", "without DPP"] <- tpNorm
  p_matrix["NDH", "NDH", "without DPP"] <-1-tpProgT2D-tpDeath-tpNorm
  p_matrix["NDH", "T2D", "without DPP"] <- tpProgT2D
  p_matrix["NDH", "Dead", "without DPP"] <- tpDeath
  p_matrix["T2D", "NDH", "without DPP"] <- tpRegress
  p_matrix["T2D", "T2D", "without DPP"] <- 1-tpRegress - (tpDeath*tpExcessDeath)
  p_matrix["T2D", "Dead", "without DPP"] <- (tpDeath)*tpExcessDeath
  p_matrix["Dead", "Dead", "without DPP"] <- 1
  
  p_matrix["Normal", "Normal", "with DPP"] <- 1 - tpProgNDH - tpDeath
  p_matrix["Normal", "NDH", "with DPP"] <- tpProgNDH
  p_matrix["Normal", "Dead", "with DPP"] <- tpDeath
  p_matrix["NDH", "Normal", "with DPP"] <- tpNorm
  p_matrix["NDH", "NDH", "with DPP"] <- 1 - (tpProgT2D*effect) - tpDeath - tpNorm
  p_matrix["NDH", "T2D", "with DPP"] <- (tpProgT2D*effect)
  p_matrix["NDH", "Dead", "with DPP"] <- tpDeath
  p_matrix["T2D", "NDH", "with DPP"] <- tpRegress
  p_matrix["T2D", "T2D", "with DPP"] <- 1 - tpRegress - (tpDeath*tpExcessDeath)
  p_matrix["T2D", "Dead", "with DPP"] <- tpDeath*tpExcessDeath
  p_matrix["Dead", "Dead", "with DPP"] <- 1
  
  return(p_matrix)
}


###########
#Making utility change with age
#state_q_matrix
state_q_matrix_cycle <- function(state_q_matrix, age)
{
  #uNorm 
  uNorm_lookup <-
    c("[30,39]" = rbeta(1,5243.210127, 647.8882745),
      "[40,49]" = rbeta(1,3749.390288, 552.6603545),
      "[50,59]" = rbeta(1,4317.392142, 733.4582261),
      "[60,69]" = rbeta(1,3275.399881, 605.1330164),
      "[70,79]" = rbeta(1,2793.5369, 592.8093771),
      "[80,120]" = rbeta(1,1048.602386, 338.6611273))
  
  
  age_grpnorm <- cut(age, breaks= c(30,39,49,59,69,79,120))
  uNorm <-uNorm_lookup[age_grpnorm]
  
  #uNDH
  uNDH_lookup <-
    c("[30,39]" = rbeta(1,18852.2607, 3585.449629),
      "[40,49]" = rbeta(1,40397.03381, 8321.148331),
      "[50,59]" = rbeta(1,79459.6746, 19951.65504),
      "[60,69]" = rbeta(1,137085.1892, 32219.07565),
      "[70,79]" = rbeta(1,164314.0338, 37615.0213),
      "[80,120]" = rbeta(1,53445.85673, 14946.63982))
  
  
  age_grpNDH <- cut(age, breaks= c(30,39,49,59,69,79,120))
  uNDH <-uNDH_lookup[age_grpNDH]
  
  #uT2D
  uT2D_lookup <-
    c("[30,39]" = rbeta(1,27.8648312, 8.05807557),
      "[40,49]" = rbeta(1,48.96668398, 23.6572448),
      "[50,59]" = rbeta(1,158.9763077, 60.70760667),
      "[60,69]" = rbeta(1,238.2666204, 102.8036175),
      "[70,79]" = rbeta(1,294.0078638, 100.8414828),
      "[80,120]" = rbeta(1,198.9254422, 97.37988722))
  
  
  age_grpT2D <- cut(age, breaks= c(30,39,49,59,69,79,120))
  uT2D <-uT2D_lookup[age_grpT2D]
  
  state_q_matrix[1,1]<- uNorm
  state_q_matrix["without DPP", "NDH"] <- uNDH
  state_q_matrix["without DPP", "T2D"] <- uT2D
  state_q_matrix["without DPP", "Dead"] <-0
  state_q_matrix["with DPP", "Normal"] <- uNorm
  state_q_matrix["with DPP", "NDH"] <- uNDH
  state_q_matrix["with DPP", "T2D"] <-uT2D
  state_q_matrix["with DPP", "Dead"]<-0
  
  
  return(state_q_matrix)
}




## Run model ----
for(n in 1:n_trials) {
  
  for (i in 1:n_treatments) {
    
    age <- Initial_age
    
    for (j in 2:n_cycles) {
      
      p_matrix <- p_matrix_cycle(p_matrix, age, j - 1)
      
      pop[, cycle = j, treatment = i] <-
        pop[, cycle = j - 1, treatment = i] %*% p_matrix[, , treatment = i]
      
      #Transition costs - none included in this model
      trans[, cycle = j, treatment = i] <-
        pop[, cycle = j - 1, treatment = i] %*% (trans_c_matrix * p_matrix[, , treatment = i])
      
      halfpop["Normal", cycle=j, treatment=i]<- 0.5*(pop["Normal", cycle=(j-1), treatment=i]+pop["Normal", cycle=j, treatment=i])
      halfpop["NDH", cycle=j, treatment=i]<- 0.5*(pop["NDH", cycle=(j-1), treatment=i]+pop["NDH", cycle=j, treatment=i])
      halfpop["T2D", cycle=j, treatment=i]<- 0.5*(pop["T2D", cycle=(j-1), treatment=i]+pop["T2D", cycle=j, treatment=i])
      halfpop["Dead", cycle=j, treatment=i]<- 0.5*(pop["Dead", cycle=(j-1), treatment=i]+pop["Dead", cycle=j, treatment=i])
      
      cycle_qalys[i, ] <-
        (state_q_matrix[treatment = i, ] %*% halfpop[, , treatment = i]) * (1/(1 + cDr))^(-1:(n_cycles-2))
      
      age <- age + 1
    }
    
    #Applying half cycle correction
    cycle_state_costs[i, ] <-
      (state_c_matrix[treatment = i, ] %*% halfpop[, , treatment = i]) * (1/(1 + cDr))^(-1:(n_cycles-2)) 
    
    cycle_costs[i, ] <- cycle_state_costs[i, ] #+ cycle_trans_costs[i, ]
    
    total_costs[i] <- sum(cycle_costs[treatment = i, -1])
    total_QALYs[i] <- sum(cycle_qalys[treatment = i, -1])
  }
  #adding in costs of DPP participation
  total_costs["with DPP"] <- total_costs["with DPP"] + dppcost_sum
  
  #adding in benefits from DPP participation
  #multipling by 0.5 as this is utilty gain - need to convert to QALY
  #agecat1: rbeta(1,35.52768105,0.153650952)
  #agecat2: rbeta(1,34.56207, 0.211612)
  #agecat3: rbeta(1,34.18377, 0.307204)
  #agecat4: rbeta(1,36.65039, 0.515016)
  #agecat5: rbeta(1,37.76731, 0.58579)
  #agecat6: rbeta(1, 34.39323, 0.371137)
  
  total_QALYs["with DPP"] <- total_QALYs["with DPP"] + 0.5*utility_sum
  
  
  costs[n, ] <- total_costs
  qalys[n,] <- total_QALYs
}


agecat4 <- data.frame(costs, qalys)
names(agecat4)[names(agecat4) == "without.DPP.1"] <- "Without.DPP.qalys.AGECAT4"
names(agecat4)[names(agecat4) == "with.DPP.1"] <- "With.DPP.qalys.AGECAT4"
names(agecat4)[names(agecat4) == "without.DPP"] <- "Without.DPP.costs.AGECAT4"
names(agecat4)[names(agecat4) == "with.DPP"] <- "With.DPP.costs.AGECAT4"
writexl::write_xlsx(agecat4, "age category4 - 267.xlsx")
save(agecat4,file="agecat4.Rda")

popagecat4 <- data.frame(pop)
writexl::write_xlsx(popagecat4, "age category4 trace.xlsx")
save(popagecat4,file="agecat4trace.Rda")


#######
# AGE CAT 5
#######
n_cohort <- 244
n_cycles <- 36
Initial_age <- 76


#NUMBER OF PSA ITERATIONS
n_trials<-10000
costs <- matrix(NA,nrow=n_trials,ncol=n_treatments, dimnames=list(NULL,t_names))
qalys <- matrix(NA,nrow=n_trials,ncol=n_treatments, dimnames=list(NULL,t_names))


#Utility of states
# uNorm <- 0.8395923 / uNDH <- 0.8085401 / uT2D <- 0.708957 / uDeath <-0 / uDPP<-0 


#Cost of states
cNorm <- rgamma(1, 44.44444, scale=45.119475)
cNDH <- rgamma(1,  44.44444, scale=50.0571)
cT2D <- rgamma(1,  44.44444, scale=99.4628)
cDeath <- 0
cDPP <- 0
#Utility of states
uNorm <- rbeta(1,27.53076303,4.883160879)
uNDH <- rbeta(1,26.33007731, 6.334144666)
uT2D <- rbeta(1, 26.52841, 10.65937)
uDeath <-0
uDPP<-0 
#Discount rates
oDr <- 0.035
cDr <- 0.035


#Transition probabilities
#normal to NDH
tpProgNDH <- 0.074202731
#dying from NGT or NDH
tpDeath <- 0.05
#Progression to T2D from NDH
tpProgT2D <- 0.0249
#From NDH to NGT
tpNorm <- 0.0795
#From T2D to NDH
tpRegress <- 0.002796084
#Additional risk of death for T2D
tpExcessDeath <- 1.6

#Effectiveness of the DPP
effect <- 0.8

#utility gains from DPP
utility <- c(1-rbeta(n_cohort,3650006.111, 56613.42414))
utility_sum <-sum(utility)

#costs of the DPP
dppcost <- c(rgamma(n_cohort,16.73964,scale=8.469118))
dppcost_sum <- sum(dppcost)


# cost of staying in state

state_c_matrix <-
  matrix(c(cNorm, cNDH, cT2D, 0,
           cNorm, cNDH + cDPP, cT2D, 0),
         byrow = TRUE,
         nrow = n_treatments,
         dimnames = list(t_names,
                         s_names))

# qaly when staying in state
state_q_matrix <-
  matrix(c(uNorm, uNDH, uT2D, 0,
           uNorm, uNDH + uDPP, uT2D, 0),
         byrow = TRUE,
         nrow = n_treatments,
         dimnames = list(t_names,
                         s_names))

# cost of moving to a state
# same for both treatments
#only cost associated with death - no cost in this model
trans_c_matrix <-
  matrix(c(0, 0, 0, cDeath,
           0, 0, 0, cDeath,
           0, 0, 0, cDeath,
           0, 0, 0, 0),
         byrow = TRUE,
         nrow = n_states,
         dimnames = list(from = s_names,
                         to = s_names))

# Transition probabilities ---- 

# time-homogeneous
#not by row - by column
p_matrix <- array(data = c(1-tpProgNDH - tpDeath, tpNorm, 0, 0,
                           tpProgNDH, 1-tpProgT2D-tpDeath-tpNorm, tpRegress, 0,
                           0, tpProgT2D,  1-tpRegress - (tpDeath*tpExcessDeath), 0,
                           tpDeath, tpDeath, (tpDeath*tpExcessDeath), 1,
                           #below would be adjusted for effect of DPP
                           1-tpProgNDH - tpDeath, tpNorm, 0, 0,
                           tpProgNDH, 1-(tpProgT2D*effect)-tpDeath-tpNorm, tpRegress, 0,
                           0, (tpProgT2D*effect),  1-tpRegress - (tpDeath*tpExcessDeath), 0,
                           tpDeath, tpDeath, (tpDeath*tpExcessDeath), 1),
                  dim = c(n_states, n_states, n_treatments),
                  dimnames = list(from = s_names,
                                  to = s_names,
                                  t_names))

# Store population output for each cycle 

# state populations
pop <- array(data = NA,
             dim = c(n_states, n_cycles, n_treatments),
             dimnames = list(state = s_names,
                             cycle = NULL,
                             treatment = t_names))

pop["Normal", cycle = 1, ] <-0
pop["NDH", cycle = 1, ] <- n_cohort
pop["T2D", cycle = 1, ] <- 0
pop["Dead", cycle = 1, ] <- 0

halfpop<-array(data=NA, dim=c(n_states,n_cycles,n_treatments), 
               dimnames = list(state= s_names, cycle=NULL, treatment= t_names))

halfpop["Normal", cycle=1, treatment="without DPP"]<- 0
halfpop["NDH", cycle=1, treatment="without DPP"]<- 0
halfpop["T2D", cycle=1, treatment="without DPP"]<- 0
halfpop["Dead", cycle=1, treatment="without DPP"]<- 0
halfpop["Normal", cycle=1, treatment="with DPP"]<- 0
halfpop["NDH", cycle=1, treatment="with DPP"]<- 0
halfpop["T2D", cycle=1, treatment="with DPP"]<- 0
halfpop["Dead", cycle=1, treatment="with DPP"]<- 0

# _arrived_ state populations
trans <- array(data = NA,
               dim = c(n_states, n_cycles, n_treatments),
               dimnames = list(state = s_names,
                               cycle = NULL,
                               treatment = t_names))

trans[, cycle = 1, ] <- 0


# Sum costs and QALYs for each cycle at a time for each treatment arm

cycle_empty_array <-
  array(NA,
        dim = c(n_treatments, n_cycles),
        dimnames = list(treatment = t_names,
                        cycle = NULL))

cycle_state_costs <- cycle_trans_costs <- cycle_qalys <- cycle_empty_array 
cycle_costs <- cycle_QALYs <- cycle_empty_array
LE <- LYs <- cycle_empty_array    # life expectancy; life-years
cycle_QALE <- cycle_empty_array   # quality-adjusted life expectancy

total_costs <- setNames(c(NA, NA), t_names)
total_QALYs <- setNames(c(NA, NA), t_names)

# #Transition probabilities
#normal to NDH
tpProgNDH <- 0.074202731
# #dying from NGT or NDH
tpDeath <- 0.05
# #Progression from NDH to T2D
tpProgT2D <- 0.0249
# #From NDH to NGT
tpNorm <- 0.0795
# #From T2D to NDH
tpRegress <- 0.002796084
# #Additional risk of death for T2D
tpExcessDeath <- 1.6
#Effectiveness of the DPP
effect <- 0.8

tpNorm <- rbeta(1,335,3709)
tpProgT2D <- rbeta(1,25.27,974.73)
tpRegress <- rbeta(1,2.8,997.2)
tpProgNDH <- rbeta(1,77.1,922.9)
tpExcessDeath <- rlnorm(1, meanlog = 0.470004, sdlog = 0.064111)

# # Time-dependent probability matrix ----
# Time dependent transition - risk of death
p_matrix_cycle <- function(p_matrix, age, cycle,
                           tpProgNDH = rbeta(1,77.1,922.9),
                           #tpProgT2D = rbeta(1,25.27,974.73),
                           tpNorm = rbeta(1,335,3709),
                           tpRegress = rbeta(1,2.8,997.2),
                           tpExcessDeath = rlnorm(1, meanlog = 0.470004, sdlog = 0.064111)
) {
  tpDeath_lookup <-
    c("[30,34]" = 0.000674,
      "[35,39]" = 0.001007,
      "[40,44]" = 0.001499,
      "[45,49]" = 0.002298,
      "[50,54]" = 0.003359,
      "[55,59]" = 0.005056,
      "[60,64]" = 0.007940,
      "[65,69]" = 0.012389,
      "[70,74]" = 0.019640,
      "[75,79]" = 0.034483,
      "[80,84]" = 0.060938,
      "[85,89]" = 0.110627,
      "[90,94]" = 0.187364,
      "[95,120]" = 0.295984)
  
  age_grp <- cut(age, breaks = c(30,34,39,44,49,54,59,64,69,74,79, 84, 89, 94, 120))
  tpDeath <- tpDeath_lookup[age_grp]
  
  
  tpProgT2D_lookup <-
    c("[30,39]" = rnorm(1,mean=0.108, sd=0.0095),
      "[40,49]" = rnorm(1,mean=0.108, sd=0.0095),
      "[50,59]" = rnorm(1,mean=0.108, sd=0.0095),
      "[60,69]" = rnorm(1,mean=0.108, sd=0.0095),
      "[70,79]" = rnorm(1,mean=0.108, sd=0.0095),
      "[80,120]" = rnorm(1,mean=0.108, sd=0.0095))
  
  age_grpT2D <- cut(age, breaks= c(30,39,49,59,69,79,120))
  tpProgT2D <-tpProgT2D_lookup[age_grpT2D]
  
  # [ includes, ( doesn't include
  effect_lookup <-
    c("(0,3]" = rlnorm(1, meanlog = -0.22314, sdlog = 0.044757),
      "(3,10]" = rnorm(1, mean = 1.000, sd = 0.00000001),
      "(10,15]"= rnorm(1, mean = 1.000, sd = 0.00000001),
      "(15,50]" = rnorm(1, mean = 1.000, sd = 0.00000001))
  
  cycle2 <- cut(cycle, breaks = c(0,3,10,15,50))
  
  effect <- effect_lookup[cycle2]
  
  #########
  
  p_matrix["Normal", "Normal", "without DPP"] <- 1 - tpProgNDH - tpDeath
  p_matrix["Normal", "NDH", "without DPP"] <- tpProgNDH 
  p_matrix["Normal", "Dead", "without DPP"] <- tpDeath
  p_matrix["NDH", "Normal", "without DPP"] <- tpNorm
  p_matrix["NDH", "NDH", "without DPP"] <-1-tpProgT2D-tpDeath-tpNorm
  p_matrix["NDH", "T2D", "without DPP"] <- tpProgT2D
  p_matrix["NDH", "Dead", "without DPP"] <- tpDeath
  p_matrix["T2D", "NDH", "without DPP"] <- tpRegress
  p_matrix["T2D", "T2D", "without DPP"] <- 1-tpRegress - (tpDeath*tpExcessDeath)
  p_matrix["T2D", "Dead", "without DPP"] <- (tpDeath)*tpExcessDeath
  p_matrix["Dead", "Dead", "without DPP"] <- 1
  
  p_matrix["Normal", "Normal", "with DPP"] <- 1 - tpProgNDH - tpDeath
  p_matrix["Normal", "NDH", "with DPP"] <- tpProgNDH
  p_matrix["Normal", "Dead", "with DPP"] <- tpDeath
  p_matrix["NDH", "Normal", "with DPP"] <- tpNorm
  p_matrix["NDH", "NDH", "with DPP"] <- 1 - (tpProgT2D*effect) - tpDeath - tpNorm
  p_matrix["NDH", "T2D", "with DPP"] <- (tpProgT2D*effect)
  p_matrix["NDH", "Dead", "with DPP"] <- tpDeath
  p_matrix["T2D", "NDH", "with DPP"] <- tpRegress
  p_matrix["T2D", "T2D", "with DPP"] <- 1 - tpRegress - (tpDeath*tpExcessDeath)
  p_matrix["T2D", "Dead", "with DPP"] <- tpDeath*tpExcessDeath
  p_matrix["Dead", "Dead", "with DPP"] <- 1
  
  return(p_matrix)
}


###########
#Making utility change with age
#state_q_matrix
state_q_matrix_cycle <- function(state_q_matrix, age)
{
  #uNorm 
  uNorm_lookup <-
    c("[30,39]" = rbeta(1,5243.210127, 647.8882745),
      "[40,49]" = rbeta(1,3749.390288, 552.6603545),
      "[50,59]" = rbeta(1,4317.392142, 733.4582261),
      "[60,69]" = rbeta(1,3275.399881, 605.1330164),
      "[70,79]" = rbeta(1,2793.5369, 592.8093771),
      "[80,120]" = rbeta(1,1048.602386, 338.6611273))
  
  
  age_grpnorm <- cut(age, breaks= c(30,39,49,59,69,79,120))
  uNorm <-uNorm_lookup[age_grpnorm]
  
  #uNDH
  uNDH_lookup <-
    c("[30,39]" = rbeta(1,18852.2607, 3585.449629),
      "[40,49]" = rbeta(1,40397.03381, 8321.148331),
      "[50,59]" = rbeta(1,79459.6746, 19951.65504),
      "[60,69]" = rbeta(1,137085.1892, 32219.07565),
      "[70,79]" = rbeta(1,164314.0338, 37615.0213),
      "[80,120]" = rbeta(1,53445.85673, 14946.63982))
  
  
  age_grpNDH <- cut(age, breaks= c(30,39,49,59,69,79,120))
  uNDH <-uNDH_lookup[age_grpNDH]
  
  #uT2D
  uT2D_lookup <-
    c("[30,39]" = rbeta(1,27.8648312, 8.05807557),
      "[40,49]" = rbeta(1,48.96668398, 23.6572448),
      "[50,59]" = rbeta(1,158.9763077, 60.70760667),
      "[60,69]" = rbeta(1,238.2666204, 102.8036175),
      "[70,79]" = rbeta(1,294.0078638, 100.8414828),
      "[80,120]" = rbeta(1,198.9254422, 97.37988722))
  
  
  age_grpT2D <- cut(age, breaks= c(30,39,49,59,69,79,120))
  uT2D <-uT2D_lookup[age_grpT2D]
  
  state_q_matrix[1,1]<- uNorm
  state_q_matrix["without DPP", "NDH"] <- uNDH
  state_q_matrix["without DPP", "T2D"] <- uT2D
  state_q_matrix["without DPP", "Dead"] <-0
  state_q_matrix["with DPP", "Normal"] <- uNorm
  state_q_matrix["with DPP", "NDH"] <- uNDH
  state_q_matrix["with DPP", "T2D"] <-uT2D
  state_q_matrix["with DPP", "Dead"]<-0
  
  
  return(state_q_matrix)
}




## Run model ----
for(n in 1:n_trials) {
  
  for (i in 1:n_treatments) {
    
    age <- Initial_age
    
    for (j in 2:n_cycles) {
      
      p_matrix <- p_matrix_cycle(p_matrix, age, j - 1)
      
      pop[, cycle = j, treatment = i] <-
        pop[, cycle = j - 1, treatment = i] %*% p_matrix[, , treatment = i]
      
      #Transition costs - none included in this model
      trans[, cycle = j, treatment = i] <-
        pop[, cycle = j - 1, treatment = i] %*% (trans_c_matrix * p_matrix[, , treatment = i])
      
      halfpop["Normal", cycle=j, treatment=i]<- 0.5*(pop["Normal", cycle=(j-1), treatment=i]+pop["Normal", cycle=j, treatment=i])
      halfpop["NDH", cycle=j, treatment=i]<- 0.5*(pop["NDH", cycle=(j-1), treatment=i]+pop["NDH", cycle=j, treatment=i])
      halfpop["T2D", cycle=j, treatment=i]<- 0.5*(pop["T2D", cycle=(j-1), treatment=i]+pop["T2D", cycle=j, treatment=i])
      halfpop["Dead", cycle=j, treatment=i]<- 0.5*(pop["Dead", cycle=(j-1), treatment=i]+pop["Dead", cycle=j, treatment=i])
      
      cycle_qalys[i, ] <-
        (state_q_matrix[treatment = i, ] %*% halfpop[, , treatment = i]) * (1/(1 + cDr))^(-1:(n_cycles-2))
      
      age <- age + 1
    }
    
    #Applying half cycle correction
    cycle_state_costs[i, ] <-
      (state_c_matrix[treatment = i, ] %*% halfpop[, , treatment = i]) * (1/(1 + cDr))^(-1:(n_cycles-2)) 
    
    cycle_costs[i, ] <- cycle_state_costs[i, ] #+ cycle_trans_costs[i, ]
    total_costs[i] <- sum(cycle_costs[treatment = i, -1])
    total_QALYs[i] <- sum(cycle_qalys[treatment = i, -1])
  }
  #adding in costs of DPP participation
  total_costs["with DPP"] <- total_costs["with DPP"] + dppcost_sum
  
  #adding in benefits from DPP participation
  #multipling by 0.5 as this is utilty gain - need to convert to QALY
  #agecat1: rbeta(1,35.52768105,0.153650952)
  #agecat2: rbeta(1,34.56207, 0.211612)
  #agecat3: rbeta(1,34.18377, 0.307204)
  #agecat4: rbeta(1,36.65039, 0.515016)
  #agecat5: rbeta(1,37.76731, 0.58579)
  #agecat6: rbeta(1, 34.39323, 0.371137)
  
  total_QALYs["with DPP"] <- total_QALYs["with DPP"] + 0.5*utility_sum
  
  
  costs[n, ] <- total_costs
  qalys[n,] <- total_QALYs
}


agecat5 <- data.frame(costs, qalys)
names(agecat5)[names(agecat5) == "without.DPP.1"] <- "Without.DPP.qalys.AGECAT5"
names(agecat5)[names(agecat5) == "with.DPP.1"] <- "With.DPP.qalys.AGECAT5"
names(agecat5)[names(agecat5) == "without.DPP"] <- "Without.DPP.costs.AGECAT5"
names(agecat5)[names(agecat5) == "with.DPP"] <- "With.DPP.costs.AGECAT5"
writexl::write_xlsx(agecat5, "age category5 - 244.xlsx")

save(agecat5,file="agecat5.Rda")

popagecat5 <- data.frame(pop)
writexl::write_xlsx(popagecat5, "age category5 trace.xlsx")
save(popagecat5,file="agecat5trace.Rda")


########################################
# AGE CAT 6
########################################
n_cohort <- 90
n_cycles <- 36
Initial_age <- 85


#NUMBER OF PSA ITERATIONS
n_trials<-10000
costs <- matrix(NA,nrow=n_trials,ncol=n_treatments, dimnames=list(NULL,t_names))
qalys <- matrix(NA,nrow=n_trials,ncol=n_treatments, dimnames=list(NULL,t_names))


#Utility of states
# uNorm <- 0.8395923 / uNDH <- 0.8085401 / uT2D <- 0.708957 / uDeath <-0 / uDPP<-0 


#Cost of states
cNorm <- rgamma(1, 44.44444, scale=45.119475)
cNDH <- rgamma(1,  44.44444, scale=50.0571)
cT2D <- rgamma(1,  44.44444, scale=99.4628)
cDeath <- 0
cDPP <- 0
#Utility of states
uNorm <- rbeta(1,27.53076303,4.883160879)
uNDH <- rbeta(1,26.33007731, 6.334144666)
uT2D <- rbeta(1, 26.52841, 10.65937)
uDeath <-0
uDPP<-0 
#Discount rates
oDr <- 0.035
cDr <- 0.035


#Transition probabilities
#normal to NDH
tpProgNDH <- 0.074202731
#dying from NGT or NDH
tpDeath <- 0.05
#Progression to T2D from NDH
tpProgT2D <- 0.0249
#From NDH to NGT
tpNorm <- 0.0795
#From T2D to NDH
tpRegress <- 0.002796084
#Additional risk of death for T2D
tpExcessDeath <- 1.6

#Effectiveness of the DPP
effect <- 0.8

#utility gains from DPP
utility <- c(1-rbeta(n_cohort,1250791.599, 13497.29538))
utility_sum <-sum(utility)

#costs of the DPP
dppcost <- c(rgamma(n_cohort,16.73964,scale=8.469118))
dppcost_sum <- sum(dppcost)

# cost of staying in state

state_c_matrix <-
  matrix(c(cNorm, cNDH, cT2D, 0,
           cNorm, cNDH + cDPP, cT2D, 0),
         byrow = TRUE,
         nrow = n_treatments,
         dimnames = list(t_names,
                         s_names))

# qaly when staying in state
state_q_matrix <-
  matrix(c(uNorm, uNDH, uT2D, 0,
           uNorm, uNDH + uDPP, uT2D, 0),
         byrow = TRUE,
         nrow = n_treatments,
         dimnames = list(t_names,
                         s_names))

# cost of moving to a state
# same for both treatments
#only cost associated with death - no cost in this model
trans_c_matrix <-
  matrix(c(0, 0, 0, cDeath,
           0, 0, 0, cDeath,
           0, 0, 0, cDeath,
           0, 0, 0, 0),
         byrow = TRUE,
         nrow = n_states,
         dimnames = list(from = s_names,
                         to = s_names))

# Transition probabilities ---- 

# time-homogeneous
#not by row - by column
p_matrix <- array(data = c(1-tpProgNDH - tpDeath, tpNorm, 0, 0,
                           tpProgNDH, 1-tpProgT2D-tpDeath-tpNorm, tpRegress, 0,
                           0, tpProgT2D,  1-tpRegress - (tpDeath*tpExcessDeath), 0,
                           tpDeath, tpDeath, (tpDeath*tpExcessDeath), 1,
                           #below would be adjusted for effect of DPP
                           1-tpProgNDH - tpDeath, tpNorm, 0, 0,
                           tpProgNDH, 1-(tpProgT2D*effect)-tpDeath-tpNorm, tpRegress, 0,
                           0, (tpProgT2D*effect),  1-tpRegress - (tpDeath*tpExcessDeath), 0,
                           tpDeath, tpDeath, (tpDeath*tpExcessDeath), 1),
                  dim = c(n_states, n_states, n_treatments),
                  dimnames = list(from = s_names,
                                  to = s_names,
                                  t_names))

# Store population output for each cycle 

# state populations
pop <- array(data = NA,
             dim = c(n_states, n_cycles, n_treatments),
             dimnames = list(state = s_names,
                             cycle = NULL,
                             treatment = t_names))

pop["Normal", cycle = 1, ] <-0
pop["NDH", cycle = 1, ] <- n_cohort
pop["T2D", cycle = 1, ] <- 0
pop["Dead", cycle = 1, ] <- 0

halfpop<-array(data=NA, dim=c(n_states,n_cycles,n_treatments), 
               dimnames = list(state= s_names, cycle=NULL, treatment= t_names))

halfpop["Normal", cycle=1, treatment="without DPP"]<- 0
halfpop["NDH", cycle=1, treatment="without DPP"]<- 0
halfpop["T2D", cycle=1, treatment="without DPP"]<- 0
halfpop["Dead", cycle=1, treatment="without DPP"]<- 0
halfpop["Normal", cycle=1, treatment="with DPP"]<- 0
halfpop["NDH", cycle=1, treatment="with DPP"]<- 0
halfpop["T2D", cycle=1, treatment="with DPP"]<- 0
halfpop["Dead", cycle=1, treatment="with DPP"]<- 0

# _arrived_ state populations
trans <- array(data = NA,
               dim = c(n_states, n_cycles, n_treatments),
               dimnames = list(state = s_names,
                               cycle = NULL,
                               treatment = t_names))

trans[, cycle = 1, ] <- 0


# Sum costs and QALYs for each cycle at a time for each treatment arm

cycle_empty_array <-
  array(NA,
        dim = c(n_treatments, n_cycles),
        dimnames = list(treatment = t_names,
                        cycle = NULL))

cycle_state_costs <- cycle_trans_costs <- cycle_qalys <- cycle_empty_array 
cycle_costs <- cycle_QALYs <- cycle_empty_array
LE <- LYs <- cycle_empty_array    # life expectancy; life-years
cycle_QALE <- cycle_empty_array   # quality-adjusted life expectancy

total_costs <- setNames(c(NA, NA), t_names)
total_QALYs <- setNames(c(NA, NA), t_names)

# #Transition probabilities
#normal to NDH
tpProgNDH <- 0.074202731
# #dying from NGT or NDH
tpDeath <- 0.05
# #Progression from NDH to T2D
tpProgT2D <- 0.0249
# #From NDH to NGT
tpNorm <- 0.0795
# #From T2D to NDH
tpRegress <- 0.002796084
# #Additional risk of death for T2D
tpExcessDeath <- 1.6
#Effectiveness of the DPP
effect <- 0.8

tpNorm <- rbeta(1,335,3709)
tpProgT2D <- rbeta(1,25.27,974.73)
tpRegress <- rbeta(1,2.8,997.2)
tpProgNDH <- rbeta(1,77.1,922.9)
tpExcessDeath <- rlnorm(1, meanlog = 0.470004, sdlog = 0.064111)

# # Time-dependent probability matrix ----
# Time dependent transition - risk of death
p_matrix_cycle <- function(p_matrix, age, cycle,
                           tpProgNDH = rbeta(1,77.1,922.9),
                           #tpProgT2D = rbeta(1,25.27,974.73),
                           tpNorm = rbeta(1,335,3709),
                           tpRegress = rbeta(1,2.8,997.2),
                           tpExcessDeath = rlnorm(1, meanlog = 0.470004, sdlog = 0.064111)
) {
  tpDeath_lookup <-
    c("[30,34]" = 0.000674,
      "[35,39]" = 0.001007,
      "[40,44]" = 0.001499,
      "[45,49]" = 0.002298,
      "[50,54]" = 0.003359,
      "[55,59]" = 0.005056,
      "[60,64]" = 0.007940,
      "[65,69]" = 0.012389,
      "[70,74]" = 0.019640,
      "[75,79]" = 0.034483,
      "[80,84]" = 0.060938,
      "[85,89]" = 0.110627,
      "[90,94]" = 0.187364,
      "[95,120]" = 0.295984)
  
  age_grp <- cut(age, breaks = c(30,34,39,44,49,54,59,64,69,74,79, 84, 89, 94, 120))
  tpDeath <- tpDeath_lookup[age_grp]
  
  tpProgT2D_lookup <-
    c("[30,39]" = rnorm(1,mean=0.108, sd=0.0095),
      "[40,49]" = rnorm(1,mean=0.108, sd=0.0095),
      "[50,59]" = rnorm(1,mean=0.108, sd=0.0095),
      "[60,69]" = rnorm(1,mean=0.108, sd=0.0095),
      "[70,79]" = rnorm(1,mean=0.108, sd=0.0095),
      "[80,120]" = rnorm(1,mean=0.108, sd=0.0095))
  
  age_grpT2D <- cut(age, breaks= c(30,39,49,59,69,79,120))
  tpProgT2D <-tpProgT2D_lookup[age_grpT2D]
  
  # [ includes, ( doesn't include
  effect_lookup <-
    c("(0,3]" = rlnorm(1, meanlog = -0.22314, sdlog = 0.044757),
      "(3,10]" = rnorm(1, mean = 1.000, sd = 0.00000001),
      "(10,15]"= rnorm(1, mean = 1.000, sd = 0.00000001),
      "(15,50]" = rnorm(1, mean = 1.000, sd = 0.00000001))
  
  cycle2 <- cut(cycle, breaks = c(0,3,10,15,50))
  
  effect <- effect_lookup[cycle2]
  
  #########
  
  p_matrix["Normal", "Normal", "without DPP"] <- 1 - tpProgNDH - tpDeath
  p_matrix["Normal", "NDH", "without DPP"] <- tpProgNDH 
  p_matrix["Normal", "Dead", "without DPP"] <- tpDeath
  p_matrix["NDH", "Normal", "without DPP"] <- tpNorm
  p_matrix["NDH", "NDH", "without DPP"] <-1-tpProgT2D-tpDeath-tpNorm
  p_matrix["NDH", "T2D", "without DPP"] <- tpProgT2D
  p_matrix["NDH", "Dead", "without DPP"] <- tpDeath
  p_matrix["T2D", "NDH", "without DPP"] <- tpRegress
  p_matrix["T2D", "T2D", "without DPP"] <- 1-tpRegress - (tpDeath*tpExcessDeath)
  p_matrix["T2D", "Dead", "without DPP"] <- (tpDeath)*tpExcessDeath
  p_matrix["Dead", "Dead", "without DPP"] <- 1
  
  p_matrix["Normal", "Normal", "with DPP"] <- 1 - tpProgNDH - tpDeath
  p_matrix["Normal", "NDH", "with DPP"] <- tpProgNDH
  p_matrix["Normal", "Dead", "with DPP"] <- tpDeath
  p_matrix["NDH", "Normal", "with DPP"] <- tpNorm
  p_matrix["NDH", "NDH", "with DPP"] <- 1 - (tpProgT2D*effect) - tpDeath - tpNorm
  p_matrix["NDH", "T2D", "with DPP"] <- (tpProgT2D*effect)
  p_matrix["NDH", "Dead", "with DPP"] <- tpDeath
  p_matrix["T2D", "NDH", "with DPP"] <- tpRegress
  p_matrix["T2D", "T2D", "with DPP"] <- 1 - tpRegress - (tpDeath*tpExcessDeath)
  p_matrix["T2D", "Dead", "with DPP"] <- tpDeath*tpExcessDeath
  p_matrix["Dead", "Dead", "with DPP"] <- 1
  
  return(p_matrix)
}


###########
#Making utility change with age
#state_q_matrix
state_q_matrix_cycle <- function(state_q_matrix, age)
{
  #uNorm 
  uNorm_lookup <-
    c("[30,39]" = rbeta(1,5243.210127, 647.8882745),
      "[40,49]" = rbeta(1,3749.390288, 552.6603545),
      "[50,59]" = rbeta(1,4317.392142, 733.4582261),
      "[60,69]" = rbeta(1,3275.399881, 605.1330164),
      "[70,79]" = rbeta(1,2793.5369, 592.8093771),
      "[80,120]" = rbeta(1,1048.602386, 338.6611273))
  
  
  age_grpnorm <- cut(age, breaks= c(30,39,49,59,69,79,120))
  uNorm <-uNorm_lookup[age_grpnorm]
  
  #uNDH
  uNDH_lookup <-
    c("[30,39]" = rbeta(1,18852.2607, 3585.449629),
      "[40,49]" = rbeta(1,40397.03381, 8321.148331),
      "[50,59]" = rbeta(1,79459.6746, 19951.65504),
      "[60,69]" = rbeta(1,137085.1892, 32219.07565),
      "[70,79]" = rbeta(1,164314.0338, 37615.0213),
      "[80,120]" = rbeta(1,53445.85673, 14946.63982))
  
  
  age_grpNDH <- cut(age, breaks= c(30,39,49,59,69,79,120))
  uNDH <-uNDH_lookup[age_grpNDH]
  
  #uT2D
  uT2D_lookup <-
    c("[30,39]" = rbeta(1,27.8648312, 8.05807557),
      "[40,49]" = rbeta(1,48.96668398, 23.6572448),
      "[50,59]" = rbeta(1,158.9763077, 60.70760667),
      "[60,69]" = rbeta(1,238.2666204, 102.8036175),
      "[70,79]" = rbeta(1,294.0078638, 100.8414828),
      "[80,120]" = rbeta(1,198.9254422, 97.37988722))
  
  
  age_grpT2D <- cut(age, breaks= c(30,39,49,59,69,79,120))
  uT2D <-uT2D_lookup[age_grpT2D]
  
  state_q_matrix[1,1]<- uNorm
  state_q_matrix["without DPP", "NDH"] <- uNDH
  state_q_matrix["without DPP", "T2D"] <- uT2D
  state_q_matrix["without DPP", "Dead"] <-0
  state_q_matrix["with DPP", "Normal"] <- uNorm
  state_q_matrix["with DPP", "NDH"] <- uNDH
  state_q_matrix["with DPP", "T2D"] <-uT2D
  state_q_matrix["with DPP", "Dead"]<-0
  
  
  return(state_q_matrix)
}




## Run model ----
for(n in 1:n_trials) {
  
  for (i in 1:n_treatments) {
    
    age <- Initial_age
    
    for (j in 2:n_cycles) {
      
      p_matrix <- p_matrix_cycle(p_matrix, age, j - 1)
      
      pop[, cycle = j, treatment = i] <-
        pop[, cycle = j - 1, treatment = i] %*% p_matrix[, , treatment = i]
      
      #Transition costs - none included in this model
      trans[, cycle = j, treatment = i] <-
        pop[, cycle = j - 1, treatment = i] %*% (trans_c_matrix * p_matrix[, , treatment = i])
      
      halfpop["Normal", cycle=j, treatment=i]<- 0.5*(pop["Normal", cycle=(j-1), treatment=i]+pop["Normal", cycle=j, treatment=i])
      halfpop["NDH", cycle=j, treatment=i]<- 0.5*(pop["NDH", cycle=(j-1), treatment=i]+pop["NDH", cycle=j, treatment=i])
      halfpop["T2D", cycle=j, treatment=i]<- 0.5*(pop["T2D", cycle=(j-1), treatment=i]+pop["T2D", cycle=j, treatment=i])
      halfpop["Dead", cycle=j, treatment=i]<- 0.5*(pop["Dead", cycle=(j-1), treatment=i]+pop["Dead", cycle=j, treatment=i])
      
      cycle_qalys[i, ] <-
        (state_q_matrix[treatment = i, ] %*% halfpop[, , treatment = i]) * (1/(1 + cDr))^(-1:(n_cycles-2))
      
      age <- age + 1
    }
    
    #Applying half cycle correction
    cycle_state_costs[i, ] <-
      (state_c_matrix[treatment = i, ] %*% halfpop[, , treatment = i]) * (1/(1 + cDr))^(-1:(n_cycles-2)) 
    
    cycle_costs[i, ] <- cycle_state_costs[i, ] #+ cycle_trans_costs[i, ]
    
    total_costs[i] <- sum(cycle_costs[treatment = i, -1])
    total_QALYs[i] <- sum(cycle_qalys[treatment = i, -1])
  }
  #adding in costs of DPP participation
  total_costs["with DPP"] <- total_costs["with DPP"] + dppcost_sum
  
  #adding in benefits from DPP participation
  #multipling by 0.5 as this is utilty gain - need to convert to QALY
  #agecat1: rbeta(1,35.52768105,0.153650952)
  #agecat2: rbeta(1,34.56207, 0.211612)
  #agecat3: rbeta(1,34.18377, 0.307204)
  #agecat4: rbeta(1,36.65039, 0.515016)
  #agecat5: rbeta(1,37.76731, 0.58579)
  #agecat6: rbeta(1, 34.39323, 0.371137)
  
  total_QALYs["with DPP"] <- total_QALYs["with DPP"] + 0.5*utility_sum
  
  
  costs[n, ] <- total_costs
  qalys[n,] <- total_QALYs
}


agecat6 <- data.frame(costs, qalys)
names(agecat6)[names(agecat6) == "without.DPP.1"] <- "Without.DPP.qalys.AGECAT6"
names(agecat6)[names(agecat6) == "with.DPP.1"] <- "With.DPP.qalys.AGECAT6"
names(agecat6)[names(agecat6) == "without.DPP"] <- "Without.DPP.costs.AGECAT6"
names(agecat6)[names(agecat6) == "with.DPP"] <- "With.DPP.costs.AGECAT6"
writexl::write_xlsx(agecat6, "age category6 - 90.xlsx")
save(agecat6,file="agecat6.Rda")

popagecat6 <- data.frame(pop)
writexl::write_xlsx(popagecat6, "age category6 trace.xlsx")
save(popagecat6,file="agecat6trace.Rda")

############################
# Model Analysis #
############################

set.seed(9009)


load("agecat6.Rda")
load("agecat5.Rda")
load("agecat4.Rda")
load("agecat3.Rda")
load("agecat2.Rda")
load("agecat1.Rda")



#willingness to pay threshold
wtp <- 20000
#creating a common variable ID so they can be merged 'trials'

n_trials<-10000
trials<-c(seq(1, n_trials, by=1))
agecat1$trials <- trials
agecat2$trials <- trials
agecat3$trials <- trials
agecat4$trials <- trials
agecat5$trials <- trials
agecat6$trials <- trials


#merging the two data frames
age12 <- merge(agecat1, agecat2, by="trials")
age123 <- merge(age12, agecat3, by="trials")
age1234 <- merge(age123, agecat4, by="trials")
age12345 <- merge(age1234, agecat5, by="trials")
whole_cohort <- merge(age12345, agecat6, by="trials")

writexl::write_xlsx(whole_cohort, "whole cohort.xlsx")
save(whole_cohort,file="whole_cohort.Rda")

#create a variable that's the sum of all the rows
#create datasets for QALYs & COSTS according to DPP / without DPP
#sum acrross the rows - using rowSums
#create a dataset from this

#With DPP - costs
DPP_costs =data.frame(whole_cohort$With.DPP.costs.AGECAT6,whole_cohort$With.DPP.costs.AGECAT5,whole_cohort$With.DPP.costs.AGECAT4,whole_cohort$With.DPP.costs.AGECAT3,whole_cohort$With.DPP.costs.AGECAT2,whole_cohort$With.DPP.costs.AGECAT1)
DPP_costs$total_costs=rowSums(DPP_costs)
#without DPP - costs
withoutDPP_costs =data.frame(whole_cohort$Without.DPP.costs.AGECAT6,whole_cohort$Without.DPP.costs.AGECAT5,whole_cohort$Without.DPP.costs.AGECAT4,whole_cohort$Without.DPP.costs.AGECAT3,whole_cohort$Without.DPP.costs.AGECAT2,whole_cohort$Without.DPP.costs.AGECAT1)
withoutDPP_costs$total_costs=rowSums(withoutDPP_costs)

#With DPP - qalys
DPP_qalys =data.frame(whole_cohort$With.DPP.qalys.AGECAT6,whole_cohort$With.DPP.qalys.AGECAT5,whole_cohort$With.DPP.qalys.AGECAT4,whole_cohort$With.DPP.qalys.AGECAT3,whole_cohort$With.DPP.qalys.AGECAT2,whole_cohort$With.DPP.qalys.AGECAT1)
DPP_qalys$total_qalys=rowSums(DPP_qalys)
#without DPP - qalys
withoutDPP_qalys =data.frame(whole_cohort$Without.DPP.qalys.AGECAT6,whole_cohort$Without.DPP.qalys.AGECAT5,whole_cohort$Without.DPP.qalys.AGECAT4,whole_cohort$Without.DPP.qalys.AGECAT3,whole_cohort$Without.DPP.qalys.AGECAT2,whole_cohort$Without.DPP.qalys.AGECAT1)
withoutDPP_qalys$total_qalys=rowSums(withoutDPP_qalys)

PSA_results =data.frame(DPP_costs$total_costs, withoutDPP_costs$total_costs, DPP_qalys$total_qalys, withoutDPP_qalys$total_qalys)
writexl::write_xlsx(PSA_results, "PSA results from mixed age cohort.xlsx")
save(PSA_results,file="PSA_results.Rda")


#creating a variable that calculates the difference between costs / qalys
PSA_results$cost_dif <- PSA_results$DPP_costs.total_costs - PSA_results$withoutDPP_costs.total_costs 
PSA_results$qaly_dif <- PSA_results$DPP_qalys.total_qalys - PSA_results$withoutDPP_qalys.total_qalys 
#calculating the incremental net monetary benefit at different WTP thresholds
PSA_results$inmb0 <- ((PSA_results$qaly_dif * 0) - PSA_results$cost_dif)/1000
PSA_results$inmb1 <- ((PSA_results$qaly_dif * 1000) - PSA_results$cost_dif)/1000
PSA_results$inmb2 <- ((PSA_results$qaly_dif * 2000) - PSA_results$cost_dif)/1000
PSA_results$inmb3 <- ((PSA_results$qaly_dif * 3000) - PSA_results$cost_dif)/1000
PSA_results$inmb4 <- ((PSA_results$qaly_dif * 4000) - PSA_results$cost_dif)/1000
PSA_results$inmb5 <- ((PSA_results$qaly_dif * 5000) - PSA_results$cost_dif)/1000
PSA_results$inmb6 <- ((PSA_results$qaly_dif * 6000) - PSA_results$cost_dif)/1000
PSA_results$inmb7 <- ((PSA_results$qaly_dif * 7000) - PSA_results$cost_dif)/1000
PSA_results$inmb8 <- ((PSA_results$qaly_dif * 8000) - PSA_results$cost_dif)/1000
PSA_results$inmb9 <- ((PSA_results$qaly_dif * 9000) - PSA_results$cost_dif)/1000
PSA_results$inmb10 <- ((PSA_results$qaly_dif * 10000) - PSA_results$cost_dif)/1000
PSA_results$inmb11 <- ((PSA_results$qaly_dif * 11000) - PSA_results$cost_dif)/1000
PSA_results$inmb12 <- ((PSA_results$qaly_dif * 12000) - PSA_results$cost_dif)/1000
PSA_results$inmb13 <- ((PSA_results$qaly_dif * 13000) - PSA_results$cost_dif)/1000
PSA_results$inmb14 <- ((PSA_results$qaly_dif * 14000) - PSA_results$cost_dif)/1000
PSA_results$inmb15 <- ((PSA_results$qaly_dif * 15000) - PSA_results$cost_dif)/1000
PSA_results$inmb16 <- ((PSA_results$qaly_dif * 16000) - PSA_results$cost_dif)/1000
PSA_results$inmb17 <- ((PSA_results$qaly_dif * 17000) - PSA_results$cost_dif)/1000
PSA_results$inmb18 <- ((PSA_results$qaly_dif * 18000) - PSA_results$cost_dif)/1000
PSA_results$inmb19 <- ((PSA_results$qaly_dif * 19000) - PSA_results$cost_dif)/1000
PSA_results$inmb20 <- ((PSA_results$qaly_dif * 20000) - PSA_results$cost_dif)/1000
PSA_results$inmb21 <- ((PSA_results$qaly_dif * 21000) - PSA_results$cost_dif)/1000
PSA_results$inmb22 <- ((PSA_results$qaly_dif * 22000) - PSA_results$cost_dif)/1000
PSA_results$inmb23 <- ((PSA_results$qaly_dif * 23000) - PSA_results$cost_dif)/1000
PSA_results$inmb24 <- ((PSA_results$qaly_dif * 24000) - PSA_results$cost_dif)/1000
PSA_results$inmb25 <- ((PSA_results$qaly_dif * 25000) - PSA_results$cost_dif)/1000
PSA_results$inmb26 <- ((PSA_results$qaly_dif * 26000) - PSA_results$cost_dif)/1000
PSA_results$inmb27 <- ((PSA_results$qaly_dif * 27000) - PSA_results$cost_dif)/1000
PSA_results$inmb28 <- ((PSA_results$qaly_dif * 28000) - PSA_results$cost_dif)/1000
PSA_results$inmb29 <- ((PSA_results$qaly_dif * 29000) - PSA_results$cost_dif)/1000
PSA_results$inmb30 <- ((PSA_results$qaly_dif * 30000) - PSA_results$cost_dif)/1000
PSA_results$inmb31 <- ((PSA_results$qaly_dif * 31000) - PSA_results$cost_dif)/1000
PSA_results$inmb32 <- ((PSA_results$qaly_dif * 32000) - PSA_results$cost_dif)/1000
PSA_results$inmb33 <- ((PSA_results$qaly_dif * 33000) - PSA_results$cost_dif)/1000
PSA_results$inmb34 <- ((PSA_results$qaly_dif * 34000) - PSA_results$cost_dif)/1000
PSA_results$inmb35 <- ((PSA_results$qaly_dif * 35000) - PSA_results$cost_dif)/1000
PSA_results$inmb36 <- ((PSA_results$qaly_dif * 36000) - PSA_results$cost_dif)/1000
PSA_results$inmb37 <- ((PSA_results$qaly_dif * 37000) - PSA_results$cost_dif)/1000
PSA_results$inmb38 <- ((PSA_results$qaly_dif * 38000) - PSA_results$cost_dif)/1000
PSA_results$inmb39 <- ((PSA_results$qaly_dif * 39000) - PSA_results$cost_dif)/1000
PSA_results$inmb40 <- ((PSA_results$qaly_dif * 40000) - PSA_results$cost_dif)/1000
PSA_results$inmb41 <- ((PSA_results$qaly_dif * 41000) - PSA_results$cost_dif)/1000
PSA_results$inmb42 <- ((PSA_results$qaly_dif * 42000) - PSA_results$cost_dif)/1000
PSA_results$inmb43 <- ((PSA_results$qaly_dif * 43000) - PSA_results$cost_dif)/1000
PSA_results$inmb44 <- ((PSA_results$qaly_dif * 44000) - PSA_results$cost_dif)/1000
PSA_results$inmb45 <- ((PSA_results$qaly_dif * 45000) - PSA_results$cost_dif)/1000

#calculating number of trials that are cost-effective given willingness to pay
wtprange<-c(seq(0,45000,1000))
per_ce<-c(rep(0,length(wtprange)))
ceac<-data.frame(wtprange, per_ce)

ceac$per_ce[1] <- (sum(PSA_results$inmb0 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[2] <- (sum(PSA_results$inmb1 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[3] <- (sum(PSA_results$inmb2 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[4] <- (sum(PSA_results$inmb3 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[5] <- (sum(PSA_results$inmb4 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[6] <- (sum(PSA_results$inmb5 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[7] <- (sum(PSA_results$inmb6 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[8] <- (sum(PSA_results$inmb7 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[9] <- (sum(PSA_results$inmb8 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[10] <- (sum(PSA_results$inmb9 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[11] <- (sum(PSA_results$inmb10 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[12] <- (sum(PSA_results$inmb11 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[13] <- (sum(PSA_results$inmb12 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[14] <- (sum(PSA_results$inmb13 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[15] <- (sum(PSA_results$inmb14 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[16] <- (sum(PSA_results$inmb15 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[17] <- (sum(PSA_results$inmb16 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[18] <- (sum(PSA_results$inmb17 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[19] <- (sum(PSA_results$inmb18 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[20] <- (sum(PSA_results$inmb19 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[21] <- (sum(PSA_results$inmb20 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[22] <- (sum(PSA_results$inmb21 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[23] <- (sum(PSA_results$inmb22 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[24] <- (sum(PSA_results$inmb23 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[25] <- (sum(PSA_results$inmb24 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[26] <- (sum(PSA_results$inmb25 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[27] <- (sum(PSA_results$inmb26 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[28] <- (sum(PSA_results$inmb27 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[29] <- (sum(PSA_results$inmb28 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[30] <- (sum(PSA_results$inmb29 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[31] <- (sum(PSA_results$inmb30 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[32] <- (sum(PSA_results$inmb31 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[33] <- (sum(PSA_results$inmb32 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[34] <- (sum(PSA_results$inmb33 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[35] <- (sum(PSA_results$inmb34 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[36] <- (sum(PSA_results$inmb35 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[37] <- (sum(PSA_results$inmb36 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[38] <- (sum(PSA_results$inmb37 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[39] <- (sum(PSA_results$inmb38 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[40] <- (sum(PSA_results$inmb39 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[41] <- (sum(PSA_results$inmb40 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[42] <- (sum(PSA_results$inmb41 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[43] <- (sum(PSA_results$inmb42 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[44] <- (sum(PSA_results$inmb43 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[45] <- (sum(PSA_results$inmb44 >= 0, na.rm=TRUE)/n_trials)
ceac$per_ce[46] <- (sum(PSA_results$inmb45 >= 0, na.rm=TRUE)/n_trials)


plot(ceac$wtprange,ceac$per_ce, type="l", ylim=c(0,1), col="black", lwd=2, xlab="Willingness to Pay (£)", ylab="Probability Cost-effective", main="Cost-effectiveness Acceptability Curve")
writexl::write_xlsx(ceac, "ceac.xlsx")

percent_CE_20 <- (sum(PSA_results$inmb20 >= 0, na.rm=TRUE)/n_trials)*100
percent_CE_30 <- (sum(PSA_results$inmb30 >= 0, na.rm=TRUE)/n_trials)*100


NE <-sum(PSA_results$cost_dif>= 0 & PSA_results$qaly_dif>= 0, na.rm=TRUE)
NW <- sum(PSA_results$cost_dif>= 0 & PSA_results$qaly_dif<= 0, na.rm=TRUE)
SE <- sum(PSA_results$cost_dif<= 0 & PSA_results$qaly_dif>= 0, na.rm=TRUE)
SW <- sum(PSA_results$cost_dif<= 0 & PSA_results$qaly_dif<= 0, na.rm=TRUE)

PSA_results$qaly_difpp <- (PSA_results$qaly_dif)/1000
PSA_results$cost_difpp <- (PSA_results$cost_dif)/1000


#exporting the results
writexl::write_xlsx(PSA_results, "incremental_nmb.xlsx")

#plotting per person incremental
#find out range of x variable
range(PSA_results$qaly_difpp)
#0.01826279 0.17007453
ceacplot <- ggplot(data=PSA_results,aes(y=cost_difpp))+
  geom_point(aes(x=qaly_difpp),color="blue")+
  geom_vline(xintercept=0,color="grey", size=1)+
  geom_hline(yintercept=0,color="grey", size=1)+
  xlim(-0.01,0.2)+
  xlab("Incremental QALYs")+
  ylab("Incremental Costs")+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5))

ceacplot + geom_abline(intercept = 0, slope = 20000, color="black", linetype="solid", size=1)

#Average costs:
mean(DPP_costs$total_costs)
#39365800
mean(withoutDPP_costs$total_costs)
#39687183

#Average QALYs:
mean(DPP_qalys$total_qalys)
#9738.393
mean(withoutDPP_qalys$total_qalys)
#9648.213

#Average incremental net benefits
#20000 WTP
mean(PSA_results$inmb20)
#2124.97
#30,000 WTP
mean(PSA_results$inmb30)
#3026.764

#Number of montecarlo simulations that are cost-effective
ceac$per_ce[21]
#100%
ceac$per_ce[31]
#100.0%


