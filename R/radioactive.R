#File: radioactive.R
#Author: Juha Karvanen
#Date: 2020-01-27
#Summary: functions used in "Estimating mean lifetime from partially observed events in nuclear physics"

library(data.table)
library(Rcpp)

#################################################################################################################
#################################################################################################################
#################################################################################################################
# 1) Functions for estimation #
###############################



##########################################################################################################
# Estimator for complete data #
###############################
# arrivals = a vector of arrival times
# departures = a vector of departure times
#
# Returns the estimated mean lifetime
##########################################################################################################
complete_estim <- function(arrivals,departures) {
  return( (sum(departures) - sum(arrivals)) / length(arrivals) )   
}


##########################################################################################################
# Subroutine blocksum for finding the blocks #
##############################################
cpp1 <- cppFunction("
  NumericVector blocksum(NumericVector x) {
    NumericVector bsum(x.length()); 
    double tempsum = 0;
    for(int i = 0; i < x.length(); i++){
      tempsum += x[i];
      if(tempsum > 1) tempsum = 1;
      bsum[i] = tempsum;
      if(tempsum > 0) tempsum = 0;
    }
    return bsum;
  }")


is_valid_thinning <- function(x_diff) {
  return( (min(cumsum(x_diff)) >= 0) )   
}  

generate_compatible <- function(x_diff, N=1) {
  n <- sum( x_diff == 1 )
  m <- sum( x_diff == -1 )
  result <- matrix(NA,N,n+m)
  for(i in 1:N) {
    arrind <- (1:(n+m))[x_diff == 1]
    x_diff_new <- x_diff
    k <- n
    while( k > m ) {
      repeat {
        newindind <- sample.int(length(arrind),1)  
        newind <- arrind[newindind]
        arrind <- arrind[-newindind]
        x_diff_new [newind] <- 0
        if( min(cumsum(x_diff_new )) >= 0 ) break
        x_diff_new[newind] <- 1
      }
      k <- k - 1
    }  
    result[i,] <- x_diff_new 
  }
  return( result )
}
# generate_compatible(c(1,-1,1,1,-1,1,1,1,-1,-1))
# aa <- generate_compatible(c(1,-1,1,1,-1,1,1,1,-1,-1,1,1,-1),100000)
# aa <- generate_compatible(c(1,-1,1,1,-1,1,1,1,-1,-1), 100000)
# aatext <- apply(aa,1,paste,collapse="")
# table(aatext)
# 
# aa <- generate_compatible(c(1,1,-1,1,1,-1), 100000)
# aatext <- apply(aa,1,paste,collapse="")
# table(aatext)
# The probabilities should be equal but they are not:
# 11-100-1 : 1/4 * 1/3 + 1/4 * 1/3 = 1/6
# 01-101-1 : 1/4 * 1/3 + 1/4 * 1/2 = 5/24


generate_compatible2 <- function(x_diff, N=1) {
  n <- sum( x_diff == 1 )
  m <- sum( x_diff == -1 )
  result <- matrix(NA,N,n+m)
  weight <- rep(1,N)
  for(i in 1:N) {
    arrind <- (1:(n+m))[x_diff == 1]
    x_diff_new <- x_diff
    k <- n
    while( k > m ) {
      repeat {
        newindind <- sample.int(length(arrind),1)  
        newind <- arrind[newindind]
        arrind <- arrind[-newindind]
        x_diff_new [newind] <- 0
        if( min(cumsum(x_diff_new )) >= 0 ) break
        x_diff_new[newind] <- 1
        weight[i] <- 0 #weight[i] * (length(arrind)/(1+length(arrind)))
      }
      k <- k - 1
    }  
    result[i,] <- x_diff_new 
  }
  return( list(result = result, weight = weight ))
}
# generate_compatible2(c(1,-1,1,1,-1,1,1,1,-1,-1))
aa <- generate_compatible2(c(1,1,-1,1,1,-1), 100000)
aatext <- apply(aa$result,1,paste,collapse="")
table(aatext)
# library(plyr)
ddply(data.frame(aatext=aatext, weight=aa$weight),.(aatext),summarise,freq=sum(weight))

##########################################################################################################
# Thinning #
############
# adlist = a list containing vectors 'arrivals' and 'departures'
# doremove = logical. Are thinned arrivals removed or only marked?
# returnblocks = logical. Is the information on the blocks returned?
#
# Returns a thinned data.table with event times (time) and types (x_diff) 
##########################################################################################################
thin <- function(adlist, doremove = T, returnblocks = T){
  departures <- adlist$departures
  arrivals <- adlist$arrivals[adlist$arrivals <= max(departures)]
  adm <- data.table(time = c(arrivals, departures), 
                    x_diff = c( rep(1, length(arrivals)), rep(-1, length(departures))))
  ram <- adm[ order(-adm$time), ]
  m <- nrow(ram)
  ram$blocksum <- blocksum(ram$x_diff)
  am <- ram[m:1,]
  if(returnblocks) {
    am$blocksumprev <- c(-1, am$blocksum[-m])
    am$blockfirst <- as.numeric( (am$blocksum >= 0) & (am$blocksumprev < 0))
    am$blockind <- cumsum(am$blockfirst)
  }
  if(doremove) {
    am <- subset(am, blocksum != 1)
  } else {  
    am$remove <- (am$blocksum == 1) 
  }
  return(am)
}


#####################################################################################################
# Subroutine for generating departure times given arrivals #
############################################################
# arrivals = a vector of arrival times
# lambda = the intensity (decay rate) parameter of exponential distribution
# p_obsdep = the probability of observing departure
# followup_time = the upper limit for the departure times 
#
# Returns a list containing vectors 'arrivals' and 'departures'
#####################################################################################################

simulate_dep <- function(arrivals, lambda, p_obsdep, followup_time)
{
  narrivals <- length(arrivals)
  lifetime <- rexp(narrivals, lambda)
  departures <- arrivals + lifetime
  departures <- departures[departures <= followup_time]
  obsdep <- as.logical(runif(length(departures)) < p_obsdep)
  return(list(arrivals = arrivals, departures = departures[obsdep]))
}


#################################################################################################################
# Bias correction by noisy binary search #
##########################################
# mu_biased = the biased estimate (raw estimate) for the mean lifetime
# arrivals = a vector of arrival times
# followup_time = the upper limit for the departure times
# p_obsdep = the probability of observing departure
# simNmin = the minimum number of repetations for each candidate estimate 
# simNmin = the maximum number of repetations for each candidate estimate 
# simalpha = the risk level: (1 -simalpha) confidence intervals are used
# tolerance = the relative tolerance (absolute difference divided by mu_biased)
# itermax = the maximum number of rounds in binary search
#
# Returns a bias-corrected estimate of the mean lifetime
#################################################################################################################
correct_bias <- function(mu_biased, arrivals, followup_time, p_obsdep=p_obsdep, 
                         simNmin = 5, simNmax=40, simalpha = 0.01, toler=0.01, itermax=100){
  
  difference <- 999 * mu_biased
  j <- 1
  mu_candidate <- 2*mu_biased
  itermax <- itermax
  lowerbound <- 0
  upperbound <- NA
  repeat {  
    lambda_candidate <- 1/mu_candidate
    estimates <- rep(NA,simNmax)
    for(i in 1:simNmin){
      adlist <- simulate_dep(arrivals = arrivals, lambda = lambda_candidate,
                             p_obsdep = p_obsdep, followup_time = followup_time)
      
      ds <- thin(adlist)
      estimates[i] <- complete_estim( departures = ds$time[ds$x_diff == -1], 
                                      arrivals = ds$time[ds$x_diff == 1])
    }
    difference <- mean(estimates, na.rm = T) - mu_biased
    if((abs(difference) / mu_biased <= toler)) return(mu_candidate) 
    semean <- sd(estimates, na.rm = T)/sqrt(simNmin)
    upperdiff <- difference + qnorm(1 - simalpha/2)*semean
    lowerdiff <- difference - qnorm(1 - simalpha/2)*semean
    simN <- simNmin + 1
    while((simN <= simNmax) & (lowerdiff < 0) & (upperdiff > 0)) {
      adlist <- simulate_dep(arrivals = arrivals, lambda = lambda_candidate,
                             p_obsdep = p_obsdep, followup_time = followup_time)
      ds <- thin(adlist)
      estimates[simN] <- complete_estim( departures = ds$time[ds$x_diff == -1], 
                                         arrivals = ds$time[ds$x_diff == 1])
      difference <- mean(estimates, na.rm = T) - mu_biased
      if((abs(difference) / mu_biased <= toler)) return(mu_candidate) 
      semean <- sd(estimates, na.rm = T)/sqrt(simN)
      upperdiff <- difference + qnorm(1 - simalpha/2)*semean
      lowerdiff <- difference - qnorm(1 - simalpha/2)*semean
      simN <- simN + 1
    }  
    if((abs(difference) / mu_biased <= toler) | (j > itermax) )
    {
      return(mu_candidate)
    }  
    if(is.na(upperbound)) {
      if(difference < 0) {
        lowerbound <- mu_candidate
        mu_candidate <- 2*mu_candidate
      }
      if(difference > 0){
        upperbound <- mu_candidate
        mu_candidate <- (lowerbound + upperbound)/2
      }
    } else {  
      if(difference < 0) {
        lowerbound <- mu_candidate
      }
      if(difference > 0){
        upperbound <- mu_candidate
      }
      mu_candidate <- (lowerbound + upperbound)/2
    }
    j <- j+1
  }
}


##################################################################################################
# Estimator for partially observed departures  (the main function) #
####################################################################
# arrtimes = a vector of arrival times
# deptimes = a vector of departure times
# followup_time = the upper limit for the departure times 
# p_obsdep = the probability of observing departure
# biascorrection = logical. Apply bias correction?
# simNmin = bias correction: the minimum number of repetations for each candidate estimate 
# simNmin = bias correction: the maximum number of repetations for each candidate estimate 
# simalpha = bias correction: the risk level: (1 -simalpha) confidence intervals are used
# tolerance = bias correction: the relative tolerance (absolute difference divided by mu_biased)
# itermax = bias correction: the maximum number of rounds in binary search
#
# Returns a point estimate for the mean lifetime
##################################################################################################
estimate <- function(arrtimes, deptimes, followup_time, p_obsdep=0.5, biascorrection=TRUE,
                     simNmin = 5, simNmax=40, simalpha = 0.01, toler=0.01, itermax=100){
  ams <- thin(list(arrivals = arrtimes, departures = deptimes))
  estimate1 <- complete_estim( departures = ams$time[ams$x_diff == -1], 
                               arrivals = ams$time[ams$x_diff == 1])
  if(biascorrection==TRUE){
    correctedEstimate <- correct_bias(mu_biased=estimate1, 
                                      arrivals = arrtimes, 
                                      p_obsdep=p_obsdep, followup_time = followup_time,
                                      simNmin = simNmin, simNmax = simNmax,
                                      simalpha = simalpha, toler = toler, itermax = itermax)
    return(list(raw_estimate = estimate1, corrected_estimate = correctedEstimate))
  } else {
    return(list(raw_estimate = estimate1, corrected_estimate = NA)) 
  }
}



#################################################################################################################
# Parametric bootstrap #
########################
# Nboot = the number of bootstrap rounds
# mu_estimate = a point estimate for the mean lifetime
# arrivals = a vector of arrival times
# followup_time = the upper limit for the departure times 
# p_obsdep = the probability of observing departure
# simNmin = bias correction: the minimum number of repetations for each candidate estimate 
# simNmin = bias correction: the maximum number of repetations for each candidate estimate 
# simalpha = bias correction: the risk level: (1 -simalpha) confidence intervals are used
# tolerance = bias correction: the relative tolerance (absolute difference divided by mu_biased)
# itermax = bias correction: the maximum number of rounds in binary search
# verbose = logical. Print information on the bootstrap round?
#
# Returns a vector of <Nboot> bootstrap estimates for the mean lifetime
#################################################################################################################
parambootstrap <- function(Nboot=1000, mu_estimate, arrivals, followup_time, p_obsdep=0.5,
                           simNmin = 5, simNmax=100, simalpha = 0.05, toler = 0.001, itermax=100, 
                           verbose = F) {
  mu_boot <- rep(NA, Nboot)
  lambda_estimate <- 1/mu_estimate
  for(j in 1:Nboot) {
    if(verbose) cat("bootstrap round: ", j, "\n")
    adlist <- simulate_dep(arrivals = arrivals, lambda = lambda_estimate,
                           p_obsdep = p_obsdep, followup_time = followup_time)
    ams <- thin(adlist)
    mu_temp <- complete_estim( departures = ams$time[ams$x_diff == -1], arrivals = ams$time[ams$x_diff == 1])
    mu_boot[j] <- correct_bias(mu_biased=mu_temp, arrivals = arrivals, 
                               p_obsdep=p_obsdep, followup_time=followup_time,
                               simNmin = simNmin, simNmax=simNmax,
                               simalpha = simalpha, toler=toler, itermax=itermax)
  }
  return(mu_boot)
}  


#################################################################################################################
#################################################################################################################
#################################################################################################################
# 2) Functions for simulation experiments #
###########################################


#################################################################################################################
# Function to simulate arrivals and departures #
################################################
# activity_time = the length of the activity period
# lambda_arr = the arrival rate
# lambda = the intensity (decay rate)
# p_obsdep = the probability of observing departure
# followup_time = the upper limit for the departure times 
#
# Returns a list containing vectors 'arrivals' and 'departures'
#################################################################################################################
simulate_arr_dep <- function(activity_time, lambda_arr, lambda, p_obsdep, followup_time)
{
  narrivals <- rpois(1, activity_time * lambda_arr)
  arrivals <- runif(narrivals, 0, activity_time)
  lifetime <- rexp(narrivals, lambda)
  departures <- arrivals + lifetime
  departures <- departures[departures <= followup_time]
  obsdepind <- as.logical(rbinom(length(departures), 1, p_obsdep))
  obsdep <- departures[obsdepind]
  return(list(arrivals = arrivals, 
              departures = obsdep, 
              ndep_during = sum(obsdep <= activity_time), 
              ndep_after = sum(obsdep > activity_time)))
}

#################################################################################################################
# Function to carry out a simulation experiments  #
###################################################
# simsetup = a data.table or data.frame with the variables: 
#             lambda_arr = the arrival rate
#             activity_time = the length of the activity period
#             followup_time = the upper limit for the departure times 
# reps = the number of simulation rounds
# seed = the seed for the random number generator
#
# Returns a data.table or data.frame with the variables:
#             lambda_arr = the arrival rate
#             activity_time = the length of the activity period
#             followup_time = the upper limit for the departure times 
#             raw_estimate = an estimate for the mean lifetime without bias correction
#             corrected_estimate = an estimate for the mean lifetime with bias correction
#################################################################################################################
simulate <- function(simsetup, reps = 1, seed = NULL) {
  if(!is.null(seed)) set.seed(seed)
  simresults <- simsetup[rep(1:nrow(simsetup), reps), ]
  simresults$raw_estimate <- NA
  simresults$corrected_estimate <- NA
  for(j in 1:nrow(simresults)) {
    adlist <- simulate_arr_dep(activity_time = simresults$activity_time[j],
                               lambda_arr = simresults$lambda_arr[j],
                               lambda = 1, p_obsdep = 0.5, 
                               followup_time = simresults$followup_time[j])
    estim <- estimate(arrtimes=adlist$arrivals, deptimes=adlist$departures, 
                      followup_time=simresults$followup_time[j], 
                      biascorrection=T,
                      simalpha = 0.05, simNmax = 100, toler = 0.001)
    simresults$raw_estimate[j] <- estim$raw_estimate
    simresults$corrected_estimate[j] <- estim$corrected_estimate
    simresults$ndep_during[j] <- adlist$ndep_during
    simresults$ndep_after[j] <- adlist$ndep_after
  }  
  return(simresults)
}  

#################################################################################################################
# Function to carry out a simulation experiment (The second simulation example in the paper)  #
###############################################################################################
# reps = the number of simulation rounds
#
# Returns a data.table or data.frame with the variables:
#             lambda_arr = the arrival rate
#             activity_time = the length of the activity period
#             followup_time = the upper limit for the departure times 
#             raw_estimate = an estimate for the mean lifetime without bias correction
#             corrected_estimate = an estimate for the mean lifetime with bias correction
simulation_example <- function(reps = 1) {
  lambda_arr <- c(0.1,0.2,0.5,1,2,5,10,20,30,40,50,60,70,80,90,100,200,500,1000)
  simsetup <- data.table(lambda_arr = rep(lambda_arr,3), 
                         activity_time = rep(c(3000,2000,1000), each = length(lambda_arr)))
  simsetup$followup_time <- simsetup$activity_time + 100
  result <- simulate(simsetup, reps = reps)
  return(result)
}
# Run the example: 
# simexample <-simulation_example(1)
# Consider using parallel computing for large simulations.



# 3) A simulation experiment to demonstrate  the inconsistency of estimator (5) in the article
simulate_one5 <- function(n = 10000, lambda_arr = 1, lambda = 1) {
  x <- rexp(n,rate = lambda)
  a <- cumsum(rexp(n, rate = lambda_arr))
  d <- a + x 
  obsdepind <- as.logical(rbinom(length(d), 1, 0.5))
  dobs <- d[obsdepind]
  return((2*sum(dobs) - sum(a))/n)
}  

simulate_estimator5 <- function(reps = 10, 
                                nvec = c(100, 1000, 10000, 100000)) {
  nn <- length(nvec)
  results <- data.frame(n = rep(nvec, reps), estim = NA)
  for(i in 1:reps) {
    for(j in 1:nn) {
      results$estim[(i-1)*nn + j] <- simulate_one5(n = nvec[j])
    }
  }  
  return(results)
}  

# results5 <- simulate_estimator5()
# with(results5, plot(n, estim, log = "x"))
