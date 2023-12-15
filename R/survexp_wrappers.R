# library(survival)
# library(data.table)
# library(tidyverse)
# library(dtplyr)

# Rate table
# building + access
# check that all ratetable variables are present in (or mapped to) the survival dataset

# Compute population survival based on the ratetable
# function: startdate, enddate, covars real, ratetable unit
# hack the nessie function from relsurv package (in order to remove the absurd print)?
# or simply use the survexp function from survival package (but I'm not quite sure this is what the function is doing)

# test = slopop
# patients_copy=copy(patients)
# survexp(formula = c(20,20) ~ 1, data=patients_copy[1:2],
#                  rmap = list(year=date_chir, age=age_at_dg, sex=sexe),method = 'individual.s', ratetable = test)
# 
# survexp(formula = ~ 1, data=patients_copy[1],
#                  rmap = list(year=date_chir, age=age_at_dg, sex=sexe),method = 'ederer',times=c(10,20), ratetable = test)
# seems I can either compute several times for a single patient or one date for several patients at a time, not both

# create some kind of class for ratetableDT to ensure time windows are contiguous, ensure DT has correct columns
# use YoB instead of age if time varying hazards

# from survexp.us: All numeric dimensions of a rate table must be in the same units.
# The survexp.us rate table contains daily hazard rates, the age cutpoints are in days, and the calendar year cutpoints are a Date.
# I'm still confused about how the implementation could know when times are given in days or anything else, and how it uses Dates accordingly...
# I will simply make the hypothesis that this is all very sloppy and that the only accepted time unit is 'day'

# S code to create a ratetable
# temp <- matrix(scan("data.smoke"), ncol=8, byrow=T)/100000smoke.rate <- c(rep(temp[,1],6), rep(temp[,2],6), temp[,3:8])attributes(smoke.rate) <- list(dim=c(7,2,2,6,3),dimnames=list(c("45-49","50-54","55-59","60-64","65-69","70-74","75-79"),c("1-20", "21+"),c("Male","Female"),c("<1", "1-2", "3-5", "6-10", "11-15", ">=16"),c("Never", "Current", "Former")),dimid=c("age", "amount", "sex", "duration", "status"),factor=c(0,1,1,0,1),cutpoints=list(c(45,50,55,60,65,70,75),NULL, NULL,c(0,1,3,6,11,16),NULL),class='ratetable')is.ratetable(smoke.rate)

# The most important point is to note that age has been rescaled.  This table contains ratesper year, whereas the US tables contained rates per day.  It is crucial thatall  of  the  time  variables  (age,  duration,  etc)  be  scaled  to  the  same  units,or the results may not be even remotely correct.  The US rate tables werecreated using days as the basic unit since year of entry will normally be ajulian date; for the smoking data years seemed more natural.
# A lot of this hastle could be removed using proper classes for times such as provided by lubridate

# For a general ratetable the only sensible thing to do might be to force to have 1 time varying column (because more than 1 does not really make sense) and force it to be -Inf + Inf if not useful

# Compute indiv cumulated hazard based on ratetable with rates depending on times


#' @title Compute expected survival for one individual at several times
#' @description A wrapper function calling survival::survexp to compute expected survival probablity of a single individual for one or several times.
#' @param individual_df A DF or DT containing a single row with a single patient's covariates
#' @param eval_times Times at which the expected survival should be evaluated. Times must be in days if using a ratetable with a Date component. 
#' @param ratetable A ratetable object 
#' @param rmap An rmap argument as the one passed to survival::survexp. This argument will only be evaluated inside survexp
#' @param fast Binary, if TRUE the function will omit some validity checks on the arguments. Should be used with caution.
#'
#' @return Returns a tibble with two columns: eval_times and surv the survival probability from the survexp object at those times. Upon failure of surv.exp because of missing data in the DF row, a vector of NA will be returned. 
#'
#' @example indiv_survprob_pch(patients[1],eval_times = c(10,20,30),ratetable = slopop,rmap=list(year=date_chir, age=age_at_dg, sex=sexe))
indiv_survprob_pch<-function(individual_df,eval_times,ratetable,rmap,fast=FALSE){
  
  # # pch for piecewise constant hazard
  if(!fast){
    
  if(!all(is.numeric(eval_times))) stop("Numeric values are expected for relative `eval_times`.")
  if(!all(eval_times>0)) stop("Relative `eval_times` must be greater than 0")
  # initially I wanted to allow the possibility to give absolute times and not only relative ones
  # however this needed to fiddle with the rmap object and find the correct date argument etc
  # in the end this is really not worth it
  # same goes for specifying a start time (s time)
  
  if(nrow(individual_df)>1){
    stop("`individual_df` should contain a single row (1 subject and its covariates)")
  }
  }
  
  # eval_times can be relative and must be expressed in days
  exp.surv<- tryCatch(
  eval(substitute(survival::survexp(formula = ~ 1, data=individual_df,
                                              rmap = rmapsub ,method = 'ederer',times=c(eval_times), ratetable = ratetable,na.action=na.omit),
                 list(rmapsub=substitute(rmap))))$surv,
  error = function(err) {
    if(err$message =="Data set has 0 rows"){
      return(rep(NA,length(eval_times)))
    }
    else{
      message("Exception caught upon calling survival::survexp() in indiv_survprob_pch().")
      message(paste0("Caught on individual_df: ",paste(individual_df,collapse = " ")))
      return(stop(err))
    }
  }
  )
  return(tibble(eval_times=eval_times, surv=exp.surv))
}


#' @title Apply individual surv prob computation to DF
#' @description This function is a wrapper around `indiv_survprob_pch`. It calls the latter function on every row of the DF in order to compute the expected survival probability for each individual row at each provided `eval_times`. The point of this is to call survexp only once per patient and compute the survival probablity at all patime points for each patient, thus avoiding redundant computation.
#' @param patientsDF a DF/DT/tibble containing the data needed to compute individual expected survival probability
#' @param eval_times The times at which survival must be evaluated. Times must be in days if using a ratetable with a Date component. Provided times can be passed as single numeric value or vector of numeric values or as unamed list of numeric values if common for all individuals. Otherwise the name of the column containing the list,vector or single value evaluation point can be passed as character string or as variable name.
#' @param ratetable A ratetable object. See ?survival::ratetable for details.
#' @param rmap An rmap argument as the one passed to survival::survexp. This argument will only be evaluated inside survexp
#'
#' @return Returns a tidy tibble with 3 columns: the row name of the individual row (as given by the call of row.names() on the patientsDF), the evaluation time and the corresponding survival probability. 
#' @export
#'
#' @examples compute_survprob_pch(patients,eval_times = c(10,20,30),ratetable = slopop,rmap=list(year=date_chir, age=age_at_dg, sex=sexe))
compute_survprob_pch<-function(patientsDF,eval_times,ratetable,rmap){
  patientsDF<-as_tibble(patientsDF)
  patientsDF<- patientsDF %>% 
    tibble::add_column(SQVVcCs1lD4R7tDVlOoVrowid=row.names(patientsDF),.name_repair = "check_unique") # Use a random col name, this will throw an error in case the column already exists
  
  if(purrr::is_scalar_character(eval_times)){
    eval_times<-patientsDF[,eval_times] # this doesn't make sense at all! FIXME
    # I might have to evaluate it at that point then?
  }
  else if(is.numeric(eval_times)){
    # Make sure it's not a date # FIXME
    # list of list of times or list of vectors, need to have an id column for this?
    if(length(eval_times)==1) {
      patientsDF<-dplyr::mutate(patientsDF,SQVVcCs1lD4R7tDVlOoVeval_times=eval_times)
    }
    else{
      old_cols<-colnames(patientsDF)
      for (t in eval_times) {
        patientsDF<-patientsDF%>% add_column("SQVVcCs1lD4R7tDVlOoVeval_times{{t}}":=t,.name_repair = "universal")
      }
      patientsDF<- patientsDF %>% pivot_longer(-all_of(old_cols),names_to = NULL ,values_to = "SQVVcCs1lD4R7tDVlOoVeval_times")
    }
 
  }
  else if(is.list(eval_times)){
    is_num<-lapply(eval_times, is.numeric)
    lens<- lapply(eval_times, length)
    is_list <-lapply(eval_times, is.list)
    
    if(all(is_num) & all(lens==1) & is.null(names(eval_times))){
      # the provided list is an unnamed list of length one numeric values
      patientsDF<-dplyr::mutate(patientsDF,SQVVcCs1lD4R7tDVlOoVeval_times=list(eval_times)) 
    }
    else if(all(names(eval_times %in% row.names(patientsDF)))){
      # TODO finish this list distribution
    }
    patientsDF<-patientsDF %>% tidyr::unnest_longer(SQVVcCs1lD4R7tDVlOoVeval_times) 
  }
  else{
    stop("`eval_times` must be a numeric vector, a list ,a string designating the dataframe column to be used.")
  }
  # TODO eventually add array of arrays
  # TODO output a period object? (such that the output contains all the information)  
  enquo_rmap<-enexpr(rmap)
    
    exp.surv<-eval(substitute(survival::survexp(formula = SQVVcCs1lD4R7tDVlOoVeval_times ~ 1, data=patientsDF,
                    rmap = rmapsub ,method = 'individual.s', ratetable = ratetable,na.action=na.exclude),
  list(rmapsub=substitute(rmap))))
    
    final_df<- patientsDF %>% select(SQVVcCs1lD4R7tDVlOoVrowid,SQVVcCs1lD4R7tDVlOoVeval_times)%>%
      mutate(expsurvs=exp.surv) %>%  
    rename(row.name=SQVVcCs1lD4R7tDVlOoVrowid,eval_time=SQVVcCs1lD4R7tDVlOoVeval_times)
  
  return(final_df)
}


dftoRatetable<-function(lambda,data=popmortality_df,qty=c('survival','mortality','yearlyHaz','dailyHaz'),rmap){
  # mortality: mortality between two time points
  # survival: conditionnal probability of survival for the time interval
  # yearlyHaz: yearly hazard rate (cumulated hazard over a year)
  # dailyHaz: daily hazard (cumulated hazard over a day)
  
  rmap['calendar_year'] # check that e
  
  # Create a ratetable containing dept information (this is incredibly and absurdly complicated)
  # adapted from: https://stackoverflow.com/questions/31064599/make-the-rate-table-for-relative-survival-analysis
  dept_ratetable_list = list()
  depts = levels(lifetable$dept)
  for (dep in depts){
    dept_ratetable_list = c(dept_ratetable_list,
                            list(transrate(lifetable[sexe=="M" & dept==dep,.(agerev,annee,yearsurv=1-mua)]%>%pivot_wider(id_cols=agerev,names_from = annee, values_from=yearsurv)%>%select(!agerev)%>%as.matrix(),
                                           lifetable[sexe=="F" & dept==dep,.(agerev,annee,yearsurv=1-mua)]%>%pivot_wider(id_cols=agerev,names_from = annee, values_from=yearsurv)%>%select(!agerev)%>%as.matrix(),
                                           yearlim=c(lifetable[dept==dep,min(annee)],lifetable[,max(annee)]),
                                           int.length=1)) 
    )
  }
  # TODO rewrite this not using transrate since transrate limits to evenly spread year intervals
  names(dept_ratetable_list)<-depts
  rt_lifetable = joinrate(dept_ratetable_list,dim.name = "dept")
  
}

# Proper alternative format to ratetable:
# DT with yearBirth added column instead of age, key on it and other covariates
# apply dumn function: return 0 if less than date s, return fraction if 

# Prepare data => cf table 11.2 from handbook of survival analysis P.231
# expand in dummy variables (cf Layla's code)
# + need for actual dates to compute population survival

# Estimation of censoring probabilities
# cf ipcw.* functions from riskRegression package
# formula for censoring dist
# should be done once for all transitions
# add one covariate for "illness" state? (is censoring different?)

# Regression
# use GLM or or geese (or wglm from risk Regression package? => but does it allow to start from time "s"?)
# allow fitting GLM's for several time points with constant time parameters
# => this means I need to compute pop survival and censoring prob for each individual before calling the GLM utility
# and store them in a table (I most probably will need some pivot_wider/longer etc magic)

# 2 options to include the net survival estimate in the computation:
# rescale the predictor and the weights inversely
# or play with offset but this will be much dependent on the link function (???)

# One must be careful to the weighting of the different time points for constant parameters
# This is briefly mentionned in Handbook of survival analysis p.228 and seems to be taken care of
# by the 'id' field of geese/geeglm p.231

# should I foresee issues in learning all parameters at the same time even when all parameters are time independent 
# (in terms of optimization)?
# MUST HAVE A TIME DEPENDENT INTERCEPT

# 

# Bootstrap: in theory censoring probability estimation should be bootstrapped too, but expected survival should be only computed once
# this means the bootstrap ~~loop~~ or function should be on top level and apply to all transitions if we want to avoid
# multiplying by the number of transitions the amount of computation time for estimating censoring distribution
# this means one needs a wrapper function (and not a loop) to avoid code duplication, since it must be called once on the complete
# dataset and then for bootstrap estimates

# Tests for time dependence of parameters, and significantly departing from 0

# Simulation function

# function to map times to last event times and vice versa
# TODO take care of t0 which may not correspond to an event time?

#' @title Map time steps to event times
#' @description Filter and map a list/vector of time steps to a list/vector of event times.
#'  This is useful for evaluating a piecewise constant function (e.g Kaplan-Meier survival, binomial regression)
#'  but limit the number of redundant computation steps, by evaluating the function only at breakpoints.
#' @param time_steps An array/vector of desired timesteps 
#' @param event_times An array/vector of actual event-times or breakpoints
#' @return An array of unique time points ordered in increasing order
#' @export FALSE
#'
#' @examples
#' .filter_times(c(1,2,4,5,6),c(1,3,5))
.filter_times<-function(time_steps,event_times){
  if(min(time_steps)<min(event_times) ){
    stop("All time steps must be >= min(event_times)")
  }
  #sort vector in decreasing order
  time_steps <- time_steps[order(time_steps,decreasing=TRUE)]
  # sort event_times in increasing order
  event_times<-event_times[order(event_times)]
  mapped_event_times <- c() 
  ind=0
  for(t in time_steps){
    # Find last event time (first backward) inferior or equal to t
    ind = purrr::detect_index(event_times,function(x){x<=t},.dir = "backward")
    mapped_event_times <- c(mapped_event_times,event_times[ind])
    event_times<-event_times[1:ind] # remove values than cannot fulfill the condition anymore (speedup?, depends on R implementation)
  }
  res <- unique(mapped_event_times)
  return(res[order(res)])
}

#FIXME unfinished function
.remap_values_to_times<-function(val,filtered_event_times,time_steps){
  # take event times values and return a vector of values mapping to time_steps
  if(length(val)!=length(filtered_event_times)){
    stop("Length of val and filtered_event_times mismatch, each val should correspond to a filtered_event_time")
  }
  if(min(time_steps)<min(filtered_event_times) | max(time_steps)>max(filtered_event_times)){
    stop("All time steps must lie in [min(event_times),max(event_times)]")
  }
}



# nessie function from the relsurv package, under GNU GPL
nessie <- function (formula = formula(data), data = parent.frame(), ratetable = relsurv::slopop, 
          times, rmap) 
{
  call <- match.call()
  if (!missing(rmap)) {
    rmap <- substitute(rmap)
  }
  na.action <- NA
  rform <- rformulate(formula, data, ratetable, na.action, 
                      rmap)
  templab <- attr(rform$Terms, "term.labels")
  if (!is.null(attr(rform$Terms, "specials")$ratetable)) 
    templab <- templab[-length(templab)]
  nameslist <- vector("list", length(templab))
  for (it in 1:length(nameslist)) {
    valuetab <- table(data[, match(templab[it], names(data))])
    nameslist[[it]] <- paste(templab[it], names(valuetab), 
                             sep = "")
  }
  names(nameslist) <- templab
  data <- rform$data
  p <- rform$m
  if (p > 0) {
    data$Xs <- my.strata(rform$X[, , drop = F], nameslist = nameslist)
  }
  else data$Xs <- rep(1, nrow(data))
  if (!missing(times)) 
    tis <- times
  else tis <- unique(sort(floor(rform$Y/365.241)))
  tis <- unique(c(0, tis))
  tisd <- tis * 365.241
  out <- NULL
  out$n <- table(data$Xs)
  out$sp <- out$strata <- NULL
  for (kt in order(names(table(data$Xs)))) {
    inx <- which(data$Xs == names(out$n)[kt])
    temp <- exp.prep(rform$R[inx, , drop = FALSE], rform$Y[inx], 
                     rform$ratetable, rform$status[inx], times = tisd, 
                     fast = FALSE)
    out$time <- c(out$time, tisd)
    out$sp <- c(out$sp, temp$sis)
    out$strata <- c(out$strata, length(tis))
    temp <- exp.prep(rform$R[inx, , drop = FALSE], rform$Y[inx], 
                     rform$ratetable, rform$status[inx], times = (seq(0, 
                                                                      100, by = 0.5) * 365.241)[-1], fast = FALSE)
    out$povp <- c(out$povp, mean(temp$sit/365.241))
  }
  names(out$strata) <- names(out$n)[order(names(table(data$Xs)))]
  if (p == 0) 
    out$strata <- NULL
  mata <- matrix(out$sp, ncol = length(tis), byrow = TRUE)
  mata <- data.frame(mata)
  mata <- cbind(mata, out$povp)
  row.names(mata) <- names(out$n)[order(names(table(data$Xs)))]
  names(mata) <- c(tis, "c.exp.surv")
  cat("\n")
  print(round(mata, 1))
  cat("\n")
  out$mata <- mata
  out$n <- as.vector(out$n)
  class(out) <- "nessie"
  invisible(out)
}

