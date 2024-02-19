# Code adapted from idmTPreg package v1.1 under GNU-GPLv2
# Original Authors: Leyla Azarang and Manuel Oviedo de la Fuente
# Adapted for net survival setting by Quentin Marcou, based on <ref>

# a relative survival logit link function for 11,12,22 transitions
# rellogit <- function(t) {
#   SL= Survival(t)
#   linkfun <- function(miu) log((miu/SL)/abs(1-(miu/SL)))
#   linkinv <- function(et)  SL*exp(et)/(1+exp(et))
#   mu.eta <- function(et) {SL*exp(et)/(1+exp(et))^2  }
#   valideta <- function(et) TRUE
#   link <- "log((miu/SL)/(1-(miu/SL)))"
#   structure(list(linkfun = linkfun, linkinv = linkinv,
#                  mu.eta = mu.eta, valideta = valideta,
#                  name = link),
#             class = "link-glm")
# }
#
# # an offset survival logit link function for 13 and 23 transitions
# offsetlogit <-function(t) {
#   dp =1- Survival(t) #death probability
#   SL = Survival(t)
#   linkfun <- function(miu) log((miu-dp)/abs(1-(miu-dp)))
#   linkinv <- function(et)  (exp(et)*SL+dp)/(1+exp(et))
#   mu.eta <- function(et) {((2*SL-1)*exp(et))/((1+exp(et))^2) }
#   valideta <- function(et) TRUE
#   link <- "log((miu-dp)/(1-(miu-dp)))"
#   structure(list(linkfun = linkfun, linkinv = linkinv,
#                  mu.eta = mu.eta, valideta = valideta,
#                  name = link),
#             class = "link-glm")
# }


#' @title Wraps the call to mod.glm.fit to handle convergence errors
#' @description
#' Tries to handle convergence and other issues upon calling the mod.glm.fit
#' function. Some errors returned by mod.glm.fit are cast to warnings. In order
#' to reduce crashes several values of `maxit` are tested if convergence issues
#' are detected. NAs will be returned in case of crash/divergence or lack of
#' informative observations.
#'
#' @param X
#' @param response
#' @param family
#' @param weights
#' @param maxit
#' @param maxmaxit
#' @param warning_str
#' @param ...
#'
#' @return A glm object
#'
#' @examples
mod.glm.fit.errorwrapper<-function(X,response,family,weights,maxit=glm.control()$maxit,maxmaxit=1000,warning_str="",...){
  result <-tryCatchLog::tryCatchLog(
    {
      #Try
      withCallingHandlers({
        mod.glm.fit2(X, response, family = family, weights = weights,start = rep(0,ncol(X)),control = glm.control(maxit = maxit))
      },
      warning=function(warn){

        if(stringr::str_detect(warn$message,"no observations informative at iteration")){
          stop(paste0("Warning caught ",warning_str,": ",warn$message," returning NA"))
        }
        else{
          warning(paste0(warning_str,warn$message))
        }
      }
      )

    },
    error=function(err) {
      # Catch
      if(err$message =="inner loop 1; cannot correct step size" ||
         err$message =="inner loop 2; cannot correct step size"){
        if(maxit*10<=maxmaxit){
          # recursively call the wrapper with a greater maxit
          return(mod.glm.fit.errorwrapper(X=X, response=response, family = family, weights = weights,
                                          maxit = maxit*10, maxmaxit = maxmaxit,warning_str = warning_str))
        }
        else{
          warning(paste0("Step size correction issue with maxmaxit reached: ",warning_str,"returning NA as result"))
          tmp<-rep(NA,dim(X)[[2L]])
          names(tmp)<-dimnames(X)[[2L]]
          return(list(coefficients=tmp,converged = FALSE))
        }
      }
      else if(stringr::str_detect(err$message,"no observations informative at iteration")){
        tmp<-rep(NA,dim(X)[[2L]])
        names(tmp)<-dimnames(X)[[2L]]
        return(list(coefficients=tmp,converged = FALSE))
      }
      else if(stringr::str_detect(err$message,"NA/NaN/Inf in 'y'")){
        # Handle an error that i do not fully understand, this seems to be linked to the divergence of the algorithm
        # with derivatives regarding parameters exploding. I do not know whether this could be linked to Layla's
        # modification to the glm.fit code
        warning("mod.glm.fit crashed with error \"",err$message,"\"",warning_str,", returning NA")
        tmp<-rep(NA,dim(X)[[2L]])
        names(tmp)<-dimnames(X)[[2L]]
        return(list(coefficients=tmp,converged = FALSE))
      }
      else{
        message(paste0("Exception caught upon calling modl.glm.fit",warning_str))
        stop(err$message)
      }
    }
  )
  return(result)
}


#' @title Check input and output of mod.glm.fit
#' @description
#' Avoids calling the glm fitting function when there are no contrasted response
#' values (all 0 or all 1). Ensures return of NA coefficient when the algorithm
#' did not converge.
#'
#'
#' @param X
#' @param response
#' @param family
#' @param weights
#' @param maxit
#' @param maxmaxit
#' @param warning_str
#' @param ...
#'
#' @return A set of glm coefficients
#'
#' @examples
mod.glm.fit.callingwrapper<-function(X,response,family,weights,maxit=glm.control()$maxit,maxmaxit=1000,warning_str="",...){
  if(any(response) & !all(response)){
    # There must be at least one event in the sample in order to learn smthg
    result<-mod.glm.fit.errorwrapper(X=X, response=response, family = family, weights = weights,
                                     maxit = maxit, maxmaxit = maxmaxit,warning_str = warning_str)
    #print(result$converged)
    if(!result$converged || result$boundary){
      result$coefficients=result$coefficients*NA #set coef to NA if the algorithm did not converge
    }
    return(coefficients(result))
  }
  else{
    # if no event, coefficients are meaningless and one should return NA
    warning("All provided responses are equal",warning_str,", cannot fit GLM, returning NA")
    tmp<-rep(NA,dim(X)[[2L]])
    names(tmp)<-dimnames(X)[[2L]]
    return(tmp)
  }
}



renewnetTPreg <-
function(formula, data, ratetable, link,rmap,time_dep_popvars=list('year','age'), s = 0, t = NULL,R = 199, by = NULL, trans)
{
# Dictionnary of used variables:
	# X: the model matrix, created from the data given the formula, model.matrix expands factors in dummy variables
	# comdata: "complete" data (no NA in any column), columns are ordered in a certain way, TODO stop creating dumb variables (ordata,comdata) and just edit the data variable
	# data: initially the whole data, then comdata subsets (data1, data2). Columns: Tt (\tifle(T),Zt(\tilde(Z)))
	# 1/2 suffix: patients in state 1 or 2 at time s
	# SfitXY: a Survfit object for the transition between state X and Y
	# ShatXY: a dataframe containing survival estimate at different times (starting at time 0)
	# SfitX: corresponds to ShatXX
	# Sfit: global censoring distribution, with a really great name!
       	# vec.tXY:

  if (missing(data))
    stop("Argument 'data' is missing with no default")
  if (!is.data.frame(data))
    stop("Argument 'data' must be a data.frame")
    if (sum(! c("id","Zt","Tt","delta1","delta", "age" ) %in% (colnames(data)))>0) # at least age should be included
        stop("data should  contain  id, Zt, Tt, delta1, delta, age variables") # TODO later on add also "sex" covariable
  if (ncol(data)<=4)
    stop("'data' must have covariables")
  if(sum(is.na(data))!=0){
    miscolumn <- sapply(data, function(x) sum(is.na(x)))
    miscolnam <- names(miscolumn[miscolumn!=0])
    if(sum(miscolumn)!=0)
      warning(sapply(miscolnam, function(x)paste(x," variable in 'data' has missing value(s)", ", ",sep="")))
  }


 # This piece of code seems to be used to make sure the formula and data columns match
 # It then builds a model.matrix object (which allows automatic expansion in dummy variables?)
 # Still I do not understand all this fiddling with getting the function call (espcially since there's no ... argument in the function)
 # why not use directly stats::model.frame(formula,data) ?
  formula <- formula # wtf?
  cl <- match.call() #get the function call (or ""command line"") string with named arguments
  mf <- match.call(expand.dots = FALSE) # same without expending ... argument, the idea is to reuse the function call, and substitute the function by model.frame
  m  <- match(c("formula", "data"), names(mf), 0L) # look for "formula" and "data" in arguments outside '...' and return their indices, return 0L (=0) if not found
  mf <- mf[c(1L, m)]
  # add options for model.frame()
  mf$drop.unused.levels <- TRUE # simplify factors and retain only used levels
  mf$na.action <- na.pass # all missing values will be retained in the model frame
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  if (is.empty.model(mt)) {
    stop("'formula' must match with 'data'")
  }
  else{ #FIXME useless else since if is a stop
    X<- model.matrix(mt, mf, contrasts)
    ind = match(colnames(X) , colnames(data))
  ind = ind[!is.na(ind)]
  covname= colnames(data)[ind]
  }

  # Filter useful columns based rmap and formula arguments
  formnames<- all.vars(formula,functions = FALSE,unique = TRUE)
  rmapnames<- all.vars(substitute(rmap),functions = FALSE,unique = TRUE)
  ordata = data[, unique(c(c("id", "Zt", "Tt", "delta1", "delta"),rmapnames, formnames ),fromLast=FALSE)]

  L.or <- nrow(ordata)
  #Remove lines with at least one missing value
  comdata <- ordata[complete.cases(ordata),]
  X<-X[complete.cases(ordata), ,drop=FALSE]

  L.com <- nrow(comdata)
  n.misobs <- L.or - L.com
  if(is.null(by)){
    by <- floor((max(comdata$Zt) - min(comdata$Zt))/quantile(comdata$Zt,0.01))
  }

  if(is.null(t)){
    t = max(comdata$Zt[comdata$delta1 == 1], na.rm = T)
  }

  if(t <= s || s<0){
    stop("argument 's' must be smaller than 't' and larger than 0")
  }
  else{
    # FIXME Why is this all contained in an "else" statement?
    if (!(link %in% c("logit", "probit", "cauchit")))
      stop( paste("binomial family do not have", "'", link, "'",  "link"))

    if(sum(comdata$delta1 < comdata$delta) != 0){
      stop("'delta' must be 0 when 'delta1' is 0")
    }
    if (!(trans %in% c("11", "12", "13","23","all")))
      stop(paste(trans, "is not a valid transition for a progressive illness-death model"))

    if(s == 0 & (trans == "23" || trans == "all" )) # TODO there's probably better way around
      stop("for the transition '23' argument 's' must be larger than 0")
    if (trans == "23" || trans == "all" ){
      X2 <- X[comdata$Zt<=s & comdata$Tt>s , ,drop=FALSE] # patients in state 2 at time s
      data2 <- comdata[comdata$Zt <= s & comdata$Tt > s,]
      Sfit23 <- summary( survfit(Surv(data2$Tt, data2$delta == 0)~ +1))
      Shat23 <- rbind(c(0,1),data.frame(time = Sfit23$time, surv = Sfit23$surv))
      Shat.function23 <- function(x){
	# A shitty function returning the last Shat$surv value known before time x
        Shatx23 <- if(length(Shat23$time) > 1)  tail(subset(Shat23, time <= x)$surv,1)
        else 1
        return(Shatx23)
      }
    }
    if (trans=="11" ||trans=="12" ||trans=="13" || trans=="all"){
      X <- X[comdata$Zt > s, ,drop=FALSE] # FIXME rename to X1 in order to remain consistent with X2/data2/... and data1
      data1 <- comdata[comdata$Zt > s,]
      Sfit1 <- summary( survfit(Surv(data1$Zt, data1$delta1 == 0)~ +1))
      Shat1 <- rbind(c(0,1), data.frame(time = Sfit1$time, surv = Sfit1$surv))
      Shat.function1 <- function(x){
	# A shitty function returning the last Shat$surv value known before time x
	Shatx1 <- if(length(Shat1$time) > 1)  tail(subset(Shat1,time <= x)$surv, 1)
        else 1
        return(Shatx1)
      }
      Sfit <- summary( survfit(Surv(data1$Tt, data1$delta == 0)~ +1)) # Note delta==0 => means we are learning the censoring distribution here
      Shat <- rbind(c(0,1), data.frame(time = Sfit$time, surv = Sfit$surv))
      Shat.function <- function(x){
	# A shitty function returning the last Shat$surv value known before time x
        Shatx <- if(length(Shat$time) > 1)  tail(subset(Shat, time <=x )$surv, 1)
        else 1
        return(Shatx)
      }
    }

    co <- vector("list", 4)
    names(co) <- c("co11", "co12", "co13", "co23")

    # Check correctness of time_dep_popvars and its interplay with rmap
    rmapsub<-substitute(rmap)
    if(is.list(time_dep_popvars) || is.character(time_dep_popvars)){
      # FIXME add check for list content being strings
      if(!all(time_dep_popvars %in% names(rmapsub)[-1])){
        stop("Names in `time_dep_popvars` do not correspond to ratetable dimension names in the `rmap` argument.")
      }

      # Now change rcall to account for s days time shift in survival computation
      for (var in time_dep_popvars){
        rmapsub[[var]]<-call('+',rmapsub[[var]],s)
      }

    }
    else if(!is.null(time_dep_popvars)){
      stop("`time_dep_popvars` must be a list or vector of strings or NULL")
    }


    # Shit part to compute expected survival

    # rellogit <- function(t,data_df) {
    #   SL <- eval(substitute(compute_survprob_pch(data_df,t-s,ratetable,rmap=rmapsubs),list(rmapsubs=rmapsub)))$expsurvs
    #   linkfun <- function(miu) log((miu/SL)/abs(1-(miu/SL)))
    #   linkinv <- function(et)  SL*exp(et)/(1+exp(et))
    #   mu.eta <- function(et) {SL*exp(et)/(1+exp(et))^2  }
    #   valideta <- function(et) TRUE
    #   link <- "log((miu/SL)/(1-(miu/SL)))"
    #   structure(list(linkfun = linkfun, linkinv = linkinv,
    #                  mu.eta = mu.eta, valideta = valideta,
    #                  name = link),
    #             class = "link-glm")
    # }
    #
    # # an offset survival logit link function for 13 and 23 transitions
    # offsetlogit <-function(t,data_df) {
    #   SL <- eval(substitute(compute_survprob_pch(data_df,t-s,ratetable,rmap= rmapsubs),list(rmapsubs=rmapsub)))$expsurvs
    #   dp <- 1- SL #death probability
    #   linkfun <- function(miu) log((miu-dp)/(1-(miu-dp)))
    #   linkinv <- function(et)  (exp(et)*SL+dp)/(1+exp(et))
    #   mu.eta <- function(et) {((2*SL-1)*exp(et))/((1+exp(et))^2) }
    #   valideta <- function(et) TRUE
    #   link <- "log((miu-dp)/(1-(miu-dp)))"
    #   structure(list(linkfun = linkfun, linkinv = linkinv,
    #                  mu.eta = mu.eta, valideta = valideta,
    #                  name = link),
    #             class = "link-glm")
    # }

# Look at the 1->1 transition
    if (trans == "11" || trans == "all") {
      print("Enter 11")
      vec.t11 <- data1$Zt
      M11 <- max(vec.t11)
      vec.t11 <- vec.t11[order(vec.t11[vec.t11 > s])]# dafuq? #FIXME
      # vec.t11 times are already guaranteed to be >s, since data1 contains only obs for Zt >s
      vec.t11 <- vec.t11[vec.t11 <= t]
      vec.t11 <- vec.t11[seq(1, length(vec.t11), by)]
      vec.t11 <- c(vec.t11, t)
      vec.t11 <- unique(vec.t11)
      L.t11 <- length(vec.t11)
      if (vec.t11[L.t11] >= M11)
        # FIXME why not just if(t>M11)
        stop(
          "for the tansition '11' the effects can not be estimated for the given 't'(large 't' returns all responses equal to 0) "
        )

      # Compute point estimate
      eta.list <-
        future.apply::future_lapply(vec.t11, function(x)
          fit_single_time_point_estimate(
            s,
            t = x,
            transition = "11",
            X = X,
            data_df = data1,
            ratetable = ratetable,
            rmapsub = rmapsub
          ))
      # Compute bootstrap
      boot.eta <-
        lapply(vec.t11, function(x){ # use single worker apply here as multiworker is used for bootstrapping
          compute_single_time_bootsraps(
            s,
            t = x,
            transition = "11",
            X = X,
            data_df = data1,
            ratetable = ratetable,
            rmapsub = rmapsub
          )})

      eta_list <- future.apply::future_lapply(eta.list, function(x)
        x[["eta"]])
      sd_list <- lapply(eta.list, function(x)
        x[["sd"]])
      coef <- do.call("bind_rows", eta_list)
      sd <- do.call("bind_rows", sd_list)

      CO <-
        list(
          transition = "11",
          formula = formula,
          time = vec.t11,
          coefficients = coef,
          SD = sd,
          LWL = coef - 1.96 * sd,
          UPL = coef + 1.96 * sd,
          p.value = 2 * pnorm(-abs(as.matrix(coef / sd)))
        )
      if (trans == "all") {
        co$co11 = CO
      }
      else {
        co <-
          list(
            "co" = CO,
            call = match.call(),
            formula = formula,
            transition = trans,
            s = s,
            t = t,
            n.misobs = n.misobs
          )
        class(co) = "TPreg"
        return(co)
      }
    }


   # Look at the 1->2 (illness) transition
    if(trans == "12" || trans == "all"){
      index <- data1$Zt < data1$Tt
      vec.t12 <- c(data1$Zt[index], data1$Tt[index])
      # QUESTION why are we using both Zt and Tt here? Only Zt accounts for individual at risk of transition 1->2
      M12 <- max(vec.t12)
      vec.t12 <- vec.t12[order(vec.t12[vec.t12 > s])] #FIXME >s already guaranteed by data1 filtering
      vec.t12 <- vec.t12[vec.t12 <= t]
      vec.t12 <- vec.t12[seq(1, length(vec.t12), by)]
      vec.t12 <- c(vec.t12, t)
      vec.t12 <- unique(vec.t12)
      L.t12 <- length(vec.t12)
      if(vec.t12[L.t12] >= M12) stop(" for the tansition '12' the effects can not be estimated for the given 't', (large 't' returns all responses equal to 0)")
      eta.list <- lapply( vec.t12, function(x){
        jumptime <- x
        res <- (data1$Zt <= jumptime & jumptime < data1$Tt)
        delta_t <- ifelse(data1$Tt <= jumptime , data1$delta, 1)
        hatG <- sapply(pmin(data1$Tt, jumptime), Shat.function)
        wei <- delta_t/hatG


        t=x
        vv <- rellogit(t,data1)
        family <- binomial(link = vv)
        eta<-mod.glm.fit.callingwrapper( X, res, family = family, weights = wei,
                                       warning_str = paste0(" for transition 1->2, s=",s," t=", jumptime), maxmaxit = 1000)
        data1 <- data1
        formula1 <- formula
        X <- X

	# Bootstrap
        r <-
          foreach(
            j = 1:R,
            .combine = rbind,
            .errorhandling = "stop",
            .inorder = FALSE
          ) %dofuture% {
            iboot <- sample(1:nrow(data1), replace=TRUE)
          boot.data <- data1[iboot, ]
          vv <- rellogit(t,boot.data)
          family <- binomial(link = vv)

          return(mod.glm.fit.callingwrapper(X[iboot, ,drop=FALSE],res[iboot],family=family,weights=wei[iboot],
                                          warning_str = paste0(" on bootstrap sample ",j," for transition 1->2, s=",s," t=", jumptime),
                                          maxmaxit = 1000))

        }
        boot.eta <- r
        boot.sd <- apply(boot.eta, 2, sd, na.rm = FALSE)
        return(list(eta=bind_rows(eta), sd=bind_rows(boot.sd)))
      })
      eta_list<-lapply(eta.list,function(x) x[["eta"]])
      sd_list<-lapply(eta.list,function(x) x[["sd"]])
      coef <- do.call("bind_rows", eta_list)
      sd <- do.call("bind_rows", sd_list)

      CO = list( transition = "12",formula=formula, time = vec.t12, coefficients = coef, SD = sd, LWL = coef - 1.96*sd,UPL=coef+1.96*sd, p.value = 2*pnorm(-abs(as.matrix(coef/sd))))
      if(trans == "all"){
        co$co12 = CO
      }
      else {
        co <- list("co" = CO, call = match.call(),formula=formula, transition = trans, s = s, t = t, n.misobs = n.misobs)
        class(co) = "TPreg"
        return(co)
      }
    }

    # Look at the 1->3 (direct death) transition
    if(trans == "13" || trans == "all"){
      vec.t13 <- data1$Tt
      M13 <- max(vec.t13)
      vec.t13 <- vec.t13[order(vec.t13[vec.t13 > s])]
      vec.t13 <- vec.t13[vec.t13 <= t]
      vec.t13 <- vec.t13[seq(1, length(vec.t13), by)]
      vec.t13 <- c(vec.t13, t)
      vec.t13 <- unique(vec.t13)
      L.t13 <- length(vec.t13)
      if(vec.t13[L.t13] >= M13) stop(" for the transition '13' the effects can not be estimated for the given 't', (large 't' returns all responses equal to 1) ")
      eta.list <- lapply( vec.t13, function(x){
        jumptime <- x
        res <- (data1$Tt <= jumptime)
        delta_t <- ifelse(data1$Tt <= jumptime , data1$delta, 1)
        hatG <- sapply(pmin(data1$Tt, jumptime), Shat.function)
        wei <- delta_t/hatG

        t<-x
        vv <- offsetlogit(t,data1)
        family <- binomial(link = vv)


        eta<-mod.glm.fit.callingwrapper( X, res, family = family, weights = wei,
                                      warning_str = paste0(" for transition 1->3, s=",s," t=", jumptime), maxmaxit = 1000)
        data1 <- data1
        formula1 <- formula
        X <- X
	# Bootstrap
        r <-
          foreach(
            j = 1:R,
            .combine = rbind,
            .errorhandling = "stop",
            .inorder = FALSE
          ) %dofuture% {
            iboot <- sample(1:nrow(data1), replace=TRUE)
          boot.data <- data1[iboot, ]
          vv <- offsetlogit(t,boot.data)
          family <- binomial(link = vv)

          return(mod.glm.fit.callingwrapper(X[iboot, ,drop=FALSE], res[iboot], family = family, weights = wei[iboot],
                                          warning_str = paste0(" on bootstrap sample ",j," for transition 1->3, s=",s," t=", jumptime),
                                          maxmaxit = 1000))

        }
        #print(r)
        boot.eta <- r
        boot.sd <- apply(boot.eta, 2, sd, na.rm = FALSE) #this will return NA if any bootstrap sample did not contain any event
        return(list(eta=bind_rows(eta), sd=bind_rows(boot.sd)))
      })
      eta_list<-lapply(eta.list,function(x) x[["eta"]])
      sd_list<-lapply(eta.list,function(x) x[["sd"]])
      coef <- do.call("bind_rows", eta_list)
      sd <- do.call("bind_rows", sd_list)

      CO <- list(transition = "13",formula=formula, time = vec.t13, coefficients = coef, SD = sd, LWL = coef - 1.96*sd, UPL = coef + 1.96*sd, p.value = 2*pnorm(-abs(as.matrix(coef/sd))))
      if(trans == "all"){
        co$co13 = CO
      }
      else {
        co <- list("co" = CO, call = match.call(),formula=formula, transition = trans, s = s, t = t, n.misobs = n.misobs)
        class(co) = "TPreg"
        return(co)
      }
    }


    # Look at the 2->3 (death after illness) transition
    if(trans == "23" || trans == "all"){
      vec.t23 <- data2$Tt
      M23 <- max(vec.t23)
      vec.t23 <- vec.t23[order(vec.t23[vec.t23 > s])]
      vec.t23 <- vec.t23[vec.t23 <= t]
      vec.t23 <- vec.t23[seq(1, length(vec.t23), by)]
      vec.t23 <- c(vec.t23,t)
      vec.t23 <- unique(vec.t23)
      L.t23 <- length(vec.t23)
      if(vec.t23[L.t23] >= M23) stop(" for the tansition '23' the effects can not be estimated for the given 't'(large 't' returns all responses equal to 1)")
      eta.list <- lapply( vec.t23, function(x){
        jumptime <- x
        res <-(data2$Tt <= jumptime)
        delta_t <- ifelse(data2$Tt <= jumptime, data2$delta, 1)
        hatG <- sapply(pmin(data2$Tt, jumptime), Shat.function23)
        wei <- delta_t/hatG

        t<-x
        vv <- offsetlogit(t,data2)
        family <- binomial(link = vv)

        eta<-mod.glm.fit.callingwrapper( X2, res, family = family, weights=wei,
                                       warning_str = paste0(" for transition 2->3, s=",s," t=", jumptime), maxmaxit = 1000)
        data2 <- data2
        formula1 <- formula
        X2 <- X2

	#Bootstrap
        r <-
          foreach(
            j = 1:R,
            .combine = rbind,
            .errorhandling = "stop",
            .inorder = FALSE
          ) %dofuture% {
            iboot <- sample(1:nrow(data2), replace=TRUE)
          boot.data <- data2[iboot, ]
          vv <- offsetlogit(t,boot.data)
          family <- binomial(link = vv)

          return(mod.glm.fit.callingwrapper(X2[iboot, ,drop=FALSE], res[iboot], family = family, weights = wei[iboot],
                                          warning_str = paste0(" on bootstrap sample ",j," for transition 1->3, s=",s," t=", jumptime),
                                          maxmaxit = 1000))
        }
        boot.eta <- r
        boot.sd <- apply(boot.eta,2, sd, na.rm = FALSE)
        return(list(eta=bind_rows(eta), sd=bind_rows(boot.sd)))
      })
      eta_list<-lapply(eta.list,function(x) x[["eta"]])
      sd_list<-lapply(eta.list,function(x) x[["sd"]])
      coef <- do.call("bind_rows", eta_list)
      sd <- do.call("bind_rows", sd_list)

      CO <- list(transition = "23", formula=formula, time = vec.t23, coefficients = coef, SD = sd, LWL = coef - 1.96*sd, UPL = coef + 1.96*sd, p.value = 2*pnorm(-abs(as.matrix(coef/sd))))
      if(trans == "all"){
        co$co23=CO
      }
      else {
        co <- list("co" = CO, call = match.call(),formula=formula, transition = trans, s = s, t = t, n.misobs=n.misobs)
        class(co)="TPreg"
        return(co)
      }
    }

    # ? edit return object properties if "all"??
    if(trans == "all"){
      co$call = match.call()
      co$formula = formula
      co$transition = "all"
      co$s = s
      co$t = t
      co$n.misobs = n.misobs
      class(co) = "TPreg"
      return(co)
    }
  }
}

estimate_censoring_dist <-
  function(s, t, X, data_df, rhs_formula = NULL) {
    # =========================================================================
    # Theory
    # =========================================================================
    # Azarang et al, 2017 define two different censoring distributions G:
    # 1->1, 1->3 and 1->2 transitions:
    # $$ G^{(s)}_x(t) =  P(C>=t|C>=s,X=x) $$
    # 2->2 and 2->3 transition:
    # $$ G^{[s]}_x(t) =  P(C>=t|C>=s,Z<=s<T, X=x) $$
    # However Z and T are assumed conditionnally independent of C given X
    # (non-informative censoring), and the reverse is true. The two definitions
    # above are inconsistent regarding this assumption
    # Thus I'll use $ G^{(s)}_x(t) = G^{[s]}_x(t) =  P(C>=t|C>=s,X=x) $ for all
    # transitions.

    # =========================================================================
    # Implementation
    # =========================================================================
    # Convert to data.table for easier row-based selection
    if (!data.table::is.data.table(data_df)) {
      data_df = data.table::as.data.table(data_df)
    }
    data_sub = data_df[!(Tt < s & delta == 0)]
    # Create survival objects, carefully handle T<s censoring
    # Tt<s & delta==1 (eq. to C<s) have been filtered out above
    surv_resp = survival::Surv(pmax(data_sub$Tt - s, 0), data_sub$delta == 0)
    # TODO allow the use of rhs_formula and other estimators
    cens_fit = survival::survfit(surv_resp ~ +1)
    # Use summary() to keep only times of censoring events
    cens_fit = summary(cens_fit, censored = FALSE)
    # Add start time (t=s) censoring probability
    cens_surv = tibble::tibble(time = cens_fit$time, surv = cens_fit$surv) %>%
      tibble::add_row(time = 0.0, surv = 1.0)
    return(cens_surv)
  }

#' Estimated survival at time t from a survfit like df.
#'
#' @param t Single numeric value
#' @param survfit_data_df Should contain at least columns `surv` and `time`.
#'
#' @return A single numeric value
#'
#'
#' @examples
get_survival_at <- function(t, survfit_data_df) {
  # Return the `surv` value for the row with lowest time greater than t.
  surv_t  = survfit_data_df %>%
    dplyr::filter(time >= t) %>%
    dplyr::slice_min(order_by = time)
  return(surv_t$surv[[1]])
}

fit_single_time_point_estimate <-
  function(s, t, transition, res, X, data_df, ratetable, rmapsub) {
    # Compute censoring weights
    cens_surv = estimate_censoring_dist(s, t, X, data_df)

    if (!data.table::is.data.table(data_df)) {
      data_df = data.table::as.data.table(data_df)
    }

    # Update censoring indicators based on considered time t
    censor_weights = NULL
    censor_indicators = NULL
    shorthand_fun = function(x) get_survival_at(x, cens_surv)
    if (transition == "11") {
      # FIXME implement data filtering based on state at time s
      data_df[Zt > t, delta1 := 1]
      censor_surv_t = sapply(pmin(data_df$Zt, t),
                             shorthand_fun)
      censor_indicators = data_df$delta1
    }
    else {
      # FIXME implement data filtering based on state at time s
      data_df[Tt > t, delta := 1]
      censor_surv_t = sapply(pmin(data_df$Tt, t),
                             shorthand_fun)
      censor_indicators = data_df$delta
    }
    censor_weights = censor_indicators / censor_surv_t

    # Create link function and family objects taking into account background mortality
    custom_link = if (transition %in% c('11', '12', '22')) {
      rellogit(
        s = s,
        t = t,
        data_df = data_df,
        ratetable = ratetable,
        rmapsub = rmapsub
      )
    }
    else {
      offsetlogit(
        s = s,
        t = t,
        data_df = data_df,
        ratetable = ratetable,
        rmapsub = rmapsub
      )
    }
    # TODO check whether the binomial variance should be adjusted too 
    custom_family = binomial(link = custom_link)

    # Extract binary response
    y <- if(transition == "11"){
      # TODO Check that censoring is correctly applied here
      # I feel like because the event is an absence of one censoring time is not
      # correctly applied
      data_df$Zt > t
    }
    else {
      stop("Transitions other than 11 not implemented")
    }

    # Fit the GLM
    eta <-
      mod.glm.fit.callingwrapper(
        X,
        y,
        family = custom_family,
        weights = censor_weights,
        warning_str = paste0(" for transition ", transition, ", s=", s, " t=", t),
        maxmaxit = 1000
      )
    return(eta)
  }

compute_single_time_bootstrap_sample <-
  function(s, t, transition, X, data_df, ratetable, rmapsub) {
    # Sample row ids with replacement
    n = nrow(data_df)
    boot_ids = sample(1:n, n, replace = TRUE)
    return(fit_single_time_point_estimate(s, t, transition, X[boot_ids,], data_df[boot_ids,], ratetable, rmapsub))
  }

compute_single_time_bootsraps <-
  function(n_boot, s, t, transition, X, data_df, ratetable, rmapsub) {
    boot_res <-
      future.apply::future_replicate(n = n_boot, expr = compute_single_time_bootstrap_sample(s, t, transition, X, data_df, ratetable, rmapsub))
    return(boot_res)
  }

summarize_single_time_bootstraps <- function(boot_res) {
  #return(quantile(boot_res, probs = c(0.025, .975)))
  return(apply(boot_res, 2, sd, na.rm = FALSE))
}

rellogit <- function(s, t, data_df, ratetable, rmapsub) {
  SL <-
    eval(substitute(
      compute_survprob_pch(data_df, t - s, ratetable, rmap = rmapsubs),
      list(rmapsubs = rmapsub)
    ))$expsurvs
  linkfun <- function(mu){
    log((mu / SL) / abs(1 - (mu / SL)))
  }
  linkinv <- function(eta){
    SL * exp(eta) / (1 + exp(eta))
  }
  mu.eta <- function(eta) {
    SL * exp(eta) / (1 + exp(eta)) ^ 2
  }
  valideta <- function(eta){
    TRUE
  }
  name <- "Relative Logit"
  return(structure(
    list(
      linkfun = linkfun,
      linkinv = linkinv,
      mu.eta = mu.eta,
      valideta = valideta,
      name = name
    ),
    class = "link-glm"
  ))
}

# an offset survival logit link function for 13 and 23 transitions
offsetlogit <- function(s, t, data_df, ratetable, rmapsub) {
  SL <-
    eval(substitute(
      compute_survprob_pch(data_df, t - s, ratetable, rmap = rmapsubs),
      list(rmapsubs = rmapsub)
    ))$expsurvs
  dp <- 1 - SL #death probability
  linkfun <- function(miu)
    log((miu - dp) / (1 - (miu - dp)))
  linkinv <- function(et)
    (exp(et) * SL + dp) / (1 + exp(et))
  mu.eta <- function(et) {
    ((2 * SL - 1) * exp(et)) / ((1 + exp(et)) ^ 2)
  }
  valideta <- function(et)
    TRUE
  link <- "log((miu-dp)/(1-(miu-dp)))"
  structure(
    list(
      linkfun = linkfun,
      linkinv = linkinv,
      mu.eta = mu.eta,
      valideta = valideta,
      name = link
    ),
    class = "link-glm"
  )
}
