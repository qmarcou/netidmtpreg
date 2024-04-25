# Code adapted from idmTPreg package v1.1 under GNU-GPLv2
# Original Authors: Leyla Azarang and Manuel Oviedo de la Fuente
# Adapted for net survival setting by Quentin Marcou, based on <ref>


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
mod.glm.fit.errorwrapper <-
  function(X,
           response,
           family,
           weights,
           maxit = glm.control()$maxit,
           maxmaxit = 1000,
           warning_str = "",
           ...) {
    result <- tryCatchLog::tryCatchLog({
      # Try
      withCallingHandlers({
        mod.glm.fit2(
          X,
          response,
          family = family,
          weights = weights,
          start = rep.int(0, times = ncol(X)),
          control = glm.control(maxit = maxit)
        )
      },
      warning = function(warn) {
        if (stringr::str_detect(warn$message, "no observations informative at iteration")) {
          stop(paste0(
            "Warning caught ",
            warning_str,
            ": ",
            warn$message,
            " returning NA"
          ))
        } else {
          warning(paste0(warning_str, warn$message))
        }
      })
    },
    error = function(err) {
      # Catch
      if (err$message == "inner loop 1; cannot correct step size" ||
          err$message == "inner loop 2; cannot correct step size") {
        if (maxit * 10 <= maxmaxit) {
          # recursively call the wrapper with a greater maxit
          return(
            mod.glm.fit.errorwrapper(
              X = X,
              response = response,
              family = family,
              weights = weights,
              maxit = maxit * 10,
              maxmaxit = maxmaxit,
              warning_str = warning_str
            )
          )
        } else {
          warning(
            paste0(
              "Step size correction issue with maxmaxit reached: ",
              warning_str,
              "returning NA as result"
            )
          )
          tmp <- rep.int(NA, times = dim(X)[[2L]])
          names(tmp) <- dimnames(X)[[2L]]
          return(list(coefficients = tmp, converged = FALSE))
        }
      } else if (stringr::str_detect(err$message, "no observations informative at iteration")) {
        tmp <- rep.int(NA, times = dim(X)[[2L]])
        names(tmp) <- dimnames(X)[[2L]]
        return(list(coefficients = tmp, converged = FALSE))
      } else if (stringr::str_detect(err$message, "NA/NaN/Inf in 'y'")) {
        # Handle an error that i do not fully understand, this seems to be linked to the divergence of the algorithm
        # with derivatives regarding parameters exploding. I do not know whether this could be linked to Layla's
        # modification to the glm.fit code
        warning(
          "mod.glm.fit crashed with error \"",
          err$message,
          "\"",
          warning_str,
          ", returning NA"
        )
        tmp <- rep.int(NA, times = dim(X)[[2L]])
        names(tmp) <- dimnames(X)[[2L]]
        return(list(coefficients = tmp, converged = FALSE))
      } else {
        message(paste0("Exception caught upon calling modl.glm.fit", warning_str))
        stop(err$message)
      }
    })
    return(result)
  }


#' @title Check input and output of mod.glm.fit
#' @description
#' Avoids calling the glm fitting function when there are no contrasted response
#' values (all 0 or all 1). Ensures return of NA coefficient when the algorithm
#' did not converge.
#'
#'
#' @param X array, must have at least 2 dimensions
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
mod.glm.fit.callingwrapper <-
  function(X,
           response,
           family,
           weights,
           maxit = glm.control()$maxit,
           maxmaxit = 1000,
           warning_str = "",
           ...) {
    if (any(response) & !all(response)) {
      # There must be at least one event in the sample in order to learn smthg
      result <- mod.glm.fit.errorwrapper(
        X = X,
        response = response,
        family = family,
        weights = weights,
        maxit = maxit,
        maxmaxit = maxmaxit,
        warning_str = warning_str
      )
      # print(result$converged)
      if (!result$converged || result$boundary) {
        result$coefficients <-
          result$coefficients * NA # set coef to NA if the algorithm did not converge
      }
      return(coefficients(result))
    } else {
      # if no event, coefficients are meaningless and one should return NA
      warning("All provided responses are equal",
              warning_str,
              ", cannot fit GLM, returning NA")
      tmp <- rep.int(NA, times = dim(X)[[2L]])
      names(tmp) <- dimnames(X)[[2L]]
      return(tmp)
    }
  }




#' Fit (net) survival Illness-Death transition probabilities
#'
#' Core function of the package, implementing a
#'
#' @param s DESCRIPTION.
#' @param t Times at which survival should be estimated, either `NULL`, scalar
#' or vector of numerics. Default to NULL. If a scalar value is given, it will
#' be interpreted as the maximum time at which survival shall be estimated
#' using a non-parametric approach. If a vector or list-like object, the values
#' will be interpreted as the list of times at which one whishes to get point
#' estimations of survivals. By default estimation will be made at the time of
#' last events before the requested times. This behavior can be altered by
#' setting `readjust_t` to `FALSE`, mostly for testing purposes. If set to
#' `NULL` (default) `t` is set as a scalar value to the last observed
#' transition time.
#' @param trans DESCRIPTION.
#' @param formula DESCRIPTION.
#' @param ratetable DESCRIPTION.
#' @param time_dep_popvars DESCRIPTION.
#' @param rmap DESCRIPTION.
#' @param data DESCRIPTION.
#' @param link DESCRIPTION.
#' @param R DESCRIPTION.
#' @param by DESCRIPTION.
#' @param readjust_t boolean, default `TRUE`. Whether to use exact `t` values
#' (`FALSE`) or adjust t values to observed transition times (`TRUE`).  
#'
#' @return RETURN_DESCRIPTION
#'
#' @references
#' Azarang, L., Scheike, T., & de Uña‐Álvarez, J. (2017). Direct modeling of
#' regression effects for transition probabilities in the progressive
#' illness–death model. **Statistics in Medicine**, 36(12), 1964-1976.
#'
#' Azarang L, Giorgi R. (2021). Estimation of covariate effects on net
#' survivals in the relative survival progressive illness-death model.
#' **Statistical Methods in Medical Research**, 30(6), 1538-1553.
#' doi:10.1177/09622802211003608
#' @examples
#' # ADD_EXAMPLES_HERE
renewnetTPreg <- function(s = 0,
                          t = NULL,
                          trans,
                          formula,
                          ratetable,
                          time_dep_popvars = NULL,
                          rmap = NULL,
                          data,
                          link = "logit",
                          R = 199,
                          by = NULL,
                          readjust_t = TRUE) {
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

    if (missing(data)) {
      stop("Argument 'data' is missing with no default")
    }
    if (!is.data.frame(data)) {
      stop("Argument 'data' must be a data.frame")
    }
    if (sum(!c("id", "Zt", "Tt", "delta1", "delta", "age") %in% (colnames(data))) > 0) { # at least age should be included
      stop("data should  contain  id, Zt, Tt, delta1, delta, age variables")
    } # TODO later on add also "sex" covariable
    if (ncol(data) <= 4) {
      stop("'data' must have covariables")
    }
    if (sum(is.na(data)) != 0) {
      miscolumn <- sapply(data, function(x) sum(is.na(x)))
      miscolnam <- names(miscolumn[miscolumn != 0])
      if (sum(miscolumn) != 0) {
        warning(sapply(miscolnam, function(x) paste(x, " variable in 'data' has missing value(s)", ", ", sep = "")))
      }
    }


    # This piece of code seems to be used to make sure the formula and data columns match
    # It then builds a model.matrix object (which allows automatic expansion in dummy variables?)
    # Still I do not understand all this fiddling with getting the function call (espcially since there's no ... argument in the function)
    # why not use directly stats::model.frame(formula,data) ?
    formula <- formula # wtf?
    cl <- match.call() # get the function call (or ""command line"") string with named arguments
    mf <- match.call(expand.dots = FALSE) # same without expending ... argument, the idea is to reuse the function call, and substitute the function by model.frame
    m <- match(c("formula", "data"), names(mf), 0L) # look for "formula" and "data" in arguments outside '...' and return their indices, return 0L (=0) if not found
    mf <- mf[c(1L, m)]
    # add options for model.frame()
    mf$drop.unused.levels <- TRUE # simplify factors and retain only used levels
    mf$na.action <- na.pass # all missing values will be retained in the model frame
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    if (is.empty.model(mt)) {
      stop("'formula' must match with 'data'")
    } else { # FIXME useless else since if is a stop
      X <- model.matrix(mt, mf, contrasts)
      ind <- match(colnames(X), colnames(data))
      ind <- ind[!is.na(ind)]
      covname <- colnames(data)[ind]
    }

    # Filter useful columns based rmap and formula arguments
    formnames <- all.vars(formula, functions = FALSE, unique = TRUE)
    rmapnames <- all.vars(substitute(rmap), functions = FALSE, unique = TRUE)
    ordata <- data[, unique(c(c("id", "Zt", "Tt", "delta1", "delta"), rmapnames, formnames), fromLast = FALSE)]

    L.or <- nrow(ordata)
    # Remove lines with at least one missing value
    comdata <- ordata[complete.cases(ordata), ]
    X <- X[complete.cases(ordata), , drop = FALSE]

    L.com <- nrow(comdata)
    n.misobs <- L.or - L.com
    if (is.null(by)) {
      by <- floor((max(comdata$Zt) - min(comdata$Zt)) / quantile(comdata$Zt, 0.01))
    }

    if (is.null(t)) {
      t <- max(comdata$Zt[comdata$delta1 == 1], na.rm = T)
    }

    if (any(t <= s) || s < 0) {
      stop("argument 's' must be smaller than 't' and larger than 0")
    }

    # FIXME Why is this all contained in an "else" statement?
    if (!(link %in% c("logit", "probit", "cauchit"))) {
      stop(paste("binomial family do not have", "'", link, "'", "link"))
    }

    if (sum(comdata$delta1 < comdata$delta) != 0) {
      stop("'delta' must be 0 when 'delta1' is 0")
    }
    if (!(trans %in% c("11", "12", "13", "23", "all"))) {
      stop(paste(trans, "is not a valid transition for a progressive illness-death model"))
    }

    if (s == 0 & (trans == "23" || trans == "all")) { # TODO there's probably better way around
      stop("for the transition '23' argument 's' must be larger than 0")
    }

    co <- vector("list", 4)
    names(co) <- c("co11", "co12", "co13", "co23")

    # Check correctness of time_dep_popvars and its interplay with rmap
    rmapsub <- substitute(rmap)
    if (is.list(time_dep_popvars) || is.character(time_dep_popvars)) {
      # FIXME add check for list content being strings
      if (!all(time_dep_popvars %in% names(rmapsub)[-1])) {
        stop("Names in `time_dep_popvars` do not correspond to ratetable dimension names in the `rmap` argument.")
      }

      # Now change rcall to account for s days time shift in survival computation
      for (var in time_dep_popvars) {
        rmapsub[[var]] <- call("+", rmapsub[[var]], s)
      }
    } else if (!is.null(time_dep_popvars)) {
      stop("`time_dep_popvars` must be a list or vector of strings or NULL")
    }

    if (trans == "all") {
      transitions <- c("11", "12", "13", "23")
    } else {
      transitions <- c(trans)
    }

    for (trans in transitions) {
      # Select individuals at risk at time s
      if (trans %in% c("11", "12", "13")) {
        X_sub <- X[comdata$Zt > s, , drop = FALSE]
        data_sub <- comdata[comdata$Zt > s, ]
      } else {
        # patients in state 2 at time s
        X_sub <- X[comdata$Zt <= s & comdata$Tt > s, , drop = FALSE]
        data_sub <- comdata[comdata$Zt <= s & comdata$Tt > s, ]
        vec.t <- data_sub$Zt
      }

        # Define time steps for non parametric estimation
        if (trans == "11") {
          vec.t <- data_sub$Zt
        } else if (trans == "12") {
          index <- data_sub$Zt < data_sub$Tt
          vec.t <- c(data_sub$Zt[index], data_sub$Tt[index])
        } else {
          vec.t <- data_sub$Tt
        }

        # Assert whether meaningful results can be obtained at requested times
        if (max(t) > max(vec.t)) {
          # FIXME handle censoring times in Zt
          stop(
            glue::glue("Effects cannot be estimated for the transition '{trans}'
            for the given 't={t}'(large 't' returns all responses
            equal to 0). Last observed transition of interest at time
            {max(vec.t)}.")
          )
        }

      if (assertthat::is.scalar(t)) {
        # Order and subsample estimation times for the non-parametric setting
        vec.t <- vec.t[vec.t < t]
        vec.t <- vec.t %>%
          unique() %>%
          sort()
        vec.t <- vec.t[seq(1, length(vec.t), by)]
        vec.t <- c(vec.t, t)
      } else {
        if (readjust_t) {
          not_implemented("vector-like t with readjust_t=TRUE is not yet
             implemented")
          # 
          vec.t <- NULL
        } else {
          vec.t <- as.numeric(t)
        }
      }

      # Compute point estimate
      print("estimate")
      coef_estimates <-
        eval(substitute( # Substitute rmap altered with s time shift
          future.apply::future_lapply(vec.t, function(x) {
            fit_single_time_point_estimate(
              s,
              t = x,
              transition = trans,
              X = X_sub,
              data_df = data_sub,
              ratetable = ratetable,
              rmap = rmapsub
            )
          }),
          list(rmapsub = substitute(rmapsub))
        ))
      # unlist the results into a single tibble instead of a list, add t
      # information as first column
      coef_estimates <-
        coef_estimates %>%
        dplyr::bind_rows() %>%
        tibble::add_column(t = vec.t, .before = 1)

      # Compute bootstrap
      print("bootstrap")
      boot_summaries <-
        eval(substitute( # Substitute rmap altered with s time shift
          # use single worker apply here as multiworker is used for bootstrapping
          lapply(vec.t, function(x) {
            compute_single_time_bootstraps( # Compute bootstraps estimates
              n_boot = R,
              s = s,
              t = x,
              transition = trans,
              X = X_sub,
              data_df = data_sub,
              ratetable = ratetable,
              rmap = rmapsub
            ) %>% # and summarize them (estimate SD, CIs)
              summarize_single_time_bootstraps()
          }),
          list(rmapsub = substitute(rmapsub))
        ))
      # unlist the results into a single tibble instead of a list, add t
      # information as first column
      boot_summaries <-
        boot_summaries %>%
        dplyr::bind_rows() %>%
        tibble::add_column(t = vec.t, .before = 1)


      # Combine point estimate and bootstrap information into a single tibble
      ## Format tibble to get one line per time point and parameter
      boot_summaries <-
        boot_summaries %>%
        tidyr::pivot_longer(
          cols = !t,
          names_sep = "_",
          names_to = c("covariate", "statistic")
        ) %>%
        tidyr::pivot_wider(names_from = "statistic", values_from = "value")
      ## Same thing for the point estimates tibble
      coef_estimates <-
        coef_estimates %>% tidyr::pivot_longer(
          cols = !t,
          names_to = "covariate",
          values_to = "estimate"
        )
      ## Join the resulting tables on time step and covariate name
      results_df <- dplyr::full_join(
        x = coef_estimates, y = boot_summaries,
        by = c("t", "covariate")
      )
      ## FIXME Retransform the result to existing format for temporary
      ## consistency with existing methods
      extract_stat_matrix <- function(stat_name) {
        ## define a closure for conciseness
        return(
          results_df %>%
            dplyr::select(t, covariate, {{ stat_name }}) %>%
            tidyr::pivot_wider(
              names_from = "covariate",
              values_from = {{ stat_name }}
            ) %>%
            dplyr::select(!t) %>%
            as.matrix()
        )
      }

      # Create returned object
      CO <-
        list(
          transition = trans,
          formula = formula,
          time = vec.t,
          coefficients = extract_stat_matrix("estimate"),
          SD = extract_stat_matrix("sd"),
          LWL = extract_stat_matrix("ci.lb"),
          UPL = extract_stat_matrix("ci.ub"),
          n.failed.boot = NULL
        )
      if (trans == "all") {
        co$co11 <- CO
      } else {
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
        class(co) <- "TPreg"
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
      data_df <- data.table::as.data.table(data_df)
    }
    data_sub <- data_df[!(Tt < s & delta == 0)]
    # Create survival objects, carefully handle T<s censoring
    # Tt<s & delta==1 (eq. to C<s) have been filtered out above
    surv_resp <- survival::Surv(pmax(data_sub$Tt - s, 0), data_sub$delta == 0)
    # TODO allow the use of rhs_formula and other estimators
    cens_fit <- survival::survfit(surv_resp ~ +1)
    # Use summary() to keep only times of censoring events
    cens_fit <- summary(cens_fit, censored = FALSE)
    # Add start time (time = s) censoring probability as first row
    cens_surv <- tibble::tibble(time = cens_fit$time, surv = cens_fit$surv) %>%
      tibble::add_row(time = 0.0, surv = 1.0, .before = 1)

    return(cens_surv)
  }

#' Estimated survival at time t from a survfit like df.
#'
#' @param t Single or vector numeric value $\geq$ 0.
#' @param survfit_data A survfit, summary.survfit or data.frame object. If
#'   data.frame, it should contain at least columns `surv` and `time`.
#' @param safe boolean, Indicate whether to enable input argument checking,
#'   default to TRUE.
#'
#' @return A survival probabilities of same length as `t`.
#'
#'
#' @examples
get_survival_at <- function(t, survfit_data, safe = TRUE) {
  # Return the `surv` value for the row with greatest time lower than t using a
  # step function
  if (is(survfit_data, "survfit")) survfit_data <- summary(survfit_data)
  if (is(survfit_data, "summary.survfit")) {
    survfit_data <- tibble::tibble(
      time = survfit_data$time,
      surv = survfit_data$surv
    ) %>%
      tibble::add_row(time = 0.0, surv = 1.0, .before = 1)
  }
  # Sanity checks
  if (safe) {
    assert_positive(t, "t", invalid_argument)
    assertr::verify(
      survfit_data,
      assertr::has_all_names("time", "surv"),
      error_fun = invalid_argument(
        "survfit_data", "contain both a `time` and `surv` column"
      )
    )
    assert_probability(survfit_data$surv, "survfit_data$surv", invalid_argument)
  }

  # Handle the no event case
  if (nrow(survfit_data) == 1) {
    if (all(t >= survfit_data$time[[1]])) {
      return(rep_len(1.0, length(t)))
    } else {
      invalid_argument(
        "survfit_data",
        "must contain several breakpoints or all `t` values must be greater
        than the provided time breakpoint"
      )
    }
  }
  # Normal case: use a stepfunc for efficient vectorised lookup
  surv_t <- stats::stepfun(x = survfit_data$time[-1], y = survfit_data$surv)

  return(surv_t(t))
}

fit_single_time_point_estimate <-
  function(s, t, transition, X, data_df, ratetable, rmap) {
    # Convert data_df to data.table for efficiency
    if (!data.table::is.data.table(data_df)) {
      data_df <- data.table::as.data.table(data_df)
    }

    # Compute censoring weights
    cens_surv <- estimate_censoring_dist(s, t, X, data_df)

    # Update censoring indicators based on considered time t
    censor_weights <- NULL
    censor_indicators <- NULL
    shorthand_fun <- function(x) get_survival_at(x, cens_surv)
    if (transition == "11") {
      # FIXME implement data filtering based on state at time s
      # Construct \Delta^{1}_{t} = 1_{min(Z,t) \leq C} in Azarang 2017
      # data_df==delta1 already implies Z \leq C, update indicator based on t
      data_df[t <= Zt & delta1 == 0, delta1 := 1]
      censor_surv_t <- shorthand_fun(pmin(data_df$Zt, t))
      censor_indicators <- data_df$delta1
    }
    else {
      # FIXME implement data filtering based on state at time s
      # \Delta_{t} = 1_{min(T,t) \leq C} in Azarang 2017, same reasoning as
      # before
      data_df[t <= Tt & delta == 0, delta := 1]
      censor_surv_t <- shorthand_fun(pmin(data_df$Tt, t))
      censor_indicators <- data_df$delta
    }
    censor_weights <- censor_indicators / censor_surv_t

    # Create link function and family objects taking into account background
    # mortality
    custom_link_constructor <- if (transition %in% c("11", "12", "22")) {
      rellogit
    } else {
      offsetlogit
    }

    custom_link <- eval(substitute(
      expr = custom_link_constructor(
        s = s,
        t = t,
        data_df = data_df,
        ratetable = ratetable,
        rmap = rmapsub # Protect rmap non standard evaluation for survexp
      ),
      env = list(rmapsub = substitute(rmap))
    ))

    # TODO check whether the binomial variance should be adjusted too
    custom_family <- binomial(link = custom_link)

    # Extract binary response
    y <- if (transition == "11") {
      data_df$Zt > t
    } else if (transition %in% c("13", "23")) {
      data_df$Tt <= t
    } else if (transition %in% c("12", "22")) {
      data_df$Zt <= t & data_df$Tt > t
      # For "22" data_df$Zt <= t is already granted by the fact that s<t
      # This is not computationnally efficient but reduces code duplication
    } else {
      stop(glue::glue("Unknown transition '{transition}'"))
    }

    # Fit the GLM
    eta <-
      mod.glm.fit.callingwrapper(
        X = X,
        response = y,
        family = custom_family,
        weights = censor_weights,
        warning_str = paste0(" for transition ", transition, ", s=", s, " t=", t),
        maxmaxit = 1000
      )
    return(eta)
  }

compute_single_time_bootstrap_sample <-
  function(s, t, transition, X, data_df, ratetable, rmap) {
    # Sample row ids with replacement
    n <- nrow(data_df)
    boot_ids <- sample(1:n, n, replace = TRUE)
    return(
      eval(substitute( # Protect rmap non standard evaluation for survexp
        fit_single_time_point_estimate(
          s = s,
          t = t,
          transition = transition,
          # BUGFIX keep X dimensionality in case there a single term in the
          # formula
          X = X[boot_ids, , drop = FALSE],
          data_df = data_df[boot_ids, ],
          ratetable = ratetable,
          rmap = rmapsub
        ),
        list(rmapsub = substitute(rmap))
      ))
    )
  }

compute_single_time_bootstraps <-
  function(n_boot, s, t, transition, X, data_df, ratetable, rmap) {
    boot_res <-
      eval(substitute( # Protect rmap non standard evaluation for survexp
        future.apply::future_replicate(
          n = n_boot,
          expr = compute_single_time_bootstrap_sample(
            s = s,
            t = t,
            transition = transition,
            X = X,
            data_df = data_df,
            ratetable = ratetable,
            rmap = rmapsub
          ),
          simplify = "array"
        ),
        list(rmapsub = substitute(rmap))
      ))

    # Return a tibble where each col is a coefficient and each row a bootstrap
    if(dim(X)[[2]]==1){
      res <- tibble::as_tibble_col(boot_res,column_name = names(boot_res)[[1]])
    }
    else{
      res <- tibble::as_tibble(t(boot_res))
    }
    return(res)
  }

summarize_single_time_bootstraps <- function(boot_res_df) {
  boot_summary <-
    boot_res_df %>% dplyr::summarise(dplyr::across(
      .cols = everything(),
      .fns = list(
        sd = ~ sd(.x, na.rm = TRUE),
        ci.lb = ~ quantile(.x, prob = 0.025, na.rm = TRUE),
        ci.ub = ~ quantile(.x, prob = .975, na.rm = TRUE),
        n.failed.boot = ~ rlang::are_na(.x) %>% sum()
      )
    ))
  return(boot_summary)
}

rellogit <- function(s, t, data_df, ratetable, rmap) {
  pop_surv <- 1.0 # Allow for crude mortality model fitting
  if(!is.null(ratetable)){
      pop_surv <-
    eval(substitute(
      compute_survprob_pch(data_df, t - s, ratetable, rmap = rmapsub),
      list(rmapsub = substitute(rmap))
    ))$expsurvs
  }

  # Response (mu) to linear predictor (eta) space
  linkfun <- function(mu) {
    log((mu / pop_surv) / abs(1 - (mu / pop_surv)))
  }

  # Linear predictor (eta) to response (mu) space
  linkinv <- function(eta) {
    pop_surv * exp(eta) / (1 + exp(eta))
  }

  # Derivative of the inverse-link function with respect to the linear
  # predictor (eta).
  mu.eta <- function(eta) {
    pop_surv * exp(eta) / (1 + exp(eta)) ^ 2
  }

  # Is the linear predictor eta within the domain of linkinv.
  valideta <- function(eta) {
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


#' An offset survival logit link function closure for 13 and 23 transitions.
#'
#' Returns a offset logit `link-glm` object taking into account population
#' mortality a the considered times s and t.
#'
#' @param s positive scalar.
#' @param t positive scalar.
#' @param data_df A `data.frame` containing population characteristics.
#' @param ratetable a `ratetable` object.
#' @param rmap see survival::survexp.fit.
#'
#' @return a `link-glm` object
#' @examples
#' # ADD_EXAMPLES_HERE
offsetlogit <- function(s, t, data_df, ratetable, rmap) {

  pop_surv <- 1.0 # Allow for crude mortality model fitting
  if (!is.null(ratetable)) {
    pop_surv <-
      eval(substitute(
        compute_survprob_pch(data_df, t - s, ratetable, rmap = rmapsub),
        list(rmapsub = substitute(rmap))
      ))$expsurvs
  }

  pop_death <- 1 - pop_surv #death probability
  # Response (mu) to linear predictor (eta) space
  linkfun <- function(mu) {
    return(log((mu - pop_death) / (1 - (mu - pop_death))))
  }

  # Linear predictor (eta) to response (mu) space
  linkinv <- function(eta) {
    return((exp(eta) * pop_surv + pop_death) / (1 + exp(eta)))
  }

  # Derivative of the inverse-link function with respect to the linear
  # predictor (eta).
  mu.eta <- function(eta) {
    return(((2 * pop_surv - 1) * exp(eta)) / ((1 + exp(eta)) ^ 2))
  }

  # Is the linear predictor eta within the domain of linkinv.
  valideta <- function(eta) {
    # assert eta are finite real numbers
    return(rlang::is_double(eta, finite = TRUE))
  }

  link <- "log((mu-pop_death)/(1-(mu-pop_death)))"
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
