#' Generate synthetic time to event data from a homogeneous Poisson process.
#'
#' Note that this is not a efficient way to model time homogeneous non recurring
#' events that may be better captured by a continuous time exponential
#' distribution from \code{\link{generate_exponential_time_to_event}}.
#'
#' @param n_individuals integer, Number of individuals to simulate
#' @param lambda float > 0, Rate of the Poisson point process
#' @param n_timesteps integer, Number of evenly spaced timesteps
#' @param recurring logical, whether event is recurring or fully absorbing.
#'
#' @return If recurring=FALSE a vector of first time to event, if recurring=TRUE
#'   a list of event times of length n_individuals.
#' @export
#'
#' @examples
generate_poisson_time_to_event <-
  function(n_individuals,
           lambda,
           n_timesteps,
           recurring = FALSE) {
    # Check arguments
    assert_count(n_individuals, "n_individuals", invalid_argument)
    assert_count(n_timesteps, "n_timesteps", invalid_argument)
    multi_assert(lambda,
                 "lambda",
                 invalid_argument,
                 c(assert_scalar, assert_positive))

    # Draw random observations with fixed time steps in a single vector
    counts <-
      stats::rpois(n = n_individuals * n_timesteps, lambda = lambda)

    # Reshape to n*m array
    counts <- array(counts, dim = c(n_individuals, n_timesteps))

    # Identify time to events
    events <- counts > 0
    events_times <-
      apply(X = events, MARGIN = 1, function(x)
        which(x)) # find non zeros
    events_times <- lapply(X = events_times, function(x)
      if (length(x) == 0)
        NA
      else
        x) # assign NA when no event is observed

    if (!recurring) {
      # Get the first time to event
      events_times <- sapply(X = events_times, min)
    }
    else{
      not_implemented("Recurring events handling is not implemented yet.")
      if (any(counts > 1)) {
        message(
          "Several events occured in a single time step and are not identifiyable, consider setting a finer grained time step."
        )
      }
    }
    return(events_times)
  }

#'Generate time to events from an exponential distribution
#'
#'This is a simple wrapper around the stats::rexp function, without added value
#'beyond consistent argument checking and format with other generating functions
#'from the package.
#'@param n_individuals integer, Number of individuals to simulate
#'@param lambda float > 0, Rate of the underlying Poisson point process
#'
#'
#'@return Atomic vector of time to (first) event.
#'@export
#'
#' @examples
generate_exponential_time_to_event <-
  function(n_individuals,
           lambda) {
    # TODO generalize to allow the use of a formula for different rates in the
    # pop (first step towards population mortality modeling)
    # Check arguments
    assert_count(n_individuals, "n_individuals", invalid_argument)
    multi_assert(lambda,
                 "lambda",
                 invalid_argument,
                 c(assert_scalar, assert_positive))
    # Generate and return random exponential time to events
    return(rexp(n = n_individuals, rate = lambda))
  }

#' Generate uncensored independent exponential illness-death data.
#'
#' Mostly for testing purposes. Illness and death times are assumed independent,
#' in particular the death rate is assumed equal in healthy and illness states.
#' This simplification creates data that are not really following an
#' Illness-Death model but a simpler semi competitive model, where death
#' prevents observation of illness times.
#'
#' @param n_individuals integer, Number of individuals to simulate
#' @param lambda_illness float > 0, Rate of the underlying Poisson point process
#'   to transition from healthy to illness
#' @param lambda_death float > 0, Rate of death, either from healthy or illness
#'   state.
#'
#' @return a tibble of iddata
#' @export
#'
#' @examples
generate_uncensored_ind_exp_idm_data <-
  function(n_individuals,
           lambda_illness,
           lambda_death) {

    # Defer parameter checking to inner generating functions
    illness_times <-
      generate_exponential_time_to_event(n_individuals = n_individuals, lambda = lambda_illness)
    death_times <-
      generate_exponential_time_to_event(n_individuals = n_individuals, lambda = lambda_death)
    healthy_sojourn_time <- pmin(illness_times, death_times)

    # FIXME Use the more generic `iddata` preparation function
    id_data <-
      tibble::tibble(
        id = seq.int(from = 1, to = n_individuals, by = 1),
        Zt = healthy_sojourn_time,
        delta1 = 1,
        # Uncensored
        Tt = death_times,
        delta = 1 # Uncensored
      )
    return(id_data)
  }

#' Apply censoring times to an existing Illness-Death dataframe.
#'
#' Update event times (Zt and Tt) and censoring indicators according to the
#' prescribed censoring time(s). Note that previous censoring information will
#' be updated where newer censoring times are provided.
#'
#' @param iddata_df An Illness-Death dataframe in `iddata` format.
#' @param censoring_times Censoring time(s), either a single value or a vector
#'  length `nrow(iddata_df)`
#'
#' @return An iddata tibble with updated censoring times
apply_iddata_censoring <- function(
    iddata_df,
    censoring_times) {
  # TODO Check that censoring time and iddata_df have same length

  # Apply censoring
  iddata_df <- data.table::as.data.table(iddata_df) %>%
    tibble::add_column(SQVVcCs1lD4R7tDVlOoV_cens_times = censoring_times)
  iddata_df[
    SQVVcCs1lD4R7tDVlOoV_cens_times < Zt,
    ":="(delta1 = 0, Zt = SQVVcCs1lD4R7tDVlOoV_cens_times)
  ]
  iddata_df[
    SQVVcCs1lD4R7tDVlOoV_cens_times < Tt,
    ":="(delta = 0, Tt = SQVVcCs1lD4R7tDVlOoV_cens_times)
  ]
  iddata_df <- tibble::as_tibble(iddata_df) %>%
    dplyr::mutate(SQVVcCs1lD4R7tDVlOoV_cens_times = NULL)
  return(iddata_df)
}

#' Update death times in an existing Illness-Death dataframe.
#'
#' Update event times (Zt and Tt) according to the prescribed death time(s).
#' Death times are only updated when the new death time is earlier than the
#' previous one for a given individual. This function cannot postpone an
#' existing death time as the parameters provide no way to know whether illness
#' occurred between old and new death time.
#'
#' @param iddata_df An Illness-Death dataframe in `iddata` format.
#' @param death_times Updated death time(s), either a single value or a vector
#'  length `nrow(iddata_df)`
#'
#' @return An iddata tibble with updated death times
apply_iddata_death <- function(
    iddata_df,
    death_times) {
  # Apply provided death times
  iddata_df <- data.table::as.data.table(iddata_df) %>%
    tibble::add_column(SQVVcCs1lD4R7tDVlOoV_death_times = death_times)
  iddata_df[
    SQVVcCs1lD4R7tDVlOoV_death_times < Zt,
    Zt := SQVVcCs1lD4R7tDVlOoV_death_times
  ]
  iddata_df[
    SQVVcCs1lD4R7tDVlOoV_death_times < Tt,
    Tt := SQVVcCs1lD4R7tDVlOoV_death_times
  ]
  iddata_df <- tibble::as_tibble(iddata_df) %>%
    dplyr::mutate(SQVVcCs1lD4R7tDVlOoV_death_times = NULL)
  return(iddata_df)
}