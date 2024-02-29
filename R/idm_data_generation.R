#' Generate synthetic time to event data from a homogeneous Poisson process.
#'
#' Note that this is not a efficient way to model time homogeneous non recurring
#'  events that may be better captured by a continuous time exponential
#'  distribution.
#'
#' @param n_individuals integer, Number of individuals to simulate
#' @param lambda float > 0, Rate of the Poisson point process
#' @param n_timesteps integer, Number of evenly spaced timesteps
#' @param recurring logical, whether event is recurring or fully absorbing.
#'
#' @return If recurring=FALSE a vector of first time to event, if recurring=TRUE
#' a list of event times of length n_individuals.
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
    multi_assert(lambda, "lambda", invalid_argument, c(assert_scalar,assert_positive))

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
