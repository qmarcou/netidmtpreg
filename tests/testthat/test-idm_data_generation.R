devtools::load_all()
testthat::test_that("Homogeneous Poisson count process generation works", {
  n_ind = 1000
  lambda = .25
  n_steps = 10
  synth_data <-
    generate_poisson_time_to_event(
      n_individuals = n_ind,
      lambda = lambda,
      n_timesteps = n_steps,
      recurring = FALSE
    )
  testthat::expect_vector(synth_data, size = n_ind)
  synth_data <-
    tibble::as_tibble_col(synth_data, column_name = "event_time")
  synth_data %>% assertr::assert(assertr::within_bounds(1, n_steps), event_time)

  # Recurring events not implemented
  testthat::expect_error(
    generate_poisson_time_to_event(
      n_individuals = n_ind,
      lambda = lambda,
      n_timesteps = n_steps,
      recurring = TRUE
    ),
    class = "not_implemented_error"
  )
  ##############################################################################
  # Function arguments checks
  ##############################################################################

  # Non negative integer expected for n_individuals
  testthat::expect_no_error(
    generate_poisson_time_to_event(
      n_individuals = 3.0,
      lambda = lambda,
      n_timesteps = n_steps,
      recurring = FALSE
    )
  )
  testthat::expect_error(
    generate_poisson_time_to_event(
      n_individuals = 3.1,
      lambda = lambda,
      n_timesteps = n_steps,
      recurring = FALSE
    ),
    class = "invalid_argument_error"
  )
  testthat::expect_error(
    generate_poisson_time_to_event(
      n_individuals = -n_ind,
      lambda = lambda,
      n_timesteps = n_steps,
      recurring = FALSE
    ),
    class = "invalid_argument_error"
  )
  # Non negative integer expected for n_timesteps
  testthat::expect_no_error(
    generate_poisson_time_to_event(
      n_individuals = n_ind,
      lambda = lambda,
      n_timesteps = 3.0,
      recurring = FALSE
    )
  )
  testthat::expect_error(
    generate_poisson_time_to_event(
      n_individuals = n_ind,
      lambda = lambda,
      n_timesteps = 3.1,
      recurring = FALSE
    ),
    class = "invalid_argument_error"
  )
  testthat::expect_error(
    generate_poisson_time_to_event(
      n_individuals = n_ind,
      lambda = lambda,
      n_timesteps = -n_steps,
      recurring = FALSE
    ),
    class = "invalid_argument_error"
  )
})

testthat::test_that("Test uncensored simplified exponential IDM data generation", {
  n_ind <- 1000
  lambda_illness <- 1
  lambda_death <- .5
  idm_synth_data <- NULL
  testthat::expect_no_error(
    idm_synth_data <- generate_uncensored_ind_exp_idm_data(
      n_individuals = n_ind,
      lambda_illness = lambda_illness,
      lambda_death = lambda_death
    )
  )
  testthat::expect_s3_class(idm_synth_data, "data.frame")
  testthat::expect_length(idm_synth_data, 5) # id, Zt, delta1, Tt, delta
  assertr::verify(idm_synth_data,
                  assertr::has_only_names(c('id', 'Zt', 'delta1', 'Tt', 'delta')))
  # TODO use `is.iddata()` function
  testthat::expect_equal(nrow(idm_synth_data), n_ind)
  assertr::assert(idm_synth_data, assertr::in_set(TRUE, allow.na = FALSE), delta, delta1)
})
