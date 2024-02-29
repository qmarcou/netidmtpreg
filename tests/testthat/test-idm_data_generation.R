devtools::load_all()
test_that("Homogeneous Poisson count process generation works", {
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
