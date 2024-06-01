test_that("Plotting TPreg objects does not throw any error", {
  # FIXME use package's data to reduce code duplication
  # Crude mortality binomial regression
  set.seed(246852389)
  n_ind <- 1e2
  l_illness <- 1.0
  l_death <- 0.1
  s_time <- 0.5
  synth_idm_data <- generate_uncensored_ind_exp_idm_data(
    n_individuals = n_ind,
    lambda_illness = l_illness,
    lambda_death = l_death
  )
  # FIXME I should not have to have those
  # Generate random age and sex labels
  synth_idm_data <-
    synth_idm_data %>% tibble::add_column(
      sex = ifelse(rbinom(n_ind, 1, prob = .5), "male", "female"),
      age = runif(n = n_ind, min = 50, max = 80)
    )
  # for (transition in c("11", "all")) {
  for (transition in c("11")) {
    estimate <-
      fit_netTPreg(~1, synth_idm_data,
        ratetable = NULL,
        rmap = NULL,
        time_dep_popvars = NULL,
        s = s_time,
        t = 3.0,
        by = n_ind / 10,
        trans = transition,
        link = "logit",
        R = 10 # Number of bootstraps
      )
    testthat::expect_no_error(
      plot(estimate)
    )

    # Check with more than just an intercept
    estimate <-
      fit_netTPreg(~sex, synth_idm_data,
        ratetable = NULL,
        rmap = NULL,
        time_dep_popvars = NULL,
        s = s_time,
        t = 3.0,
        by = n_ind / 10,
        trans = transition,
        link = "logit",
        R = 10 # Number of bootstraps
      )
    testthat::expect_no_error(
      plot(estimate)
    )
  }
})
