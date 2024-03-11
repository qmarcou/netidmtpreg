devtools::load_all()
testthat::test_that("Test convergence error wrapper for mod.glm.fit.", {
  # TODO
  testthat::skip("not implemented")
})

testthat::test_that("Test calling wrapper for mod.glm.fit.", {
  # Completely random data
  n <- 100
  dim_x <- 3
  y <- rbinom(n = 100, size = 1, prob = .5)
  x <- rnorm(dim_x * n) %>% array(dim = c(n, dim_x))
  w <- runif(n, min = 0, max = 1)
  testthat::expect_no_error(mod.glm.fit.callingwrapper(
    X = x,
    response = y,
    weights = w,
    family = binomial()
  ))
  # Uninformative response vector returns NA coefficients
  beta <- mod.glm.fit.callingwrapper(
    X = x,
    response = rep.int(0, times = n),
    weights = w,
    family = binomial()
  )
  testthat::expect_true(all(rlang::are_na(beta)))


  beta <- mod.glm.fit.callingwrapper(
    X = x,
    response = rep.int(1, times = n),
    weights = w,
    family = binomial()
  )
  testthat::expect_true(all(rlang::are_na(beta)))
})

testthat::test_that("IDM crude survival model Fitting", {
  # Test single dimensionnal model.matrix (e.g intercept only formula ~ 1)
  # Fixed by feb8378ef0a17d1aca30f2fda55ec57c77711e64
  # TODO test crude mortality estimates correctness
  testthat::skip("not implemented")
  # Crude mortality binomial regression
  synth_idm_data <- generate_uncensored_ind_exp_idm_data(
    n_individuals = 1e4,
    lambda_illness = 1.0,
    lambda_death = 0.1
  )
  for (transition in c("all", "11", "12", "22", "13", "23")) {
    renewnetTPreg(~1, synth_idm_data, ratetable = NULL, s = 0, t = 10, by = 1, trans = transition)
  }
})

testthat::test_that("IDM Net survival model Fitting", {
  # Check that the model runs even with nonsense population information
  n_ind <- 1e4
  synth_idm_data <- generate_uncensored_ind_exp_idm_data(
    n_individuals = n_ind,
    lambda_illness = 1.0,
    lambda_death = 0.1
  )
  # Generate random age and sex labels
  synth_idm_data <-
    synth_idm_data %>% tibble::add_column(
      sex = ifelse(rbinom(n_ind, 1, prob = .5), "male", "female"),
      age = runif(n = n_ind, min = 50, max = 80)
    )
  # Generate random start of follow up dates
  # FIXME a date before 1940 or after 2012 (limits of uspop ratetable) is
  # extremely unlikely with these parameters but not impossible.
  synth_idm_data <-
    synth_idm_data %>% tibble::add_column(start_date = as.Date.numeric(
      x = rnorm(n = n_ind, mean = 0, sd = 1e2),
      origin = as.Date("15/06/1976", "%d/%m/%Y")
    ))
  for (transition in c("11")) {
    renewnetTPreg(
      formula = ~1,
      synth_idm_data,
      # Use a standard ratetable
      ratetable = survival::survexp.us,
      rmap = list(
        age = age,
        sex = sex,
        year = start_date
      ),
      time_dep_popvars = list("age", "year"),
      s = 0,
      t = 1.5,
      by = n_ind / 2,
      trans = transition,
      link = "logit",
      R = 2 # Number of bootstraps
    )
  }
  testthat::skip("not implemented")
})

testthat::test_that("Test single time point estimation", {
  # TODO alleviate code duplication?
  n_ind <- 1e4
  synth_idm_data <- generate_uncensored_ind_exp_idm_data(
    n_individuals = n_ind,
    lambda_illness = 1.0,
    lambda_death = 0.1
  )
  # Generate random age and sex labels
  synth_idm_data <-
    synth_idm_data %>% tibble::add_column(
      sex = ifelse(rbinom(n_ind, 1, prob = .5), "male", "female"),
      age = runif(n = n_ind, min = 50, max = 80)
    )
  # Generate random start of follow up dates
  # FIXME a date before 1940 or after 2012 (limits of uspop ratetable) is
  # extremely unlikely with these parameters but not impossible.
  synth_idm_data <-
    synth_idm_data %>% tibble::add_column(start_date = as.Date.numeric(
      x = rnorm(n = n_ind, mean = 0, sd = 1e2),
      origin = as.Date("15/06/1976", "%d/%m/%Y")
    ))
  # Fit a GLM estimate
  testthat::expect_no_error(
    fit_single_time_point_estimate(
      s = 0,
      # No time correction needed in rmap for s=0
      t = 1.5,
      transition = "11",
      X = array(1, dim = c(n_ind, 1)),
      # Intercept only model.matrix
      data_df = synth_idm_data,
      ratetable = survival::survexp.us,
      rmap = list(year = start_date, sex = sex, age = age)
    )
  )
})

testthat::test_that("Check censoring dist fitting and prediction", {
  # Check that censoring fitting and prediction works without censoring
  n_ind <- 1e4
  synth_idm_data <- generate_uncensored_ind_exp_idm_data(
    n_individuals = n_ind,
    lambda_illness = 1.0,
    lambda_death = 0.1
  )
  cens_fit <-
    estimate_censoring_dist(
      s = 0,
      t = 1.5,
      X = NULL,
      data_df = synth_idm_data
    )
  # get censoring estimate at boundary
  assertthat::are_equal(get_survival_at(1.5, cens_fit), 1.0)
  testthat::skip("not implemented")
  # TODO Add actual censoring
})
