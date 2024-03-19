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

devtools::dev_mode(on = TRUE)
devtools::install_local()
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
  for(session_type in c("sequential", "multisession")){
    future::plan(session_type)
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
    future::plan(future::sequential)
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
  # Add two censoring events to check survfit_df ordering
  synth_idm_data[as.integer(c(n_ind/2,n_ind/3)), "delta"] <- 0
  cens_fit <-
    estimate_censoring_dist(
      s = 0,
      t = 1.5,
      X = NULL,
      data_df = synth_idm_data
    )
  # get censoring estimate at boundary
  testthat::expect_no_error(get_survival_at(1.5, cens_fit))

  testthat::skip("not implemented")
  # TODO Add actual censoring
})

testthat::test_that("Test get_survival_at function", {
  survfit_df <- tibble::tibble(
    time = c(0, 1, 2, 3, 4, 5, 6),
    surv = c(1.0, .85, .75, .60, .55, .40, .20)
  )
  # define a closure for testing
  shorthand_fun <- function(x) get_survival_at(x, survfit_df)

  # Test normal scalar usage
  testthat::expect_equal(1.0, shorthand_fun(0.0))
  testthat::expect_equal(1.0, shorthand_fun(.5))
  testthat::expect_equal(.75, shorthand_fun(2.0))
  testthat::expect_equal(.75, shorthand_fun(3.0 - 1e-5))
  testthat::expect_equal(.4, shorthand_fun(5))
  testthat::expect_equal(.20, shorthand_fun(6.0))
  testthat::expect_equal(.20, shorthand_fun(7.0))

  # Test vectorised input
  testthat::expect_equal(
    c(1.0, 1.0, .75, .75, .4, .20, .20),
    shorthand_fun(c(0.0, .5, 2.0, 3.0 - 1e-5, 5, 6.0, 7.0))
  )

  # Test edge cases
  ## slightly negative value
  testthat::expect_error(shorthand_fun(-1e-10), class = "invalid_argument_error")
  testthat::expect_error(shorthand_fun(c(1, 2, 3, -1e-10, 5, 6)), class = "invalid_argument_error")
  ## Infinite positive value
  testthat::expect_equal(.20, shorthand_fun(Inf))
  ## Single breakpoint survfit_df
  survfit_df_2 <- tibble::tibble(
    time = c(0),
    surv = c(1.0)
  )
  shorthand_fun_2 <- function(x) get_survival_at(x, survfit_df_2)
  testthat::expect_equal(c(1.0, 1.0, 1.0), shorthand_fun_2(c(0, 1, Inf)))
  testthat::expect_error(shorthand_fun_2(-1e-10), class = "invalid_argument_error")
  testthat::expect_error(shorthand_fun_2(c(-1e-10, 1)), class = "invalid_argument_error")
  survfit_df_3 <- tibble::tibble(
    time = c(1.0), # not defined on 0.0 this time
    surv = c(1.0)
  )
  shorthand_fun_3 <- function(x) get_survival_at(x, survfit_df_3)
  testthat::expect_equal(c(1.0, 1.0), shorthand_fun_3(c(1, Inf)))
  testthat::expect_error(shorthand_fun_3(1 - 1e-10), class = "invalid_argument_error")
  # survival probability greater than 1
  survfit_df_4 <- tibble::tibble(
    time = c(0, 1),
    surv = c(1.1, .85)
  )
  testthat::expect_error(get_survival_at(0.0, survfit_df_4), class = "invalid_argument_error")
  # negative survival probability
  survfit_df_5 <- tibble::tibble(
    time = c(0, 1),
    surv = c(1.0, -1e-5)
  )
  testthat::expect_error(get_survival_at(0.0, survfit_df_5), class = "invalid_argument_error")
  # Missing columns
  testthat::expect_error(get_survival_at(
    0.0,
    tibble::tibble(time = c(.0))
  ), class = "invalid_argument_error")
  testthat::expect_error(get_survival_at(
    0.0,
    tibble::tibble(surv = c(1.0))
  ), class = "invalid_argument_error")
  # Works with extra columns
  testthat::expect_no_error(get_survival_at(
    0.0,
    survfit_df_2 %>% tibble::add_column(x = c("foo"))
  ))
})

testthat::test_that("Test summarize_single_time_bootstraps", {
  mock_boot_df <-
    tibble::tibble("(Intercept)" = rnorm(1e5,mean = 1, sd = 1),
                   "X1" = rnorm(1e5, mean = 0, sd = 1))
  testthat::expect_no_error(boot_summary <- summarize_single_time_bootstraps(mock_boot_df))
  assertr::verify(
    boot_summary,
    assertr::has_only_names(
      "(Intercept)_sd",
      "(Intercept)_ci.lb",
      "(Intercept)_ci.ub",
      "(Intercept)_n.failed.boot",
      "X1_sd",
      "X1_ci.lb",
      "X1_ci.ub",
      "X1_n.failed.boot"
    )
  )
  # TODO check values
  #assertr::verify()
})
