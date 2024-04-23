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

# Crude survival testing
testthat::test_that("IDM crude survival model estimates", {
  # Test single dimensional model.matrix (e.g intercept only formula ~ 1)
  # Fixed by feb8378ef0a17d1aca30f2fda55ec57c77711e64
  # Crude mortality binomial regression
  set.seed(246852389)
  n_ind <- 1e5
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
  estimates <- list()
  for (transition in c("11", "12", "13", "23")) {
    estimates[transition] <-
      renewnetTPreg(~1, synth_idm_data,
        ratetable = NULL,
        rmap = NULL,
        time_dep_popvars = NULL,
        s = s_time,
        t = 3.0,
        by = n_ind / 10,
        trans = transition,
        link = "logit",
        R = 1 # Number of bootstraps
      )
  }

  # Check agreement with generating parameters
  expected_tp <- collections::dict()
  ## Check 11 transition estimates
  expected_tp$set("11", pexp(estimates[["11"]]$time - s_time,
    rate = (l_illness + l_death),
    lower = FALSE # P(T>t)
  ))

  ## Check 13 transition estimates
  # In this simple test death rate is the same starting from state 1 or 2
  # therefore state 2 does not influence 13 transition
  expected_tp$set("13", pexp(estimates[["13"]]$time - s_time,
    rate = l_death,
    lower = TRUE # P(T<=t)
  ))

  ## Check 12 transition estimates: 12 transition and no death
  p_illness <- pexp(estimates[["12"]]$time - s_time,
    rate = l_illness,
    lower = TRUE # P(T<=t)
  )
  p_no_death <- pexp(estimates[["12"]]$time - s_time,
    rate = l_death,
    lower = FALSE # P(T>t)
  )
  expected_tp$set("12", p_illness * p_no_death)

  ## Check 23 transition estimates
  expected_tp$set("23", pexp(estimates[["23"]]$time - s_time,
    rate = l_death,
    lower = TRUE # P(T<=t)
  ))

  ## Check 22 transition estimates: survival probability
  expected_tp$set("22", pexp(estimates[["22"]]$time - s_time,
    rate = l_death,
    lower = FALSE # P(T<t)
  ))

  for (transition in c("11", "12", "13", "23")) {
    tp_val <- expected_tp$get(transition)

    testthat::expect_equal(
      object = as.numeric(expit(estimates[[!!transition]]$coefficients)),
      expected = tp_val,
      tolerance = .05 # ~5% relative difference tolerance 
    )
  }
})

testthat::test_that("IDM crude survival regression gives correct
                      estimates in absence of covariate and censoring 
                      and exponentially distributed events", {

                      })

testthat::test_that("IDM Net survival model fitting runs inside futures", {
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
  for (session_type in c("sequential", "multisession")){
    future::plan(session_type)
    for (transition in c("all", "11", "12", "13", "23")) {
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
        s = .5,
        t = 1.5,
        by = n_ind / 2,
        trans = transition,
        link = "logit",
        R = 2 # Number of bootstraps
      )
    }
    future::plan(future::sequential)
  }
})

testthat::test_that("IDM Net survival fitting arguments handling", {
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

  # vector-like t with adjusted_t gives is estimated at expected times
  t <- c(1.3, 1.4, 1.5)
  estim <- renewnetTPreg(
    s = .5,
    t = t,
    trans = "11",
    formula = ~1,
    data = synth_idm_data,
    ratetable = NULL,
    R = 2, # Number of bootstraps,
    readjust_t = FALSE
  )
  testthat::expect_equal(object = estim$co$time, expected = t)
  testthat::expect_error( # large t in unadjusted vector of times
    estim <- renewnetTPreg(
      s = .5,
      t = c(t, 1e5),
      trans = "11",
      formula = ~1,
      data = synth_idm_data,
      ratetable = NULL,
      R = 2, # Number of bootstraps,
      readjust_t = FALSE
    )
  )
  testthat::skip("not implemented")
  # trans = "all" returns all transitions
})

testthat::test_that(
  "Estimated IDM Net survival and crude survival without
population mortality are equal",
  {
    # Check that the model runs even with nonsense population information
    n_ind <- 1e5
    l_illness <- 1.0
    l_death <- 1.0
    l_pop_death <- 1.0
    s_time <- 0
    synth_idm_data <- generate_uncensored_ind_exp_idm_data(
      n_individuals = n_ind,
      lambda_illness = l_illness,
      lambda_death = l_death
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

    # Generate population mortality assuming equal constant population rate
    population_death_times <- generate_exponential_time_to_event(
      n_individuals = n_ind,
      lambda = l_pop_death
    )
    # create a corresponding ratetable object
    const_ratetable <- survival::survexp.us
    const_ratetable[] <- l_pop_death
    # Update death time accordingly to create an observed crude survival dataset
    crude_synth_idm_data <- apply_iddata_death(
      synth_idm_data,
      population_death_times
    )

    for (transition in c("11")) {
      # FIXME: need to align times
      # net_truth <- renewnetTPreg(
      #   formula = ~1,
      #   synth_idm_data,
      #   # Use a standard ratetable
      #   ratetable = NULL,
      #   rmap = list(
      #     age = age,
      #     sex = sex,
      #     year = start_date
      #   ),
      #   time_dep_popvars = list("age", "year"),
      #   s = 0,
      #   t = 1.5,
      #   by = n_ind / 2,
      #   trans = transition,
      #   link = "logit",
      #   R = 1 # Number of bootstraps
      # )
      net_estimated <- renewnetTPreg(
        formula = ~1,
        crude_synth_idm_data,
        # Use a standard ratetable
        ratetable = const_ratetable,
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
        R = 1 # Number of bootstraps
      )
      # FIXME: this cannot work, times are not aligned
      # testthat::expect_equal(
      #   object = net_estimated$co$coefficients,
      #   expected = net_truth$co$coefficients,
      #   tolerance = .01
      # )
      # Compute expected coefficients using survfit
      net_survfit <- survival::survfit(
        survival::Surv(time = pmin(Zt, Tt), event = delta) ~ 1,
        data = synth_idm_data
      )
      net_surv_probs <- get_survival_at(net_estimated$co$time, net_survfit)
      # net_surv_probs <- pexp(net_estimated$co$time,
      #     rate = (l_illness + l_death),
      #     lower = FALSE # P(T>t)
      #   )
      testthat::expect_equal(
        object = as.numeric(expit(net_estimated$co$coefficients)),
        expected = net_surv_probs,
        tolerance = .01
      )
    }

    testthat::skip("not implemented")
    # Generate population death times using a ratetable
  }
)

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

  # Test vectorized input
  testthat::expect_equal(
    c(1.0, 1.0, .75, .75, .4, .20, .20),
    shorthand_fun(c(0.0, .5, 2.0, 3.0 - 1e-5, 5, 6.0, 7.0))
  )

  # Test edge cases
  ## slightly negative value
  testthat::expect_error(
    shorthand_fun(-1e-10),
    class = "invalid_argument_error"
  )
  testthat::expect_error(
    shorthand_fun(c(1, 2, 3, -1e-10, 5, 6)),
    class = "invalid_argument_error"
  )
  ## Infinite positive value
  testthat::expect_equal(.20, shorthand_fun(Inf))
  ## Single breakpoint survfit_df
  survfit_df_2 <- tibble::tibble(
    time = c(0),
    surv = c(1.0)
  )
  shorthand_fun_2 <- function(x) get_survival_at(x, survfit_df_2)
  testthat::expect_equal(c(1.0, 1.0, 1.0), shorthand_fun_2(c(0, 1, Inf)))
  testthat::expect_error(
    shorthand_fun_2(-1e-10),
    class = "invalid_argument_error"
  )
  testthat::expect_error(
    shorthand_fun_2(c(-1e-10, 1)),
    class = "invalid_argument_error"
  )
  survfit_df_3 <- tibble::tibble(
    time = c(1.0), # not defined on 0.0 this time
    surv = c(1.0)
  )
  shorthand_fun_3 <- function(x) get_survival_at(x, survfit_df_3)
  testthat::expect_equal(c(1.0, 1.0), shorthand_fun_3(c(1, Inf)))
  testthat::expect_error(
    shorthand_fun_3(1 - 1e-10),
    class = "invalid_argument_error"
  )
  # survival probability greater than 1
  survfit_df_4 <- tibble::tibble(
    time = c(0, 1),
    surv = c(1.1, .85)
  )
  testthat::expect_error(
    get_survival_at(0.0, survfit_df_4),
    class = "invalid_argument_error"
  )
  # negative survival probability
  survfit_df_5 <- tibble::tibble(
    time = c(0, 1),
    surv = c(1.0, -1e-5)
  )
  testthat::expect_error(
    get_survival_at(0.0, survfit_df_5),
    class = "invalid_argument_error"
  )
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
  # Test with survfit object
  survfit_obj <- survival::survfit(
    survival::Surv(c(1, 2, 3, 4), c(1, 1, 1, 1)) ~ 1
  )
  testthat::expect_equal(
    get_survival_at(c(0, 1.5, 5), survfit_obj),
    c(1.0, .75, .0)
  )
  # Test summary.survfit object
  sum_survfit <- summary(survfit_obj)
  testthat::expect_equal(
    get_survival_at(c(0, 1.5, 5), survfit_obj),
    c(1.0, .75, .0)
  )
})

testthat::test_that("Test summarize_single_time_bootstraps", {
  mock_boot_df <-
    tibble::tibble("(Intercept)" = rnorm(1e5,mean = 1, sd = 1),
                   "X1" = rnorm(1e5, mean = 0, sd = 1))
  testthat::expect_no_error(
    boot_summary <- summarize_single_time_bootstraps(mock_boot_df)
  )
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
