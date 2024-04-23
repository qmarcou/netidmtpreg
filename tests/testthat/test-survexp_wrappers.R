test_that("Test compute_survprob_pch", {
  data_df <- survival::jasa
  rt <- survival::survexp.us
  eval_times <- c(10, 20, 30)
  # Check several rmap argument styles following examples from survexp doc
  ## data dependent style (survexp example)
  testthat::expect_no_error(
    exp_surv_df <- compute_survprob_pch(
      data_df,
      eval_times = eval_times,
      ratetable = rt,
      rmap = list(
        sex = "male",
        year = accept.dt,
        age = (accept.dt - birth.dt)
      )
    )
  )
  ### Assert returns probabilities
  exp_surv_df %>% assertr::assert(assertr::within_bounds(0.0,1.0), expsurvs)
  ### Assert eval_times respected
  exp_surv_df %>% assertr::assert(assertr::in_set(eval_times), eval_time)
  testthat::expect_equal(nrow(exp_surv_df), nrow(data_df)*length(eval_times))

  ## hardcoded style
  testthat::expect_no_error(
    exp_surv_df <- compute_survprob_pch(
      data_df,
      eval_times = eval_times,
      ratetable = rt,
      rmap = list(
        sex = "male",
        year = accept.dt,
        age = 65 * 365 # age in days
      )
    )
  )
  testthat::expect_no_error(
    exp_surv_df <- compute_survprob_pch(
      data_df,
      eval_times = eval_times,
      ratetable = rt,
      rmap = list(
        sex = "male",
        year = as.Date("1970","%Y"),
        age =  (accept.dt - birth.dt) # age as a difftime object
      )
    )
  )
  testthat::expect_no_error(
    exp_surv_df <- compute_survprob_pch(
      data_df,
      eval_times = eval_times,
      ratetable = rt,
      rmap = list(
        sex = "male",
        year = 1970,
        age = 65 * 365 # age in days
      )
    )
  )

  testthat::skip("This test fails due to a recylcling problem with difftime(?).")
  testthat::expect_no_error(
    exp_surv_df <- compute_survprob_pch(
      data_df,
      eval_times = eval_times,
      ratetable = rt,
      rmap = list(
        sex = "male",
        year = accept.dt,
        age = as.difftime(65*365, units = "days") # age as a difftime object
      )
    )
  )

  testthat::skip("These tests fail with segfault from survexp, I do not understand why.")
  # FIXME This test is just a combination of the two above
  testthat::expect_no_error(exp_surv_df_2 <- compute_survprob_pch(
    data_df,
    eval_times = eval_times,
    ratetable = rt,
    rmap = list(
      sex = "male",
      year = as.Date("1970","%Y"), # Today's day and month but in 1970
      age = as.difftime(65*365,units = "days") # age in days
    )
  ))
  testthat::expect_no_error(exp_surv_df_2 <- compute_survprob_pch(
    data_df,
    eval_times = eval_times,
    ratetable = rt,
    rmap = list(
      sex = "male",
      year = as.Date("1970","%Y"), # Today's day and month but in 1970
      age = 65*365 # age in days
    )
  ))
})
