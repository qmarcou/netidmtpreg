# Transform a TPreg's coefficients related data to a single tibble
as_tibble.TPreg <- function(x, tmpsep = ".-_-.") {
  # Get covariate names
  var_names <- all.vars(x$co$formula)
  # Check whether formula implies an intercept
  # check existence of an intercept https://stackoverflow.com/questions/50073814/create-expression-test-to-see-if-as-formulax-1-in-r-contains-an-intercept
  if (attr(stats::terms(x$co$formula), "intercept")) {
    var_names <- c("(Intercept)", var_names)
  }

  # Gather all information in one messy tibble with suffixes
  raw_tib <- tibble::tibble(t = x$co$time)
  for (qty in c("coefficients", "SD", "LWL", "UPL", "p.value")) {
    raw_tib <- dplyr::bind_cols(
      raw_tib,
      tibble::as_tibble(x$co[[qty]]) %>%
        dplyr::rename_with(.fn = function(x) {
          return(paste0(x, tmpsep, qty))
        }, .cols = dplyr::everything())
    )
  }

  # Now make tidy tibbles for each variable and combine them
  tidy_tib <- tibble::tibble()
  for (varname in var_names) {
    tmp_tib <- dplyr::bind_cols(
      tibble::tibble(t = x$co$time, covar = varname), raw_tib %>%
        dplyr::select(dplyr::starts_with(varname))
    )
    tidy_tib <- dplyr::bind_rows(
      tidy_tib,
      tmp_tib %>%
        tidyr::pivot_longer(
          cols = !c("t", "covar"),
          names_to = c("covar.val", "qty"),
          names_prefix = varname,
          names_sep = tmpsep
        )
    )
  }
  tidy_tib <- tidy_tib %>%
    tidyr::pivot_wider(values_from = value, names_from = qty)
  return(tidy_tib)
}

combine_TPreg_tidy <- function(TPregs_objects, model_names = NULL) {
  if (is.null(names(TPregs_objects)) & is.null(model_names)) {
    stop("You must provide names for the TPreg objects, either via a names
      object or via the model_names argument")
  }
  if (rlang::is_null(model_names)) {
    model_names <- names(TPregs_objects)
  }
  TPregs_objects <- as.list(TPregs_objects)
  names(TPregs_objects) <- model_names
  TPregs_objects <- lapply(TPregs_objects, as_tibble.TPreg)
  final_tib <- tibble::tibble(
    model.name = names(TPregs_objects),
    tidyTPreg = TPregs_objects
  ) %>% tidyr::unnest(tidyTPreg)
  return(final_tib)
}
