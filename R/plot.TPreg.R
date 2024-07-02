#' @importFrom ggplot2 autolayer
#' @export
autolayer.TPreg <-
  function(object, ..., model = NULL) {
    # https://stackoverflow.com/a/7099056
    # https://ggplot2.tidyverse.org/reference/automatic_plotting.html
    # https://stackoverflow.com/a/47491926

    # Create tiddy tibble
    tidy_tpreg <- tidy.TPreg(object)
    # Add model name to
    if (!is.null(model)) {
      tidy_tpreg <- tidy_tpreg %>% tibble::add_column(model = model)
    }

    # Create the ggplot object
    ggplot_layer <- geom_lineci(data = tidy_tpreg, ...)
    return(ggplot_layer)
  }

#' @importFrom ggplot2 autoplot
#' @export
autoplot.TPreg <-
  function(object, ..., model = NULL) {
    # Workaround to avoid cutting out ribbon confidence intervals
    # https://stackoverflow.com/a/38777929
    # should be mentioned in README/vignette

    tidy_tpreg <- tidy.TPreg(object)

    ggplot_obj <- if (is.null(model)) {
      ggplot2::ggplot(
        data = tidy_tpreg,
        mapping = ggplot2::aes(color = covar.val)
      )
    } else {
      ggplot2::ggplot(
        data = tidy_tpreg,
        mapping = ggplot2::aes(color = covar.val, linetype = model)
      )
    }

    ggplot_obj <- ggplot_obj +
      autolayer.TPreg(object = object, ...) +
      ggplot2::geom_hline(yintercept = 0, color = "red") +
      ggplot2::geom_vline(xintercept = object$s) +
      ggplot2::facet_wrap(~covar)
    return(ggplot_obj)
  }


#' @export
plot.TPreg <-
  function(x, ...) {
    autoplot.TPreg(object = x, ...)
  }


geom_lineci <- function(mapping = NULL, data = NULL,
                        position = "identity",
                        ...,
                        fill = "#888888",
                        alpha = 0.2,
                        size = 0) {
  # Some references:
  # https://stackoverflow.com/a/32616695
  # https://github.com/tidyverse/ggplot2/blob/main/R/geom-smooth.R
  # https://ggplot2-book.org/extensions#sec-new-geoms
  list(
    ggplot2::geom_ribbon(
      data = data,
      mapping = ggplot2::aes(x = t, ymin = LWL, ymax = UPL),
      fill = fill,
      alpha = alpha,
      size = size,
      ...
    ),
    ggplot2::geom_line(
      data = data,
      mapping = ggplot2::aes(x = t, y = coefficients, ...)
    )
  )
}
