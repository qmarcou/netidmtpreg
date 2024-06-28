autolayer.TPreg <- function(tpreg_obj, ...) {
  # https://stackoverflow.com/a/7099056
  # https://ggplot2.tidyverse.org/reference/automatic_plotting.html
  # https://stackoverflow.com/a/47491926

  # BUG: cannot directly combine geoms this way without a background plot
  # need to create a custom geom to bundele them:
  # https://ggplot2-book.org/extensions#sec-combining-multiple-geoms
  tidy_tpreg <- tidy.TPreg(tpreg_obj)
  ggplot_layer <-  ggplot2::geom_ribbon(
    data = tidy_tpreg,
    mapping = ggplot2::aes(x = t, ymin = LWL, ymax = UPL),
    alpha = 0.2, size = 0
  ) +
  ggplot2::geom_line(
    data = tidy_tpreg,
    mapping = ggplot2::aes(x = t, y = coefficients)
    )
  return(ggplot_layer)
}

autoplot.TPreg <- function(tpreg_obj, ...) {
  my_aes <- ggplot2::aes(
    color = ggplot2::sym("covar.val"),
    fill = ggplot2::sym("covar.val")
  )
  tidy_tpreg <- tidy.TPreg(tpreg_obj)
  ggplot_obj <- ggplot2::ggplot(data = tidy_tpreg, mapping = my_aes) +
    autolayer.TPreg(tpreg_obj = tpreg_obj, ...) +
    ggplot2::geom_hline(yintercept = 0, color = "red") +
    ggplot2::geom_vline(xintercept = s) +
    ggplot2::facet_wrap(~covar)
  return(ggplot_obj)
}


#' @export
plot.TPreg <-
  function(x, covar, rug = TRUE, main, ylab, xlab, Ylim, ...) {
    return(x %>% as_tibble.TPreg() %>% plot_TPregs())
  }

plot_TPregs <- function(TPregs_tidy, s = NULL) {
  my_aes <- list(
    color = ggplot2::sym("covar.val"),
    fill = ggplot2::sym("covar.val")
  )
  if ("model.name" %in% colnames(TPregs_tidy)) {
    my_aes[["linetype"]] <- ggplot2::sym("model.name")
  }
  my_aes <- ggpubr::create_aes(my_aes)
  ggplot_obj <- ggplot2::ggplot(data = TPregs_tidy, mapping = my_aes) +
    ggplot2::geom_ribbon(
      ggplot2::aes(x = t, ymin = LWL, ymax = UPL),
      alpha = 0.2, size = 0
    ) +
    ggplot2::geom_line(ggplot2::aes(x = t, y = coefficients)) +
    ggplot2::geom_hline(yintercept = 0, color = "red") +
    ggplot2::geom_vline(xintercept = s) +
    ggplot2::facet_wrap(~covar)
  return(ggplot_obj)
}
# problem of shared colorscale between faceted plots, cf https://stackoverflow.com/questions/3805029/different-legends-and-fill-colours-for-facetted-ggplot, see bottom answer using ggnewscale

