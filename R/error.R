#' Custom exception for invalid arguments
#'
#' Derived from https://adv-r.hadley.nz/conditions.html#custom-conditions
#' @param arg string, faulty argument name
#' @param must string, description of the expected nature of the argument
#' @param not string, optional, description of the received argument
#'
#' @return None
#' @export
#'
#' @examples
invalid_argument <- function(arg, must, not = NULL) {
  msg <- glue::glue("`{arg}` must {must}")
  if (!is.null(not)) {
    not <- typeof(not)
    msg <- glue::glue("{msg}; not {not}.")
  }

  rlang::abort(
    class = "invalid_argument_error",
    message = msg,
    arg = arg,
    must = must,
    not = not
  )
}

not_implemented <- function(msg){
  rlang::abort(
    message = msg,
    class = "not_implemented_error"
  )
}

assert_scalar <- function(x, x_name, abort_class) {
  if (!assertthat::is.scalar(x))
    abort_class(arg = x_name, must = "be a scalar value.")
}

assert_positive <- function(x, x_name, abort_class) {
  if (any(x < 0))
    abort_class(arg = x_name, must = "be a positive number.")
}

assert_probability <- function(x, x_name, abort_class) {
  if (any(x < 0 | x >1))
    abort_class(arg = x_name, must = "be in [0.0, 1.0].")
}

assert_count <- function(x, x_name, abort_class) {
  if (!assertthat::is.count(x))
    abort_class(arg = x_name, must = "be a scalar positive integer (count) number.")
}

assert_integerish <- function(x, x_name, abort_class) {
  if (!rlang::is_integerish(x))
    abort_class(arg = x_name, must = "be a positive number.")
}

#' Perform multiple assertions raising a common aborting class.
#'
#' @param x any, object to be tested
#' @param x_name string, name of the object to be tested
#' @param assert_list a vector or list of asserting functions
#' @param abort_class callable, the custom abort class
#'
#' @return None
#' @export
#'
#' @examples
multi_assert <- function(x, x_name, abort_class, assert_list) {
  for (asserting_func in assert_list) {
    assertthat::has_args(asserting_func, c("x", "x_name", "abort_class"), exact = TRUE)
    asserting_func(x, x_name, abort_class)
  }
}
