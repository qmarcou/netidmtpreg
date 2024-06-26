% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/error.R
\name{multi_assert}
\alias{multi_assert}
\title{Perform multiple assertions raising a common aborting class.}
\usage{
multi_assert(x, x_name, abort_class, assert_list)
}
\arguments{
\item{x}{any, object to be tested}

\item{x_name}{string, name of the object to be tested}

\item{abort_class}{callable, the custom abort class}

\item{assert_list}{a vector or list of asserting functions}
}
\value{
None
}
\description{
Perform multiple assertions raising a common aborting class.
}
