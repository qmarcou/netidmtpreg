% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_netTPreg.R
\name{rellogit}
\alias{rellogit}
\title{An relative survival logit link function closure.}
\usage{
rellogit(s, t, data_df, ratetable, rmap)
}
\arguments{
\item{s}{positive scalar.}

\item{t}{positive scalar.}

\item{data_df}{A \code{data.frame} containing population characteristics.}

\item{ratetable}{a \code{ratetable} object.}

\item{rmap}{see survival::survexp.fit.}
}
\value{
a \code{link-glm} object
}
\description{
Returns a offset logit \code{link-glm} object taking into account population
mortality a the considered times s and t. This relative logit function
arises when decomposing observed transition probabilities into net
transitions and population survival for all transitions where death (3) is
not the target state.
}
