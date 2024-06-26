% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survexp_wrappers.R
\name{compute_survprob_pch}
\alias{compute_survprob_pch}
\title{Apply individual surv prob computation to DF}
\usage{
compute_survprob_pch(patientsDF, eval_times, ratetable, rmap)
}
\arguments{
\item{patientsDF}{a DF/DT/tibble containing the data needed to compute
individual expected survival probability}

\item{eval_times}{The times at which survival must be evaluated. Times must
be in days if using a ratetable with a Date component. Provided times can
be passed as single numeric value or vector of numeric values or as unamed
list of numeric values if common for all individuals. Otherwise the name of
the column containing the list,vector or single value evaluation point can
be passed as character string or as variable name.}

\item{ratetable}{A ratetable object. See ?survival::ratetable for details.}

\item{rmap}{An rmap argument as the one passed to survival::survexp. This
argument will only be evaluated inside survexp}
}
\value{
Returns a tidy tibble with 3 columns: the row name of the individual
row (as given by the call of row.names() on the patientsDF), the evaluation
time and the corresponding survival probability.
}
\description{
This function is a wrapper around \code{indiv_survprob_pch}. It calls
the latter function on every row of the DF in order to compute the expected
survival probability for each individual row at each provided \code{eval_times}.
The point of this is to call survexp only once per patient and compute the
survival probablity at all patime points for each patient, thus avoiding
redundant computation.
}
\examples{
compute_survprob_pch(patients,eval_times = c(10,20,30),ratetable = slopop,rmap=list(year=date_chir, age=age_at_dg, sex=sexe))
}
