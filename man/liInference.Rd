% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/liInference.R
\name{liInference}
\alias{liInference}
\title{Perform hyperparameter inference assuming the model given in Li et al. (2022) on a cycle length dataset.}
\usage{
liInference(
  Y,
  cluster,
  S = 10,
  startingParams = c(kappa = 180, gamma = 6, alpha = 2, beta = 20)
)
}
\arguments{
\item{Y}{A vector of observed cycle lengths.}

\item{cluster}{A vector indicating the individual cluster/group membership for each observation Y.}

\item{S}{Maximum number of possible skipped cycles (see Li et al. for details).}

\item{startingParams}{A vector of starting values for hyperparameters
(default values from Li et al.).}
}
\value{
A list containing the results of hyperparameter inference.
}
\description{
This function performs hyperparameter inference on a given dataset of individuals
and their tracked cycles, assuming the model specified in Li et al. (2022).
Default starting values for hyperparameters and optimization tuning parameters
are those given in Li et al.
}
\references{
Li, Kathy, et al. "A predictive model for next cycle start date that accounts for adherence in menstrual self-tracking." Journal of the American Medical Informatics Association 29.1 (2022): 3-11.
}
