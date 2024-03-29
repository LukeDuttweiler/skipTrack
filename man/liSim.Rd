% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simTrackData.R
\name{liSim}
\alias{liSim}
\title{Simulate user tracked menstrual cycle data for an individual using the li method.}
\usage{
liSim(i, skipProb, maxCycles)
}
\arguments{
\item{i}{Individual identifier. Character, numeric or integer.}

\item{skipProb}{Vector, ignored for this method.}

\item{maxCycles}{Integer, Maximum possible number of true cycles per tracked cycle.}
}
\value{
A data.frame with columns 'Individual', 'TrackedCycles',
'NumTrue', 'Mean', 'SkipProb'.
}
\description{
This function generates synthetic data for user tracked menstrual cycles for a
single individual using the li method.
}
\seealso{
\code{\link{simTrackData}}
}
