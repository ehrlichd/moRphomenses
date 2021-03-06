% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/align.R
\name{mm_ArrayData}
\alias{mm_ArrayData}
\title{Array Data}
\usage{
mm_ArrayData(
  ObsIDs,
  ObsDays,
  ObsValue,
  ObsMid,
  StartDay = 1,
  EndDay = NULL,
  ScaleTo,
  ScaleToMid = NULL
)
}
\arguments{
\item{ObsIDs}{A vector that contains indiviaul IDs repeated for muliple days of collection.}

\item{ObsDays}{A vector that contains information on time, IE Day 1, Day 2, Day 3. Note: this vector should include integers, continuous data might produce unintended results.}

\item{ObsValue}{A vector containing the variable sampled.}

\item{ObsMid}{A vector containng the midpoint day for each individual. Note: ObsMid must have the same number of observations as unique Individuals.}

\item{StartDay}{Default is starting at ObsDay 1, can specify other values to subsample.}

\item{EndDay}{If NULL (default), the highest ObsDay is used for each individual.}

\item{ScaleTo}{Integer. Number of days to up/down sample observations to using \code{\link{mm_GetInterval}}.}

\item{ScaleToMid}{If NULL (default) 0 will be centered using mm_interval.}
}
\value{
Returns a 3D array of data to be analyzed with individuals in the 3rd dimension.
}
\description{
Construct a ragged array (containing missing data) of a specified length (up/down sampling individuals to fit).
}
