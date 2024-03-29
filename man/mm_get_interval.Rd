% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/align.R
\name{mm_get_interval}
\alias{mm_get_interval}
\title{Create equallly spaced intervals.}
\usage{
mm_get_interval(days, day0 = NULL)
}
\arguments{
\item{days}{The number of days(divisions) fit between -1 and 1 (inclusive)}

\item{day0}{If NULL (default), the median integer will be calculated. This produces symmetrical ranges when days = odd number. Can be specified for asymmetric ranges.}
}
\description{
Create a sequence from -1:1 of specified length. MIDpoint (day0) can be
}
\examples{
mm_get_interval(15) ## Symmetrical sequence from -1 to 1 with 0 in the middle.
mm_get_interval(15, day0 = 8) ## The same sequence, explicitly specifying the midpoint

mm_get_interval(15, day0 = 3) ## 15 divisions with an asymmetric distribution.

}
