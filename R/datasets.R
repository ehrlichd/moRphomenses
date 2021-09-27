#' Sample hormone dataset
#'
#' A sample new dataset generated from parameters of n=98 individuals. This dataset is ready for shape analysis and analogous to UNS treatment in Ehrlich et al. Array is easily scaled to replicate GMS or (recommended) MMS.
#'
#' @name mm_data_array
#' @docType data
#' @format An array with dimensions 28 x 2 x 60.
#' \describe{
#'   \item{x}{aligned cycles (range -1:1)}
#'   \item{y}{E1G values (ng/ml)}
#'   \item{z}{individuals(n = 60)}
#'   ...
#' }
#'
"mm_data_array"

#' Sample hormone classifiers
#'
#' Sample classifiers to be paired with sample array. This table contains 60 rows to match the 60 individuals across the third dimension of the array
#' @name mm_data_short
#' @docType data
#' @format A data.frame with 60 rows and 7 columns.
#' \describe{
#'   \item{cycL}{Overall cycle length)}
#'   \item{folL}{follicular phase length)}
#'   \item{lutL}{luteal phase length)}
#'   \item{avgE1G}{average E1G across cycle}
#'   \item{ovDay}{day of ovulation}
#'   \item{grp}{phenotypic subgroup as inferred by Ehrlich et al}
#'   \item{ind}{individual id}
#'   ...
#' }
#'
"mm_data_short"

#' Sample hormone data - long format
#'
#' A sample new dataset generated from parameters of n=98 individuals. Data are provided in "long format" with some values (eg cycle length) repeated alongside daily E1G observations. Data require subsequent processing prior to analysis.
#' @name mm_data_long
#' @docType data
#' @format A data.frame with 1759 rows and 10 columns.
#' \describe{
#'
#'
#'
#'   \item{ovAdjDay}{Cycle Day, aligned by day of ovulation==0}
#'   \item{E1G}{daily value of E1G (nn/ml)}
#'   \item{cycL}{Overall cycle length)}
#'   \item{folL}{follicular phase length)}
#'   \item{lutL}{luteal phase length)}
#'   \item{avgE1G}{average E1G across cycle}
#'   \item{ovDay}{day of ovulation}
#'   \item{grp}{phenotypic subgroup as inferred by Ehrlich et al}
#'   \item{ind}{individual id}
#'   ...
#' }
#'
"mm_data_long"
