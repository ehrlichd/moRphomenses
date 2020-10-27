#' Split data table into a 3D array and a data table of unique classifiers
#'
#' @param rawData A dataframe or table in long format. IE rows contain the same individual multiple times from different sampling events
#'
#' @param ID A column index or name that contains unique IDs of individuals
#' @param  classifiers A list of column indeces or names that contain classifiers to be used in subsequent analysis
#' @param day A vector containing sample day. Intended to be used with integer/factor data. A continuous sample day variable might produce unintended results.
#' @param var A vector containing the variable sampled.
#'
#' @return Returns a data table of n individuals with additional classifier data (EG, groupID, age, sex) as well as a 3D array of data to be analyzed. Data array is in the form [StudyDay x Value X Indivual]
#'  @export

mm_split <- function(rawData, ID, classifiers = list(), day, var){
  id_levs <- as.factor(rawData[ID])
  n <- length(levels(id_levs))


  class <- data.frame("IDs" = id_levs, rawData[classifiers])

  sub <- list()
  rDay <- data.frame("IDs" = id_levs, minDay=numeric(), maxDay=numeric(), maxLength=numeric())
  for (i in 1:n){
    sub[[i]] <- droplevels(rawData[,c(ID, day, var)])
    rDay$minDay[i] <- min(sub[[i]])[2]
    rDay$maxDay[i] <- max(sub[[i]])[2]

  }
  rDay$maxLength <- rDay$minDay - rDay$maxDay

  arr <- array(dim = )
}
