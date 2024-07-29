#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#



#' Launch mm_Explorer
#'
#' @param dataset The dataset to be analyzed
#' @export

mm_Explorer <- function(dataset=NULL){
  if(is.null(dataset)){
    dataset <- mm_data
  }

  ## I don't expect this to work, but give it a try...otherwise,
  ## CUT AND PASTE THE FINAL VERSION FROM app.R

  runApp(appDir = system.file("app.R", package = "moRphomenses"))


}
