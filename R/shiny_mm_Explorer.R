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
#' @export

mm_Explorer <- function(){

  requireNamespace("shiny", quietly = TRUE)

  # app <- shiny::shinyApp(ui = ui, server = server)
  #
  # shiny::runApp(app)

  # ## system.file not finding the file
  # ## This appears to work if
  shiny::runApp(appDir = system.file("mm_Explorer", package = "moRphomenses"))


}
