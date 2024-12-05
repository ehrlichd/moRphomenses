

#' Launch mm_Explorer
#'
#' @export

mm_Explorer <- function(){

  requireNamespace("shiny", quietly = TRUE)

  app_dir <- system.file("mm_Explorer_local", package = "moRphomenses")

  shiny::runApp(app_dir, launch.browser = TRUE,display.mode = "normal")


}
