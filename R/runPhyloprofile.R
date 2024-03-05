#' Run PhyloCellulase app
#' @export
#' @return A shiny application - GUI version of PhyloProfile
#' @import BiocStyle
#' @import DT
#' @importFrom colourpicker colourInput
#' @import energy
#' @import shinyBS
#' @import shinycssloaders
#' @rawNamespace import(RCurl, except = reset)
#' @rawNamespace import(shinyjs, except = colourInput)

runPhyloCellulase <- function(){
    appDir <- system.file("PhyloProfile", package = "PhyloCellulase")
    if (appDir == "") {
        stop(
            "Could not find apps directory. Try re-installing `PhyloCellulase`.",
            call = FALSE
        )
    }

    shiny::runApp(
        appDir,
        launch.browser = TRUE,
        display.mode = "normal"
    )
}
