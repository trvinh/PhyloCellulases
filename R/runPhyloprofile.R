#' Run PhyloCellulase app
#' @export
#' @return A shiny application - GUI version of PhyloProfile
#' @import BiocStyle
#' @import DT
#' @importFrom colourpicker colourInput
#' @import energy
#' @import shinyBS
#' @import shinycssloaders
#' @import PhyloProfile
#' @import ape
#' @import data.table
#' @import ggplot2
#' @rawNamespace import(RCurl, except = reset)
#' @rawNamespace import(shinyjs, except = colourInput)

runPhylopPCD <- function(){
    appDir <- system.file("PhyloProfile", package = "PhylopPCD")
    if (appDir == "") {
        stop(
            "Could not find apps directory. Try re-installing `PhylopPCD`.",
            call = FALSE
        )
    }

    shiny::runApp(
        appDir,
        launch.browser = TRUE,
        display.mode = "normal"
    )
}
