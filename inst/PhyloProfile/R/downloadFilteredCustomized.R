#' Download filtered data from customized profile
#' @param data MAIN data for downloading (from module "downloadFilteredMain.R"
#' @param fasta fasta sequences (from reactive fn "customizedFastaDownload")
#' @param inSeq selected sequences in customized profile (from input$inSeq)
#' @param inTaxa selected taxa in customized profile (from input$inTaxa)
#' @return data of customized profile for downloading
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

downloadFilteredCustomizedUI <- function(id) {
    ns <- NS(id)

    tabPanel(
        "Customized data",
        column(
            12,
            strong(
                em(
                    "NOTE: Depending on your choice in Download filtered data
                    -> Main data, either all or only representative sequences
                    will be downloaded!"
                ),
                style = "color:red"
            ),
            hr()
        ),
        column(
            12,
            DT::dataTableOutput(ns("filteredCustomData"))
        ),
        column(
            4,
            downloadButton(ns("downloadCustomData"),
                           "Download customized data")
        ),
        # column(
        #     4,
        #     downloadButton(ns("downloadCustomFasta"),
        #                    "Download FASTA sequences"),
        #     uiOutput(ns("downloadCustomFasta.ui"))
        # ),
        column(
            4,
            downloadButton(ns("downloadCustomLong"),
                           "Download data as PhyloProfile input format")
        )
    )
}

downloadFilteredCustomized <- function(
    input, output, session, data, fasta, inSeq, inTaxa, var1, var2
){
    # filtered data for downloading (Customized Profile) -----------------------
    downloadCustomData <- reactive({
        req(inSeq())
        req(inTaxa())
        data <- as.data.frame(data())
        # get subset of data according to selected genes/taxa
        if (!is.null(inSeq()) | !is.null(inTaxa())) {
            if (inSeq()[1] != "all" & inTaxa()[1] == "all") {
                # select data for selected sequences only
                customData <- subset(data, geneID %in% inSeq())
            } else if (inSeq()[1] == "all" & inTaxa()[1] != "all") {
                # select data for selected taxa only
                customData <- subset(data, fullName %in% inTaxa())
            } else if (inSeq()[1] != "all" & inTaxa()[1] != "all") {
                # select data for selected sequences and taxa
                customData <- subset(data, geneID %in% inSeq()
                                      & fullName %in% inTaxa())
            } else {
                customData <- data
            }
        } else {
            customData <- data
        }
        
        # filter data 
        if (length(var1()) == 1) {
            customData <- customData[customData$var1 >= var1()[1], ]
        } else {
            if (max(customData$var1) <= var1()[2]) {
                customData <- customData[
                    customData$var1 >= var1()[1] & customData$var1 <= var1()[2], 
                ]
            }
        }
        if (!all(is.na(customData$var2))) {
            if (length(var2()) == 1) {
                customData <- customData[customData$var2 >= var2()[1], ]
            } else {
                customData <- customData[
                    customData$var2 >= var2()[1] & customData$var2 <= var2()[2], 
                ]
            }
        } else {
            customData$var2 <- 0
        }
        customData <- customData[!is.na(customData$geneID), ]
        
        # return data
        customData <- customData[, c("geneID",
                               "orthoID",
                               "fullName",
                               "ncbiID",
                               "var1",
                               "var2")]
        colnames(customData) <- c("geneID", "orthoID", "fullName", "ncbiID", "FAS_F", "FAS_B")
        
        customData <- as.matrix(customData)
        return(customData)
    })

    # download data ------------------------------------------------------------
    output$downloadCustomData <- downloadHandler(
        filename = function(){
            c("customFilteredData.out")
        },
        content = function(file){
            dataOut <- downloadCustomData()
            write.table(dataOut, file,
                        sep = "\t",
                        row.names = FALSE,
                        quote = FALSE)
        }
    )

    # render download data table -----------------------------------------------
    output$filteredCustomData <- DT::renderDataTable({
        data <- downloadCustomData()
        data
    })

    # download FASTA -----------------------------------------------------------
    output$downloadCustomFasta <- downloadHandler(
        filename = function(){
            c("customFilteredSeq.fa")
        },
        content = function(file){
            fastaOutDf <- fasta()
            write.table(fastaOutDf,
                        file,
                        sep = "\t",
                        col.names = FALSE,
                        row.names = FALSE,
                        quote = FALSE)
        }
    )

    # download data as long format ---------------------------------------------
    downloadCustomDataLong <- reactive({
        downloadCustomData <- downloadCustomData()
        return(downloadCustomData[,c(1,4,2,5,6)])
    })

    output$downloadCustomLong <- downloadHandler(
        filename = function(){
            c("customFilteredData.phyloprofile")
        },
        content = function(file){
            dataOut <- downloadCustomDataLong()
            write.table(dataOut, file,
                        sep = "\t",
                        row.names = FALSE,
                        quote = FALSE)
        }
    )

    return(downloadCustomData)
}
