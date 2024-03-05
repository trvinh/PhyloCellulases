#' Import function files
sourceFiles = list.files( path = "R", pattern = "*.R$", full.names = TRUE)
lapply(sourceFiles, source, .GlobalEnv)
library(PhyloProfile)

#' set size limit for input (9999mb)
options(
    shiny.maxRequestSize = 9999 * 1024 ^ 2 # size limit for input 9999mb
)

#' MAIN SERVER =================================================================
shinyServer(function(input, output, session) {
    # Automatically stop a Shiny app when closing the browser tab
    session$allowReconnect(TRUE)

    # ========================= DOWNLOAD INPUT FILES  ==========================
    # observe({
    #     fileExist <- file.exists("data/ribi.phyloprofile")
    #     if (fileExist == FALSE) {
    #         msg <- paste0(
    #             "Please wait while phyloprofile data are being downloaded!!!"
    #         )
    #         createAlert(
    #             session, "fileExistMsgUI", "fileExistMsg", title = "",
    #             content = msg,
    #             append = FALSE
    #         )
    #         download.file(
    #             "https://applbio.biologie.uni-frankfurt.de/download/RibosomeBiogenesis/PP_RibosomeBiogenesis/FINAL_297_HsaSce.phyloprofile",
    #             destfile = "data/ribi.phyloprofile",
    #             method = "libcurl"
    #         )
    #     } else closeAlert(session, "fileExistMsg")
    # })

    observe({
        fileExist <- file.exists("data/ribi.fasta")
        if (fileExist == FALSE) {
            msg <- paste0(
                "Please wait while fasta data are being downloaded!!!"
            )
            createAlert(
                session, "fileExistMsgUI", "fileExistMsg", title = "",
                content = msg,
                append = FALSE
            )
            download.file(
                "https://applbio.biologie.uni-frankfurt.de/download/Cellulases/PhyloCellulase/rbf.fasta",
                destfile = "data/ribi.fasta",
                method = "libcurl"
            )
        } else closeAlert(session, "fileExistMsg")
    })

    observe({
        fileExist <- dir.exists("data/domains")
        if (fileExist == FALSE) {
            msg <- paste0(
                "Please wait while domain data are being downloaded!!!"
            )
            createAlert(
                session, "fileExistMsgUI", "fileExistMsg", title = "",
                content = msg,
                append = FALSE
            )
            download.file(
                "https://applbio.biologie.uni-frankfurt.de/download/Cellulases/PhyloCellulase/domains.tar.gz",
                destfile = "data/domains.tar.gz",
                method = "libcurl"
            )
            print("Extracting domain files...")
            untar("data/domains.tar.gz", exdir = "data/domains")
            file.remove("data/domains.tar.gz")
            print("Done!")
        } else closeAlert(session, "fileExistMsg")
    })

    observe({
        fileExist <- file.exists("data/filteredDfspecies.rds")
        if (fileExist == FALSE) {
            msg <- paste0(
                "Please wait while RDS data are being downloaded!!!"
            )
            createAlert(
                session, "fileExistMsgUI", "fileExistMsg", title = "",
                content = msg,
                append = FALSE
            )
            download.file(
                "https://applbio.biologie.uni-frankfurt.de/download/Cellulases/PhyloCellulase/rds.tar.gz",
                destfile = "data/rds.tar.gz",
                method = "libcurl"
            )
            print("Extracting rds files...")
            untar("data/rds.tar.gz", exdir = "data")
            file.remove("data/rds.tar.gz")
            print("Done!")
        } else closeAlert(session, "fileExistMsg")
    })

    output$warningMsg <- renderUI({
        msg <- paste(
            "<p><span style=\"color: #ff0000;\"><em><strong>Due to the large",
            "amount of data, it could take a while for rendering the plot!",
            "</strong></em></span></p>"
        )
        HTML(msg)
    })


    # ========================= INITIAL PARAMETERS  ============================
    mainInput <- "data/ribi.phyloprofile"
    var1ID <- "FAS_F"
    var2ID <- "FAS_B"
    getRefspec <- function(rankSelect) {
        ribiDf <- data.frame(
            "name" = c(
                "Saccharomyces cerevisiae S288C",
                "Saccharomyces cerevisiae",
                "Saccharomyces",
                "Saccharomycetaceae",
                "Saccharomycetales",
                "Saccharomycetes",
                "Ascomycota",
                "Fungi",
                "Eukaryota"
            ),
            "rank" = c(
                "strain",
                "species",
                "genus",
                "family",
                "order",
                "class",
                "phylum",
                "kingdom",
                "superkingdom"
            )
        )
        return(ribiDf$name[ribiDf$rank == rankSelect])
    }

    var1AggregateBy <- "max"
    var2AggregateBy <- "max"
    var1Relation <- "protein"
    var2Relation <- "protein"

    # =========================== RENDER FILTER SLIDEBARS ======================

    # * render filter slidebars for Main plot ----------------------------------
    output$var1Cutoff.ui <- renderUI({
        createSliderCutoff(
            "var1", paste(var1ID, "cutoff:"), 0.0, 1.0, var1ID
        )
    })

    output$var2Cutoff.ui <- renderUI({
        createSliderCutoff(
            "var2", paste(var2ID, "cutoff:"), 0.0, 1.0, var2ID
        )
    })

    output$percentCutoff.ui <- renderUI({
        createSliderCutoff(
            "percent", "% of present taxa:", 0.0, 1.0, "percent"
        )
    })

    # * render filter slidebars for Customized plot ----------------------------
    output$var1Filter.ui <- renderUI({
        req(input$var1)
        createSliderCutoff(
            "var1cus",
            paste(var1ID, "cutoff:"),
            input$var1[1], input$var1[2], var1ID
        )
    })

    output$var2Filter.ui <- renderUI({
        req(input$var2)
        createSliderCutoff(
            "var2cus",
            paste(var2ID, "cutoff:"),
            input$var2[1], input$var2[2], var2ID
        )
    })

    output$percentFilter.ui <- renderUI({
        req(input$percent)
        createSliderCutoff(
            "percent2",
            "% of present taxa:",
            input$percent[1], input$percent[2], "percent"
        )
    })

    output$coorthologFilter.ui <- renderUI({
        numericInput(
            "coortholog2",
            "Max co-orthologs",
            min = 1,
            max = 999,
            step = 1,
            value = input$coortholog,
            width = 150
        )
    })

    # * update value for filter slidebars of Main Plot -------------------------
    # ** based on customized profile
    observe({
        newVar1 <- input$var1cus
        updateSliderCutoff(
            session,
            "var1", paste(var1ID, "cutoff:"), newVar1, var1ID
        )
    })

    observe({
        newVar2 <- input$var2cus
        updateSliderCutoff(
            session,
            "var2", paste(var2ID, "cutoff:"), newVar2, var2ID
        )
    })

    observe({
        newPercent <- input$percent2
        updateSliderCutoff(
            session,
            "percent", "% of present taxa:", newPercent, "percent"
        )
    })

    observe({
        newCoortholog <- input$coortholog2
        updateNumericInput(
            session,
            "coortholog",
            value = newCoortholog
        )
    })

    # * reset cutoffs of Main plot ---------------------------------------------
    observeEvent(input$resetMain, {
        shinyjs::reset("var1")
        shinyjs::reset("var2")
        shinyjs::reset("percent")
        shinyjs::reset("coortholog")
    })

    # * reset cutoffs of Customized plot ---------------------------------------
    observeEvent(input$resetSelected, {
        shinyjs::reset("var1")
        shinyjs::reset("var2")
        shinyjs::reset("percent")
        shinyjs::reset("coortholog")
    })

    # ====================== PROCESSING INPUT DATA =============================
    # * convert main input file in any format into long format dataframe -------
    getMainInput <- reactive({
        withProgress(message = 'Reading main input...', value = 0.5, {
            # # inFile <- system.file(
            # #     "extdata", "ribi/ribi.phyloprofile",
            # #     package="PhyloCellulase"
            # # )
            # inFile <- "data/ribi.phyloprofile"
            # longDataframe <- createLongMatrix(inFile)
            #
            # # convert geneID, ncbiID and orthoID into factor and
            # # var1, var2 into numeric
            # for (i in seq_len(3)) {
            #     longDataframe[, i] <- as.factor(longDataframe[, i])
            # }
            # if (ncol(longDataframe) > 3) {
            #     for (j in seq(4, ncol(longDataframe))){
            #         longDataframe[,j] <- suppressWarnings(
            #             as.numeric(as.character(longDataframe[,j]))
            #         )
            #     }
            # }
            #
            # # remove duplicated lines
            # longDataframe <- longDataframe[!duplicated(longDataframe),]
            # # update number of genes to plot based on input
            # if (nlevels(as.factor(longDataframe$geneID)) <= 1500) {
            #     updateNumericInput(
            #         session,
            #         "endIndex", value = nlevels(as.factor(longDataframe$geneID))
            #     )
            # }
            # # return
            # saveRDS(longDataframe, file = "data/mainInput.rds")
            longDataframe <- readRDS("data/mainInput.rds")
            return(longDataframe)
        })
    })

    # * parse domain info into data frame --------------------------------------
    getDomainInformation <- reactive({
        withProgress(message = 'Reading domain input...', value = 0.5, {
            if (input$tabs == "Main profile") {
                # info = groupID,orthoID,supertaxon,mVar1,%spec,var2
                info <- mainpointInfo()
            } else if (input$tabs == "Customized profile") {
                info <- selectedpointInfo()
            }
            domainDf <- parseDomainInput(
                info[1],
                "data/domains/", #input$domainPath,
                "folder"
            )
            return(domainDf)
        })
    })

    # * get ID list of input taxa from main input ------------------------------
    inputTaxonID <- reactive({
        withProgress(message = 'Getting input taxon IDs...', value = 0.5, {
            longDataframe <- getMainInput()
            inputTaxa <- getInputTaxaID(longDataframe)
            return(inputTaxa)
        })
    })

    # * get NAME list of all (super)taxa ---------------------------------------
    inputTaxonName <- reactive({
        req(getMainInput())
        if (input$rankSelect == "") return()
        withProgress(message = 'Getting input taxon names...', value = 0.5, {
            inputTaxaName <- PhyloCellulase::getInputTaxaNameCr(
              input$rankSelect, inputTaxonID()
            )
            return(inputTaxaName)
        })
    })

    # * sort taxonomy data of input taxa ---------------------------------------
    sortedtaxaList <- reactive({
        # withProgress(message = 'Sorting input taxa...', value = 0.5, {
        #     # get input taxonomy tree
        #     # treeIn <- system.file(
        #     #     "extdata", "ribi/ribi.nwk",
        #     #     package="PhyloCellulase"
        #     # )
        #     # inputTaxaTree <- read.tree(file = treeIn)
        #
        #     # sort taxonomy matrix based on selected refTaxon
        #     sortedOut <- PhyloCellulase::sortInputTaxaCr(
        #         taxonIDs = inputTaxonID(),
        #         rankName = input$rankSelect,
        #         refTaxon = getRefspec(input$rankSelect),
        #         taxaTree = NULL #inputTaxaTree
        #     )
        #     # return
        #     print(getwd())
        #     saveRDS(sortedOut, file = paste0("data/sortedInput",input$rankSelect,".rds"))
        #     return(sortedOut)
        # })
        sortedOut <- readRDS(paste0("data/sortedInput",input$rankSelect,".rds"))
        return(sortedOut)
    })

    # * count taxa for each supertaxon -----------------------------------------
    getCountTaxa <- reactive({
        taxaCount <- plyr::count(sortedtaxaList(), "supertaxon")
        return(taxaCount)
    })

    # * get subset data for plotting (default 30 genes if > 50 genes) ----------
    # preData <- reactive({
    #     longDataframe <- getMainInput()
    #     req(longDataframe)
    #     # isolate start and end gene index
    #     input$updateBtn
    #     # if (input$autoUpdate == TRUE) {
    #     #     startIndex <- input$stIndex
    #     #     endIndex <- input$endIndex
    #     # } else {
    #     #     startIndex <- isolate(input$stIndex)
    #     #     endIndex <- isolate(input$endIndex)
    #     # }
    #     #
    #     # if (is.na(endIndex)) endIndex <- 1000
    #     startIndex <- 1
    #     endIndex <- nlevels(as.factor(longDataframe$geneID))
    #     withProgress(message = 'Subseting data...', value = 0.5, {
    #         longDataframe <- unsortID(longDataframe, FALSE)
    #         listIn <- input$list
    #         if (!is.null(listIn)) {
    #             list <- read.table(file = listIn$datapath, header = FALSE)
    #             listGeneOri <- list$V1
    #             if (startIndex <= length(listGeneOri)) {
    #                 listGene <- listGeneOri[listGeneOri[startIndex:endIndex]]
    #             } else listGene <- listGeneOri
    #             data <- longDataframe[longDataframe$geneID %in% listGene, ]
    #         } else {
    #             subsetID <-
    #                 levels(longDataframe$geneID)[startIndex:endIndex]
    #             data <- longDataframe[longDataframe$geneID %in% subsetID, ]
    #         }
    #
    #         if (ncol(data) < 5) {
    #             for (i in seq_len(5 - ncol(data))) {
    #                 data[paste0("newVar", i)] <- 1
    #             }
    #         }
    #
    #         # return preData
    #         if (nrow(data) == 0) return()
    #         colnames(data) <- c("geneID", "ncbiID", "orthoID", "var1", "var2")
    #         return(data)
    #     })
    # })

    # * creating main dataframe for subset taxa (in species/strain level) ------
    getFullData <- reactive({
        # req(preData())
        req(getCountTaxa())
        req(sortedtaxaList())
        # {
        #     input$plotCustom
        #     input$updateBtn
        # }
        # withProgress(message = 'Parsing profile data...', value = 0.5, {
        #     if (input$autoUpdate == TRUE) {
        #         coorthologCutoffMax <- input$coortholog
        #     } else {
        #         coorthologCutoffMax <- isolate(input$coortholog)
        #     }
        #     fullMdData <- parseInfoProfile(
        #         inputDf = preData(),
        #         sortedInputTaxa = sortedtaxaList(),
        #         taxaCount = getCountTaxa(),
        #         coorthoCOMax = coorthologCutoffMax
        #     )
        #     saveRDS(fullMdData, file = paste0("data/fullMdData",input$rankSelect,".rds"))
        #     return(fullMdData)
        # })
        fullMdData <- readRDS(paste0("data/fullMdData",input$rankSelect,".rds"))
        return(fullMdData)
    })

    # * filter full data -------------------------------------------------------
    filteredDataHeat <- reactive({
        {
            input$plotCustom
            input$updateBtn
        }
        # check input file
        filein <- mainInput
        req(filein)
        req(input$rankSelect)
        withProgress(message = 'Creating data for plotting...', value = 0.5, {
            # get all cutoffs
            if (input$autoUpdate == TRUE) {
                percentCutoff <- input$percent
                coorthologCutoffMax <- input$coortholog
                var1Cutoff <- input$var1
                var2Cutoff <- input$var2
                colorByGroup <- FALSE #input$colorByGroup
            } else {
                percentCutoff <- isolate(input$percent)
                coorthologCutoffMax <- isolate(input$coortholog)
                var1Cutoff <- isolate(input$var1)
                var2Cutoff <- isolate(input$var2)
                colorByGroup <- FALSE #isolate(input$colorByGroup)
            }

            if (
                # 4==5
                (is.null(var1Cutoff) && is.null(var2Cutoff) &&
                is.null(percentCutoff) && coorthologCutoffMax == 999) ||
                (var1Cutoff[1] == 0.0 && var1Cutoff[2] == 1.0 &&
                var2Cutoff[1] == 0.0 && var2Cutoff[2] == 1.0 &&
                percentCutoff[1] == 0.0 && percentCutoff[2] == 1.0 &&
                coorthologCutoffMax == 999)
            ) {
                filteredDf <- readRDS(paste0("data/filteredDf",input$rankSelect,".rds"))
                return(filteredDf)
            } else {
                # get selected supertaxon name
                split <- strsplit(as.character(getRefspec(input$rankSelect)), "_")
                inSelect <- as.character(split[[1]][1])

                # get gene categories
                inputCatDt <- NULL

                # create data for heatmap plotting
                filteredDf <- filterProfileData(
                    DF = getFullData(),
                    taxaCount = getCountTaxa(),
                    refTaxon = inSelect,
                    percentCutoff,
                    coorthologCutoffMax,
                    var1Cutoff,
                    var2Cutoff,
                    var1Relation,
                    var2Relation,
                    groupByCat = colorByGroup,
                    catDt = inputCatDt,
                    var1AggregateBy = var1AggregateBy,
                    var2AggregateBy = var2AggregateBy
                )
                # saveRDS(filteredDf, file = paste0("data/filteredDf",input$rankSelect,".rds"))
                return(filteredDf)
            }
        })
    })

    # * heatmap data input -----------------------------------------------------
    dataHeat <- reactive({
        req(filteredDataHeat())
        dataHeat <- PhyloCellulase::reduceProfileCr(filteredDataHeat())
        return(dataHeat)
    })

    # =========================== MAIN PROFILE TAB =============================
    # ** warning message -------------------------------------------------------
    observe({
        desc = paste(
            "<p>Rendering the plot may take some time due to large amount of data!</p>",
            "<p>For small monitors, the preset plot size (8800x5000) may throw",
            "an error. In that case, please reduce the plot size and press the",
            "\'Update plot\' button.</p>"
        )

        if (input$tabs == "Main profile") {
            createAlert(
                session, "warningUI", "warningPlot",
                title = "", content = desc, append = FALSE
            )
        }
    })

    # * get total number of genes ----------------------------------------------
    output$totalGeneNumber.ui <- renderUI({
        geneList <- getMainInput()
        out <- as.list(levels(factor(geneList$geneID)))
        listIn <- input$list
        if (!is.null(listIn)) {
            list <- read.table(file = listIn$datapath, header = FALSE)
            out <- as.list(list$V1)
        }
        if (length(out) > 0) {
            strong(paste0("Total number of genes:  ", length(out)))
        }
    })

    # * get list of taxa for highlighting --------------------------------------
    output$highlightTaxonUI <- renderUI({
        choice <- inputTaxonName()
        out <- as.list(levels(factor(choice$fullName)))
        out <- append("none", out)

        selectInput("taxonHighlight", "Select (super)taxon to highlight:",
                    out, selected = out[1])
    })

    # * get list of genes for highlighting -------------------------------------
    output$highlightGeneUI <- renderUI({
        geneList <- dataHeat()
        out <- as.list(levels(factor(geneList$geneID)))
        out <- append("none", out)
        selectInput(
            "geneHighlight", "Select gene to highlight:", out, selected = out[1]
        )
    })

    # * update plot size based on input ----------------------------------------
    observe({
        longDataframe <- getMainInput()
        req(longDataframe)
        req(input$rankSelect)
        if (input$autoSizing) {
            inputSuperTaxon <- inputTaxonName()
            nrTaxa <- nlevels(as.factor(inputSuperTaxon$fullName))
            nrGene <- nlevels(as.factor(longDataframe$geneID))
            # adapte to axis type
            if (input$xAxis == "taxa") {
                h <- nrGene
                w <- nrTaxa
            } else {
                w <- nrGene
                h <- nrTaxa
            }
            # adapt to dot zoom factor
            if (input$dotZoom < -0.5){
                hv <- (200 + 12 * h) * (1 + input$dotZoom) + 500
                wv <- (200 + 12 * w) * (1 + input$dotZoom) + 500
            }  else if ((input$dotZoom < 0)) {
                hv <- (200 + 12 * h) * (1 + input$dotZoom) + 200
                wv <- (200 + 12 * w) * (1 + input$dotZoom) + 200
            } else {
                hv <- (200 + 12 * h) * (1 + input$dotZoom)
                wv <- (200 + 12 * w) * (1 + input$dotZoom)
            }
            # minimum size
            if (hv < 300) hv <- 300
            if (wv < 300) wv <- 300
            # update plot size based on number of genes/taxa
            hv <- hv + 300
            wv <- wv + 300

            if (input$rankSelect == "species") {
                if (input$xAxis == "taxa") {
                    if (hv > 4136) hv <- 4136
                    if (wv > 21000) wv <- 21000
                } else {
                    if (wv > 4136) wv <- 4136
                    if (hv > 21000) hv <- 21000
                }
            }

            if (h <= 20) {
                updateSelectInput(
                    session, "mainLegend",
                    label = "Legend position:",
                    choices = list("Right" = "right",
                                   "Left" = "left",
                                   "Top" = "top",
                                   "Bottom" = "bottom",
                                   "Hide" = "none"),
                    selected = "top"
                )
                updateNumericInput(
                    session,
                    "width", value = wv  + 50
                )
            } else if (h <= 30) {
                updateNumericInput(
                    session,
                    "width", value = wv + 50
                )
            } else {
                updateNumericInput(
                    session,
                    "width", value = wv
                )
            }
            updateNumericInput(
                session,
                "height", value = hv
            )
        }
    })

    # observe({
    #     if (input$rankSelect == "species") {
    #         updateSliderInput(
    #             session,
    #             "dotZoom", "",
    #             min = -1,
    #             max = 3,
    #             step = 0.1,
    #             value = -0.2
    #         )
    #         updateNumericInput(
    #             session,
    #             "xSize",
    #             "X-axis label size (px)",
    #             min = 3,
    #             max = 99,
    #             step = 1,
    #             value = 12
    #         )
    #         updateNumericInput(
    #             session,
    #             "ySize",
    #             "Y-axis label size (px)",
    #             min = 3,
    #             max = 99,
    #             step = 1,
    #             value = 12
    #         )
    #     }
    # })

    # * reset configuration windows of Main plot -------------------------------
    observeEvent(input$resetMainConfig, {
        shinyjs::reset("xSize")
        shinyjs::reset("ySize")
        shinyjs::reset("legendSize")
        shinyjs::reset("xAngle")
        shinyjs::reset("dotZoom")
    })

    # * close configuration windows of Main plot -------------------------------
    observeEvent(input$applyMainConfig, {
        toggleModal(session, "mainPlotConfigBs", toggle = "close")
    })

    # * parameters for the main profile plot -----------------------------------
    getParameterInputMain <- reactive({
        input$updateBtn
        if (input$autoUpdate == TRUE) {
            inputPara <- list(
                "xAxis" = input$xAxis,
                "var1ID" = var1ID,
                "var2ID"  = var2ID,
                "midVar1" = input$midVar1,
                "midVar2" = input$midVar2,
                "lowColorVar1" =  input$lowColorVar1,
                "midColorVar1" =  input$midColorVar1,
                "highColorVar1" = input$highColorVar1,
                "lowColorVar2" = input$lowColorVar2,
                "midColorVar2" =  input$midColorVar2,
                "highColorVar2" = input$highColorVar2,
                "paraColor" = input$paraColor,
                "xSize" = input$xSize,
                "ySize" = input$ySize,
                "legendSize" = input$legendSize,
                "mainLegend" = input$mainLegend,
                "dotZoom" = input$dotZoom,
                "xAngle" = input$xAngle,
                "guideline" = 1,
                "width" = input$width,
                "height" = input$height,
                "colorByGroup" = FALSE #input$colorByGroup
            )
        } else {
            inputPara <- isolate(
                list(
                    "xAxis" = input$xAxis,
                    "var1ID" = var1ID,
                    "var2ID"  = var2ID,
                    "midVar1" = input$midVar1,
                    "midVar2" = input$midVar2,
                    "lowColorVar1" =  input$lowColorVar1,
                    "midColorVar1" =  input$midColorVar1,
                    "highColorVar1" = input$highColorVar1,
                    "lowColorVar2" = input$lowColorVar2,
                    "midColorVar2" =  input$midColorVar2,
                    "highColorVar2" = input$highColorVar2,
                    "paraColor" = input$paraColor,
                    "xSize" = input$xSize,
                    "ySize" = input$ySize,
                    "legendSize" = input$legendSize,
                    "mainLegend" = input$mainLegend,
                    "dotZoom" = input$dotZoom,
                    "xAngle" = input$xAngle,
                    "guideline" = 1,
                    "width" = input$width,
                    "height" = input$height,
                    "colorByGroup" = FALSE #input$colorByGroup
                )
            )
        }
        return(inputPara)
    })

    # * render dot size to dotSizeInfo ---------------------------------------
    output$dotSizeInfo <- renderUI({
        dataHeat <- dataHeat()
        dataHeat$presSpec[dataHeat$presSpec == 0] <- NA
        presentVl <- dataHeat$presSpec[!is.na(dataHeat$presSpec)]

        minDot <- (floor(min(presentVl) * 10) / 10 * 5) * (1 + input$dotZoom)
        maxDot <- (floor(max(presentVl) * 10) / 10 * 5) * (1 + input$dotZoom)

        em(paste0("current point's size: ", minDot, " - ", maxDot))
    })

    # * plot main profile ------------------------------------------------------
    mainpointInfo <- callModule(
        createProfilePlot, "mainProfile",
        data = dataHeat,
        # clusteredDataHeat = clusteredDataHeat,
        # applyCluster = reactive(input$applyCluster),
        parameters = getParameterInputMain,
        inSeq = reactive(input$inSeq),
        inTaxa = reactive(input$inTaxa),
        rankSelect = reactive(input$rankSelect),
        inSelect = reactive(getRefspec(input$rankSelect)),
        taxonHighlight = reactive(input$taxonHighlight),
        geneHighlight = reactive(input$geneHighlight),
        typeProfile = reactive("mainProfile")
    )

    # ======================== CUSTOMIZED PROFILE TAB ==========================

    # * get list of all sequence IDs for customized profile -----
    output$geneIn <- renderUI({
        filein <- mainInput
        fileCustom <- input$customFile
        data <- getFullData()
        outAll <- c("all", as.list(levels(factor(data$geneID))))
        if (!is.null(fileCustom)) {
            customList <- read.table(
                file = fileCustom$datapath, header = FALSE
            )
            customList$V1 <- as.factor(customList$V1)
            outAll <- as.list(levels(customList$V1))
        }
        if (outAll[1] == "all") {
            createSelectGene("inSeq", outAll, "all")
        } else {
            createSelectGene("inSeq", outAll, outAll)
        }
    })

    # * get list of all taxa for customized profile ----------------------------
    output$taxaIn <- renderUI({
        filein <- mainInput
        if (is.null(filein)) return(selectInput("inTaxa", "", "all"))
        choice <- inputTaxonName()
        out <- c("all", as.list(levels(factor(choice$fullName))))
        selectInput("inTaxa", "",
                    out,
                    selected = out[1],
                    multiple = TRUE,
                    selectize = FALSE)
    })

    # * check if all genes and all species are selected ------------------------
    output$sameProfile <- reactive({
        if (length(input$inSeq[1]) == 0) return(FALSE)
        else {
            if (input$inSeq[1] == "all" & input$inTaxa[1] == "all") return(TRUE)
        }
    })
    outputOptions(output, "sameProfile", suspendWhenHidden = FALSE)

    # * update customized plot size based on input -----------------------------
    observe({
        longDataframe <- getMainInput()
        req(longDataframe)
        req(input$inTaxa)
        req(input$inSeq)
        if (input$selectedAutoSizing) {
            nrTaxa <- length(input$inTaxa)
            nrGene <- length(input$inSeq)
            if (input$inTaxa[1] == "all") {
                inputSuperTaxon <- inputTaxonName()
                nrTaxa <- nlevels(as.factor(inputSuperTaxon$fullName))
            }
            if (input$inSeq[1] == "all") {
                nrGene <- nlevels(as.factor(longDataframe$geneID))
            }
            # adapte to axis type
            if (input$xAxisSelected == "taxa") {
                h <- nrGene
                w <- nrTaxa
            } else {
                w <- nrGene
                h <- nrTaxa
            }
            # adapt to dot zoom factor
            if (input$dotZoomSelect < -0.5){
                hv <- (200 + 12 * h) * (1 + input$dotZoomSelect) + 500
                wv <- (200 + 12 * w) * (1 + input$dotZoomSelect) + 500
            }  else if ((input$dotZoomSelect < 0)) {
                hv <- (200 + 12 * h) * (1 + input$dotZoomSelect) + 200
                wv <- (200 + 12 * w) * (1 + input$dotZoomSelect) + 200
            } else {
                hv <- (200 + 12 * h) * (1 + input$dotZoomSelect)
                wv <- (200 + 12 * w) * (1 + input$dotZoomSelect)
            }
            # minimum size
            if (hv < 300) hv <- 300
            if (wv < 300) wv <- 300
            # update plot size based on number of genes/taxa
            hv <- hv + 300
            wv <- wv + 300
            if (h <= 20) {
                updateSelectInput(
                    session, "selectedLegend",
                    label = "Legend position:",
                    choices = list("Right" = "right",
                                   "Left" = "left",
                                   "Top" = "top",
                                   "Bottom" = "bottom",
                                   "Hide" = "none"),
                    selected = "top"
                )
                updateNumericInput(
                    session,
                    "selectedWidth", value = wv  + 50
                )
            } else if (h <= 30) {
                updateNumericInput(
                    session,
                    "selectedWidth", value = wv + 50
                )
            } else {
                updateNumericInput(
                    session,
                    "selectedWidth", value = wv
                )
            }
            updateNumericInput(
                session,
                "selectedHeight", value = hv
            )
        }
    })

    # * reset configuration windows of Customized plot -------------------------
    observeEvent(input$resetSelectedConfig, {
        shinyjs::reset("xSizeSelect")
        shinyjs::reset("ySizeSelect")
        shinyjs::reset("legendSizeSelect")
        shinyjs::reset("xAngleSelect")
        shinyjs::reset("dotZoomSelect")
    })

    # ** close configuration windows of Customized plot ------------------------
    observeEvent(input$applySelectedConfig, {
        toggleModal(session, "selectedPlotConfigBs", toggle = "close")
    })

    # * parameters for the customized profile plot -----------------------------
    getParameterInputCustomized <- reactive({
        input$plotCustom
        if (input$autoUpdateSelected == TRUE) {
            inputPara <- list(
                "xAxis" = input$xAxisSelected,
                "var1ID" = var1ID,
                "var2ID"  = var2ID,
                "midVar1" = input$midVar1,
                "midVar2" = input$midVar2,
                "lowColorVar1" =  input$lowColorVar1,
                "midColorVar1" =  input$midColorVar1,
                "highColorVar1" = input$highColorVar1,
                "lowColorVar2" = input$lowColorVar2,
                "midColorVar2" =  input$midColorVar2,
                "highColorVar2" = input$highColorVar2,
                "paraColor" = input$paraColor,
                "xSize" = input$xSizeSelect,
                "ySize" = input$ySizeSelect,
                "legendSize" = input$legendSizeSelect,
                "mainLegend" = input$selectedLegend,
                "dotZoom" = input$dotZoomSelect,
                "xAngle" = input$xAngleSelect,
                "guideline" = 0,
                "width" = input$selectedWidth,
                "height" = input$selectedHeight,
                "colorByGroup" = FALSE #input$colorByGroup
            )
        } else {
            inputPara <- isolate(
                list(
                    "xAxis" = input$xAxisSelected,
                    "var1ID" = var1ID,
                    "var2ID"  = var2ID,
                    "midVar1" = input$midVar1,
                    "midVar2" = input$midVar2,
                    "lowColorVar1" =  input$lowColorVar1,
                    "midColorVar1" =  input$midColorVar1,
                    "highColorVar1" = input$highColorVar1,
                    "lowColorVar2" = input$lowColorVar2,
                    "midColorVar2" =  input$midColorVar2,
                    "highColorVar2" = input$highColorVar2,
                    "paraColor" = input$paraColor,
                    "xSize" = input$xSizeSelect,
                    "ySize" = input$ySizeSelect,
                    "legendSize" = input$legendSizeSelect,
                    "mainLegend" = input$selectedLegend,
                    "dotZoom" = input$dotZoomSelect,
                    "xAngle" = input$xAngleSelect,
                    "guideline" = 0,
                    "width" = input$selectedWidth,
                    "height" = input$selectedHeight,
                    "colorByGroup" = FALSE #input$colorByGroup
                )
            )
        }
        return(inputPara)
    })

    # * plot customized profile ------------------------------------------------
    selectedpointInfo <- callModule(
        createProfilePlot, "customizedProfile",
        data = dataHeat,
        # clusteredDataHeat = clusteredDataHeat,
        # applyCluster = reactive(input$applyCluster),
        parameters = getParameterInputCustomized,
        inSeq = reactive(input$inSeq),
        inTaxa = reactive(input$inTaxa),
        rankSelect = reactive(input$rankSelect),
        inSelect = reactive(getRefspec(input$rankSelect)),
        taxonHighlight = reactive("none"),
        geneHighlight = reactive("none"),
        typeProfile = reactive("customizedProfile")
    )

    # ============================== POINT INFO ================================

    # * get status of pointInfo for activating Detailed Plot button -----------
    output$pointInfoStatus <- reactive({
        if (input$tabs == "Main profile") {
            # info contains groupID,orthoID,supertaxon,mVar1,%spec,var2
            info <- mainpointInfo()
        } else if (input$tabs == "Customized profile") {
            info <- selectedpointInfo()
        } else info <- NULL
        return(is.null(info))
    })
    outputOptions(output, "pointInfoStatus", suspendWhenHidden = FALSE)

    # * show info into "point's info" box --------------------------------------
    output$pointInfo <- renderText({
        # GET INFO BASED ON CURRENT TAB
        if (input$tabs == "Main profile") {
            # info contains groupID,orthoID,supertaxon,mVar1,%spec,var2
            info <- mainpointInfo()
        } else if (input$tabs == "Customized profile") {
            info <- selectedpointInfo()
        } else return()

        req(info)
        orthoID <- info[2]

        if (is.na(orthoID)) return()
        else {
            a <- toString(paste("Seed-ID:", info[1]))
            b <- toString(paste0(
                "Hit-ID: ", orthoID,
                " (", info[3], ")"
            ))
            c <- ""
            if (var1ID != "") {
                c <- toString(paste(
                    var1AggregateBy, var1ID, ":", info[4]
                ))
            }
            d <- ""
            if (var2ID != "") {
                d <- toString(paste(
                    var2AggregateBy, var2ID, ":", info[6]
                ))
            }
            e <- toString(paste("% present taxa:", info[5]))
            paste(a, b, c, d, e, sep = "\n")
        }
    })

    # ============================= DETAILED PLOT ==============================
    # * data for detailed plot -------------------------------------------------
    detailPlotDt <- reactive({
        # GET INFO BASED ON CURRENT TAB
        if (input$tabs == "Main profile") {
            # info contains groupID,orthoID,supertaxon,mVar1,%spec,var2
            info <- mainpointInfo()
        } else if (input$tabs == "Customized profile") {
            info <- selectedpointInfo()
        }

        req(info)
        withProgress(message = 'Getting data for detailed plot...', value=0.5, {
            ### get refspec name
            split <- strsplit(as.character(getRefspec(input$rankSelect)), "_")
            inSelect <- as.character(split[[1]][1])

            ### get info for present taxa in selected supertaxon (1)
            fullDf <- getFullData()
            ### filter data if needed
            if  (input$detailedFilter == TRUE) {
                fullDf <- filteredDataHeat()
                if (info[3] == inSelect) {
                    fullDf <- fullDf[
                        fullDf$var1 >= input$var1[1]
                        & fullDf$var1 <= input$var1[2],
                    ]
                    fullDf <- fullDf[
                        fullDf$var2 >= input$var2[1]
                        & fullDf$var2 <= input$var2[2],
                    ]
                }
                updateCheckboxInput(
                    session, "detailedRemoveNA", value = TRUE
                )
            }
            plotTaxon <- unique(
                fullDf$supertaxon[grep(info[3], fullDf$supertaxon)]
            )
            plotGeneID <- info[1]
            selDf <- fullDf[fullDf$geneID == plotGeneID
                            & fullDf$supertaxon == plotTaxon, ]
            ### get all taxa of this supertaxon (2)
            allTaxaDf <- sortedtaxaList()
            allTaxaDf <- allTaxaDf[allTaxaDf$supertaxon == plotTaxon,
                                   c("abbrName", "fullName")]

            ### merge (1) and (2) together
            joinedDf <- merge(selDf, allTaxaDf, by = c("abbrName"), all.y =TRUE)
            joinedDf <- subset(
                joinedDf,
                select = c(
                    "abbrName", "fullName.y", "geneID", "orthoID", "var1","var2"
                )
            )
            names(joinedDf)[names(joinedDf) == "fullName.y"] <- "fullName"

            # replace var1/var2 as NA for all "NA orthologs"
            joinedDf$var1[is.na(joinedDf$orthoID)] <- NA
            joinedDf$var2[is.na(joinedDf$orthoID)] <- NA

            # remove NA orthologs if required
            if (input$detailedRemoveNA == TRUE) {
                joinedDf <- joinedDf[!is.na(joinedDf$orthoID), ]
            }

            ### return data for detailed plot
            return(joinedDf)
        })
    })

    # * render detailed plot ---------------------------------------------------
    pointInfoDetail <- callModule(
        createDetailedPlot, "detailedPlot",
        data = detailPlotDt,
        var1ID = reactive(var1ID),
        var2ID = reactive(var2ID),
        detailedText = reactive(input$detailedText),
        detailedHeight = reactive(input$detailedHeight)
    )

    # * render database links --------------------------------------------------
    output$dbLink <- renderUI({
        info <- pointInfoDetail() # info = seedID, orthoID, var1
        req(info)
        seqID <- toString(info[2])
        tmp <- as.list(strsplit(seqID, "\\|")[[1]])
        linkText <- ""
        # get taxon ID
        taxon <- tmp[[2]]
        taxId <- as.list(strsplit(taxon, "@")[[1]])[[2]]
        taxUrl <- paste0(
            "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=",
            taxId
        )
        linkText <- paste0(
            "<p><a href='", taxUrl, "' target='_blank'>",
            "NCBI taxonomy entry for <strong>", taxId , "</strong></a></p>"
        )
        # get protein ID
        protId <- tmp[[3]]
        uniprotUrl <- paste0("https://www.uniprot.org/uniprot/", protId)
        ncbiUrl <- paste0("https://www.ncbi.nlm.nih.gov/protein/", protId)
        if (RCurl::url.exists(uniprotUrl)) {
            linkText <- paste0(
                linkText, "<p><a href='", uniprotUrl, "' target='_blank'>",
                "UniProt entry for <strong>", protId, "</strong></a></p>"
            )
        } else if (RCurl::url.exists(ncbiUrl)) {
            linkText <- paste0(
                linkText, "<p><a href='", ncbiUrl, "' target='_blank'>",
                "NCBI protein entry for <strong>", protId, "</strong></a></p>"
            )
        }
        # render links
        linkText <- paste0(
            linkText,
            "<p><em><strong>Disclaimer:</strong> ",
            "External links are automatically generated and may point to ",
            "a wrong target (see <a ",
            "href=\"https://github.com/BIONF/PhyloProfile/wiki/FAQ",
            "#wrong-info-from-public-databases\" ",
            "target=\"_blank\">FAQ</a>)</em></p>"
        )
        HTML(linkText)
    })

    # * render FASTA sequence --------------------------------------------------
    output$fasta <- renderText({
        info <- pointInfoDetail() # info = seedID, orthoID, var1
        req(info)
        seqID <- toString(info[2])
        fastain <- "data/ribi.fasta"
        # fastain <- system.file(
        #     "extdata", "ribi/ribi.fasta",
        #     package="PhyloCellulase"
        # )
        fastaOut <- getFastaFromFile(seqID, fastain)
        return(paste(fastaOut[1]))
    })

    # ======================== FEATURE ARCHITECTURE PLOT =======================
    # * render domain plot -----------------------------------------------------
    observeEvent(input$doDomainPlot, {
        callModule(
            createArchitecturePlot, "archiPlot",
            pointInfo = pointInfoDetail,
            domainInfo = getDomainInformation,
            labelArchiSize = reactive(input$labelArchiSize),
            titleArchiSize = reactive(input$titleArchiSize),
            archiHeight = reactive(input$archiHeight),
            archiWidth = reactive(input$archiWidth)
        )
    })

    # ======================== FILTERED DATA DOWNLOADING =======================

    # * for main profile =======================================================
    mainFastaDownload <- reactive({
        downloadDf <- as.data.frame(downloadData())
        seqIDs <- downloadDf$orthoID
        fastain <- "data/ribi.fasta"
        # fastain <- system.file(
        #     "extdata", "ribi/ribi.fasta",
        #     package="PhyloCellulase"
        # )
        mainFastaOut <- getFastaFromFile(seqIDs, fastain)
        return(mainFastaOut)
    })

    downloadData <- callModule(
        downloadFilteredMain,
        "filteredMainDownload",
        data = getFullData,
        taxaCount = getCountTaxa,
        fasta = mainFastaDownload,
        var1ID = reactive(var1ID),
        var2ID = reactive(var2ID),
        var1 = reactive(input$var1),
        var2 = reactive(input$var2),
        percent = reactive(input$percent)
    )

    # * for customized profile =================================================
    customizedFastaDownload <- reactive({
        downloadDf <- as.data.frame(downloadCustomData())
        seqIDs <- downloadDf$orthoID
        fastain <- "data/ribi.fasta"
        # fastain <- system.file(
        #     "extdata", "ribi/ribi.fasta",
        #     package="PhyloCellulase"
        # )
        fastaOutDf <- getFastaFromFile(seqIDs, fastain)
        return(fastaOutDf)
    })

    downloadCustomData <- callModule(
        downloadFilteredCustomized,
        "filteredCustomizedDownload",
        data = downloadData,
        fasta = customizedFastaDownload,
        inSeq = reactive(input$inSeq),
        inTaxa = reactive(input$inTaxa)
    )
})
