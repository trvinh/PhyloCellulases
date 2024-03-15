#' Import function files
sourceFiles = list.files( path = "R", pattern = "*.R$", full.names = TRUE)
lapply(sourceFiles, source, .GlobalEnv)
library(PhyloProfile)

#' set size limit for input (9999mb)
options(
    scipen = 999,
    shiny.maxRequestSize = 9999 * 1024 ^ 2 # size limit for input 9999mb
)

#' MAIN SERVER =================================================================
shinyServer(function(input, output, session) {
    homePath = c(wd='~/')
    # Automatically stop a Shiny app when closing the browser tab
    session$allowReconnect(TRUE)

    nameFullFile <- paste0(
        getwd(), "/data/preProcessedTaxonomy.txt"
    )
    currentNCBIinfo <- NULL
    if (file.exists(nameFullFile))
        currentNCBIinfo <- as.data.frame(data.table::fread(nameFullFile))

    # =========================== INITIAL CHECKING  ============================
    # * check for internet connection ------------------------------------------
    observe({
        if (hasInternet() == FALSE) toggleState("demoData")
    })

    output$noInternetMsg <- renderUI({
        if (hasInternet() == FALSE && input$demoData != "preCalcDt") {
            strong(
                em("Internet connection is required for using demo data!"),
                style = "color:red"
            )
        } else return()
    })
    
    # ========================== PRE-DEFINED VARIABLES =========================
    fasta_file <- "data/cellulase_all.extended.fa"
    
    var1ID <- "FAS_F"
    var2ID <- "FAS_B"
    seedSource <- "ncbi"
    orthoSource <- "ncbi"
    var1AggregateBy <- "max"
    var2AggregateBy <- "max"
    var1Relation <- "protein"
    var2Relation <- "protein"
    seqIdFormat <- 1 # bionf format
    separator <- 1 # |
    applyCluster <- TRUE
    clusterMethod <- "complete"
    taxDBpath <- "/Users/vinh/projects/fdog_manuscript/cellulase/new/pp_rhiso"
    
    # * reset profile plot colors ----------------------------------------------
    observeEvent(input$defaultColorVar2, {
        shinyjs::reset("lowColorVar2")
        shinyjs::reset("midColorVar2")
        shinyjs::reset("highColorVar2")
        shinyjs::reset("midVar2")
    })

    observeEvent(input$defaultColorVar1, {
        shinyjs::reset("lowColorVar1")
        shinyjs::reset("midColorVar1")
        shinyjs::reset("highColorVar1")
        shinyjs::reset("midVar1")
    })

    observeEvent(input$defaultColorPara, {
        shinyjs::reset("paraColor")
    })

    # * render list of taxonomy ranks ------------------------------------------
    output$rankSelect <- renderUI({
        selectInput(
            "rankSelect", label = h5("Select taxonomy rank:"),
            choices = getTaxonomyRanks(),
            selected = "class"
        )
    })

    # # * render list of (super)taxa ---------------------------------------------
    # observe({
    #     choice <- inputTaxonName()
    #     choice$fullName <- as.factor(choice$fullName)
    # 
    #     if (input$demoData == "arthropoda") {
    #         hellemDf <- data.frame(
    #             "name" = c(
    #                 "Drosophila melanogaster",
    #                 "Drosophila melanogaster",
    #                 "Drosophila",
    #                 "Drosophilidae",
    #                 "Diptera",
    #                 "Insecta",
    #                 "Arthropoda",
    #                 "Metazoa",
    #                 "Eukaryota"
    #             ),
    #             "rank" = c(
    #                 "strain",
    #                 "species",
    #                 "genus",
    #                 "family",
    #                 "order",
    #                 "class",
    #                 "phylum",
    #                 "kingdom",
    #                 "superkingdom"
    #             )
    #         )
    #         rankName <- input$rankSelect
    # 
    #         updateSelectizeInput(
    #             session, "inSelect", "", server = TRUE,
    #             choices = as.list(levels(choice$fullName)),
    #             selected = hellemDf$name[hellemDf$rank == rankName]
    #         )
    #     } else if (input$demoData == "ampk-tor") {
    #         humanDf <- data.frame(
    #             "name" = c(
    #                 "Homo sapiens",
    #                 "Homo sapiens",
    #                 "Homo",
    #                 "Hominidae",
    #                 "Primates",
    #                 "Mammalia",
    #                 "Chordata",
    #                 "Metazoa",
    #                 "Eukaryota"
    #             ),
    #             "rank" = c(
    #                 "strain",
    #                 "species",
    #                 "genus",
    #                 "family",
    #                 "order",
    #                 "class",
    #                 "phylum",
    #                 "kingdom",
    #                 "superkingdom"
    #             )
    #         )
    #         rankName <- input$rankSelect
    # 
    #         updateSelectizeInput(
    #             session, "inSelect", "", server = TRUE,
    #             choices = as.list(levels(choice$fullName)),
    #             selected = humanDf$name[humanDf$rank == rankName]
    #         )
    #     } else if (input$demoData == "preCalcDt") {
    #         refspec <- levels(choice$fullName)[1]
    #         if (!is.null(i_refspec)) {
    #             if (i_refspec %in% levels(choice$fullName))
    #                 refspec <- i_refspec
    #         }
    #         updateSelectizeInput(
    #             session, "inSelect", "", server = TRUE,
    #             choices = as.list(levels(choice$fullName)),
    #             selected = refspec
    #         )
    #     } else {
    #         if (length(choice$fullName) > 0) {
    #             updateSelectizeInput(
    #                 session, "inSelect", "", server = TRUE,
    #                 choices = as.list(levels(choice$fullName)),
    #                 selected = levels(choice$fullName)[1]
    #             )
    #         }
    #     }
    # })

    # * enable/disable update plot button --------------------------------------
    observe({
        if (input$autoUpdate == TRUE) {
            updateButton(session, "updateBtn", disabled = TRUE)
        } else {
            updateButton(session, "updateBtn", disabled = FALSE)
        }
    })

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
        createSliderCutoff(
            "var1cus", paste(var1ID, "cutoff:"), 0.0, 1.0, var1ID
        )
    })

    output$var2Filter.ui <- renderUI({
        createSliderCutoff(
            "var2cus", paste(var2ID, "cutoff:"), 0.0, 1.0, var2ID
        )
    })

    output$percentFilter.ui <- renderUI({
        createSliderCutoff(
            "percent2", "% of present taxa:", 0.0, 1.0, "percent"
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
    
    # * reset cutoffs of Main plot ---------------------------------------------
    observeEvent(input$resetMain, {
        shinyjs::reset("var1")
        shinyjs::reset("var2")
        shinyjs::reset("percent")
        shinyjs::reset("coortholog")
    })

    # * reset cutoffs of Customized plot ---------------------------------------
    observeEvent(input$resetSelected, {
        shinyjs::reset("var1cus")
        shinyjs::reset("var2cus")
        shinyjs::reset("percent2")
        shinyjs::reset("coortholog2")
    })

    # ====================== PROCESSING INPUT DATA =============================

    # * convert main input file in any format into long format dataframe -------
    getMainInput <- reactive({
        withProgress(message = 'Reading main input...', value = 0.5, {
            longDataframe <- readRDS("data/mainInput.rds")
            return(longDataframe)
        })
    })

    # * parse domain info into data frame --------------------------------------
    getDomainInformation <- reactive({
        withProgress(message = 'Reading domain input...', value = 0.5, {
            if (input$tabs == "Main profile") {
                # info contains groupID, orthoID,
                # supertaxon, mVar1, %spec, var2
                info <- mainpointInfo()
            } else if (input$tabs == "Customized profile") {
                info <- selectedpointInfo()
            }
            req(info[1])
            domainDf <- parseDomainInput(
                info[1],
                "data/domains/",
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
        })
    })

    # * get NAME list of all (super)taxa ---------------------------------------
    inputTaxonName <- reactive({
        req(input$rankSelect)
        req(getMainInput())
        if (input$rankSelect == "") return()
        withProgress(message = 'Getting input taxon names...', value = 0.5, {
            inputTaxaName <- getInputTaxaName(
                input$rankSelect, inputTaxonID(), taxDBpath
            )
            return(inputTaxaName)
        })
    })

    # * sort taxonomy data of input taxa ---------------------------------------
    sortedtaxaList <- reactive({
        req(input$rankSelect)
        sortedOut <- readRDS(paste0("data/sortedInput",input$rankSelect,".rds"))
        return(sortedOut)
    })

    # * count taxa for each supertaxon -----------------------------------------
    getCountTaxa <- reactive({
        taxaCount <- plyr::count(sortedtaxaList(), "supertaxon")
        return(taxaCount)
    })

    # * creating main dataframe for subset taxa (in species/strain level) ------
    getFullData <- reactive({
        req(getCountTaxa())
        req(sortedtaxaList())
        fullMdData <- readRDS(paste0("data/fullMdData",input$rankSelect,".rds"))
    })

    # * filter full data -------------------------------------------------------
    filteredDataHeat <- reactive({
        req(input$rankSelect)
        {
            input$plotCustom
            input$updateBtn
        }
        
        # check input file
        # filein <- input$mainInput
        # if (
        #     input$demoData == "arthropoda" | input$demoData == "ampk-tor" |
        #     input$demoData == "preCalcDt"
        # ) {
        #     filein <- 1
        # }
        # req(filein)
        withProgress(message = 'Creating data for plotting...', value = 0.5, {
            # get all cutoffs
            if (input$autoUpdate == TRUE) {
                percentCutoff <- input$percent
                coorthologCutoffMax <- input$coortholog
                var1Cutoff <- input$var1
                var2Cutoff <- input$var2
                colorByGroup <- input$colorByGroup
            } else {
                percentCutoff <- isolate(input$percent)
                coorthologCutoffMax <- isolate(input$coortholog)
                var1Cutoff <- isolate(input$var1)
                var2Cutoff <- isolate(input$var2)
                colorByGroup <- isolate(input$colorByGroup)
            }

            # get selected supertaxon name
            split <- strsplit(as.character(input$inSelect), "_")
            inSelect <- as.character(split[[1]][1])

            # # get gene categories
            # inputCatDt <- NULL
            # if (length(colorByGroup) > 0 && colorByGroup == TRUE) {
            #     # get gene category
            #     geneCategoryFile <- input$geneCategory
            #     if (!is.null(geneCategoryFile)) {
            #         inputCatDt <- read.table(
            #             file = geneCategoryFile$datapath,
            #             sep = "\t",
            #             header = FALSE,
            #             check.names = FALSE,
            #             comment.char = "",
            #             fill = TRUE
            #         )
            #         colnames(inputCatDt) <- c("geneID","group")
            #     } else if (!is.null(i_geneCategory)){
            #         inputCatDt <- read.table(
            #             file = i_geneCategory,
            #             sep = "\t",
            #             header = FALSE,
            #             check.names = FALSE,
            #             comment.char = "",
            #             fill = TRUE
            #         )
            #         colnames(inputCatDt) <- c("geneID","group")
            #     } else inputCatDt <- NULL
            # } else colorByGroup = FALSE
            # 
            # # create data for heatmap plotting
            # filteredDf <- filterProfileData(
            #     DF = getFullData(),
            #     taxaCount = getCountTaxa(),
            #     refTaxon = inSelect,
            #     percentCutoff,
            #     coorthologCutoffMax,
            #     var1Cutoff,
            #     var2Cutoff,
            #     var1Relation,
            #     var2Relation,
            #     groupByCat = colorByGroup,
            #     catDt = inputCatDt,
            #     var1AggregateBy = var1AggregateBy,
            #     var2AggregateBy = var2AggregateBy
            # )
            # saveRDS(filteredDf, file = paste0("data/filteredDf",input$rankSelect,".rds"))
            filteredDf <- readRDS(paste0("data/filteredDf",input$rankSelect,".rds"))
            ### NEED TO APPLY FILTERED HERERE!!!!
            return(filteredDf)
        })
    })

    # * heatmap data input -----------------------------------------------------
    dataHeat <- reactive({
        req(filteredDataHeat())
        dataHeat <- reduceProfile(filteredDataHeat())
        return(dataHeat)
    })

    # * clustered heatmap data -------------------------------------------------
    clusteredDataHeat <- reactive({
        dataHeat <- dataHeat()

        if (nlevels(as.factor(dataHeat$geneID)) > 1) {
            withProgress(message = 'Clustering profile data...', value = 0.5, {
                dat <- readRDS(paste0("data/data4cluster",input$rankSelect,".rds"))
                # do clustering based on distance matrix
                row.order <- hclust(
                    getDistanceMatrixProfiles(), method = clusterMethod
                )$order

                # re-order distance matrix accoring to clustering
                datNew <- dat[row.order, ] #col.order

                # return clustered gene ID list
                clusteredGeneIDs <- as.factor(row.names(datNew))

                # sort original data according to clusteredGeneIDs
                dataHeat$geneID <- factor(dataHeat$geneID, levels=clusteredGeneIDs)

                dataHeat <- dataHeat[!is.na(dataHeat$geneID),]
                return(dataHeat)
            })
        } else return(dataHeat)
        return(dataHeat)
    })

    # * get list of all input (super)taxa and their ncbi IDs -------------------
    allInputTaxa <- reactive({
        return()
    })

    # =========================== MAIN PROFILE TAB =============================

    # * render popup for selecting rank and return list of subset taxa ---------
    mainTaxaName <- callModule(
        selectTaxonRank,
        "selectTaxonRankMain",
        rankSelect = reactive(input$rankSelect),
        inputTaxonID = inputTaxonID,
        taxDB = reactive(taxDBpath)
    )

    # * get list of taxa for highlighting --------------------------------------
    output$taxonHighlight.ui <- renderUI({
        choice <- inputTaxonName()
        out <- levels(factor(choice$fullName))
        if (input$applyMainTaxa == TRUE) {
            out <- mainTaxaName()
            selectizeInput(
                "taxonHighlight", "", out, selected = out, multiple = TRUE,
                options = list(placeholder = 'none')
            )
        } else {
            selectizeInput(
                "taxonHighlight", "", out, multiple = TRUE,
                options = list(placeholder = 'none')
            )
        }
    })

    # * get list of genes for highlighting -------------------------------------
   observe({
        geneList <- dataHeat()
        out <- levels(factor(geneList$geneID))

        if (!(is.null(input$geneHighlightFile))) {
            fileHighlight <- input$geneHighlightFile
            highlightList <- read.table(
                file = fileHighlight$datapath, header = FALSE
            )
            updateSelectizeInput(
                session, "geneHighlight", server = TRUE,
                choices = out, selected = intersect(out, highlightList$V1)
            )
        } else {
            updateSelectizeInput(
                session, "geneHighlight", server = TRUE,
                choices = out
            )
        }
    })

    # * disable/enable highlighing orthologs having the same ID ----------------
    observe({
        longDataframe <- getMainInput()
        req(longDataframe)
        req(input$rankSelect)
        lowestRank <- getLowestRank(longDataframe, taxDBpath)
        if (!(lowestRank == input$rankSelect))
            shinyjs::disable("colorByOrthoID")
        else
            shinyjs::enable("colorByOrthoID")
    })

    # * render list of superRanks for adding vertical lines --------------------
    output$superRankSelect.ui <- renderUI({
        allRanks <- getTaxonomyRanks()
        selectInput(
            "superRankSelect", label = "Display taxonomic labels for:",
            choices = c(allRanks[!(allRanks %in% c("strain"))]),
            selected = ""
        )
    })

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

        colorByGroup <- FALSE # input$colorByGroup
        # get category colors
        catColors <- NULL

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
                "guideline" = 0,
                "width" = input$width,
                "height" = input$height,
                "colorByGroup" = colorByGroup,
                "catColors" = catColors,
                "colorByOrthoID" = input$colorByOrthoID,
                "groupLabelSize" = input$groupLabelSize,
                "groupLabelDist" = input$groupLabelDist,
                "groupLabelAngle" = input$groupLabelAngle
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
                    "guideline" = 0,
                    "width" = input$width,
                    "height" = input$height,
                    "colorByGroup" = colorByGroup,
                    "catColors" = catColors,
                    "colorByOrthoID" = input$colorByOrthoID,
                    "groupLabelSize" = input$groupLabelSize,
                    "groupLabelDist" = input$groupLabelDist,
                    "groupLabelAngle" = input$groupLabelAngle
                )
            )
        }
        return(inputPara)
    })

    # * render dot size to dotSizeInfo ---------------------------------------
    output$dotSizeInfo <- renderUI({
        # req(v$doPlot)

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
        clusteredDataHeat = clusteredDataHeat,
        applyCluster = reactive(applyCluster),
        parameters = getParameterInputMain,
        inSeq = reactive(input$inSeq),
        inTaxa = reactive(input$inTaxa),
        rankSelect = reactive(input$rankSelect),
        inSelect = reactive(input$inSelect),
        taxonHighlight = reactive(input$taxonHighlight),
        geneHighlight = reactive(input$geneHighlight),
        typeProfile = reactive("mainProfile"),
        taxDB = reactive(taxDBpath),
        superRank = reactive(input$superRankSelect),
        allTaxa = allInputTaxa
    )

    # ======================== CUSTOMIZED PROFILE TAB ==========================

    # * get list of all sequence IDs for customized profile -----
    output$cusGene.ui <- renderUI({
        # filein <- input$mainInput
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
            selectizeInput(
                "inSeq", "", outAll, selected = "all", multiple = TRUE
            )
        } else {
            selectizeInput(
                "inSeq", "", outAll, selected = outAll, multiple = TRUE
            )
        }
    })

    # * render popup for selecting rank and return list of subset taxa ---------
    output$cusRankSelect.ui <- renderUI({
        selectInput(
            "cusRankSelect", label = h5("Select taxonomy rank:"),
            choices = append("none", getTaxonomyRanks()),
            selected = "none"
        )
    })
    
    # get (super)taxa based on selected rank
    cusSuperTaxa <- reactive({
        rankSelectCus <- input$cusRankSelect
        
        if (length(rankSelectCus) == 0 || rankSelectCus == "none") return()
        else {
            # load list of unsorted taxa
            Dt <- getTaxonomyMatrix(taxDBpath, TRUE, inputTaxonID())
            # load list of taxon name
            nameList <- getNameList(taxDBpath)
            # get rank name from rankSelect
            rankName <- rankSelectCus
            
            choice <- data.frame(
                ncbiID = unlist(Dt[rankName]), stringsAsFactors = FALSE
            )
            choice <- merge(choice, nameList, by = "ncbiID", all = FALSE)
            return(choice)
        }
    })
    
    # render list of possible taxon names from getTaxaCus()
    output$cusSelectTaxonRank.ui <- renderUI({
        req(input$cusRankSelect)
        if (input$cusRankSelect == "none") return()
        choice <- cusSuperTaxa()
        choice$fullName <- as.factor(choice$fullName)
        selectInput(
            "cusTaxonSelect",
            label = h5("Choose (super)taxon of interest:"),
            as.list(levels(choice$fullName)),
            levels(choice$fullName)[1]
        )
    })

    # * get list of all taxa for customized profile ----------------------------
    cusTaxa <- reactive({
        req(input$cusRankSelect)
        taxaSelectCus <- input$cusTaxonSelect
        rankName <- input$cusRankSelect
        if (taxaSelectCus == "") return()
        
        # load list of unsorted taxa
        Dt <- getTaxonomyMatrix(taxDBpath, TRUE, inputTaxonID())
        
        # get ID of selected (super)taxon from input$taxaSelectCus
        taxaList <- getNameList(taxDBpath)
        superID <- taxaList$ncbiID[taxaList$fullName == taxaSelectCus
                                   & taxaList$rank %in% c(rankName, "norank")]

        # from that ID, get list of all taxa for main selected taxon
        mainRankName <- "strain" #rankSelect()
        customizedtaxaID <-
            levels(as.factor(Dt[mainRankName][Dt[rankName] == superID, ]))
        
        cusTaxaName <-
            taxaList$fullName[taxaList$ncbiID %in% customizedtaxaID]
        return(cusTaxaName)
    })
    
    output$cusTaxa.ui <- renderUI({
        input$plotCustom
        out <- isolate(cusTaxa())
        selectizeInput("inTaxa","Selected taxa",out, selected = out, multiple = TRUE)
    })

    # * render list of superRanks for adding vertical lines --------------------
    output$cusSuperRankSelect.ui <- renderUI({
        allRanks <- getTaxonomyRanks()
        selectInput(
            "cusSuperRankSelect", label = "Display taxonomic labels for:",
            choices = c(allRanks[!(allRanks %in% c("strain"))]),
            selected = ""
        )
    })

    # * check if all genes and all species are selected ------------------------
    output$sameProfile <- reactive({
        # if (v$doPlot == FALSE) return(FALSE)
        if (length(input$inSeq[1]) == 0) return(FALSE)
        else {
            if (length(input$inSeq) == 0 || length(input$inTaxa) == 0)
                return(TRUE)
            if ("all" %in% input$inSeq & "all" %in% input$inTaxa) return(TRUE)
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
            if (nrTaxa < 10000 && nrGene < 10000) {
                if ("all" %in% input$inTaxa) {
                    inputSuperTaxon <- inputTaxonName()
                    nrTaxa <- nlevels(as.factor(inputSuperTaxon$fullName))
                }
                if ("all" %in% input$inSeq) {
                    nrGene <- 320 # input$endIndex
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

        colorByGroup <- FALSE #input$colorByGroup
        # get category colors
        catColors <- NULL

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
                "colorByGroup" = colorByGroup,
                "catColors" = catColors,
                "colorByOrthoID" = input$colorByOrthoID,
                "groupLabelSize" = input$groupLabelSizeSelect,
                "groupLabelDist" = input$groupLabelDistSelect,
                "groupLabelAngle" = input$groupLabelAngleSelect
            )
        )
        return(inputPara)
    })

    # * plot customized profile ------------------------------------------------
    cusFilteredDataHeat <- reactive ({
        filteredDf <- readRDS("data/filteredDfstrain.rds")
        ### NEED TO APPLY FILTERED HERERE!!!!
        return(filteredDf)
    })
    cusDataHeat <- reactive({
        return(reduceProfile(cusFilteredDataHeat()))
    })
    cusClusteredDataHeat <- reactive({
        cusDataHeat <- cusDataHeat()
        if (nlevels(as.factor(cusDataHeat$geneID)) > 1) {
            withProgress(message = 'Clustering profile data...', value = 0.5, {
                # dat <- getProfiles()
                # saveRDS(dat, file = paste0("data/data4cluster",input$rankSelect,".rds"))
                dat <- readRDS(paste0("data/data4cluster",input$rankSelect,".rds"))
                # do clustering based on distance matrix
                row.order <- hclust(
                    getDistanceMatrixProfiles(), method = clusterMethod
                )$order
                
                # re-order distance matrix accoring to clustering
                datNew <- dat[row.order, ] #col.order
                
                # return clustered gene ID list
                clusteredGeneIDs <- as.factor(row.names(datNew))
                
                # sort original data according to clusteredGeneIDs
                cusDataHeat$geneID <- factor(cusDataHeat$geneID, levels=clusteredGeneIDs)
                
                cusDataHeat <- cusDataHeat[!is.na(cusDataHeat$geneID),]
                return(cusDataHeat)
            })
        } else return(cusDataHeat)
    })
    
    selectedpointInfo <- callModule(
        createProfilePlot, "customizedProfile",
        data = cusDataHeat, #dataHeat,
        clusteredDataHeat = cusClusteredDataHeat, #clusteredDataHeat,
        applyCluster = reactive(applyCluster),
        parameters = getParameterInputCustomized,
        inSeq = reactive(input$inSeq),
        inTaxa = reactive(input$inTaxa),
        rankSelect = reactive(input$rankSelect),
        inSelect = reactive(input$inSelect),
        taxonHighlight = reactive("none"),
        geneHighlight = reactive("none"),
        typeProfile = reactive("customizedProfile"),
        taxDB = reactive(taxDBpath),
        superRank = reactive(input$cusSuperRankSelect),
        allTaxa = allInputTaxa
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
        orthoID <- info[[2]]
        if (length(info[[2]]) > 1) orthoID <- paste0(info[[2]][1], ",...")

        if (is.na(orthoID)) return()
        else {
            a <- toString(paste("Seed-ID:", info[[1]]))
            b <- toString(paste0(
                "Hit-ID: ", orthoID
            ))
            c <- ""
            if (var1ID != "") {
                c <- toString(paste(
                    var1AggregateBy, var1ID, ":", info[[5]]
                ))
            }
            d <- ""
            if (var2ID != "") {
                d <- toString(paste(
                    var2AggregateBy, var2ID, ":", info[[7]]
                ))
            }

            if (info[[10]] == "Y") {
                if (info[[3]] == 1) {
                    s <- toString(
                        paste0(info[[3]], " ortholog in ", info[[4]], ":")
                    )
                } else {
                    s <- toString(
                        paste0(info[[3]], " co-orthologs in ", info[[4]], ":")
                    )
                }
                e <- ""
            } else {
                s <- toString(paste0("Best ortholog in ", info[[4]], ":"))
                e <- toString(
                    paste0(
                        "% present taxa: ", info[[6]], " (", info[[8]], " out of ",
                        info[[9]],  ")", collapse = ""
                    )
                )
            }
            paste(a, s, b, c, d, e, sep = "\n")
        }
    })

    # ============================= DETAILED PLOT ==============================
    # * data for detailed plot -------------------------------------------------
    detailPlotDt <- reactive({
        # req(v$doPlot)
        req(input$detailedBtn)
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
            split <- strsplit(as.character(input$inSelect), "_")
            inSelect <- as.character(split[[1]][1])

            ### get info for present taxa in selected supertaxon (1)
            fullDf <- getFullData()
            if (input$tabs == "Customized profile") {
                fullDf <- readRDS(paste0("data/fullMdDatastrain.rds"))
            }
            ### filter data if needed
            if  (input$detailedFilter == TRUE) {
                fullDf <- filteredDataHeat()
                if (info[[4]] == inSelect) {
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
                fullDf$supertaxon[grep(paste0(info[[4]],"$"), fullDf$supertaxon)]
            )
            plotGeneID <- info[[1]]
            selDf <- fullDf[fullDf$geneID == plotGeneID
                            & fullDf$supertaxon == plotTaxon, ]
            ### get all taxa of this supertaxon (2)
            allTaxaDf <- sortedtaxaList()
            if (input$tabs == "Customized profile") {
                allTaxaDf <- readRDS(paste0("data/sortedInputstrain.rds"))
            }
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
    parseProId <- function(protId, separator, seqIdFormat) {
        tmp <- ""
        if (separator == 1) {
            if (grepl("\\|", protId))
                tmp <- as.list(strsplit(protId, "\\|")[[1]])
        } else if (separator == 2) {
            if (grepl("@", protId))
                tmp <- as.list(strsplit(protId, "@")[[1]])
        } else if (separator == 3) {
            if (grepl("#", protId))
                tmp <- as.list(strsplit(protId, "#")[[1]])
        } else if (separator == 4) {
            if (grepl(";", protId))
                tmp <- as.list(strsplit(protId, ";")[[1]])
        }
        if (length(tmp) > 1) {
            if (seqIdFormat == 3) {
                protId <- tmp[[length(tmp)]]
            } else if (seqIdFormat == 4) {
                protId <- tmp[[1]]
            } else if (seqIdFormat == 1) {
                if (length(tmp) >= 3) {
                    protId <- tmp[[3]]
                } else {
                    warning("Wrong ID format was set!")
                    return(NULL)
                }
            }
        }
        return(protId)
    }

    output$dbLink.ui <- renderUI({
        info <- pointInfoDetail() # info = seedID, orthoID, var1, var2, ncbiID
        req(info)
        linkText <- ""
        # get seed ID
        seedId <- toString(info[1])
        separators <- c("|", "@", "#", ";")
        if (grepl(paste(separators, collapse = "|"), seedId)) {
            seedId <- parseProId(seedId, separator, seqIdFormat)
        }
        if (length(seedId) > 0) {
            if (length(strsplit(seedId, "_")[[1]]) > 3) {
                seedIdMod <- paste0("XP_", gsub("^.*_", "", seedId))
            } else {
                seedIdMod <- gsub("^.*_", "", seedId)
            }
            linkText <- paste0(linkText, createDBlink(seedIdMod, "NCBI"))
        }
        # get ortho ID
        protId <- toString(info[2])
        if (grepl(paste(separators, collapse = "|"), protId)) {
            protId <- parseProId(protId, separator, seqIdFormat)
        }
        if (length(protId) > 0) {
            linkText <- paste0(linkText, createDBlink(protId, "NCBI"))
        }
        # get taxon ID
        taxId <- gsub("ncbi", "", info[5])
        taxHierarchy <- PhyloProfile:::getTaxHierarchy(taxId, currentNCBIinfo)
        taxUrls <- paste(taxHierarchy[[1]]$link, collapse = ", ")
        taxUrls <- gsub("<p>", "", taxUrls)
        taxUrls <- gsub("</p>", "", taxUrls)
        linkText <- paste0(
            linkText, "<p><strong>NCBI taxonomy: </strong>", taxUrls, "</p>"
        )

        # render links
        linkText <- paste0(
            linkText,
            "<p><em><strong>Disclaimer:</strong> ",
            "External links are automatically generated and may point to ",
            "a wrong target. Please adapt the sequence ID format according ",
            "to your data in the Input and Settings tab (see <a ",
            "href=\"https://github.com/BIONF/PhyloProfile/wiki/FAQ",
            "#wrong-info-from-public-databases\" ",
            "target=\"_blank\">FAQ</a>)</em></p>"
        )
        HTML(linkText)
    })

    # * render FASTA sequence --------------------------------------------------
    output$fasta <- renderText({
        # req(v$doPlot)
        info <- pointInfoDetail() # info = seedID, orthoID, var1
        req(info)
        seqID <- toString(info[2])
        fastaOut <- getFastaFromFile(seqID, fasta_file)
        return(paste(fastaOut[1]))
    })

    # ======================== FEATURE ARCHITECTURE PLOT =======================
    # * check if seed and orthoID are specified and get domain file/path -------
    observe({
        infoTmp <- c()
        if (input$tabs == "Main profile") {
            infoTmp <- mainpointInfo()
        } else if (input$tabs == "Customized profile") {
            infoTmp <-selectedpointInfo()
        }
        if (!is.null(infoTmp)) {
            tmp <- getDomainFile()
        }
    })

    getDomainFile <- reactive({
        # get lowest rank
        # activate doDomainPlotMain if either working on the lowest rank
        # or only 1 ortholog present
        longDataframe <- getMainInput()
        req(longDataframe)
        req(input$rankSelect)
        lowestRank <- getLowestRank(longDataframe, taxDBpath)

        # get info from POINT INFO box
        info <- c()
        infoTmp <- c()
        if (input$tabs == "Main profile") {
            # info contains groupID,orthoID,supertaxon,mVar1,%spec,var2
            infoTmp <- mainpointInfo()
        } else if (input$tabs == "Customized profile") {
            infoTmp <- selectedpointInfo()
        }
        # only when no co-ortholog exists
        if (length(infoTmp[[2]]) == 1) {
            info <- c(infoTmp[[1]], as.character(infoTmp[[2]]))
        } else {
            # else, get info from detailed plot
            updateButton(session, "doDomainPlotMain", disabled = TRUE)
            if (!is.null(pointInfoDetail())) {
                info <- pointInfoDetail() # info = seedID, orthoID, var1
            }
        }

        if (is.null(info)) {
            updateButton(session, "doDomainPlot", disabled = TRUE)
            updateButton(session, "doDomainPlotMain", disabled = TRUE)
            return("noSelectHit")
        } else {
            req(info)
            domainDf <- parseDomainInput(
                info[1], "data/domains/", "folder"
            )
            if (length(domainDf) == 1) {
                if (domainDf == "noSelectHit" |
                    domainDf == "noFileInFolder") {
                    updateButton(
                        session, "doDomainPlot", disabled = TRUE
                    )
                    updateButton(
                        session, "doDomainPlotMain", disabled = TRUE
                    )
                    return(domainDf)
                } else {
                    updateButton(
                        session, "doDomainPlot", disabled = FALSE
                    )
                    if (lowestRank == input$rankSelect || infoTmp[[8]][1] == 1)
                        updateButton(session, "doDomainPlotMain", disabled = FALSE)
                    else
                        updateButton(session, "doDomainPlotMain", disabled = TRUE)
                }
            } else {
                updateButton(
                    session, "doDomainPlot", disabled = FALSE
                )
                if (lowestRank == input$rankSelect || infoTmp[[8]][1] == 1)
                    updateButton(session, "doDomainPlotMain", disabled = FALSE)
                else
                    updateButton(session, "doDomainPlotMain", disabled = TRUE)
            }
        }
        return(info)
    })

    # * check domain file ------------------------------------------------------
    output$checkDomainFiles <- renderUI({
        fileDomain <- getDomainFile()
        if (length(fileDomain) == 1) {
            if (fileDomain == "noFileInput") {
                em("Domain file not provided!!")
            } else if (fileDomain == "noFileInFolder") {
                msg <- paste0(
                    "<p><em>Domain file not found!! </em></p>
                <p><em>Please make sure that file name has to be in this format:
                <strong>&lt;seedID&gt;.extension</strong>, where extension is
                limited to <strong>txt</strong>, <strong>csv</strong>,
                <strong>list</strong>, <strong>domains</strong> or
                <strong>architecture</strong>.</em></p>"
                )
                HTML(msg)
            } else if (fileDomain == "noSelectHit") {
                em("Please select one ortholog sequence!!")
            }
        }
    })

    # * render domain plot -----------------------------------------------------
    observeEvent(input$doDomainPlot, {
        callModule(
            createArchitecturePlot, "archiPlot",
            pointInfo = getDomainFile,
            domainInfo = getDomainInformation,
            labelArchiSize = reactive(input$labelArchiSize),
            titleArchiSize = reactive(input$titleArchiSize),
            archiHeight = reactive(input$archiHeight),
            archiWidth = reactive(input$archiWidth),
            seqIdFormat = reactive(seqIdFormat),
            currentNCBIinfo = reactive(currentNCBIinfo)
        )
    })
    observeEvent(input$doDomainPlotMain, {
        callModule(
            createArchitecturePlot, "archiPlotMain",
            pointInfo = getDomainFile,
            domainInfo = getDomainInformation,
            labelArchiSize = reactive(input$labelArchiSizeM),
            titleArchiSize = reactive(input$titleArchiSizeM),
            archiHeight = reactive(input$archiHeightM),
            archiWidth = reactive(input$archiWidthM),
            seqIdFormat = reactive(seqIdFormat),
            currentNCBIinfo = reactive(currentNCBIinfo)
        )
    })

    # ======================== FILTERED DATA DOWNLOADING =======================

    # * for main profile =======================================================
    mainFastaDownload <- reactive({
        downloadDf <- as.data.frame(downloadData())
        seqIDs <- downloadDf$orthoID

        if (input$demoData == "none") {
            filein <- input$mainInput
            inputType <- checkInputValidity(filein$datapath)
            # get fata from oma
            if (inputType == "oma") {
                allOmaDf <- finalOmaDf()
                filteredDownloadDf <- as.data.frame(downloadData())
                filteredOmaDf <-
                    subset(allOmaDf,
                           allOmaDf$orthoID %in% filteredDownloadDf$orthoID &
                               allOmaDf$seed %in% filteredDownloadDf$geneID)
                mainFastaOut <- getAllFastaOma(filteredOmaDf)
            }
            # get fasta from main input
            else if (inputType == "fasta") {
                seqIDMod <- paste0(
                    as.character(downloadDf$geneID), "|",
                    as.character(downloadDf$ncbiID), "|",
                    as.character(downloadDf$orthoID)
                )
                mainFastaOut <- getFastaFromFasInput(
                    seqIDMod, file = filein$datapath
                )
            } else {
                # get from concaternated file
                if (input$inputType == "Concatenated fasta file") {
                    fastain <- input$concatFasta
                    mainFastaOut <- getFastaFromFile(seqIDs, fastain$datapath)
                }
                # get from folder
                else {
                    mainFastaOut <- getFastaFromFolder(
                        seqIDs,
                        input$path,
                        input$dirFormat,
                        input$fileExt,
                        input$idFormat
                    )
                }
            }
        } else {
            if (input$demoData == "preCalcDt") {
                if (!is.null(i_fastaInput))
                    mainFastaOut <- getFastaFromFile(seqIDs, i_fastaInput)
            } else {
                # get fasta from demo online data
                mainFastaOut <- getFastaDemo(seqIDs, demoData = input$demoData)
            }
        }
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

        if (input$demoData == "none") {
            filein <- input$mainInput
            inputType <- checkInputValidity(filein$datapath)
            # get fata from oma
            if (inputType == "oma") {
                allOmaDf <- finalOmaDf()
                filteredDownloadDf <- as.data.frame(downloadCustomData())
                filteredOmaDf <-
                    subset(allOmaDf,
                           allOmaDf$orthoID %in% filteredDownloadDf$orthoID &
                               allOmaDf$seed %in% filteredDownloadDf$geneID)
                fastaOutDf <- getAllFastaOma(filteredOmaDf)
            }
            # get fasta from main input
            else if (inputType == "fasta") {
                seqIDMod <- paste0(
                    as.character(downloadDf$geneID), "|",
                    as.character(downloadDf$ncbiID), "|",
                    as.character(downloadDf$orthoID)
                )
                fastaOutDf <- getFastaFromFasInput(
                    seqIDMod, file = filein$datapath
                )
            } else {
                # get from concaternated file
                if (input$inputType == "Concatenated fasta file") {
                    fastain <- input$concatFasta
                    fastaOutDf <- getFastaFromFile(seqIDs, fastain$datapath)
                }
                # get from folder
                else {
                    fastaOutDf <- getFastaFromFolder(
                        seqIDs,
                        input$path,
                        input$dirFormat,
                        input$fileExt,
                        input$idFormat
                    )
                }
            }
        } else {
            if (input$demoData == "preCalcDt") {
                if (!is.null(i_fastaInput))
                    fastaOutDf <- getFastaFromFile(seqIDs, i_fastaInput)
            } else {
                # get fasta from demo online data
                fastaOutDf <- getFastaDemo(seqIDs, demoData = input$demoData)
            }
        }
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

    # ============================ ANALYSIS FUNCTIONS ==========================

    # * PROFILE CLUSTERING =====================================================

    # ** calculate distance matrix ---------------------------------------------
    getDistanceMatrixProfiles <- reactive({
        withProgress(message = 'Calculating distance matrix...', value = 0.5, {
            # req(getProfiles())
            # # if (is.null(input$distMethod))
            # #     distMethod <- i_distMethod
            # # else distMethod <- input$distMethod
            # distMethod <- "euclidean"
            # distanceMatrix <- getDistanceMatrix(getProfiles(), distMethod)
            # saveRDS(distanceMatrix, file = paste0("data/distanceMatrix",input$rankSelect,".rds"))
            distanceMatrix <- readRDS(paste0("data/distanceMatrix",input$rankSelect,".rds"))
            return(distanceMatrix)
        })
    })

    # * DISTRIBUTION ANALYSIS ==================================================
    # ** description for distribution analysis function ------------------------
    observe({
        desc = paste(
            "Plot the distributions of the values incurred by the integrated
            information layers."
        )

        if (input$tabs == "Distribution analysis") {
            createAlert(
                session, "descDistributionUI", "descDistribution",
                content = desc, append = FALSE
            )
        }
    })

    # ** list of available variables for distribution plot ---------------------
    output$selected.distribution <- renderUI({
        if (nchar(var1ID) == 0 & nchar(var2ID) == 0) {
            varList <- "% present taxa"
        } else if (nchar(var1ID) == 0 & nchar(var2ID) > 0) {
            varList <- as.list(c(var2ID, "% present taxa"))
        } else if (nchar(var1ID) > 0 & nchar(var2ID) == 0) {
            varList <- as.list(c(var1ID, "% present taxa"))
        } else {
            varList <- as.list(c(var1ID, var2ID, "% present taxa"))
        }

        selectInput(
            "selectedDist", "Choose variable to plot:", varList, varList[1]
        )
    })

    # ** var1 / var2 distribution data -----------------------------------------
    distributionDf <- reactive({
        # req(v$doPlot)
        withProgress(message = 'Getting data for analyzing...', value = 0.5, {
            splitDt <- createVariableDistributionData(
                getMainInput(), input$var1, input$var2
            )
            # filter data base on customized plot (if chosen)
            if (input$dataset.distribution == "Customized data") {
                req(input$inSeq)
                splitDt <- createVariableDistributionDataSubset(
                    getFullData(),
                    splitDt,
                    input$inSeq,
                    input$inTaxa
                )
            }
            # return dt
            return(splitDt)
        })
    })

    # ** render distribution plots ---------------------------------------------
    observe({
        # req(v$doPlot)
        req(input$selectedDist)

        if (input$selectedDist == "% present taxa") {
            callModule(
                analyzeDistribution, "distPlot",
                data = reactive(
                    createPercentageDistributionData(
                        getMainInput(), input$rankSelect, taxDBpath
                    )
                ),
                varID = reactive(input$selectedDist),
                varType = reactive("presSpec"),
                percent = reactive(input$percent),
                distTextSize = reactive(input$distTextSize),
                distWidth = reactive(input$distWidth)
            )
        } else {
            if (input$selectedDist == var1ID) {
                callModule(
                    analyzeDistribution, "distPlot",
                    data = distributionDf,
                    varID = reactive(input$selectedDist),
                    varType = reactive("var1"),
                    percent = reactive(input$percent),
                    distTextSize = reactive(input$distTextSize),
                    distWidth = reactive(input$distWidth)
                )
            } else if (input$selectedDist == var2ID) {
                callModule(
                    analyzeDistribution, "distPlot",
                    data = distributionDf,
                    varID = reactive(input$selectedDist),
                    varType = reactive("var2"),
                    percent = reactive(input$percent),
                    distTextSize = reactive(input$distTextSize),
                    distWidth = reactive(input$distWidth)
                )
            }
        }
    })
})
