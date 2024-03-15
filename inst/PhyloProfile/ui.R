#' Import function files
sourceFiles = list.files(path = "R", pattern = "*.R$", full.names = TRUE)
lapply(sourceFiles, source, .GlobalEnv)

#' MAIN UI ====================================================================
shinyUI(
    fluidPage(
        includeCSS("www/custom.css"),
        tags$style(type = "text/css", "body {padding-top: 80px;}"),
        useShinyjs(),

        # Application title
        titlePanel("", windowTitle = "PhyloProfile Cellulases"),

        # TOP WELLPANEL FOR PLOT CONFIGURATION ---------------------------------
        conditionalPanel(
            condition = "input.tabs=='Main profile'",
            wellPanel(
                fluidRow(
                    column(
                        2,
                        radioButtons(
                            inputId = "xAxis",
                            label = "Choose type of x-axis:",
                            choices = list("taxa", "genes"),
                            selected = "taxa",
                            inline = TRUE
                        ),
                        hr(),
                        checkboxInput(
                            "autoUpdate",
                            strong(em("Auto update plot")),
                            value = FALSE,
                            width = NULL
                        ),
                        checkboxInput(
                            "keepOrder",
                            strong(em("Retain gene order")),
                            value = TRUE,
                            width = NULL
                        ),
                        bsPopover(
                            "keepOrder",
                            "",
                            "Do no change gene order while filtering data",
                            "bottom"
                        )
                    ),
                    column(
                        1,
                        createPlotSize("width", "Width (px)", 6000)
                        # ,
                        # checkboxInput(
                        #     "autoSizing",
                        #     strong(em("Auto sizing")),
                        #     value = FALSE,
                        #     width = NULL
                        # )
                    ),
                    column(
                        1, createPlotSize("height", "Height (px)", 3000),
                        actionButton("mainPlotConfig", "Appearance")
                    ),
                    column(
                        2, uiOutput("var1Cutoff.ui")
                    ),
                    column(
                        2, uiOutput("var2Cutoff.ui")
                    ),
                    column(
                        2, uiOutput("percentCutoff.ui")
                    ),
                    column(
                        2,
                        numericInput(
                            "coortholog",
                            "Max co-orthologs",
                            min = 1,
                            max = 999,
                            step = 1,
                            value = 999,
                            width = 150
                        ),
                        bsButton(
                            "resetMain",
                            "Reset cutoffs",
                            style = "danger",
                            icon = icon("backward")
                        )
                    )
                )
            )
        ),

        conditionalPanel(
            condition = "input.tabs=='Customized profile'",
            wellPanel(
                fluidRow(
                    column(
                        2,
                        radioButtons(
                            inputId = "xAxisSelected",
                            label = "Choose type of x-axis:",
                            choices = list("taxa", "genes"),
                            selected = "taxa",
                            inline = TRUE
                        ),
                        hr(),
                    ),
                    column(
                        1,
                        createPlotSize("selectedWidth", "Width (px)", 6000),
                        checkboxInput(
                            "selectedAutoSizing",
                            strong(em("Auto sizing")),
                            value = FALSE,
                            width = NULL
                        )
                    ),
                    column(
                        1, createPlotSize("selectedHeight", "Height (px)", 3000),
                        actionButton("selectedPlotConfig", "Appearance")
                    ),
                    column(
                        2, uiOutput("var1Filter.ui")
                    ),
                    column(
                        2, uiOutput("var2Filter.ui")
                    ),
                    column(
                        2, uiOutput("percentFilter.ui")
                    ),
                    column(
                        2,
                        uiOutput("coorthologFilter.ui"),
                        bsButton(
                            "resetSelected",
                            "Reset cutoffs",
                            style = "danger",
                            icon = icon("backward")
                        )
                    )
                )
            )
        ),

        # MAIN NARVARPAGE TABS -------------------------------------------------
        navbarPage(
            em(strong("PhyloProfile Cellulases")),
            id = "tabs",
            collapsible = TRUE,
            inverse = TRUE,
            fluid = TRUE,
            position = "fixed-top",

            # MAIN PROFILE TAB =================================================
            tabPanel(
                "Main profile",
                sidebarLayout(
                    # * sidebar panel for profile highlight --------------------
                    sidebarPanel(
                        width = 4,
                        uiOutput("rankSelect"),
                        hr(),
                        column(
                            12,
                            style = "padding:0px;",
                            strong("Select gene to highlight:")
                        ),
                        column(
                            12,
                            fluidRow(
                                column(
                                    8,
                                    style = "padding:0px;",
                                    selectizeInput(
                                        "geneHighlight","", NULL, multiple=TRUE, 
                                        options = list(placeholder = 'none')
                                    ),
                                    bsPopover(
                                        "geneHighlight",
                                        "",
                                        "Select gene to highlight",
                                        "right"
                                    )
                                ),
                                column(
                                    4,
                                    fileInput(
                                        "geneHighlightFile", "", width = "100%"
                                    )
                                )
                            )
                        ),
                        column(
                            12,
                            style = "padding:0px;",
                            strong(
                                "Select (super)taxon to highlight:"
                            )
                        ),
                        column(
                            8,
                            style = "padding:0px;",
                            uiOutput("taxonHighlight.ui")
                        ),
                        column(
                            4,
                            h3(""),
                            bsButton("taxonHighlightBrowse", "Browse...")
                        ),
                        column(
                            12,
                            # checkboxInput(
                            #     "colorByGroup",
                            #     strong("Highlight genes by categories"),
                            #     value = FALSE
                            # ),
                            checkboxInput(
                                "colorByOrthoID",
                                strong("Highlight duplicated ortholog IDs"),
                                value = FALSE
                            ),
                            bsPopover(
                                "colorByOrthoID",
                                "",
                                paste(
                                    "Please check in the Clustering profiles",
                                    "function, if the profiles are clustered",
                                    "using ortho IDs"
                                ),
                                "bottom"
                            ),
                            hr()
                        ),
                        uiOutput("superRankSelect.ui"),
                        hr(),
                        bsButton(
                            "updateBtn", "Update plot", style = "warning",
                            icon("sync"), disabled = FALSE
                        )
                    ),
                    # * main panel for profile plot ----------------------------
                    mainPanel(
                        createProfilePlotUI("mainProfile")
                        # conditionalPanel(
                        #     condition = "input.do > 0",
                        #     createProfilePlotUI("mainProfile")
                        # )
                    )
                )
            ),

            # CUSTOMIZED PROFILE TAB ===========================================
            tabPanel(
                "Customized profile",
                sidebarLayout(
                    # * sidebar panel for subseting data -----------------------
                    sidebarPanel(
                        width = 4,
                        column(
                            12,
                            style = "padding:0px;",
                            strong(
                                "Select (super)taxon/(super)taxa of interest:"
                            )
                        ),
                        column(
                            12,
                            style = "padding:0px;",
                            uiOutput("cusRankSelect.ui"),
                            # selectTaxonRankUI("selectTaxonRank")
                        ),
                        column(
                            12,
                            style = "padding:0px;",
                            uiOutput("cusSelectTaxonRank.ui"),
                            strong(h5("Choose (super)taxon of interest:")),
                            selectizeInput(
                                "inSelect", "", choices = NULL, selected = NULL
                            )
                        ),
                        column(
                            12,
                            style = "padding:0px;",
                            uiOutput("cusTaxa.ui")
                        ),
                        column(
                            12,
                            style = "padding:0px;",
                            strong("Select sequence(s) of interest:")
                        ),
                        
                        column(
                            12,
                            fluidRow(
                                column(
                                    8,
                                    style = "padding:0px;",
                                    uiOutput("cusGene.ui")
                                ),
                                column(
                                    4,
                                    fileInput("customFile", "", width = "100%")
                                )
                            )
                        ),
                        
                        uiOutput("cusSuperRankSelect.ui"),

                        h5(""),
                        bsButton(
                            "plotCustom",
                            "Do plot",
                            style = "warning",
                            icon("sync")
                        )
                    ),

                    # * main panel for customized profile plot -----------------
                    mainPanel(
                        conditionalPanel(
                            condition = "output.sameProfile == true",
                            h4(
                                "Please select subset of genes and/
                                or taxa for customized profile!"
                            )
                        ),
                        createProfilePlotUI("customizedProfile")
                    )
                )
            ),

            # FUNCTION TAB =====================================================
            navbarMenu(
                "Function",
                # * Distribution analysis --------------------------------------
                tabPanel(
                    "Distribution analysis",
                    h4(strong("Distribution analysis")),
                    bsAlert("descDistributionUI"),

                    wellPanel(
                        fluidRow(
                            column(
                                2,
                                selectInput(
                                    "dataset.distribution", "Select data",
                                    choices = c("Main data", "Customized data"),
                                    selected = "Main data"
                                ),
                                uiOutput("selected.distribution")
                            ),
                            column(
                                2, uiOutput("var1Dist.ui")
                            ),
                            column(
                                2, uiOutput("var2Dist.ui")
                            ),
                            column(
                                2, uiOutput("percentDist.ui")
                            ),
                            column(
                                2,
                                createTextSize(
                                    "distTextSize", "Label size", 12, 100
                                )
                            ),
                            column(
                                2,
                                createPlotSize(
                                    "distWidth", "Width (px)", 600
                                )
                            )
                        )
                    ),
                    analyzeDistributionUI("distPlot")
                )
            ),

            # DATA DOWNLOAD TAB ================================================
            navbarMenu(
                "Export data",
                # * Export data ------------------------------------------------
                downloadFilteredMainUI("filteredMainDownload"),
                downloadFilteredCustomizedUI("filteredCustomizedDownload"),
                
           ),

            # HELP TAB =========================================================
            navbarMenu(
                "Help",
                tabPanel(
                    a(
                        "Wiki",
                        href = "https://github.com/BIONF/PhyloProfile/wiki",
                        target = "_blank"
                    )
                ),
                tabPanel(
                    a(
                        "About",
                        href = "https://BIONF.github.io/PhyloProfile/",
                        target = "_blank"
                    )
                )
            )
        ),

        # LIST OF POP-UP WINDOWS ===============================================

        # * popup for plotting detailed plot -----------------------------------
        bsModal(
            "modalBs",
            "Detailed plot",
            "detailedBtn",
            size = "large",
            fluidRow(
                column(
                    2, createPlotSize("detailedHeight", "Height (px)", 100)
                ),
                column(
                    3, createTextSize("detailedText", "Text size (px)", 12, 150)
                ),
                column(
                    7,
                    checkboxInput(
                        "detailedRemoveNA",
                        strong("Hide taxa that have no ortholog (NAs)",
                               style = "color:red"),
                        value = FALSE
                    ),
                    checkboxInput(
                        "detailedFilter",
                        strong("Apply filters",
                               style = "color:red"),
                        value = FALSE
                    )
                )
            ),
            hr(),
            createDetailedPlotUI("detailedPlot"),
            bsButton(
                "doDomainPlot", "Show domain architecture", disabled = TRUE
            ),
            uiOutput("checkDomainFiles"),
            br(),
            h4("Sequence:"),
            verbatimTextOutput("fasta"),
            br(),
            h4("Links:"),
            uiOutput("dbLink.ui")
        ),

        # * popup for plotting domain architecture plot ------------------------
        bsModal(
            "plotArchi",
            "Domain architecture",
            "doDomainPlot",
            size = "large",
            fluidRow(
                column(
                    2, createPlotSize("archiHeight", "Plot height(px)", 400)
                ),
                column(
                    2, createPlotSize("archiWidth", "Plot width(px)", 800)
                ),
                column(
                    2,
                    createTextSize("titleArchiSize", "Title size(px)", 14, 150)
                ),
                column(
                    2,
                    createTextSize("labelArchiSize","X-axis size(px)",12,150)
                )
            ),
            createArchitecturePlotUI("archiPlot")
        ),
        bsModal(
            "plotArchiFromMain",
            "Domain architecture",
            "doDomainPlotMain",
            size = "large",
            fluidRow(
                column(
                    2, createPlotSize("archiHeightM", "Plot height(px)", 400)
                ),
                column(
                    2, createPlotSize("archiWidthM", "Plot width(px)", 800)
                ),
                column(
                    2,
                    createTextSize("titleArchiSizeM", "Title size(px)", 14, 150)
                ),
                column(
                    2,
                    createTextSize("labelArchiSizeM","X-axis size(px)",12,150)
                )
            ),
            createArchitecturePlotUI("archiPlotMain")
        ),

        # * popup for setting Main plot configurations -------------------------
        bsModal(
            "mainPlotConfigBs",
            "Plot appearance configuration",
            "mainPlotConfig",
            size = "small",
            column(
                6, createTextSize("xSize", "X-axis label size (px)", 14, 100)
            ),
            column(
                6, createTextSize("ySize", "Y-axis label size (px)", 14, 100)
            ),
            column(
                6,
                createTextSize("legendSize", "Legend label size (px)", 8, 150)
            ),
            column(
                6,
                selectInput(
                    "mainLegend", label = "Legend position:",
                    choices = list("Right" = "right",
                                   "Left" = "left",
                                   "Top" = "top",
                                   "Bottom" = "bottom",
                                   "Hide" = "none"),
                    selected = "right",
                    width = 150
                )
            ),
            column(
                6, 
                createTextSize("groupLabelSize", "Group label size (px)",7,100)
            ),
            column(
                6,
                numericInput(
                    "groupLabelDist", "Height for group label",
                    min = 0, max = 100, step = 1, value = 7, width = 100
                )
            ),
            column(
                12,
                HTML("<strong>Angle for taxonomic group label</strong>:<br>"),
                sliderInput(
                    "groupLabelAngle",
                    "",
                    min = 0,
                    max = 90,
                    step = 10,
                    value = 90,
                    width = 250
                ),
                br()
            ),
            column(
                12,
                HTML("<strong>Angle for x-axis label</strong>:<br>"),
                sliderInput(
                    "xAngle",
                    "",
                    min = 0,
                    max = 90,
                    step = 10,
                    value = 60,
                    width = 250
                ),
                br()
            ),
            column(
                12,
                HTML("<strong>Zooming factor (α) for dots on
                    profile</strong>:<br>"),
                sliderInput(
                    "dotZoom", "",
                    min = -1,
                    max = 3,
                    step = 0.1,
                    value = 0,
                    width = 250
                ),
                HTML("<em>dot size = (1+α)*defaultSize<br>defaultSize
                    =[0:5]</em>"),
                uiOutput("dotSizeInfo"),
                br()
            ),
            br(),
            strong(h4("Color configuration:")),
            actionButton(
                "setColor", "Change colors",
                style = "padding:4px; font-size:100%"
            ),
            hr(),
            bsButton("resetMainConfig", "Reset", style = "danger"),
            bsButton("applyMainConfig", "Done", style = "warning")
        ),
        
        # * popup for setting plot colors (profiles) ---------------------------
        bsModal(
            "color",
            "Set colors for profile",
            "setColor",
            size = "small",
            colourpicker::colourInput(
                "lowColorVar1",
                "Low variable 1 (dot)",
                value = "#FF8C00"
            ),
            colourpicker::colourInput(
                "midColorVar1",
                "Mid variable 1 (dot)",
                value = "#40ABCF"
            ),
            colourpicker::colourInput(
                "highColorVar1",
                "High variable 1 (dot)",
                value = "#164294"
            ),
            numericInput(
                "midVar1",
                "Mitpoint varriable 1",
                min = 0,
                max = 1,
                step = 0.01,
                value = 0.5
            ),
            actionButton(
                "defaultColorVar1",
                "Default",
                style = "padding:4px; font-size:100%"
            ),
            hr(),
            colourpicker::colourInput(
                "lowColorVar2",
                "Low variable 2 (background)",
                value = "#CC8D8D"
            ),
            colourpicker::colourInput(
                "midColorVar2",
                "Mid variable 2 (background)",
                value = "#FFFFFF"
            ),
            colourpicker::colourInput(
                "highColorVar2",
                "High variable 2 (background)",
                value = "#616587"
            ),
            numericInput(
                "midVar2",
                "Mitpoint varriable 2",
                min = 0,
                max = 1,
                step = 0.01,
                value = 1
            ),
            actionButton(
                "defaultColorVar2",
                "Default",
                style = "padding:4px; font-size:100%"
            ),
            hr(),
            colourpicker::colourInput(
                "paraColor",
                "Color for inparalogs",
                value = "#07d000"
            ),
            actionButton(
                "defaultColorPara",
                "Default",
                style = "padding:4px; font-size:100%"
            )
        ),

        # * popup for setting Customized plot configurations -------------------
        bsModal(
            "selectedPlotConfigBs",
            "Plot appearance configuration",
            "selectedPlotConfig",
            size = "small",
            column(
                6,
                createTextSize("xSizeSelect", "X-axis label size (px)", 14, 100)
            ),
            column(
                6,
                createTextSize("ySizeSelect", "Y-axis label size (px)", 14, 100)
            ),

            column(
                6,
                createTextSize("legendSizeSelect", "Legend label size (px)",
                               8, 150)
            ),
            column(
                6,
                selectInput(
                    "selectedLegend", label = "Legend position:",
                    choices = list("Right" = "right",
                                   "Left" = "left",
                                   "Top" = "top",
                                   "Bottom" = "bottom",
                                   "Hide" = "none"),
                    selected = "right",
                    width = 150
                )
            ),
            column(
                6, 
                createTextSize(
                    "groupLabelSizeSelect", "Group label size (px)", 7, 100
                )
            ),
            column(
                6,
                numericInput(
                    "groupLabelDistSelect", "Height for group label",
                    min = 0, max = 100, step = 1, value = 3, width = 100
                )
            ),
            column(
                12,
                HTML("<strong>Angle for taxonomic group label</strong>:<br>"),
                sliderInput(
                    "groupLabelAngleSelect",
                    "",
                    min = 0,
                    max = 90,
                    step = 10,
                    value = 90,
                    width = 250
                ),
                br()
            ),
            column(
                12,
                HTML("<strong>Angle for x-axis label</strong>:<br>"),
                sliderInput(
                    "xAngleSelect", "",
                    min = 0,
                    max = 90,
                    step = 10,
                    value = 60,
                    width = 250
                ),
                br()
            ),
            column(
                12,
                HTML("<strong>Zooming factor (α) for dots on
                     profile</strong>:<br>"),
                sliderInput(
                    "dotZoomSelect", "",
                    min = -1,
                    max = 3,
                    step = 0.1,
                    value = 0,
                    width = 250
                ),
                HTML("<em>dot size = (1+α)*defaultSize<br>
                     defaultSize=[0:5]</em>"),
                uiOutput("dotSizeInfoSelect"),
                br()
            ),
            br(),
            hr(),
            bsButton("resetSelectedConfig", "Reset", style = "danger"),
            bsButton("applySelectedConfig", "Done", style = "warning")
        ),

        # * popup for select taxa on Main Profile ------------------------
        bsModal(
            "highlight",
            "Select taxon/taxa of interest",
            "taxonHighlightBrowse",
            size = "small",
            selectTaxonRankUI("selectTaxonRankMain"),
            checkboxInput(
                "applyMainTaxa",
                strong("Apply to main profile",
                       style = "color:red"),
                value = FALSE
            )
        ),

        # * popup for select taxa on Core gene finding -------------------------
        bsModal(
            "browseTaxaCoreBs",
            "Select taxon/taxa of interest",
            "browseTaxaCore",
            size = "small",
            selectTaxonRankUI("selectTaxonRankCore"),
            checkboxInput(
                "applyCoreTaxa",
                strong("Apply", style = "color:red"),
                value = FALSE
            )
        ),

        # POINT INFO BOX =======================================================
        conditionalPanel(
            condition =
                "input.tabs=='Main profile' ||
                input.tabs=='Customized profile'",
            absolutePanel(
                bottom = 5, left = 30,
                fixed = TRUE,
                draggable = TRUE,
                h5("Point's info:"),
                verbatimTextOutput("pointInfo"),
                conditionalPanel(
                    condition = "output.pointInfoStatus == 0",
                    bsButton(
                        "detailedBtn",
                        "Detailed plot",
                        style = "success",
                        disabled = FALSE
                    ),
                    bsButton(
                        "doDomainPlotMain",
                        "Domain plot",
                        style = "success",
                        disabled = TRUE
                    )
                ),
                style = "opacity: 0.80"
            )
        )
    )
)
