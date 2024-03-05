#' Import function files
sourceFiles = list.files(path = "R", pattern = "*.R$", full.names = TRUE)
lapply(sourceFiles, source, .GlobalEnv)

#' MAIN UI ====================================================================
shinyUI(
    fluidPage(

        tags$style(type = "text/css", "body {padding-top: 80px;}"),
        useShinyjs(),

        # Application title
        titlePanel("", windowTitle = "PhyloCellulase"),

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
                            "autoSizing",
                            strong(em("Auto sizing")),
                            value = FALSE,
                            width = NULL
                        )
                    ),
                    column(
                        1,
                        createPlotSize("width", "Width (px)", 8800),
                        actionButton(
                            "setColor", "Change colors",
                            style = "padding:4px; font-size:100%"
                        )
                    ),
                    column(
                        1, createPlotSize("height", "Height (px)", 5000),
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
                        checkboxInput(
                            "autoUpdateSelected",
                            strong(em("Auto update plot")),
                            value = TRUE,
                            width = NULL
                        ),
                        checkboxInput(
                            "selectedAutoSizing",
                            strong(em("Auto sizing")),
                            value = TRUE,
                            width = NULL
                        )
                    ),
                    column(
                        1,
                        createPlotSize("selectedWidth", "Width (px)", 600)
                    ),
                    column(
                        1, createPlotSize("selectedHeight", "Height (px)", 600),
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
            em(strong("PhyloCellulase")),
            id = "tabs",
            collapsible = TRUE,
            inverse = TRUE,
            fluid = TRUE,
            position = "fixed-top",

            # MAIN PROFILE TAB =================================================
            tabPanel(
                "Main profile",
                bsAlert("warningUI"),
                sidebarLayout(
                    # * sidebar panel for profile highlight --------------------
                    sidebarPanel(
                        strong(("Select working taxonomy rank:")),
                        selectInput(
                            "rankSelect", label = "",
                            choices = c(
                                #"strain",
                                "species",
                                "genus",
                                "family",
                                "order",
                                "class",
                                "phylum",
                                "kingdom",
                                "superkingdom"
                            ), #PhyloProfile::getTaxonomyRanks(),
                            selected = "family"
                        ),
                        uiOutput("totalGeneNumber.ui"),
                        uiOutput("highlightGeneUI"),
                        bsPopover(
                            "highlightGeneUI",
                            "",
                            "Select gene to highlight",
                            "bottom"
                        ),
                        br(),
                        uiOutput("highlightTaxonUI"),

                        # column(
                        #     4,
                        #     numericInput(
                        #         "stIndex",
                        #         "Show from:",
                        #         min = 1,
                        #         max = 1600,
                        #         value = 1,
                        #         width = 100
                        #     ),
                        #     style = "padding:0px;"
                        # ),
                        #
                        # column(
                        #     4,
                        #     numericInput(
                        #         "endIndex",
                        #         "...to:",
                        #         min = 1,
                        #         max = 1600,
                        #         value = 1000,
                        #         width = 100
                        #     ),
                        #     style = "padding:0px;"
                        # ),

                        conditionalPanel(
                            condition = "input.autoUpdate == false",
                            bsButton(
                                "updateBtn",
                                "Update plot",
                                style = "warning",
                                icon("refresh")
                            )
                            # ,
                            # bsPopover(
                            #     "updateBtn",
                            #     "",
                            #     paste0(
                            #         "Due to the large amount of data, it will ",
                            #         "take a while for plotting!"
                            #     ),
                            #     "bottom"
                            # )
                        ),
                        width = 2
                    ),
                    # * main panel for profile plot ----------------------------
                    mainPanel(
                        createProfilePlotUI("mainProfile"),
                        width = 10
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
                            strong("Select sequence(s) of interest:")
                        ),

                        column(
                            12,
                            fluidRow(
                                column(
                                    8,
                                    style = "padding:0px;",
                                    uiOutput("geneIn")
                                ),
                                column(
                                    4,
                                    fileInput("customFile", "", width = "100%")
                                )
                            )
                        ),

                        column(
                            12,
                            style = "padding:0px;",
                            strong(
                                "Select (super)taxon/(super)taxa of interest:"
                            )
                        ),
                        column(
                            12,
                            fluidRow(
                                column(
                                    8,
                                    style = "padding:0px;",
                                    uiOutput("taxaIn")
                                )
                            )
                        ),

                        h5(""),
                        conditionalPanel(
                            condition = "input.autoUpdateSelected == false",
                            bsButton(
                                "plotCustom",
                                "Update plot",
                                style = "warning",
                                icon("refresh")
                            )
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
                        createProfilePlotUI("customizedProfile"),
                        width = 8
                    )
                )
            ),

            # DATA DOWNLOAD TAB ================================================
            navbarMenu(
                "Download filtered data",
                downloadFilteredMainUI("filteredMainDownload"),
                downloadFilteredCustomizedUI("filteredCustomizedDownload")
            ),

            # HELP TAB =========================================================
            navbarMenu(
                "Help",
                tabPanel(
                    a(
                        "About",
                        href = "https://BIONF.github.io/PhyloCellulase/",
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
                "doDomainPlot", "Show domain architecture", disabled = FALSE
            ),
            uiOutput("checkDomainFiles"),
            br(),
            h4("Sequence:"),
            verbatimTextOutput("fasta"),
            br(),
            h4("Links:"),
            uiOutput("dbLink")
        ),

        # * popup for plotting domain architecture plot ------------------------
        bsModal(
            "plotArchi",
            "Domain architecture",
            "doDomainPlot",
            size = "large",
            fluidRow(
                column(
                    2, createPlotSize("archiHeight", "Plot height(px)", 200)
                ),
                column(
                    2, createPlotSize("archiWidth", "Plot width(px)", 800)
                ),
                column(
                    2,
                    createTextSize("titleArchiSize", "Title size(px)", 11, 150)
                ),
                column(
                    2,
                    createTextSize("labelArchiSize", "SeqID size(px)", 11, 150)
                )
            ),
            uiOutput("test.ui"),
            createArchitecturePlotUI("archiPlot")
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
                value = "#40ABCF" #"#164294"
            ),
            colourpicker::colourInput(
                "highColorVar1",
                "High variable 1 (dot)",
                value = "#164294"# "steelblue"
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

        # * popup for setting Main plot configurations -------------------------
        bsModal(
            "mainPlotConfigBs",
            "Plot appearance configuration",
            "mainPlotConfig",
            size = "small",
            strong(em(
                paste(
                    "Note: due to the size limit, the font size and dot",
                    "zooming factor may not affect the plot on SPECIES level!"
                )
            )),
            hr(),
            column(
                6, createTextSize("xSize", "X-axis label size (px)", 12, 100)
            ),
            column(
                6, createTextSize("ySize", "Y-axis label size (px)", 12, 100)
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
                    value = -0.3,
                    width = 250
                ),
                HTML("<em>dot size = (1+α)*defaultSize<br>defaultSize
                    =[0:5]</em>"),
                uiOutput("dotSizeInfo"),
                br()
            ),
            br(),
            hr(),
            bsButton("resetMainConfig", "Reset", style = "danger"),
            bsButton("applyMainConfig", "Done", style = "warning")
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
                    )
                ),
                style = "opacity: 0.80"
            )
        )
    )
)
