if (!require("shinydashboard", quietly = TRUE))
  install.packages("shinydashboard")
if (!require("shinyWidgets", quietly = TRUE))
  install.packages("shinyWidgets")
if (!require("shinydashboardPlus", quietly = TRUE))
  install.packages("shinydashboardPlus")
if (!require("shinybusy", quietly = TRUE))
  install.packages("shinybusy")
if (!require("shinyBS", quietly = TRUE))
  install.packages("shinyBS")
if (!require("shinyjs", quietly = TRUE))
  install.packages("shinyjs")

if (!require("periscope", quietly = TRUE))
  install.packages("periscope")
if (!require("splines", quietly = TRUE))
  install.packages("splines")
if (!require("INLA", quietly = TRUE))
  install.packages("INLA")
if (!require("inlabru", quietly = TRUE))
  install.packages("inlabru")
if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
if (!require("lattice", quietly = TRUE))
  install.packages("lattice")
if (!require("rintrojs", quietly = TRUE))
  install.packages("rintrojs")
if (!require("patchwork", quietly = TRUE))
  install.packages("patchwork")
if (!require("viridis", quietly = TRUE))
  install.packages("viridis")
if (!require("rgeos", quietly = TRUE))
  install.packages("rgeos")

library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinydashboardPlus)
library(shinyBS)
library(shinyjs)
library(periscope)
library(splines)
library(INLA)
library(inlabru)
library(ggplot2)
library(dplyr)
library(lattice)
library(rintrojs)
library(patchwork)
library(viridis)
library(rgeos)

# library(DT)
# library(data.table)
# library(canvasXpress)
# library(htmlwidgets)
# library(colourpicker)
# library(plotly)
# library(tools)
# library(processx)
# library(fields)
# library(dipsaus)

#add_reset_button("./App_Abundance_Modelling/")

# Some functions to the UI

# js_test <- "
#     $(document).ready(function() {
#       $('#Equation').qtip({
#         overwrite: true,
#         content: {
#           text: $('#tooltip')
#         },
#         position: {
#           my: 'top left',
#           at: 'bottom right'
#         },
#         show: {
#           ready: false
#         },
#         hide: {
#           event: 'unfocus'
#         },
#         style: {
#           classes: 'qtip-youtube qtip-rounded'
#         },
#         events: {
#           blur: function(event, api) {
#             api.elements.tooltip.hide();
#           }
#         }
#       });
#     });
# "

# section below allows in-line LaTeX via $ in mathjax.


# modify_stop_propagation <- function(x) {
#     x$children[[1]]$attribs$onclick = "event.stopPropagation()"
#     x
# }
# 
# convertMenuItem <- function(mi,tabName) {
#     mi$children[[1]]$attribs['data-toggle']="tab"
#     mi$children[[1]]$attribs['data-value'] = tabName
#     if(length(mi$attribs$class)>0 && mi$attribs$class=="treeview"){
#         mi$attribs$class=NULL
#     }
#     mi
# }

# Define UI for application

header <- dashboardHeader(
    title = div(img(src = "logoUV-cropped.svg", height = 45), span(HTML("Spatial R-INLA"), style={'padding-left: 15px'})),
    titleWidth = 400

    #Dropdown Menu
    # dropdownMenu(type = "messages",
    #              messageItem(
    #                  from = "Sales Dept",
    #                  message = "Sales are steady this month."
    #              ),
    #              messageItem(
    #                  from = "New User",
    #                  message = "How do I register?",
    #                  icon = icon("question"),
    #                  time = "13:45"
    #              ),
    #              messageItem(
    #                  from = "Support",
    #                  message = "The new server is ready.",
    #                  icon = icon("life-ring"),
    #                  time = "2014-12-01"
    #              )),
    # 
    # ,dropdownMenu(type = "notifications",
    #              notificationItem(
    #                  text = "5 new users today",
    #                  icon("users")
    #              ),
    #              notificationItem(
    #                  text = "12 items delivered",
    #                  icon("truck"),
    #                  status = "success"
    #              ),
    #              notificationItem(
    #                  text = "Server load at 86%",
    #                  icon = icon("exclamation-triangle"),
    #                  status = "warning"
    #              )),
    # 
    # dropdownMenu(type = "tasks", badgeStatus = "success")
)
## Sidebar content
sidebar <- dashboardSidebar(minified = FALSE, width=200,
    sidebarMenu(id = "sidebarID",
        menuItem("Introduction", tabName = "introduction", icon = icon("exclamation-circle")),
        menuItem("Data Simulation", tabName = "datasimulation", icon = icon("database")),
        menuItem("Upload Data", tabName = "dataloading", icon = icon("copy", lib = "glyphicon")),
        menuItem("Model Analysis", tabName = "modelanalysis", icon = icon("chart-line"), expandedName = "MAexpand",
                 menuSubItem("Independent Model", tabName="modiid"),
                 menuSubItem("Preferential Model", tabName="modpref"))
        # menuItem("Feedback Model", tabName = "feedbackmodelanalysis", icon = icon("chart-line"), expandedName = "FMexpand",
        #          menuSubItem("Independent Model", tabName="feedmodiid"),
        #          menuSubItem("Preferential Model", tabName="feedmodpref"))
    )
)

## Body content
body <- dashboardBody(
    tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
    ),
    introjsUI(),
    shinybusy::add_busy_spinner(spin = "atom", margins=c(20,20), position="bottom-right", ),# fading-circle
    # useShinyjs(),
    tabItems(
        # Introducton tab content
        tabItem(tabName = "introduction"

        ),
        
        # Map and sample simulation tab content
        tabItem(tabName = "datasimulation",
                fluidRow(
                    column(
                        style = 'padding-bottom:15px',
                        width = 12,
                        introBox(
                            bsButton("parametersMapSim", 
                                     label = "Map parameters",
                                     type = "toggle", 
                                     value = TRUE,
                                     icon = icon("drafting-compass"), 
                                     style = "primary"),
                            bsButton("parametersSampleSim", 
                                     label = "Sample parameters", 
                                     type = "toggle", 
                                     value = TRUE,
                                     icon = icon("spinner", class = "spinner-box"), 
                                     style = "primary"),
                            bsButton("infoSim", 
                                     label = "Summary information", 
                                     type = "toggle", 
                                     value = FALSE,
                                     icon = icon("info"), 
                                     style = "warning")
                        )
                    )
                ),
                fluidRow(
                    conditionalPanel(
                        condition = "input.parametersMapSim|input.parametersSampleSim",
                        box(width = 3,
                            conditionalPanel(
                                condition = "input.parametersMapSim",
                                box(width = 12, status = "info",
                                    title = "Map parameters",
                                    solidHeader = TRUE,
                                    collapsible = TRUE,
                                    actionButton("makeSim", 
                                             label = "Simulate", 
                                             #type = "action",
                                             #block = TRUE,
                                             width="100%",
                                             #style = "unite",
                                             icon = icon("jedi-order")),
                                    sliderInput("limlattice",
                                                label = "Analysis area:",
                                                min = -5,
                                                max = 5,
                                                value = c(0,1)),
                                    sliderInput("lengthlattice",
                                                label = "Vector dimension (grid):",
                                                min=10,
                                                max=200,
                                                value = 100),
                                    sliderInput("dimmap",
                                                label = "Map resolution:",
                                                min=50,
                                                max=300,
                                                value = 150),
                                    numericInput("seedGlobal",
                                                 label="Global seed:",
                                                 value=1),
                                    box(width = 12, title="Spatial parameters",
                                        status="info", solidHeader=TRUE, collapsible=TRUE,
                                    numericInput("seedSP",
                                                 label="Seed (Spatial effect):",
                                                 value=1),
                                    numericInput("sigma.matern", #width = '45%',
                                                 label=HTML(paste0("Matérn stdev. (\U03C3", tags$sub("0"),"):")),
                                                 value=0.5), 
                                    numericInput("range.matern", #width = '45%',
                                                 label=HTML(paste0("Matérn range (\U03C1", tags$sub("0"),"):")),
                                                 value=0.2)),
                                    numericInput("var",
                                                 label=HTML(paste0("Distribution variance (\U03C3", tags$sup("2"),"):")),
                                                 value=0.2),
                                    textInput("bathymetry.formula",
                                              label="Bathymetry formula (x, y)",
                                              value="0.5 + (x-min(x)) + 0.8*(y-mean(y))**2 + rnorm(length(x), 0, 0.01)",
                                              placeholder="sqrt(abs(x+y))+x**2+log(y-min(y)+1)"),
                                    selectInput("typeBathyeff", label="Types of bathymetry effects",
                                                choices=c("Lineal"="lin", "Random Walk (1d)"="rw1",  
                                                          "Random Walk (2d)"="rw2", 
                                                          "Custom function"="customfunction")
                                    ),
                                    conditionalPanel(condition = "input.typeBathyeff=='rw1'", 
                                                     box(status="info", title="Parameters RW1", width=12, solidHeader = TRUE,
                                                         numericInput("prec.rw1", label="Precision (RW1)", value=10, min=0),
                                                         numericInput("init.rw1", label="Init value (RW1)", value=0.2),
                                                         numericInput("nknots.rw1", label="Number of knots (RW1)", value=60),
                                                         textInput("minmax.effect.rw1",
                                                                   label="Extreme values",
                                                                   value="0, 5",
                                                                   placeholder="sqrt(abs(x+y))+x**2+log(y-min(y)+1)"))
                                    ),
                                    conditionalPanel(condition = "input.typeBathyeff=='rw2'", 
                                                     box(status="info", title="Parameters RW2", width=12, solidHeader = TRUE,
                                                         numericInput("prec.rw2", label="Precision (RW2)", value=10, min=0),
                                                         numericInput("init.rw2", label="Init value (RW2)", value=0.2),
                                                         numericInput("nknots.rw2", label="Number of knots (RW2)", value=60),
                                                         textInput("minmax.effect.rw2",
                                                                   label="Extreme values",
                                                                   value="0, 5",
                                                                   placeholder="sqrt(abs(x+y))+x**2+log(y-min(y)+1)"))
                                    ),
                                    conditionalPanel(condition = "input.typeBathyeff=='customfunction'",
                                                     box(status="info", title="Custom function", width=12, solidHeader = TRUE,
                                                         textInput("formula.eff.bathymetry",
                                                                   label="Bathymetry effect formula (b)",
                                                                   value="0.3 + 0.5*b + 1.5*(b-mean(b))**3 + rnorm(length(b), 0, 0.01)",
                                                                   placeholder="0.1+b+b**2"))
                                    ),
                                    textInput("beta",
                                              label=HTML(paste0("Predictor coefficients (\U03B2", tags$sub("0"),
                                                                ", \U03B2",tags$sub("1"),", ...", "):")),
                                              value="0,1",
                                              placeholder = "0,1"),
                                    selectInput("datadistributionSim", label="Data distribution",
                                                choices=c("Gamma distribution" = "gamma", "Gaussian distribution" = "gaussian")
                                    )
                                )),
                            conditionalPanel(
                                condition = "input.parametersSampleSim",
                                box(width = 12, status = "info",
                                    title = "Sample parameters",
                                    solidHeader = TRUE,
                                    collapsible = TRUE,
                                    footer = "The (iid) are the parameters to the independent sampling process, 
                                    and the (ps) are the parameters to the preferential sampling process.",
                                    actionButton("makeSample",
                                                 label = "Sample",
                                                 width="100%",
                                                 icon = icon("jedi-order")),
                                    numericInput("seedSampleR",
                                                 label="Seed (iid):",
                                                 value=1),
                                    numericInput("niid.samples", "Number of samples (iid)",
                                                 value = 60),
                                    numericInput("seedSampleP",
                                                 label="Seed (ps):",
                                                 value=1),
                                    numericInput("nps.samples", "Number of samples (ps)",
                                                 value = 60),
                                    numericInput("r.scale", "Scale factor (ps)",
                                                 value = 8)
                                ))
                            )
                    ),
                    #plot2_hover <- hoverOpts(id = "examplePlot2_hover"),
                    conditionalPanel(condition="input.makeSim>=1", 
                    box(id = "resultbox", width=9, collapsible=FALSE,
                    box(id = "mesh.speff.box", width = 12, status = "info",
                        title = "Mesh and Spatial effect", solidHeader = TRUE,
                        collapsible = TRUE, collapsed = FALSE,
                        column(width=6,
                               downloadFileButton(id="ggplotMeshSim",
                                                  downloadtypes = c("png", "csv", "txt"),
                                                  hovertext = "Download image and data"),
                               downloadablePlotUI(id = "ggplotMeshSim",
                                                  list("png", "txt", "tsv"),
                                                  btn_valign = "bottom",
                                                  btn_halign = "right",
                                                  btn_overlap = TRUE)),
                        column(width=6,
                               downloadFileButton(id = "ggplotSpSim",
                                                  downloadtypes = c("png", "csv", "txt"),
                                                  hovertext = "Download image and data"),
                               downloadablePlotUI(id = "ggplotSpSim",
                                                  list("png", "txt", "tsv"),
                                                  btn_valign = "bottom",
                                                  btn_halign = "right",
                                                  btn_overlap = TRUE))
                        ),
                    box(id = "bathymetry.bath.effect", width = 12, status = "info",
                        title = "Bathymetric chart and its effect", solidHeader = TRUE,
                        collapsible = TRUE, collapsed = FALSE,
                        column(width=6,
                               downloadFileButton(id="ggplotBatChartSim",
                                                  downloadtypes = c("png", "csv", "txt"),
                                                  hovertext = "Download image and data"),
                               downloadablePlotUI(id = "ggplotBatChartSim",
                                                  list("png", "txt", "tsv"),
                                                  btn_valign = "bottom",
                                                  btn_halign = "right",
                                                  btn_overlap = TRUE)),
                        column(width=6,
                               downloadFileButton(id="ggplotBatEffSim",
                                                  downloadtypes = c("png", "csv", "txt"),
                                                  hovertext = "Download image and data"),
                               downloadablePlotUI(id = "ggplotBatEffSim",
                                                  list("png", "txt", "tsv"),
                                                  btn_valign = "bottom",
                                                  btn_halign = "right",
                                                  btn_overlap = TRUE))
                        ),
                    box(id = "abundance.batheffect.maps", width = 12, status = "info",
                        title = "Bathymetric effect chart and Abundance map", solidHeader = TRUE,
                        collapsible = TRUE, collapsed = FALSE,
                        column(width=6,
                               downloadFileButton(id="ggplotBatEffChartSim",
                                                  downloadtypes = c("png", "csv", "txt"),
                                                  hovertext = "Download image and data"),
                               downloadablePlotUI(id = "ggplotBatEffChartSim",
                                                  list("png", "txt", "tsv"),
                                                  btn_valign = "bottom",
                                                  btn_halign = "right",
                                                  btn_overlap = TRUE)),
                        column(width=6,
                               downloadFileButton(id="ggplotAbundanceSim",
                                                  downloadtypes = c("png", "csv", "txt"),
                                                  hovertext = "Download image and data"),
                               downloadablePlotUI(id = "ggplotAbundanceSim",
                                                  list("png", "txt", "tsv"),
                                                  btn_valign = "bottom",
                                                  btn_halign = "right",
                                                  btn_overlap = TRUE))
                    ),
                    conditionalPanel(condition="input.makeSample>=1",
                    box(id = "sample.maps", width = 12, status = "info",
                        title = "Independent and preferential samples", solidHeader = TRUE,
                        collapsible = TRUE, collapsed = FALSE,
                        column(width=6,
                               downloadFileButton(id="ggplotRandomSampleSim",
                                                  downloadtypes = c("png", "csv", "txt"),
                                                  hovertext = "Download image and data"),
                               downloadablePlotUI(id="ggplotRandomSampleSim",
                                                  list("png", "txt", "tsv"),
                                                  btn_valign = "bottom",
                                                  btn_halign = "right",
                                                  btn_overlap = TRUE)),
                        column(width=6,
                               downloadFileButton(id="ggplotPrefSampleSim",
                                                  downloadtypes = c("png", "csv", "txt"),
                                                  hovertext = "Download image and data"),
                               downloadablePlotUI(id = "ggplotPrefSampleSim",
                                                  list("png", "txt", "tsv"),
                                                  btn_valign = "bottom",
                                                  btn_halign = "right",
                                                  btn_overlap = TRUE))
                    ))
                    ))
                )
        ),
        
        # Own data loading tab content
        tabItem(tabName = "dataloading",
                fluidRow(
                    column(
                        style = 'padding-bottom:15px',
                        width = 12,
                        introBox(
                            bsButton("LoadingPanel", 
                                     label = "Loading Panel",
                                     type = "toggle", 
                                     value = TRUE,
                                     icon = icon("file-export"), 
                                     style = "primary"),
                            bsButton("infoSimLoad", 
                                     label = "Summary information", 
                                     type = "toggle", 
                                     value = FALSE,
                                     icon = icon("info"), 
                                     style = "warning")
                        )
                    )
                ),
                fluidRow(
                    conditionalPanel(condition="input.LoadingPanel",
                    box(width=4,
                    fileInput("file.uploadData",
                              label="Load a Data Frame (observations):",
                              placeholder = "No file selected (.csv or .rds)"),
                    fileInput("file.uploadDataBathymetryRaster",
                              label="Load a Bathymetric raster Data Frame:",
                              placeholder = "No file selected (.csv or .rds)")
                )),
                conditionalPanel(condition="output.fileUploadedSample|output.fileUploadedRaster",
                box(width=8, id="TablePanelRaster",
                    # tableOutput("content")
                    conditionalPanel("output.fileUploadedSample",
                    column(width=6,
                    downloadableTableUI("table.read.sample",
                                        downloadtypes=c("csv", "tsv")
                                        ))),
                    conditionalPanel("output.fileUploadedRaster",
                    column(width=6,
                           downloadableTableUI("table.read.bathymetry.raster",
                                               downloadtypes=c("csv", "tsv")
                           ))),
                    conditionalPanel("output.fileUploadedSample",
                    column(width=6,
                           downloadFileButton(id="quilt.plot.sample.bathymetry",
                                              downloadtypes = c("png"),
                                              hovertext = "Download image"),
                           downloadablePlotUI(id = "quilt.plot.sample.bathymetry",
                                              list("png"),
                                              btn_valign = "bottom",
                                              btn_halign = "right",
                                              btn_overlap = TRUE)),
                    column(width=6,
                           downloadFileButton(id="quilt.plot.sample.abundance",
                                              downloadtypes = c("png"),
                                              hovertext = "Download image"),
                           downloadablePlotUI(id = "quilt.plot.sample.abundance",
                                              list("png"),
                                              btn_valign = "bottom",
                                              btn_halign = "right",
                                              btn_overlap = TRUE))),
                    conditionalPanel("output.fileUploadedRaster",
                    column(width=12,
                    downloadFileButton(id="quilt.plot.bathymetry",
                                       downloadtypes = c("png"),
                                       hovertext = "Download image"),
                    downloadablePlotUI(id = "quilt.plot.bathymetry",
                                       list("png"),
                                       btn_valign = "bottom",
                                       btn_halign = "right",
                                       btn_overlap = TRUE)))
                )))
        ),

        # Independent and Random Model
        tabItem(tabName="modiid",withMathJax(),
                fluidRow(
                column(
                    style = 'padding-bottom:15px',
                    width = 12,
                introBox(
                bsButton("FitRandomPanel", 
                         label = "Fit Panel",
                         type = "toggle", 
                         value = TRUE,
                         icon = icon("drafting-compass"), 
                         style = "primary"),
                bsButton("infoRandomFit", 
                         label = "Summary information", 
                         type = "toggle", 
                         value = FALSE,
                         icon = icon("info"), 
                         style = "warning")
                ))),
                fluidRow(
                conditionalPanel(condition="input.FitRandomPanel",
                box(width=3, 
                actionButton("fitIndRandom", label="Fit Model", width="100%",
                             icon=icon("jedi-order")),
                selectInput("DataSimulatedLoaded", label="Select data to analyse",
                            choices=list("Simulated"="sim", "Loaded"="load"),
                            selected="sim"),
                selectInput("BathymetryRasterSPDE", label="Bathymetry to data prediction",
                            choices=list("Raster"="raster", "Solve as SPDE"="solvebathy"),
                            selected="solvebathy"),
                conditionalPanel(condition="input.BathymetryRasterSPDE=='raster'",
                selectInput("rasterBatymetryPred", label="Options to prediction",
                            selected="SPDEraster",
                            choices=list("Prediction grid through raster (SPDE)"="SPDEraster",
                                         # "Use raster to make a prediction grid through the A (INLA) matrix 
                                         # (if and only if the raster is a regular grid)"="smoothraster",
                                         "Raster as prediction locations"="rasterpred")
                            ),
                conditionalPanel(condition="input.rasterBatymetryPred=='SPDEraster'",
                                 sliderInput("SPDErasterInput", "Predictive grid dimensionality:",
                                             min=50, max=200, value=100)),
                conditionalPanel(condition="input.rasterBatymetryPred=='smoothraster'",
                                 sliderInput("smoothrasterInput", "Predictive grid dimensionality:",
                                             min=50, max=200, value=100))
                ),
                conditionalPanel(condition="input.BathymetryRasterSPDE=='solvebathy'",
                sliderInput("dimrandompred",
                            label = "Predictive grid dimensionality:",
                            min = 50,
                            max = 200,
                            value = 100)),
                sliderInput("dimrandommap",
                            label= "Predictive map resolution (spatial effect):",
                            min = 100,
                            max = 300,
                            value = 150),
                box(id="RandomPriorBeta", width=12, title="Prior parameters (fixed effects)",
                    status="info", solidHeader=TRUE, collapsible=TRUE,
                box(title="Bathymetry", width=12, status="info", solidHeader=TRUE, collapsible=TRUE,
                selectInput("RPModelBathy", label="Bathymetric model",
                            choices=list("Lineal"="lin", "B-spline"="bs", 
                                         "Random Walk (1d)"="rw1", "Random Walk (2d)"="rw2",
                                         "SPDE (1d)"="spde"),
                            selected="lin"),
                conditionalPanel(condition="input.RPModelBathy=='lin'",
                                 radioGroupButtons(inputId = "autocustomLinBathy", 
                                   choices = list("Auto"="auto", "Custom"="custom"),
                                   status = "primary"),
                                 conditionalPanel(condition="input.autocustomLinBathy=='custom'",
                                 numericInput("RandomMeanPriorLinBathymetry",
                                           label=HTML("Mean prior") , value=0),
                                 numericInput("RandomPrecPrioLinrBathymetry",
                                           label=HTML("Precision prior"), value=0.001))),
                conditionalPanel(condition="input.RPModelBathy=='bs'",
                     numericInput("bsRandomnKnots", label="Number of knots",
                                  value=6, min=1, step=1)),
                conditionalPanel(condition="input.RPModelBathy=='rw1'",
                                 numericInput("rw1RandomnKnots", label="Number of knots",
                                              value=6, min=1, step=1),
                                 radioGroupButtons(inputId = "autocustomRandomRw1", 
                                                   label="Prior RW1",
                                                   choices = list("Default"="default", "Custom"="custom"),
                                                   status = "primary"),
                                 conditionalPanel(condition="input.autocustomRandomRw1=='custom'",
                                                  textInput("RandomPrecRw1", label="Prior loggamma (mean, prec) to prec (rw1)",
                                                            value="1,0.001", placeholder="1,0.001"))),
                conditionalPanel(condition="input.RPModelBathy=='rw2'",
                                 numericInput("rw2RandomnKnots", label="Number of knots",
                                              value=6, min=1, step=1),
                                 radioGroupButtons(inputId = "autocustomRandomRw2", 
                                                   label="Prior RW2",
                                                   choices = list("Default"="default", "Custom"="custom"),
                                                   status = "primary"),
                                 conditionalPanel(condition="input.autocustomRandomRw2=='custom'",
                                 textInput("RandomPrecRw2", label="Prior precision (mean, prec) to RW (rw2)",
                                           value="1,0.001", placeholder="1,0.001"))),
                conditionalPanel(condition="input.RPModelBathy=='spde'",
                                 numericInput("spdeRandomnKnots", label="Number of knots",
                                              value=6, min=1, step=1),
                                 radioGroupButtons(inputId = "autocustomRandomSpde", 
                                                   label="Prior SPDE",
                                                   choices = list("Auto"="auto", "Custom"="custom"),
                                                   status = "primary"),
                                 conditionalPanel(condition="input.autocustomRandomSpde=='custom'",
                                                  textInput("RandomTheta1Spde", 
                                                            label=HTML(paste0("Prior ", "\U03F4", tags$sub("1"), " (mean, prec)")),
                                                            value="0,0.01", placeholder="0,0.01"),
                                                  textInput("RandomTheta2Spde", 
                                                            label=HTML(paste0("Prior ", "\U03F4", tags$sub("2"), " (mean, prec)")),
                                                            value="0,0.01", placeholder="0,0.01")))
                ),
                selectInput("optionRPB", label="Interceptor",#"Prior options (another linear covariates, e.g. Interceptor)",
                            choices=list("Default"="default", "Custom"="custom"),
                            selected="default"),
                conditionalPanel(condition="input.optionRPB=='custom'",
                textInput("RandomMeanPriorBeta",
                          label=HTML("Mean prior vector (\U03B2)") , value="0", placeholder = "0"),
                textInput("RandomPrecPriorB",
                          label=HTML("Precision prior vector (\U03B2)"), value="0.001", placeholder = "0.001"))
                ),
                box(id="RandomPriorSpatial", width=12, title="Penalized priors (spatial effects)",
                    status="info", solidHeader=TRUE, collapsible=TRUE,
                    selectInput("optionRPS", label="Prior options",
                                choices=list("Auto"="auto", "Custom"="custom"),
                                selected="auto"),
                    conditionalPanel(condition="input.optionRPS=='custom'",
                             textInput("RandomPriorRange",
                                       label=HTML(paste0(p(HTML(paste0("Range Matérn ", "(\U03C1", tags$sub("0"), ", p", tags$sub("1"), "):" ))), p("\\( \\mathbf{P}(\\boldsymbol\\rho<\\boldsymbol\\rho_o)=\\mathbf{p_1} \\)"))), 
                                       value="0.5,0.5", placeholder = "diff(limlattice)/2,0.5"),
                             textInput("RandomPriorStdev",
                                       label=HTML(paste0(p(HTML(paste0("Stdev. Matérn ", "(\U03C3", tags$sub("0"), ", p", tags$sub("2"), "):" ))), p("\\( \\mathbf{P}(\\boldsymbol\\sigma<\\boldsymbol\\sigma_o)=\\mathbf{p_2} \\)"))), 
                                       value="0.5,0.5", placeholder = "diff(limlattice)/2,0.5")
                             )),
                box(id="RandomFamily", width=12, title="Family and its hyperparameters",
                    status="info", solidHeader=TRUE, collapsible=TRUE,
                    selectInput("SelectRandomFamily", label="Family distribution",
                                choices=list("Gamma"="gamma","Gaussian"="gaussian"),
                                selected="gamma"),
                    radioGroupButtons(inputId = "autocustomRandomFamily", 
                                      label="Hyperparameters",
                                      choices = list("Default"="default", "Custom"="custom"),
                                      status = "primary"),
                    conditionalPanel(condition="input.autocustomRandomFamily=='custom'",
                                     textInput("RandomFamilyHyper", label="Log-Gamma parameters",
                                               value="0.01, 0.01", placeholder = "0.01, 0.01"))
                    ),
                box(id="advancedINLARandom", width=12,
                    title="Advanced INLA config.",
                    collapsible=TRUE, collapsed=TRUE,
                    status = "info", solidHeader=TRUE,
                    selectInput("strategyapproxINLARandom",
                                 label = "Aproximation strategy",
                                 choices=list("Auto"="auto", "Gaussian"="gaussian",
                                              "Simplified Laplace"="simplified.laplace",
                                              "Laplace"="laplace", "Adaptive"="adaptive"),
                                  selected = "auto"),
                    selectInput("strategyintINLARandom",
                                label = "Integration strategy",
                                choices=list("Auto"="auto", "Central compositive design"="ccd",
                                             "Grid strategy"="grid",
                                             "Empirical Bayes"="eb"),
                                selected = "auto"),
                    radioGroupButtons(inputId = "autocustomRandomMode", 
                                      label="Hyperparameters modes",
                                      choices = list("Default"="default", "Custom"="custom"),
                                      status = "primary"),
                    conditionalPanel(condition="input.autocustomRandomMode=='custom'",
                                     textInput("RandomModeHyper", label="Values",
                                               value="0,0", placeholder = "0,0"))
                    
                    ))
                ),
                conditionalPanel(condition="input.fitIndRandom>=1",
                box(id="plotstablesFitRandom", width=9,
                box(title="Abundance predictive maps", solidHeader=TRUE,
                    collapsible=TRUE, width=12, status="info",
                column(width=12,
                    downloadFileButton(id="ggplotAbundanceMeanMedianStedevFit",
                                       downloadtypes = c("png", "csv", "txt"),
                                       hovertext = "Download image and data"),
                    downloadablePlotUI(id = "ggplotAbundanceMeanMedianStedevFit",
                                       downloadtypes = c("png", "csv", "txt"),
                                       btn_valign = "bottom",
                                       btn_halign = "right",
                                       btn_overlap = TRUE))),
                box(title="Posterior maps of the spatial effect", solidHeader=TRUE,
                    collapsible=TRUE, width=12, status="info",
                column(width=12,
                       downloadFileButton(id="ggplotSpatialMeanMedianStdev",
                                          downloadtypes = c("png", "csv", "txt"),
                                          hovertext = "Download image and data"),
                       downloadablePlotUI(id = "ggplotSpatialMeanMedianStdev",
                                          downloadtypes = c("png", "csv", "txt"),
                                          btn_valign = "bottom",
                                          btn_halign = "right",
                                          btn_overlap = TRUE))),
                box(title="Density functions of the fixed effects", solidHeader=TRUE,
                    collapsible=TRUE, width=12, status="info",
                column(width=12,
                       downloadFileButton(id="ggplotFixParamFit",
                                          downloadtypes = c("png", "csv", "txt"),
                                          hovertext = "Download image and data"),
                       downloadablePlotUI(id = "ggplotFixParamFit",
                                          downloadtypes = c("png", "csv", "txt"),
                                          btn_valign = "bottom",
                                          btn_halign = "right",
                                          btn_overlap = TRUE))),
                box(title="Density functions of the hyperparameters", solidHeader=TRUE,
                    collapsible=TRUE, width=12, status="info",
                column(width=12,
                       downloadFileButton(id="ggplotHyperParamFit",
                                          downloadtypes = c("png", "csv", "txt"),
                                          hovertext = "Download image and data"),
                       downloadablePlotUI(id = "ggplotHyperParamFit",
                                          downloadtypes = c("png", "csv", "txt"),
                                          btn_valign = "bottom",
                                          btn_halign = "right",
                                          btn_overlap = TRUE))),
                box(title="Tables of fixed effects, hyperparameters and internal hyperparameters", solidHeader=TRUE,
                    collapsible=TRUE, width=12, status="info",
                column(width=6,
                       downloadableTableUI("tableRandomFixedPar", downloadtypes=c("csv", "tsv"))
                       ),
                column(width=6,
                       downloadableTableUI("tableRandomHyperPar", downloadtypes=c("csv", "tsv"))
                       ),
                column(width=12,
                       downloadableTableUI("tableRandomInternalHyperPar", downloadtypes=c("csv", "tsv"))),
                column(width=8,
                       downloadableTableUI("dataCPOtable", downloadtypes=c("csv", "tsv"))),
                column(width=4,
                       downloadableTableUI("dataDICtable", downloadtypes=c("csv", "tsv"))
                ))
                ))
        )),
        tabItem(tabName="modpref", withMathJax(),
                fluidRow(
                  column(
                    style = 'padding-bottom:15px',
                    width = 12,
                    introBox(
                      bsButton("FitPrefPanel", 
                               label = "Fit Panel",
                               type = "toggle", 
                               value = TRUE,
                               icon = icon("drafting-compass"), 
                               style = "primary"),
                      bsButton("infoPrefFit", 
                               label = "Summary information", 
                               type = "toggle", 
                               value = FALSE,
                               icon = icon("info"), 
                               style = "warning")
                    ))),
                fluidRow(
                  conditionalPanel(condition="input.FitPrefPanel",
                                   box(width=3, 
                                       actionButton("fitPref", label="Fit Model", width="100%",
                                                    icon=icon("jedi-order")),
                                       selectInput("PrefDataSimulatedLoaded", label="Select data to analyse",
                                                   choices=list("Simulated"="sim", "Loaded"="load"),
                                                   selected="sim"),
                                       selectInput("PrefBathymetryRasterSPDE", label="Bathymetry to data prediction",
                                                   choices=list("Raster"="raster", "Solve as SPDE"="solvebathy"),
                                                   selected="solvebathy"),
                                       conditionalPanel(condition="input.PrefBathymetryRasterSPDE=='raster'",
                                            selectInput("PrefrasterBatymetryPred", label="Options to prediction",
                                                        selected="SPDEraster",
                                                        choices=list("Prediction grid through raster (SPDE)"="SPDEraster") #,
                                                                     # "Use raster to make a prediction grid through the A (INLA) matrix 
                                                                     # (if and only if the raster is a regular grid)"="smoothraster",
                                                                     # "Raster as prediction locations"="rasterpred")
                                            ),
                                              conditionalPanel(condition="input.PrefrasterBatymetryPred=='SPDEraster'",
                                                               sliderInput("SPDErasterInput", "Predictive grid dimensionality:",
                                                                           min=50, max=200, value=100)),
                                              conditionalPanel(condition="input.PrefrasterBatymetryPred=='smoothraster'",
                                                               sliderInput("smoothrasterInput", "Predictive grid dimensionality:",
                                                                           min=50, max=200, value=100))
                                       ),
                                       conditionalPanel(condition="input.PrefBathymetryRasterSPDE=='solvebathy'",
                                                        sliderInput("dimrandompred",
                                                                    label = "Predictive grid dimensionality:",
                                                                    min = 50,
                                                                    max = 200,
                                                                    value = 100)),
                                       sliderInput("dimprefmap",
                                                   label= "Predictive map resolution (spatial effect):",
                                                   min = 100,
                                                   max = 300,
                                                   value = 150),
                                       box(id="PrefPointProcessMesh", width=12, title="Custom Mesh", closable = FALSE,
                                           status="warning", collapsible = TRUE, collapsed = TRUE, solidHeader = TRUE,
                                           materialSwitch(inputId="PPPCustomMesh", label="Allow custom mesh", 
                                                          value=FALSE, status="primary"),
                                           numericInput(inputId="PrefPPMeshQloc", label="Qloc",
                                                        width="100%", value=0.03, min=0, max=1, step=0.01)
                                           
                                       ),
                                       box(id="PointPrefProcces", width=12, title="Point Process", closable = FALSE,
                                         status="warning", collapsible = TRUE, collapsed = TRUE, solidHeader = TRUE,
                                         actionButton(inputId="FitPointProcess", label="Fit point process",
                                                      width="100%", icon=icon("jedi-order"))
                                       ),
                                       box(id="PrefPriorBeta", width=12, title="Prior parameters (fixed effects)",
                                           status="info", solidHeader=TRUE, collapsible=TRUE,
                                           box(title="Bathymetry", width=12, status="info", solidHeader=TRUE, collapsible=TRUE,
                                               selectInput("PPModelBathy", label="Bathymetric model",
                                                           choices=list("Lineal"="lin", #"B-spline"="bs", 
                                                                        "Random Walk (1d)"="rw1", "Random Walk (2d)"="rw2",
                                                                        "SPDE (1d)"="spde"),
                                                           selected="lin"),
                                               conditionalPanel(condition="input.PPModelBathy=='lin'",
                                                                radioGroupButtons(inputId = "autocustomPrefLinBathy", 
                                                                                  choices = list("Auto"="auto", "Custom"="custom"),
                                                                                  status = "primary"),
                                                                conditionalPanel(condition="input.autocustomPrefLinBathy=='custom'",
                                                                                 numericInput("PrefMeanPriorLinBathymetry",
                                                                                              label=HTML("Mean prior") , value=0),
                                                                                 numericInput("PrefPrecPrioLinrBathymetry",
                                                                                              label=HTML("Precision prior"), value=0.001))),
                                               # conditionalPanel(condition="input.PPModelBathy=='bs'",
                                               #                  numericInput("bsPrefnKnots", label="Number of knots",
                                               #                               value=6, min=1, step=1)),
                                               conditionalPanel(condition="input.PPModelBathy=='rw1'",
                                                                numericInput("rw1PrefnKnots", label="Number of knots",
                                                                             value=6, min=1, step=1),
                                                                radioGroupButtons(inputId = "autocustomPrefRw1", 
                                                                                  label="Prior RW1",
                                                                                  choices = list("Default"="default", "Custom"="custom"),
                                                                                  status = "primary"),
                                                                conditionalPanel(condition="input.autocustomPrefRw1=='custom'",
                                                                                 textInput("PrefPrecRw1", label="Prior loggamma (mean, prec) to prec (rw1)",
                                                                                           value="1,0.001", placeholder="1,0.001"))),
                                               conditionalPanel(condition="input.PPModelBathy=='rw2'",
                                                                numericInput("rw2PrefnKnots", label="Number of knots",
                                                                             value=6, min=1, step=1),
                                                                radioGroupButtons(inputId = "autocustomPrefRw2", 
                                                                                  label="Prior RW2",
                                                                                  choices = list("Default"="default", "Custom"="custom"),
                                                                                  status = "primary"),
                                                                conditionalPanel(condition="input.autocustomPrefRw2=='custom'",
                                                                                 textInput("PrefPrecRw2", label="Prior precision (mean, prec) to RW (rw2)",
                                                                                           value="1,0.001", placeholder="1,0.001"))),
                                               conditionalPanel(condition="input.PPModelBathy=='spde'",
                                                                numericInput("spdePrefnKnots", label="Number of knots",
                                                                             value=6, min=1, step=1),
                                                                radioGroupButtons(inputId = "autocustomPrefSpde", 
                                                                                  label="Prior SPDE",
                                                                                  choices = list("Auto"="auto", "Custom"="custom"),
                                                                                  status = "primary"),
                                                                conditionalPanel(condition="input.autocustomPrefSpde=='custom'",
                                                                                 textInput("PrefTheta1Spde", 
                                                                                           label=HTML(paste0("Prior ", "\U03F4", tags$sub("1"), " (mean, prec)")),
                                                                                           value="0,0.01", placeholder="0,0.01"),
                                                                                 textInput("PrefTheta2Spde", 
                                                                                           label=HTML(paste0("Prior ", "\U03F4", tags$sub("2"), " (mean, prec)")),
                                                                                           value="0,0.01", placeholder="0,0.01")))
                                           ),
                                           selectInput("optionPPB", label="Interceptor",#"Prior options (another linear covariates, e.g. Interceptor)",
                                                       choices=list("Default"="default", "Custom"="custom"),
                                                       selected="default"),
                                           conditionalPanel(condition="input.optionPPB=='custom'",
                                                            textInput("PrefMeanPriorBeta",
                                                                      label=HTML("Mean prior vector (\U03B2)") , value="0", placeholder = "0"),
                                                            textInput("PrefPrecPriorB",
                                                                      label=HTML("Precision prior vector (\U03B2)"), value="0.001", placeholder = "0.001"))
                                       ),
                                       box(id="PrefPriorSpatial", width=12, title="Penalized priors (spatial effects)",
                                           status="info", solidHeader=TRUE, collapsible=TRUE,
                                           selectInput("optionPPS", label="Prior options",
                                                       choices=list("Auto"="auto", "Custom"="custom"),
                                                       selected="auto"),
                                           conditionalPanel(condition="input.optionPPS=='custom'",
                                                            textInput("PrefPriorRange",
                                                                      label=HTML(paste0(p(HTML(paste0("Range Matérn ", "(\U03C1", tags$sub("0"), ", p", tags$sub("1"), "):" ))), p("\\( \\mathbf{P}(\\boldsymbol\\rho<\\boldsymbol\\rho_o)=\\mathbf{p_1} \\)"))), 
                                                                      value="0.5,0.5", placeholder = "diff(limlattice)/2,0.5"),
                                                            textInput("PrefPriorStdev",
                                                                      label=HTML(paste0(p(HTML(paste0("Stdev. Matérn ", "(\U03C3", tags$sub("0"), ", p", tags$sub("2"), "):" ))), p("\\( \\mathbf{P}(\\boldsymbol\\sigma<\\boldsymbol\\sigma_o)=\\mathbf{p_2} \\)"))), 
                                                                      value="0.5,0.5", placeholder = "diff(limlattice)/2,0.5")
                                           )),
                                       box(id="PrefFamily", width=12, title="Family and its hyperparameters",
                                           status="info", solidHeader=TRUE, collapsible=TRUE,
                                           selectInput("SelectPrefFamily", label="Family distribution",
                                                       choices=list("Gamma"="gamma","Gaussian"="gaussian"),
                                                       selected="gamma"),
                                           radioGroupButtons(inputId = "autocustomPrefFamily", 
                                                             label="Hyperparameters",
                                                             choices = list("Default"="default", "Custom"="custom"),
                                                             status = "primary"),
                                           conditionalPanel(condition="input.autocustomPrefFamily=='custom'",
                                                            textInput("PrefFamilyHyper", label="Log-Gamma parameters",
                                                                      value="0.01, 0.01", placeholder = "0.01, 0.01"))
                                       ),
                                       box(id="advancedINLAPref", width=12,
                                           title="Advanced INLA config.",
                                           collapsible=TRUE, collapsed=TRUE,
                                           status = "info", solidHeader=TRUE,
                                           selectInput("strategyapproxINLAPref",
                                                       label = "Aproximation strategy",
                                                       choices=list("Auto"="auto", "Gaussian"="gaussian",
                                                                    "Simplified Laplace"="simplified.laplace",
                                                                    "Laplace"="laplace", "Adaptive"="adaptive"),
                                                       selected = "auto"),
                                           selectInput("strategyintINLAPref",
                                                       label = "Integration strategy",
                                                       choices=list("Auto"="auto", "Central compositive design"="ccd",
                                                                    "Grid strategy"="grid",
                                                                    "Empirical Bayes"="eb"),
                                                       selected = "auto"),
                                           bsButton("PrefBetaCopy", label = "Prior to Beta-copy",
                                                    type = "toggle", 
                                                    value = FALSE),
                                           conditionalPanel(condition="input.PrefBetaCopy",
                                                            textInput("PrefpriorBetacopy", label="Values",
                                                                      value="0,0.001", placeholder = "0,0.001")),
                                           radioGroupButtons(inputId = "autocustomPrefMode", 
                                                             label="Hyperparameters modes",
                                                             choices = list("Default"="default", "Custom"="custom"),
                                                             status = "primary"),
                                           conditionalPanel(condition="input.autocustomPrefMode=='custom'",
                                                            textInput("PrefModeHyper", label="Values",
                                                                      value="0,0", placeholder = "0,0"))
                                           
                                       ))
                  ),
                  conditionalPanel(condition="input.PPPCustomMesh",
                                   box(title="Mesh (Preferential and Point Process)", width=9, solidHeader=TRUE, status="info",
                                   downloadFileButton(id="ggplotPPPMesh",
                                                      downloadtypes = c("png"),
                                                      hovertext = "Download image and data"),
                                   downloadablePlotUI(id = "ggplotPPPMesh",
                                                      list("png"),
                                                      btn_valign = "bottom",
                                                      btn_halign = "right",
                                                      btn_overlap = TRUE))
                                   ),
                  conditionalPanel(condition="input.FitPointProcess>=1",
                  box(title="Point process results", width=9, collapsible=TRUE, solidHeader=TRUE, status="primary",
                                           column(width=6,
                                                  downloadableTableUI("tablePointProcessFixedPar", downloadtypes=c("csv", "tsv"))
                                           ),
                                           column(width=6,
                                                  downloadableTableUI("tablePointProcessHyperPar", downloadtypes=c("csv", "tsv"))
                                           ))
                      ),
                  conditionalPanel(condition="input.fitPref>=1",
                                   box(title="Preferential model results", width=9, collapsible=TRUE, solidHeader=TRUE, status="primary",
                                       box(title="Abundance predictive maps", solidHeader=TRUE,
                                           collapsible=TRUE, width=12, status="info",
                                           column(width=12,
                                                  downloadFileButton(id="ggplotPrefAbundanceMeanMedianStedevFit",
                                                                     downloadtypes = c("png", "csv", "txt"),
                                                                     hovertext = "Download image and data"),
                                                  downloadablePlotUI(id = "ggplotPrefAbundanceMeanMedianStedevFit",
                                                                     downloadtypes = c("png", "csv", "txt"),
                                                                     btn_valign = "bottom",
                                                                     btn_halign = "right",
                                                                     btn_overlap = TRUE))),
                                       box(title="Posterior maps of the spatial effect", solidHeader=TRUE,
                                           collapsible=TRUE, width=12, status="info",
                                           column(width=12,
                                                  downloadFileButton(id="ggplotPrefSpatialMeanMedianStdev",
                                                                     downloadtypes = c("png", "csv", "txt"),
                                                                     hovertext = "Download image and data"),
                                                  downloadablePlotUI(id = "ggplotPrefSpatialMeanMedianStdev",
                                                                     downloadtypes = c("png", "csv", "txt"),
                                                                     btn_valign = "bottom",
                                                                     btn_halign = "right",
                                                                     btn_overlap = TRUE))),
                                       box(title="Density functions of the fixed effects", solidHeader=TRUE,
                                           collapsible=TRUE, width=12, status="info",
                                           column(width=12,
                                                  downloadFileButton(id="ggplotPrefFixParamFit",
                                                                     downloadtypes = c("png", "csv", "txt"),
                                                                     hovertext = "Download image and data"),
                                                  downloadablePlotUI(id = "ggplotPrefFixParamFit",
                                                                     downloadtypes = c("png", "csv", "txt"),
                                                                     btn_valign = "bottom",
                                                                     btn_halign = "right",
                                                                     btn_overlap = TRUE))),
                                       box(title="Density functions of the hyperparameters", solidHeader=TRUE,
                                           collapsible=TRUE, width=12, status="info",
                                           column(width=12,
                                                  downloadFileButton(id="ggplotPrefHyperParamFit",
                                                                     downloadtypes = c("png", "csv", "txt"),
                                                                     hovertext = "Download image and data"),
                                                  downloadablePlotUI(id = "ggplotPrefHyperParamFit",
                                                                     downloadtypes = c("png", "csv", "txt"),
                                                                     btn_valign = "bottom",
                                                                     btn_halign = "right",
                                                                     btn_overlap = TRUE))),
                                       box(title="Tables of fixed effects, hyperparameters and internal hyperparameters", solidHeader=TRUE,
                                           collapsible=TRUE, width=12, status="info",
                                           column(width=6,
                                                  downloadableTableUI("tablePrefModelFixedPar", downloadtypes=c("csv", "tsv"))
                                           ),
                                           column(width=6,
                                                  downloadableTableUI("tablePrefModelHyperPar", downloadtypes=c("csv", "tsv"))
                                           ),
                                           column(width=12,
                                                  downloadableTableUI("tablePrefModelInternalHyperPar", downloadtypes=c("csv", "tsv"))),
                                           column(width=8,
                                                  downloadableTableUI("dataPrefCPOtable", downloadtypes=c("csv", "tsv"))),
                                           column(width=4,
                                                  downloadableTableUI("dataPrefDICtable", downloadtypes=c("csv", "tsv"))
                                           )) 
                                       ))
                )#,
        # tabItem(tabName="feedmodiid", "Independent Feedback Model"
        # ),
        # tabItem(tabName="feedmodpref", "Preferential Feedback Model "
        # )
    )
))

controlbar <-  dashboardControlbar(
    collapsed = TRUE,
    
    tags$head(
      tags$style(
        "
.body{
    background: #333;
    height: 100vh;
    display: flex;
    justify-content: center;
    align-items: center;
}

.button{
    border: none;
    background: linear-gradient(to right, #E32636, #BE0032);
    color: #fff;
    padding: .5em 1em;
    border-radius: 0.3em;
    font-size: 1em;
    cursor: pointer;
    position: relative;
    overflow: hidden;
}

.button::before{
    content: '';
    height: 160%;
    width: 40px;
    background: rgba(255,255,255,.3);
    position: absolute;
    top: -10px;
    transform: translateX(-80px) rotate(20deg);
    transition: all .5s;
}

.button:hover::before{
    transform: translateX(230px) rotate(20deg);
}
"
      )
    ),
    br(),
    column(width = 12,
           actionButton("EnhancementIrrGrid", label="Enhance irregular grid images", class = "button")),
    br(), br(), br(),
    column(width=12,
           materialSwitch(inputId = "RefineInterpolationIrrGrid", label = "Alternative enhancement", value = FALSE, status = "primary"))
           # switchInput(inputId="RefineInterpolationIrrGrid", label="Refine", size = "small")
           
           
    # bsButton("parameasdassadtersMapSim", 
    #          label = "Map parameters",
    #          type = "toggle", 
    #          value = TRUE,
    #          class="button")
)
    
    

    
    # , tags$button(
    #     id = 'close',
    #     type = "button",
    #     class = "btn action-button",
    #     onclick = "setTimeout(function(){window.close();},500);",  # close browser
    #     "Close window")
    # )

dashboardPage(skin = "red", header, sidebar, body, controlbar)
