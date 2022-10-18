
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

if (!require("utils", quietly = TRUE))
  install.packages("utils")

if (!require("INLA", quietly = TRUE)){
  # R.version <- readline(prompt = "Are you running this app with the last version of R [y/n]? ")
  writeLines("The first time, you may need to restart the application when the packages are finished installing. You would notice if the process still on but the app window is closed.")
  R.version <- menu(choices = c("Yes", "No"), title = paste0("Your R version is ", version$major, ".", version$minor, ". Are you running this app on the last version of R?"))
  if(R.version==1){
    install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
  } else { 
    install.packages("remotes")
    remotes::install_version("INLA", version="22.05.03",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
    }
}
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
library(lattice)
library(rintrojs)
library(patchwork)
library(viridis)
library(rgeos)
library(dplyr)

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
        tabItem(tabName = "introduction",
                fluidRow(
                  column(
                    # style= 'padding-bottom:15px',
                    width = 12,
                    mainPanel(
                      h1("Introduction", style="color: DarkRed; text-align: center; text-decoration: underline; font-weight: bold"),
                      withMathJax(),
                      p("The funtionality of this shiny app is to solve two kinds of spatial models using the Bayesian approach and the INLA methodology. 
                        Bayesian approach is one statistic paradigm in which the parameters are considered as random variables, so they will have distribution functions estimated though Bayes' theorem. 
                        The INLA methodology is a deterministic approximation for hierarchical bayesian models with an underlayng Gaussian field structure for the parameters.
                        The two kinds of models to resolve are (i) a geostatistical model and (ii) a specific joint model; which are now briefly described. 
                        Although they will be explained in more detail in their corresponding sections, where you will find a", span("Summary", style="font-style: italic"), "to show such information."),
                      p("The structure of the geoestatistical model for the response variable \\(y\\) is the following, in which the likelihood function \\(f(\\cdot)\\) of observations \\(y\\) is related to the predictor by means of the link function \\(g(\\cdot)\\):"),
                      helpText("\\begin{equation}
                                \\begin{array}{c}
                               y_i \\sim f(y_i|\\boldsymbol\\theta),\\\\
                               g(E(y_i)) = g(\\mu_i) = \\beta_0 + \\mathbf{X}_i\\boldsymbol\\beta + u_i.
                               \\end{array}
                               \\end{equation}"),
                      p(paste0("Where \\(\\beta_0\\) is the intercept, ", "\\(\\mathbf{X}\\) ", "are the explanatory variables and ", 
                               "\\(\\boldsymbol\\beta\\) the parameters, ", 
                               "meanwhile ", "\\(u_i\\) ", "is the spatial effect.")),
                      p("The specific joint model we are using combines two likelihoods, one related to the geostatistical data and another related to the locations. 
                        Between the two likelihoods we have one or more shared components (for us only the spatial term is shared), 
                        what implies that their estimates would be done though both likelihoods. 
                        The structure of this joint model is:"),
                      helpText("\\begin{equation}
                                \\begin{array}{c}
                               y_i \\sim f(y_i|\\boldsymbol\\theta),\\\\
                               p_i \\sim Poisson(p_i|\\boldsymbol\\theta),\\\\
                               g(E(y_i)) = g(\\mu_i)= \\beta_0 + \\mathbf{X}_i\\boldsymbol\\beta + u_i,\\\\
                               \\log(E(\\lambda_i)) = \\log(\\lambda_i)= \\beta_0' + \\mathbf{X}_i \\boldsymbol\\beta' + \\alpha \\cdot u_i.
                               \\end{array}
                               \\end{equation}"),
                      p("Where the intercept for the geostatistical layer is \\(\\beta_0\\) and 
                        \\(\\beta_0'\\) is the one for the point process layer,
                        \\(\\mathbf{X}\\) are the explanatory variables, \\(\\boldsymbol\\beta\\) the parameters of the geostatistical layer and
                        \\(\\boldsymbol\\beta'\\) the parameters of the point process layer. Meanwhile \\(u_i\\) is the spatial effect and
                        \\(\\alpha\\) is a internal shared parameter, a scale term."),
                      br(),
                      p("Therefore, the application allows to simulate an example data to test the inference results controlling its parameters and also allows data loading and modeling, 
                        according to the model structures showed previosly. The simulation is done in the", span("Simulation Data", style = "font-style: italic"),  
                        "section, in which is possible to control the parameters that regulates the simulation of the covariables, intercept, spatial effect and the geostatistical data output. 
                        In te next section", span("(Model Analysis)", style="font-style: italic"), "is possible to shape the inference through some buttons implemented in the UI, which make it possible to control the basic prior parameters for intercept, explanatory variables, spatial effect and the likelihood function, 
                        as well as for specific INLA controll elements. These controll items are defined in", span("Advance INLA config.", style="font-style: italic"), "where we can define the modes for the hyperparameters, the approximation type and the integration strategy."),
                      h1("Glosary", style="color: DarkRed; text-align: center; text-decoration: underline; font-weight: bold"),
                      br(),
                      p("Since this application is intended for all kind of users, a glossary is included containing a brief selection of terms, with some explanations, as they are used along the whole application.
                        This glosary has two ordered levels, one by field and the second ordered by the terms belonging to that specific field."),
                      tags$div(
                        HTML("
                        <ul> <strong>Bayesian inference:</strong>
                        <ol>
                        <li> <i>Bayes' theorem (conditional probability distribution):</i> </li>
                        <li> <i>Hierarchical model:</i> </li>
                        <li> <i>Likelihood:</i> </li>
                        <li> <i>Link function:</i> </li>
                        <li> <i>Marginal distribution(s):</i> </li>
                        <li> <i>Mean:</i> </li>
                        <li> <i>Mode:</i> </li>
                        <li> <i>Parameter(s):</i> </li>
                        <li> <i>Posterior distribution:</i> </li>
                        <li> <i>Posterior predictive distribution:</i> </li>
                        <li> <i>Prior distribution:</i> </li>
                        <li> <i>Prior predictive distribution:</i> </li>
                        </ol>
                        </ul>
                        
                        <ul> <strong>INLA basic elements:</strong>
                        <ol>
                        <li> <i>Aproximation strategy:</i> </li>
                        <li> <i>Integration strategy:</i> </li>
                        <li> <i>Latent field:</i> </li>
                        <li> <i>Gaussian Field (GMRF):</i> </li>
                        <li> <i>Gaussian Markov Random Field (GMRF):</i> </li>
                        </ol>
                        </ul>
                        
                        <ul> <strong>INLA spatial approach:</strong>
                        <ol>
                        <li> <i>Finite Element Method (FEM):</i> </li>
                        <li> <i>Matérn covariance function:</i> </li>
                        <li> <i>Mesh:</i> </li>
                        <li> <i>Stochastic Partial Differential Equation (SPDE)</i>: </li>
                        </ol>
                        </ul>
                             
                             "
                             ))
                    )
                  )
                )

        ),
        
        # Map and sample simulation tab content
        tabItem(tabName = "datasimulation",
                conditionalPanel(
                  condition = "input.infoSim",
                fluidRow(
                  column(
                    # style= 'padding-bottom:15px',
                    width = 12,
                    mainPanel(
                      h1("Summary", style="color: DarkRed; text-align: center; text-decoration: underline; font-weight: bold"),
                      withMathJax(),
                      p("In this section we will explain in some detail the structure and code used for data simulation. The main data \\(y\\) is reated to 
                        the geostatistical likelihood, which is built from the model structure data: intercept, bathymetry (one numerical covariate) and spatial data. "),
                      p("As we have showed previosly the model structure is:"),
                      helpText("\\begin{equation}
                                \\begin{array}{c}
                               y_i \\sim f(y_i|\\boldsymbol\\theta),\\\\
                               g(\\mu_i)= \\beta_0 +\\mathbf{X}_i\\boldsymbol\\beta + u_i.
                               \\end{array}
                               \\end{equation}"),
                      p("But it must be clear that since the app is built for spatial modeling the data is set in a subregion \\(\\mathcal{D}\\subset\\mathbb{R}^2\\). 
                        Therefore the first elements to configure are related to the study region specifications:"),
                      tags$div(
                        HTML("<ol>
                                <li> Analysis area: a square region governed by the user-specified upper and lower limits.</li>
                                <li> Vector dimension: the grid dimension for spatial simulation. </li>
                                <li> Map resolution: the final grid dimensión resulting from a <i>Finite Elemente Method</i> (FEM) projection of the first grid simulation, which means a map resolution increase. </li>
                             </ol>")
                      ),
                      tags$div(
                        HTML("Subsequently , two items are available to control the randomness of the stochastic process: <i>Global seed</i> and <i> Spatial effect seed</i>. 
                             The first rules the randomness of the whole geostatistic process simulation, meanwhile the latter controls the specific randomness of the spatial effect simulation.")),
                      br(),
                      tags$div(
                        HTML(
                          "With respect to the simulation of the spatial effect, given that it is characterized by the \\(u_i \\sim N(\\mathbf{0}, \\boldsymbol\\Sigma) \\) expression, 
                          where the variace matrix is determined by a Matérn covariance function \\(\\Sigma=\\sigma \\cdot C_{\\nu}(d, \\rho)\\) with two parameters: spatial range \\(\\rho\\) and standard marginal deviation \\(\\sigma\\). Two numerical inputs are availabled in the UI to stablish the values for such parameters.")
                        ),
                      br(),
                      tags$div(
                        HTML(
                          "For the covariate we can define the formula (<i>Bathymetric formula</i>) and choose the effect type:
                          <ol> 
                          <li> Linear: the value of its effect will be set in <i> Predictor coefficients</i>.</li>
                          <li> Random walk (1o): we can configure the first order random walk \\(\\Delta x_i= x_i - x_{i-1}=N(0, \\tau)\\) through the following parameters: an initial value \\(x_0\\), the precision \\(\\tau\\), the number of knots for simulate the random walk, and the limit values which can't be traspassing by any \\(x_i\\) of the random walk.</li>
                          <li> Random walk (2o): to structure the second order random walk \\(\\Delta^2 x_i= x_i - 2\\cdot x_{i-1} + x_{i-2} = N(0, \\tau)\\) we must define an initial value \\(x_0\\), the precision \\(\\tau\\), the number of knots for simulate the random walk, and the limit values which can't be traspassing by any \\(x_i\\) of the random walk.</li>
                          <li> Custom function: it's a continous univariate function that configure the relation beetween the covariate values and its effect, \\(y=f(x)\\).</li>
                          </ol>
                          Since the linear effect and the custom function effect are a simple covariate transformation, they do not need any further explanation. But it could be done for the random walk to get some insight. 
                          The random walk is a stochastic process, such its implementation over the simulate covariate was design in three steps: 
                          (i) first, the <n> knots are located at equally distances between the minimum and maximum values of the covariate, 
                          (ii) a linear interpolation between every knot is dont by taking every simulated value of the covariate and 
                          (iii) finally, it is checked if any resulting value exceeds the indicated limits for the effect, so that if it does, the two previous steps are repeated until the condition is met. The <i>predictor coefficients</i> allow to define the itercept and the linear effect, if this was chosen for the covariate.</i>")
                        ),
                      br(),
                      p("Our response variable is \\(y\\), a variable for the geostatistical process implies a continuous variable for the likelihood function or data distribution function. 
                        Then a Gamma and Gaussian function could be used by means of its corresponding UI button, allowing to set the distribution for data output:"),
                      helpText("$$f(\\cdot)=\\{N(\\cdot) \\veebar Gamma(\\cdot)\\}.$$"),
                      p("We must also indicate the value of the variance associated with the likelihood distribution."),
                      br(),
                      tags$div(
                        HTML("Finally, we explain the section dedicated to the simulation of both <i>independent</i> (iid) and <i>preferential</i> sampling.
                        For each of which an option is available to specify the seed of the stochastic sampling process and the number of samples desired.
                             Independent sampling is simply performed using the <code>sample</code> function on the raster resulting from the simulation of the geostatistical data. 
                             Meanwhile, since the specific joint model assumes than the geostatistical and the point process are linked, ")
                      ),
                      helpText("\\begin{equation}
                                \\begin{array}{c}
                               y_i \\sim f(y_i|\\boldsymbol\\theta),\\\\
                               p_i \\sim Poisson(p_i|\\boldsymbol\\theta),\\\\
                               g(\\mu_i)= \\mathbf{X}_i\\boldsymbol\\beta + u_i,\\\\
                               \\log(\\lambda_i)= \\mathbf{X}_i \\boldsymbol\\beta' + \\alpha \\cdot u_i;
                               \\end{array}
                               \\end{equation}"),
                      tags$div(
                        HTML("such that the sampling is preferred by the geoestatistical data. Therefore, for the preferential sampling we will have that the probability associated to each 'point' of the raster is:")
                      ),
                      helpText("\\begin{equation}
                               p_i=\\frac{\\exp[r\\cdot y_i]}{\\sum_i\\exp[r \\cdot y_i]\\cdot S_i}  ) \\propto \\exp[ r\\cdot y_i].
                               \\end{equation}"),
                      tags$div(
                        HTML("Where \\( r \\) is the preferential degree, how much the geostatistical values influed in the sampling, \\(y_i\\) is the value of the observed variable and \\( S_i \\) is surface related to the raster pixel \\( i \\). 
                             Then, \\(\\mathbf{p}=\\{p_1,p_2,...,p_n\\}\\) is the vector we will use in the argument <i>prob</i> of <code>sample(x, size, prob)</code> for the preferential sampling.")
                      ),
                      br()
                    )
                  ))
                ),
                fluidRow(
                    column(
                        style = 'padding-bottom:15px',
                        width = 12,
                        introBox(
                            bsButton("parametersMapSim", 
                                     label = "Geostatistical parameters",
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
                                    title = "Geostatistical parameters",
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
                                                 label="Spatial effect seed:",
                                                 value=1),
                                    numericInput("sigma.matern", #width = '45%',
                                                 label=HTML(paste0("Matérn stdev. (\U03C3", tags$sub("0"),"):")),
                                                 value=0.5), 
                                    numericInput("range.matern", #width = '45%',
                                                 label=HTML(paste0("Matérn range (\U03C1", tags$sub("0"),"):")),
                                                 value=0.2)),
                                    textInput("bathymetry.formula",
                                              label="Bathymetry formula (x, y)",
                                              value="0.5 + (x-min(x)) + 0.8*(y-mean(y))**2 + rnorm(length(x), 0, 0.01)",
                                              placeholder="sqrt(abs(x+y))+x**2+log(y-min(y)+1)"),
                                    selectInput("typeBathyeff", label="Types of bathymetry effects",
                                                choices=c("Linear"="lin", "Random Walk (1o)"="rw1",  
                                                          "Random Walk (2o)"="rw2", 
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
                                    numericInput("var",
                                                 label=HTML(paste0("Distribution variance (\U03C3", tags$sup("2"),"):")),
                                                 value=0.2),
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
                                    footer = "The (iid) are the parameters for the independent sampling, 
                                    and the (ps) are the parameters for the preferential sampling.",
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
                        title = "Bathymetric map and its effect", solidHeader = TRUE,
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
                        title = "Bathymetric effect map and Abundance map", solidHeader = TRUE,
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
                conditionalPanel(
                  condition = "input.infoSimLoad",
                  fluidRow(
                    column(
                      # style= 'padding-bottom:15px',
                      width = 12,
                      mainPanel(
                        h1("Summary", style="color: DarkRed; text-align: center; text-decoration: underline; font-weight: bold"),
                        withMathJax(),
                        tags$div(
                          HTML("In this section we briefly explain how data reading works in this app. For this purpose, there are two elements in the UI: 
                               <ol>
                               <li><strong>Load Analysis Data Frame:</strong> it allows to read the loaded data frame, but the two first columns must be the \\((x,y)\\) coordinates, the third the observational variable \\(y\\) and the others the explanatory variables \\(\\mathbf(X)\\).
                               <li><strong>Load Covariates Data Frame:</strong> it allows to read a data frame realated to some covariate and using it for the prediction grid. The two first columns must be the \\((x,y)\\) and the third the covariate.</li>
                              </ol>
                               Once any data frame is loaded its table and a plot will appear in the interface.  
                               ")
                        ),
                        br()
                      )
                    ))
                ),    
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
                              label="Load Analysis Data Frame:",
                              placeholder = "No file selected (.csv or .rds)"),
                    fileInput("file.uploadDataBathymetryRaster",
                              label="Load Covariates Data Frame:",
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
                conditionalPanel(
                  condition = "input.infoRandomFit",
                  fluidRow(
                    column(
                      # style= 'padding-bottom:15px',
                      width = 12,
                      mainPanel(
                        h1("Summary", style="color: DarkRed; text-align: center; text-decoration: underline; font-weight: bold"),
                        withMathJax(),
                        tags$div(
                          HTML("")
                        ),
                        br()
                      )
                    ))
                ),
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
                selectInput("DataSimulatedLoaded", label="Select data to be analyzed",
                            choices=list("Simulated"="sim", "Loaded"="load"),
                            selected="sim"),
                selectInput("BathymetryRasterSPDE", label="Bathymetry for data prediction",
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
                box(id="RandomMesh", width=12, title="Custom Mesh", closable = FALSE,
                    status="warning", collapsible = TRUE, collapsed = TRUE, solidHeader = TRUE,
                    materialSwitch(inputId="RandomCustomMesh", label="Allow custom mesh", 
                                   value=FALSE, status="primary"),
                    numericInput(inputId="RandomMeshQloc", label="Qloc",
                                 width="100%", value=0.03, min=0, max=1, step=0.01)
                    
                ),
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
                selectInput("optionRPB", label="Intercept",#"Prior options (another linear covariates, e.g. Intercept)",
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
                conditionalPanel(condition="input.RandomCustomMesh",
                                 box(title="Independent Model Mesh", width=9, solidHeader=TRUE, status="info",
                                     downloadFileButton(id="ggplotRandomMesh",
                                                        downloadtypes = c("png"),
                                                        hovertext = "Download image and data"),
                                     downloadablePlotUI(id = "ggplotRandomMesh",
                                                        list("png"),
                                                        btn_valign = "bottom",
                                                        btn_halign = "right",
                                                        btn_overlap = TRUE))
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
                conditionalPanel(
                  condition = "input.infoPrefFit",
                  fluidRow(
                    column(
                      # style= 'padding-bottom:15px',
                      width = 12,
                      mainPanel(
                        h1("Summary", style="color: DarkRed; text-align: center; text-decoration: underline; font-weight: bold"),
                        withMathJax(),
                        tags$div(
                          HTML("")
                        ),
                        br()
                      )
                    ))
                ),
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
                                       selectInput("PrefDataSimulatedLoaded", label="Select data to be analyzed",
                                                   choices=list("Simulated"="sim", "Loaded"="load"),
                                                   selected="sim"),
                                       selectInput("PrefBathymetryRasterSPDE", label="Bathymetry for data prediction",
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
                                           selectInput("optionPPB", label="Intercept",#"Prior options (another linear covariates, e.g. Intercept)",
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
                                           # bsButton("PrefBetaCopy", label = "Prior for Beta-copy",
                                           #          type = "toggle", 
                                           #          value = FALSE),
                                           # conditionalPanel(condition="input.PrefBetaCopy",
                                           #                  textInput("PrefpriorBetacopy", label="Values",
                                           #                            value="0,0.001", placeholder = "0,0.001")),
                                           textInput("PrefpriorBetacopy", label="Prior Beta-copy values",
                                                     value="0,0.001", placeholder = "0,0.001"),
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
