if (!require("shinydashboard", quietly = TRUE)) {
  install.packages("shinydashboard")
}
if (!require("shinyWidgets", quietly = TRUE)) {
  install.packages("shinyWidgets")
}
if (!require("shinydashboardPlus", quietly = TRUE)) {
  install.packages("shinydashboardPlus")
}
if (!require("shinybusy", quietly = TRUE)) {
  install.packages("shinybusy")
}
if (!require("shinyBS", quietly = TRUE)) {
  install.packages("shinyBS")
}
if (!require("shinyjs", quietly = TRUE)) {
  install.packages("shinyjs")
}

if (!require("periscope", quietly = TRUE)) {
  install.packages("periscope")
}

if (!require("utils", quietly = TRUE)) {
  install.packages("utils")
}

if (!require("stringr", quietly = TRUE)) {
  install.packages("stringr")
}

if (!require("INLA", quietly = TRUE)) {
  
}
if (!require("inlabru", quietly = TRUE)) {
  install.packages("inlabru")
}
if (!require("fmesher", quietly = TRUE)) {
  remotes::install_github("inlabru-org/fmesher", ref = "stable")
}
if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!require("ggtext", quietly = TRUE)) {
  install.packages("ggtext")
}
 

if (!require("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!require("lattice", quietly = TRUE)) {
  install.packages("lattice")
}
if (!require("rintrojs", quietly = TRUE)) {
  install.packages("rintrojs")
}
if (!require("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}
if (!require("viridis", quietly = TRUE)) {
  install.packages("viridis")
}
if (!require("sp", quietly = TRUE)) {
  install.packages("sp")
}
if (!require("sf", quietly = TRUE)) {
  install.packages("sf")
}
if (!require("gridExtra", quietly = TRUE)) {
  install.packages("gridExtra")
}

library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinydashboardPlus)
library(shinyBS)
library(shinyjs)
library(periscope)
library(INLA)
library(inlabru)
library(fmesher)
library(ggplot2)
library(ggtext)
library(lattice)
library(rintrojs)
library(patchwork)
library(viridis)
library(sp)
library(sf)
library(stringr)
library(dplyr)
library(gridExtra)

header <- dashboardHeader(
  tags$li(
    class = "dropdown",
    tags$style(".main-header {max-height: 200px}"),
    tags$style(".main-header .logo {height: 115px; padding-top: 15px}")
  ),
  title = div(img(src = "logoUV-cropped.svg", height = 85), span(HTML("BAYSPINS"), style = {
    "padding-left: 15px"
  })),
  titleWidth = 550
)
## Sidebar content
sidebar <- dashboardSidebar(
  tags$style(".left-side, .main-sidebar {padding-top: 115px}"),
  minified = FALSE, width = 200,
  sidebarMenu(
    id = "sidebarID",
    menuItem("Introduction",
      tabName = "introduction", icon = icon("exclamation-circle"),
      menuSubItem("Presentation", tabName = "presentation"),
      menuSubItem("Theoretical framework", tabName = "theoryframe"),
      menuSubItem("Glosary", tabName = "glosary")
    ),
    menuItem("Data Simulation", tabName = "datasimulation", icon = icon("database")),
    menuItem("Upload Data", tabName = "dataloading", icon = icon("copy", lib = "glyphicon")),
    menuItem("Model Fitting",
      tabName = "modelanalysis", icon = icon("chart-line"), expandedName = "MAexpand",
      menuSubItem("Independent Model", tabName = "modind"),
      menuSubItem("LGCP Model", tabName = "modlgcp"),
      menuSubItem("Preferential Model", tabName = "modpref"),
      menuSubItem("Mixture Model", tabName = "modmixture")
    )
  )
)

## Body content
body <- dashboardBody(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
  ),
  introjsUI(),
  shinybusy::add_busy_spinner(spin = "atom", margins = c(20, 20), position = "bottom-right", ), # fading-circle
  # useShinyjs(),
  tabItems(
    # UI_Presentation_Section ----
    tabItem(
      tabName = "presentation",
      fluidRow(
        box(
          width = 12, title = h2(strong("Presentation"), align = "center"), status = "primary", solidHeader = TRUE,
          p("Species Distribution Models (SDMs) are a statistical tool that has seen a substantial expansion in its 
                      implementation over the last two decades. Along with their widespread use, the complexity of the data analysed 
                      and the structures of the models used have increased. However, the INLA approach, which is increasingly 
                      used in the field of ecological sciences, has not yet been integrated into an application that can synthesise 
                      the complexity of its code into a user-friendly interface for continuous spatial modelling."),
          p("To overcome this shortcoming, we present a novel application that allows the use of the INLA 
                      methodology for those users who are not very experienced or for those users with experience who prefer a tool 
                      that allows them to carry out an initial analysis quickly, avoiding the process of writing code. 
            The modelling structures that have been implemented in the application are: geostatistical models, log-Gaussian Cox process models, 
            preferential models and mixture models."),
          p("Although this application was intended to provide a simple and versatile tool for the SDM domain, it can be used in any other domain 
            where geostatistical, LGCP, preferential or mixture models are required. Examples of other areas where the application can be used are 
            environmental sciences, spatial econometrics, geospatial health, etc."),
          htmlOutput("Presentation1", fill = TRUE)
        )
      )
    ),
    # Theory tab content
    tabItem(
      tabName = "theoryframe",
      fluidRow(
        column(
          width = 12,
          box(title = h1("Introduction", style = "color: White; text-align: center; text-decoration: underline; font-weight: bold", align = "center"), 
              width=12, status = "primary", solidHeader = TRUE,
          # h1("Introduction", style = "color: DarkRed; text-align: center; text-decoration: underline; font-weight: bold"),
          withMathJax(),
          p("This Shiny app is designed to employ the Bayesian approach and the INLA methodology for solving four distinct spatial models. 
          The Bayesian approach is a statistical paradigm where parameters are treated as random variables, and their distribution functions are estimated through Bayes' theorem."),
          p("The INLA methodology, on the other hand, serves as a deterministic approximation for hierarchical Bayesian models, incorporating an underlying Gaussian field structure for the parameters. 
            The app aims to address the following two types of models:"),
          HTML("
                    
                    <ol>
                    <li> <strong><i>(i) Independent Model (Geostatistical model) </i></strong></li>
                    <li> <strong><i>(i) Log-Gaussian Cox Process </i></strong></li>
                    <li> <strong><i>(i) Preferential Model </i></strong></li>
                    <li> <strong><i>(i) Mixture Model </i></strong></li>
                    </ol>"),
          p("Brief descriptions of these models are provided below."),
          h3("1. Independent model"),
          p("The structure of the geoestatistical model for the response variable \\(y\\) is the following, in which the likelihood function \\(f(\\cdot)\\) of observations \\(y\\) is related to the predictor by means of the link function \\(g(\\cdot)\\):"),
          helpText("\\begin{equation}
                           \\begin{array}{c}
                           y_i \\sim f(y_i|\\boldsymbol\\theta),\\\\
                           g(E(y_i)) = g(\\mu_i) = \\beta_0 + \\mathbf{X}_i\\boldsymbol\\beta + \\sum_k f_k(z_{ik}|\\boldsymbol\\theta_k) + u_i.
                           \\end{array}
                           \\end{equation}"),
          p(paste0(
            "Where \\(\\beta_0\\) is the intercept, ", "\\(\\mathbf{X}\\) ", "are the explanatory variables, ",
            "\\(\\boldsymbol\\beta\\) the parameters, ", "\\(f_k(z_{k})\\)", " are smothing terms ",
            "and ", "\\(u_i\\) ", "is the spatial effect."
          )),
          h3("2. Log-Gaussian Cox process (LGCP)"),
          p(paste0("The second model we consider is the log-Gaussian Cox Process (LGCP) model. This model facilitates the assessment of the process that generates the locations of a sample.
            The specific positions of observed events are influenced by an underlying spatial process, typically represented through an intensity function, ", "\\lambda(s_i)",
            ". This function denotes the average event count per spatial unit and can be shaped by various factors, including covariates and other random effects, 
            whether structured or unstructured:
            ")),
          helpText("\\begin{equation}
                           \\begin{array}{c}
                           s_i \\sim LGCP(s_i|\\boldsymbol\\theta),\\\\
                           \\log(E(\\lambda_i)) = \\log(\\lambda_i)= \\beta_0 + \\mathbf{X}_i \\boldsymbol\\beta + \\sum_k f_k(z_{ik}|\\boldsymbol\\theta_k) + \\alpha \\cdot u_i.
                           \\end{array}
                           \\end{equation}"),
          p(paste0("where we have defined the latent structure for the intensity function with linear ", "\\(\\boldsymbol\\beta\\)", ", non-linear ", "\\(f_k(z_{k})\\)", 
                   " effects and one spatially structured random effect ", "\\(u_i\\)", "." )),
          h3("3. Preferential model"),
          p("The preferential model combines two likelihoods, one related to the geostatistical data (response variable or mark values) and another related to the poin pattern of the locations.
          Between the two likelihoods we have one or more shared components (for us only the spatial term is shared), 
          what implies that their estimates would be done though both likelihoods. The structure of the model can be written as:"),
          helpText("\\begin{equation}
                            \\begin{array}{c}
                           y_i \\sim f(y_i|\\boldsymbol\\theta),\\\\
                           s_i \\sim LGCP(s_i|\\boldsymbol\\theta),\\\\
                           g(E(y_i)) = g(\\mu_i)= \\beta_0 + \\mathbf{X}_i\\boldsymbol\\beta + \\sum_k f_k(z_{ik}|\\boldsymbol\\theta_k) + u_i,\\\\
                           \\log(E(\\lambda_i)) = \\log(\\lambda_i)= \\beta_0' + \\mathbf{X}_i \\boldsymbol\\beta' + \\sum_k f_k(z_{ik}) + \\alpha \\cdot u_i.
                           \\end{array}
                           \\end{equation}"),
          p("Where the intercept for the geostatistical layer is \\(\\beta_0\\) and 
                    \\(\\beta_0'\\) is the one for the point process layer,
                    \\(\\mathbf{X}\\) are the explanatory variables, \\(\\boldsymbol\\beta\\) the parameters of the geostatistical layer and
                    \\(\\boldsymbol\\beta'\\) the parameters of the point process layer. Meanwhile \\(u_i\\) is the spatial effect and
                    \\(\\alpha\\) is a internal shared parameter, a scale term."),
          h3("4. Mixture model"),
          p(paste("In the Mixture Model we assume that the wole data is a mixture of different data sets coming from several sampling structures.", 
                  "Hence, this model integrates a geostatistical process for analyzing the values of the marks along with multiple Log-Gaussian Cox Processes (LGCP) 
                  for each dataset originating from distinct preferential sampling processes:")),
          helpText("\\begin{equation}
                            \\begin{array}{c}
                           y_i \\sim f(y_i|\\boldsymbol\\theta),\\\\
                           s^k_i \\sim LGCP(s^k_i|\\boldsymbol\\theta^k),\\\\
                           g(E(y_i)) = g(\\mu_i)= \\beta_0 + \\mathbf{X}_i\\boldsymbol\\beta + \\sum_j f_j(z_{ij}|\\boldsymbol\\theta_j) + u_i,\\\\
                           \\log(E(\\lambda^k_i)) = \\log(\\lambda^k_i)= \\beta^k_0 + \\mathbf{X}_i \\boldsymbol\\beta^k + \\sum_j f_j(z_{ij}|\\boldsymbol\\theta^k_j) + \\alpha^k \\cdot u_i.
                           \\end{array}
                           \\end{equation}"),
          p("Where the parameters have the same meaning as in the preferential model, but with the particularity that in the LGCP process we have the label \\(k\\) 
            to differentiate the inference parameters for the different LGCPs associated to the corresponding group of samples."),
          br(),
          p(
            "Therefore, the application allows to simulate an example data to test the inference results controlling its parameters and also allows data loading and modeling, 
                    according to the model structures showed previosly. The simulation is done in the", span("Simulation Data", style = "font-style: italic"),
            "section, in which is possible to control the parameters that regulates the simulation of the covariables, intercept, spatial effect and the geostatistical data output. 
                    In te next section", span("(Model Fitting)", style = "font-style: italic"), "is possible to shape the inference through some buttons implemented in the UI, which make it possible to control the basic prior parameters for intercept, explanatory variables, spatial effect and the likelihood function, 
                    as well as for specific INLA controll elements. These controll items are defined in", span("Advance INLA config.", style = "font-style: italic"), "where we can define the modes for the hyperparameters, the approximation type and the integration strategy."
          )),
        )
      )
    ),
    tabItem(
      tabName = "glosary",
      fluidRow(
        column(
          width = 12,
          box(
            width = 12, title = h1("Glosary", style = "color: White; text-align: center; text-decoration: underline; font-weight: bold", align = "center"),
            status = "primary", solidHeader = TRUE,
          p("Since this application is intended for all kind of users, a glossary is included containing a brief selection of terms, with some explanations, as they are used along the whole application.
                    This  has two ordered levels, one by field and the second ordered by the terms belonging to that specific field."),
          tags$div(
            HTML("
                    <ul> <h3><strong>Bayesian inference:</strong></h3>
                    <ol>
                    <li> <strong><i>Bayes' theorem (conditional probability distribution):</i></strong> A fundamental theorem in probability theory that describes the probability of an event based on prior knowledge of conditions that might be related to the event. It is named after the Reverend Thomas Bayes and provides a way to update probabilities as new evidence or information becomes available. </li>
                    <li> <strong><i>Hierarchical model:</i></strong> A statistical model that includes nested levels of parameters or random effects. In hierarchical modeling, data are assumed to have a structure that allows for variation at multiple levels, and it is particularly useful when dealing with grouped, clustered or high complex data.</li>
                    <li> <strong><i>Model parameter(s):</i></strong> The parameters in a statistical model are the coefficients or constants that define the model. They are the values that need to be estimated from the data, and they influence the shape and characteristics of the model. </li>
                    <li> <strong><i>Likelihood:</i></strong> The likelihood function represents the probability of observing the given data given a particular set of parameter values in a statistical model. It is a key component in Bayesian and frequentist statistical inference. </li>
                    <li> <strong><i>Link function:</i></strong> In generalized linear models (GLM), the link function describes the relationship between the linear predictor and the mean of the distribution function. It connects the linear combination of predictors to the parameters of the distribution. </li>
                    <li> <strong><i>Marginal distribution(s):</i></strong> The distribution of one or more variables in a subset of a larger set of variables, obtained by integrating or summing over the other variables. It provides information about the probability distribution of a subset of variables, ignoring the values of the other variables. </li>
                    <li> <strong><i>Mean:</i></strong> A measure of central tendency that represents the average value of a set of data points. It is calculated by summing up all values and dividing by the number of observations. </li>
                    <li> <strong><i>Mode:</i></strong> The mode of a distribution is the value that occurs most frequently. In a dataset, it is the data point with the highest frequency. </li>
                    <li> <strong><i>Prior distribution:</i></strong> The prior distribution represents the initial beliefs or knowledge about the parameters before observing the data. It is an essential component in Bayesian statistics and is updated to the posterior distribution through Bayes' theorem. </li>
                    <li> <strong><i>Prior predictive distribution:</i></strong> The prior predictive distribution incorporates uncertainty from the prior distribution to generate predictions before observing any data. It represents the range of possible outcomes based on the initial beliefs about the parameters. </li>
                    <li> <strong><i>Posterior distribution:</i></strong> In Bayesian statistics, the posterior distribution is the updated probability distribution of the parameters after taking into consideration the observed data and the prior distribution. It is obtained using Bayes' theorem. </li>
                    <li> <strong><i>Posterior predictive distribution:</i></strong> This distribution combines the information from the likelihood and the posterior distribution to make predictions for future observations. It reflects the uncertainty in predictions by considering both the uncertainty in parameters and the variability in the data. </li>
                    </ol>
                    </ul>
                    
                    <ul> <h3><strong>INLA basic elements:</strong></h3>
                    <ol>
                    <li> <strong><i>Aproximation strategy:</i></strong> INLA provides three distinct methods for approximating non-Gaussian data. These methods include a straightforward Gaussian approximation, a simplified Laplace approach—which is essentially a Gaussian approximation refined through spline correction—and the Laplace approximation, involving the utilization of two nested Gaussian approximations. </li>
                    <li> <strong><i>Integration strategy:</i></strong> INLA provides three different methods for the integration of the hyperparameters when the mode set to 'classic'. These three strategies are the grid, the center composite design and empirical bayes.</li>
                    <li> <strong><i>Latent field:</i></strong> A latent field refers to an unobserved, underlying random field that is part of the statistical model. In INLA, the latent field represents the unobservable values or processes of interest. </li>
                    <li> <strong><i>Gaussian Field (GF):</i></strong> A Gaussian field is a type of random field where the distribution of the values at any finite collection of points follows a multivariate normal (Gaussian) distribution. Gaussian fields are commonly used in spatial statistics and geostatistics to model spatial variability or dependence. </li>
                    <li> <strong><i>Gaussian Markov Random Field (GMRF):</i></strong> A Gaussian Markov Random Field (GMRF) is a type of random field where the conditional distribution of each random variable, given the values of all other variables, is Gaussian and the variables are associated with a graph structure. GMRFs are used in spatial statistics and Bayesian modeling to represent spatial dependence in a computationally efficient manner, often making use of sparse precision matrices. The 'Markov' aspect implies that the conditional independence structure is defined by a Markov property on the underlying graph. </li>
                    </ol>
                    </ul>
                    
                    <ul> <h3><strong>INLA spatial approach:</strong></h3>
                    <ol>
                    <li> <strong><i>Finite Element Method (FEM):</i></strong> A numerical technique for finding approximate solutions to boundary value problems for partial differential equations. It involves subdividing a complex system into smaller, simpler parts called finite elements. The behavior of each element is then described by a set of mathematical equations, and the solutions are combined to approximate the behavior of the entire system. </li>
                    <li> <strong><i>Matèrn covariance function:</i></strong> A mathematical function commonly used in spatial statistics and geostatistics to model the correlation or covariance between spatial data points. The Matérn covariance function is characterized by a smoothness parameter that influences the smoothness of the resulting spatial process. </li>
                    <li> <strong><i>Mesh:</i></strong> In the context of the Integrated Nested Laplace Approximation (INLA), a mesh refers to a discretization of the spatial domain. The mesh is used to represent the spatial structure of the data and is essential for the computational aspects of INLA. </li>
                    <li> <strong><i>Stochastic Partial Differential Equation (SPDE):</i></strong> A type of partial differential equation where one or more of the parameters or coefficients are subject to random variations. SPDEs are often used in the modeling of random fields and spatial processes. They provide a framework for describing the behavior of systems that exhibit both spatial variability and randomness. The use of SPDEs is common in Bayesian spatial statistics and geostatistics. </li>
                    </ol>
                    </ul>
                         
                         ")
            )
          )
        )
      )
    ),

    # UI_Simulation_Section ----
    tabItem(
      tabName = "datasimulation",
      conditionalPanel(
        condition = "input.infoSim",
        fluidRow(
          column(
            width = 12,
            mainPanel(
              box(title = h1("Summary", style = "color: White; text-align: center; text-decoration: underline; font-weight: bold", align = "center"), 
                  width=12, status = "primary", solidHeader = TRUE,
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
                             The first rules the randomness of the whole geostatistic process simulation, meanwhile the latter controls the specific randomness of the spatial effect simulation.")
                ),
                br(),
                tags$div(
                  HTML(
                    "With respect to the simulation of the spatial effect, given that it is characterized by the \\(u_i \\sim N(\\mathbf{0}, \\boldsymbol\\Sigma) \\) expression, 
                          where the variace matrix is determined by a Matérn covariance function \\(\\Sigma=\\sigma \\cdot C_{\\nu}(d, \\rho)\\) with two parameters: spatial range \\(\\rho\\) and standard marginal deviation \\(\\sigma\\). Two numerical inputs are availabled in the UI to stablish the values for such parameters."
                    )
                ),
                br(),
                tags$div(
                  HTML(
                    "For the covariate we can define the formula (<i>Covariate formula</i>) and choose the effect type:
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
                          (iii) finally, it is checked if any resulting value exceeds the indicated limits for the effect, so that if it does, the two previous steps are repeated until the condition is met. The <i>predictor coefficients</i> allow to define the itercept and the linear effect, if this was chosen for the covariate.</i>"
                  )
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
                )
              ),
              br()
            )
          )
        )
      ),
      fluidRow(
        column(
          style = "padding-bottom:15px",
          width = 12,
          introBox(
            bsButton("parametersMapSim",
              label = "Geostatistical parameters",
              type = "toggle",
              value = TRUE,
              icon = icon("drafting-compass"),
              style = "primary"
            ),
            bsButton("parametersSampleSim",
              label = "Sample parameters",
              type = "toggle",
              value = TRUE,
              icon = icon("spinner", class = "spinner-box"),
              style = "primary"
            ),
            bsButton("infoSim",
              label = "Summary information",
              type = "toggle",
              value = FALSE,
              icon = icon("info"),
              style = "warning"
            )
          )
        )
      ),
      fluidRow(
        conditionalPanel(
          condition = "input.parametersMapSim|input.parametersSampleSim",
          box(
            width = 3,
            conditionalPanel(
              condition = "input.parametersMapSim",
              box(
                width = 12, status = "info",
                title = "Geostatistical parameters",
                solidHeader = TRUE,
                collapsible = TRUE,
                actionButton("makeSim",
                  label = "Simulate",
                  # type = "action",
                  # block = TRUE,
                  width = "100%",
                  # style = "unite",
                  icon = icon("jedi-order")
                ),
                sliderInput("limlattice",
                  label = "Analysis area:",
                  min = -5,
                  max = 5,
                  value = c(0, 1)
                ),
                sliderInput("lengthlattice",
                  label = "Vector dimension (grid):",
                  min = 10,
                  max = 200,
                  value = 40
                ),
                sliderInput("dimmap",
                  label = "Map resolution:",
                  min = 50,
                  max = 300,
                  value = 150
                ),
                numericInput("seedGlobal",
                  label = "Global seed:",
                  value = 1
                ),
                box(
                  width = 12, title = "Spatial parameters",
                  status = "info", solidHeader = TRUE, collapsible = TRUE,
                  numericInput("seedSP",
                    label = "Spatial effect seed:",
                    value = 1
                  ),
                  numericInput("sigma.matern", # width = '45%',
                    label = HTML(paste0("Matérn stdev. (\U03C3", tags$sub("0"), "):")),
                    value = 0.5
                  ),
                  numericInput("range.matern", # width = '45%',
                    label = HTML(paste0("Matérn range (\U03C1", tags$sub("0"), "):")),
                    value = 0.2
                  )
                ),
                textInput("bathymetry.formula",
                  label = "Covariate formula (x, y)",
                  value = "0.5 + (x-min(x)) + 0.8*(y-mean(y))**2 + rnorm(length(x), 0, 0.01)",
                  placeholder = "sqrt(abs(x+y))+x**2+log(y-min(y)+1)"
                ),
                selectInput("typeBathyeff",
                  label = "Types of covariate effects",
                  choices = c(
                    "Linear" = "lin", "Random Walk (1o)" = "rw1",
                    "Random Walk (2o)" = "rw2",
                    "Custom function" = "customfunction"
                  )
                ),
                conditionalPanel(
                  condition = "input.typeBathyeff=='rw1'",
                  box(
                    status = "info", title = "Parameters RW1", width = 12, solidHeader = TRUE,
                    numericInput("prec.rw1", label = "Precision (RW1)", value = 10, min = 0),
                    numericInput("init.rw1", label = "Init value (RW1)", value = 0.2),
                    numericInput("nknots.rw1", label = "Number of knots (RW1)", value = 60),
                    textInput("minmax.effect.rw1",
                      label = "Extreme values",
                      value = "0, 5",
                      placeholder = "sqrt(abs(x+y))+x**2+log(y-min(y)+1)"
                    )
                  )
                ),
                conditionalPanel(
                  condition = "input.typeBathyeff=='rw2'",
                  box(
                    status = "info", title = "Parameters RW2", width = 12, solidHeader = TRUE,
                    numericInput("prec.rw2", label = "Precision (RW2)", value = 10, min = 0),
                    numericInput("init.rw2", label = "Init value (RW2)", value = 0.2),
                    numericInput("nknots.rw2", label = "Number of knots (RW2)", value = 60),
                    textInput("minmax.effect.rw2",
                      label = "Extreme values",
                      value = "0, 5",
                      placeholder = "sqrt(abs(x+y))+x**2+log(y-min(y)+1)"
                    )
                  )
                ),
                conditionalPanel(
                  condition = "input.typeBathyeff=='customfunction'",
                  box(
                    status = "info", title = "Custom function", width = 12, solidHeader = TRUE,
                    textInput("formula.eff.bathymetry",
                      label = "Bathymetry effect formula (b)",
                      value = "0.3 + 0.5*b + 1.5*(b-mean(b))**3 + rnorm(length(b), 0, 0.01)",
                      placeholder = "0.1+b+b**2"
                    )
                  )
                ),
                textInput("beta",
                  label = HTML(paste0(
                    "Predictor coefficients (\U03B2", tags$sub("0"),
                    ", \U03B2", tags$sub("1"), ", ...", "):"
                  )),
                  value = "0,1",
                  placeholder = "0,1"
                ),
                numericInput("var",
                  label = HTML(paste0("Distribution variance (\U03C3", tags$sup("2"), "):")),
                  value = 0.2
                ),
                selectInput("datadistributionSim",
                  label = "Data distribution",
                  choices = c("Gamma" = "gamma", "Gaussian" = "gaussian") #, "Bernoulli" = "bernoulli")
                )
              )
            ),
            conditionalPanel(
              condition = "input.parametersSampleSim",
              box(
                width = 12, status = "info",
                title = "Sample parameters",
                solidHeader = TRUE,
                collapsible = TRUE,
                footer = "The (iid) are the parameters for the independent sampling, 
                                    and the (ps) are the parameters for the preferential sampling.",
                actionButton("makeSample",
                  label = "Sample",
                  width = "100%",
                  icon = icon("jedi-order")
                ),
                numericInput("seedSampleR",
                  label = "Seed (iid):",
                  value = 1
                ),
                numericInput("niid.samples", "Number of samples (iid)",
                  value = 120
                ),
                numericInput("seedSampleP",
                  label = "Seed (ps):",
                  value = 1
                ),
                numericInput("nps.samples", "Number of samples (ps)",
                  value = 120
                ),
                numericInput("r.scale", "Scale factor (ps)",
                  value = 7
                )
              )
            ),
            conditionalPanel(
              condition = "input.parametersSampleSim",
              box(
                width = 12, status = "info",
                title = "Mixture sample parameters",
                solidHeader = TRUE,
                collapsible = TRUE,
                footer = "We can simulate two preferential samples or simulate one as an independent sample by seting the scale value equal to zero.",
                actionButton("makeSampleMixture",
                             label = "Sample",
                             width = "100%",
                             icon = icon("jedi-order")
                ),
                numericInput("seedSample.mixture1",
                             label = "Seed (m1):",
                             value = 1
                ),
                numericInput("n.mixture1.samples", "Number of samples (m1)",
                             value = 60
                ),
                numericInput("r.scale.mixture1", "Scale factor (m1)",
                             value = -7
                ),
                numericInput("seedSample.mixture2",
                             label = "Seed (m2):",
                             value = 1
                ),
                numericInput("n.mixture2.samples", "Number of samples (m2)",
                             value = 60
                ),
                numericInput("r.scale.mixture2", "Scale factor (m2)",
                             value = 7
                )
              )
            )
          )
        ),
        conditionalPanel(
          condition = "input.makeSim>=1",
          box(
            id = "resultbox", width = 9, collapsible = FALSE,
            box(
              id = "mesh.speff.box", width = 12, status = "info",
              title = "Mesh and Spatial effect", solidHeader = TRUE,
              collapsible = TRUE, collapsed = FALSE,
              column(
                width = 6,
                downloadFileButton(
                  id = "ggplotMeshSim",
                  downloadtypes = c("png", "csv", "txt"),
                  hovertext = "Download image and data"
                ),
                downloadablePlotUI(
                  id = "ggplotMeshSim",
                  list("png", "txt", "tsv"),
                  btn_valign = "bottom",
                  btn_halign = "right",
                  btn_overlap = TRUE
                )
              ),
              column(
                width = 6,
                downloadFileButton(
                  id = "ggplotSpSim",
                  downloadtypes = c("png", "csv", "txt"),
                  hovertext = "Download image and data"
                ),
                downloadablePlotUI(
                  id = "ggplotSpSim",
                  list("png", "txt", "tsv"),
                  btn_valign = "bottom",
                  btn_halign = "right",
                  btn_overlap = TRUE
                )
              )
            ),
            box(
              id = "bathymetry.bath.effect", width = 12, status = "info",
              title = "Covariate map and its effect", solidHeader = TRUE,
              collapsible = TRUE, collapsed = FALSE,
              column(
                width = 6,
                downloadFileButton(
                  id = "ggplotBatChartSim",
                  downloadtypes = c("png", "csv", "txt"),
                  hovertext = "Download image and data"
                ),
                downloadablePlotUI(
                  id = "ggplotBatChartSim",
                  list("png", "txt", "tsv"),
                  btn_valign = "bottom",
                  btn_halign = "right",
                  btn_overlap = TRUE
                )
              ),
              column(
                width = 6,
                downloadFileButton(
                  id = "ggplotBatEffSim",
                  downloadtypes = c("png", "csv", "txt"),
                  hovertext = "Download image and data"
                ),
                downloadablePlotUI(
                  id = "ggplotBatEffSim",
                  list("png", "txt", "tsv"),
                  btn_valign = "bottom",
                  btn_halign = "right",
                  btn_overlap = TRUE
                )
              )
            ),
            box(
              id = "abundance.batheffect.maps", width = 12, status = "info",
              title = "Covariate effect map and Response map", solidHeader = TRUE,
              collapsible = TRUE, collapsed = FALSE,
              column(
                width = 6,
                downloadFileButton(
                  id = "ggplotBatEffChartSim",
                  downloadtypes = c("png", "csv", "txt"),
                  hovertext = "Download image and data"
                ),
                downloadablePlotUI(
                  id = "ggplotBatEffChartSim",
                  list("png", "txt", "tsv"),
                  btn_valign = "bottom",
                  btn_halign = "right",
                  btn_overlap = TRUE
                )
              ),
              column(
                width = 6,
                downloadFileButton(
                  id = "ggplotAbundanceSim",
                  downloadtypes = c("png", "csv", "txt"),
                  hovertext = "Download image and data"
                ),
                downloadablePlotUI(
                  id = "ggplotAbundanceSim",
                  list("png", "txt", "tsv"),
                  btn_valign = "bottom",
                  btn_halign = "right",
                  btn_overlap = TRUE
                )
              )
            ),
            conditionalPanel(
              condition = "input.makeSample>=1",
              box(
                id = "sample.maps", width = 12, status = "info",
                title = "Independent and preferential samples", solidHeader = TRUE,
                collapsible = TRUE, collapsed = FALSE,
                column(
                  width = 6,
                  downloadFileButton(
                    id = "ggplotRandomSampleSim",
                    downloadtypes = c("png", "csv", "txt"),
                    hovertext = "Download image and data"
                  ),
                  downloadablePlotUI(
                    id = "ggplotRandomSampleSim",
                    list("png", "txt", "tsv"),
                    btn_valign = "bottom",
                    btn_halign = "right",
                    btn_overlap = TRUE
                  )
                ),
                column(
                  width = 6,
                  downloadFileButton(
                    id = "ggplotPrefSampleSim",
                    downloadtypes = c("png", "csv", "txt"),
                    hovertext = "Download image and data"
                  ),
                  downloadablePlotUI(
                    id = "ggplotPrefSampleSim",
                    list("png", "txt", "tsv"),
                    btn_valign = "bottom",
                    btn_halign = "right",
                    btn_overlap = TRUE
                  )
                )
              )
            ),
            conditionalPanel(
              condition = "input.makeSampleMixture>=1",
              box(
                id = "sample.mixture.maps", width = 12, status = "info",
                title = "Mixture samples", solidHeader = TRUE,
                collapsible = TRUE, collapsed = FALSE,
                column(
                  width = 12,
                  downloadFileButton(
                    id = "ggplotMixtureSampleSim",
                    downloadtypes = c("png", "csv", "txt"),
                    hovertext = "Download image and data"
                  ),
                  downloadablePlotUI(
                    id = "ggplotMixtureSampleSim",
                    list("png", "txt", "tsv"),
                    btn_valign = "bottom",
                    btn_halign = "right",
                    btn_overlap = TRUE
                  )
                )
              )
            )
          )
        )
      )
    ),

    # UI_LoadingData_Section ----
    
    tabItem(
      tabName = "dataloading",
      conditionalPanel(
        condition = "input.infoSimLoad",
        fluidRow(
          column(
            # style= 'padding-bottom:15px',
            width = 12,
            mainPanel(
              box(title = h1("Summary", style = "color: White; text-align: center; text-decoration: underline; font-weight: bold", align = "center"), 
                  width=12, status = "primary", solidHeader = TRUE,
              withMathJax(),
              tags$div(
                HTML("In this section we briefly explain how data reading works in this app. For this purpose, there are two elements in the UI: 
                               <ol>
                               <li><strong>Main Analysis Data Frame:</strong> it allows to read the loaded data frame, but the two first columns must be the \\((x,y)\\) coordinates, the third the observational variable \\(y\\) and the others the explanatory variables \\(\\mathbf(X)\\).
                               <li><strong>Auxiliary Covariate Data Frame:</strong> it allows to read a data frame realated to some covariate and using it for the prediction grid. The two first columns must be the \\((x,y)\\) and the third the covariate.</li>
                              </ol>
                               Once any data frame is loaded its table and a plot will appear in the interface.  
                               ")
              )),
              br()
            )
          )
        )
      ),
      fluidRow(
        column(
          style = "padding-bottom:15px",
          width = 12,
          introBox(
            bsButton("LoadingPanel",
              label = "Loading Panel",
              type = "toggle",
              value = TRUE,
              icon = icon("file-export"),
              style = "primary"
            ),
            bsButton("infoSimLoad",
              label = "Summary information",
              type = "toggle",
              value = FALSE,
              icon = icon("info"),
              style = "warning"
            )
          )
        )
      ),
      fluidRow(
        conditionalPanel(
          condition = "input.LoadingPanel",
          box(
            width = 4,
            fileInput("file.uploadData",
              label = "Main Analysis Data Frame:",
              placeholder = "No file selected (.csv or .rds)"
            ),
            # sliderInput("AspectRatioSample", label="Aspect Ratio", min=1, max=10, value=2),
            fileInput("file.uploadDataRaster",
              label = "Auxiliary Covariate Data Frame:",
              placeholder = "No file selected (.csv or .rds)"
            )
          )
        ),
        conditionalPanel(
          condition = "output.fileUploadedSample|output.fileUploadedRaster",
          box(
            width = 8, id = "TablePanelRaster",
            # tableOutput("content")
            conditionalPanel(
              "output.fileUploadedSample",
              column(
                width = 12,
                downloadableTableUI("table.read.sample",
                  downloadtypes = c("csv", "tsv")
                )
              )
            ),
            conditionalPanel(
              "output.fileUploadedSample",
              column(
                width = 12,
                downloadFileButton(
                  id = "quilt.plot.sample",
                  downloadtypes = c("png"),
                  hovertext = "Download image"
                ),
                plotOutput(outputId = "main_plot2", width = "100%", height = "auto")
                # downloadablePlotUI(id = "quilt.plot.sample",
                #                    list("png"),
                #                    btn_valign = "bottom",
                #                    btn_halign = "right",
                #                    btn_overlap = TRUE)
              )
            ),
            conditionalPanel(
              "output.fileUploadedRaster",
              column(
                width = 12,
                downloadableTableUI("table.read.raster",
                  downloadtypes = c("csv", "tsv")
                )
              )
            ),
            conditionalPanel(
              "output.fileUploadedRaster",
              column(
                width = 12,
                downloadFileButton(
                  id = "quilt.plot.raster",
                  downloadtypes = c("png"),
                  hovertext = "Download image"
                ),
                plotOutput(outputId = "main_plot3", width = "100%", height = "auto")
              )
            )
          )
        )
      )
    ),

    # UI_Independent_Section ----
    tabItem(
      tabName = "modind", withMathJax(),
      fluidRow(
        box(
          width = 3,
          actionButton("fitInd", label = "Fit Model", width = "100%", icon = icon("jedi-order")), br(), br(),
          box(id = "SaveIndModel", width = 12, title = "Save model data and code", closable = FALSE,
              status = "info", collapsible = TRUE, collapsed = TRUE, solidHeader = TRUE,
              prettyToggle(
                inputId = "saveInd",
                label_on = "Save", 
                label_off = "Do not save"
              ),
              textInput(inputId = "IndSavePath", label = "Absolute path to save the model", value = "", placeholder = "C:/example/example/"),
              textInput(inputId = "IndSaveName", label = "Filename", value = "Ind_example.RDS", placeholder = "Ind_example.RDS")
          ),
          selectInput("IndDataSimulatedLoaded",
                      label = "Select data to be analyzed",
                      choices = list("Simulated" = "sim", "Loaded" = "load"), selected = "sim"
          ),
          selectInput("IndRasterSPDE",
                      label = "Use auxilary data",
                      choices = list("Yes" = "raster", "No" = "solvecov"), selected = "solvecov"
          ),
          conditionalPanel(
            condition = "input.IndRasterSPDE=='raster'",
            selectInput("IndRasterPred",
                        label = "Auxiliary data functionality",
                        selected = "SPDEraster",
                        choices = list("Additional data" = "SPDEraster",
                                       "Fixing prediction points" = "rasterpred")
            )
          ),
          conditionalPanel(
            condition = "input.IndRasterPred=='SPDEraster' || input.IndRasterSPDE=='solvecov' || input.IndDataSimulatedLoaded=='sim'",
            sliderInput("IndSPDErasterInput", "Predictive grid dimensionality:",
                        min = 50, max = 200, value = 100
            )
          ),
          sliderInput("dimIndmap",
                      label = "Predictive map resolution (spatial effect):",
                      min = 100, max = 300, value = 150
          ),
          box(
            id = "IndMeshCustom", width = 12, title = "Custom Mesh", closable = FALSE,
            status = "warning", collapsible = TRUE, collapsed = TRUE, solidHeader = TRUE,
            materialSwitch(
              inputId = "IndCustomMesh", label = "Allow custom mesh",
              value = FALSE, status = "primary"
            ),
            conditionalPanel(
              condition = "input.IndCustomMesh",
              actionButton("buildIndMesh", label = "Build Mesh", class = "button"),
              radioGroupButtons(
                inputId = "boundaryIndMesh",
                label = "Outer Boundary Mesh",
                choices = c("Default", "Non-Convex" = "IndMeshnonconvex", "Custom" = "IndMeshcustomboundary"),
                individual = TRUE,
                checkIcon = list(
                  yes = tags$i(class = "fa fa-circle", style = "color: steelblue"),
                  no = tags$i(class = "fa fa-circle-o", style = "color: steelblue")
                )
              ),
              conditionalPanel(
                condition = "input.boundaryIndMesh=='IndMeshnonconvex'",
                numericInput(inputId = "curvatureIndMesh", label = "Curvature Radious", value = 0),
                numericInput(inputId = "resolutionIndMesh", label = "Resolution", value = 40, min = 40)
              ),
              conditionalPanel(
                condition = "input.boundaryIndMesh=='IndMeshcustomboundary'",
                fileInput(
                  inputId = "shapefileIndMesh", label = "Outer Boundary",
                  accept = c(".csv", ".rds"), multiple = TRUE
                )
              ),
              radioGroupButtons(
                inputId = "interiorIndMesh",
                label = "Inner Boundary Mesh",
                choices = c("Default", "Non-Convex" = "interiorIndMeshnonconvex", "Custom" = "interiorIndMeshcustomboundary"),
                individual = TRUE,
                checkIcon = list(
                  yes = tags$i(class = "fa fa-circle", style = "color: steelblue"),
                  no = tags$i(class = "fa fa-circle-o", style = "color: steelblue")
                )
              ),
              conditionalPanel(
                condition = "input.interiorIndMesh=='interiorIndMeshnonconvex'",
                numericInput(inputId = "interiorcurvatureIndMesh", label = "Curvature Radious", value = 0),
                numericInput(inputId = "interiorresolutionIndMesh", label = "Resolution", value = 40, min = 40)
              ),
              conditionalPanel(
                condition = "input.interiorIndMesh=='interiorIndMeshcustomboundary'",
                fileInput(
                  inputId = "interiorshapefileIndMesh", label = "Inner Boundary",
                  accept = c(".csv", ".rds"), multiple = TRUE
                )
              ),
              radioGroupButtons(
                inputId = "selectionIndMesh",
                label = "Mesh density",
                choices = c("Quantile location" = "qlocation", "Edge length" = "edgelength"),
                individual = TRUE,
                checkIcon = list(
                  yes = tags$i(class = "fa fa-circle", style = "color: steelblue"),
                  no = tags$i(class = "fa fa-circle-o", style = "color: steelblue")
                )
              ),
              conditionalPanel(
                condition = "input.selectionIndMesh=='qlocation'",
                numericInput(
                  inputId = "IndMeshQloc", label = "Quantile Location",
                  width = "100%", value = 0.03, min = 0, max = 1, step = 0.01
                )
              ),
              conditionalPanel(
                condition = "input.selectionIndMesh=='edgelength'",
                numericInput(inputId = "EdgeLengthIndMesh", label = "Edge Length", width = "100%", value = 1, min = 0)
              )
            )
          ),
          box(
            id = "IndFamily", width = 12, title = "Family and its hyperparameters",
            status = "info", solidHeader = TRUE, collapsible = TRUE,
            selectInput("SelectIndFamily",
                        label = "Family distribution",
                        choices = list("Bernoulli" = "bernoulli", "Binomial" = "binomial", "Poisson" = "poisson", "Beta" = "beta", "Gaussian" = "gaussian", "Gamma" = "gamma", "Exponential" = "exponential"),
                        selected = "gaussian"
            ),
            conditionalPanel(
              condition = "input.SelectIndFamily=='beta' || input.SelectIndFamily=='gaussian' || input.SelectIndFamily=='gamma'",
              radioGroupButtons(
                inputId = "autocustomIndFamily",
                label = "Hyperparameters",
                choices = list("Default" = "default", "Custom" = "custom"),
                status = "primary"
              ),
              conditionalPanel(
                condition = "input.autocustomIndFamily=='custom'",
                selectInput(inputId="IndFamilyPriorKind",
                            label="Prior distribution",
                            choices = list("Log-Gamma"="loggamma", "PC-prior"="pc", "Uniform"="unif", "Flat Uniform"="unifflat"),
                            selected="loggamma"
                ),
                textInput("IndFamilyHyper",
                          label = "Prior parameters",
                          value = "0.01, 0.01", placeholder = "0.01, 0.01"
                )
              )
            )
          ),
          box(
            id = "IndResponseVariable", width = 12, title = "Response variable",
            status = "info", solidHeader = TRUE, collapsible = TRUE,
            uiOutput("UserResponseInd")
          ),
          box(
            id = "IndDefaultComponents", width = 12, title = "Intercept and Explanatory Variables",
            status = "info", solidHeader = TRUE, collapsible = TRUE,
            checkboxGroupInput(
              inputId = "DefaultComponentsInd",
              label = "Default Components",
              choices = c("Intercept", "Spatial Effect"),
              selected = c("Intercept", "Spatial Effect")
            ),
            uiOutput("checkBoxIndDataFrame"),
            uiOutput("SelectIndEffectCov"),
            conditionalPanel(
              condition = "input.IndDataSimulatedLoaded=='sim'",
              uiOutput("SelectSimIndEffectCov1")
            ),
            conditionalPanel(
              condition = "input.IndDataSimulatedLoaded=='load'",
              uiOutput("SelectLoadIndEffectCovFactorPred")
            )
          ),
          conditionalPanel(
            condition = "input.DefaultComponentsInd.includes('Spatial Effect')",
            box(
              id = "IndPriorSpatial", width = 12, title = "Spatial Effect",
              status = "info", solidHeader = TRUE, collapsible = TRUE,
              box(
                title = "Joit prior preview", width = 12, status = "info", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                popify(textInput(inputId = "Indrangebasepriorprev", label = "Base Prior Range", value = "0.01,5,200,0.5,0,1", width = NULL, placeholder = "0.01,5,200,0.5,0,1"),
                       title = "Base Prior Range",
                       content = "The <b>meaning of the values</b> according to the order of presentation are: the lower value of the axis, the upper value, the number of breaks, the value for rho0, the mean of the normal and the standard deviation of the normal.",
                       placement = "right"
                ),
                popify(textInput(inputId = "Indsigmabasepriorprev", label = "Base Prior Sigma", value = "0.01,5,200,1,0,1", width = NULL, placeholder = "0.01,5,200,1,0,1"),
                       title = "Base Prior Stdev. Deviation",
                       content = "The <b>meaning of the values</b> according to the order of presentation are: the lower value of the axis, the upper value, the number of breaks, the value for rho0, the mean of the normal and the standard deviation of the normal.",
                       placement = "right"
                ),
                popify(textInput(inputId = "Indrangepcpriorprev", label = "PC-Prior Range", value = "0.01,5,200,0.5,0.5", width = NULL, placeholder = "0.01,5,200,0.5,0.5"),
                       title = "PC-Prior Range",
                       content = "The <b>meaning of the values</b> according to the order of presentation are: the lower value of the axis, the upper value, the number of breaks, the value for rho0, the mean of the normal and the standard deviation of the normal.",
                       placement = "right"
                ),
                popify(textInput(inputId = "Indsigmapcpriorprev", label = "PC-Prior Sigma", value = "0.01,5,200,1,0.5", width = NULL, placeholder = "0.01,5,200,1,0.5"),
                       title = "PC-Prior Stdev. Deviation",
                       content = "The <b>meaning of the values</b> according to the order of presentation are: the lower value of the axis, the upper value, the number of breaks, the value for rho0, the mean of the normal and the standard deviation of the normal.",
                       placement = "right"
                ),
                actionButton(inputId="Indpreviewpriordistributions", label = "Preview Joint Prior Distributions", class = "button")
              ),
              radioGroupButtons(inputId = "KindPriorSpatialEffectInd", label = "Kind of Prior distribution",
                                choices = c("PC-priors" = "PC.prior", "Base Priors" = "Base"), status = "success", justified = TRUE),
              selectInput("optionIndS",
                          label = "Prior options",
                          choices = list("Auto" = "auto", "Custom" = "custom"),
                          selected = "auto"
              ),
              conditionalPanel(
                condition = "input.optionIndS=='custom'&input.KindPriorSpatialEffectInd=='PC.prior'",
                textInput("IndPriorRangePC",
                          label = HTML(paste0(p(HTML(paste0("Range Matérn ", "(\U03C1", tags$sub("0"), ", p", tags$sub("1"), "):"))), p("\\( \\mathbf{P}(\\boldsymbol\\rho<\\boldsymbol\\rho_o)=\\mathbf{p_1} \\)"))),
                          value = "0.5,0.5", placeholder = "diff(limlattice)/2,0.5"
                ),
                textInput("IndPriorStdevPC",
                          label = HTML(paste0(p(HTML(paste0("Stdev. Matérn ", "(\U03C3", tags$sub("0"), ", p", tags$sub("2"), "):"))), p("\\( \\mathbf{P}(\\boldsymbol\\sigma>\\boldsymbol\\sigma_o)=\\mathbf{p_2} \\)"))),
                          value = "0.5,0.5", placeholder = "diff(limlattice)/2,0.5"
                )
              ),
              conditionalPanel(
                condition = "input.optionIndS=='custom'&input.KindPriorSpatialEffectInd=='Base'",
                textInput("IndPriorRangeBase",
                          label = HTML(paste0(p(HTML(paste0("Range Matérn ", "(\U03C1", tags$sub("0"), ", \U03BC", tags$sub("1"), ", \U03C3", tags$sub("1"), "):"))), p("\\(  \\)"))),
                          value = "0.5,0.5,1.0", placeholder = "diff(limlattice)/2,0.5"
                ),
                textInput("IndPriorStdevBase",
                          label = HTML(paste0(p(HTML(paste0("Stdev. Matérn ", "(\U03C3", tags$sub("0"), ", \U03BC", tags$sub("2"), ", \U03C3", tags$sub("2"), "):"))), p("\\(  \\)"))),
                          value = "0.5,0.5,1.0", placeholder = "diff(limlattice)/2,0.5"
                )
              )
            )
          ),
          box(
            id = "advancedINLAInd", width = 12,
            title = "Advanced INLA config.",
            collapsible = TRUE, collapsed = TRUE,
            status = "info", solidHeader = TRUE,
            selectInput("INLAModeInd",
                        label = "INLA Mode",
                        choices = list("Default (Compact)" = "compact",
                                       "Classic" = "classic",
                                       "Experimental" = "experimental")
            ),
            conditionalPanel(
              condition="input.INLAModeInd=='classic'",
              selectInput("strategyapproxINLAInd",
                          label = "Aproximation strategy",
                          choices = list(
                            "Auto" = "auto", "Gaussian" = "gaussian",
                            "Simplified Laplace" = "simplified.laplace",
                            "Laplace" = "laplace", "Adaptive" = "adaptive"
                          ),
                          selected = "auto"
              ),
              selectInput("strategyintINLAInd",
                          label = "Integration strategy",
                          choices = list(
                            "Auto" = "auto", "Central compositive design" = "ccd",
                            "Grid strategy" = "grid",
                            "Empirical Bayes" = "eb"
                          ),
                          selected = "auto"
              )
            ),
            textInput("IndpriorBetacopy",
                      label = "Prior Beta-copy values",
                      value = "0,0.001", placeholder = "0,0.001"
            ),
            radioGroupButtons(
              inputId = "autocustomIndMode",
              label = "Hyperparameters modes",
              choices = list("Default" = "default", "Custom" = "custom"),
              status = "primary"
            ),
            conditionalPanel(
              condition = "input.autocustomIndMode=='custom'",
              textInput("IndModeHyper",
                        label = "Values",
                        value = "0,0", placeholder = "0,0"
              )
            )
          )
        ),
        conditionalPanel(
          condition = "input.IndCustomMesh",
          box(
            title = "Mesh (Ind Model)", width = 9, solidHeader = TRUE, status = "info",
            downloadFileButton(
              id = "ggplotIndMesh",
              downloadtypes = c("png"),
              hovertext = "Download image and data"
            ),
            downloadablePlotUI(
              id = "ggplotIndMesh",
              list("png"),
              btn_valign = "bottom",
              btn_halign = "right",
              btn_overlap = TRUE
            )
          )
        ),
        conditionalPanel(
          condition = "input.Indpreviewpriordistributions>=1",
          box(
            title = "Preview Joint Prior Distributions (Spatial effects)", width = 9, solidHeader = TRUE, collapsible = TRUE, status = "info",
            downloadFileButton(
              id = "IndPreviewJointPriorPlot",
              downloadtypes = c("png"),
              hovertext = "Download image"
            ),
            downloadablePlotUI(
              id = "IndPreviewJointPriorPlot",
              list("png"),
              btn_valign = "bottom",
              btn_halign = "right",
              btn_overlap = TRUE
            )
          )
        ), 
        ## UI_Ind_Outputs_SubSection ====
        conditionalPanel(
          condition = "input.fitInd>=1",
          box(
            title = "Independent model results", width = 9, collapsible = TRUE, solidHeader = TRUE, status = "primary",
            box(
              title = "Response predictive maps", solidHeader = TRUE,
              collapsible = TRUE, width = 12, status = "info",
              column(
                width = 12,
                downloadFileButton(
                  id = "ggplotIndAbundanceMeanMedianStdevFit",
                  downloadtypes = c("png", "csv", "txt"),
                  hovertext = "Download image and data"
                ),
                downloadablePlotUI(
                  id = "ggplotIndAbundanceMeanMedianStdevFit",
                  downloadtypes = c("png", "csv", "txt"),
                  btn_valign = "bottom",
                  btn_halign = "right",
                  btn_overlap = TRUE
                )
              )
            ),
            box(
              title = "Linear predictor maps", solidHeader = TRUE,
              collapsible = TRUE, width = 12, status = "info",
              column(
                width = 12,
                downloadFileButton(
                  id = "ggplotIndPredictorMeanMedianStdevFit",
                  downloadtypes = c("png", "csv", "txt"),
                  hovertext = "Download image and data"
                ),
                downloadablePlotUI(
                  id = "ggplotIndPredictorMeanMedianStdevFit",
                  downloadtypes = c("png", "csv", "txt"),
                  btn_valign = "bottom",
                  btn_halign = "right",
                  btn_overlap = TRUE
                )
              )
            ),
            box(
              title = "Posterior maps of the spatial effect", solidHeader = TRUE,
              collapsible = TRUE, width = 12, status = "info",
              column(
                width = 12,
                downloadFileButton(
                  id = "ggplotIndSpatialMeanMedianStdev",
                  downloadtypes = c("png", "csv", "txt"),
                  hovertext = "Download image and data"
                ),
                downloadablePlotUI(
                  id = "ggplotIndSpatialMeanMedianStdev",
                  downloadtypes = c("png", "csv", "txt"),
                  btn_valign = "bottom",
                  btn_halign = "right",
                  btn_overlap = TRUE
                )
              )
            ),
            box(
              title = "Density functions of the fixed effects", solidHeader = TRUE,
              collapsible = TRUE, width = 12, status = "info",
              column(
                width = 12,
                downloadFileButton(
                  id = "ggplotIndFixParamFit",
                  downloadtypes = c("png", "csv", "txt"),
                  hovertext = "Download image and data"
                ),
                downloadablePlotUI(
                  id = "ggplotIndFixParamFit",
                  downloadtypes = c("png", "csv", "txt"),
                  btn_valign = "bottom",
                  btn_halign = "right",
                  btn_overlap = TRUE
                )
              )
            ),
            box(
              title = "Density functions of the hyperparameters", solidHeader = TRUE,
              collapsible = TRUE, width = 12, status = "info",
              column(
                width = 12,
                downloadFileButton(
                  id = "ggplotIndHyperParamFit",
                  downloadtypes = c("png", "csv", "txt"),
                  hovertext = "Download image and data"
                ),
                downloadablePlotUI(
                  id = "ggplotIndHyperParamFit",
                  downloadtypes = c("png", "csv", "txt"),
                  btn_valign = "bottom",
                  btn_halign = "right",
                  btn_overlap = TRUE
                )
              )
            ),
            box(
              title = "Tables of fixed effects, hyperparameters and internal hyperparameters", solidHeader = TRUE,
              collapsible = TRUE, width = 12, status = "info",
              column(
                width = 6,
                downloadableTableUI("tableIndModelFixedPar", downloadtypes = c("csv", "tsv"))
              ),
              column(
                width = 6,
                downloadableTableUI("tableIndModelHyperPar", downloadtypes = c("csv", "tsv"))
              ),
              column(
                width = 12,
                downloadableTableUI("tableIndModelInternalHyperPar", downloadtypes = c("csv", "tsv"))
              ),
              column(
                width = 12,
                downloadableTableUI("dataIndCPOtable", downloadtypes = c("csv", "tsv"))
              ),
              column(
                width = 6,
                downloadableTableUI("dataIndDICtable", downloadtypes = c("csv", "tsv"))
              ),
              column(
                width = 6,
                downloadableTableUI("dataIndWAICtable", downloadtypes = c("csv", "tsv"))
              )
            )
          )
        )
      )
    ),
    
    # UI_LGCP_Section ----
    tabItem(
      tabName = "modlgcp", withMathJax(),
      fluidRow(
        box(
          width = 3,
          actionButton("fitLgcp", label = "Fit Model", width = "100%", icon = icon("jedi-order")), br(), br(),
          box(id = "SaveLgcpModel", width = 12, title = "Save model data and code", closable = FALSE,
              status = "info", collapsible = TRUE, collapsed = TRUE, solidHeader = TRUE,
              prettyToggle(
                inputId = "saveLgcp",
                label_on = "Save", 
                label_off = "Do not save"
              ),
              textInput(inputId = "LgcpSavePath", label = "Absolute path to save the model", value = "", placeholder = "C:/example/example/"),
              textInput(inputId = "LgcpSaveName", label = "Filename", value = "Lgcp_example.RDS", placeholder = "Lgcp_example.RDS")
          ),
          selectInput("LgcpDataSimulatedLoaded",
                      label = "Select data to be analyzed",
                      choices = list("Simulated" = "sim", "Loaded" = "load"), selected = "sim"
          ),
          selectInput("LgcpRasterSPDE",
                      label = "Use auxilary data",
                      choices = list("Yes" = "raster", "No" = "solvecov"), selected = "solvecov"
          ),
          conditionalPanel(
            condition = "input.LgcpRasterSPDE=='raster'",
            selectInput("LgcpRasterPred",
                        label = "Auxiliary data functionality",
                        selected = "SPDEraster",
                        choices = list("Additional data" = "SPDEraster",
                                       "Fixing prediction points" = "rasterpred")
            )
          ),
          conditionalPanel(
            condition = "input.LgcpRasterPred=='SPDEraster' || input.LgcpRasterSPDE=='solvecov' || input.LgcpDataSimulatedLoaded=='sim'",
            sliderInput("LgcpSPDErasterInput", "Predictive grid dimensionality:",
                        min = 50, max = 200, value = 100
            )
          ),
          sliderInput("dimLgcpmap",
                      label = "Predictive map resolution (spatial effect):",
                      min = 100, max = 300, value = 150
          ),
          box(
            id = "LgcpMeshCustom", width = 12, title = "Custom Mesh", closable = FALSE,
            status = "warning", collapsible = TRUE, collapsed = TRUE, solidHeader = TRUE,
            materialSwitch(
              inputId = "LgcpCustomMesh", label = "Allow custom mesh",
              value = FALSE, status = "primary"
            ),
            conditionalPanel(
              condition = "input.LgcpCustomMesh",
              actionButton("buildLgcpMesh", label = "Build Mesh", class = "button"),
              radioGroupButtons(
                inputId = "boundaryLgcpMesh",
                label = "Outer Boundary Mesh",
                choices = c("Default", "Non-Convex" = "LgcpMeshnonconvex", "Custom" = "LgcpMeshcustomboundary"),
                individual = TRUE,
                checkIcon = list(
                  yes = tags$i(class = "fa fa-circle", style = "color: steelblue"),
                  no = tags$i(class = "fa fa-circle-o", style = "color: steelblue")
                )
              ),
              conditionalPanel(
                condition = "input.boundaryLgcpMesh=='LgcpMeshnonconvex'",
                numericInput(inputId = "curvatureLgcpMesh", label = "Curvature Radious", value = 0),
                numericInput(inputId = "resolutionLgcpMesh", label = "Resolution", value = 40, min = 40)
              ),
              conditionalPanel(
                condition = "input.boundaryLgcpMesh=='LgcpMeshcustomboundary'",
                fileInput(
                  inputId = "shapefileLgcpMesh", label = "Outer Boundary",
                  accept = c(".csv", ".rds"), multiple = TRUE
                )
              ),
              radioGroupButtons(
                inputId = "interiorLgcpMesh",
                label = "Inner Boundary Mesh",
                choices = c("Default", "Non-Convex" = "interiorLgcpMeshnonconvex", "Custom" = "interiorLgcpMeshcustomboundary"),
                individual = TRUE,
                checkIcon = list(
                  yes = tags$i(class = "fa fa-circle", style = "color: steelblue"),
                  no = tags$i(class = "fa fa-circle-o", style = "color: steelblue")
                )
              ),
              conditionalPanel(
                condition = "input.interiorLgcpMesh=='interiorLgcpMeshnonconvex'",
                numericInput(inputId = "interiorcurvatureLgcpMesh", label = "Curvature Radious", value = 0),
                numericInput(inputId = "interiorresolutionLgcpMesh", label = "Resolution", value = 40, min = 40)
              ),
              conditionalPanel(
                condition = "input.interiorLgcpMesh=='interiorLgcpMeshcustomboundary'",
                fileInput(
                  inputId = "interiorshapefileLgcpMesh", label = "Inner Boundary",
                  accept = c(".csv", ".rds"), multiple = TRUE
                )
              ),
              radioGroupButtons(
                inputId = "selectionLgcpMesh",
                label = "Mesh density",
                choices = c("Quantile location" = "qlocation", "Edge length" = "edgelength"),
                individual = TRUE,
                checkIcon = list(
                  yes = tags$i(class = "fa fa-circle", style = "color: steelblue"),
                  no = tags$i(class = "fa fa-circle-o", style = "color: steelblue")
                )
              ),
              conditionalPanel(
                condition = "input.selectionLgcpMesh=='qlocation'",
                numericInput(
                  inputId = "LgcpMeshQloc", label = "Quantile Location",
                  width = "100%", value = 0.03, min = 0, max = 1, step = 0.01
                )
              ),
              conditionalPanel(
                condition = "input.selectionLgcpMesh=='edgelength'",
                numericInput(inputId = "EdgeLengthLgcpMesh", label = "Edge Length", width = "100%", value = 1, min = 0)
              )
            )
          ),
          box(
            id = "LgcpFamily", width = 12, title = "Family and its hyperparameters",
            status = "info", solidHeader = TRUE, collapsible = TRUE,
            selectInput("SelectLgcpFamily",
                        label = "Family distribution",
                        choices = list("Poisson" = "poisson"),
                        selected = "poisson"
            ),
            conditionalPanel(
              condition = "input.SelectLgcpFamily=='beta' || input.SelectLgcpFamily=='gaussian' || input.SelectLgcpFamily=='gamma'",
              radioGroupButtons(
                inputId = "autocustomLgcpFamily",
                label = "Hyperparameters",
                choices = list("Default" = "default", "Custom" = "custom"),
                status = "primary"
              ),
              conditionalPanel(
                condition = "input.autocustomLgcpFamily=='custom'",
                selectInput(inputId="LgcpFamilyPriorKind",
                            label="Prior distribution",
                            choices = list("Log-Gamma"="loggamma", "PC-prior"="pc", "Uniform"="unif", "Flat Uniform"="unifflat"),
                            selected="loggamma"
                ),
                textInput("LgcpFamilyHyper",
                          label = "Prior parameters",
                          value = "0.01, 0.01", placeholder = "0.01, 0.01"
                )
              )
            )
          ),
          box(
            id = "LgcpDefaultComponents", width = 12, title = "Intercept and Explanatory Variables",
            status = "info", solidHeader = TRUE, collapsible = TRUE,
            checkboxGroupInput(
              inputId = "DefaultComponentsLgcp",
              label = "Default Components",
              choices = c("Intercept", "Spatial Effect"),
              selected = c("Intercept", "Spatial Effect")
            ),
            uiOutput("checkBoxLgcpDataFrame"),
            uiOutput("SelectLgcpEffectCov"),
            conditionalPanel(
              condition = "input.LgcpDataSimulatedLoaded=='sim'",
              uiOutput("SelectSimLgcpEffectCov1")
            ),
            conditionalPanel(
              condition = "input.LgcpDataSimulatedLoaded=='load'",
              uiOutput("SelectLoadLgcpEffectCovFactorPred")
            )
          ),
          conditionalPanel(
            condition = "input.DefaultComponentsLgcp.includes('Spatial Effect')",
            box(
              id = "LgcpPriorSpatial", width = 12, title = "Spatial Effect",
              status = "info", solidHeader = TRUE, collapsible = TRUE,
              box(
                title = "Joit prior preview", width = 12, status = "info", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                popify(textInput(inputId = "Lgcprangebasepriorprev", label = "Base Prior Range", value = "0.01,5,200,0.5,0,1", width = NULL, placeholder = "0.01,5,200,0.5,0,1"),
                       title = "Base Prior Range",
                       content = "The <b>meaning of the values</b> according to the order of presentation are: the lower value of the axis, the upper value, the number of breaks, the value for rho0, the mean of the normal and the standard deviation of the normal.",
                       placement = "right"
                ),
                popify(textInput(inputId = "Lgcpsigmabasepriorprev", label = "Base Prior Sigma", value = "0.01,5,200,1,0,1", width = NULL, placeholder = "0.01,5,200,1,0,1"),
                       title = "Base Prior Stdev. Deviation",
                       content = "The <b>meaning of the values</b> according to the order of presentation are: the lower value of the axis, the upper value, the number of breaks, the value for rho0, the mean of the normal and the standard deviation of the normal.",
                       placement = "right"
                ),
                popify(textInput(inputId = "Lgcprangepcpriorprev", label = "PC-Prior Range", value = "0.01,5,200,0.5,0.5", width = NULL, placeholder = "0.01,5,200,0.5,0.5"),
                       title = "PC-Prior Range",
                       content = "The <b>meaning of the values</b> according to the order of presentation are: the lower value of the axis, the upper value, the number of breaks, the value for rho0, the mean of the normal and the standard deviation of the normal.",
                       placement = "right"
                ),
                popify(textInput(inputId = "Lgcpsigmapcpriorprev", label = "PC-Prior Sigma", value = "0.01,5,200,1,0.5", width = NULL, placeholder = "0.01,5,200,1,0.5"),
                       title = "PC-Prior Stdev. Deviation",
                       content = "The <b>meaning of the values</b> according to the order of presentation are: the lower value of the axis, the upper value, the number of breaks, the value for rho0, the mean of the normal and the standard deviation of the normal.",
                       placement = "right"
                ),
                actionButton(inputId="Lgcppreviewpriordistributions", label = "Preview Joint Prior Distributions", class = "button")
              ),
              radioGroupButtons(inputId = "KindPriorSpatialEffectLgcp", label = "Kind of Prior distribution",
                                choices = c("PC-priors" = "PC.prior", "Base Priors" = "Base"), status = "success", justified = TRUE),
              selectInput("optionLgcpS",
                          label = "Prior options",
                          choices = list("Auto" = "auto", "Custom" = "custom"),
                          selected = "auto"
              ),
              conditionalPanel(
                condition = "input.optionLgcpS=='custom'&input.KindPriorSpatialEffectLgcp=='PC.prior'",
                textInput("LgcpPriorRangePC",
                          label = HTML(paste0(p(HTML(paste0("Range Matérn ", "(\U03C1", tags$sub("0"), ", p", tags$sub("1"), "):"))), p("\\( \\mathbf{P}(\\boldsymbol\\rho<\\boldsymbol\\rho_o)=\\mathbf{p_1} \\)"))),
                          value = "0.5,0.5", placeholder = "diff(limlattice)/2,0.5"
                ),
                textInput("LgcpPriorStdevPC",
                          label = HTML(paste0(p(HTML(paste0("Stdev. Matérn ", "(\U03C3", tags$sub("0"), ", p", tags$sub("2"), "):"))), p("\\( \\mathbf{P}(\\boldsymbol\\sigma>\\boldsymbol\\sigma_o)=\\mathbf{p_2} \\)"))),
                          value = "0.5,0.5", placeholder = "diff(limlattice)/2,0.5"
                )
              ),
              conditionalPanel(
                condition = "input.optionLgcpS=='custom'&input.KindPriorSpatialEffectLgcp=='Base'",
                textInput("LgcpPriorRangeBase",
                          label = HTML(paste0(p(HTML(paste0("Range Matérn ", "(\U03C1", tags$sub("0"), ", \U03BC", tags$sub("1"), ", \U03C3", tags$sub("1"), "):"))), p("\\(  \\)"))),
                          value = "0.5,0.5,1.0", placeholder = "diff(limlattice)/2,0.5"
                ),
                textInput("LgcpPriorStdevBase",
                          label = HTML(paste0(p(HTML(paste0("Stdev. Matérn ", "(\U03C3", tags$sub("0"), ", \U03BC", tags$sub("2"), ", \U03C3", tags$sub("2"), "):"))), p("\\(  \\)"))),
                          value = "0.5,0.5,1.0", placeholder = "diff(limlattice)/2,0.5"
                )
              )
            )
          ),
          box(
            id = "advancedINLALgcp", width = 12,
            title = "Advanced INLA config.",
            collapsible = TRUE, collapsed = TRUE,
            status = "info", solidHeader = TRUE,
            selectInput("INLAModeLgcp",
                        label = "INLA Mode",
                        choices = list("Default (Compact)" = "compact",
                                       "Classic" = "classic",
                                       "Experimental" = "experimental")
            ),
            conditionalPanel(
              condition="input.INLAModeLgcp=='classic'",
              selectInput("strategyapproxINLALgcp",
                          label = "Aproximation strategy",
                          choices = list(
                            "Auto" = "auto", "Gaussian" = "gaussian",
                            "Simplified Laplace" = "simplified.laplace",
                            "Laplace" = "laplace", "Adaptive" = "adaptive"
                          ),
                          selected = "auto"
              ),
              selectInput("strategyintINLALgcp",
                          label = "Integration strategy",
                          choices = list(
                            "Auto" = "auto", "Central compositive design" = "ccd",
                            "Grid strategy" = "grid",
                            "Empirical Bayes" = "eb"
                          ),
                          selected = "auto"
              )
            ),
            textInput("LgcppriorBetacopy",
                      label = "Prior Beta-copy values",
                      value = "0,0.001", placeholder = "0,0.001"
            ),
            radioGroupButtons(
              inputId = "autocustomLgcpMode",
              label = "Hyperparameters modes",
              choices = list("Default" = "default", "Custom" = "custom"),
              status = "primary"
            ),
            conditionalPanel(
              condition = "input.autocustomLgcpMode=='custom'",
              textInput("LgcpModeHyper",
                        label = "Values",
                        value = "0,0", placeholder = "0,0"
              )
            )
          )
        ),
        conditionalPanel(
          condition = "input.LgcpCustomMesh",
          box(
            title = "Mesh (Lgcp Model)", width = 9, solidHeader = TRUE, status = "info",
            downloadFileButton(
              id = "ggplotLgcpMesh",
              downloadtypes = c("png"),
              hovertext = "Download image and data"
            ),
            downloadablePlotUI(
              id = "ggplotLgcpMesh",
              list("png"),
              btn_valign = "bottom",
              btn_halign = "right",
              btn_overlap = TRUE
            )
          )
        ),
        conditionalPanel(
          condition = "input.Lgcppreviewpriordistributions>=1",
          box(
            title = "Preview Joint Prior Distributions (Spatial effects)", width = 9, solidHeader = TRUE, collapsible = TRUE, status = "info",
            downloadFileButton(
              id = "LgcpPreviewJointPriorPlot",
              downloadtypes = c("png"),
              hovertext = "Download image"
            ),
            downloadablePlotUI(
              id = "LgcpPreviewJointPriorPlot",
              list("png"),
              btn_valign = "bottom",
              btn_halign = "right",
              btn_overlap = TRUE
            )
          )
        ), 
        ## UI_LGCP_Outputs_SubSection ====
        conditionalPanel(
          condition = "input.fitLgcp>=1",
          box(
            title = "LGCP model results", width = 9, collapsible = TRUE, solidHeader = TRUE, status = "primary",
            box(
              title = "Intensity predictive maps", solidHeader = TRUE,
              collapsible = TRUE, width = 12, status = "info",
              column(
                width = 12,
                downloadFileButton(
                  id = "ggplotLgcpAbundanceMeanMedianStdevFit",
                  downloadtypes = c("png", "csv", "txt"),
                  hovertext = "Download image and data"
                ),
                downloadablePlotUI(
                  id = "ggplotLgcpAbundanceMeanMedianStdevFit",
                  downloadtypes = c("png", "csv", "txt"),
                  btn_valign = "bottom",
                  btn_halign = "right",
                  btn_overlap = TRUE
                )
              )
            ),
            box(
              title = "Log-Intesity predictive maps", solidHeader = TRUE,
              collapsible = TRUE, width = 12, status = "info",
              column(
                width = 12,
                downloadFileButton(
                  id = "ggplotLgcpPredictorMeanMedianStdevFit",
                  downloadtypes = c("png", "csv", "txt"),
                  hovertext = "Download image and data"
                ),
                downloadablePlotUI(
                  id = "ggplotLgcpPredictorMeanMedianStdevFit",
                  downloadtypes = c("png", "csv", "txt"),
                  btn_valign = "bottom",
                  btn_halign = "right",
                  btn_overlap = TRUE
                )
              )
            ),
            box(
              title = "Posterior maps of the spatial effect", solidHeader = TRUE,
              collapsible = TRUE, width = 12, status = "info",
              column(
                width = 12,
                downloadFileButton(
                  id = "ggplotLgcpSpatialMeanMedianStdev",
                  downloadtypes = c("png", "csv", "txt"),
                  hovertext = "Download image and data"
                ),
                downloadablePlotUI(
                  id = "ggplotLgcpSpatialMeanMedianStdev",
                  downloadtypes = c("png", "csv", "txt"),
                  btn_valign = "bottom",
                  btn_halign = "right",
                  btn_overlap = TRUE
                )
              )
            ),
            box(
              title = "Density functions of the fixed effects", solidHeader = TRUE,
              collapsible = TRUE, width = 12, status = "info",
              column(
                width = 12,
                downloadFileButton(
                  id = "ggplotLgcpFixParamFit",
                  downloadtypes = c("png", "csv", "txt"),
                  hovertext = "Download image and data"
                ),
                downloadablePlotUI(
                  id = "ggplotLgcpFixParamFit",
                  downloadtypes = c("png", "csv", "txt"),
                  btn_valign = "bottom",
                  btn_halign = "right",
                  btn_overlap = TRUE
                )
              )
            ),
            box(
              title = "Density functions of the hyperparameters", solidHeader = TRUE,
              collapsible = TRUE, width = 12, status = "info",
              column(
                width = 12,
                downloadFileButton(
                  id = "ggplotLgcpHyperParamFit",
                  downloadtypes = c("png", "csv", "txt"),
                  hovertext = "Download image and data"
                ),
                downloadablePlotUI(
                  id = "ggplotLgcpHyperParamFit",
                  downloadtypes = c("png", "csv", "txt"),
                  btn_valign = "bottom",
                  btn_halign = "right",
                  btn_overlap = TRUE
                )
              )
            ),
            box(
              title = "Tables of fixed effects, hyperparameters and internal hyperparameters", solidHeader = TRUE,
              collapsible = TRUE, width = 12, status = "info",
              column(
                width = 6,
                downloadableTableUI("tableLgcpModelFixedPar", downloadtypes = c("csv", "tsv"))
              ),
              column(
                width = 6,
                downloadableTableUI("tableLgcpModelHyperPar", downloadtypes = c("csv", "tsv"))
              ),
              column(
                width = 12,
                downloadableTableUI("tableLgcpModelInternalHyperPar", downloadtypes = c("csv", "tsv"))
              ),
              column(
                width = 12,
                downloadableTableUI("dataLgcpCPOtable", downloadtypes = c("csv", "tsv"))
              ),
              column(
                width = 6,
                downloadableTableUI("dataLgcpDICtable", downloadtypes = c("csv", "tsv"))
              ),
              column(
                width = 6,
                downloadableTableUI("dataLgcpWAICtable", downloadtypes = c("csv", "tsv"))
              )
            )
          )
        )
      )
    ),
    
    
    # UI_Preferential_Section ----
    tabItem(
      tabName = "modpref", withMathJax(),
      fluidRow(
        box(
          width = 3,
          actionButton("fitPref", label = "Fit Model", width = "100%", icon = icon("jedi-order")), br(), br(),
          box(id = "SavePrefModel", width = 12, title = "Save model data and code", closable = FALSE,
              status = "info", collapsible = TRUE, collapsed = TRUE, solidHeader = TRUE,
              prettyToggle(
                inputId = "savePref",
                label_on = "Save", 
                label_off = "Do not save"
              ),
              textInput(inputId = "PrefSavePath", label = "Absolute path to save the model", value = "", placeholder = "C:/example/example/"),
              textInput(inputId = "PrefSaveName", label = "Filename", value = "Pref_example.RDS", placeholder = "Pref_example.RDS")
              ),
          selectInput("PrefDataSimulatedLoaded",
                      label = "Select data to be analyzed",
                      choices = list("Simulated" = "sim", "Loaded" = "load"), selected = "sim"
          ),
          selectInput("PrefRasterSPDE",
                      label = "Use auxilary data",
                      choices = list("Yes" = "raster", "No" = "solvecov"), selected = "solvecov"
          ),
          conditionalPanel(
            condition = "input.PrefRasterSPDE=='raster'",
            selectInput("PrefRasterPred",
                        label = "Auxiliary data functionality",
                        selected = "SPDEraster",
                        choices = list("Additional data" = "SPDEraster",
                                       "Fixing prediction points" = "rasterpred")
            )
          ),
          conditionalPanel(
            condition = "input.PrefRasterPred=='SPDEraster' || input.PrefRasterSPDE=='solvecov' || input.PrefDataSimulatedLoaded=='sim'",
            sliderInput("PrefSPDErasterInput", "Predictive grid dimensionality:",
                        min = 50, max = 200, value = 100
            )
          ),
          sliderInput("dimPrefmap",
                      label = "Predictive map resolution (spatial effect):",
                      min = 100, max = 300, value = 150
          ),
          box(
            id = "PrefMeshCustom", width = 12, title = "Custom Mesh", closable = FALSE,
            status = "warning", collapsible = TRUE, collapsed = TRUE, solidHeader = TRUE,
            materialSwitch(
              inputId = "PrefCustomMesh", label = "Allow custom mesh",
              value = FALSE, status = "primary"
            ),
            conditionalPanel(
              condition = "input.PrefCustomMesh",
              actionButton("buildPrefMesh", label = "Build Mesh", class = "button"),
              radioGroupButtons(
                inputId = "boundaryPrefMesh",
                label = "Outer Boundary Mesh",
                choices = c("Default", "Non-Convex" = "PrefMeshnonconvex", "Custom" = "PrefMeshcustomboundary"),
                individual = TRUE,
                checkIcon = list(
                  yes = tags$i(class = "fa fa-circle", style = "color: steelblue"),
                  no = tags$i(class = "fa fa-circle-o", style = "color: steelblue")
                )
              ),
              conditionalPanel(
                condition = "input.boundaryPrefMesh=='PrefMeshnonconvex'",
                numericInput(inputId = "curvaturePrefMesh", label = "Curvature Radious", value = 0),
                numericInput(inputId = "resolutionPrefMesh", label = "Resolution", value = 40, min = 40)
              ),
              conditionalPanel(
                condition = "input.boundaryPrefMesh=='PrefMeshcustomboundary'",
                fileInput(
                  inputId = "shapefilePrefMesh", label = "Outer Boundary",
                  accept = c(".csv", ".rds"), multiple = TRUE
                )
              ),
              radioGroupButtons(
                inputId = "interiorPrefMesh",
                label = "Inner Boundary Mesh",
                choices = c("Default", "Non-Convex" = "interiorPrefMeshnonconvex", "Custom" = "interiorPrefMeshcustomboundary"),
                individual = TRUE,
                checkIcon = list(
                  yes = tags$i(class = "fa fa-circle", style = "color: steelblue"),
                  no = tags$i(class = "fa fa-circle-o", style = "color: steelblue")
                )
              ),
              conditionalPanel(
                condition = "input.interiorPrefMesh=='interiorPrefMeshnonconvex'",
                numericInput(inputId = "interiorcurvaturePrefMesh", label = "Curvature Radious", value = 0),
                numericInput(inputId = "interiorresolutionPrefMesh", label = "Resolution", value = 40, min = 40)
              ),
              conditionalPanel(
                condition = "input.interiorPrefMesh=='interiorPrefMeshcustomboundary'",
                fileInput(
                  inputId = "interiorshapefilePrefMesh", label = "Inner Boundary",
                  accept = c(".csv", ".rds"), multiple = TRUE
                )
              ),
              radioGroupButtons(
                inputId = "selectionPrefMesh",
                label = "Mesh density",
                choices = c("Quantile location" = "qlocation", "Edge length" = "edgelength"),
                individual = TRUE,
                checkIcon = list(
                  yes = tags$i(class = "fa fa-circle", style = "color: steelblue"),
                  no = tags$i(class = "fa fa-circle-o", style = "color: steelblue")
                )
              ),
              conditionalPanel(
                condition = "input.selectionPrefMesh=='qlocation'",
                numericInput(
                  inputId = "PrefMeshQloc", label = "Quantile Location",
                  width = "100%", value = 0.03, min = 0, max = 1, step = 0.01
                )
              ),
              conditionalPanel(
                condition = "input.selectionPrefMesh=='edgelength'",
                numericInput(inputId = "EdgeLengthPrefMesh", label = "Edge Length", width = "100%", value = 1, min = 0)
              )
            )
          ),
          box(
            id = "PrefFamily", width = 12, title = "Family and its hyperparameters",
            status = "info", solidHeader = TRUE, collapsible = TRUE,
            selectInput("SelectPrefFamily",
                        label = "Family distribution",
                        choices = list("Bernoulli" = "bernoulli", "Binomial" = "binomial", "Poisson" = "poisson", "Beta" = "beta", "Gaussian" = "gaussian", "Gamma" = "gamma", "Exponential" = "exponential"),
                        selected = "gaussian"
            ),
            conditionalPanel(
              condition = "input.SelectPrefFamily=='beta' || input.SelectPrefFamily=='gaussian' || input.SelectPrefFamily=='gamma'",
              radioGroupButtons(
                inputId = "autocustomPrefFamily",
                label = "Hyperparameters",
                choices = list("Default" = "default", "Custom" = "custom"),
                status = "primary"
              ),
              conditionalPanel(
                condition = "input.autocustomPrefFamily=='custom'",
                selectInput(inputId="PrefFamilyPriorKind",
                            label="Prior distribution",
                            choices = list("Log-Gamma"="loggamma", "PC-prior"="pc", "Uniform"="unif", "Flat Uniform"="unifflat"),
                            selected="loggamma"
                ),
                textInput("PrefFamilyHyper",
                          label = "Prior parameters",
                          value = "0.01, 0.01", placeholder = "0.01, 0.01"
                )
              )
            )
          ),
          box(
            id = "PrefResponseVariable", width = 12, title = "Response variable",
            status = "info", solidHeader = TRUE, collapsible = TRUE,
            uiOutput("UserResponsePref")
          ),
          box(
            id = "PrefDefaultComponents", width = 12, title = "Intercept and Explanatory Variables",
            status = "info", solidHeader = TRUE, collapsible = TRUE,
            checkboxGroupInput(
              inputId = "DefaultComponentsPref",
              label = "Default Components",
              choices = c("Intercept", "Spatial Effect"),
              selected = c("Intercept", "Spatial Effect")
            ),
            uiOutput("checkBoxPrefDataFrame"),
            uiOutput("SelectPrefEffectCov"),
            uiOutput("checkBoxPrefSharing"),
            conditionalPanel(
              condition = "input.PrefDataSimulatedLoaded=='sim'",
              uiOutput("SelectSimPrefEffectCov1")
            ),
            conditionalPanel(
              condition = "input.PrefDataSimulatedLoaded=='load'",
              uiOutput("SelectLoadPrefEffectCovFactorPred")
            )
          ),
          conditionalPanel(
            condition = "input.DefaultComponentsPref.includes('Spatial Effect')",
            box(
              id = "PrefPriorSpatial", width = 12, title = "Spatial Effect",
              status = "info", solidHeader = TRUE, collapsible = TRUE,
              box(
                title = "Joit prior preview", width = 12, status = "info", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                popify(textInput(inputId = "Prefrangebasepriorprev", label = "Base Prior Range", value = "0.01,5,200,0.5,0,1", width = NULL, placeholder = "0.01,5,200,0.5,0,1"),
                       title = "Base Prior Range",
                       content = "The <b>meaning of the values</b> according to the order of presentation are: the lower value of the axis, the upper value, the number of breaks, the value for rho0, the mean of the normal and the standard deviation of the normal.",
                       placement = "right"
                ),
                popify(textInput(inputId = "Prefsigmabasepriorprev", label = "Base Prior Sigma", value = "0.01,5,200,1,0,1", width = NULL, placeholder = "0.01,5,200,1,0,1"),
                       title = "Base Prior Stdev. Deviation",
                       content = "The <b>meaning of the values</b> according to the order of presentation are: the lower value of the axis, the upper value, the number of breaks, the value for rho0, the mean of the normal and the standard deviation of the normal.",
                       placement = "right"
                ),
                popify(textInput(inputId = "Prefrangepcpriorprev", label = "PC-Prior Range", value = "0.01,5,200,0.5,0.5", width = NULL, placeholder = "0.01,5,200,0.5,0.5"),
                       title = "PC-Prior Range",
                       content = "The <b>meaning of the values</b> according to the order of presentation are: the lower value of the axis, the upper value, the number of breaks, the value for rho0, the mean of the normal and the standard deviation of the normal.",
                       placement = "right"
                ),
                popify(textInput(inputId = "Prefsigmapcpriorprev", label = "PC-Prior Sigma", value = "0.01,5,200,1,0.5", width = NULL, placeholder = "0.01,5,200,1,0.5"),
                       title = "PC-Prior Stdev. Deviation",
                       content = "The <b>meaning of the values</b> according to the order of presentation are: the lower value of the axis, the upper value, the number of breaks, the value for rho0, the mean of the normal and the standard deviation of the normal.",
                       placement = "right"
                ),
                actionButton(inputId="Prefpreviewpriordistributions", label = "Preview Joint Prior Distributions", class = "button")
              ),
              radioGroupButtons(inputId = "KindPriorSpatialEffectPref", label = "Kind of Prior distribution",
                                choices = c("PC-priors" = "PC.prior", "Base Priors" = "Base"), status = "success", justified = TRUE),
              selectInput("optionPrefS",
                          label = "Prior options",
                          choices = list("Auto" = "auto", "Custom" = "custom"),
                          selected = "auto"
              ),
              conditionalPanel(
                condition = "input.optionPrefS=='custom'&input.KindPriorSpatialEffectPref=='PC.prior'",
                textInput("PrefPriorRangePC",
                          label = HTML(paste0(p(HTML(paste0("Range Matérn ", "(\U03C1", tags$sub("0"), ", p", tags$sub("1"), "):"))), p("\\( \\mathbf{P}(\\boldsymbol\\rho<\\boldsymbol\\rho_o)=\\mathbf{p_1} \\)"))),
                          value = "0.5,0.5", placeholder = "diff(limlattice)/2,0.5"
                ),
                textInput("PrefPriorStdevPC",
                          label = HTML(paste0(p(HTML(paste0("Stdev. Matérn ", "(\U03C3", tags$sub("0"), ", p", tags$sub("2"), "):"))), p("\\( \\mathbf{P}(\\boldsymbol\\sigma>\\boldsymbol\\sigma_o)=\\mathbf{p_2} \\)"))),
                          value = "0.5,0.5", placeholder = "diff(limlattice)/2,0.5"
                )
              ),
              conditionalPanel(
                condition = "input.optionPrefS=='custom'&input.KindPriorSpatialEffectPref=='Base'",
                textInput("PrefPriorRangeBase",
                          label = HTML(paste0(p(HTML(paste0("Range Matérn ", "(\U03C1", tags$sub("0"), ", \U03BC", tags$sub("1"), ", \U03C3", tags$sub("1"), "):"))), p("\\(  \\)"))),
                          value = "0.5,0.5,1.0", placeholder = "diff(limlattice)/2,0.5"
                ),
                textInput("PrefPriorStdevBase",
                          label = HTML(paste0(p(HTML(paste0("Stdev. Matérn ", "(\U03C3", tags$sub("0"), ", \U03BC", tags$sub("2"), ", \U03C3", tags$sub("2"), "):"))), p("\\(  \\)"))),
                          value = "0.5,0.5,1.0", placeholder = "diff(limlattice)/2,0.5"
                )
              )
            )
          ),
          box(
            id = "advancedINLAPref", width = 12,
            title = "Advanced INLA config.",
            collapsible = TRUE, collapsed = TRUE,
            status = "info", solidHeader = TRUE,
            selectInput("INLAModePref",
                        label = "INLA Mode",
                        choices = list("Default (Compact)" = "compact",
                                       "Classic" = "classic",
                                       "Experimental" = "experimental")
            ),
            conditionalPanel(
              condition="input.INLAModePref=='classic'",
              selectInput("strategyapproxINLAPref",
                          label = "Aproximation strategy",
                          choices = list(
                            "Auto" = "auto", "Gaussian" = "gaussian",
                            "Simplified Laplace" = "simplified.laplace",
                            "Laplace" = "laplace", "Adaptive" = "adaptive"
                          ),
                          selected = "auto"
              ),
              selectInput("strategyintINLAPref",
                          label = "Integration strategy",
                          choices = list(
                            "Auto" = "auto", "Central compositive design" = "ccd",
                            "Grid strategy" = "grid",
                            "Empirical Bayes" = "eb"
                          ),
                          selected = "auto"
              )
            ),
            textInput("PrefpriorBetacopy",
                      label = "Prior Beta-copy values",
                      value = "0,0.001", placeholder = "0,0.001"
            ),
            radioGroupButtons(
              inputId = "autocustomPrefMode",
              label = "Hyperparameters modes",
              choices = list("Default" = "default", "Custom" = "custom"),
              status = "primary"
            ),
            conditionalPanel(
              condition = "input.autocustomPrefMode=='custom'",
              textInput("PrefModeHyper",
                        label = "Values",
                        value = "0,0", placeholder = "0,0"
              )
            )
          )
        ),
        conditionalPanel(
          condition = "input.PrefCustomMesh",
          box(
            title = "Mesh (Pref Model)", width = 9, solidHeader = TRUE, status = "info",
            downloadFileButton(
              id = "ggplotPrefMesh",
              downloadtypes = c("png"),
              hovertext = "Download image and data"
            ),
            downloadablePlotUI(
              id = "ggplotPrefMesh",
              list("png"),
              btn_valign = "bottom",
              btn_halign = "right",
              btn_overlap = TRUE
            )
          )
        ),
        conditionalPanel(
          condition = "input.Prefpreviewpriordistributions>=1",
          box(
            title = "Preview Joint Prior Distributions (Spatial effects)", width = 9, solidHeader = TRUE, collapsible = TRUE, status = "info",
            downloadFileButton(
              id = "PrefPreviewJointPriorPlot",
              downloadtypes = c("png"),
              hovertext = "Download image"
            ),
            downloadablePlotUI(
              id = "PrefPreviewJointPriorPlot",
              list("png"),
              btn_valign = "bottom",
              btn_halign = "right",
              btn_overlap = TRUE
            )
          )
        ), 
        ## UI_Pref_Outputs_SubSection ====
        conditionalPanel(
          condition = "input.fitPref>=1",
          box(
            title = "Preferential model results", width = 9, collapsible = TRUE, solidHeader = TRUE, status = "primary",
            box(
              title = "Response predictive maps", solidHeader = TRUE,
              collapsible = TRUE, width = 12, status = "info",
              column(
                width = 12,
                downloadFileButton(
                  id = "ggplotPrefAbundanceMeanMedianStdevFit",
                  downloadtypes = c("png", "csv", "txt"),
                  hovertext = "Download image and data"
                ),
                downloadablePlotUI(
                  id = "ggplotPrefAbundanceMeanMedianStdevFit",
                  downloadtypes = c("png", "csv", "txt"),
                  btn_valign = "bottom",
                  btn_halign = "right",
                  btn_overlap = TRUE
                )
              )
            ),
            box(
              title = "Linear predictor maps", solidHeader = TRUE,
              collapsible = TRUE, width = 12, status = "info",
              column(
                width = 12,
                downloadFileButton(
                  id = "ggplotPrefPredictorMeanMedianStdevFit",
                  downloadtypes = c("png", "csv", "txt"),
                  hovertext = "Download image and data"
                ),
                downloadablePlotUI(
                  id = "ggplotPrefPredictorMeanMedianStdevFit",
                  downloadtypes = c("png", "csv", "txt"),
                  btn_valign = "bottom",
                  btn_halign = "right",
                  btn_overlap = TRUE
                )
              )
            ),
            box(
              title = "Posterior maps of the spatial effect", solidHeader = TRUE,
              collapsible = TRUE, width = 12, status = "info",
              column(
                width = 12,
                downloadFileButton(
                  id = "ggplotPrefSpatialMeanMedianStdev",
                  downloadtypes = c("png", "csv", "txt"),
                  hovertext = "Download image and data"
                ),
                downloadablePlotUI(
                  id = "ggplotPrefSpatialMeanMedianStdev",
                  downloadtypes = c("png", "csv", "txt"),
                  btn_valign = "bottom",
                  btn_halign = "right",
                  btn_overlap = TRUE
                )
              )
            ),
            box(
              title = "Density functions of the fixed effects", solidHeader = TRUE,
              collapsible = TRUE, width = 12, status = "info",
              column(
                width = 12,
                downloadFileButton(
                  id = "ggplotPrefFixParamFit",
                  downloadtypes = c("png", "csv", "txt"),
                  hovertext = "Download image and data"
                ),
                downloadablePlotUI(
                  id = "ggplotPrefFixParamFit",
                  downloadtypes = c("png", "csv", "txt"),
                  btn_valign = "bottom",
                  btn_halign = "right",
                  btn_overlap = TRUE
                )
              )
            ),
            box(
              title = "Density functions of the hyperparameters", solidHeader = TRUE,
              collapsible = TRUE, width = 12, status = "info",
              column(
                width = 12,
                downloadFileButton(
                  id = "ggplotPrefHyperParamFit",
                  downloadtypes = c("png", "csv", "txt"),
                  hovertext = "Download image and data"
                ),
                downloadablePlotUI(
                  id = "ggplotPrefHyperParamFit",
                  downloadtypes = c("png", "csv", "txt"),
                  btn_valign = "bottom",
                  btn_halign = "right",
                  btn_overlap = TRUE
                )
              )
            ),
            box(
              title = "Tables of fixed effects, hyperparameters and internal hyperparameters", solidHeader = TRUE,
              collapsible = TRUE, width = 12, status = "info",
              column(
                width = 6,
                downloadableTableUI("tablePrefModelFixedPar", downloadtypes = c("csv", "tsv"))
              ),
              column(
                width = 6,
                downloadableTableUI("tablePrefModelHyperPar", downloadtypes = c("csv", "tsv"))
              ),
              column(
                width = 12,
                downloadableTableUI("tablePrefModelInternalHyperPar", downloadtypes = c("csv", "tsv"))
              ),
              column(
                width = 12,
                downloadableTableUI("dataPrefCPOtable", downloadtypes = c("csv", "tsv"))
              ),
              column(
                width = 6,
                downloadableTableUI("dataPrefDICtable", downloadtypes = c("csv", "tsv"))
              ),
              column(
                width = 6,
                downloadableTableUI("dataPrefWAICtable", downloadtypes = c("csv", "tsv"))
              )
            )
          )
        )
      )
    ),
    
    # UI_Mixture_Section ----
    tabItem(
      tabName = "modmixture", withMathJax(),
      fluidRow(
        box(
          width = 3,
          actionButton("fitMixture", label = "Fit Model", width = "100%", icon = icon("jedi-order")), br(), br(),
          box(id = "SaveMixtureModel", width = 12, title = "Save model data and code", closable = FALSE,
              status = "info", collapsible = TRUE, collapsed = TRUE, solidHeader = TRUE,
              prettyToggle(
                inputId = "saveMixture",
                label_on = "Save", 
                label_off = "Do not save"
              ),
              textInput(inputId = "MixtureSavePath", label = "Absolute path to save the model", value = "", placeholder = "C:/example/example/"),
              textInput(inputId = "MixtureSaveName", label = "Filename", value = "Mixture_example.RDS", placeholder = "Mixture_example.RDS")
          ),
          selectInput("MixtureDataSimulatedLoaded",
                      label = "Select data to be analyzed",
                      choices = list("Simulated" = "sim", "Loaded" = "load"), selected = "sim"
          ),
          selectInput("MixtureRasterSPDE",
                      label = "Use auxilary data",
                      choices = list("Yes" = "raster", "No" = "solvecov"), selected = "solvecov"
          ),
          conditionalPanel(
            condition = "input.MixtureRasterSPDE=='raster'",
            selectInput("MixtureRasterPred",
                        label = "Auxiliary data functionality",
                        selected = "SPDEraster",
                        choices = list("Additional data" = "SPDEraster",
                                       "Fixing prediction points" = "rasterpred")
            )
          ),
          conditionalPanel(
            condition = "input.MixtureRasterPred=='SPDEraster' || input.MixtureRasterSPDE=='solvecov' || input.MixtureDataSimulatedLoaded=='sim'",
            sliderInput("MixtureSPDErasterInput", "Predictive grid dimensionality:",
                        min = 50, max = 200, value = 100
            )
          ),
          sliderInput("dimMixturemap",
                      label = "Predictive map resolution (spatial effect):",
                      min = 100, max = 300, value = 150
          ),
          box(
            id = "MixtureMeshCustom", width = 12, title = "Custom Mesh", closable = FALSE,
            status = "warning", collapsible = TRUE, collapsed = TRUE, solidHeader = TRUE,
            materialSwitch(
              inputId = "MixtureCustomMesh", label = "Allow custom mesh",
              value = FALSE, status = "primary"
            ),
            conditionalPanel(
              condition = "input.MixtureCustomMesh",
              actionButton("buildMixtureMesh", label = "Build Mesh", class = "button"),
              radioGroupButtons(
                inputId = "boundaryMixtureMesh",
                label = "Outer Boundary Mesh",
                choices = c("Default", "Non-Convex" = "MixtureMeshnonconvex", "Custom" = "MixtureMeshcustomboundary"),
                individual = TRUE,
                checkIcon = list(
                  yes = tags$i(class = "fa fa-circle", style = "color: steelblue"),
                  no = tags$i(class = "fa fa-circle-o", style = "color: steelblue")
                )
              ),
              conditionalPanel(
                condition = "input.boundaryMixtureMesh=='MixtureMeshnonconvex'",
                numericInput(inputId = "curvatureMixtureMesh", label = "Curvature Radious", value = 0),
                numericInput(inputId = "resolutionMixtureMesh", label = "Resolution", value = 40, min = 40)
              ),
              conditionalPanel(
                condition = "input.boundaryMixtureMesh=='MixtureMeshcustomboundary'",
                fileInput(
                  inputId = "shapefileMixtureMesh", label = "Outer Boundary",
                  accept = c(".csv", ".rds"), multiple = TRUE
                )
              ),
              radioGroupButtons(
                inputId = "interiorMixtureMesh",
                label = "Inner Boundary Mesh",
                choices = c("Default", "Non-Convex" = "interiorMixtureMeshnonconvex", "Custom" = "interiorMixtureMeshcustomboundary"),
                individual = TRUE,
                checkIcon = list(
                  yes = tags$i(class = "fa fa-circle", style = "color: steelblue"),
                  no = tags$i(class = "fa fa-circle-o", style = "color: steelblue")
                )
              ),
              conditionalPanel(
                condition = "input.interiorMixtureMesh=='interiorMixtureMeshnonconvex'",
                numericInput(inputId = "interiorcurvatureMixtureMesh", label = "Curvature Radious", value = 0),
                numericInput(inputId = "interiorresolutionMixtureMesh", label = "Resolution", value = 40, min = 40)
              ),
              conditionalPanel(
                condition = "input.interiorMixtureMesh=='interiorMixtureMeshcustomboundary'",
                fileInput(
                  inputId = "interiorshapefileMixtureMesh", label = "Inner Boundary",
                  accept = c(".csv", ".rds"), multiple = TRUE
                )
              ),
              radioGroupButtons(
                inputId = "selectionMixtureMesh",
                label = "Mesh density",
                choices = c("Quantile location" = "qlocation", "Edge length" = "edgelength"),
                individual = TRUE,
                checkIcon = list(
                  yes = tags$i(class = "fa fa-circle", style = "color: steelblue"),
                  no = tags$i(class = "fa fa-circle-o", style = "color: steelblue")
                )
              ),
              conditionalPanel(
                condition = "input.selectionMixtureMesh=='qlocation'",
                numericInput(
                  inputId = "MixtureMeshQloc", label = "Quantile Location",
                  width = "100%", value = 0.03, min = 0, max = 1, step = 0.01
                )
              ),
              conditionalPanel(
                condition = "input.selectionMixtureMesh=='edgelength'",
                numericInput(inputId = "EdgeLengthMixtureMesh", label = "Edge Length", width = "100%", value = 1, min = 0)
              )
            )
          ),
          box(
            id = "MixtureFamily", width = 12, title = "Family and its hyperparameters",
            status = "info", solidHeader = TRUE, collapsible = TRUE,
            selectInput("SelectMixtureFamily",
                        label = "Family distribution",
                        choices = list("Bernoulli" = "bernoulli", "Binomial" = "binomial", "Poisson" = "poisson", "Beta" = "beta", "Gaussian" = "gaussian", "Gamma" = "gamma", "Exponential" = "exponential"),
                        selected = "gaussian"
            ),
            conditionalPanel(
              condition = "input.SelectMixtureFamily=='beta' || input.SelectMixtureFamily=='gaussian' || input.SelectMixtureFamily=='gamma'",
              radioGroupButtons(
                inputId = "autocustomMixtureFamily",
                label = "Hyperparameters",
                choices = list("Default" = "default", "Custom" = "custom"),
                status = "primary"
              ),
              conditionalPanel(
                condition = "input.autocustomMixtureFamily=='custom'",
                selectInput(inputId="MixtureFamilyPriorKind",
                            label="Prior distribution",
                            choices = list("Log-Gamma"="loggamma", "PC-prior"="pc", "Uniform"="unif", "Flat Uniform"="unifflat"),
                            selected="loggamma"
                ),
                textInput("MixtureFamilyHyper",
                          label = "Prior parameters",
                          value = "0.01, 0.01", placeholder = "0.01, 0.01"
                )
              )
            )
          ),
          box(
            id = "MixtureResponseVariable", width = 12, title = "Response variable",
            status = "info", solidHeader = TRUE, collapsible = TRUE,
            uiOutput("UserResponseMixture")
          ),
          box(
            id = "MixtureDefaultComponents", width = 12, title = "Intercept and Explanatory Variables",
            status = "info", solidHeader = TRUE, collapsible = TRUE,
            checkboxGroupInput(
              inputId = "DefaultComponentsMixture",
              label = "Default Components",
              choices = c("Intercept", "Spatial Effect"),
              selected = c("Intercept", "Spatial Effect")
            ),
            uiOutput("checkBoxMixtureDataFrame"),
            uiOutput("checkBoxSelectMixtureDependents"),
            uiOutput("SelectMixtureEffectCov"),
            uiOutput("checkBoxMixtureSharing"),
            conditionalPanel(
              condition = "input.MixtureDataSimulatedLoaded=='sim'",
              uiOutput("SelectSimMixtureEffectCov1")
            ),
            conditionalPanel(
              condition = "input.MixtureDataSimulatedLoaded=='load'",
              uiOutput("SelectLoadMixtureEffectCovFactorPred")
            )
          ),
          conditionalPanel(
            condition = "input.DefaultComponentsMixture.includes('Spatial Effect')",
            box(
              id = "MixturePriorSpatial", width = 12, title = "Spatial Effect",
              status = "info", solidHeader = TRUE, collapsible = TRUE,
              box(
                title = "Joit prior preview", width = 12, status = "info", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                popify(textInput(inputId = "Mixturerangebasepriorprev", label = "Base Prior Range", value = "0.01,5,200,0.5,0,1", width = NULL, placeholder = "0.01,5,200,0.5,0,1"),
                       title = "Base Prior Range",
                       content = "The <b>meaning of the values</b> according to the order of presentation are: the lower value of the axis, the upper value, the number of breaks, the value for rho0, the mean of the normal and the standard deviation of the normal.",
                       placement = "right"
                ),
                popify(textInput(inputId = "Mixturesigmabasepriorprev", label = "Base Prior Sigma", value = "0.01,5,200,1,0,1", width = NULL, placeholder = "0.01,5,200,1,0,1"),
                       title = "Base Prior Stdev. Deviation",
                       content = "The <b>meaning of the values</b> according to the order of presentation are: the lower value of the axis, the upper value, the number of breaks, the value for rho0, the mean of the normal and the standard deviation of the normal.",
                       placement = "right"
                ),
                popify(textInput(inputId = "Mixturerangepcpriorprev", label = "PC-Prior Range", value = "0.01,5,200,0.5,0.5", width = NULL, placeholder = "0.01,5,200,0.5,0.5"),
                       title = "PC-Prior Range",
                       content = "The <b>meaning of the values</b> according to the order of presentation are: the lower value of the axis, the upper value, the number of breaks, the value for rho0, the mean of the normal and the standard deviation of the normal.",
                       placement = "right"
                ),
                popify(textInput(inputId = "Mixturesigmapcpriorprev", label = "PC-Prior Sigma", value = "0.01,5,200,1,0.5", width = NULL, placeholder = "0.01,5,200,1,0.5"),
                       title = "PC-Prior Stdev. Deviation",
                       content = "The <b>meaning of the values</b> according to the order of presentation are: the lower value of the axis, the upper value, the number of breaks, the value for rho0, the mean of the normal and the standard deviation of the normal.",
                       placement = "right"
                ),
                actionButton(inputId="Mixturepreviewpriordistributions", label = "Preview Joint Prior Distributions", class = "button")
              ),
              radioGroupButtons(inputId = "KindPriorSpatialEffectMixture", label = "Kind of Prior distribution", 
                                choices = c("PC-priors" = "PC.prior", "Base Priors" = "Base"), status = "success", justified = TRUE),
              selectInput("optionMixtureS",
                          label = "Prior options",
                          choices = list("Auto" = "auto", "Custom" = "custom"),
                          selected = "auto"
              ),
              conditionalPanel(
                condition = "input.optionMixtureS=='custom'&input.KindPriorSpatialEffectMixture=='PC.prior'",
                textInput("MixturePriorRangePC",
                          label = HTML(paste0(p(HTML(paste0("Range Matérn ", "(\U03C1", tags$sub("0"), ", p", tags$sub("1"), "):"))), p("\\( \\mathbf{P}(\\boldsymbol\\rho<\\boldsymbol\\rho_o)=\\mathbf{p_1} \\)"))),
                          value = "0.5,0.5", placeholder = "diff(limlattice)/2,0.5"
                ),
                textInput("MixturePriorStdevPC",
                          label = HTML(paste0(p(HTML(paste0("Stdev. Matérn ", "(\U03C3", tags$sub("0"), ", p", tags$sub("2"), "):"))), p("\\( \\mathbf{P}(\\boldsymbol\\sigma>\\boldsymbol\\sigma_o)=\\mathbf{p_2} \\)"))),
                          value = "0.5,0.5", placeholder = "diff(limlattice)/2,0.5"
                )
              ),
              conditionalPanel(
                condition = "input.optionMixtureS=='custom'&input.KindPriorSpatialEffectMixture=='Base'",
                textInput("MixturePriorRangeBase",
                          label = HTML(paste0(p(HTML(paste0("Range Matérn ", "(\U03C1", tags$sub("0"), ", \U03BC", tags$sub("1"), ", \U03C3", tags$sub("1"), "):"))), p("\\(  \\)"))),
                          value = "0.5,0.5,1.0", placeholder = "diff(limlattice)/2,0.5"
                ),
                textInput("MixturePriorStdevBase",
                          label = HTML(paste0(p(HTML(paste0("Stdev. Matérn ", "(\U03C3", tags$sub("0"), ", \U03BC", tags$sub("2"), ", \U03C3", tags$sub("2"), "):"))), p("\\(  \\)"))),
                          value = "0.5,0.5,1.0", placeholder = "diff(limlattice)/2,0.5"
                )
              )
            )
          ),
          box(
            id = "advancedINLAMixture", width = 12,
            title = "Advanced INLA config.",
            collapsible = TRUE, collapsed = TRUE,
            status = "info", solidHeader = TRUE,
            selectInput("INLAModeMixture",
                        label = "INLA Mode",
                        choices = list("Default (Compact)" = "compact",
                                       "Classic" = "classic",
                                       "Experimental" = "experimental")
            ),
            conditionalPanel(
              condition="input.INLAModeMixture=='classic'",
              selectInput("strategyapproxINLAMixture",
                          label = "Aproximation strategy",
                          choices = list(
                            "Auto" = "auto", "Gaussian" = "gaussian",
                            "Simplified Laplace" = "simplified.laplace",
                            "Laplace" = "laplace", "Adaptive" = "adaptive"
                          ),
                          selected = "auto"
              ),
              selectInput("strategyintINLAMixture",
                          label = "Integration strategy",
                          choices = list(
                            "Auto" = "auto", "Central compositive design" = "ccd",
                            "Grid strategy" = "grid",
                            "Empirical Bayes" = "eb"
                          ),
                          selected = "auto"
              )
            ),
            textInput("MixturepriorBetacopy",
                      label = "Prior Beta-copy values",
                      value = "0,0.001", placeholder = "0,0.001"
            ),
            radioGroupButtons(
              inputId = "autocustomMixtureMode",
              label = "Hyperparameters modes",
              choices = list("Default" = "default", "Custom" = "custom"),
              status = "primary"
            ),
            conditionalPanel(
              condition = "input.autocustomMixtureMode=='custom'",
              textInput("MixtureModeHyper",
                        label = "Values",
                        value = "0,0", placeholder = "0,0"
              )
            )
          )
        ),
        conditionalPanel(
          condition = "input.MixtureCustomMesh",
          box(
            title = "Mesh (Mixture Model)", width = 9, solidHeader = TRUE, status = "info",
            downloadFileButton(
              id = "ggplotMixtureMesh",
              downloadtypes = c("png"),
              hovertext = "Download image and data"
            ),
            downloadablePlotUI(
              id = "ggplotMixtureMesh",
              list("png"),
              btn_valign = "bottom",
              btn_halign = "right",
              btn_overlap = TRUE
            )
          )
        ),
        conditionalPanel(
          condition = "input.Mixturepreviewpriordistributions>=1",
          box(
            title = "Preview Joint Prior Distributions (Spatial effects)", width = 9, solidHeader = TRUE, collapsible = TRUE, status = "info",
            downloadFileButton(
              id = "MixturePreviewJointPriorPlot",
              downloadtypes = c("png"),
              hovertext = "Download image"
            ),
            downloadablePlotUI(
              id = "MixturePreviewJointPriorPlot",
              list("png"),
              btn_valign = "bottom",
              btn_halign = "right",
              btn_overlap = TRUE
            )
          )
        ), ## UI_Mixture_Outputs_SubSection ====
        conditionalPanel(
          condition = "input.fitMixture>=1",
          box(
            title = "Mixture model results", width = 9, collapsible = TRUE, solidHeader = TRUE, status = "primary",
            box(
              title = "Response predictive maps", solidHeader = TRUE,
              collapsible = TRUE, width = 12, status = "info",
              column(
                width = 12,
                downloadFileButton(
                  id = "ggplotMixtureAbundanceMeanMedianStdevFit",
                  downloadtypes = c("png", "csv", "txt"),
                  hovertext = "Download image and data"
                ),
                downloadablePlotUI(
                  id = "ggplotMixtureAbundanceMeanMedianStdevFit",
                  downloadtypes = c("png", "csv", "txt"),
                  btn_valign = "bottom",
                  btn_halign = "right",
                  btn_overlap = TRUE
                )
              )
            ),
            box(
              title = "Linear predictor maps", solidHeader = TRUE,
              collapsible = TRUE, width = 12, status = "info",
              column(
                width = 12,
                downloadFileButton(
                  id = "ggplotMixturePredictorMeanMedianStdevFit",
                  downloadtypes = c("png", "csv", "txt"),
                  hovertext = "Download image and data"
                ),
                downloadablePlotUI(
                  id = "ggplotMixturePredictorMeanMedianStdevFit",
                  downloadtypes = c("png", "csv", "txt"),
                  btn_valign = "bottom",
                  btn_halign = "right",
                  btn_overlap = TRUE
                )
              )
            ),
            box(
              title = "Posterior maps of the spatial effect", solidHeader = TRUE,
              collapsible = TRUE, width = 12, status = "info",
              column(
                width = 12,
                downloadFileButton(
                  id = "ggplotMixtureSpatialMeanMedianStdev",
                  downloadtypes = c("png", "csv", "txt"),
                  hovertext = "Download image and data"
                ),
                downloadablePlotUI(
                  id = "ggplotMixtureSpatialMeanMedianStdev",
                  downloadtypes = c("png", "csv", "txt"),
                  btn_valign = "bottom",
                  btn_halign = "right",
                  btn_overlap = TRUE
                )
              )
            ),
            box(
              title = "Density functions of the fixed effects", solidHeader = TRUE,
              collapsible = TRUE, width = 12, status = "info",
              column(
                width = 12,
                downloadFileButton(
                  id = "ggplotMixtureFixParamFit",
                  downloadtypes = c("png", "csv", "txt"),
                  hovertext = "Download image and data"
                ),
                downloadablePlotUI(
                  id = "ggplotMixtureFixParamFit",
                  downloadtypes = c("png", "csv", "txt"),
                  btn_valign = "bottom",
                  btn_halign = "right",
                  btn_overlap = TRUE
                )
              )
            ),
            box(
              title = "Density functions of the hyperparameters", solidHeader = TRUE,
              collapsible = TRUE, width = 12, status = "info",
              column(
                width = 12,
                downloadFileButton(
                  id = "ggplotMixtureHyperParamFit",
                  downloadtypes = c("png", "csv", "txt"),
                  hovertext = "Download image and data"
                ),
                downloadablePlotUI(
                  id = "ggplotMixtureHyperParamFit",
                  downloadtypes = c("png", "csv", "txt"),
                  btn_valign = "bottom",
                  btn_halign = "right",
                  btn_overlap = TRUE
                )
              )
            ),
            box(
              title = "Tables of fixed effects, hyperparameters and internal hyperparameters", solidHeader = TRUE,
              collapsible = TRUE, width = 12, status = "info",
              column(
                width = 6,
                downloadableTableUI("tableMixtureModelFixedPar", downloadtypes = c("csv", "tsv"))
              ),
              column(
                width = 6,
                downloadableTableUI("tableMixtureModelHyperPar", downloadtypes = c("csv", "tsv"))
              ),
              column(
                width = 12,
                downloadableTableUI("tableMixtureModelInternalHyperPar", downloadtypes = c("csv", "tsv"))
              ),
              column(
                width = 12,
                downloadableTableUI("dataMixtureCPOtable", downloadtypes = c("csv", "tsv"))
              ),
              column(
                width = 6,
                downloadableTableUI("dataMixtureDICtable", downloadtypes = c("csv", "tsv"))
              ),
              column(
                width = 6,
                downloadableTableUI("dataMixtureWAICtable", downloadtypes = c("csv", "tsv"))
              )
            )
          )
        )
      )
    )
    #  -----
  )
)
controlbar <- dashboardControlbar(
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
  br(), br(), br(), br(),
  column(
    width = 12,
    actionButton("EnhancementIrrGrid", label = "Enhance irregular grid images", class = "button")
  ),
  br(), br(), br(),
  column(
    width = 12,
    materialSwitch(inputId = "RefineInterpolationIrrGrid", label = "Alternative enhancement", value = FALSE, status = "primary")
  )
)
dashboardPage(skin = "red", header, sidebar, body, controlbar)
