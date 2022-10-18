# Shiny App for Spatial Modeling (RStudio & R-INLA)

This application is designed for the resolution of spatial models using the Bayesian paradigm and the INLA methodology. Therefore, we will briefly explain its functionality and its main sections. In addition, in this *README* file will be expose the theoretical foundaments behind  the application. Which means, a summary of Bayesian inference and the foundaments of the INLA methodology.

To run this app you can use the function `runGitHub("ShinyAppSpatialModelFeedback", "MarioFigueiraP")` from `library(shiny)`.

# App dependencies

As mentioned above, this application was built to solve spatial models, for wich it needs several packages, although these dependencies would be installed automatically when running the application itself. If the automatic process fails you can try to install these packages manually, which are:

```
install.packages("shiny")
install.packages("shinydashboard")
install.packages("shinyWidgets")
install.packages("shinydashboardPlus")
install.packages("shinyBS")
install.packages("shinyjs")
install.packages("periscope")
install.packages("splines")
install.packages("INLA",
  repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), 
  dep=TRUE) # It is a core package for the app.
install.packages("inlabru")
install.packages("ggplot2")
install.packages("lattice")
install.packages("rintrojs")
install.packages("patchwork")
install.packages("viridis")
install.packages("rgeos")
install.packages("dplyr")
```

The INLA package installation could give some problems or errors, in such case it is desirable to visit the [INLA home page](https://www.r-inla.org/), where the installation is explained in some detail and many FAQ are answered.

# Main application sections

The app is made up of four major blocks:

1. Introduction.
2. Data Simulation.
3. Upload Data.
4. Model Analysis.

*Data Simulation* and *Upload Data* are sequentially linked to *Model Anaysis*
