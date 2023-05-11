library(shinydashboard)
library(plotly)

scan.num <- 20

shinyUI(
  pageWithSidebar(
    headerPanel("Spectra processing parameters investigation"),
    sidebarPanel(
      width = 3,
      h3("Input spectra file"),
      shiny::selectInput("input.select.file",
                         "Select spectra file",
                         "Loading files from DB",
                         multiple = FALSE,
                         selectize = FALSE),
      h3("Processing parameters"),
      shiny::selectInput("input.select.mode",
                         "Ion acquisition mode",
                         "Select spectra file",
                         multiple = FALSE,
                         selectize = FALSE),
      numericInput("input.scan.num", "MS scan number", value = scan.num),
      numericInput("input.scan.num.l", "MS scan position to the left", value = 10),
      numericInput("input.scan.num.r", "MS scan position to the right", value = 10),
      numericInput("input.halfWindowSize", "halfWindowSize", value = 9),
      numericInput('input.tol.align', "Mass peak tolerance align", value = 2e-4),
      numericInput('input.tol.bin.peaks', "Mass peak tolerance binPeaks", value = 2e-4)
      ),
    mainPanel(
      h3("Resulted spectra"),
      plotlyOutput("pid.raw", height = "2100px")
    )
  )
)
