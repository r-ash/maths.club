f <- function(x) {
  exp(-(x^2)) * 0.5 + exp(-((x - 1)^2)) * 0.3 + exp(-((x + 2)^2)) * 0.2
}

ui <- bslib::page_sidebar(
  theme = bslib::bs_theme(bootswatch = "minty"),
  sidebar = bslib::sidebar(
    shiny::numericInput("initial", "Start point", 0),
    shiny::numericInput("nsamples", "Number of samples", 10000),
    shiny::numericInput("sd", "Proposal dist sd", 0.3, step = 0.1)
  ),
  plotOutput("samples")
)

server <- function(input, output, session) {
  sample_data <- shiny::reactive({
    maths.club:::metropolis(f, input$initial, input$nsamples, input$sd)
  })
  x_min <- shiny::reactive({
    if (input$initial < -5) {
      input$initial
    } else {
      -5
    }
  })
  x_max <- shiny::reactive({
    if (input$initial > 5) {
      input$initial
    } else {
      5
    }
  })

  output$samples <- shiny::renderPlot({
    maths.club:::plot_both(sample_data(), f, x_min(), x_max())
  })
}

shiny::shinyApp(ui, server)
