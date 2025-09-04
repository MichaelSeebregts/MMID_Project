
library(shiny)
library(tidyverse)
library(glue)
library(bslib)
library(shinycssloaders) 

source(file.path("SBRMIC008_MMID_App", "MainModel.R"))

# Helper: pull the deSolve output matrix out of whatever run_model() returns
get_out_matrix <- function(res) {
  if (is.null(res)) return(NULL)
  # Case 1: run_model() returns the matrix directly
  if (is.matrix(res)) return(res)
  # Case 2: run_model() returns a list — take the first matrix-looking element
  if (is.list(res)) {
    i <- which(vapply(res, function(x) is.matrix(x) && ncol(x) >= 2, logical(1)))
    if (length(i)) return(res[[i[1]]])
  }
  stop("Unexpected shape from run_model(); expected a matrix or a list containing one.")
}

# Series extractor that uses the model's time grid (first column of the out matrix)
get_series <- function(out, vname, patch) {
  j <- match(vname, var_names)
  if (is.na(j)) stop("Unknown variable name: ", vname)
  out[, 1 + varind[j, patch], drop = TRUE]
}


ui <- navbarPage(
  id = "topnav",
  title = "Malaria Model",
  theme = bs_theme(version = 5, bootswatch = "flatly",
                   base_font = bslib::font_google("Inter")),
  
  # 1) HOME (landing) ---------------------------------------------------------
  tabPanel(title = "Home", value = "home",
           bslib::layout_columns(col_widths = c(7,5), gap = "24px",
                                 bslib::card(
                                   bslib::card_header("Welcome"),
                                   tags$p("Use this app to explore the malaria model across patches."),
                                   tags$ul(
                                     tags$li("Adjust λᵥ (vivax) and λ_F (falciparum)"),
                                     tags$li("Run the model and visualize trajectories"),
                                     tags$li("Compare Vivax components vs Falciparum")
                                   ),
                                   actionButton("go_to_model", "Start exploring", class = "btn-primary")
                                 ),
                                 bslib::card(
                                   bslib::card_header("Tips"),
                                   tags$p("Keep 'Parameters.xlsx' in the app directory."),
                                   tags$p("Heavy computations (transition setup) occur once at startup.")
                                 )
           )
  ),
  
  # 2) MODEL tab (your existing app) -----------------------------------------
  tabPanel(title = "Model", value = "model",
           # No sidebar here: just a tabset with cards
           tabsetPanel(id = "plot_tabs",
                       tabPanel("Overview",
                                bslib::card(
                                  bslib::card_header("Overview trajectories"),
                                  plotOutput("trajPlot", height = "420px"),
                                  verbatimTextOutput("summaryText")
                                )
                       ),
                       tabPanel("Vivax breakdown",
                                bslib::card(
                                  bslib::card_header("Vivax components"),
                                  plotOutput("vivaxPlot", height = "420px")
                                )
                       ),
                       tabPanel("Falciparum breakdown",
                                bslib::card(
                                  bslib::card_header("Falciparum components"),
                                  plotOutput("falciPlot", height = "420px")
                                )
                       )
           )
  ),
  
  tabPanel(title = "Parameters", value = "parameters",
           # Put inputs directly in the content area; split into two nice cards
           bslib::layout_columns(col_widths = c(6, 6), gap = "24px",
                                 
                                 bslib::card(
                                   bslib::card_header("Transmission parameters"),
                                   sliderInput("lambda_v", "Vivax force of infection (λᵥ)",
                                               min = 0, max = 2, value = 0.8, step = 0.01),
                                   sliderInput("lambda_F", "Falciparum force of infection (λ_F)",
                                               min = 0, max = 2, value = 0.8, step = 0.01)
                                 ),
                                 
                                 bslib::card(
                                   bslib::card_header("Simulation settings"),
                                   numericInput("patch", "Patch to plot", value = 5, min = 1, max = N, step = 1),
                                   div(class = "d-flex gap-2",
                                       actionButton("go", "Run model", class = "btn btn-primary"),
                                       actionButton("go_to_model2", "View results", class = "btn btn-outline-secondary")
                                   ),
                                   helpText("Tip: Heavy computations run once at startup.")
                                 )
                                 
           )
  ),
  
  tabPanel(title = "About", value = "about",
           bslib::card(
             bslib::card_header("About this app"),
             tags$p("Short blurb about the model, data sources, and authors.")
           )
  ),
  
  # Default selected tab (Home)
  selected = "home"
)

# --- SERVER --------------------------------------------------
server <- function(input, output, session) {
  # keep patch bounded by N from MainModel.R
  observe({ updateNumericInput(session, "patch", max = N) })
  
  # simple nav helpers for your Home buttons
  observeEvent(input$go_to_model,  updateTabsetPanel(session, "topnav", selected = "parameters"))
  observeEvent(input$go_to_model2, updateTabsetPanel(session, "topnav", selected = "model"))
  
  # Run model when user clicks "Run"
  results <- eventReactive(input$go, {
    scen <- list(lambda_v = input$lambda_v, lambda_F = input$lambda_F)
    run_model(parameters, scen, time, initcondrun)
  }, ignoreInit = TRUE)
  
  # ---- Overview plot ----
  output$trajPlot <- renderPlot({
    req(results())
    out <- get_out_matrix(results())
    req(out)
    p <- input$patch
    
    S_ch       <- get_series(out, "S_ch", p)
    Iv_RDT_TP  <- get_series(out, "Iv_RDT_TP_ch", p)
    Iv_RDT_FN  <- get_series(out, "Iv_RDT_FN_ch", p)
    Iv_M_TP    <- get_series(out, "Iv_M_TP_ch",   p)
    Iv_M_FN    <- get_series(out, "Iv_M_FN_ch",   p)
    Av         <- get_series(out, "Av_ch",        p)
    
    IF_RDT_TP  <- get_series(out, "IF_RDT_TP_ch", p)
    IF_RDT_FN  <- get_series(out, "IF_RDT_FN_ch", p)
    IF_M_TP    <- get_series(out, "IF_M_TP_ch",   p)
    IF_M_FN    <- get_series(out, "IF_M_FN_ch",   p)
    AF         <- get_series(out, "AF_ch",        p)
    
    totalVivax  <- Iv_RDT_TP + Iv_RDT_FN + Iv_M_TP + Iv_M_FN + Av
    totalFalcip <- IF_RDT_TP + IF_RDT_FN + IF_M_TP + IF_M_FN + AF
    
    df <- tibble(
      time = out[, 1],
      `S_ch` = S_ch,
      `Vivax (total)` = totalVivax,
      `Falciparum (total)` = totalFalcip
    ) |> pivot_longer(-time, names_to = "series", values_to = "value")
    
    ggplot(df, aes(time, value, color = series)) +
      geom_line(linewidth = 1) +
      labs(
        x = "Year", y = "Population", color = NULL,
        title = paste("Patch", p),
        subtitle = glue("lambda_v = {input$lambda_v}, lambda_F = {input$lambda_F}")
      ) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "top")
  })
  
  # ---- Vivax breakdown ----
  output$vivaxPlot <- renderPlot({
    req(results())
    out <- get_out_matrix(results()); req(out)
    p <- input$patch
    
    parts <- list(
      `Vivax: Av`        = "Av_ch",
      `Vivax: Iv_RDT_TP` = "Iv_RDT_TP_ch",
      `Vivax: Iv_RDT_FN` = "Iv_RDT_FN_ch",
      `Vivax: Iv_M_TP`   = "Iv_M_TP_ch",
      `Vivax: Iv_M_FN`   = "Iv_M_FN_ch"
    )
    
    df <- purrr::map(parts, ~ get_series(out, .x, p)) |>
      list2DF() |>
      mutate(time = out[, 1]) |>
      pivot_longer(-time, names_to = "series", values_to = "value")
    
    ggplot(df, aes(time, value, color = series)) +
      geom_line(linewidth = 1) +
      labs(x = "Year", y = "Population", color = NULL, title = paste("Vivax (Patch", p, ")")) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "top")
  })
  
  # ---- Falciparum breakdown ----
  output$falciPlot <- renderPlot({
    req(results())
    out <- get_out_matrix(results()); req(out)
    p <- input$patch
    
    parts <- list(
      `Falci: AF`        = "AF_ch",
      `Falci: IF_RDT_TP` = "IF_RDT_TP_ch",
      `Falci: IF_RDT_FN` = "IF_RDT_FN_ch",
      `Falci: IF_M_TP`   = "IF_M_TP_ch",
      `Falci: IF_M_FN`   = "IF_M_FN_ch"
    )
    
    df <- purrr::map(parts, ~ get_series(out, .x, p)) |>
      list2DF() |>
      mutate(time = out[, 1]) |>
      pivot_longer(-time, names_to = "series", values_to = "value")
    
    ggplot(df, aes(time, value, color = series)) +
      geom_line(linewidth = 1) +
      labs(x = "Year", y = "Population", color = NULL, title = paste("Falciparum (Patch", p, ")")) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "top")
  })
}



shinyApp(ui, server)
