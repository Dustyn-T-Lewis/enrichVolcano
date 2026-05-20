# enrichVolcano interactive GUI. Launch with enrichVolcano::ev_app().
# Tune input, databases, and thresholds, then render the volcano-in-ring live.

library(shiny)
library(enrichVolcano)

reg <- tryCatch(list_databases(), error = function(e) NULL)
db_choices <- if (!is.null(reg)) stats::setNames(reg$name, reg$display_name) else
  c("Hallmark" = "hallmark", "Reactome" = "reactome", "GO BP" = "go_bp")

species_choices <- c("human", "mouse", "rat", "zebrafish", "fly", "yeast", "pig")

ui <- fluidPage(
  titlePanel("enrichVolcano — volcano-in-ring builder"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      radioButtons("source", "Data source",
                   c("Bundled YvO example" = "example", "Upload CSV" = "upload")),
      conditionalPanel(
        "input.source == 'upload'",
        fileInput("file", "Differential abundance CSV", accept = ".csv"),
        helpText("Any limma / proteoDA / DEP / DESeq2 / edgeR / wide-suffix table.")
      ),
      selectInput("contrast", "Contrast", choices = NULL),
      selectInput("species", "Species", choices = species_choices),
      selectizeInput("databases", "Databases", choices = db_choices,
                     selected = c("hallmark", "reactome", "go_bp"),
                     multiple = TRUE),
      sliderInput("enrich_padj", "Pathway FDR cutoff", 0.01, 0.5, 0.05, 0.01),
      sliderInput("max_terms", "Max pathways shown", 4, 25, 14, 1),
      fluidRow(
        column(6, numericInput("min_size", "Min set size", 15, 1, 500)),
        column(6, numericInput("max_size", "Max set size", 500, 10, 5000))
      ),
      selectInput("p_method", "Volcano y-axis",
                  c("pi-score (Eq.2)" = "pi_eq2", "pi-score (Eq.1)" = "pi_eq1",
                    "raw p" = "raw_p", "adjusted p" = "adj_p")),
      sliderInput("dedup_cutoff", "Jaccard dedup cutoff", 0.1, 0.9, 0.5, 0.05),
      actionButton("build", "Build figure", class = "btn-primary"),
      tags$hr(),
      downloadButton("dl_pdf", "Download PDF"),
      downloadButton("dl_png", "Download PNG")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Figure", plotOutput("plot", height = "640px")),
        tabPanel("Enriched pathways", tableOutput("tbl")),
        tabPanel("Input summary", verbatimTextOutput("summary"))
      )
    )
  )
)

server <- function(input, output, session) {

  raw_data <- reactive({
    if (input$source == "example") {
      fp <- system.file("extdata", "examples", "yvo_tidy.rds",
                        package = "enrichVolcano")
      validate(need(nzchar(fp), "Bundled example not found."))
      readRDS(fp)
    } else {
      req(input$file)
      utils::read.csv(input$file$datapath, check.names = FALSE,
                      stringsAsFactors = FALSE)
    }
  })

  validated <- reactive({
    d <- raw_data()
    tryCatch(ev_validate(d),
             error = function(e) {
               validate(need(FALSE, paste("Could not read input:",
                                           conditionMessage(e))))
             })
  })

  observeEvent(validated(), {
    ctrs <- unique(validated()$contrast)
    updateSelectInput(session, "contrast", choices = ctrs, selected = ctrs[1])
  })

  output$summary <- renderText({
    v <- validated()
    paste0(
      "Detected format: ", attr(v, "ev_source"), "\n",
      "Proteins: ", length(unique(v$gene)), "\n",
      "Contrasts: ", length(unique(v$contrast)), "\n  ",
      paste(unique(v$contrast), collapse = "\n  ")
    )
  })

  result <- eventReactive(input$build, {
    req(input$contrast, length(input$databases) > 0)
    withProgress(message = "Running enrichment...", {
      enrich_volcano(
        raw_data(), contrast = input$contrast, species = input$species,
        databases = input$databases, p_method = input$p_method,
        enrich_mode = "fgsea", enrich_padj = input$enrich_padj,
        dedup = list(method = "jaccard", cutoff = input$dedup_cutoff,
                     scope = "within_db"),
        ring = list(max_terms = input$max_terms, order_by = "padj",
                    magnitude = "neg_log_padj", color = "nes")
      )
    })
  })

  output$plot <- renderPlot({
    print(result())
  }, res = 96)

  output$tbl <- renderTable({
    enr <- attr(result(), "ev_data")$enrichment
    req(nrow(enr) > 0)
    enr <- enr[order(enr$padj), c("database", "pathway", "NES", "padj", "size")]
    utils::head(enr, 40)
  })

  out_plot <- function(file, device, w = 9, h = 9) {
    ggplot2::ggsave(file, result(), device = device, width = w, height = h)
  }
  output$dl_pdf <- downloadHandler(
    filename = function() paste0("volcano_ring_", input$contrast, ".pdf"),
    content = function(file) out_plot(file, "pdf")
  )
  output$dl_png <- downloadHandler(
    filename = function() paste0("volcano_ring_", input$contrast, ".png"),
    content = function(file) out_plot(file, "png")
  )
}

shinyApp(ui, server)
