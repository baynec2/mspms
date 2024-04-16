library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(mspms)
library(dplyr)

ui <- dashboardPage(#skin = "midnight",
    dashboardHeader(title = "MSP-MS"),
    dashboardSidebar(
      # Defining the sidebar menu
      sidebarMenu(
        # Making the tabs
        menuItem("About", tabName = "about", icon = icon("magnifying-glass-chart")),
        menuItem("Output", tabName = "output", icon = icon("database")),
        menuItem("Stats", tabName = "stats", icon = icon("star-of-life")),
        menuItem("DataViz", tabName = "viz", icon = icon("chart-simple"))
      ),
      #Putting the place to load the files
      fileInput("upload1", "Upload PEAKS LFQ file",accept = ".csv"),
      fileInput("upload2", "Upload PEAKS ID file",accept = ".csv"),
      fileInput("upload3", "Upload design matrix",accept = ".csv")
      ),

    dashboardBody(
      ### About page, put instructions here ##
      tabItems(
        tabItem(tabName = "about",
               htmlOutput("about")),
        tabItem(tabName = "output",
                h1("normalized data"),
                fluidRow(
                  box(DT::DTOutput('processed_data'),width = 12),
                downloadButton(
                  "downloadData",
                  label = "Download"
                ))),
        tabItem(tabName = "stats",
                fluidRow(
                  box(DT::DTOutput('ttest_results'),width = 6),
                  box(DT::DTOutput('anova_results'),width = 6)),
                downloadButton(
                  "downloadttest",
                  label = "Download t-test"
                  ),
                downloadButton(
                  "downloadanova",
                  label = "Download anova"
                )
                ),
        tabItem(tabName = "viz",
                fluidRow(
                  box(plotOutput("plot1"),width = 6),
                  box(plotOutput("plot2"),width = 6)
                ),
                fluidRow(
                  box(plotOutput("plot3"),width = 12)

                )
                )
        )
    )
)



# Define server logic
server <- function(input, output){

  # Rendering the about page
  output$about <- renderUI({
    includeHTML("./about.html")})


  # Reading in the design matrix
  design_matrix = reactive({
    readr::read_csv(input$upload3$datapath)
  })

  # Processing the data normalization data
  processed_data = reactive({
        prepared_data =  mspms::prepare_for_normalyzer(input$upload1$datapath,input$upload2$datapath)
        # mspms workflow
        design_matrix = design_matrix()
        normalyzed_data = mspms::normalyze(prepared_data,design_matrix)
        outliers = mspms::handle_outliers(normalyzed_data,design_matrix)
        imputed = mspms::impute(outliers)
        joined_with_library = mspms::join_with_library(imputed)
        final_data = mspms::add_cleavages(joined_with_library)
        final_data
        })


  #Downloading processed data
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("processed_data", ".csv", sep = "")
    },
    content = function(file){
      readr::write_csv(processed_data(), file)
    }
  )

  output$processed_data = DT::renderDT(
    processed_data(),
    options = list(scrollX = TRUE))

  # Preparing for stats

  prepared_for_stats = reactive({
      processed_data() %>%
      mspms::prepare_for_stats(design_matrix())
  })

  anova = reactive({
    prepared_for_stats() %>%
    mspms::mspms_anova()
  })


  #Downloading anova data
  output$downloadanova <- downloadHandler(
    filename = function() {
      paste("anova_results", ".csv", sep = "")
    },
    content = function(file){
      readr::write_csv(anova(), file)
    }
  )

  output$anova_results = DT::renderDT(
    anova(),
    options = list(scrollX = TRUE))

  # Doing the t-tests

  ttests = reactive({
    prepared_for_stats() %>%
    mspms::mspms_t_tests()
  })


  #Downloading processed data
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("ttest_results", ".csv", sep = "")
    },
    content = function(file){
      readr::write_csv(ttests(), file)
    }
  )

  output$ttest_results = DT::renderDT(
    ttests(),
    options = list(scrollX = TRUE))


  #doing some data visualizations.

  # Volcano plots
  log2fct = reactive({
     prepared_for_stats() %>%
     mspms::log2fc_t_test()
   })

  output$plot1 = renderPlot({

    prepared_for_stats() %>%
      mspms::plot_pca()

    })


  output$plot2 = plotly::renderPlotly({

    prepared_for_stats() %>%
      mspms::plot_heatmap()



    })

  output$plot3 = renderPlot({

    log2fct() %>%
      ggplot(aes(x = log2fc,y = -log10(p.adj)))+
      geom_point()+
      geom_hline(yintercept = -log10(0.05),linetype = "dashed",color = "red")+
      geom_vline(xintercept = 3, linetype = "dashed",color = "red")+
      geom_vline(xintercept = -3, linetype = "dashed",color = "red")+
      theme_minimal()+
      labs(x = "Log2 Fold Change",y = "-log10(p value)")+
      facet_wrap(~comparison,scales = "free")

  })



  }








# Run the application
shinyApp(ui = ui, server = server)
