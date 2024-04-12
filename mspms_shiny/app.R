library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(mspms)

ui <- dashboardPage(skin = "midnight",
    dashboardHeader(title = "MSP-MS"),
    dashboardSidebar(
      sidebarMenu(
        # Making the tabs
        menuItem("About", tabName = "about", icon = icon("magnifying-glass-chart")),
        menuItem("QC Plots", tabName = "QC", icon = icon("chart-simple"))
      ),
       fileInput("upload1", "Upload PEAKS LFQ file",accept = ".csv"),
       fileInput("upload2", "Upload PEAKS ID file",accept = ".csv"),
       downloadButton(
        "downloadData",
        label = "Download",
        icon = shiny::icon("download")
      )
    ),

    dashboardBody(
      ### About page, put instructions here ##
      tabItems(
        tabItem(tabName = "about",
                fluidPage(
                  tags$iframe(src = './about.html',
                              width = '100%', height = '800px',
                              frameborder = 0, scrolling = 'auto'
                  )
                )),
        tabItem(tabName = "QC",
                h1("Load_files"),
                fluidRow(
                  box()
                  )),
                  )
                )

        )
    )
)

# Define server logic
server <- function(input, output){


prepared_data = reactive({
    req(input$upload1)
    req(input$upload2)

  # mspms workflow
  # Prepare the data for normalyzer analysis
  prepared_data = mspms::prepare_for_normalyzer(input$upload1$datapath,input$upload2$datapath)})
  # Extracting design matrix
  design_matrix = mspms::extract_design_matrix(prepared_data)
  normalyzed_data = mspms::normalyze(prepared_data,design_matrix)
  outliers = mspms::handle_outliers(normalyzed_data)
  imputed = mspms::impute(outliers)
  joined_with_library = mspms::join_with_library(imputed)


  final_data = mspms::add_cleavages(joined_with_library)

  #Downloading processed data
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("processed_data", ".csv", sep = "")
    },
    content = function(file) {
      readr::write_csv(final_data, file)
    }
  )





}


# Run the application
shinyApp(ui = ui, server = server)
