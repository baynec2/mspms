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
        menuItem("Output", tabName = "output", icon = icon("chart-simple"))
      ),
      #Putting the place to load the files
       fileInput("upload1", "Upload PEAKS LFQ file",accept = ".csv"),
       fileInput("upload2", "Upload PEAKS ID file",accept = ".csv")
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
                )))
                  )
                  )
    )


# Define server logic
server <- function(input, output){

  processed_data = reactive({
        prepared_data =  mspms::prepare_for_normalyzer(input$upload1$datapath,input$upload2$datapath)
        # mspms workflow
        design_matrix = mspms::extract_design_matrix(prepared_data)
        normalyzed_data = mspms::normalyze(prepared_data,design_matrix)
        outliers = mspms::handle_outliers(normalyzed_data)
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


  output$about <- renderUI({
    includeHTML("./about.html")
  })
  }






# Run the application
shinyApp(ui = ui, server = server)
