## Author: Taylor Falk
## tfalk@bu.edu
## BU BF591
## Assignment 7

# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(bslib)
library(ggplot2)
library(colourpicker) # you might need to install this

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("BF591 Assignment 7"),
  h4("To use this application, download the CSV deseq_res.csv from the data directory of this app's repository."),
  
  sidebarLayout(
    sidebarPanel(
      # Add elements for user inputs (e.g., file upload, sliders) here
      fileInput("file1", "Load differential expression results"),
      HTML("A volcano plot can be generated with <b>\"log2 fold-change\"</b> on the x-axis and <b>\"p-adjusted\"</b> on the y-axis."),
      radioButtons("xAxisRadio", 
                   label = "Choose the column for the x-axis",
                   choices = list("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"), 
                   selected = "log2FoldChange"),
      radioButtons("yAxisRadio", 
                   label = "Choose the column for the y-axis",
                   choices = list("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"), 
                   selected = "padj"),
      colourInput("color1", "Base point color", value = "#22577A"),
      colourInput("color2", "Highlight point color", value = "#FFCF56"),
      sliderInput("magnitudeSlider", "Select the magnitude of the p adjusted coloring:", min = -300, max = 0, value = -150),
      actionButton("submit_btn", "Plot", icon("car"), width = "100%")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Plot", plotOutput("volcano", height = "675px")),
        tabPanel("Table", tableOutput("table"))
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    #' load_Data
    #'
    #' @details Okay this one is a little weird but bear with me here. This is 
    #' still a "function", but it will take no arguments. The `reactive({})` bit 
    #' says "if any of my inputs (as in, input$...) are changed, run me again". 
    #' This is useful when a user clicks a new button or loads a new file. In 
    #' our case, look for the uploaded file's datapath argument and load it with 
    #' read.csv. Return this data frame in the normal return() style.
    load_data <- reactive({
      data <- read.csv(input$file1$datapath)
      return(data)
    })
    
    #' Volcano plot
    #'
    #' @param dataf The loaded data frame.
    #' @param x_name The column name to plot on the x-axis
    #' @param y_name The column name to plot on the y-axis
    #' @param slider A negative integer value representing the magnitude of
    #' p-adjusted values to color. Most of our data will be between -1 and -300.
    #' @param color1 One of the colors for the points.
    #' @param color2 The other colors for the points. Hexadecimal strings: "#CDC4B5"
    #'
    #' @return A ggplot object of a volcano plot
    #' @details I bet you're tired of these plots by now. Me too, don't worry.
    #' This is _just_ a normal function. No reactivity, no bells, no whistles. 
    #' Write a normal volcano plot using geom_point, and integrate all the above 
    #' values into it as shown in the example app. The testing script will treat 
    #' this as a normal function.
    #' 
    #' !!sym() may be required to access column names in ggplot aes().
    #'
    #' @examples volcano_plot(df, "log2fc", "padj", -100, "blue", "taupe")
    volcano_plot <-
        function(dataf, x_name, y_name, slider, color1, color2) {
          plot <- ggplot(dataf, aes_string(x = x_name, y = sprintf("-log10(%s) + 1e-10", y_name), color = "factor(padj < 10^slider)")) +
            geom_point(aes(), shape = 19) +
            scale_color_manual(values = c("FALSE" = color1, "TRUE" = color2, "NA" = "gray"), name = sprintf("padj < 1 × 10^%d", slider)) +
            labs(x = x_name, y = sprintf("-log10(%s)", y_name)) +
            theme_minimal() +
            theme(legend.position = "bottom", legend.direction = "horizontal") 
          
          return(plot)
        }
    
    #' Draw and filter table
    #'
    #' @param dataf Data frame loaded by load_data()
    #' @param slider Negative number, typically from the slider input.
    #'
    #' @return Data frame filtered to p-adjusted values that are less than 
    #' 1 * 10^slider, columns for p-value and p-adjusted value have more digits 
    #' displayed.
    #' @details Same as above, this function is a standard R function. Tests will 
    #' evaluate it normally. Not only does this function filter the data frame to 
    #' rows that are above the slider magnitude, it should also change the format 
    #' of the p-value columns to display more digits. This is so that it looks 
    #' better when displayed on the web page. I would suggest the function 
    #' `formatC()`
    #'
    #' @examples draw_table(deseq_df, -210)
    #'    X  baseMean     log2FC     lfcSE      stat       pvalue         padj
    #'gene1 11690.780   9.852926 0.2644650  37.25607 8.45125e-304 1.54472e-299
    #'gene2  3550.435  -6.183714 0.1792708 -34.49369 9.97262e-261 9.11398e-257
    draw_table <- function(dataf, slider) {
      # remove rows with NAs
      dataf <- na.omit(dataf)
      # filter the data frame based on the slider magnitude
      filtered_data <- dataf[dataf$padj < 10^slider, ]
      # format p-value columns to display more digits
      filtered_data$pvalue <- formatC(filtered_data$pvalue, format = "e", digits = 3)
      filtered_data$padj <- formatC(filtered_data$padj, format = "e", digits = 3)

      return(filtered_data)
    }

    #' These outputs aren't really functions, so they don't get a full skeleton,
    #' but use the renderPlot() and renderTable() functions to return() a plot
    #' or table object, and those will be displayed in your application.
    output$volcano <- renderPlot({
      req(input$submit_btn)
      data <- load_data()
      plot <- volcano_plot(data, input$xAxisRadio, input$yAxisRadio, input$magnitudeSlider, input$color1, input$color2)
      return(plot)
    })

    # Same here, just return the table as you want to see it in the web page
    output$table <- renderTable({
      req(input$submit_btn)
      data <- load_data()
      filtered_table <- draw_table(data, input$magnitudeSlider)
      
      return(filtered_table)
    })
}

# Run the application
shinyApp(ui = ui, server = server)
