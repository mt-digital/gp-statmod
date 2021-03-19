#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplotify)

source('util.R')
source('model.R')

# Define UI for application that draws a histogram
ui <- fluidPage(

    titlePanel("Group polarization counterexample generator."),

    h1(getQueryString()[["caseStudy"]]),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
                     # h1("yo"),
            htmlOutput("text"),
        ),

        # Show a plot of the generated distribution
        mainPanel(
            # Initialize parameter controls defaulted to a Schkade result.
             numericInput("latentMean",
                        "Hypothesized latent mean:",
                        min = -1,
                        max = 12,
                        step = 0.1,
                        value = 5.5),
            numericInput("observedPreDelibMean",
                        "Reported, target pre-deliberation mean",
                        min = 1,
                        max = 10,
                        step = 0.1,
                        value = 9.2),
            numericInput("observedPostDelibMean",
                        "Reported, target post-deliberation mean",
                        min = 1,
                        max = 10,
                        step = 0.1,
                        value = 9.4),
            numericInput("minBinValue", "Minimum opinion bin value", 1, step = 1),
            numericInput("maxBinValue", "Maximum opinion bin value", 10, step = 1),
            htmlOutput("caseStudy"),
            br(),
            h4("Large-N model exact calculations"),
            fluidRow(
                column(6,
                    tags$ul(tags$li("Latent Mean: XX"), tags$li("YOOOO"))
                ),
                column(6,
                    h5("pre:"),
                    plotOutput("preBarplot"),
                    h5("post:"),
                    plotOutput("postBarplot")
                )
             ),
            br(),
            h5("Confirmatory simulations for empirical N"),
           # TODO: polts n stuff
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

    # Sets up case study if there is one.
    caseStudy <- 
        reactive({
            # Check if there is a caseStudy query parameter, designed to be
            # passed when user clicks on nav element in sidebar.
            if (!is.null(getQueryString()$caseStudy))
            {
                caseStudy <- getQueryString()$caseStudy
                setCaseStudyValues(caseStudy, session)

                return (caseStudy)
            }
            else
            {
                # Defaulting to Schkade result for now.
                return ("Schkade, et al., (2010) - Liberals on Affirmative Action")
            }
        })

    output$caseStudy <- renderText({caseStudy()})

    kVec <- reactive({input$minBinValue:input$maxBinValue})

    # Use a guess of 1.5 for latentSD.
    latentPreSD <- reactive({
        solveForLatentSD(kVec(), input$latentMean, input$observedPreDelibMean, 1.5)
    })
    latentPostSD <- reactive({
        solveForLatentSD(kVec(), input$latentMean, input$observedPostDelibMean, 1.5)
    })

    probVec <- reactive({makeProbVec(kVec(), input$latentMean, latentPreSD())})
    output$preBarplot <- renderPlot({
        barplot(probVec())
    })

    probVecPost <- reactive({
        makeProbVec(kVec(), input$latentMean, latentPostSD())
    })
    output$postBarplot <- renderPlot({
        barplot(probVecPost())
    })

    observe({
        val <- input$minBinValue
        updateSliderInput(session, "latentMean", min = val)
    })

    observe({
        val <- input$maxBinValue
        updateSliderInput(session, "latentMean", max = val)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
