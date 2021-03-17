#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

source('util.R')

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
           # printOutput(
           plotOutput("preBarplot"),
           plotOutput("postBarplot")
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

    # output$standardDeviationReport <- renderText

    output$preBarplot <- renderPlot({
        largeNBarplot(input, "pre")
    })

    output$post <- renderPlot({
        largeNBarplot(input, "post")
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
