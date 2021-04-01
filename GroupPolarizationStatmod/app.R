#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(ggplot2)
library(ggplotify)
library(reshape2)
library(shiny)

source("model.R")
source("numerical.R")


# Define UI for application that draws a histogram
ui <- function(request) { fluidPage(

    titlePanel("Group polarization counterexample generator."),

    h1(getQueryString()[["caseStudy"]]),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
                     # h1("yo"),
            htmlOutput("text"),
            width = 3
        ),

        # Show a plot of the generated distribution
        mainPanel(
        # Initialize parameter controls defaulted to a Schkade result.
            fluidRow(
                column(6, 
                    numericInput("minBinValue", "Minimum opinion bin value", 
                                 1, step = 1),
                    numericInput("maxBinValue", "Maximum opinion bin value", 
                                 10, step = 1),
                    fluidRow(
                        column(6,
                            numericInput("step", "Hillclimb solver step size", 
                                         0.01, step = 0.001)
                        ),
                        column(6,
                            numericInput("successTol", "Hillclimb success threshold", 
                                         1e-2)
                        )
                    )
                ),
                column(6,
                    numericInput("latentMean",
                               "Hypothesized latent mean:",
                               # min = -1,
                               # max = 12,
                               step = 0.1,
                               value = 5.5),
                    numericInput("observedPreDelibMean",
                                "Reported, target pre-deliberation mean",
                                # min = 1,
                                # max = 10,
                                step = 0.1,
                                value = 9.2),
                    numericInput("observedPostDelibMean",
                                "Reported, target post-deliberation mean",
                                # min = 1,
                                # max = 10,
                                step = 0.1,
                                value = 9.4)
                ),
            ),
            # htmlOutput("caseStudy"),
            br(),
            h4("Large-N model exact calculations"),
            plotOutput("largeNBarplot"),
            br(),
            h4("Confirmatory simulations for empirical N"),
           # TODO: polts n stuff
        )
    )
)
}

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
    latentPreSDResult <- reactive({
        solveForLatentSD(kVec(), input$latentMean, input$observedPreDelibMean, 1.5,
                         step = input$step)
    })
    latentPostSDResult <- reactive({
        solveForLatentSD(kVec(), input$latentMean, input$observedPostDelibMean, 1.5,
                         step = input$step)
    })

    probVecPre <- reactive({
        makeProbVec(kVec(), input$latentMean, latentPreSDResult()[1])
    })

    output$preBarplot <- renderPlot({
        barplot(probVec())
    })

    probVecPost <- reactive({
        # 
        makeProbVec(kVec(), input$latentMean, latentPostSDResult()[1])
    })
    output$postBarplot <- renderPlot({
        barplot(probVecPost())
    })

    largeNPlotDF <- reactive({
        melt(data.frame(probVecPre = probVecPre(), 
                        probVecPost = probVecPost(), 
                        kVec = kVec()), 
            id.vars = 'kVec')
    })

    sdObsPre <- reactive({ sdObs(kVec(), probVecPre()) })
    sdObsPost <- reactive({ sdObs(kVec(), probVecPost()) })
    # Build the title based on the solved-for latent pre- and post-deliberation
    # opinion SDs, and report the squared error and number of steps required.
    plotTitle <- reactive({

        fmtVal <- function(val) { format(val, digits=3) }
        title <-  paste(
            # "Solved ",
            "latPreSD=", fmtVal(latentPreSDResult()[1]), 
            ", latPostSD=", fmtVal(latentPostSDResult()[1]), "\n",
            # "latPreSD=", fmtVal(latentPreSDResult()[1]^0.5), 
            # ", latPostSD=", fmtVal(latentPostSDResult()[1]^0.5), "\n",
            "simObsPreSD=", fmtVal(sdObsPre()), 
            ", simObsPostSD=", fmtVal(sdObsPost()), 
            ".\nSq. Error: pre = ", fmtVal(latentPreSDResult()[2]),
            ", post = ", fmtVal(latentPostSDResult()[2]),
            sep = ""
        )

        # If the pre- or post-SD result reports a failure to find a solution
        # then report this in the main figure.
        hillclimbSuccess <- (
            # The second element in the hillclimb results is the squared error.
            latentPreSDResult()[2] < input$successTol && latentPostSDResult()[2] < input$successTol
        )
        if (!hillclimbSuccess)
        {
            title <- paste(title, "(FAILED TO FIND SD SATISFYING CONSTRAINTS)", 
                           sep = "\n")
        }

        return (title)
    })
    output$largeNBarplot <- renderPlot({

        ggplot(largeNPlotDF(), aes(x = kVec, y = value, fill=variable)) + 
            geom_bar(stat='identity', position='dodge') +
            scale_x_continuous(breaks = kVec()) +
            xlab('Reported opinion') + ylab('Frequency') +
            ggtitle(plotTitle())
    })

    observe({
        reactiveValuesToList(input)
        session$doBookmark()
    })

    onBookmarked(function(url) { updateQueryString(url) })
}

# Run the application 
shinyApp(ui = ui, server = server, enableBookmarking = "url")
