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
library(readxl)
library(shiny)

source("model.R")
source("numerical.R")

STUDIES <- read_excel("data/CaseStudies_ProtoDB-TEST.xlsx", "ArticleStudiesMinimal")
# STUDIES <- read_excel("data/CaseStudies_ProtoDB.xlsx", "ArticleStudiesMinimal")
FULL_TAGS <- sort(paste(STUDIES$ArticleTag, STUDIES$TreatmentTag, sep=" - "))

# Define UI for application that draws a histogram
ui <- function(request) { fluidPage(

    tags$head(tags$style(HTML("
        #plausible-div {
            text-align: center;
        }
        #save-btn-div {
            text-align: center;
        }
        #saveBtn {
            font-size: 25px;   
            width: 50%;
            color: white;
            background-color: dodgerblue;
        }
    
    "))),
  
    titlePanel("Group polarization counterexample generator."),

    h3(textOutput("treatmentTag")),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          selectInput("treatmentTag", "Choose a study within an article:",
                      choices = FULL_TAGS, size=25, selectize = FALSE),
            width = 3
        ),

        # Show a plot of the generated distribution
        mainPanel(
        # Initialize parameter controls defaulted to a Schkade result.
            fluidRow(
                column(6,
                    numericInput("ObservedMeanPre",
                                "Reported, target pre-deliberation mean",
                                step = 0.01,
                                value = 9.2),
                    numericInput("ObservedMeanPost",
                                "Reported, target post-deliberation mean",
                                step = 0.01,
                                value = 9.4),
                    numericInput("LatentMean",
                               "Hypothesized latent mean:",
                               step = 0.01,
                               value = 5.5),
                ),
                column(6, 
                    numericInput("MinBinValue", 
                                 "Minimum opinion bin value", 
                                 1, step = 1),
                    numericInput("MaxBinValue", 
                                 "Maximum opinion bin value", 
                                 10, step = 1),
                    fluidRow(
                        column(6,
                            numericInput("step", 
                                         "Hillclimb solver step size", 
                                         0.01, step = 0.01)
                        ),
                        column(6,
                            numericInput("successTol", 
                                         "Hillclimb success threshold", 
                                         0.05, step = 0.005)
                        )
                    )
                ),
                # div(
                #     id = "plausible-div",
                #     checkboxInput("Plausible", "Plausibly false positive?", FALSE)
                # ),
            ),

            # htmlOutput("caseStudy"),
            br(),
            h4("Large-N model exact calculations"),
            plotOutput("largeNBarplot"),
            br(),
           # TODO: SAVE button 
        )
    ),
    # print(tags),
    tags$head(tags$script(src = "message-handler.js")),
    div(id = "save-btn-div", 
        actionButton("saveBtn", "Save")
    )
)
}


########## SERVER ########## 

# Define server logic required to draw a histogram
server <- function(input, output, session) {

    STUDIES <- read_excel("data/CaseStudies_ProtoDB-TEST.xlsx", "ArticleStudiesMinimal")
    FULL_TAGS <- sort(paste(STUDIES$ArticleTag, STUDIES$TreatmentTag, sep=" - "))

    # output$caseStudy <- renderText({caseStudy()})
    treatmentSplit <- reactive({strsplit(input$treatmentTag, " - ")})

    article <- reactive({treatmentSplit()[[1]][1]})
    treatment <- reactive({treatmentSplit()[[1]][2]})
    treatmentRow <- reactive({STUDIES[STUDIES$TreatmentTag == treatment(), ]})
    
    # print({treatmentRow()})

    kVec <- reactive({treatmentRow()$MinBinValue:treatmentRow()$MaxBinValue})
    output$MinBinValue <- reactive({treatmentRow()$MinBinValue})
    output$MaxBinValue <- reactive({treatmentRow()$MaxBinValue})
    output$treatmentTag <- renderText({
        paste("Currently analyzing treatment:", input$treatmentTag)
    })

    # Use a guess of 1.5 for latentSD.
    # solveForLatentSD returns three values: xcurr, ycurr, its
    latentPreSDResult <- reactive({
        solveForLatentSD(kVec(), input$LatentMean, input$ObservedMeanPre, 1.5,
                         step = input$step)
    })
    latentPostSDResult <- reactive({
        solveForLatentSD(kVec(), input$LatentMean, input$ObservedMeanPost, 1.5,
                         step = input$step)
    })

    probVecPre <- reactive({
        makeProbVec(kVec(), input$LatentMean, latentPreSDResult()[1])
    })

    output$preBarplot <- renderPlot({
        barplot(probVec())
    })

    probVecPost <- reactive({
        # 
        makeProbVec(kVec(), input$LatentMean, latentPostSDResult()[1])
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

    # observe({
    #     reactiveValuesToList(input)
    #     session$doBookmark()
    # })
    observeEvent(input$treatmentTag, {
        print(input$treatmentTag)
        treatmentSplit <- strsplit(input$treatmentTag, " - ")
        print(treatmentSplit)
        article <- treatmentSplit[[1]][1]
        treatment <- treatmentSplit[[1]][2]
        print(article)
        print(treatment)
        # print(STUDIES[STUDIES$TreatmentTag == treatment, ])
        treatmentRow <- STUDIES[STUDIES$TreatmentTag == treatment, ]
        print(treatmentRow)
        updateNumericInput(session, 
                           "LatentMean", 
                           value = treatmentRow$LatentMean)
        updateNumericInput(session, 
                           "ObservedMeanPre", 
                           value = treatmentRow$ObservedMeanPre)
        updateNumericInput(session, 
                           "ObservedMeanPost", 
                           value = treatmentRow$ObservedMeanPost)
        updateNumericInput(session, 
                           "MinBinValue", 
                           value = treatmentRow$MinBinValue)
        updateNumericInput(session, 
                           "MaxBinValue", 
                           value = treatmentRow$MaxBinValue)
        updateNumericInput(session, 
                           "step", 
                           value = treatmentRow$HillclimbStepSize)
        updateNumericInput(session, 
                           "successTol", 
                           value = treatmentRow$HillclimbSuccessThreshold)
        updateCheckboxInput(session,
                            "Plausible",
                            value = treatmentRow$Plausible)
        print("")
        print("*YOO whatup*")
    })
    
    observeEvent(input$saveBtn, {
        treatmentSplit <- strsplit(input$treatmentTag, " - ")
        print(treatmentSplit)
        article <- treatmentSplit[[1]][1]
        treatment <- treatmentSplit[[1]][2]
        print(article)
        print(treatment)
        # print(STUDIES[STUDIES$TreatmentTag == treatment, ])
        treatmentRow <- STUDIES[STUDIES$TreatmentTag == treatment, ]
        treatmentRow <- STUDIES[STUDIES$TreatmentTag == treatment, ]
        treatmentRow$ObservedMeanPre <- input$ObservedMeanPre
        treatmentRow$ObservedMeanPost <- input$ObservedMeanPost
        treatmentRow$LatentMean <- input$LatentMean
        print(treatmentRow)
    })

    onBookmarked(function(url) { updateQueryString(url) })
}

# Run the application 
shinyApp(ui = ui, server = server, enableBookmarking = "url")
