#setwd("G:\\Shared drives\\Ehsan PhD work\\Codes\\Git\\Iterative-removals-IC-cost\\")
source("CRTVarAdj_func.R", local=TRUE)
source("IterRemCst_func.R", local=TRUE)
source("FigGenDf_func.R", local=TRUE)

library("shiny")
library("shinyBS")
library("ggplot2")
library("reshape2")
library("plyr")
library("swCRTdesign")
library("matrixcalc")
library("scales")
library("tidyverse")
library("shinythemes")
library("Matrix")
library("plotly")
library("RColorBrewer")
library("tidyr")
library("memoise")

ui <- fluidPage(
    
    tags$head(includeHTML(("google-analytics.html"))),
  
    titlePanel(h1("Cost-efficient incomplete stepped wedge designs",
                  h2("Using an iterative approach"),
                  h3(""))),
    sidebarLayout(
        sidebarPanel(
            h4("Trial configuration"),
            sliderInput(inputId = "Tp", label = "Number of periods:",
                        min = 5, max = 20, value = 6,  step = 1,ticks = FALSE),
            numericInput("m",
                         "Number of subjects in each cluster-period, m:",
                         min = 1,
                         max=1000,
                         step = 1,
                         value = 7),
            numericInput("N",
                         "Number of clusters per each sequence, N:",
                         min = 1,
                         max=1000,
                         step = 1,
                         value = 7),
            hr(),
            
            h4("Correlation parameters"),
            
            numericInput("rho0", "Intra-cluster correlation",
                         min = 0,
                         max=1,
                         step = 0.001,
                         value = 0.05),
            radioButtons("type", label = ("Allow for decay correlation"),
                         choices = list("Yes" = 1, "No" = 0), selected = 1),
            numericInput("r",
                         "Cluster auto-correlation",
                         min = 0,
                         max=1,
                         step = 0.05,
                         value = 0.95),
            hr(),
            
            h4("Effect Size"),
            
            sliderInput("effsiz", "Effect size:",
                        min = 0.05, max = 1.0,
                        value = 0.27, step = 0.01,ticks = FALSE),
            
            hr(),
            
            h4("Minimal Acceptable Power"),
            
            numericInput("accept_pwr", "Acceptable power (%):",
                         min = 0, max = 100,
                         value = 80, step = 1),
            
            hr(),
            
            h4("Cost components"),
            
            # Input: Interval, cost per cluster
            numericInput("c", "Cost per cluster",
                         min = 100, max = 100000,
                         value = 2500, step = 100),
            
            # Input: Interval, cost per subject receiving intervention 
            numericInput("p", "Cost per subject receiving intervention ",
                         min = 0, max = 5000,
                         value = 140, step = 10),
            
            # Input: Interval, cost per subject receiving control
            numericInput("pprim", "Cost per subject receiving control",
                         min = 0, max = 5000,
                         value = 80, step = 10),
            
            # Input: Interval, Cost of restarting an intervention condition
            numericInput("g", "Restart cost (intervention)",
                         min = 0, max = 50000,
                         value = 230, step = 50),
            bsTooltip("g","Cost of restarting an intervention condition","right"),
            # Input: Interval, Cost of restarting a control condition
            numericInput("gprim", "Restart cost (control)",
                         min = 0, max = 50000,
                         value = 0, step = 50),
            bsTooltip("gprim","Cost of restarting a control condition","right"),
            # clicks the button
            actionButton("update", "Update"),
        ),
        mainPanel(
            tabsetPanel(
                tabPanel("The iterative removal of cluster-period cells",
                         uiOutput("plotheader1a"), uiOutput("plotheader1b"),
                         plotlyOutput("varREMplot"),
                         textOutput("ICremtext")
                         #textOutput("errorText")
                ),
                tabPanel("Relative cost efficiency",
                         uiOutput("plotheader6a"), uiOutput("plotheader6b"),
                         plotlyOutput("RCEplot"),
                         textOutput("ICRCEtext")
                ),
                tabPanel("Power",
                         uiOutput("plotheader4a"), uiOutput("plotheader4b"),
                         plotlyOutput("Powplot"),
                         textOutput("ICpowtext")
                ),
                tabPanel("Variance",
                         uiOutput("plotheader2a"), uiOutput("plotheader2b"),
                         plotlyOutput("Varsplot"),
                         textOutput("ICvartext")
                ),
                tabPanel("Precision loss",
                         uiOutput("plotheader3a"), uiOutput("plotheader3b"),
                         plotlyOutput("Prelossplot"),
                         textOutput("ICRvartext")
                ),
                
                tabPanel("Total cost",
                         uiOutput("plotheader5a"), uiOutput("plotheader5b"),
                         plotlyOutput("Cstplot"),
                         textOutput("ICcsttext")
                ),
                tabPanel("Cost vs Variance",
                         uiOutput("plotheader7a"), uiOutput("plotheader7b"),
                         plotlyOutput("CstVarplot"),
                         textOutput("ICcstVartext")
                ),
                tabPanel("Contact us",
                         verbatimTextOutput("text")
                )
            )
        )
    )
)
server <- function(input, output, session) {
    output$ICremtext <- renderText({
        "Information content of progressively reduced stepped wedge designs"
    })
    output$ICvartext <- renderText({
        "..."
    })
    output$ICRvartext <- renderText({
        ""
    })
    output$ICPowtext <- renderText({
        "..."
    })
    output$ICCsttext <- renderText({
        "..."
    })
    output$ICRCEtext <- renderText({
        ""
    })
    output$ICcstVartext <- renderText({
        "..."
    })
    output$text <- renderText({
        paste(
            "For more information please contact us:",
            "Monash University",
            "School of Public Health and Preventive Medicine",
            "553 St Kilda Road", 
            "Melbourne VIC 3004", 
            "Australia",
            "ehsan.rezaeidarzi@monash.edu", sep="\n")
    })
    
    values <- reactiveValues(
        Tp = 6,
        m = 7,
        N=7,
        rho0 = 0.05,
        r=0.95,
        type=1,
        effsiz=0.27,
        accept_pwr=80,
        c=2500,
        p=140,
        pprim=80,
        g=230,
        gprim=0
    )
    
    observeEvent(input$update, {
        values$Tp <- input$Tp
        values$m <- input$m
        values$N <- input$N
        values$rho0 <- input$rho0
        values$r <- input$r
        values$type <- input$type
        values$effsiz <- input$effsiz
        values$accept_pwr<- input$accept_pwr
        values$c <- input$c
        values$p <- input$p
        values$pprim <- input$pprim
        values$g <- input$g
        values$gprim <- input$gprim
    })
    output$plotheader1a <- eventReactive(input$update, {
        header1a()
    })
    output$plotheader1b <- eventReactive(input$update, {
        header1b()
    })
    header1a <- renderPrint({
        tags$h3("The iterative removal of cells with low-information content")
    })
    header1b <- renderPrint({
        tags$h4("Updating the information content of remaining cells, and iterating")
    })
    output$plotheader2a <- eventReactive(input$update, {
        header2a()
    })
    output$plotheader2b <- eventReactive(input$update, {
        header2b()
    })
    header2a <- renderPrint({
        tags$h3("Variance of treatment effect estimator")
    })
    header2b <- renderPrint({
        tags$h4("Number of iterations")
    })
    output$plotheader3a <- eventReactive(input$update, {
        header3a()
    })
    output$plotheader3b <- eventReactive(input$update, {
        header3b()
    })
    header3a <- renderPrint({
        tags$h3("Precision loss compared to complete design")
    })
    header3b <- renderPrint({
        tags$h4("")
    })
    output$plotheader4a <- eventReactive(input$update, {
        header4a()
    })
    output$plotheader4b <- eventReactive(input$update, {
        header4b()
    })
    header4a <- renderPrint({
        tags$h3(paste0("Power to detect effect size of ", input$effsiz))
    })
    header4b <- renderPrint({
        tags$h4("Number of iterations")
    })
    output$plotheader5a <- eventReactive(input$update, {
        header5a()
    })
    output$plotheader5b <- eventReactive(input$update, {
        header5b()
    })
    header5a <- renderPrint({
        tags$h3("Trial cost by each design")
    })
    header5b <- renderPrint({
        tags$h4("")
    })
    output$plotheader6a <- eventReactive(input$update, {
        header6a()
    })
    output$plotheader6b <- eventReactive(input$update, {
        header6b()
    })
    header6a <- renderPrint({
        tags$h3("Relative cost efficiency = Cost efficiency of the progressively reduced designs / Cost efficiency of the complete design")
    })
    header6b <- renderPrint({
        tags$h4("Cost efficiency = Precision of the treatment effect / Total study cost")
    })
    output$plotheader7a <- eventReactive(input$update, {
        header7a()
    })
    output$plotheader7b <- eventReactive(input$update, {
        header7b()
    })
    header7a <- renderPrint({
        tags$h3("")
    })
    header7b <- renderPrint({
        tags$h4("")
    })
    
   
    ##option 2: fast    
    output$varREMplot <- renderPlotly({
        
        Xdlist <- list()  #design matrix
        Xdlist[[1]] <- SWdesmat(values$Tp)
        Inipow<- pow(CRTVarGeneralAdj(Xdlist[[1]],values$m,values$rho0,values$r,values$type)/values$N,values$effsiz,siglevel=0.05)*100
        
        tryCatch({
            FigRes <- FigGenDf(values$Tp, values$N, values$m, values$rho0, values$r, values$type, values$c, values$p, values$pprim, values$g, values$gprim, values$effsiz,values$accept_pwr)
          
            color_palette <- colorRampPalette(brewer.pal(8, "YlOrRd"))(length(unique(FigRes[[2]]$value)) - 2)

            p1 <- plot_ly(FigRes[[2]], x = ~as.numeric(Period), xgap = 4, y = ~as.factor(Sequence), ygap = 4, frame = ~iter,
                          z = ~value, type = 'heatmap', colors = color_palette,
                          hoverinfo = "text", showscale = FALSE, hoverlabel = list(bordercolor = NULL, font = list(size = 14)),
                          hovertext = ~ifelse(value < 100,
                                              paste("Xdvalue:", Xdvalue, "<br>InfCont:", value),
                                              paste("Xdvalue:", Xdvalue, "<br>InfCont:", 'IC cannot be calculated'))) %>%
                layout(plot_bgcolor = "gray", showlegend = FALSE,
                       xaxis = list(title = "Period", titlefont = list(size = 18), showline = TRUE,
                                    tickmode = "auto", tickfont = list(size = 16), nticks = 6, ticks = "inside",
                                    mirror = TRUE, showgrid = FALSE),
                       yaxis = list(title = "Sequence", titlefont = list(size = 18), tickfont = list(size = 16), autorange = "reversed",
                                    mirror = TRUE, showline = TRUE, showgrid = FALSE))
            
            p1
        }, error = function(e) {
            # If an error occurs, display a text message instead of the plot
            message <- paste0("Error: Unable to generate the plot.\n Please adjust the effect size or reduce the minimal acceptable power limit\n
            The initial power for effect size of", " ",values$effsiz," ","is"," ",format(round(Inipow,2),2),"%")
            plot_ly() %>%
                layout(title = message, font = list(size = 10, color = "Red"), margin = list(t = 50))
        })
    })
    output$Varsplot <- renderPlotly({
        FigRes <- FigGenDf(values$Tp, values$N, values$m, values$rho0, values$r, values$type, values$c, values$p, values$pprim, values$g, values$gprim, values$effsiz, values$accept_pwr)
        p2 <- generatePlotly(FigRes, ~iter, ~variance, "Iteration", "Variance",
                             ~paste("Iteration:", iter, "<br>Variance:", round(variance, 3)))
        p2
    })
    
    output$Prelossplot <- renderPlotly({
        FigRes <- FigGenDf(values$Tp, values$N, values$m, values$rho0, values$r, values$type, values$c, values$p, values$pprim, values$g, values$gprim, values$effsiz, values$accept_pwr)
        p3 <- generatePlotly(FigRes, ~iter, ~Preloss, "Iteration", "Precision loss (%)",
                             ~paste("Iteration:", iter, "<br>Preloss:", format(round(Preloss, 2), 2), "%"))
        p3
    })
    output$Powplot <- renderPlotly({
        FigRes <- FigGenDf(values$Tp, values$N, values$m, values$rho0, values$r, values$type, values$c, values$p, values$pprim, values$g, values$gprim, values$effsiz, values$accept_pwr)
        p4 <- generatePlotly(FigRes, ~iter, ~power, "Iteration", "Power (%)",
                             ~paste("Iteration:", iter, "<br>Power:", format(round(power, 2), 2), "%"))
        p4
    })
    
    output$Cstplot <- renderPlotly({
        FigRes <- FigGenDf(values$Tp, values$N, values$m, values$rho0, values$r, values$type, values$c, values$p, values$pprim, values$g, values$gprim, values$effsiz, values$accept_pwr)
        p5 <- generatePlotly(FigRes, ~iter, ~cost, "Iteration", "Cost ($)",
                             ~paste("Iteration:", iter, "<br>Cost:", "$", format(cost, big.mark = ",", scientific = FALSE)))
        p5
    })
    output$CstVarplot<-renderPlotly({

        FigRes <-  FigGenDf(values$Tp, values$N, values$m, values$rho0, values$r, values$type, values$c, values$p, values$pprim, values$g, values$gprim,values$effsiz,values$accept_pwr)
        p6 <- generatePlotly(FigRes, ~variance, ~cost, "Variance", "Cost ($)",
                             ~paste("Iteration:", iter, "<br>Cost:","$",format(cost,big.mark=",",scientific=FALSE) ,
                               "<br>Variance:",  round(variance,4),
                               "<br>Power:",  format(round(power,2),2),"%"))
        p6
    })
    #Relative cost efficiency plot
    output$RCEplot<-renderPlotly({

        FigRes <-  FigGenDf(values$Tp, values$N, values$m, values$rho0, values$r, values$type, values$c, values$p, values$pprim, values$g, values$gprim,values$effsiz,values$accept_pwr)
        p7 <- generatePlotly(FigRes, ~iter, ~RCE, "Iteration", "Relative Cost Efficiency (RCE)",
                             ~paste("Iteration:", iter, "<br>Cost:","$",format(cost,big.mark=",",scientific=FALSE)
                                         ,"<br>RCE", round(RCE,4)))
        p7

    })
}
    # Run the application
    shinyApp(ui = ui, server = server)
    
    
    
    