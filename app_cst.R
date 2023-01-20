#setwd("G:\\Shared drives\\Ehsan PhD work\\Codes\\Git\\Iterative-removals-IC-cost\\")
source("CRTVarAdj_func.R", local=TRUE)
source("IterRemCst_func.R", local=TRUE)
#source("IterRemoval_func.R", local=TRUE)

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

ui <- fluidPage(
    titlePanel(h1("Information content of progressively reduced stepped wedge designs",h2(""),h3(""))),
    sidebarLayout(
        
        sidebarPanel(
            h4("Trial configuration"),
            
            sliderInput(inputId = "Tp", label = "Number of periods:",
                        min = 5, max = 20, value = 5,  step = 1,ticks = FALSE),
            numericInput("m",
                         "Number of subjects in each cluster-period, m:",
                         min = 1,
                         max=1000,
                         step = 1,
                         value = 100),
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
            
            sliderInput("effsize", "Effect size:",
                        min = 0.05, max = 1.0,
                        value = 0.2, step = 0.05,ticks = FALSE),
            
            hr(),
            
            h4("Cost components"),
            
            # Input: Interval, cost per cluster
            numericInput("c", "Cost per cluster",
                         min = 100, max = 100000,
                         value = 2500, step = 100),
            
            # Input: Interval, cost per subject receiving intervention 
            numericInput("s", "Cost per subject receiving intervention ",
                         min = 0, max = 5000,
                         value = 100, step = 50),
            
            # Input: Interval, cost per subject receiving control
            numericInput("sprim", "Cost per subject receiving control",
                         min = 0, max = 5000,
                         value = 100, step = 50),
            
            # Input: Interval, Cost of restarting an intervention condition
            numericInput("R", "Restart cost (intervention)",
                         min = 0, max = 50000,
                         value = 250, step = 50),
            bsTooltip("R","Cost of restarting an intervention condition","right"),
            # Input: Interval, Cost of restarting a control condition
            numericInput("Rprim", "Restart cost (control)",
                         min = 0, max = 50000,
                         value = 250, step = 50),
            bsTooltip("Rprim","Cost of restarting a control condition","right"),
            # clicks the button
            actionButton("update", "Update"),
        ),
        mainPanel(
            tabsetPanel(
                tabPanel("The iterative removal of cluster-period cells",
                         uiOutput("plotheader1a"), uiOutput("plotheader1b"),
                         plotlyOutput("varREMplot"),
                         textOutput("ICremtext")
                ),
                tabPanel("Variance",
                         uiOutput("plotheader2a"), uiOutput("plotheader2b"),
                         plotlyOutput("Varsplot"),
                         textOutput("ICvartext")
                ),
                tabPanel("Efficiency loss",
                         uiOutput("plotheader3a"), uiOutput("plotheader3b"),
                         plotlyOutput("Efflossplot"),
                         textOutput("ICRvartext")
                ),
                tabPanel("Power",
                         uiOutput("plotheader4a"), uiOutput("plotheader4b"),
                         plotlyOutput("Powplot"),
                         textOutput("ICpowtext")
                ),
                tabPanel("Trial cost",
                         uiOutput("plotheader5a"), uiOutput("plotheader5b"),
                         plotlyOutput("Cstplot"),
                         textOutput("ICcsttext")
                ),
                tabPanel("Design choice",
                           uiOutput("plotheader6a"), uiOutput("plotheader6b"),
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
    output$ICremtext <- renderText({
        "..."
    })
    output$ICvartext <- renderText({
        "..."
    })
    output$ICRvartext <- renderText({
        "..."
    })
    output$ICPowtext <- renderText({
        "..."
    })
    output$ICCsttext <- renderText({
        "..."
    })
    output$ICstVartext <- renderText({
        "..."
    })
    
    values <- reactiveValues(
        Tp = 5,
        m = 100,
        rho0 = 0.05,
        r=0.95,
        type=1,
        effsize=0.2,
        c=2500,
        s=50,
        sprim=50,
        R=250,
        Rprim=250
    )
    
    observeEvent(input$update, {
        values$Tp <- input$Tp
        values$m <- input$m
        values$rho0 <- input$rho0
        values$r <- input$r
        values$type <- input$type
        values$effsize <- input$effsize
        values$c <- input$c
        values$s <- input$s
        values$sprim <- input$sprim
        values$R <- input$R
        values$Rprim <- input$Rprim
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
        tags$h3("Efficiency loss compared to complete design")
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
        tags$h3(paste0("Power to detect effect size of ", input$effsize))
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
        tags$h3("")
    })
    header6b <- renderPrint({
        tags$h4("")
    })
    
    output$varREMplot<-renderPlotly({
        
        Tp=values$Tp
        K=values$Tp-1
        m=values$m
        rho0=values$rho0
        r=values$r
        type=values$type
        c=values$c
        s=values$s
        sprim=values$sprim
        R=values$R
        Rprim=values$R
        
        #Put parameters here
        IterRes<- IterRemove(Tp,m,rho0,r,type,c,s,sprim,R,Rprim)
        melted_varmatexcl<- melt(IterRes[[1]])
        melted_desmatexcl<- melt(IterRes[[2]])
        
        names(melted_desmatexcl)[names(melted_desmatexcl)=="value"] <- "Xdvalue"
        melted_varmatexcl_t<- jointdataset <- merge(melted_varmatexcl, melted_desmatexcl, by = c('Var1','Var2','L1'))
        melted_varmatexcl_t$value<-round(melted_varmatexcl_t$value, 4)
        
        color_palette <-colorRampPalette(brewer.pal(8, "YlOrRd"))(length(unique(melted_varmatexcl_t$value))-2)
        
        
        names(melted_varmatexcl_t)[names(melted_varmatexcl_t)=="Var1"] <- "Sequence"
        names(melted_varmatexcl_t)[names(melted_varmatexcl_t)=="Var2"] <- "Period"
        names(melted_varmatexcl_t)[names(melted_varmatexcl_t)=="L1"] <- "iter"
        
        Xdes <- SWdesmat(Tp)
        varmatall<- IterRes[[3]]
        cvec<- IterRes[[4]]
        
        Tp <- ncol(Xdes)
        K  <- nrow(Xdes)
        
        #power
        iter=1:length(varmatall)
        df=as.data.frame(varmatall)
        
        pow <- function(vars, effsize, siglevel=0.05){
            z <- qnorm(siglevel/2)
            pow <- pnorm(z + sqrt(1/vars)*effsize)
            return(pow)
        }
        
        powdf <- function(df, effsize, siglevel=0.05){
            powvals <- apply(df, MARGIN=2, pow, effsize, siglevel)
            powdf <- data.frame(iter, df$varmatall,powvals*100)
            colnames(powdf) <- c("iter","variance","power")
            return(powdf)
        }
        res <- powdf(df,values$effsize)
        res$r <- r
        
        
        res <- cbind(res,res$variance[1]/res$variance,(1-(res$variance[1]/res$variance))*100,cvec)
        colnames(res) <- c("iter","variance","power","r","Rvariance","Effloss","cost")
        
        
        melted_varmatexcl_t <- merge(res, melted_varmatexcl_t, by = "iter", all = TRUE)
        
        ##option 2: fast
        p1<-plot_ly(melted_varmatexcl_t, x = ~as.numeric(Period), xgap = 4,y = ~as.factor(Sequence),ygap =4, frame = ~iter,
                    z=~value,type = 'heatmap', colors=color_palette,
                    hoverinfo="text",showscale = FALSE,hoverlabel=list(bordercolor=NULL, font=list(size=14)),
                    hovertext=~ifelse(value<100,
                                      paste("Xdvalue:",Xdvalue,"<br>InfCont:",value),
                                      paste("Xdvalue:",Xdvalue,"<br>InfCont:",'IC cannot be calculated')))%>%
            layout(plot_bgcolor="gray",showlegend = FALSE,
                   xaxis=list(title="Period", titlefont=list(size=18), showline=TRUE,
                              tickmode="auto", tickfont=list(size=16), nticks=6, ticks="inside",
                              mirror=TRUE,showgrid=FALSE),
                   yaxis=list(title="Sequence", titlefont=list(size=18), tickfont=list(size=16), autorange="reversed",
                              mirror=TRUE, showline=TRUE,showgrid=FALSE))
        
        output$Varsplot<-renderPlotly({
            
            p <- plot_ly(res, height=500, width=800, x=~iter, y=~variance, name="Variance", type="scatter",
                         mode="lines", hoverinfo="text", hoverlabel=list(bordercolor=NULL, font=list(size=16)),
                         text=~paste("Iteration:", iter, "<br>Variance:", round(variance, 3)),
                         line=list(color="#F8766D", width=4, dash="dash"))%>%
                layout(xaxis=list(title="Iteration", titlefont=list(size=18), showline=TRUE,
                                  tickmode="auto", tickfont=list(size=16), nticks=6, ticks="inside",
                                  mirror=TRUE, showgrid=FALSE),
                       yaxis=list(title="Variance", titlefont=list(size=18), tickfont=list(size=16),
                                  mirror=TRUE, showline=TRUE),
                       legend=list(orientation="h", xanchor="center", yanchor="bottom", x=0.5, y=-0.5, font=list(size=16)),
                       margin=list(l=100, r=40))
            print(p)
        })
        
        output$Efflossplot<-renderPlotly({
            
            p <- plot_ly(res, height=500, width=800, x=~iter, y=~Effloss, name="Effloss", type="scatter",
                         mode="lines", hoverinfo="text", hoverlabel=list(bordercolor=NULL, font=list(size=16)),
                         text=~paste("Iteration:", iter, "<br>Effloss:", format(round(Effloss,2),2),"%"),
                         line=list(color="#F8766D", width=4, dash="dash"))%>%
                layout(xaxis=list(title="Iteration", titlefont=list(size=18), showline=TRUE,
                                  tickmode="auto", tickfont=list(size=16), nticks=6, ticks="inside",
                                  mirror=TRUE, showgrid=FALSE),
                       yaxis=list(title="Efficiency loss (%)", titlefont=list(size=18), tickfont=list(size=16),
                                  mirror=TRUE, showline=TRUE),
                       legend=list(orientation="h", xanchor="center", yanchor="bottom", x=0.5, y=-0.5, font=list(size=16)),
                       margin=list(l=100, r=40))
            print(p)
        })
        
        output$Powplot<-renderPlotly({
            p <- plot_ly(res, height=500, width=800, x=~iter, y=~power, name="Power", type="scatter",
                         mode="lines", hoverinfo="text", hoverlabel=list(bordercolor=NULL, font=list(size=16)),
                         text=~paste("Iteration:", iter, "<br>Power:",  format(round(power,2),2),"%"),
                         line=list(color="#00BA38", width=4, dash="dash"))%>%
                layout(xaxis=list(title="Iteration", titlefont=list(size=18), showline=TRUE,
                                  tickmode="auto", tickfont=list(size=16), nticks=6, ticks="inside",
                                  mirror=TRUE, showgrid=FALSE),
                       yaxis=list(title="Power (%)", titlefont=list(size=18), tickfont=list(size=16),
                                  mirror=TRUE, showline=TRUE),
                       legend=list(orientation="h", xanchor="center", yanchor="bottom", x=0.5, y=-0.5, font=list(size=16)),
                       margin=list(l=100, r=60))
            print(p)
        })
        
        output$Cstplot<-renderPlotly({
            p <- plot_ly(res, height=500, width=800, x=~iter, y=~cost, name="cost", type="scatter",
                         mode="lines", hoverinfo="text", hoverlabel=list(bordercolor=NULL, font=list(size=16)),
                         text=~paste("Iteration:", iter, "<br>Cost:","$",format(cost,big.mark=",",scientific=FALSE)),
                         line=list(color="#00BA38", width=4, dash="dash"))%>%
                layout(xaxis=list(title="Iteration", titlefont=list(size=18), showline=TRUE,
                                  tickmode="auto", tickfont=list(size=16), nticks=6, ticks="inside",
                                  mirror=TRUE, showgrid=FALSE),
                       yaxis=list(title="Cost ($)", titlefont=list(size=18), tickfont=list(size=16),
                                  mirror=TRUE, showline=TRUE),
                       legend=list(orientation="h", xanchor="center", yanchor="bottom", x=0.5, y=-0.5, font=list(size=16)),
                       margin=list(l=100, r=60))
            print(p)
        })
        
        output$CstVarplot<-renderPlotly({
                        p <- plot_ly(res, height=500, width=800)%>%
                            add_markers(x=~variance, y=~cost,type="scatter",
                                        mode = 'markers', hoverinfo="text", hoverlabel=list(bordercolor=NULL, font=list(size=16)),
                                        text=~paste("Iteration:", iter, "<br>Cost:","$",format(cost,big.mark=",",scientific=FALSE) ,
                                                    "<br>Variance:",  round(variance,4)),showlegend = F)%>%
                layout(xaxis=list(title="Variance", titlefont=list(size=18), showline=TRUE,
                                  tickmode="auto", tickfont=list(size=16), nticks=6, ticks="inside",
                                  mirror=TRUE, showgrid=FALSE),
                       yaxis=list(title="Cost ($)", titlefont=list(size=18), tickfont=list(size=16),
                                  mirror=TRUE, showline=TRUE),
                       legend=list(orientation="h", xanchor="center", yanchor="bottom", x=0.5, y=-0.5, font=list(size=16)),
                       margin=list(l=100, r=60))
        })

        
        print(p1)
    })
    
}

# Run the application
shinyApp(ui = ui, server = server)

