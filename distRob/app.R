#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(dashboardthemes)

library(abind)
library(viridis)
library(ggplot2)
library(english)

library(glmnet)
source("../R/SDForest.r")

# Define UI for application that draws a histogram
ui <- shinyUI(
  dashboardPage(
    #skin = 'blue',
    
    dashboardHeader(title = "hello"),
    
    dashboardSidebar(
      actionButton("show_dim", label = "Dimensions"),
      conditionalPanel(
        condition = "output.show_dim",
        numericInput("n_train", label = "Number of Training Samples", value = 500, min = 1, max = 10000),
        numericInput("n_test", label = "Number of Test Samples", value = 500, min = 1, max = 10000),            
        numericInput("p", label = "Number of Features", value = 50, min = 2, max = 10000),
        numericInput("q", label = "Dimension of hidden Confounding", value = 10, min = 1, max = 10000),
        numericInput("r", label = "Dimension of observed Anchor", value = 1, min = 1, max = 10000), 
        numericInput("n_caus", label = "Number of Causal Features", value = 1, min = 1, max = 50),
        ),

      actionButton("data_distribution", label = "Data Distribution"),
      conditionalPanel(
        condition = "output.data_distribution",
        numericInput("var_X", label = "Variance of Features", value = 1, min = 0, max = 10000),
        numericInput("var_H", label = "Variance of hidden Confounding", value = 2, min = 0, max = 10000),
        numericInput("var_A", label = "Variance of Anchor", value = 1, min = 0, max = 10000),
        numericInput("var_Y", label = "Variance of Response", value = 0.01, min = 0, max = 10000)
      ),
      
      actionButton("effect_strength", label = "Effect Strength"),
      conditionalPanel(
        condition = "output.effect_strength",
        numericInput("beta_strength", label = "Strength of Causal Features", value = 1, min = 0, max = 10000),
        numericInput("delta_strength", label = "Strength of Confounding", value = 1, min = 0, max = 10000),
        numericInput("gamma_strength", label = "Strength of Confounding Features", value = 1, min = 0, max = 10000),
        numericInput("alpha_1_strength", label = "Strength of Anchor Features", value = 1, min = 0, max = 10000),
        numericInput("alpha_2_strength", label = "Strength of Anchor", value = 1, min = 0, max = 10000),
        numericInput("alpha_3_strength", label = "Strength of Confounding Anchor", value = 0, min = 0, max = 10000)
      ),
      
      actionButton("shift", label = "Distribution Shift"),
      conditionalPanel(
        condition = "output.shift",
        numericInput("A_shift", label = "Shift of Anchor", value = 1, min = 0, max = 10000), 
        numericInput("A_gamma", label = "gamma Parameter", value = 1, min = 0, max = 10000)
      )
    ),
    
    dashboardBody(
      shinyDashboardThemes(
        theme = "grey_light"
      ),
      
      tags$head(tags$style(HTML('.content-wrapper { overflow: auto !important;}'))),
      
      
      tabBox(#title = 'Binding ELISA',
        tabPanel(title = 'Positive Control',
                 fluidRow(box(plotOutput("controlplot", height = 250), width = 12))
        ),
        
        width = 12
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output, session) {
  #reactive variables
  rv <- reactiveValues(show_dim = FALSE)
  
  observe(updateSliderInput(session, "n_caus", max = input$p))
  
  observeEvent(rv$plates, {
    updateCheckboxGroupInput(session, inputId = "chosen_antigens",
                             choices = dimnames(rv$plates)[[1]]
    )
  })

  output$show_dim <- reactive({
    return(input$show_dim %% 2)
  })
  outputOptions(output, 'show_dim', suspendWhenHidden=FALSE)

  output$data_distribution <- reactive({
    return(input$data_distribution %% 2)
  })
  outputOptions(output, 'data_distribution', suspendWhenHidden=FALSE)

  output$effect_strength <- reactive({
    return(input$effect_strength %% 2)
  })
  outputOptions(output, 'effect_strength', suspendWhenHidden=FALSE)

  output$shift <- reactive({
    return(input$shift %% 2)
  })
  outputOptions(output, 'shift', suspendWhenHidden=FALSE)


  
  output$controlplot <- renderPlot({
    n_train <- input$n_train
    n_test <- input$n_test
    
    r <- input$r
    q <- input$q
    p <- input$p

    n_caus <- min(input$n_caus, input$p)

    var_X <- input$var_X
    var_H <- input$var_H
    var_A <- input$var_A
    var_Y <- input$var_Y

    beta_strength <- input$beta_strength
    delta_strength <- input$delta_strength
    gamma_strength <- input$gamma_strength
    alpha_1_strength <- input$alpha_1_strength
    alpha_2_strength <- input$alpha_2_strength
    alpha_3_strength <- input$alpha_3_strength

    A_shift <- input$A_shift

    beta <- matrix(0, nrow = p)
    beta[1:n_caus, ] <- rnorm(n_caus, 0, beta_strength)

    delta <- matrix(rnorm(q, 0, delta_strength), nrow = q)
    gamma <- matrix(rnorm(q * p, 0, gamma_strength), nrow = q)
    alpha_1 <- matrix(rnorm(r * p, 0, alpha_1_strength), nrow = r)
    alpha_2 <- matrix(rnorm(r, 0, alpha_2_strength), nrow = r)
    alpha_3 <- matrix(rnorm(r * q, 0, alpha_3_strength), nrow = r)

    A <- matrix(rnorm(n_train * r, 0, var_A), nrow = n_train)
    H <- A %*% alpha_3 + matrix(rnorm(n_train * q, 0, var_H), nrow = n_train)
    X <- A %*% alpha_1 + H %*% gamma + matrix(rnorm(n_train * p, 0, var_X), nrow = n_train)
    Y <- X %*% beta + H %*% delta + A %*% alpha_2 + rnorm(n_train, 0, var_Y)

    mean(abs(X %*% beta))
    mean(abs(H %*% delta))
    mean(abs(A %*% alpha_2))
    mean(abs(rnorm(n_train, 0, var_Y)))

    A_test <- matrix(rnorm(n_test * r, 0, var_A), nrow = n_test) + A_shift
    H_test <- A_test %*% alpha_3 + matrix(rnorm(n_test * q, 0, var_H), nrow = n_test)
    X_test <- A_test %*% alpha_1 + H_test %*% gamma + matrix(rnorm(n_test * p, 0, var_X), nrow = n_test)
    Y_test <- X_test %*% beta + H_test %*% delta + A_test %*% alpha_2 + rnorm(n_test, 0, var_Y)

    W <- get_W(A, gamma = input$A_gamma)
    Q <- get_Q(X, type = 'trim')



    fit_lin_cv <- cv.glmnet(x = X, y = Y, alpha = 1, nfolds = 10)
    fit_lin <- glmnet(x = X, y = Y, alpha = 1, lambda = fit_lin_cv$lambda.1se)

    fit_linA_cv <- cv.glmnet(x = W %*% X, y = W %*% Y, alpha = 1, nfolds = 10)
    fit_linA <- glmnet(x = W %*% X, y = W %*% Y, alpha = 1, lambda = fit_linA_cv$lambda.1se)


    fit_linQ_cv <- cv.glmnet(x = Q %*% X, y = Q %*% Y, alpha = 1, nfolds = 10)

    fit_linQ <- glmnet(x = Q %*% X, y = Q %*% Y, alpha = 1, lambda = fit_linQ_cv$lambda.1se)

    #fit_linA <- lm.fit(x = W %*% X, y = W %*% Y)$coefficients
    #fit_lin <- lm.fit(x = X, y = Y)$coefficients
    #fit_linQ <- lm.fit(x = Q %*% X, y = Q %*% Y)$coefficients

    coefA <- coef(fit_linA)[-1, 1]
    coef <- coef(fit_lin)[-1, 1]
    coefQ <- coef(fit_linQ)[-1, 1]


    pred_linA <- predict(fit_linA, newx = X_test)[, 1]
    pred_lin <- predict(fit_lin, newx = X_test)[, 1]
    pred_linQ <- predict(fit_linQ, newx = X_test)[, 1]


    coef_col <- rep('black', p)
    coef_col[1:n_caus] <- '#0ea10e'

    ceof_df <- data.frame(beta = abs(beta), fit_linA = abs(coefA), 
      fit_lin = abs(coef), fit_linQ = abs(coefQ))
    #plot(ceof_df, col = coef_col, pch = 20, cex = 1)


    par(mfrow = c(1, 2))
    plot(X_test[, 1], Y_test)
    points(X_test[, 1], pred_lin, col = 'blue', pch = 20, cex = 1)
    points(X_test[, 1], pred_linA, col = 'red', pch = 20, cex = 0.8)
    points(X_test[, 1], pred_linQ, col = '#be24b7', pch = 20, cex = 0.5)
    points(X_test[, 1], X_test %*% beta, col = 'green', pch = 20, cex = 0.3)
  })
  

  rv$barplot <- reactive({    
    
    n <- length(input$chosen_antigens)
    if(input$plate_type == 'Antigens'){
      dilutions <- c('undil', paste('1:', input$dilution, sep = ''), paste('1:', input$dilution**2, sep = ''), paste('1:', input$dilution**3, sep = ''))
    }else{
      dilutions <- c('1:100', '1:400', '1:1600', '1:6400')
    }    
    
    if(n != 0){
      barplot_data <- data.frame(c(rep(dilutions[1], n), rep(dilutions[2], n), rep(dilutions[3], n), rep(dilutions[4], n)), 
                                 rep(input$chosen_antigens, 4), as.vector(rv$plates[input$chosen_antigens, , input$antibody]))
      names(barplot_data) <- c('concentration', 'antigen', 'absorbance')
      barplot_data$concentration <- factor(barplot_data$concentration, levels = dilutions)
      barplot_data$antigen <- factor(barplot_data$antigen, levels = input$chosen_antigens)
      
      
      ggplot(barplot_data, aes(fill=antigen, y=absorbance, x=concentration)) + 
        geom_bar(position="dodge", stat="identity", width = 0.35 * n**0.5) +
        scale_fill_viridis(discrete = T) +
        theme(panel.grid = element_line(size = (0.2), colour="grey"), plot.title = element_text(size=20)) + 
        ggtitle(input$antibody_name,) + 
        labs(y = expression("OD"["450 nm"]), x = 'dilution factor') + 
        ylim(0, max(c(1.3, barplot_data$absorbance)) + 0.2) +
        labs(fill = input$plate_type)
    }
    
  })
  
})

# Run the application 
shinyApp(ui = ui, server = server)
