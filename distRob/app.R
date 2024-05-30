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

library(igraph)
library(ggplot2)
library(plotly)
library(tidyr)
library(glmnet)
source("SDForest.r")

# Define UI for application that draws a histogram
ui <- shinyUI(
  dashboardPage(
    #skin = 'blue',
    
    dashboardHeader(title = ''),
    
    dashboardSidebar(
      actionButton("show_dim", label = "Dimensions"),
      conditionalPanel(
        condition = "output.show_dim",
        sliderInput("n_train", label = "Number of Training Samples", value = 500, min = 1, max = 10000),
        sliderInput("n_test", label = "Number of Test Samples", value = 500, min = 1, max = 1000),            
        sliderInput("p", label = "Number of Features", value = 500, min = 2, max = 1000),
        sliderInput("q", label = "Dimension of hidden Confounding", value = 2, min = 1, max = 1000),
        sliderInput("r", label = "Dimension of observed Anchor", value = 4, min = 1, max = 1000), 
        sliderInput("n_caus", label = "Number of Causal Features", value = 1, min = 1, max = 50),
        ),

      #actionButton("data_distribution", label = "Data Distribution"),
      #conditionalPanel(
      #  condition = "output.data_distribution",
      #  numericInput("var_X", label = "Variance of Features", value = 1, min = 0, max = 10000),
      #  numericInput("var_H", label = "Variance of hidden Confounding", value = 2, min = 0, max = 10000),
      #  numericInput("var_A", label = "Variance of Anchor", value = 1, min = 0, max = 10000),
      #  numericInput("var_Y", label = "Variance of Response", value = 0.01, min = 0, max = 10000)
      #),
      
      actionButton("effect_strength", label = "Effect Strength"),
      conditionalPanel(
        condition = "output.effect_strength",
        sliderInput("var_Y", label = "Variance of Response", value = 0.1, min = 0, max = 10, step = 0.01),
        sliderInput("beta_strength", label = "beta", value = 1, min = 0, max = 10, step = 0.001),
        sliderInput("delta_strength", label = "delta variance", value = 1, min = 0, max = 10, step = 0.01),
        sliderInput("gamma_strength", label = "gamma variance", value = 1, min = 0, max = 10, step = 0.01),
        sliderInput("alpha_1_strength", label = "alpha 1 variance", value = 1, min = 0, max = 10, step = 0.01),
        sliderInput("alpha_2_strength", label = "alpha 2 variance", value = 1, min = 0, max = 10, step = 0.01),
        sliderInput("alpha_3_strength", label = "alpha 3 variance", value = 0, min = 0, max = 10, step = 0.01)
      ),
      
      actionButton("shift", label = "Distribution Shift"),
      conditionalPanel(
        condition = "output.shift",
        sliderInput("A_shift", label = "Shift of Anchor", value = 5, min = 0, max = 100, step = 0.01), 
        sliderInput("A_gamma", label = "gamma Parameter", value = 10, min = 0, max = 100, step = 0.01)
      ),
      textOutput('sim')
    ),
    
    
    dashboardBody(
      shinyDashboardThemes(
        theme = "grey_light"
      ),
      
      tags$head(tags$style(HTML('.content-wrapper { overflow: auto !important;}'))),
      
      
      tabBox(#title = 'Binding ELISA',
        tabPanel(title = 'Linear Regression',
          fluidRow(box(actionButton("run", "Run simulation"), width = 12)),
          fluidRow(box(plotOutput("graph"), width = 6), 
          box(plotlyOutput("perf"), width = 6)),
          fluidRow(box(plotlyOutput("sim_res"), width = 12))
          #fluidRow(box(plotOutput("perf"), width = 12))
        ),
        
        width = 12
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output, session) {
  #reactive variables
  rv <- reactiveValues(perf = ggplot(), sim_res = ggplot(), beta = NULL, delta = NULL, gamma = NULL, alpha_1 = NULL, 
  alpha_2 = NULL, alpha_3 = NULL, A = NULL, H = NULL, X = NULL, Y = NULL, A_test = NULL, 
  H_test = NULL, X_test = NULL, Y_test = NULL)
  
  
  observe(updateSliderInput(session, "n_caus", max = input$p))
  

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
  outputOptions(output, 'shift', suspendWhenHidden = FALSE)

  output$sim <- renderText({
    print(input$run)
    n_train <- input$n_train
    n_test <- input$n_test
    
    r <- input$r
    q <- input$q
    p <- input$p

    n_caus <- min(input$n_caus, input$p)

    var_X <- 1
    var_H <- 1
    var_A <- 1
    var_Y <- input$var_Y

    beta_strength <- input$beta_strength
    delta_strength <- input$delta_strength
    gamma_strength <- input$gamma_strength
    alpha_1_strength <- input$alpha_1_strength
    alpha_2_strength <- input$alpha_2_strength
    alpha_3_strength <- input$alpha_3_strength

    A_shift <- input$A_shift

    beta <- matrix(0, nrow = p)
    #beta[1:n_caus, ] <- rnorm(n_caus, 0, sqrt(beta_strength))
    beta[1:n_caus] <- beta_strength

    delta <- matrix(rnorm(q, 0, sqrt(delta_strength)), nrow = q)
    gamma <- matrix(rnorm(q * p, 0, sqrt(gamma_strength)), nrow = q)
    alpha_2 <- matrix(rnorm(r * p, 0, sqrt(alpha_2_strength)), nrow = r)
    alpha_3 <- matrix(rnorm(r, 0, sqrt(alpha_3_strength)), nrow = r)
    alpha_1 <- matrix(rnorm(r * q, 0, sqrt(alpha_1_strength)), nrow = r)

    A <- matrix(rnorm(n_train * r, 0, sqrt(var_A)), nrow = n_train)
    H <- A %*% alpha_1 + matrix(rnorm(n_train * q, 0, sqrt(var_H)), nrow = n_train)
    X <- A %*% alpha_2 + H %*% gamma + matrix(rnorm(n_train * p, 0, sqrt(var_X)), nrow = n_train)
    Y <- X %*% beta + H %*% delta + A %*% alpha_3 + rnorm(n_train, 0, sqrt(var_Y))

    A_test <- matrix(rnorm(n_test * r, 0, sqrt(var_A)), nrow = n_test) + A_shift
    H_test <- A_test %*% alpha_1 + matrix(rnorm(n_test * q, 0, sqrt(var_H)), nrow = n_test)
    X_test <- A_test %*% alpha_2 + H_test %*% gamma + matrix(rnorm(n_test * p, 0, sqrt(var_X)), nrow = n_test)
    Y_test <- X_test %*% beta + H_test %*% delta + A_test %*% alpha_3 + rnorm(n_test, 0, sqrt(var_Y))


    rv$beta <- beta
    rv$delta <- delta
    rv$gamma <- gamma
    rv$alpha_1 <- alpha_1
    rv$alpha_2 <- alpha_2
    rv$alpha_3 <- alpha_3
    rv$A <- A
    rv$H <- H
    rv$X <- X
    rv$Y <- Y
    rv$A_test <- A_test
    rv$H_test <- H_test
    rv$X_test <- X_test
    rv$Y_test <- Y_test
    return(as.character(r))
  })

  observeEvent(input$run, {
    A <- rv$A
    H <- rv$H
    X <- rv$X
    Y <- rv$Y
    A_test <- rv$A_test
    H_test <- rv$H_test
    X_test <- rv$X_test
    Y_test <- rv$Y_test

    p <- input$p

    W <- get_W(A, gamma = input$A_gamma)
    Q <- get_Q(X, type = 'trim')

    fit_lin_cv <- cv.glmnet(x = X, y = Y, alpha = 1, nfolds = 10)
    fit_lin <- glmnet(x = X, y = Y, alpha = 1, lambda = fit_lin_cv$lambda.1se)

    fit_linA_cv <- cv.glmnet(x = W %*% X, y = W %*% Y, alpha = 1, nfolds = 10)
    fit_linA <- glmnet(x = W %*% X, y = W %*% Y, alpha = 1, lambda = fit_linA_cv$lambda.1se)

    fit_linQ_cv <- cv.glmnet(x = Q %*% X, y = Q %*% Y, alpha = 1, nfolds = 10)
    fit_linQ <- glmnet(x = Q %*% X, y = Q %*% Y, alpha = 1, lambda = fit_linQ_cv$lambda.1se)

    coefA <- coef(fit_linA)[-1, 1]
    coef <- coef(fit_lin)[-1, 1]
    coefQ <- coef(fit_linQ)[-1, 1]

    pred_linA <- predict(fit_linA, newx = X_test)[, 1]
    pred_lin <- predict(fit_lin, newx = X_test)[, 1]
    pred_linQ <- predict(fit_linQ, newx = X_test)[, 1]

    coef_col <- rep('black', p)
    coef_col[1:input$n_caus] <- '#0ea10e'

    ceof_df <- data.frame(beta = abs(rv$beta), fit_linA = abs(coefA), 
      fit_lin = abs(coef), fit_linQ = abs(coefQ))

    df_res <- data.frame(X_1 = X_test[, 1], Y = Y_test, 
      Anchor = pred_linA, lasso = pred_lin, trim = pred_linQ, 
      True = X_test %*% rv$beta)
    
    df_res <- gather(df_res, key = 'Method', value = 'Y', -X_1)

    col_method <- c(Y = 'black', Anchor = 'red', lasso = 'blue', True = 'green', trim = 'purple')
    df_res$size <- 1
    df_res$size[df_res$Method == 'Anchor'] <- 1.5

    rv$sim_res <- ggplot(df_res, aes(x = X_1, y = Y, color = Method, pch = Method)) + 
      geom_point(size = df_res$size) + 
      theme_bw() + 
      xlab('X1') +
      ylab('Test Y') +
      theme(legend.position = 'bottom', legend.title = element_blank()) +
      scale_color_manual(values = col_method)


    perf <- data.frame(lasso = (Y_test - pred_lin)**2, Anchor = (Y_test - pred_linA)**2, 
      True = (Y_test - X_test %*% rv$beta)**2, trim = (Y_test - pred_linQ)**2)
    
    perf <- gather(perf, key = 'Method', value = 'MSE')
    rv$perf <- ggplot(perf, aes(x = Method, y = MSE, fill = Method)) + 
      geom_boxplot() + 
      theme_bw() + 
      ylab('Test MSE') + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      #theme(legend.position = 'none') + 
      scale_fill_manual(values = col_method)
      
  })

  
  output$sim_res <- renderPlotly({
    ggplotly(rv$sim_res)
  })

  output$perf <- renderPlotly({
    ggplotly(rv$perf)
  })

  output$graph <- renderPlot({
    strengths <- c(mean(abs(rv$X %*% rv$beta)), mean(abs(rv$H %*% rv$delta)), 
      mean(abs(rv$H %*% rv$gamma)), mean(abs(rv$A %*% rv$alpha_2)), 
      mean(abs(rv$A %*% rv$alpha_1)), mean(abs(rv$A %*% rv$alpha_3)))

    effect_present <- as.logical(strengths)
    label_col <- c('black', 'black', 'black', 'black', 'black', 'black')
    label_col[!effect_present] <- 'white'

    strengths <- strengths - min(strengths)
    strengths <- strengths / max(strengths) * 3 + 1

    dims <- c(input$p, input$q, 1, input$r)
    dims <- dims - min(dims)
    dims <- dims / max(dims) * 40 + 20 

    nodes <- data.frame(name = c('X', 'H', 'Y', 'A'), 
                        size = dims, 
                        color = c('#2ad4f1', 'white', '#32f77d', '#df6c6c'),
                        label.color = c('black', 'black', 'black', 'black'), 
                        label.cex = 1.5)

    edges <- data.frame(from = c('X', 'H', 'H', 'A', 'A', 'A'), 
                        to = c('Y', 'Y', 'X', 'X', 'H', 'Y'), 
                        width = strengths, 
                        label = c('beta\n', 'delta\n', 'gamma\n', 'alpha 2\n', 'alpha 1\n', 'alpha 3\n'),
                        color = effect_present,
                        arrow.size = 1.5, 
                        label.color = label_col,
                        label.cex = 1.5, 
                        label.dist = 1.5)

    g <- graph_from_data_frame(edges, directed = T, vertices = nodes)

    set.seed(6)
    coords <- layout.reingold.tilford(barabasi.game(4))
    plot(g, layout = coords, margin = 0)
  })
  
})

# Run the application 
shinyApp(ui = ui, server = server)
