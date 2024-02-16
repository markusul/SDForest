true_function <- function(beta, js){
    res <- list(beta = beta, js = js)
    class(res) <- 'true_function'
    return(res)
}

predict.true_function <- function(object, newdata){
    f_X <- apply(newdata, 1, function(x) f_four(x, object$beta, object$js))
    return(f_X)
}

plotDep <- function(object, n_examples = 19){
  ggdep <- ggplot2::ggplot() + ggplot2::theme_bw()
  preds <- object$preds
  x_seq <- object$x_seq
  
  sample_examples <- sample(1:dim(preds)[2], n_examples)
  for(i in sample_examples){
      pred_data <- data.frame(x = x_seq, y = preds[, i])
      ggdep <- ggdep + ggplot2::geom_line(data = pred_data, ggplot2::aes(x = x, y = y), col = 'grey')
  }

  ggdep <- ggdep + ggplot2::geom_line(data = data.frame(x = x_seq, y = object$preds_mean), 
                    ggplot2::aes(x = x, y = y), col = '#08cbba', linewidth = 1.5)
  ggdep <- ggdep + ggplot2::geom_point(data = data.frame(x = object$xj, y = -5), 
                    ggplot2::aes(x = x, y = y), col = 'black', size = 1,shape = 108)
  ggdep <- ggdep + ggplot2::ylab('f(x)')
  ggdep <- ggdep + ggplot2::xlim(quantile(object$xj, 0.05), quantile(object$xj, 0.95))
  if(is.character(object$j)){
    ggdep <- ggdep + ggplot2::xlab(object$j)
  }else{
    ggdep <- ggdep + ggplot2::xlab(paste('x', object$j, sep = ''))
  }
  ggdep + ggplot2::ylim(-5, 6)
}

ranger_fun <- function(object){
    res <- list(model = object)
    class(res) <- 'ranger_fun'
    return(res)
}
predict.ranger_fun <- function(object, newdata){
    return(predict(object$model, newdata)$predictions)
}

source('R/SDForest.r')
library(gridExtra)
library(ggplot2)
library(ranger)

# default experiment
# load data
load('simulation_study/results/default_szenario.RData')
set.seed(2024)

fit2 <- ranger(x = data.frame(data$X), y = data$Y, num.trees = 100, 
  importance = 'impurity', mtry = floor(0.9 * ncol(data$X)))

true_f <- true_function(data$beta, data$j)

dep_f_1 <- condDependence(true_f, data$j[1], data$X)
dep_f_2 <- condDependence(true_f, data$j[2], data$X)
dep_f_3 <- condDependence(true_f, data$j[3], data$X)
dep_f_4 <- condDependence(true_f, data$j[4], data$X)

data$j
sort(fit$var_importance, decreasing = T)[1:6]
sort(fit2$variable.importance, decreasing = T)[1:6]

plot(fit$var_importance / max(fit$var_importance), col = 'blue', 
  ylim = c(0, 1), xlab = 'Variable', ylab = 'Variable importance')
points(fit2$variable.importance / max(fit2$variable.importance), 
  col = 'red', pch = 2)

grid.arrange(plotDep(dep_f_1), plotDep(dep_f_2), 
  plotDep(dep_f_3), plotDep(dep_f_4), ncol = 2, 
  top = 'Conditional dependence of the true function')

grid.arrange(plotDep(dep_1), plotDep(dep_2), 
  plotDep(dep_3), plotDep(dep_4), ncol = 2, 
  top = 'Conditional dependence of the SDForest')


ranger_fit <- ranger_fun(fit2)

dep_r_1 <- condDependence(ranger_fit, data$j[1], data.frame(data$X))
dep_r_2 <- condDependence(ranger_fit, data$j[2], data.frame(data$X))
dep_r_3 <- condDependence(ranger_fit, data$j[3], data.frame(data$X))
dep_r_4 <- condDependence(ranger_fit, data$j[4], data.frame(data$X))


grid.arrange(plotDep(dep_r_1), plotDep(dep_r_2), 
  plotDep(dep_r_3), plotDep(dep_r_4), ncol = 2, 
  top = 'Conditional dependence of the ranger')


gg_regpath <- ggplot()
for(i in 1:ncol(reg_path$varImp_path)){
  gg_regpath <- gg_regpath + geom_line(data = data.frame(x = reg_path$cp, 
    y = reg_path$varImp_path[, i]), aes(x = x, y = y), 
    col = if(i %in% data$j)'#d11010' else 'grey')
}
gg_regpath <- gg_regpath + theme_bw() + xlab('Complexity parameter') + 
  ylab('Variable importance') + ggtitle('Variable importance path')

gg_regpath

gg_stablepath <- ggplot()
for(i in 1:ncol(stable_path$varImp_path)){
  gg_stablepath <- gg_stablepath + geom_line(data = data.frame(x = stable_path$cp, 
    y = stable_path$varImp_path[, i]), aes(x = x, y = y), 
    col = if(i %in% data$j)'#d11010' else 'grey')
}
gg_stablepath <- gg_stablepath + theme_bw() + xlab('Complexity parameter: cp') + 
  ylab('Variable importance') + ggtitle('Stability selection path')

gg_stablepath


plotOOB(reg_path)

cp_min <- reg_path$cp[which.min(reg_path$loss_path[, 'oob SDE'])]
fit_pruned <- prune(fit, cp_min)

sort(fit_pruned$var_importance, decreasing = T)[1:6]

dep_p_1 <- condDependence(fit_pruned, data$j[1])
dep_p_2 <- condDependence(fit_pruned, data$j[2])
dep_p_3 <- condDependence(fit_pruned, data$j[3])
dep_p_4 <- condDependence(fit_pruned, data$j[4])

grid.arrange(plotDep(dep_p_1), plotDep(dep_p_2), 
  plotDep(dep_p_3), plotDep(dep_p_4), ncol = 2, 
  top = 'Conditional dependence of the pruned SDForest')


f_hat <- predict(fit, data.frame(data_test$X))
f_hat2 <- predict(fit2, data.frame(data_test$X))





Q <- get_Q(data_test$X, 'trim')
f_mse <- mean((data_test$f_X - f_hat)^2)
SDE <- mean((Q %*% data_test$Y - Q %*% f_hat)^2)


fit$oob_SDloss
f_mse
SDE
