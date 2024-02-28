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

##### default experiment #####

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
  ylim = c(0, 1), xlab = 'Variable', ylab = 'Variable importance', pch = 3)
points(fit2$variable.importance / max(fit2$variable.importance), 
  col = 'red', pch = 2)
points(data$j, rep(1, length(data$j)), col = 'green', pch = 20)

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

##### Performance depending on the dimensions #####

library(ggplot2)
library(tidyr)
load('simulation_study/results/perf_n.RData')

perf_n <- do.call(rbind, perf_n)
perf_n <- data.frame(perf_n, rownames(perf_n), row.names = NULL)
names(perf_n) <- c(n_seq, 'method')

perf_n <- gather(perf_n, n, error, -method)

ggplot(perf_n, aes(x = n, y = error, col = method)) + 
  geom_boxplot() + theme_bw() + xlab('Number of training samples') + 
  ylab('Mean squared error') + ggtitle('Performance of SDForest and ranger')


load('simulation_study/results/perf_p.RData')

perf_p <- do.call(rbind, perf_p)
perf_p <- data.frame(perf_p, rownames(perf_p), row.names = NULL)
names(perf_p) <- c(p_seq, 'method')

perf_p <- gather(perf_p, p, error, -method)

ggplot(perf_p, aes(x = p, y = error, col = method)) + 
  geom_boxplot() + theme_bw() + xlab('Number of features') + 
  ylab('Mean squared error') + ggtitle('Performance of SDForest and ranger')


load('simulation_study/results/perf_q.RData')

perf_q <- do.call(rbind, perf_q)
perf_q <- data.frame(perf_q, rownames(perf_q), row.names = NULL)
names(perf_q) <- c(q_seq, 'method')

perf_q <- gather(perf_q, q, error, -method)

ggplot(perf_q, aes(x = q, y = error, col = method)) + 
  geom_boxplot() + theme_bw() + xlab('Number of confounders') + 
  ylab('Mean squared error') + ggtitle('Performance of SDForest and ranger')


##### Regularization performance #####

load('simulation_study/results/regularization_performance.RData')

res_reg_mean <- data.frame(apply(simplify2array(res_reg), 1:2, mean))
res_reg_u <- data.frame(apply(simplify2array(res_reg), 1:2, quantile, prob = 0.95))
res_reg_l <- data.frame(apply(simplify2array(res_reg), 1:2, quantile, prob = 0.05))

res_reg_mean <- gather(res_reg_mean, key = 'type', value = 'mean', -cp)

res_reg_u <- gather(res_reg_u, key = 'type', value = 'u', -cp)

res_reg_l <- gather(res_reg_l, key = 'type', value = 'l', -cp)

res <- merge(res_reg_mean, res_reg_u)
res <- merge(res, res_reg_l)

library(ggplot2)

ggplot(res, aes(x = cp, y = mean)) + 
  geom_line(aes(col = type)) + 
  geom_ribbon(aes(ymin = l, ymax = u, fill = type), alpha = 0.2) + 
  theme_bw() + xlab('Complexity parameter') + ylab('Mean squared error') + 
  ggtitle('Regularization performance of SDForest')





