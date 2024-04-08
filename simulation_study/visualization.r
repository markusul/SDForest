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
library(ggsci)
library(ggpubr)

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
imp_1 <- fit$var_importance / max(fit$var_importance)
imp_2 <- fit2$variable.importance / max(fit2$variable.importance)

sort(imp_1, decreasing = T)[1:6]
sort(imp_2, decreasing = T)[1:6]

plot(imp_1, col = 'blue', 
  ylim = c(0, 1), xlab = 'Variable', ylab = 'Variable importance', pch = 3)
points(imp_2, col = 'red', pch = 2)
points(data$j, rep(1, length(data$j)), col = 'green', pch = 20)

true_imp <- rep('blue', length(imp_1))
true_imp[data$j] <- 'green'
plot(log(imp_1), log(imp_2), col = true_imp, pch = 20,
  xlab = 'Variable importance SDForest', ylab = 'Variable importance ranger')


ggdep1 <- plotDep(dep_1) + 
  geom_line(aes(x = dep_f_1$x_seq, y = dep_f_1$preds_mean, col = 'red'), linewidth = 1.5) + 
  ggplot2::labs(col = "") + 
  ggplot2::scale_color_manual(values = c("red"), labels = c("True Function"))

ggdep2 <- plotDep(dep_2) + 
  geom_line(aes(x = dep_f_2$x_seq, y = dep_f_2$preds_mean, col = 'red'), linewidth = 1.5) + 
  ggplot2::labs(col = "") + 
  ggplot2::scale_color_manual(values = c("red"), labels = c("True Function"))

ggdep3 <- plotDep(dep_3) +
  geom_line(aes(x = dep_f_3$x_seq, y = dep_f_3$preds_mean, col = 'red'), linewidth = 1.5) + 
  ggplot2::labs(col = "") + 
  ggplot2::scale_color_manual(values = c("red"), labels = c("True Function"))

ggdep4 <- plotDep(dep_4) +
  geom_line(aes(x = dep_f_4$x_seq, y = dep_f_4$preds_mean, col = 'red'), linewidth = 1.5) + 
  ggplot2::labs(col = "") + 
  ggplot2::scale_color_manual(values = c("red"), labels = c("True Function"))

ggarrange(ggdep1, ggdep2, ggdep3, ggdep4,
  ncol = 2, nrow = 2, common.legend = T, legend = 'bottom')

ranger_fit <- ranger_fun(fit2)

dep_r_1 <- condDependence(ranger_fit, data$j[1], data.frame(data$X))
dep_r_2 <- condDependence(ranger_fit, data$j[2], data.frame(data$X))
dep_r_3 <- condDependence(ranger_fit, data$j[3], data.frame(data$X))
dep_r_4 <- condDependence(ranger_fit, data$j[4], data.frame(data$X))

ggdep1_r <- plotDep(dep_r_1) + 
  geom_line(aes(x = dep_f_1$x_seq, y = dep_f_1$preds_mean, col = 'red'), linewidth = 1.5) + 
  ggplot2::labs(col = "") + 
  ggplot2::scale_color_manual(values = c("red"), labels = c("True Function"))

ggdep2_r <- plotDep(dep_r_2) +
  geom_line(aes(x = dep_f_2$x_seq, y = dep_f_2$preds_mean, col = 'red'), linewidth = 1.5) + 
  ggplot2::labs(col = "") + 
  ggplot2::scale_color_manual(values = c("red"), labels = c("True Function"))

ggdep3_r <- plotDep(dep_r_3) + 
  geom_line(aes(x = dep_f_3$x_seq, y = dep_f_3$preds_mean, col = 'red'), linewidth = 1.5) + 
  ggplot2::labs(col = "") + 
  ggplot2::scale_color_manual(values = c("red"), labels = c("True Function"))

ggdep4_r <- plotDep(dep_r_4) +
  geom_line(aes(x = dep_f_4$x_seq, y = dep_f_4$preds_mean, col = 'red'), linewidth = 1.5) + 
  ggplot2::labs(col = "") + 
  ggplot2::scale_color_manual(values = c("red"), labels = c("True Function"))

ggarrange(ggdep1_r, ggdep2_r, ggdep3_r, ggdep4_r,
  ncol = 2, nrow = 2, common.legend = T, legend = 'bottom')

gg_regpath <- ggplot()
for(i in 1:ncol(reg_path$varImp_path)){
  gg_regpath <- gg_regpath + geom_line(data = data.frame(x = reg_path$cp, 
    y = reg_path$varImp_path[, i]), aes(x = x, y = y), 
    col = if(i %in% data$j)'#d11010' else 'grey')
}
gg_regpath <- gg_regpath + theme_bw() + xlab('Regularization: cp') + 
  ylab('Variable importance') + ggtitle('Variable importance path')

gg_regpath

gg_stablepath <- ggplot()
for(i in 1:ncol(stable_path$varImp_path)){
  gg_stablepath <- gg_stablepath + geom_line(data = data.frame(x = stable_path$cp, 
    y = stable_path$varImp_path[, i]), aes(x = x, y = y), 
    col = if(i %in% data$j)'#d11010' else 'grey')
}
gg_stablepath <- gg_stablepath + theme_bw() + xlab('Regularization: cp') + 
  ylab('Variable importance') + ggtitle('Stability selection path')

gg_stablepath

##### Performance depending on the dimensions #####
agg_fun <- function(x){
  mean(x**2)
}


load_perf <- function(path, agg_fun){
  load(path)
  perf <- lapply(perf, function(n) cbind(t(sapply(n, function(x)sapply(x, agg_fun))), seq = seq))

  perf <- do.call(rbind, perf)
  perf <- data.frame(perf)
  perf$seq <- as.factor(perf$seq)

  perf <- gather(perf, 'method', error, -seq)
  perf
}

library(ggplot2)
library(tidyr)

files <- list.files('simulation_study/results/perf_n')
length(files)
perf_n <- lapply(paste0('simulation_study/results/perf_n/', files), 
  load_perf, agg_fun = agg_fun)

perf_n <- do.call(rbind, perf_n)

gg_n <- ggplot(perf_n, aes(x = seq, y = error, fill = method)) + 
  geom_boxplot(outlier.size = 0.4) + theme_bw() + xlab('Number of training samples') + 
  ylab('Mean squared error') + scale_fill_tron()

gg_n
ggsave(filename = "simulation_study/figures/n.jpeg", plot = gg_n, width = 6, height = 4)

files <- list.files('simulation_study/results/perf_p')
length(files)

perf_p <- lapply(paste0('simulation_study/results/perf_p/', files), 
  load_perf, agg_fun = agg_fun)

perf_p <- do.call(rbind, perf_p)

gg_p <- ggplot(perf_p, aes(x = seq, y = error, fill = method)) + 
  geom_boxplot(outlier.size = 0.4) + theme_bw() + xlab('Number of covariates') + 
  ylab('Mean squared error') + scale_fill_tron()

gg_p
ggsave(filename = "simulation_study/figures/p.jpeg", plot = gg_p, width = 6, height = 4)



files <- list.files('simulation_study/results/perf_q')
length(files)

perf_q <- lapply(paste0('simulation_study/results/perf_q/', files), 
  load_perf, agg_fun = agg_fun)

perf_q <- do.call(rbind, perf_q)

gg_q <- ggplot(perf_q, aes(x = seq, y = error, fill = method)) + 
  geom_boxplot(outlier.size = 0.4) + theme_bw() + xlab('Number of confounders') + 
  ylab('Mean error') + scale_fill_tron()
gg_q

ggsave(filename = "simulation_study/figures/q.jpeg", plot = gg_q, width = 6, height = 4)

files <- list.files('simulation_study/results/perf_max')
length(files)

perf_max <- lapply(paste0('simulation_study/results/perf_max/', files), 
  load_perf, agg_fun = agg_fun)

perf_max <- do.call(rbind, perf_max)

gg_max <- ggplot(perf_max, aes(x = seq, y = error, fill = method)) + 
  geom_boxplot(outlier.size = 0.4) + theme_bw() + xlab('Subsample size') + 
  ylab('Mean error') + scale_fill_tron()
gg_max

ggsave(filename = "simulation_study/figures/max.jpeg", plot = gg_max, width = 6, height = 4)

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


gg_reg <- ggplot(res, aes(x = cp, y = mean)) + 
  geom_line(aes(col = type, linetype = type)) + 
  geom_ribbon(aes(ymin = l, ymax = u, fill = type), alpha = 0.2) + 
  theme_bw() + xlab('Regularization: cp') + ylab('Mean squared error') +
  guides(fill = guide_legend(title = NULL), linetype = guide_legend(title = NULL), 
    col = guide_legend(title = NULL))

gg_reg
ggsave(filename = "simulation_study/figures/reg.jpeg", plot = gg_reg, width = 6, height = 4)


#### Sparsity performance ####

files <- list.files('simulation_study/results/perf_eff')
length(files)
perf_eff <- lapply(paste0('simulation_study/results/perf_eff/', files), 
  load_perf, agg_fun = agg_fun)

perf_eff <- do.call(rbind, perf_eff)

gg_eff <- ggplot(perf_eff, aes(x = seq, y = error, fill = method)) + 
  geom_boxplot(outlier.size = 0.4) + theme_bw() + xlab("Number of affected covariates") + 
  ylab('Mean squared error') + scale_fill_tron()

gg_eff
ggsave(filename = "simulation_study/figures/eff.jpeg", plot = gg_eff, width = 6, height = 4)
