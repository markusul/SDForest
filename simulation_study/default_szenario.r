source("R/SDForest.r")

p <- 400
n <- 400

set.seed(2024)
data <- simulate_data_nonlinear(20, p, n, 4)

fit <- SDForest(x = data$X, y = data$Y, cp = 0, multicore = T)
print('fit done')

reg_path <- regPath(fit, oob = T, multicore = T)
print('reg_path done')

stable_path <- stabilitySelection(fit, multicore = T)
print('stable_path done')

dep_1 <- condDependence(fit, data$j[1], multicore = T)
dep_2 <- condDependence(fit, data$j[2], multicore = T)
dep_3 <- condDependence(fit, data$j[3], multicore = T)
dep_4 <- condDependence(fit, data$j[4], multicore = T)
print('dep done')

save(fit, reg_path, stable_path, dep_1, dep_2, dep_3, dep_4, file = "simulation_study/results/default_szenario.RData")
print('save done')
