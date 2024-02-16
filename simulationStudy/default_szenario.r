source("R/SDForest.r")

p <- 400
n <- 400

set.seed(2024)
data <- simulate_data_nonlinear(20, p, n, 4)

fit <- SDForest(x = data$X, y = data$Y, cp = 0, multicore = F)

reg_path <- regPath(fit, oob = T, multicore = F)
stable_path <- stabilitySelection(fit, multicore = F)

dep_1 <- condDependence(fit, data$j[1])
dep_2 <- condDependence(fit, data$j[2])
dep_3 <- condDependence(fit, data$j[3])
dep_4 <- condDependence(fit, data$j[4])

save(fit, reg_path, stable_path, dep_1, dep_2, dep_3, dep_4, file = "simulationStudy/results/default_szenario.RData")

