source('simulation_study/utils.r')

p <- 500
q <- 20
n_test <- 500

N_rep <- 2

n_seq <- seq(100, 1000, 200)
n_seq <- c(100, 1000)


print('start')
start <- Sys.time()
perf_n <- lapply(1:N_rep, function(i) lapply(n_seq, function(n) performance_measure(n, p, q, n_test, eff = NULL)))
save(perf_n, n_seq, file = "simulation_study/results/perf_n.RData")
print('n done')
print(Sys.time() - start)