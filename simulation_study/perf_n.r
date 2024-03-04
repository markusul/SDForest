source('simulation_study/utils.r')

p <- 500
q <- 20
n_test <- 500

N_rep <- 10

seq <- seq(100, 1000, 200)

print('start')
start <- Sys.time()
perf <- lapply(1:N_rep, function(i) lapply(seq, function(n) performance_measure(n, p, q, n_test, eff = NULL)))
save(perf, seq, file = gsub('[ :]', '_', paste("simulation_study/results/perf_n/", date(), '.RData', sep='')))
print('n done')
print(Sys.time() - start)