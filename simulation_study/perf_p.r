source('simulation_study/utils.r')

n <- 500
q <- 20
n_test <- 500

N_rep <- 10

seq <- seq(100, 1000, 200)

print('start')
start <- Sys.time()
perf <- lapply(1:N_rep, function(i) lapply(seq, function(p) performance_measure(n, p, q, n_test, eff = NULL)))
save(perf, seq, file = gsub('[ :]', '_', paste("simulation_study/results/perf_p/", date(), '.RData', sep='')))
print('p done')
print(Sys.time() - start)