source('simulation_study/utils.r')

n <- 500
p <- 500

n_test <- 500

N_rep <- 10

seq <- seq(0, 100, 20)
#seq <- c(0, 100)


print('start')
start <- Sys.time()
perf <- lapply(1:N_rep, function(i) lapply(seq, function(q) performance_measure(n, p, q, n_test, eff = NULL)))
save(perf, seq, file = gsub('[ :]', '_', paste("simulation_study/results/perf_q/", date(), '.RData', sep='')))
print('q done')
print(Sys.time() - start)
