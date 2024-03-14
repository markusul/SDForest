source('simulation_study/utils.r')

n <- 500
p <- 500
q <- 20
n_test <- 500

N_rep <- 10

seq <- seq(10, 500, 50)
#eff_seq <- c(0, 499)


print('start')
start <- Sys.time()
perf <- lapply(1:N_rep, function(i) lapply(seq, function(max) performance_measure(n, p, q, n_test, eff = NULL, max)))
save(perf, seq, file = gsub('[ :]', '_', paste("simulation_study/results/perf_max/", date(), '.RData', sep='')))
print('done')
print(Sys.time() - start)