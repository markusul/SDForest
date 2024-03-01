source('simulation_study/utils.r')

n <- 500
p <- 500
q <- 20
n_test <- 500

N_rep <- 2

eff_seq <- seq(0, 499, 50)
eff_seq <- c(0, 499)


print('start')
start <- Sys.time()
perf_eff <- lapply(1:N_rep, function(i) lapply(eff_seq, function(eff) performance_measure(n, p, q, n_test, eff)))
save(perf_eff, eff_seq, file = "simulation_study/results/perf_eff.RData")
print('done')
print(Sys.time() - start)