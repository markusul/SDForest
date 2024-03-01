source('simulation_study/utils.r')

n <- 500
q <- 20
n_test <- 500

N_rep <- 2

p_seq <- seq(100, 1000, 200)
p_seq <- c(100, 1000)

print('start')
start <- Sys.time()
perf_p <- lapply(1:N_rep, function(i) lapply(p_seq, function(p) performance_measure(n, p, q, n_test, eff = NULL)))
save(perf_p, p_seq, file = "simulation_study/results/perf_p.RData")
print('p done')
print(Sys.time() - start)