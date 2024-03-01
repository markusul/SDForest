source('simulation_study/utils.r')

n <- 500
p <- 500

n_test <- 500

N_rep <- 2

q_seq <- seq(0, 100, 20)
q_seq <- c(0, 100)


print('start')
start <- Sys.time()
perf_q <- lapply(1:N_rep, function(i) lapply(q_seq, function(q) performance_measure(n, p, q, n_test, eff = NULL)))
save(perf_q, q_seq, file = "simulation_study/results/perf_q.RData")
print('q done')
print(Sys.time() - start)
