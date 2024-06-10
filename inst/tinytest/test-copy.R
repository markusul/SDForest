set.seed(1)
X <- matrix(rnorm(50 * 20), nrow = 50)
Y <- rnorm(50)

tree <- SDTree(x = X, y = Y)
pruned_tree <- prune(copy(tree), 0.2)
expect_false(all.equal(tree$predictions, pruned_tree$predictions))