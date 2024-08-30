set.seed(1)

X <- matrix(rnorm(50 * 20), nrow = 50)
Y <- rnorm(50)

tree <- SDTree(x = X, y = Y, Q_type = 'no_deconfounding', 
               cp = 0, min_sample = 5)

# does min sample work
expect_true(min(table(tree$predictions)) >= 5)
rpart_tree <- rpart::rpart(y ~ ., data.frame(y = Y, X), 
                           control = rpart::rpart.control(minbucket = 5, 
                                                          cp = 0, 
                                                          minsplit = 10, 
                                                          xval = 0))
pruned_tree <- prune(copy(tree), 0.1)
pruned_rpart_tree <- rpart::prune(rpart_tree, 0.1)

# predict = predictions
expect_equal(tree$predictions, as.vector(predict(tree, data.frame(X))))
# changes in model due to pruning and copy of tree before pruning
expect_false(all(tree$predictions == as.vector(predict(pruned_tree, data.frame(X)))))
# equality of tree and rpart tree (checked using predictions)
expect_equal(tree$predictions, as.vector(predict(rpart_tree)))
# equality of pruned tree and pruned rpart tree (checked using predictions)
expect_equal(as.vector(predict(pruned_tree, data.frame(X))), predict(pruned_rpart_tree))

partDependence(tree, 1, X, subSample = 10)

