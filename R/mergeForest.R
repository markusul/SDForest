#' Merge two forests
#' 
#' This function merges two forests. 
#' The trees are combined and the variable importance is 
#' calculated as a weighted average of the two forests. 
#' If the forests are trained on the same data, 
#' the predictions and oob_predictions are combined as well.
#' @author Markus Ulmer
#' @param fit1 first /code{SDForest} object
#' @param fit2 second /code{SDForest} object
#' @return merged /code{SDForest} object
#' @export
mergeForest <- function(fit1, fit2){
  if(any(fit1$var_names != fit2$var_names)) 
    stop('forest must be trained using the same covariates')

  len_1 <- length(fit1$forest)
  len_2 <- length(fit2$forest)
  len_new <- len_1 + len_2

  # combine trees
  fit1$forest <- c(fit1$forest, fit2$forest)

  # weighted average of variable importance
  fit1$var_importance <- (fit1$var_importance * len_1 + 
    fit2$var_importance * len_2) / len_new

  # check if both forests are trained on the same data
  # if so, combine predictions and oob_predictions
  if(all(dim(fit1$X) == dim(fit2$X)) && 
     all(fit1$X == fit2$X, fit1$Y == fit2$Y, fit1$Q == fit2$Q, 
         !is.null(fit1$X), !is.null(fit1$Y), !is.null(fit1$Q))){
    
    # weighted average of predictions
    fit1$predictions <- (fit1$predictions * length(fit1$forest) + 
      fit2$predictions * length(fit2$forest)) / len_new

    oob_weights_1 <- sapply(fit1$oob_ind, length)
    oob_weights_2 <- sapply(fit2$oob_ind, length)
    oob_weights_new <- oob_weights_1 + oob_weights_2

    oob_predictions_1 <- fit1$oob_predictions * oob_weights_1
    oob_predictions_2 <- fit2$oob_predictions * oob_weights_2

    fit1$oob_predictions <- rowSums(cbind(oob_predictions_1, oob_predictions_2), 
                                    na.rm = T)
    fit1$oob_predictions <- fit1$oob_predictions / oob_weights_new
    
    # combine what trees are available for oob_i
    fit1$oob_ind <- lapply(1:length(fit1$Y), function(i){
      c(fit1$oob_ind[[i]], fit2$oob_ind[[i]] + len_1)
    })

    # calculate oob_loss and oob_SDloss
    fit1$oob_SDloss <- loss(fit1$Q %*% fit1$Y, fit1$Q %*% fit1$oob_predictions)
    fit1$oob_loss <- loss(fit1$Y, fit1$oob_predictions)
  }else{
    warning('forests migth be trained on different data')
    fit1$X <- NULL
    fit1$Y <- NULL
    fit1$Q <- NULL
    fit1$predictions <- NULL
    fit1$oob_prediction <- NULL
    fit1$oob_ind <- NULL
    fit1$oob_loss <- NULL
    fit1$oob_SDloss <- NULL
  }
  
  fit1
}
