#' @title <DeepCirCode_test>
#'
#' @param x_test a 3D array
#' @param y_test a 2D array
#' @param species a string to be "human", or "mouse", or "fly"
#' @return NULL
#' @description This function test DeepCirCode model with test case
#' @export

DeepCirCode_test <- function(x_test,y_test,species){
  set.seed(123);
  model <- keras_model_sequential();
  model %>%

    layer_conv_1d(filters = 128, kernel_size = 12, strides = 1,padding = "valid",input_shape = c(200,4), activation = "relu") %>%
    layer_dropout(0)%>%

    layer_conv_1d(filters = 128, kernel_size = 6, strides = 1,activation = 'relu') %>%
    layer_dropout(0.5) %>%

    layer_max_pooling_1d(pool_size =4, strides =4)%>%
    layer_dropout(0.5)%>%

    layer_flatten() %>%
    layer_dense(units = 2, activation = "softmax") %>%

    summary(model)

  model %>% compile(
    loss = "binary_crossentropy",
    optimizer = optimizer_rmsprop(),
    metrics = c("accuracy")
  )

  ### load the saved best model
  cat("\n\nloading the best DeepCirCode model...\n")

  # load the model saved in the DeepCirCode package
  modelName <- paste("DeepCirCode_bestmodel","_",species,".hdf5", sep = "")
  path <- system.file(modelName, package ="DeepCirCode")

  model <- load_model_hdf5(filepath = path)
  cat("\n\nDeepCirCode model loaded :D\n")

  ### generate predictions
  cat("\n\nPredicting test set...\n")
  prediction_class <- model %>% predict_classes(x_test)
  label_class <- y_test[ ,2]
  CM <- table(label_class, prediction_class)
  cat("\n\nconfusion matrix:\n")
  print(CM);

  ### performance
  cat("\n\nDeepCirCode performance evaluation:\n")
  result <- model %>% evaluate(x_test, y_test, verbose = 0)
  print(result)

  ### ROC
  classes <- model %>% predict_proba(x_test);
  pred <- prediction(classes[,2],y_test[,2])

  # draw ROCR curve
  perf <- performance(pred, "tpr","fpr");
  plot(perf, col = "blue", lty =3, lwd=3, main = "ROC Curve ");
  abline(a = 0, b = 1, lty = 3)

  # print the auc vaules
  auc <- performance(pred, "auc")
  cat("\n\nROC_AUC value is:\n")
  print(auc@y.values)
}

