library(testthat)
library(nba)
library(rpart)
library(randomForest)
library(e1071)
library(xgboost)

context("pred_probs range(0,1)")

test_that("GLM predicted probabilities from models are between 0-1",{
  # load data
  data("mtcars")
  # fit glm
  glm_fit <- glm(am ~ mpg, mtcars, family = binomial("logit"))
  # extract fitted value
  pred <- glm_fit$fitted.values
  # expect the maximum of the predicted value to be less than 1 and minimmum to be greater than 1
  expect_lt(max(pred),1)
  expect_gt(min(pred),0)
})

test_that("rpart predicted probabilities from models are between 0-1",{
  # load data
  data("mtcars")
  # fit glm
  rpart_fit <- rpart(am ~ mpg, mtcars)
  # extract fitted value
  pred <- predict(rpart_fit, mtcars)
  # expect the maximum of the predicted value to be less than 1 and minimmum to be greater than 1
  expect_lt(max(pred),1)
  expect_gt(min(pred),0)
})

