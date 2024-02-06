test_that("Risk computaion by DR function", {
risk <-  DR(5:10, "Pouillot", 1:11)
  expect_true(is.vector(risk)==FALSE)
})
