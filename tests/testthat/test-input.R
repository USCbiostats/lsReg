test_that("lsreg_input", {
  expect_error(lsReg(), "No base model specified")
  mdl <- 1
  expect_error(lsReg(basemdl = mdl), "Base model is not of class glm")
  mdlfile <- system.file("extdata", "linmodel.rds", package = "lsReg")
  mdl <- readRDS(mdlfile)

  expect_error(lsReg(basemdl = mdl),
               "colstoadd not specified")
  expect_error(lsReg(basemdl = mdl, colstoadd = "A"),
               "colstoadd must be an integer value")
  expect_error(lsReg(basemdl = mdl, colstoadd = 1:2),
               "colstoadd must be a single integer value")
  expect_error(lsReg(basemdl = mdl, colstoadd = 1.1),
               "colstoadd must be an integer value")
  expect_error(lsReg(basemdl = mdl, colstoadd = -1),
               "colstoadd must be a postive integer")
  
  mdl$family$family <- "test"
  expect_error(lsReg(basemdl = mdl),
               "Family of model must be guassian, bionomial or poisson")
})
