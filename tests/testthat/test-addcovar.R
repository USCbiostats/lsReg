datafile <- system.file("extdata", "simulated_data.rds", package = "lsReg")
dat <- readRDS(datafile)
mdl <- glm(ylin ~ x1 + x2, data = dat, family = gaussian)
mem <- lsReg(mdl, colstoadd = 1, testtype = "wald")

test_that("addcovar rejects non-lsregmem input", {
  expect_equal(addcovar(list(), as.matrix(dat[, "z5", drop = FALSE])), 1)
  expect_equal(addcovar(1,     as.matrix(dat[, "z5", drop = FALSE])), 1)
})

test_that("addcovar runs without error on valid input", {
  expect_no_error(addcovar(mem, as.matrix(dat[, "z5", drop = FALSE])))
  expect_no_error(addcovar(mem, as.matrix(dat[, "z1", drop = FALSE])))
})
