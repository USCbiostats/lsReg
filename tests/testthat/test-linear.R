datafile <- system.file("extdata", "simulated_data.rds", package = "lsReg")
dat  <- readRDS(datafile)
base <- glm(ylin ~ x1 + x2, data = dat, family = gaussian)
xr_z5 <- as.matrix(dat[, "z5", drop = FALSE])
xr_z1 <- as.matrix(dat[, "z1", drop = FALSE])

full_z5 <- glm(ylin ~ x1 + x2 + z5, data = dat, family = gaussian)
full_z1 <- glm(ylin ~ x1 + x2 + z1, data = dat, family = gaussian)

test_that("linear Wald test statistic matches glm", {
  mem <- lsReg(base, colstoadd = 1, testtype = "wald")
  addcovar(mem, xr_z5)
  expect_equal(mem$testvalue[1, 1],
               coef(summary(full_z5))["z5", "t value"], tolerance = 1e-5)
  addcovar(mem, xr_z1)
  expect_equal(mem$testvalue[1, 1],
               coef(summary(full_z1))["z1", "t value"], tolerance = 1e-5)
})

test_that("linear Wald parameter estimate matches glm", {
  mem <- lsReg(base, colstoadd = 1, testtype = "wald")
  addcovar(mem, xr_z5)
  expect_equal(mem$fitdata$betab[1], as.numeric(coef(full_z5)["z5"]), tolerance = 1e-5)
  addcovar(mem, xr_z1)
  expect_equal(mem$fitdata$betab[1], as.numeric(coef(full_z1)["z1"]), tolerance = 1e-5)
})

test_that("linear LRT statistic matches glm log-likelihood ratio", {
  mem <- lsReg(base, colstoadd = 1, testtype = "lrt")
  addcovar(mem, xr_z5)
  lrt_z5 <- as.numeric(-2 * (logLik(base) - logLik(full_z5)))
  expect_equal(mem$testvalue[1], lrt_z5, tolerance = 1e-5)
  addcovar(mem, xr_z1)
  lrt_z1 <- as.numeric(-2 * (logLik(base) - logLik(full_z1)))
  expect_equal(mem$testvalue[1], lrt_z1, tolerance = 1e-5)
})

test_that("linear score test statistic matches statmod::glm.scoretest", {
  mem <- lsReg(base, colstoadd = 1, testtype = "score")
  addcovar(mem, xr_z5)
  expect_equal(mem$testvalue[1, 1],
               statmod::glm.scoretest(base, dat[, "z5"]), tolerance = 1e-5)
  addcovar(mem, xr_z1)
  expect_equal(mem$testvalue[1, 1],
               statmod::glm.scoretest(base, dat[, "z1"]), tolerance = 1e-5)
})

test_that("linear robust score test statistic is finite", {
  mem <- lsReg(base, colstoadd = 1, testtype = "robustscore")
  addcovar(mem, xr_z5)
  expect_true(is.finite(mem$testvalue[1, 1]))
  addcovar(mem, xr_z1)
  expect_true(is.finite(mem$testvalue[1, 1]))
})

test_that("linear robust Wald test statistic is finite", {
  mem <- lsReg(base, colstoadd = 1, testtype = "robustwald")
  addcovar(mem, xr_z5)
  expect_true(is.finite(mem$testvalue[1, 1]))
  addcovar(mem, xr_z1)
  expect_true(is.finite(mem$testvalue[1, 1]))
})

test_that("linear robust statistics are similar to standard under correct model", {
  mem_wald   <- lsReg(base, colstoadd = 1, testtype = "wald")
  mem_rwald  <- lsReg(base, colstoadd = 1, testtype = "robustwald")
  mem_score  <- lsReg(base, colstoadd = 1, testtype = "score")
  mem_rscore <- lsReg(base, colstoadd = 1, testtype = "robustscore")
  addcovar(mem_wald,   xr_z5)
  addcovar(mem_rwald,  xr_z5)
  addcovar(mem_score,  xr_z5)
  addcovar(mem_rscore, xr_z5)
  expect_equal(mem_wald$testvalue[1,1],  mem_rwald$testvalue[1,1],  tolerance = 0.1)
  expect_equal(mem_score$testvalue[1,1], mem_rscore$testvalue[1,1], tolerance = 0.1)
})
