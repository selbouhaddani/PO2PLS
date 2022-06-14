test_that("All functions", {
  expect_error(generate_params(20,21,2,1,1,0.1), NA)
  parms <- generate_params(20,21,2,1,1,0.1)
  expect_error(generate_data(50, parms), NA)
  Dat <- generate_data(50, parms)
  expect_error(PO2PLS(Dat$X, Dat$Y, 2, 1, 1, 2), NA)
  fit <- PO2PLS(Dat$X, Dat$Y, 2, 1, 1, 2)
  expect_error(E_step_slow(Dat$X, Dat$Y, parms), NA)
  expect_equal(PO2PLS:::E_step(Dat$X,Dat$Y,parms), PO2PLS:::E_step_slow(Dat$X,Dat$Y,parms))
  expect_error(PO2PLS(Dat$X, Dat$Y, 2, 1, 1, 2,init_param = "u"), NA)
  expect_error(jitter_params(parms), NA)
  expect_error(PO2PLS(Dat$X, Dat$Y, 2, 1, 1, 2, init_param = "r" ), NA)
  expect_error(print(summary(PO2PLS(Dat$X, Dat$Y, 2, 1, 1, 2, init_param = "r" ))), NA)
  expect_error(print(PO2PLS(Dat$X, Dat$Y, 2, 1, 1, 2, init_param = "r" )), NA)
  expect_error(LRT(PO2PLS(Dat$X, Dat$Y, 2, 1, 1, 2, init_param = "r" ), PO2PLS(Dat$X, Dat$Y, 1, 1, 1, 2, init_param = "r" )), NA)
  expect_error(variances_inner.po2m(PO2PLS(Dat$X, Dat$Y, 2, 1, 1, 2), Dat$X, Dat$Y), NA)
  expect_error(bootstrap_inner.po2m(fit, Dat$X, Dat$Y, 1, 2, steps=2), NA)
})





