library(testthat)
library(mirt)
library(numDeriv)
options(error = utils::recover, browserNLdisabled=TRUE)

i.count <- 8
d <- rnorm(i.count)
data <- simdata(matrix(rlnorm(i.count, sdlog=.5), ncol=1),
                matrix(c(d+1, d), ncol=2),
                200, 'graded')
spoint <- c(1, 1, -.5)

deriv <- genD(function(param) {
  tech <- list(debugderiv=list(param=1:length(param), value=param))
  fit <- mirt(data, 1, rep('graded',i.count), D=1, calcNull=FALSE, technical=tech)
  fit@MstepLogLik
}, spoint, method.args=list(eps=0.01, d=0.01, r=2))

tech <- list(debugderiv=list(param=1:length(spoint), value=spoint))
fit <- mirt(data, 1, rep('graded',i.count), D=1, calcNull=FALSE, technical=tech)

if (0) {
  for (ix in 1:i.count) {
    ii <- extract.item(fit, ix)
    print(ii@par)
  }
}

if (1) {
  ii <- extract.item(fit, 1)

  np <- length(spoint)
  expect_equal(-ii@gradient, deriv$D[1:np], tolerance=1e-4)
  print(-ii@gradient)
  
  hess <- matrix(NA, nrow=np, ncol=np)
  dx <- np+1
  for (hr in 1:np) {
    hess[hr,1:hr] <- deriv$D[dx:(dx+hr-1)]
    dx <- dx + hr
  }
  ii@hessian[is.na(hess)] <- NA
  expect_equal(hess, -ii@hessian, tolerance=1e-4)
  print(hess)
}
