context('createItem')

test_that('old2PL', {
    name <- 'old2PL'
    par <- c(a = .5, b = -2)
    est <- c(TRUE, TRUE)
    P.old2PL <- function(par,Theta,ncat){
        a <- par[1]
        b <- par[2]
        P1 <- 1 / (1 + exp(-1.702*a*(Theta - b)))
        cbind(1-P1, P1)
    }
    lbound <- c(-Inf, -Inf)
    ubound <- c(Inf, Inf)
    
    x <- createItem(name, par=par, est=est, lbound=lbound, ubound=ubound, P=P.old2PL)
    
    #So, let's estimate it!
    dat <- expand.table(LSAT7)
    sv <- mirt(dat, 1, c(rep('2PL',4), 'old2PL'), customItems=list(old2PL=x), pars = 'values', verbose=FALSE)    
    expect_is(sv, 'data.frame')          
    mod <- mirt(dat, 1, c(rep('2PL',4), 'old2PL'), customItems=list(old2PL=x), verbose=FALSE)
    expect_is(mod, 'ConfirmatoryClass')          
    expect_is(coef(mod), 'list')
    mod2 <- confmirt(dat, 1, c(rep('2PL',4), 'old2PL'), customItems=list(old2PL=x), verbose = FALSE)
    expect_is(mod2, 'ConfirmatoryClass')          
    expect_is(coef(mod2), 'list')    
    
    #' #nonlinear
    name <- 'nonlin'
    par <- c(a1 = .5, a2 = .1, d = 0)
    est <- c(TRUE, TRUE, TRUE)
    P.nonlin <- function(par,Theta,ncat){
      a1 <- par[1]
      a2 <- par[2] 
      d <- par[3]
      P1 <- 1 / (1 + exp(-1.702*(a1*Theta + a2*Theta^2 + d)))
      cbind(1-P1, P1)
    } 
     
    x2 <- createItem(name, par=par, est=est, P=P.nonlin)    
    mod <- mirt(dat, 1, c(rep('2PL',4), 'nonlin'), customItems=list(nonlin=x2), verbose=FALSE)
    expect_is(mod, 'ConfirmatoryClass')          
    expect_is(coef(mod), 'list')    
})


