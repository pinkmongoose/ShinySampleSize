#Sampling.r
#Darren Green
#14/01/2021

#pop <- list(N=50, R=5) # population size, expected number of reactors
#test <- list(sens=0.95, spec=0.98) # test parameters
#target <- list(sens=0.95, spec=0.95) # intended herd values

dbetabinom <- function(k,N,eta,theta) {
  eta <- eta + 1e-10
  theta <- theta + 1e-10
  #distribution of +units out of a population of size N based on beta distributed prevalence
  x <- exp(
    lgamma(N+1)+lgamma(k+eta)+lgamma(N-k+theta)+lgamma(eta+theta)
    -lgamma(k+1)-lgamma(N-k+1)-lgamma(N+eta+theta)-lgamma(eta)-lgamma(theta)
  )
  return(ifelse(is.na(x),0,x))
}

sdbeta <- function(a,b) sqrt(a)*sqrt(b)/(a+b)/sqrt(a+b+1)
sdbeta2 <- function(p,n) sqrt(p*(1-p)/(n+1))

dtails <- function(d,x) {
  n <- length(d)
  up <- array(0,n)
  down <- up
  for (i in 1:n) up[i] <- sum(d[1:i])
  lo <- n
  for (i in 1:n) {
    if (up[i]>=x) {lo <- i-1; break}
  }
  for (i in n:1) down[i] <- sum(d[i:n])
  hi <- 1
  for (i in n:1) {
    if (down[i]>=x) {hi <- i-1; break}
  }
  return (c(lo,hi))
}

dt <- function(n.pos, n.neg, test) {
  # distribution of + tests out of n.pos/n.neg +/- samples.
  t.convolution <- array(0, n.pos + n.neg + 1)
  if (test$Testtype=="BetaBinomial") {
    t.pos <- dbetabinom(0:n.pos, n.pos, test$beta*test$betaN, (1-test$beta)*test$betaN)
    t.neg <- dbetabinom(0:n.neg, n.neg, (1-test$alpha)*test$alphaN, test$alpha*test$alphaN)
    pos.tails <- dtails(t.pos, 1-test$pConv)
    neg.tails <- dtails(t.neg, 1-test$pConv)
    lo.pos <- pos.tails[1]
    hi.pos <- pos.tails[2]
    lo.neg <- neg.tails[1]
    hi.neg <- neg.tails[2]
  } else {
    t.pos <- dbinom(0:n.pos, n.pos, test$beta)
    t.neg <- dbinom(0:n.neg, n.neg, 1 - test$alpha)
    lo.pos <- qbinom(1-test$pConv, n.pos, test$beta)
    hi.pos <- qbinom(test$pConv, n.pos, test$beta)
    lo.neg <- qbinom(1-test$pConv, n.neg, 1 - test$alpha)
    hi.neg <- qbinom(test$pConv, n.neg, 1 - test$alpha)    
  }
  for (i in lo.pos:hi.pos)
    for (j in lo.neg:hi.neg) {
      t.convolution[i + j + 1] <- t.convolution[i + j + 1] + t.pos[i + 1] * t.neg[j + 1]
    }
  return(t.convolution)
}

dT <- function(test, pop, n) {
  # distribution of + tests out of sampled population.
  T <- array(0, n + 1)
  if (test$Sampletype=="Hypergeometric") {
    dn <- dhyper(0:n, pop$R, pop$N - pop$R, n)
  } else {
    dn <- dbinom(0:n, n, pop$R/pop$N)
  }
  # distribution of + samples out of sampled pop
  for (m in 0:n) {
    dtm <- dt(m, n - m, test)
    T <- T + dtm * dn[m + 1]
  }
  return(T)
}

herd.test <- function(test, pop, n, cutpoint) {
  dTn <- dT(test, pop, n)
  B <- sum(dTn[(cutpoint + 2):(n + 1)])
  dT0n <- dbinom(0:n, n, 1 - test$alpha)
  A <- 1 - sum(dT0n[(cutpoint + 2):(n + 1)])
  return(list(B=B, A=A))
}

herd.test.all <- function(test, pop, n) {
  dTn <- dT(test, pop, n)
  dT0n <- dbinom(0:n, n, 1 - test$alpha)
  B <- sapply(0:(n - 1), function(cutpoint) {
    sum(dTn[(cutpoint + 2):(n + 1)])
  })
  A <- sapply(0:(n - 1), function(cutpoint) {
    1 - sum(dT0n[(cutpoint + 2):(n + 1)])
  })
  return(cbind(B,A))
}

DrawROC <- function() {
  if (!length(D$H))
    plot.new()
  else {
    plot(
#      1 - D$H$A, D$H$B,
      NULL,
      type = "o",
      xlim = c(0, 1), ylim = c(0, 1),
      xlab = "1-specificity", ylab = "sensitivity",
      main = paste("ROC for n =", D$Hn),
      cex.lab=1.5, cex.axis=1.5, cex.main=1.5
    )
    rect(0,0,1-D$HAreq,1,col="#ddffdd",border=NA)
    rect(1-D$HAreq,D$HBreq,1,1,col="#ddddff",border=NA)
    rect(0,D$HBreq,1-D$HAreq,1,col="#ddffff",border=NA)
    lines(y = c(D$HBreq, D$HBreq), x = c(0, 1), col = "blue")
    lines(x = c(1 - D$HAreq, 1 - D$HAreq), y = c(0, 1), col = "blue")
    lines(y = c(0, 1), x = c(0, 1), col = "red")
    lines(x=1-D$H$A, y=D$H$B)
    points(x=1-D$H$A, y=D$H$B)
    lines(x = c(1-D$H$A[1],1), y=c(D$H$B[1],1), lty=3)        
    text(1 - D$H$A, D$H$B, labels = 0:(D$Hn - 1), cex = 1.5, pos = 4, col = "black")

  }
}

imperfect.testing <- function(test, pop, target, batch=FALSE) {
  warn <- ""
  lo.n <- 1
  hi.n <- pop$N
  n <- 1
  cutpoint <- 0
  repeat {
    h <- herd.test(test, pop, n, cutpoint)
    if (batch==F) setProgress(value = n)
    if (h$B < target$Breq) {
      lo.n <- n
    } else {
      hi.n <- n
    }
    n <- ceiling(exp(0.5 * (log(lo.n) + log(hi.n))) - 1e-5)
    if (hi.n - lo.n < 2) {
      h <- herd.test(test, pop, n, cutpoint)
      if (h$B < target$Breq) {
        warn <- "Cannot achieve desired sensitivity by sampling every unit."
        warning(warn)
        res <- c(0, 0, 0, 0, pop$N, pop$R)
        D$results <- res
        D$warn <- warn
        D$H <- NULL
        return(res)
      }
      if (h$A < target$Areq) {
        cutpoint <- cutpoint + 1
        if (cutpoint > pop$N) {
          warn <- "Cannot achieve desired specificity by sampling every unit."
          warning(warn)
          res <- c(0, 0, 0, 0, pop$N, pop$R)
          D$results <- res
          D$H <- NULL
          return(res)
        }
        hi.n <- pop$N
      } else {
        res <- c(n, cutpoint, h$B, h$A, pop$N, pop$R)
        D$results <- res
        D$H <- as.data.frame(herd.test.all(test, pop, n))
        D$HBreq <- target$Breq
        D$HAreq <- target$Areq
        D$Hn <- n
        D$warn <- warn
        #        cat(sum(tdist(test, pop, n)))
        return(res)
      }
    }
  }
}

OutOfRange <- function(x, min, max, name) {
  if (x > max)
    showModal(modalDialog(title = "Error", paste(name, "is too high.")))
  else if (x < min)
    showModal(modalDialog(title = "Error", paste(name, "is too low.")))
  else
    return(F)
  return(T)
}

RunModel <- function() {
  D$err <- F
  if (OutOfRange(input$N, 1, 25000, "Population size")) {
    D$err <- T
    return()
  }
  if (OutOfRange(input$R, 0, input$N, "Number of reactors")) {
    D$err <- T
    return()
  }
  if (OutOfRange(input$beta, 0, 1, "Test sensitivity")) {
    D$err <- T
    return()
  }
  if (OutOfRange(input$alpha, 0, 1, "Test specificity")) {
    D$err <- T
    return()
  }
  if (OutOfRange(input$Breq, 0, 0.98, "Herd sensitivity")) {
    D$err <- T
    return()
  }
  if (OutOfRange(input$Areq, 0, 0.98, "Herd specificity")) {
    D$err <- T
    return()
  }
  if (OutOfRange(input$betaN, 1, Inf, "Sensitivity sample size")) {
    D$err <- T
    return()
  }
  if (OutOfRange(input$alphaN, 1, Inf, "Specificity sample size")) {
    D$err <- T
    return()
  }
  if (OutOfRange(input$pConv, 0.99, 1, "Coverage for convolution")) {
    D$err <- T
    return()
  }
  
  pop <- list(N = input$N, R = input$R)
  pop$N <- floor(pop$N)
  if (pop$R < 1) {
    pop$R = floor(pop$R*pop$N)
  } else pop$R = floor(pop$R)
  if (pop$R/pop$N * input$beta < 1 - input$alpha) {
    showModal(modalDialog(title = "Error", paste("Pre-run check: Specificity is too low.")))
    D$err <- T
    return()
  }
  imperfect.testing(
    list(beta = input$beta, alpha = input$alpha, Testtype = input$Testtype,
         betaN = input$betaN, alphaN = input$alphaN, pConv = input$pConv, Sampletype=input$Sampletype),
    pop,
    list(Breq = input$Breq, Areq = input$Areq),
    batch=F
  )
}
