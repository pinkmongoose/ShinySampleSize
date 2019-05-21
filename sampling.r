

#pop <- list(N=50, R=5) # population size, expected number of reactors
#test <- list(sens=0.95, spec=0.98) # test parameters
#target <- list(sens=0.95, spec=0.95) # intended herd values


sdist <- function(vp, vn, test) { # distribution of + tests out of vp/vn +/- samples.
#  if (vp+vn > 100) {
#      return(dsinib(vp+vn, c(vp,vn), c(test$sens,test$spec)))
  s <- array(0,vp+vn+1)
  dp <- dbinom(0:vp, vp, test$sens)
  dn <- dbinom(0:vn, vn, 1-test$spec)
  lop <- qbinom(0.001, vp, test$sens)
  hop <- qbinom(0.999, vp, test$sens)
  lon <- qbinom(0.001, vn, 1-test$spec)
  hon <- qbinom(0.999, vn, 1-test$spec)
  for (p in lop:hop)
    for (n in lon:hon) {
      s[p+n+1] <- s[p+n+1] + dp[p+1]*dn[n+1]
    }
  return(s)
}

tdist <- function(test, pop, n) { # distribution of + tests out of sampled population.
  t <- array(0,n+1)
  rdist <- dhyper(0:n, pop$R, pop$N-pop$R, n) # distribution of + samples out of sampled pop
  for (p in 0:n) {
    s <- sdist(p, n-p, test)
    t <- t + s * rdist[p+1]
  }
  return(t)
}

herd.test <- function(test, pop, n, cutpoint) {
  t.sens <- tdist(test, pop, n)
  sens <- sum(t.sens[(cutpoint+2):(n+1)])
  t.spec <- dbinom(0:n, n, 1-test$spec)
  spec <- 1-sum(t.spec[(cutpoint+2):(n+1)])
  return(list(sens=sens,spec=spec))
}

imperfect.testing <- function(test, pop, target) {
  warn <- ""
  lo.n <- 1
  hi.n <- pop$N
  n <- 1
  cutpoint <- 0
  repeat {
    h <- herd.test(test, pop, n, cutpoint)
    setProgress(value=n)
    cat(n, cutpoint, h$sens, h$spec, "\n")
    if (h$sens < target$sens) {
      lo.n <- n
    } else {
      hi.n <- n
    }
    n <- ceiling(exp(0.5*(log(lo.n)+log(hi.n))))
    if (hi.n - lo.n < 2) {
      h <- herd.test(test, pop, n, cutpoint)
      if (h$sens < target$sens) {
        warn <- "Cannot achieve desired sensitivity by sampling every unit."
        warning(warn)
        res <- c(0,0,0,0,pop$N,pop$R)
        D$results <- res
        D$warn <- warn
        return(res)
      }
      if (h$spec < target$spec) {
        cutpoint <- cutpoint + 1
        if (cutpoint>pop$N) {
          warn <- "Cannot achieve desired specificity by sampling every unit."
          warning(warn)
          res <- c(0,0,0,0,pop$N,pop$R)
          D$results <- res
          return(res)
        }
        hi.n <- pop$N
      } else {
        res <- c(n,cutpoint,h$sens,h$spec,pop$N,pop$R)
        D$results <- res
        D$warn <- warn
        return(res)
      }      
    }
  }
}

OutOfRange <- function(x,min,max,name) {
  if (x>max) showModal(modalDialog(title="Error",paste(name,"is too high.")))
  else if (x<min) showModal(modalDialog(title="Error",paste(name,"is too low.")))
  else return(F)
  return(T)
}

RunModel <- function() {
  D$err <- F
  if (OutOfRange(input$N,1,10000,"Population size")) D$err<-T
  if (OutOfRange(input$R,1,input$N,"Number of reactors")) D$err<-T
  if (OutOfRange(input$TSens,0,1,"Test sensitivity")) D$err<-T
  if (OutOfRange(input$TSpec,0,1,"Test specificity")) D$err<-T
  if (OutOfRange(input$HSens,0,1,"Herd sensitivity")) D$err<-T
  if (OutOfRange(input$HSpec,0,1,"Herd specificity")) D$err<-T
  if (D$err) return()
  imperfect.testing(list(sens=input$TSens,spec=input$TSpec),
                    list(N=floor(input$N),R=floor(input$R)),
                    list(sens=input$HSens,spec=input$HSpec))
}

st<-Sys.time()
sdist(100,900,test)
et<-Sys.time()
print(et-st)

