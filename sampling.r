#Sampling.r
#Darren Green
#21/05/2019

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

herd.test.all <- function(test, pop, n) {
  t.sens <- tdist(test, pop, n)
  t.spec <- dbinom(0:n, n, 1-test$spec)
  sens <- sapply(0:(n-1), function(cutpoint) {
    sum(t.sens[(cutpoint+2):(n+1)])
  })
  spec <- sapply(0:(n-1), function(cutpoint) {
    1-sum(t.spec[(cutpoint+2):(n+1)])
  })
  return(cbind(sens,spec))
}

DrawROC <- function() {
  if (!length(D$H)) plot.new()
  else {
    plot(1-D$H$spec,D$H$sens,type="o",xlim=c(0,1),ylim=c(0,1),xlab="1-specificity",ylab="sensitivity",main=paste("ROC for n =",D$Hn))
    text(1-D$H$spec,D$H$sens,labels=0:(D$Hn-1),cex=0.7,pos=3)
    lines(x=c(1-D$Hspec,1-D$Hspec),y=c(0,1),col="red")
    lines(y=c(D$Hsens,D$Hsens),x=c(0,1),col="red")
    lines(y=c(0,1),x=c(0,1),col="blue")
  }
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
    # cat(n, cutpoint, h$sens, h$spec, "\n")
    if (h$sens < target$sens) {
      lo.n <- n
    } else {
      hi.n <- n
    }
    n <- ceiling(exp(0.5*(log(lo.n)+log(hi.n)))-1e-5)
    if (hi.n - lo.n < 2) {
      h <- herd.test(test, pop, n, cutpoint)
      if (h$sens < target$sens) {
        warn <- "Cannot achieve desired sensitivity by sampling every unit."
        warning(warn)
        res <- c(0,0,0,0,pop$N,pop$R)
        D$results <- res
        D$warn <- warn
        D$H <- NULL
        return(res)
      }
      if (h$spec < target$spec) {
        cutpoint <- cutpoint + 1
        if (cutpoint>pop$N) {
          warn <- "Cannot achieve desired specificity by sampling every unit."
          warning(warn)
          res <- c(0,0,0,0,pop$N,pop$R)
          D$results <- res
          D$H <- NULL
          return(res)
        }
        hi.n <- pop$N
      } else {
        res <- c(n,cutpoint,h$sens,h$spec,pop$N,pop$R)
        D$results <- res
        D$H <- as.data.frame(herd.test.all(test, pop, n))
        D$Hsens <- target$sens
        D$Hspec <- target$spec
        D$Hn <- n
        D$warn <- warn
        cat(sum(tdist(test,pop,n)))
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
  if (OutOfRange(input$N,1,25000,"Population size")) {D$err<-T; return()}
  if (OutOfRange(input$R,0,input$N,"Number of reactors")) {D$err<-T; return()}
  if (OutOfRange(input$TSens,0,1,"Test sensitivity")) {D$err<-T; return()}
  if (OutOfRange(input$TSpec,0,1,"Test specificity")) {D$err<-T; return()}
  if (OutOfRange(input$HSens,0,0.98,"Herd sensitivity")) {D$err<-T; return()}
  if (OutOfRange(input$HSpec,0,0.98,"Herd specificity")) {D$err<-T; return()}
  pop <- list(N=input$N,R=input$R)
  if (pop$R<1) {
    pop$R<-floor(pop$R*pop$N)
    pop$N<-floor(pop$N)
  } else {
    pop$R<-floor(pop$R)
    pop$N<-floor(pop$N)
  }
  if (pop$R/pop$N*input$TSens < 1-input$TSpec) {
    showModal(modalDialog(title="Error",paste("Specificity is too low.")))
    D$err<-T; return()
  }
  imperfect.testing(list(sens=input$TSens,spec=input$TSpec),
                    pop,
                    list(sens=input$HSens,spec=input$HSpec))
}
