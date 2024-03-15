#' mq_pred
#'
#' @title Small Area predictor under M-quantile model
#'
#' @description This functions is used to predict the small area mean under a M-quantile  approach.
#'
#' @usage mq_pred(ys, xs, barX, barXr, area.s, Ni, Q = sort(c(seq(0.006, 0.99, 0.045), 0.5, 0.994, 0.01, 0.02, 0.96, 0.98)), tol = 1e-06, maxit = 100, k=1.345)
#'
#' @param ys the response: a vector of length the number of rows of xs.
#' @param xs a matrix or data frame containing the explanatory variables including the intercept.
#' @param barX m*(p+1) matrix that contains the population means of each of the p auxiliary variables for the m domains. The first column is the intercept.
#' @param barXr m*(p+1) matrix that contains the population means for out of sample units of each of the p auxiliary variables for the m domains. The first column is the intercept.
#' @param barXs m*(p+1) matrix that contains the population means for observed units of each of the p auxiliary variables for the m domains. The first column is the intercept.
#' @param area.s n*1 vector (same size as y in formula) with domain codes.
#' @param Ni vector of domain population sizes.
#' @param Q the quantile(s) to be estimated, this is generally a number strictly between 0 and 1.
#' @param tol the accuracy for the stopping criterion.
#' @param maxit the limit on the number of algorithm iterations.
#' @param k tuning constant used for Huber influence function for the estimation of the coefficients.
#'
#' @import
#' MASS
#'
#' @details In this function the Huber influence function is used by default and the MAD is used by default for estimating the scale.
#'
#' @return Return a list with the following elements:
#' \item{coefficients}{a named vector of coefficients at Qi}
#' \item{MQ}{a named vector of MQ predictor at small area level}
#' \item{MQCD}{a named vector of MQ Chambers and Dunstan predictor at small area level}
#' \item{Qi}{the estimated quantile(s)}
#'
#' @references Chambers, R. and Tzavidis, N. (2006) M-quantile models for small area estimation. Biometrika, 93, 255-268.
#' @references Lahiri, P. and Salvati, N. (2022). A Nested Error Regression Model with High Dimensional Parameter for Small Area Estimation. arXiv:2201.10235, https://doi.org/10.48550/arXiv.2201.10235.
#'
#' @author Nicola Salvati
#'
#' @examples
#' library(pps)
#'
#' m=20
#' sigma2u=3
#' sigma2e<-rep(6,m)
#' ar<-seq(1,m)
#'
#' ni=rep(10,m)
#' Ni=rep(100,m)
#' N=sum(Ni)
#' n=sum(ni)
#'
#' u<-NULL
#' u<-rnorm(m,0,sqrt(sigma2u))
#' u<-rep(u,Ni)
#' e<-NULL
#' for (i in 1:m)e <- c(e,rnorm(Ni[i], 0, sqrt(sigma2e[i])))
#' gr=rep(1:m,each=100)
#' ar=unique(gr)
#' X=matrix(c(rlnorm(N,log(4.5)-0.5,0.5)),nrow=N,ncol=1)
#' beta.tmp<-c(rep(5,N/2),rep(-5,N/2))
#' y=10+beta.tmp*X+u+e
#' pop.matrix<-cbind(y,X,gr)
#' pop<-as.data.frame(pop.matrix)
#' names(pop)<-c("y","x","area")
#' Yibar=tapply(pop[,1],pop[,3],mean)
#'
#' s<-stratsrs(pop$area,ni)
#' x.s=pop[s,]$x
#' y.s=pop[s,]$y
#' regioncode.s=pop[s,]$area
#' pop.r<-pop[-s,]
#'
#' P=sort(c(seq(0.006,0.99,0.045),0.5,0.994,0.01,0.02,0.96,0.98))
#' XMean<-cbind(1,tapply(pop$x,pop$area,mean))
#' XsMean<-cbind(1,tapply(pop$x[s],pop$area[s],mean))
#' XrMean<-cbind(1,tapply(pop$x[-s],pop$area[-s],mean))
#' est<-mq_pred(ys=y.s, xs=cbind(1,x.s), barX=XMean, barXr=XrMean, barXs=XsMean, area.s=regioncode.s, Ni=Ni, Q = P, tol = 1e-06, maxit = 100, k=1.345)
#'
#'
#' @export
mq_pred<- function(ys, xs, barX, barXr, barXs, area.s, Ni, Q = sort(c(seq(0.006, 0.99, 0.045), 0.5, 0.994, 0.01, 0.02, 0.96, 0.98)), tol = 1e-06, maxit = 100, k=1.345)
{
  #ys: response variable
  #xs: design matrix without intercept
  #barX: matrix of average x for all units in the population with intercept
  #barXr: matrix of average x for out of sample units with intercept
  #barXs: matrix of average x for sampled units with intercept
  #area.s: area code for sampled unit
  #Ni: populazion size for each small area
  #Q the quantile(s) to be estimated, this is generally a number strictly between 0 and 1
  #tol the accuracy for the stopping criterion
  #maxit the limit on the number of algorithm iterations
  #k: tuning costant for M-quantile estimation
  require(MASS)
  #Finding the Quantile Orders by Linear Interpolation
  # Assumes that "zerovalinter" function has been already loaded
  "gridfitinter"<-function(y,expectile,Q)
    # computing of the expectile-order of each observation of y by interpolation
  {
    nq<-length(Q)
    diff <- y %*% t(as.matrix(rep(1, nq))) - expectile
    vectordest <- apply(diff, 1, zerovalinter,Q)
    #print(vectordest)
    #qord<-list(ord=c(vectordest))
    #qord
  }
  # COMPUTING OF THE QUANTILE-ORDERS
  "zerovalinter"<-function(y, x)
  {
    if(min(y) > 0) {
      xmin <- x[y == min(y)]
      if(length(xmin) > 0)
        xmin <- xmin[length(xmin)]
      xzero <- xmin
    }

    else {
      if(max(y) < 0) {
        xmin <- x[y == max(y)]
        if(length(xmin) > 0)
          xmin <- xmin[1]
        xzero <- xmin
      }
      else {
        y1 <- min(y[y > 0])
        if(length(y1) > 0)
          y1 <- y1[length(y1)]
        y2 <- max(y[y < 0])
        if(length(y2) > 0)
          y2 <- y2[1]
        x1 <- x[y == y1]
        if(length(x1) > 0)
          x1 <- x1[length(x1)]
        x2 <- x[y == y2]
        if(length(x2) > 0)
          x2 <- x2[1]
        xzero <- (x2 * y1 - x1 * y2)/(y1 - y2)
        xmin <- x1
        if(abs(y2) < y1)
          xmin <- x2
      }
    }
    resu <-  xzero
    resu
  }

  #QRLM function
  QRLM <-function (x, y, case.weights = rep(1, nrow(x)),k=1.345, var.weights = rep(1, nrow(x)), ..., w = rep(1, nrow(x)), init = "ls", psi = psi.huber, scale.est = c("MAD", "Huber", "proposal 2"), k2 = 1.345, method = c("M", "MM"), maxit = 20, acc = 1e-04, test.vec = "resid", q = 0.5)
  { require(MASS)
    irls.delta <- function(old, new) sqrt(sum((old - new)^2)/max(1e-20, sum(old^2)))
    irls.rrxwr <- function(x, w, r) {
      w <- sqrt(w)
      max(abs((matrix(r*w,1,length(r)) %*% x)/sqrt(matrix(w,1,length(r)) %*% (x^2))))/sqrt(sum(w*r^2))
    }
    method <- match.arg(method)
    nmx <- deparse(substitute(x))
    if (is.null(dim(x))) {
      x <- as.matrix(x)
      colnames(x) <- nmx
    }
    else x <- as.matrix(x)
    if (is.null(colnames(x)))
      colnames(x) <- paste("X", seq(ncol(x)), sep = "")
    if (qr(x)$rank < ncol(x))
      stop("x is singular: singular fits are not implemented in rlm")
    if (!(any(test.vec == c("resid", "coef", "w", "NULL")) || is.null(test.vec)))
      stop("invalid testvec")
    if (length(var.weights) != nrow(x))
      stop("Length of var.weights must equal number of observations")
    if (any(var.weights < 0))
      stop("Negative var.weights value")
    if (length(case.weights) != nrow(x))
      stop("Length of case.weights must equal number of observations")
    w <- (w * case.weights)/var.weights
    if (method == "M") {
      scale.est <- match.arg(scale.est)
      if (!is.function(psi))
        psi <- get(psi, mode = "function")
      arguments <- list(...)
      if (length(arguments)) {
        pm <- pmatch(names(arguments), names(formals(psi)), nomatch = 0)
        if (any(pm == 0))
          warning(paste("some of ... do not match"))
        pm <- names(arguments)[pm > 0]
        formals(psi)[pm] <- unlist(arguments[pm])
      }
      if (is.character(init)) {
        if (init == "ls")
          temp <- lm.wfit(x, y, w, method = "qr")
        else if (init == "lts")
          temp <- lqs.default(x, y, intercept = FALSE, nsamp = 200)
        else stop("init method is unknown")
        coef <- temp$coef
        resid <- temp$resid
      }
      else {
        if (is.list(init))
          coef <- init$coef
        else coef <- init
        resid <- y - x %*% coef
      }
    }
    else if (method == "MM") {
      scale.est <- "MM"
      temp <- lqs.default(x, y, intercept = FALSE, method = "S", k0 = 1.548)
      coef <- temp$coef
      resid <- temp$resid
      psi <- psi.bisquare
      if (length(arguments <- list(...)))
        if (match("c", names(arguments), nomatch = FALSE)) {
          c0 <- arguments$c
          if (c0 > 1.548) {
            psi$c <- c0
          }
          else warning("c must be at least 1.548 and has been ignored")
        }
      scale <- temp$scale
    }
    else stop("method is unknown")
    done <- FALSE
    conv <- NULL
    n1 <- nrow(x) - ncol(x)
    if (scale.est != "MM")
      scale <- mad(resid/sqrt(var.weights), 0)
    qest <- matrix(0, nrow = ncol(x), ncol = length(q))
    qwt <- matrix(0, nrow = nrow(x), ncol = length(q))
    qfit <- matrix(0, nrow = nrow(x), ncol = length(q))
    qres <- matrix(0, nrow = nrow(x), ncol = length(q))
    qvar.matrix <- array(rep(0,ncol(x)*ncol(x)),dim=c(ncol(x),ncol(x),length(q)))
    qscale <- NULL
    for(i in 1:length(q)) {
      for (iiter in 1:maxit) {
        if (!is.null(test.vec))
          testpv <- get(test.vec)
        if (scale.est != "MM") {
          if (scale.est == "MAD")
            scale <- median(abs(resid/sqrt(var.weights)))/0.6745
          else {gamma<- 4*k2^2*(1-pnorm(k2))*((1-q[i])^2+q[i]^2) - 4*k2*dnorm(k2)*((1-q[i])^2+q[i]^2) + 4*(1-q[i])^2*(pnorm(0)-(1-pnorm(k2))) + 4*q[i]^2*(pnorm(k2)-pnorm(0))
          scale <- sqrt(sum(pmin(resid^2/var.weights,(k2*scale)^2))/(n1*gamma))
          }
          if (scale == 0) {
            done <- TRUE
            break
          }
        }
        w <- psi(resid/(scale * sqrt(var.weights)),k=k) * case.weights
        ww <- 2 * (1 - q[i]) * w
        ww[resid > 0] <- 2 * q[i] * w[resid > 0]
        w <- ww
        temp <- lm.wfit(x, y, w, method = "qr")
        coef <- temp$coef
        resid <- temp$residuals
        if (!is.null(test.vec))
          convi <- irls.delta(testpv, get(test.vec))
        else convi <- irls.rrxwr(x, wmod, resid)
        conv <- c(conv, convi)
        done <- (convi <= acc)
        if (done)
          break
      }
      if (!done)
        warning(paste("rlm failed to converge in", maxit, "steps at q = ", q[i]))
      qest[, i] <- coef
      qscale[i]<-scale
      qwt[, i] <- w
      qfit[, i] <- temp$fitted.values
      qres[,i] <- resid

      tmp.res.mq<-qres[,i]/qscale[i]
      Epsi2<-(sum((qwt[,i]*tmp.res.mq)^2)/(nrow(x)-ncol(x)))
      Epsi<-(1/qscale[i])*(sum(2*(q[i]*(0<=tmp.res.mq & tmp.res.mq<= k)+(1-q[i])*(-k <=tmp.res.mq & tmp.res.mq<0)))/nrow(x))
      qvar.matrix[,,i]<- (((Epsi2)/Epsi^2)*solve(t(x)%*%x))

    }
    list(fitted.values = qfit, residuals = qres, q.values = q, q.weights = qwt, coef= qest,qscale=qscale,var.beta=qvar.matrix)
  }


  n<-length(ys)
  ni<-table(area.s)
  m<-length(ni)
  p<-ncol(xs)
  area_id<-sort(unique(area.s))
  MQ.Naive<-NULL
  MQ.CD<-NULL
  DIRECT<-tapply(ys,area.s,mean)
  #step 1: MQ function to get the q_ij
  mod<-QRLM(y=ys,x=xs,q=Q,maxit=maxit,acc=tol,k=k)
  qo<-matrix(c(gridfitinter(ys,mod$fitted.values,mod$q.values)),nrow=n,ncol=1)
  qmat<-matrix(c(qo,area.s),nrow=n,ncol=2)
  qmean<-tapply(qmat[,1],qmat[,2],mean)
  ob1<-QRLM(y=ys,x=xs,q=qmean,maxit=maxit,acc=tol,k=k)

  #Step 2: small area estimation
  for (i in 1:m)
  {
    MQ.Naive[i]<-(1/Ni[i])*(sum(ys[area.s==area_id[i]])+(Ni[i]-ni[i])*(barXr[i,]%*%ob1$coef[,i]))
    MQ.CD[i]<-DIRECT[i]+(XMean[i,-c(1)]-XsMean[i,-c(1)])*ob1$coef[2:p,i]
  }

  list(MQ=MQ.Naive, MQCD=MQ.CD)

}
