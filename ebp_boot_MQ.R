#' ebp_boot_MQ
#'
#' @title Bootstrap method for estimating the rmse of the MQ-based predictor.
#'
#' @description This function produces an estimation of the rmse of the MQ-based predictor using a bootstrap approach according the varying regression coefficients and variance components sigma2ei.
#'
#' @usage ebp_boot_MQ<-function(xs,barX,area.s,betaq,sigma2u,sigma2e,Ni,ni,B=50,numCores=2)
#'
#' @param xs a matrix or data frame containing the explanatory variables including the intercept.
#' @param barX m*(p+1) matrix that contains the population means of each of the p auxiliary variables for the m domains. The first column is the intercept.
#' @param area.s n*1 vector (same size as y in formula) with domain codes.
#' @param betaq a named vector of estimated coefficients obtained by ebp_mq function.
#' @param sigma2u the estimated sigma2u obtained by ebp_mq function.
#' @param sigma2e a named vector of estimated sigma2ei obtained by ebp_mq function.
#' @param Ni vector of domain population sizes.
#' @param ni vector of domain sample sizes.
#' @param B number of bootstrap iterations.
#' @param numCores number of cores used for paralleling the bootstrap procedure.
#' @param k_b tuning constant used for Huber influence function for the estimation of the coefficients.
#'
#' @import
#' doSNOW
#' tcltk
#' MASS
#' foreach
#' iterators
#'
#' @details In this function the results of ebp_nerhd are used for creating the bootstrap population.
#'
#' @return Return a list with the following elements:
#' \item{RMSE_MQ_BOOT}{the bootstrap RMSE estimates for each small area}
#'
#' @references Lahiri, P. and Salvati, N. (2022). A Nested Error Regression Model with High Dimensional Parameter for Small Area Estimation. arXiv:2201.10235, https://doi.org/10.48550/arXiv.2201.10235.
#' @references Chambers, R. and Tzavidis, N. (2006) M-quantile models for small area estimation. Biometrika, 93, 255-268.
#'
#' @author Nicola Salvati
#'
#' @seealso ebp_boot
#' @seealso ebp_boot_bhf
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
#'
#' P=sort(c(seq(0.006,0.99,0.045),0.5,0.994,0.01,0.02,0.96,0.98))
#' XMean<-cbind(1,tapply(pop$x,pop$area,mean))
#' dta<-data.frame(ys=y.s,xs=x.s,area.s=regioncode.s)
#' est<-ebp_nerhd(ys~1+xs,random=~1|area.s,est.method="global",barX=XMean,Ni=Ni,Q=P,method="LBP",tol=1e-06,maxit=100,k_b=1.345,k_sigma_u=1.345,k_sigma_e=100,fig="TRUE",uij.error.test="BHF",sigma2e.test="Raudenbush_Bryk",data=dta)
#' boot<-ebp_boot_MQ(xs=cbind(1,x.s),barX=XMean,area.s=regioncode.s,betaq=est$coefficients,sigma2u=est$est.sigma2u,sigma2e=est$est.sigma2e,Ni=Ni,ni=ni,B=10,numCores=2,k_b=1.345)
#'
#' @export
ebp_boot_MQ<-
  function(xs,barX,area.s,betaq,sigma2u,sigma2e,Ni,ni,B=50,numCores=2,k_b=100)
  {#start function
    #xs: design matrix for the sample with the intercept
    #betaq: it is the matrix of regression coefficients
    #barX: the average value of X at population level
    #Ni: population size in each area
    #ni: sample size in each area
    #sigma2u: vector of variance components
    #sigma2e: vector of variance components
    #B: number of bootstrap iterations
    #k_b=100: the tuning constant for betas
    #area.s: small area id
    #numCores: number of cores to use in the parallel
    require(doSNOW)
    require(tcltk)
    require(foreach)
    require(iterators)

    n<-sum(ni)
    m<-length(ni)
    p<-ncol(xs)
    area_id<-sort(unique(area.s))

    fun.boot<-function(h){
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
        #theta <- 2 * pnorm(k2) - 1
        #gamma <- theta + k2^2 * (1 - theta) - 2 * k2 * dnorm(k2)

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

      assign("xs",xs,pos=1)
      assign("barX",barX,pos=1)
      assign("betaq",betaq,pos=1)
      assign("sigma2u",sigma2u,pos=1)
      assign("sigma2e",sigma2e,pos=1)
      assign("Ni",Ni,pos=1)
      assign("ni",ni,pos=1)
      Yibar.boot<-NULL
      MQ.boot<-NULL

      u.boot<-NULL
      e.boot<-NULL
      for (i in 1:m)u.boot<-c(u.boot,rnorm(1,0,sqrt(sigma2u)))
      for (i in 1:m)e.boot <-c(e.boot,rnorm(ni[i], 0, sqrt(sigma2e[i])))
      y.boot<-NULL
      for (i in 1:m)y.boot=c(y.boot,xs[area.s==area_id[i],]%*%betaq[i,]+u.boot[i]+e.boot[area.s==area_id[i]])
      for (i in 1:m){
        Yibar.boot[i]<-barX[i,]%*%betaq[i,]+u.boot[i]
      }
      #Estimation procedure
      x.s.boot=xs
      y.s.boot=y.boot
      #EBLUP BHF

      mod.mq<-QRLM(x=x.s.boot,y=y.s.boot,q=sort(c(seq(0.006,0.99,0.045),0.5,0.994,0.01,0.02,0.96,0.98)),k=k_b)
      qo<-matrix(c(gridfitinter(y.s.boot,mod.mq$fitted.values,mod.mq$q.values)),nrow=n,ncol=1)
      qmat<-matrix(c(qo,area.s),nrow=n,ncol=2)
      mqo<-aggregate(qmat[,1],list(d2=qmat[,2]),mean)[,2]
      mod1.mq<-QRLM(x=x.s.boot,y=y.s.boot,q=mqo,k=k_b)

      for (i in 1:m){
        MQ.boot[i]<-barX[i,]%*%mod1.mq$coef[,i]
      }

      list(Yibar.boot=Yibar.boot,MQ.boot=MQ.boot)
    }#end fun.boot

    numCores <- numCores
    cl <- makeSOCKcluster(numCores)
    registerDoSNOW(cl)
    pb <- txtProgressBar(max=B, style=3)
    progress <- function(nnn) setTxtProgressBar(pb, nnn)
    opts <- list(progress=progress)
    out = foreach(h =1:B, .options.snow=opts) %dopar% {
      #sink('log.txt', append=TRUE)
      Sys.sleep(1)
      set.seed(h)
      fun.boot(h)
    }
    close(pb)
    stopCluster(cl)

    Yibar.boot<-matrix(0,B,m)
    MQ.boot<-matrix(0,B,m)

    for (ii in 1:B)
    {
      Yibar.boot[ii,]<-out[[ii]]$Yibar.boot
      MQ.boot[ii,]<-out[[ii]]$MQ.boot
    }
    a_phi.ebp<-NULL
    a_phi.ebp<-(sqrt(apply((MQ.boot-Yibar.boot)^2,2,mean)))
    list(RMSE_MQ_BOOT=a_phi.ebp)
    #end function
  }
