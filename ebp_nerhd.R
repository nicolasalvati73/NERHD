#' ebp_nerhd
#'
#' @title Empirical Best Predictor under NERHD model
#'
#' @description This functions is used to obtain the EBP fitting a varying regression coefficients and random error variance models based on NERHD approach for small area estimation.
#'
#' @usage ebp_nerhd(formula,random=~1 barX, est.method = c("local", "global"), Ni, Q = sort(c(seq(0.006, 0.99, 0.045), 0.5, 0.994, 0.01, 0.02, 0.96, 0.98)), method = c("LBP", "Naive"), tol = 1e-06, maxit = 100, k_b = 100, k_sigma_u = 100, k_sigma_e = 100, fig = "TRUE", uij.error.test = c("BHF", "CS"), sigma2e.test = c("Raudenbush_Bryk", "Bartlett"))
#'
#' @param formula an object of class formula of the form y ~ x1 + x2 + ... + xp for fixed coefficients.
#' @param random the random part of the model defining the groups.
#' @param barX m*(p+1) matrix that contains the population means of each of the p auxiliary variables for the m domains. The first column is the intercept.
#' @param area.s n*1 vector (same size as y in formula) with domain codes.
#' @param est.method estimation method for sigma2ei. It could be local or global.
#' @param Ni vector of domain population sizes.
#' @param Q the quantile(s) to be estimated, this is generally a number strictly between 0 and 1.
#' @param method a character string. It is the method for computiong the M-quantile coefficients at area level. It could be LBP or Naive. LBP uses Empirical Linear Best (ELB) predictor, while Naive is the traditional method suggested by Chambers & Tzavidis (2006) for computing the M-quantile coefficients at area level by averaging the unit level M-quantile coefficients. Default is LBP.
#' @param tol the accuracy for the stopping criterion.
#' @param maxit the limit on the number of algorithm iterations.
#' @param k_b tuning constant used for Huber influence function for the estimation of the coefficients.
#' @param k_sigma_u tuning constant used for Huber influence function for the estimation of the sigma2u.
#' @param k_sigma_e tuning constant used for Huber influence function for the estimation of the sigma2ei.
#' @param fig logical. If TRUE the qqplots are returned.
#' @param uij.error.test statistic to test the normality of u_ij errors: if BHF the test is run as in in Battese et al. (1988) otherwise if it is CS the test is carry out as in in Calvin and Sedransk (1991).
#' @param sigma2e.test statistic to test equality of the sigma2ei between areas: it could be  Raudenbush and Bryk's test (2002) or Bartlett's test (Snedecor, George W. and Cochran, William G. (1989)).
#' @param data a data frame containing the y, the x variables and the grouping variables
#'
#' @import
#' sae
#' Matrix
#' lme4
#' MASS
#' statmod
#'
#' @details In this function the Huber influence function is used by default and the MAD is used by default for estimating the scale.
#'
#' @return Return a list with the following elements:
#' \item{coefficients}{a named vector of coefficients}
#' \item{EBP}{a named vector of ebp predictor at small area level}
#' \item{est.sigma2u}{the value of the estimated sigma2u}
#' \item{est.sigma2e}{a named vector of estimated sigma2ei}
#' \item{Bi}{a named vector of shrinkage factor}
#' \item{Qi}{the estimated quantile(s)}
#' \item{RMSE.EBP}{a named vector of estimated sqrt of g1}
#' \item{betas.test}{the results of the test on the non-stationarity of coefficients}
#' \item{test.sigma2e}{the result of the test on the non-stationarity of variance components sigma2ei}
#' \item{test.shapiro}{the result of the shapiro test on the u_ij residuals}
#'
#' @references Battese, G., Harter, R. and Fuller, W. (1988) An error-components model for prediction of county crop areas using survey and satellite data. J. Am. Statist. Ass., 83, 28-36.
#' @references Calvin, J. A. and Sedransk, J. (1991). Bayesian and frequentist predictive inference for the patternsof care studies.Journal of the American Statistical Association, 86, 36-48.
#' @references Chambers, R. and Tzavidis, N. (2006) M-quantile models for small area estimation. Biometrika, 93, 255-268.
#' @references Lahiri, P. and Salvati, N. (2022). A Nested Error Regression Model with High Dimensional Parameter for Small Area Estimation. arXiv:2201.10235, https://doi.org/10.48550/arXiv.2201.10235
#' @references Raudenbush, S.W., & Bryk, A.S. (2002). Hierarchical Linear Models: Applications and data analysis methods. (2nd ed.). Thousand Oaks, CA: Sage Publications.
#' @references Snedecor, George W. and Cochran, William G. (1989), Statistical Methods, Eighth Edition, Iowa State University Press.
#'
#' @author Nicola Salvati
#'
#' @seealso ebp_mcjack
#' @seealso ebp_boot
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
#' dta=data.frame(ys=y.s,xs=x.s,area=regioncode.s)
#' est<-ebp_nerhd(formula=ys~1+xs,random=~1|area,est.method="global",barX=XMean,Ni=Ni,Q=P,method="LBP",tol=1e-06,maxit=100,k_b=1.345,k_sigma_u=1.345,k_sigma_e=100,fig="TRUE",uij.error.test="BHF",sigma2e.test="Raudenbush_Bryk",data=dta)
#'
#'
#' @export
ebp_nerhd <-function(formula,random=~1,barX,est.method=c("local","global"),Ni,Q=sort(c(seq(0.006,0.99,0.045),0.5,0.994,0.01,0.02,0.96,0.98)),method=c("LBP","Naive"),tol=1e-06,maxit=100,k_b=100,k_sigma_u=100,k_sigma_e=100,fig="TRUE",uij.error.test=c("BHF","CS"),sigma2e.test=c("Raudenbush_Bryk","Bartlett"),data)
{ #barX: the average value of X at population level
  #Ni: population size in each area
  #Q: the quantiles for fitting the MQ
  #tol=1e-06
  #maxit=100
  #est.method: it is the estimation method for sigma2e. It could be local or global
  #k_b=100: the tuning constant for betas
  #k_sigma_u=100: the tuning constant for sigma2u
  #k_sigma_e=100: the tuning constant for sigma2e
  #method: Naive, traditional method to estimate Qi. LBP is the new proposal based on linear best predictor
  #fig: show the plots if it is TRUE (default)
  #uij.error.test: BHF in Battese et al. (1988); CS in Calvin and Sedransk (1991)
  #sigma2e.test: test on the difference between sigma2ei. Raudenbush_Bryk or Bartlett
  require(sae)
  require(Matrix)
  require(lme4)
  require(MASS)
  require(statmod)
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
  mqre_l<-function(qtl,y,x,group,tol=1e-06,maxit=100,k=100,k_sigma_u=100,k_sigma_e=100)
  {
    ###############functions################

    my.psi<-function(u,q,k){
      sm<-median(abs(u))/0.6745
      w <- psi.huber(u/sm,k)
      ww <- 2 * (1 - q) * w
      ww[u> 0] <- 2 * q * w[u > 0]
      w <- ww
      w*u
    }

    der.psi<-function(u,q,k){
      sm<-median(abs(u))/0.6745
      u<-u/sm
      der.psi <- ifelse(abs(u) <= k, 1, 0)
      ww <- 2 * (1 - q) * der.psi
      ww[u> 0] <- 2 * q * der.psi[u > 0]
      w <- ww
      w
    }
    ########################starting data###################
    n=length(y)
    ni=table(group)
    m=length(ni)
    areanumber=m
    p=ncol(x)
    gg=unique(group)
    z<-model.matrix(y~as.factor(group)-1)

    #STEP 1 (Definition of the starting values)
    w=cbind(x,z)
    fit.H=lm(y~w)
    e.1.H=residuals(fit.H)
    sigma.e.hat.H=sum(e.1.H^2)/(n-(qr(w)$rank))

    xx=x
    fit2.H=lm(y~xx)
    e.2.H=residuals(fit2.H)
    A.H= sum(e.2.H^2)
    xx=as.matrix(xx)
    B.H=sigma.e.hat.H*(n-qr(xx)$rank)
    o.H=diag(n)-(xx%*%(ginv(t(xx)%*%xx)%*%t(xx)))
    C.H=sum(diag(t(z)%*%o.H%*%z))
    sigma.v.hat.H=(A.H-B.H)/C.H
    sigmasq0= sigma.e.hat.H #initial values of sigma_sq_e#
    sigmasq0.v=sigma.v.hat.H #initial values sigma_sq_V#

    if(dim(x)[2]==1){ls <- lm(y ~ x[,1]-1)}else{
      ls <- lm(y ~ x[,-1])}
    beta1<-c(ls$coefficients)
    estsigma2u<-sigmasq0.v
    estsigma2e<-sigmasq0
    sigma1G.old<-NULL
    sigma1E.old<-NULL
    sigma1E.old<-rep(sigmasq0,areanumber)
    sigma1G.old<-c(sigmasq0.v)
    sigma1G.est<-NULL
    sigma1E.est<-NULL

    iter<-0
    ZZ<-z%*%t(z)
    beta.old<-matrix(beta1,areanumber,p,byrow=T)
    beta.est<-matrix(0,areanumber,p)
    var.beta<-array(rep(0,ncol(x)*ncol(x)),dim=c(ncol(x),ncol(x),length(qtl)))

    #STEP 2: estimation of betas and sigma2
    while(TRUE)
    { residuals<-NULL
    for (i in 1:length(qtl))
    {
      #STEP 1 Estimation of beta
      xbeta <- c(x %*% beta.old[i,])
      R <- diag(rep(sigma1E.old[i], n),n,n)
      G <- diag(rep(sigma1G.old, areanumber),areanumber,areanumber)
      V <- R + z %*% G %*% t(z)
      V.inv<-chol2inv(chol(V))
      U <- diag(diag(V),n,n)
      U.inv <- chol2inv(chol(U))
      ystar<-c(as.numeric(sqrt(U.inv) %*% y))
      xstar<-matrix(as.numeric(sqrt(U.inv) %*% x),n,p)
      mod<-QRLM(y=ystar,x=xstar,q=qtl[i],maxit=maxit,acc=tol,k=k)
      beta.est[i,]<-mod$coef
      var.beta[,,i]<-mod$var.beta
      r <- c(sqrt(U.inv) %*% (y - xbeta))
      ri<-as.vector(r)[group==i]
      #mod.rlmerE<-rlmer(r~1+(1|group),rho.sigma.e = psi2propII(smoothPsi, k = k_sigma),
      #   rho.sigma.b = psi2propII(smoothPsi, k = k_sigma))
      psi.r <- my.psi(as.vector(ri), k = k_sigma_e, q=qtl[i])
      V.inv.i<-V.inv[group==i,group==i]
      U.i<-U[group==i,group==i]
      ZZ.i<-ZZ[group==i,group==i]
      if(ni[i]>1){qq <- V.inv.i %*% sqrt(U.i) %*% psi.r
      KK5<-sqrt(U.i) %*% V.inv.i
      KK0<-t(psi.r)%*% KK5
      A1=KK0%*%qq
      A2=KK0%*%ZZ.i%*%qq
      A <- matrix(c(A1[1],A2[1]), nrow = 2, ncol=1)}

      if(ni[i]==1){qq <- V.inv.i * sqrt(U.i) * psi.r
      KK5<-sqrt(U.i)*V.inv.i
      KK0<-t(psi.r)*KK5
      A1=KK0*qq
      A2=KK0*ZZ.i*qq
      A <- matrix(c(A1[1],A2[1]), nrow = 2, ncol=1)}

      const<- 4*k_sigma_e^2*(1-pnorm(k_sigma_e))*((1-qtl[i])^2+qtl[i]^2) - 4*k_sigma_e*dnorm(k_sigma_e)*((1-qtl[i])^2+qtl[i]^2) + 4*(1-qtl[i])^2*(pnorm(0)-(1-pnorm(k_sigma_e))) + 4*qtl[i]^2*(pnorm(k_sigma_e)-pnorm(0))
      K2_value<- const

      if(ni[i]>1){
        K2<-diag(K2_value, ni[i])
        KK1<-K2%*%V.inv.i
        KK2<-KK1%*%V.inv.i
        KK3<-KK1%*%ZZ.i%*%V.inv.i
        t1=sum(diag(KK2))
        t2=sum(diag(KK2%*%ZZ.i))
        t3=sum(diag(KK3))
        t4=sum(diag(KK3%*%ZZ.i))
      }

      if(ni[i]==1){
        K2<-diag(K2_value, ni[i])
        KK1<-K2*V.inv.i
        KK2<-KK1*V.inv.i
        KK3<-KK1*ZZ.i*V.inv.i
        t1=KK2
        t2=KK2*ZZ.i
        t3=KK3
        t4=KK3*ZZ.i}

      T1<- matrix(c(t1,t3,t2,t4), nrow = 2, ncol=2)
      sigma=chol2inv(chol(T1))%*%A
      sigma1E.est[i]<-sigma[1,1]
      if(p==1)residuals<-c(residuals,(y[group==i]-x[group==i,1]*beta.est[i,1]))
      if(p==2)residuals<-c(residuals,(y[group==i]-x[group==i,2]*beta.est[i,2]))
      if(p>2)residuals<-c(residuals,(y[group==i]-x[group==i,-1]%*%beta.est[i,-1]))
    }#end beta estimation

    residuals<-residuals-mean(residuals)
    R <- diag(rep(sigma1E.old, ni),n,n)
    G <- diag(rep(sigma1G.old, areanumber),areanumber,areanumber)
    V <- R + z %*% G %*% t(z)
    V.inv<-chol2inv(chol(V))
    U <- diag(diag(V))
    U.inv <- chol2inv(chol(U))
    r <- c(sqrt(U.inv) %*% residuals)
    psi.r <- my.psi(as.vector(r), k = k_sigma_u, q=0.5)
    qq <- V.inv %*% sqrt(U) %*% psi.r
    KK5<-sqrt(U) %*% V.inv
    KK0<-t(psi.r)%*% KK5
    A1=KK0%*%qq
    A2=KK0%*%ZZ%*%qq
    A <- matrix(c(A1[1],A2[1]), nrow = 2, ncol=1)
    const<- 4*k_sigma_u^2*(1-pnorm(k_sigma_u))*((1-0.5)^2+0.5^2) - 4*k_sigma_u*dnorm(k_sigma_u)*((1-0.5)^2+0.5^2) + 4*(1-0.5)^2*(pnorm(0)-(1-pnorm(k_sigma_u))) + 4*0.5^2*(pnorm(k_sigma_u)-pnorm(0))
    K2_value<- const
    K2<-diag(K2_value, n)
    KK1<-K2%*%V.inv
    KK2<-KK1%*%V.inv
    KK3<-KK1%*%ZZ%*%V.inv
    t1=sum(diag(KK2))
    t2=sum(diag(KK2%*%ZZ))
    t3=sum(diag(KK3))
    t4=sum(diag(KK3%*%ZZ))
    T1<- matrix(c(t1,t3,t2,t4), nrow = 2, ncol=2)
    sigma=chol2inv(chol(T1))%*%A
    sigma1G.est<-sigma[2,1]
    if(p>1)B0<-mean(beta.est[,1])
    if(p==1)B0<-0

    #STEP 3: computation variance components sigma^2_gamma

    if(p>1)fehler <- sum((beta.old[,-1] -beta.est[,-1])^2)+ sum((sigma1G.old -sigma1G.est)^2+sum((sigma1E.old -sigma1E.est)^2))
    if(p==1) fehler <- sum((beta.old -beta.est)^2)+ sum((sigma1G.old -sigma1G.est)^2+sum((sigma1E.old -sigma1E.est)^2))
    ifelse (iter < maxit, ifelse(fehler > tol, {iter=iter+1}, {break}), {break})
    beta.old <- beta.est
    sigma1G.old<-sigma1G.est
    sigma1E.old<-sigma1E.est
    } #end while, estimation procedure

    if(sigma1G.est<0){lme_fit<-lmer(y~-1+x+(1|group))
    sigma1G.est<-summary(lme_fit)$varcor$group[1]}
    if(p>1)coefficients =cbind(rep(B0,areanumber),beta.est[,-1])
    if(p==1)coefficients =beta.est
    list(coefficients =coefficients,sigma2u=(sigma1G.est),sigma2e=(sigma1E.est),varBeta=var.beta,quantile=qtl)

  }

  mqre_g<-function(qtl,y,x,group,tol=1e-06,maxit=100,k=100,k_sigma_u=100,k_sigma_e=100)
  {
    ###############functions################

    my.psi<-function(u,q,k){
      sm<-median(abs(u))/0.6745
      w <- psi.huber(u/sm,k)
      ww <- 2 * (1 - q) * w
      ww[u> 0] <- 2 * q * w[u > 0]
      w <- ww
      w*u
    }

    der.psi<-function(u,q,k){
      sm<-median(abs(u))/0.6745
      u<-u/sm
      der.psi <- ifelse(abs(u) <= k, 1, 0)
      ww <- 2 * (1 - q) * der.psi
      ww[u> 0] <- 2 * q * der.psi[u > 0]
      w <- ww
      w
    }
    ########################starting data###################
    n=length(y)
    ni=table(group)
    m=length(ni)
    areanumber=m
    p=ncol(x)
    gg=unique(group)
    z<-model.matrix(y~as.factor(group)-1)

    #STEP 1 (Definition of the starting values)
    w=cbind(x,z)
    fit.H=lm(y~w)
    e.1.H=residuals(fit.H)
    sigma.e.hat.H=sum(e.1.H^2)/(n-(qr(w)$rank))

    xx=x
    fit2.H=lm(y~xx)
    e.2.H=residuals(fit2.H)
    A.H= sum(e.2.H^2)
    xx=as.matrix(xx)
    B.H=sigma.e.hat.H*(n-qr(xx)$rank)
    o.H=diag(n)-(xx%*%(ginv(t(xx)%*%xx)%*%t(xx)))
    C.H=sum(diag(t(z)%*%o.H%*%z))
    sigma.v.hat.H=(A.H-B.H)/C.H
    sigmasq0= sigma.e.hat.H #initial values of sigma_sq_e#
    sigmasq0.v=sigma.v.hat.H #initial values sigma_sq_V#

    if(dim(x)[2]==1){ls <- lm(y ~ x[,1]-1)}else{
      ls <- lm(y ~ x[,-1])}
    beta1<-c(ls$coefficients)
    estsigma2u<-sigmasq0.v
    estsigma2e<-sigmasq0
    sigma1G.old<-NULL
    sigma1E.old<-NULL
    sigma1E.old<-rep(sigmasq0,areanumber)
    sigma1G.old<-c(sigmasq0.v)
    sigma1G.est<-NULL
    sigma1E.est<-NULL

    iter<-0
    ZZ<-z%*%t(z)
    beta.old<-matrix(beta1,areanumber,p,byrow=T)
    beta.est<-matrix(0,areanumber,p)
    var.beta<-array(rep(0,ncol(x)*ncol(x)),dim=c(ncol(x),ncol(x),length(qtl)))

    #STEP 2: estimation of betas and sigma2
    while(TRUE)
    { residuals<-NULL
    for (i in 1:length(qtl))
    {
      #STEP 1 Estimation of beta
      xbeta <- c(x %*% beta.old[i,])
      R <- diag(rep(sigma1E.old[i], n),n,n)
      G <- diag(rep(sigma1G.old, areanumber),areanumber,areanumber)
      V <- R + z %*% G %*% t(z)
      V.inv<-chol2inv(chol(V))
      U <- diag(diag(V))
      U.inv <- chol2inv(chol(U))
      ystar<-c(as.numeric(sqrt(U.inv) %*% y))
      xstar<-matrix(as.numeric(sqrt(U.inv) %*% x),n,p)
      mod<-QRLM(y=ystar,x=xstar,q=qtl[i],maxit=maxit,acc=tol,k=k)
      beta.est[i,]<-mod$coef
      var.beta[,,i]<-mod$var.beta
      r <- c(sqrt(U.inv) %*% (y - xbeta))
      #mod.rlmerE<-rlmer(r~1+(1|group),rho.sigma.e = psi2propII(smoothPsi, k = k_sigma),
      #   rho.sigma.b = psi2propII(smoothPsi, k = k_sigma))
      psi.r <- my.psi(as.vector(r), k = k_sigma_e, q=qtl[i])
      qq <- V.inv %*% sqrt(U) %*% psi.r
      KK5<-sqrt(U) %*% V.inv
      KK0<-t(psi.r)%*% KK5
      A1=KK0%*%qq
      A2=KK0%*%ZZ%*%qq
      A <- matrix(c(A1[1],A2[1]), nrow = 2, ncol=1)
      const<- 4*k_sigma_e^2*(1-pnorm(k_sigma_e))*((1-qtl[i])^2+qtl[i]^2) - 4*k_sigma_e*dnorm(k_sigma_e)*((1-qtl[i])^2+qtl[i]^2) + 4*(1-qtl[i])^2*(pnorm(0)-(1-pnorm(k_sigma_e))) + 4*qtl[i]^2*(pnorm(k_sigma_e)-pnorm(0))
      K2_value<- const
      K2<-diag(K2_value, n)
      KK1<-K2%*%V.inv
      KK2<-KK1%*%V.inv
      KK3<-KK1%*%ZZ%*%V.inv
      t1=sum(diag(KK2))
      t2=sum(diag(KK2%*%ZZ))
      t3=sum(diag(KK3))
      t4=sum(diag(KK3%*%ZZ))
      T1<- matrix(c(t1,t3,t2,t4), nrow = 2, ncol=2)
      sigma=chol2inv(chol(T1))%*%A
      sigma1E.est[i]<-sigma[1,1]
      if(p==1)residuals<-c(residuals,(y[group==i]-x[group==i,1]*beta.est[i,1]))
      if(p==2)residuals<-c(residuals,(y[group==i]-x[group==i,2]*beta.est[i,2]))
      if(p>2)residuals<-c(residuals,(y[group==i]-x[group==i,-1]%*%beta.est[i,-1]))
    }#end beta estimation

    #STEP 3: computation variance components sigma^2_gamma
    residuals<-residuals-mean(residuals)
    R <- diag(rep(sigma1E.old, ni),n,n)
    G <- diag(rep(sigma1G.old, areanumber),areanumber,areanumber)
    V <- R + z %*% G %*% t(z)
    V.inv<-chol2inv(chol(V))
    U <- diag(diag(V))
    U.inv <- chol2inv(chol(U))
    r <- c(sqrt(U.inv) %*% residuals)
    psi.r <- my.psi(as.vector(r), k = k_sigma_u, q=0.5)
    qq <- V.inv %*% sqrt(U) %*% psi.r
    KK5<-sqrt(U) %*% V.inv
    KK0<-t(psi.r)%*% KK5
    A1=KK0%*%qq
    A2=KK0%*%ZZ%*%qq
    A <- matrix(c(A1[1],A2[1]), nrow = 2, ncol=1)
    const<- 4*k_sigma_u^2*(1-pnorm(k_sigma_u))*((1-0.5)^2+0.5^2) - 4*k_sigma_u*dnorm(k_sigma_u)*((1-0.5)^2+0.5^2) + 4*(1-0.5)^2*(pnorm(0)-(1-pnorm(k_sigma_u))) + 4*0.5^2*(pnorm(k_sigma_u)-pnorm(0))
    K2_value<- const
    K2<-diag(K2_value, n)
    KK1<-K2%*%V.inv
    KK2<-KK1%*%V.inv
    KK3<-KK1%*%ZZ%*%V.inv
    t1=sum(diag(KK2))
    t2=sum(diag(KK2%*%ZZ))
    t3=sum(diag(KK3))
    t4=sum(diag(KK3%*%ZZ))
    T1<- matrix(c(t1,t3,t2,t4), nrow = 2, ncol=2)
    sigma=chol2inv(chol(T1))%*%A
    sigma1G.est<-sigma[2,1]
    if(p==1)B0<-0
    if(p>1)B0<-mean(beta.est[,1])

    if(p>1)fehler <- sum((beta.old[,-1] -beta.est[,-1])^2)+ sum((sigma1G.old -sigma1G.est)^2+sum((sigma1E.old -sigma1E.est)^2))
    if(p==1)fehler <- sum((beta.old[,1] -beta.est[,1])^2)+ sum((sigma1G.old -sigma1G.est)^2+sum((sigma1E.old -sigma1E.est)^2))
    ifelse (iter < maxit, ifelse(fehler > tol, {iter=iter+1}, {break}), {break})
    beta.old <- beta.est
    sigma1G.old<-sigma1G.est
    sigma1E.old<-sigma1E.est
    } #end while, estimation procedure

    if(sigma1G.est<0){lme_fit<-lmer(y~-1+x+(1|group))
    sigma1G.est<-summary(lme_fit)$varcor$group[1]}
    if(p>1)coefficients =cbind(rep(B0,areanumber),beta.est[,-1])
    if(p==1)coefficients =beta.est
    list(coefficients =coefficients,sigma2u=(sigma1G.est),sigma2e=(sigma1E.est),varBeta=var.beta,quantile=qtl)

  }
  if (!is.data.frame(data))
    stop("`data' must be a data frame")
  model.f<-model.frame(formula,data=data)
  ys<-as.numeric(model.response(model.f))
  xs<-model.matrix(formula, data)[, , drop = FALSE]
  rform <- random
  mform <- strsplit(as.character(random)[2], "\\|")[[1]]
  mform <- gsub(" ", "", mform)
  mform1 <- mform[1]
  mform2 <- mform[2]
  if (!(mform2 %in% names(data))) {
    stop("The specified clustering variable is not contained in the data frame.")
  }
  formula.group=formula(paste("~",mform2))
  area.s=model.matrix(formula.group, data)[, -1, drop = FALSE]
  n<-length(ys)
  ni<-table(area.s)
  m<-length(ni)
  p<-ncol(xs)
  area_id<-sort(unique(area.s))
  #step 1: MQ function to get the q_ij
  mod<-QRLM(y=ys,x=xs,q=Q,maxit=maxit,acc=tol,k=k_b)
  qo<-matrix(c(gridfitinter(ys,mod$fitted.values,mod$q.values)),nrow=n,ncol=1)
  qmat<-matrix(c(qo,area.s),nrow=n,ncol=2)
  #step 2: LBP to get the Q_i
  if(method=="LBP")
  {tmp.Xmean<-rep(1,m)
  tmp.Popn<-Ni
  Xmean <- data.frame(area_id, tmp.Xmean)
  Popn  <- data.frame(area_id, tmp.Popn)
  tmp.dta<-data.frame(qmat[,1],area.s)
  names(tmp.dta)<-c("qij","area.s")
  LBP<-eblupBHF(formula=qij~1, dom=area.s, meanxpop=Xmean, popnsize=Popn,method = "REML", data=tmp.dta)
  Qi<-LBP$eblup$eblup
  Qi<-ifelse(Qi>1,0.999,Qi)
  Qi<-ifelse(Qi<0,0.001,Qi)
  }
  if(method=="Naive"){Qi<-tapply(qmat[,1],qmat[,2],mean)}
  #step 3: estimate Beta_Q_i, sigma^2_v_i and sigma^2_e_i
  if(est.method=="local")out1<-mqre_l(qtl=Qi,y=ys,x=xs,group=area.s,tol=tol,maxit=maxit,k=k_b,k_sigma_u=k_sigma_u,k_sigma_e=k_sigma_e)
  if(est.method=="global")out1<-mqre_g(qtl=Qi,y=ys,x=xs,group=area.s,tol=tol,maxit=maxit,k=k_b,k_sigma_u=k_sigma_u,k_sigma_e=k_sigma_e)
  #step 4: Computation MQ_BP and sqrt(g1)
  MQ_BP<-NULL
  MSE.MQ_BP<-NULL
  g1<-rep(0,m)
  Bi<-NULL
  est.sigma2u<-out1$sigma2u
  est.sigma2e<-NULL
  ranef<-NULL
  coefficients<-matrix(0,m,p)
  for (i in 1:m)
  {
    if(ni[i]>1 & p>1)barx=apply(xs[area.s==area_id[i],],2,mean)
    if(ni[i]>1 & p==1)barx=mean(xs[area.s==area_id[i]])
    if(ni[i]==1)barx=xs[area.s==area_id[i],]
    if(ni[i]==1 & p==1)barx=xs[area.s==area_id[i]]
    bary=mean(ys[area.s==area_id[i]])
    Bi[i]<-out1$sigma2e[i]/ni[i]/(out1$sigma2e[i]/ni[i]+out1$sigma2u)
    est.sigma2e[i]<-out1$sigma2e[i]
    coefficients[i,]<-out1$coefficients[i,]
    if(p>1)MQ_BP[i]<-barX[i,]%*%out1$coefficients[i,]+(1-Bi[i])*(bary-barx%*%out1$coefficients[i,])
    if(p==1)MQ_BP[i]<-barX[i,2]%*%out1$coefficients[i,]+(1-Bi[i])*(bary-barx%*%out1$coefficients[i,])
    ranef[i]<-(1-Bi[i])*(bary-barx%*%out1$coefficients[i,])
    g1[i]<-as.numeric((out1$sigma2u*out1$sigma2e[i]/ni[i])/(out1$sigma2e[i]/ni[i]+out1$sigma2u))
    MSE.MQ_BP[i]<-((1-ni[i]/Ni[i])^2)*g1[i]
  }

  #step 5: diagnostic analysis as in Battese et al 1988
  lme_fit<-lmer(ys~-1+xs+(1|area.s))
  lme.sigma2u<-summary(lme_fit)$varcor$area.s[1]
  lme.sigma2e<-summary(lme_fit)$sigma^2
  BetaFix<-matrix(c(summary(lme_fit)$coefficients[,1]),p,1)

  if(uij.error.test=="BHF")
  {
    #trasformed residuals
    hat_uij<-NULL
    hat_uij1<-NULL
    for (i in 1:m)
    {
      bary=mean(ys[area.s==area_id[i]])
      if(ni[i]>1 & p>1)barx=apply(xs[area.s==area_id[i],],2,mean)
      if(ni[i]>1 & p==1)barx=mean(xs[area.s==area_id[i]])
      if(ni[i]==1 & p>1)barx=xs[area.s==area_id[i],]
      if(ni[i]==1 & p==1)barx=xs[area.s==area_id[i]]

      alphai<-1-sqrt(out1$sigma2e[i]/(out1$sigma2e[i]+ni[i]*out1$sigma2u))
      alphai1<-1-sqrt(lme.sigma2e/(lme.sigma2e+ni[i]*lme.sigma2u))
      tmp.hat_uij<-NULL
      tmp.hat_uij1<-NULL
      yi<-ys[area.s==area_id[i]]
      xsi<-xs[area.s==area_id[i],]

      for (j in 1:ni[i])
      {
        if(ni[i]==1){tmp.hat_uij1[j]<-(yi[j]-alphai1*bary)-(xsi-alphai1*barx)%*%BetaFix;
        tmp.hat_uij[j]<-(yi[j]-alphai*bary)-(xsi-alphai*barx)%*%out1$coefficients[i,]}
        if(ni[i]>1 & p>1){tmp.hat_uij1[j]<-(yi[j]-alphai1*bary)-(xsi[j,]-alphai1*barx)%*%BetaFix;
        tmp.hat_uij[j]<-(yi[j]-alphai*bary)-(xsi[j,]-alphai*barx)%*%out1$coefficients[i,]}
        if(ni[i]>1 & p==1){tmp.hat_uij1[j]<-(yi[j]-alphai1*bary)-(xsi[j]-alphai1*barx)%*%BetaFix;
        tmp.hat_uij[j]<-(yi[j]-alphai*bary)-(xsi[j]-alphai*barx)%*%out1$coefficients[i,]}
      }
      hat_uij<-c(hat_uij,tmp.hat_uij)
      hat_uij1<-c(hat_uij1,tmp.hat_uij1)
    }
 }


  if(uij.error.test=="CS")
  {
    hat_uij<-NULL
    hat_uij1<-NULL
    for (i in 1:m)
    {
      #Varying coef model
      Ri<-ifelse(ni[i]>1,diag(rep(out1$sigma2e[i], ni[i])),out1$sigma2e[i])

      Gi <-ifelse(ni[i]>1,matrix(out1$sigma2u,ni[i],ni[i]),out1$sigma2u)
      Vi <- Gi+Ri
      sqrt.Vi.inv<-sqrt(diag(chol2inv(chol(Vi))))
      yi<-ys[area.s==area_id[i]]
      xsi<-xs[area.s==area_id[i],]
      if(p>1)resi<-yi-xsi%*%out1$coefficients[i,]
      if(p==1)resi<-yi-xsi*out1$coefficients[i,]
      hat_uij<-c(hat_uij,c(sqrt.Vi.inv*resi))

      #lme model
      R <- ifelse(ni[i]>1,diag(rep(lme.sigma2e, ni[i])),lme.sigma2e)
      G <- ifelse(ni[i]>1,matrix(lme.sigma2u,ni[i],ni[i]),lme.sigma2u)
      V <- G+R
      sqrt.V.inv<-sqrt(diag(chol2inv(chol(V))))
      resi.lme<-yi-xsi%*%BetaFix
      hat_uij1<-c(hat_uij1,c(sqrt.V.inv*resi.lme))

    }

  }

  r<-(hat_uij-mean(hat_uij))/(sd(hat_uij))
  r1<-(hat_uij1-mean(hat_uij1))/(sd(hat_uij1))

  #standardized random effects according Lange Ryan
  std.gamma<-NULL
  weights.cdf<-NULL
  for (i in 1:m)
  {
    yi<-ys[area.s==area_id[i]]
    if(ni[i]>1 & p>1)barx=apply(xs[area.s==area_id[i],],2,mean)
    if(ni[i]>1 & p==1)barx=mean(xs[area.s==area_id[i]])
    if(ni[i]==1)barx=xs[area.s==area_id[i],]
    if(ni[i]==1 & p==1)barx=xs[area.s==area_id[i]]
    std.gamma[i]<-(mean(yi)-barx%*%out1$coefficients[i,])/sqrt(out1$sigma2e[i]/ni[i]+out1$sigma2u)
    weights.cdf[i]<-(ni[i]*out1$sigma2u^2)/(ni[i]*out1$sigma2u+out1$sigma2e[i])
  }

  tmp.std.gamma<-order(std.gamma)
  std.gamma<-std.gamma[tmp.std.gamma]
  weights.cdf<-weights.cdf[tmp.std.gamma]
  F.star<-cumsum(weights.cdf)/sum(weights.cdf)

  if(fig=="TRUE"){
    par(mfrow=c(1,2))
    qqnorm(r,ylab="standardized residuals",xlab="standard normal",main=expression('u'[ij]),cex=0.7,pch=3)
    grid(5,5)
    qqline(r, col=2,lty="dashed")
    qqplot(qnorm(F.star),std.gamma,ylab="standardized random effects",xlab="standard normal",main="random effects",cex=0.7,pch=3)
    grid(5,5)
    qqline(std.gamma, col=2,lty="dashed")
    #qqnorm(r1,ylab="residuals",xlab="standard normal",main="EBLUP",cex=0.7,pch=3)
    #grid(5,5)
    #qqline(r1, col=2,lty="dashed")

  }
  shapiro<-shapiro.test(hat_uij)
  shapiro.re<-shapiro.test(std.gamma)
  shapiro.lme<-shapiro.test(hat_uij1)

  #test coefficients J. A. Hausman Vol. 46, No. 6 (Nov., 1978), pp. 1251-1271
  if(p>1)VarBetaFix<-matrix(as.numeric(summary(lme_fit)$vcov),p,p,byrow=T)[-1,-1]
  if(p==1)VarBetaFix<-matrix(as.numeric(summary(lme_fit)$vcov),p,p,byrow=T)
  Fs<-NULL
  for (i in 1:m)
  {
    if(p>1)VarBeta<-out1$varBeta[-1,-1,i]
    if(p==1)VarBeta<-out1$varBeta[,,i]
    if(p>2)Fs[i]<-as.numeric(t(out1$coefficients[i,][-1]-BetaFix[-1,])%*%solve(VarBeta-VarBetaFix)%*%(out1$coefficients[i,][-1]-BetaFix[-1,]))
    if(p==2)Fs[i]<-as.numeric(t(out1$coefficients[i,][-1]-BetaFix[-1,])*(1/(VarBeta-VarBetaFix))*(out1$coefficients[i,][-1]-BetaFix[-1,]))
    if(p==1)Fs[i]<-as.numeric(t(out1$coefficients[i,]-BetaFix)*(1/(VarBeta-VarBetaFix))*(out1$coefficients[i,]-BetaFix))
    }

  p_value_F<-round(1-pchisq(sum(Fs),((p-1)*m)),4)
  betas.test<-list(statistic = c(Delta = sum(Fs)), p.value = p_value_F,
                   method = "Beta Hausman's test",data.name = "Regression Coefficients")
  class(betas.test) <- "htest"
  #Raudenbush and Bryk, 2002 test. Possible alternative Bartlett's test
  if(sigma2e.test=="Raudenbush_Bryk"){
    di<-NULL
    tmp.den<-ni-p
    tmp.den<-ifelse(tmp.den>0,tmp.den,1)
    for (i in 1:m)
    {
      di[i]<- (log(out1$sigma2e[i])-(sum(tmp.den*log(out1$sigma2e))/sum(tmp.den)))/sqrt(2/tmp.den[i])
    }
    H<-sum(di^2)
    p_value_sigma2e<-round(1-pchisq(H,(m-1)),4)
    test.sigma2e<-list(statistic = c(H = H), p.value = p_value_sigma2e,
                       method = "Raudenbush-Bryk's test", data.name = "sigma2e")
    class(test.sigma2e) <- "htest"
  }

  if(sigma2e.test=="Bartlett"){
    tmp.den<-ni-p
    tmp.den<-ifelse(tmp.den>0,tmp.den,1)
    Sp<-sum(tmp.den*out1$sigma2e)/sum(tmp.den)
    H<-(sum(tmp.den)*log(Sp)-(sum(tmp.den*log(out1$sigma2e))))/(1+(1/(3*(p*m-1)))*(sum((1/(tmp.den)))-(1/(sum(ni)-p*m))))
    p_value_sigma2e<-round(1-pchisq(H,(p*m-1)),4)
    test.sigma2e<-list(statistic = c(H = H), p.value = p_value_sigma2e,
                       method = "Bartlett's test",data.name = "sigma2e")
    class(test.sigma2e) <- "htest"
  }

  #Zewotir and Galpin (2007)

  resid<-NULL
  for (i in 1:m)
  {
    yi<-ys[area.s==area_id[i]]
    xsi<-xs[area.s==area_id[i],]
    if(p>1)resi<-yi-xsi%*%out1$coefficients[i,]-ranef[i]
    if(p==1)resi<-yi-xsi*out1$coefficients[i,]-ranef[i]
    resid<-c(resid,resi)
  }
  R <- diag(rep(out1$sigma2e, ni),n,n)
  G <- diag(rep(out1$sigma2u, m),m,m)
  Z<-model.matrix(ys~as.factor(area.s)-1)
  V <- Z%*%G%*%t(Z)+R
  Vinv <- solve(V)
  P <- Vinv-Vinv%*%xs%*%solve(t(xs)%*%Vinv%*%xs)%*%t(xs)%*%Vinv
  S <- R%*%P

  resid2s<-resid^2/sum(resid^2)
  cutoff<-(1-2*p/n)*(1-1/(lme.sigma2e/lme.sigma2u+n/(2*m)))
  cutoffresid2s<-mean(resid2s)+2*sd(resid2s)
  if(fig=="TRUE"){
    par(mfrow=c(1,1))
    plot(resid2s,diag(S),main="",xlab=expression(frac('e'[i]^2,paste('e'^t,'e'))),ylab=expression('s'[ii]),cex=0.7,pch=3)
    abline(h=cutoff)
    abline(v=cutoffresid2s)
  }

  #list of results
  list(coefficients=coefficients,est.sigma2u=est.sigma2u,est.sigma2e=est.sigma2e,EBP=MQ_BP, Bi=Bi,ranef=ranef,Qi=Qi,RMSE.EBP=sqrt(MSE.MQ_BP),betas.test=betas.test,test.sigma2e=test.sigma2e,test.shapiro=shapiro,test.shapiro.lme=shapiro.lme,test.shapiro.re=shapiro.re,method=est.method)

  #coefficients: estimated regression coefficients
  #est.sigma2u, est.sigma2e: estimated variance components
  #EBP: predictions at small area level
  #Bi: shrinkage factor at area level
  #Qi: area level M-quantile coefficients
  #RMSE.EBP: estimated RMSE as sqrt of g1
  #betas.test: result of the test on the varying coefficients (Hausman's test)
  #test.sigma2e: result of the test on the varying sigma2e (Raudenbush-Bryk's test or Bartlett's test)
  #test.shapiro: results of the shapiro test on the u_ij for the proposed model
  #test.shapiro.lme: results of the shapiro test on the u_ij for the lme model
  #test.shapiro.re: results of the shapiro test on the random effects for the proposed model
}
