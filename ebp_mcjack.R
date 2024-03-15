#' ebp_mcjack
#'
#' @title The McJack estimator of the rmse of ebp predictor.
#'
#' @description The McJack provides the estimate of the root mse for each samll area.
#'
#' @usage ebp_mcjack(ys, xs, area.s, betaq, sigma2u, sigma2e, est.method = c("local", "global"), method=c("LBP","Naive"), barX, Qi, Ni, B = 50, tol = 1e-06, maxit = 100, k_b = 100, k_sigma_u = 100, k_sigma_e = 100, numCores = 2)
#'
#' @param ys the response: a vector of length the number of rows of xs.
#' @param xs a matrix or data frame containing the explanatory variables including the intercept.
#' @param area.s n*1 vector (same size as y in formula) with domain codes.
#' @param betaq a named vector of estimated coefficients obtained by nerhd function.
#' @param sigma2u the estimated sigma2u obtained by nerhd function.
#' @param sigma2e a named vector of estimated sigma2ei obtained by nerhd function.
#' @param est.method estimation method for sigma2ei. It could be local or global.
#' @param method a character string. It is the method for computing the M-quantile coefficients at area level. It could be LBP or Naive. LBP uses Empirical Linear Best (ELB) predictor, while Naive is the traditional method suggested by Chambers & Tzavidis (2006) for computing the M-quantile coefficients at area level by averaging the unit level M-quantile coefficients. Default is LBP.
#' @param barX m*(p+1) matrix that contains the population means of each of the p auxiliary variables for the m domains. The first column is the intercept.
#' @param Qi Estimated M-quantile coefficients at area level from nerhd function.
#' @param Ni vector of domain population sizes.
#' @param B number of bootstrap iterations.
#' @param tol the accuracy for the stopping criterion.
#' @param maxit the limit on the number of algorithm iterations.
#' @param k_b tuning constant used for Huber influence function for the estimation of the coefficients.
#' @param k_sigma_u tuning constant used for Huber influence function for the estimation of the sigma2u.
#' @param k_sigma_e tuning constant used for Huber influence function for the estimation of the sigma2ei.
#' @param numCores number of cores used for paralleling the bootstrap procedure.
#'
#' @import
#' MASS
#'
#' @return RMSE_EBP_MJ is a vector that returns the estimated rmse by McJack procedure for each area.
#'
#' @references Lahiri, P. and Salvati, N. (2022). A Nested Error Regression Model with High Dimensional Parameter for Small Area Estimation. arXiv:2201.10235, https://doi.org/10.48550/arXiv.2201.10235.
#' @references Jiang, J., Lahiri, P., and Nguyen, T. (2018). A unified monte-carlo jackknife for small area estimation after model selection. Annals of Mathematical Sciences and Applications, 3 (2):405–438
#'
#' @author Nicola Salvati
#'
#' @seealso ebp_nerhd
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
#' dta<-data.frame(ys=y.s,xs=x.s,area.s=regioncode.s)
#' est<-ebp_nerhd(ys~1+xs,random=~1|area.s,est.method="global",barX=XMean,Ni=Ni,Q=P,method="LBP",tol=1e-06,maxit=100,k_b=1.345,k_sigma_u=1.345,k_sigma_e=100,fig="TRUE",uij.error.test="BHF",sigma2e.test="Raudenbush_Bryk",data=dta)
#' mcjack<-ebp_mcjack(ys=y.s,xs=cbind(1,x.s),area.s=regioncode.s,betaq=est$coefficients,sigma2u=est$est.sigma2u,sigma2e=est$est.sigma2e,est.method="global",method="LBP",barX=XMean,Qi=est$Qi,Ni=Ni,B=10,tol=1e-06,maxit=100,k_b=1.345,k_sigma_u=1.345,k_sigma_e=100,numCores=2)
#'
#'
#' @export
ebp_mcjack<-
  function(ys,xs,area.s,betaq,sigma2u,sigma2e,est.method=c("local","global"),method=c("LBP","Naive"),barX,Qi,Ni,B=50,tol=1e-06,maxit=100,k_b=100,k_sigma_u=100,k_sigma_e=100,numCores=2)
  {
    #ys: sampled response variable
    #xs: design matrix for the sample with the intercept
    #betaq: it is the matrix of regression coefficients
    #barX: the average value of X at population level
    #Ni: population size in each area
    #ni: sample size in each area
    #sigma2u: vector of variance components
    #sigma2e: vector of variance components
    #B: number of bootstrap iterations
    #tol=1e-06
    #maxit=100
    #k_b=100: the tuning constant for betas
    #k_sigma_u=100: the tuning constant for sigma2u
    #k_sigma_e=100: the tuning constant for sigma2e
    #numCores: number of cores to use in the parallel
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
    {
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
    mqre_mcjack_l <-
      function(qtl,y,x,group,group.out,tol=1e-06,maxit=100,k=100,k_sigma_u=100,k_sigma_e=100)
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
        gg=unique(group)
        tmp.ni<-NULL
        tmp.ni[gg]<-ni
        tmp.ni[group.out]<-0
        ni<-tmp.ni
        areanumber=length(qtl)
        p=ncol(x)


        #z=matrix(0,n,m)
        z<-model.matrix(y~as.factor(group)-1)
        #kk=0
        #for (j in 1:m){for(i in 1:ni[j]){
        # kk=kk+1
        #  z[kk,j]=1}}

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

        if(dim(x)[2]==1){ls <- lm(y ~ 1)}else{
          ls <- lm(y ~ x[,-1])}
        beta1<-c(ls$coefficients)
        estsigma2u<-sigmasq0.v
        estsigma2e<-sigmasq0
        sigma1G.old<-NULL
        sigma1E.old<-NULL
        sigma1E.old<-rep(sigmasq0,length(qtl))
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
        for (i in gg)
        {
          #STEP 1 Estimation of beta
          xbeta <- c(x %*% beta.old[i,])
          R <- diag(rep(sigma1E.old[i], n))
          G <- diag(rep(sigma1G.old, m))
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
          ri<-as.vector(r)[group==i]
          #mod.rlmerE<-rlmer(r~1+(1|group),rho.sigma.e = psi2propII(smoothPsi, k = k_sigma),
          #   rho.sigma.b = psi2propII(smoothPsi, k = k_sigma))
          psi.r <- my.psi(as.vector(ri), k = k_sigma_e, q=qtl[i])
          V.inv.i<-V.inv[group==i,group==i]
          U.i<-U[group==i,group==i]
          ZZ.i<-ZZ[group==i,group==i]
          if(ni[i]>1){
            qq <- V.inv.i %*% sqrt(U.i) %*% psi.r
            KK5<-sqrt(U.i) %*% V.inv.i
            KK0<-t(psi.r)%*% KK5
            A1=KK0%*%qq
            A2=KK0%*%ZZ.i%*%qq
          }
          if(ni[i]==1){
            qq <- V.inv.i * sqrt(U.i) * psi.r
            KK5<-sqrt(U.i) * V.inv.i
            KK0<-t(psi.r)* KK5
            A1=KK0*qq
            A2=KK0*ZZ.i*qq
          }
          A <- matrix(c(A1[1],A2[1]), nrow = 2, ncol=1)
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
            t4=KK3*ZZ.i
          }

          T1<- matrix(c(t1,t3,t2,t4), nrow = 2, ncol=2)
          sigma=chol2inv(chol(T1))%*%A
          sigma1E.est[i]<-sigma[1,1]

          if(p==2)residuals<-c(residuals,(y[group==i]-x[group==i,2]*beta.est[i,2]))
          if(p>2)residuals<-c(residuals,(y[group==i]-x[group==i,-1]%*%beta.est[i,-1]))
        }#end beta estimation

        #STEP 3: computation variance components sigma^2_gamma
        residuals<-residuals-mean(residuals)
        R <- diag(rep(sigma1E.old[-group.out], ni[-group.out]))
        G <- diag(rep(sigma1G.old, m))
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
        B0<-mean(beta.est[,1])

        fehler <- sum((beta.old[-group.out,-1] -beta.est[-group.out,-1])^2)+ sum((sigma1G.old -sigma1G.est)^2+sum((sigma1E.old[-group.out] -sigma1E.est[-group.out])^2))
        ifelse (iter < maxit, ifelse(fehler > tol, {iter=iter+1}, {break}), {break})
        beta.old <- beta.est
        sigma1G.old<-sigma1G.est
        sigma1E.old<-sigma1E.est
        } #end while, estimation procedure


        if(sigma1G.est<0){lme_fit<-lmer(y~1+x[,-1]+(1|group))
        sigma1G.est<-summary(lme_fit)$varcor$group[1]}
        list(coefficients =cbind(rep(B0,areanumber),beta.est[,-1]),sigma2u=(sigma1G.est),sigma2e=(sigma1E.est),varBeta=var.beta,quantile=qtl)

      }

    mqre_mcjack_g <-
      function(qtl,y,x,group,group.out,tol=1e-06,maxit=100,k=100,k_sigma_u=100,k_sigma_e=100)
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
        areanumber=length(qtl)
        p=ncol(x)
        gg=unique(group)

        #z=matrix(0,n,m)
        z<-model.matrix(y~as.factor(group)-1)
        #kk=0
        #for (j in 1:m){for(i in 1:ni[j]){
        # kk=kk+1
        #  z[kk,j]=1}}

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

        if(dim(x)[2]==1){ls <- lm(y ~ 1)}else{
          ls <- lm(y ~ x[,-1])}
        beta1<-c(ls$coefficients)
        estsigma2u<-sigmasq0.v
        estsigma2e<-sigmasq0
        sigma1G.old<-NULL
        sigma1E.old<-NULL
        sigma1E.old<-rep(sigmasq0,length(qtl))
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
        for (i in 1:areanumber)
        {
          #STEP 1 Estimation of beta
          xbeta <- c(x %*% beta.old[i,])
          R <- diag(rep(sigma1E.old[i], n))
          G <- diag(rep(sigma1G.old, m))
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

          if(p==2)residuals<-c(residuals,(y[group==i]-x[group==i,2]*beta.est[i,2]))
          if(p>2)residuals<-c(residuals,(y[group==i]-x[group==i,-1]%*%beta.est[i,-1]))
        }#end beta estimation

        #STEP 3: computation variance components sigma^2_gamma
        residuals<-residuals-mean(residuals)
        R <- diag(rep(sigma1E.old[-group.out],ni))
        G <- diag(rep(sigma1G.old, m))
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
        B0<-mean(beta.est[,1])

        fehler <- sum((beta.old[,-1] -beta.est[,-1])^2)+ sum((sigma1G.old -sigma1G.est)^2+sum((sigma1E.old -sigma1E.est)^2))
        ifelse (iter < maxit, ifelse(fehler > tol, {iter=iter+1}, {break}), {break})
        beta.old <- beta.est
        sigma1G.old<-sigma1G.est
        sigma1E.old<-sigma1E.est
        } #end while, estimation procedure


        if(sigma1G.est<0){lme_fit<-lmer(y~1+x[,-1]+(1|group))
        sigma1G.est<-summary(lme_fit)$varcor$group[1]}
        list(coefficients =cbind(rep(B0,areanumber),beta.est[,-1]),sigma2u=(sigma1G.est),sigma2e=(sigma1E.est),varBeta=var.beta,quantile=qtl)

      }

    n<-length(ys)
    ni<-table(area.s)
    m<-length(ni)
    p<-ncol(xs)
    ar<-sort(unique(area.s))

    BP_LB_tilde_aiphi<-NULL
    mod<-ebpmq_boot_mj(xs=xs,barX=barX,area.s=area.s,betaq=betaq,est.method=est.method,method=method,sigma2u=sigma2u,sigma2e=sigma2e,Ni=Ni,ni=ni,B=B,tol=tol,maxit=maxit,k_b=k_b,k_sigma_u=k_sigma_u,k_sigma_e=k_sigma_e,numCores=numCores)
    BP_LB_tilde_aiphi<-mod$a_phi.ebp

    #Jackknife
    BP_LB_tilde_aiphi_l<-matrix(0,m,m)
    for (i in 1:m)
    {
      xJack<-xs[area.s!=ar[i],]
      yJack<-ys[area.s!=ar[i]]
      groupJack<-area.s[area.s!=ar[i]]
      if(est.method=="local")out1<-mqre_mcjack_l(qtl=Qi,y=yJack,x=xJack,group=groupJack,group.out=ar[i],tol=tol,maxit=maxit,k=k_b,k_sigma_u=k_sigma_u,k_sigma_e=k_sigma_e)
      if(est.method=="global")out1<-mqre_mcjack_g(qtl=Qi,y=yJack,x=xJack,group=groupJack,group.out=ar[i],tol=tol,maxit=maxit,k=k_b,k_sigma_u=k_sigma_u,k_sigma_e=k_sigma_e)
      est.sigma2u<-NULL
      est.sigma2u<-out1$sigma2u
      est.sigma2e<-NULL
      coefficients<-matrix(0,m,p)
      for (j in 1:m)
      {
        est.sigma2e[j]<-out1$sigma2e[j]
        coefficients[j,]<-out1$coefficients[j,]
      }
      if(est.method=="local"){est.sigma2e[i]<-sigma2e[i]
      coefficients[i,]<-betaq[i,]}

      modJack<-ebpmq_boot_mj(xs=xs,barX=barX,area.s=area.s,est.method=est.method,method=method,betaq=coefficients,sigma2u=est.sigma2u,sigma2e=est.sigma2e,Ni=Ni,ni=ni,B=B,tol=tol,maxit=maxit,k_b=k_b,k_sigma_u=k_sigma_u,k_sigma_e=k_sigma_e,numCores=numCores)
      BP_LB_tilde_aiphi_l[i,]<-modJack$a_phi.ebp
      print(paste("McJack iteration n°",i))
    }
    BP_LB_hat_aphi<-NULL

    for (i in 1:m)
    {
      BP_LB_hat_aphi[i]<-BP_LB_tilde_aiphi[i]-((m-1)/(m))*sum((BP_LB_tilde_aiphi_l[,i]-BP_LB_tilde_aiphi[i]))
    }

    RMSE_BP_LB<-exp(BP_LB_hat_aphi)
    list(RMSE_EBP_MJ=RMSE_BP_LB)
  }
