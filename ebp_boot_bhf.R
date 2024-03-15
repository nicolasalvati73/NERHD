#' ebp_boot_bhf
#'
#' @title Bootstrap method for estimating the rmse of the EBLUP predictor under NERHD model.
#'
#' @description This function produces an estimation of the rmse of the EBLUP predictor using a bootstrap approach according the varying regression coefficients and variance components sigma2ei estimated by NERHD model.
#'
#' @usage ebp_boot_bhf<-function(xs,barX,area.s,betaq,sigma2u,sigma2e,Ni,ni,B=50,numCores=2)
#'
#' @param xs a matrix or data frame containing the explanatory variables including the intercept.
#' @param barX m*(p+1) matrix that contains the population means of each of the p auxiliary variables for the m domains. The first column is the intercept.
#' @param area.s n*1 vector (same size as y in formula) with domain codes.
#' @param betaq a named vector of estimated coefficients obtained by ebp function.
#' @param sigma2u the estimated sigma2u obtained by ebp function.
#' @param sigma2e a named vector of estimated sigma2ei obtained by ebp function.
#' @param Ni vector of domain population sizes.
#' @param ni vector of domain sample sizes.
#' @param B number of bootstrap iterations.
#' @param numCores number of cores used for paralleling the bootstrap procedure.
#'
#' @import
#' doSNOW
#' tcltk
#' lme4
#' foreach
#' iterators
#'
#' @details In this function the results of ebp_nerhd are used for creating the bootstrap population.
#'
#' @return Return a list with the following elements:
#' \item{RMSE_EBLUP_BOOT}{the bootstrap RMSE estimates for each small area}
#'
#' @references Battese, G., Harter, R. and Fuller, W. (1988) An error-components model for prediction of county crop areas using survey and satellite data. J. Am. Statist. Ass., 83, 28-36.
#' @references Lahiri, P. and Salvati, N. (2022). A Nested Error Regression Model with High Dimensional Parameter for Small Area Estimation. arXiv:2201.10235, https://doi.org/10.48550/arXiv.2201.10235.
#'
#' @author Nicola Salvati
#'
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
#' boot<-ebp_boot_bhf(xs=cbind(1,x.s),barX=XMean,area.s=regioncode.s,betaq=est$coefficients,sigma2u=est$est.sigma2u,sigma2e=est$est.sigma2e,Ni=Ni,ni=ni,B=10,numCores=2)
#'
#'
#' @export
ebp_boot_bhf<-
  function(xs,barX,area.s,betaq,sigma2u,sigma2e,Ni,ni,B=50,numCores=2)
  { #start function
    #xs: design matrix for the sample with the intercept
    #betaq: it is the matrix of regression coefficients
    #barX: the average value of X at population level
    #Ni: population size in each area
    #ni: sample size in each area
    #sigma2u: vector of variance components
    #sigma2e: vector of variance components
    #B: number of bootstrap iterations
    #area.s: small area id
    #numCores: number of cores to use in the parallel
    require(doSNOW)
    require(tcltk)
    require(foreach)
    require(iterators)

    m<-length(ni)
    p<-ncol(xs)
    area_id<-sort(unique(area.s))

    fun.boot<-function(h){
      require(lme4)
      assign("xs",xs,pos=1)
      assign("barX",barX,pos=1)
      assign("betaq",betaq,pos=1)
      assign("sigma2u",sigma2u,pos=1)
      assign("sigma2e",sigma2e,pos=1)
      assign("Ni",Ni,pos=1)
      assign("ni",ni,pos=1)
      Yibar.boot<-NULL
      EBLUP.boot<-NULL

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
      mod.lme<-lmer(y.s.boot~-1+x.s.boot+(1|area.s))
      est.ui<-as.numeric(ranef(mod.lme)$area.s$'(Intercept)')
      est.coef<-fixef(mod.lme)
      for (i in 1:m){
        EBLUP.boot[i]<-barX[i,]%*%est.coef+est.ui[i]
      }

      list(Yibar.boot=Yibar.boot,EBLUP.boot=EBLUP.boot)
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
    EBLUP.boot<-matrix(0,B,m)

    for (ii in 1:B)
    {
      Yibar.boot[ii,]<-out[[ii]]$Yibar.boot
      EBLUP.boot[ii,]<-out[[ii]]$EBLUP.boot
    }
    a_phi.ebp<-NULL
    a_phi.ebp<-(sqrt(apply((EBLUP.boot-Yibar.boot)^2,2,mean)))
    list(RMSE_EBLUP_BOOT=a_phi.ebp)
    #end function
  }
