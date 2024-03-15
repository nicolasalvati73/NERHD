#' obp_ner
#'
#' @title Observed Best Predictor under NER model
#'
#' @description This functions is used to pobtain the OBP under the NER model.
#'
#' @usage obp_ner(ys,xs,area.s,ni,Ni,barX,barXs,interval=c(0.001,1000))
#'
#' @param ys the response: a vector of length the number of rows of xs.
#' @param xs a matrix or data frame containing the explanatory variables including the intercept.
#' @param barX m*(p+1) matrix that contains the population means of each of the p auxiliary variables for the m domains. The first column is the intercept.
#' @param barXs m*(p+1) matrix that contains the sample means for each of the p auxiliary variables for the m domains. The first column is the intercept.
#' @param area.s n*1 vector (same size as y in formula) with domain codes.
#' @param Ni vector of domain population sizes.
#' @param ni vector of domain sample sizes.
#' @param interval the min and the max values of gamma.
#'
#' @import
#' MASS
#'
#' @details In this function the it is developed the algorithm for estimating the OBP under NER model.
#'
#' @return Return a list with the following elements:
#' \item{coefficients}{a named vector of coefficients}
#' \item{OPB}{a named vector of OBP at small area level}
#' \item{est.gamma}{the value of the estimated gamma}
#'
#' @references Jiang, J., Nguyen, T. and Rao, J.S. (2011). Best predictive small area estimation. Journal of the American Statistical Association, 106, 494, 732-745.
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
#'
#' XMean<-cbind(1,tapply(pop$x,pop$area,mean))
#' XsMean<-cbind(1,tapply(pop$x[s],pop$area[s],mean))
#' est<-obp_ner(ys=y.s,xs=cbind(1,x.s),area.s=regioncode.s,ni=ni,Ni=Ni,barX=XMean,barXs=XsMean,interval=c(0.001,1000))
#'
#'
#' @export
obp_ner<-function(ys,xs,area.s,ni,Ni,barX,barXs,interval=c(0.001,1000))
{
  require(MASS)
  ar<-sort(unique(area.s))
  p<-ncol(xs)
  m<-length(ar)
  fi=ni/Ni

  yi<- function(i){ matrix(c(ys[area.s==ar[i]]),ni[i],1)}
  xi<- function(i){xs[area.s==ar[i],]}
  yibar<- rep(0,m)
  yibar<-as.numeric(tapply(ys,area.s,mean))

  #parameter estimation
  mu.i.2.hat<-rep(0,m)
  for(i in 1:m)mu.i.2.hat[i]=mean(yi(i)^2) - ((Ni[i]-1)/(Ni[i]*(ni[i]-1)))* sum((yi(i)-yibar[i])^2)

  Q.gamma=function(gamma){
    tmp.Gammai = ((ni*gamma)/(1+ni*gamma))
    tmp.G=diag((fi+(1-fi)*tmp.Gammai))
    tmp.G.2=diag((fi+(1-fi)*tmp.Gammai)^2)

    tmp.H= t(barX-tmp.G%*%barXs)%*%(barX-tmp.G%*%barXs)
    tmp.L= (diag(1,m)-2*tmp.G)%*%barX +(tmp.G.2%*%barXs)

    Q<-t(yibar)%*%(tmp.G.2-tmp.L%*%ginv(tmp.H)%*%t(tmp.L))%*%yibar+t(rep(1,m))%*%(diag(1,m)-2*tmp.G)%*%mu.i.2.hat
    Q
  }
  fgamma=optimise(Q.gamma,interval=interval)
  est.gamma=fgamma$minimum
  Gammai = ((ni*est.gamma)/(1+ni*est.gamma))
  G=diag((fi+(1-fi)*Gammai))
  G.2=diag((fi+(1-fi)*Gammai)^2)
  H= t(barX-G%*%barXs)%*%(barX-G%*%barXs)
  L=(G.2%*%barXs)+(diag(1,m)-2*G)%*%barX
  beta=c(ginv(H)%*%t(L)%*%yibar)

  mu.i<-barX%*%beta+(fi+(1-fi)*(ni*est.gamma/(1+ni*est.gamma)))*(yibar-barXs%*%beta)
  list(OBP=mu.i,coefficients=beta,est.gamma=est.gamma)
}
