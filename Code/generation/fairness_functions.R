
##' plotcal() 
##' Produces a set of points for a calibration plot, and optionally plots them.
##' Uses either a binning method or a kernel method to determine height of points. 
##' In both method, divides unit interval into subintervals [0,1/n], [1/n,2/n), ... [(n-1)/n,1). 
##' For bin \[a,b\)
##'    x co-ordinate is (a+b)/2
##'   For binning method
##'    y co_ordinate is mean({y:predicted value for y is in [a,b)})
##'   For kernel method
##'    y co-ordinate is weighted mean of all y values with the weight of value yi given by dnorm(y-yi,sd=kernel_sd)
##'
##' @param y class labels, 0/1 or logical
##' @param ypred predictions Pr(Y=1), numeric vector
##' @param n number of subintervals/points
##' @param kernel set to TRUE to use kernel method
##' @param kernel_sd kernel width for kernel method; see above
##' @param conf include a confidence interval; alpha, c0, c2 are only relevant if conf=TRUE
##' @param alpha return a pointwise confidence envolope for conservative 1-alpha confidence interval
##' @param c0 for computing maximum bias; assume true covariance function is of the form a0+ a1x + a2x^2, with |a0|<c0, |a2|<c2 (c1 does not matter)
##' @param c2 for computing maximum bias; assume true covariance function is of the form a0+ a1x + a2x^2, with |a0|<c0, |a2|<c2 (c1 does not matter)
##' @param plot set to FALSE to suppress plot
##' @param ... further parameters passed to plot()
##' @return n x 2 matrix containing co-ordinates of points on the curve.
plotcal=function(y,ypred,n=10,kernel=F,kernel_sd=0.05,alpha=0.05,c0=0,c2=0.1,plot=TRUE,conf=TRUE,...)  {
  if (!kernel) {
   ycal=rep(0,n); xcal=ycal; ncal=ycal
   xup=rep(0,n); xdown=rep(0,n)
	 for (i in 1:n) {
		  sub=which((ypred> (i-1)/n) & (ypred< i/n))
 	  	if (length(sub)>5) ycal[i]=mean(y[sub]) else ycal[i]=NA # Disclosure control
 	  	if (length(sub)>5) xcal[i]=mean(ypred[sub]) else xcal[i]=NA # Disclosure control
 	  	if (length(sub)>5) ncal[i]=length(sub)
 	  	if (conf) {
 	  	  xse=sqrt(ycal[i]*(1-ycal[i])/length(sub))
 	  	  xup[i]=ycal[i] -qnorm(alpha/2)*xse
 	  	  xdown[i]=ycal[i] + qnorm(alpha/2)*xse
 	  	}
	 }
	 if (plot) plot(xcal,ycal,...)
	 return(list(x=xcal,y=ycal,n=ncal,upper=xup,lower=xdown))
  } else {  # use kernel method with given sd
   ypredsub=seq(0,1,length=n)
   kern=function(x,y) dnorm(x-y,sd=kernel_sd)
   wt=outer(ypredsub,ypred,kern)
   x1=(wt %*% ypred)
   y1=(wt %*% y)
   csub=ypredsub*y1/x1
   csub=pmax(pmin(csub,1),0)
   
   nsub=rep(0,length(csub))
   for (cc in 1:(length(csub)-1)) {
     sub=which(ypred>= csub[cc] & ypred < csub[cc+1])
     if (length(sub)>5) nsub[cc]=length(sub)
   }

   if (plot) plot(ypredsub,csub,...)
   
   if (conf) {
   # Confidence intervals
   wts=ypredsub*wt/as.vector(x1) # now csub= wts %*% y
   yvar=(wts^2 %*% as.vector(ypred*(1-ypred)))
   
   # Max bias
   bias=rowSums(outer(ypredsub,ypred,function(x,y) kern(x,y)*(c0*(y-x) + c2*(x^2 * y - y^2 * x))))/x1
   
   upper=csub - qnorm(alpha/2)*sqrt(yvar) + bias
   lower=csub + qnorm(alpha/2)*sqrt(yvar) - bias
   
   return(list(x=ypredsub,y=csub,n=nsub,lower=lower,upper=upper))
   } else {
   return(cbind(ypredsub,csub,nsub))
   }
  }
}


##' getroc() 
##' Comprehensive plotting function for receiver-operator characteristic curve. Also calculates AUROC and standard error. 
##' 
##' Rather than returning points corresponding to every cutoff, only returns a representative sample of equally-spaced points along the curve.
##'
##' SE of AUROC with no CV structure is from Hanley and McNeil 1982. SE of AUROC with CV folds is from LeDell et al 2012
##'
##' Does not plot anything. Object can be plotted in a default way.
##'
##' @param y class labels, 0/1 or logical
##' @param ypred predictions Pr(Y=1), numeric vector
##' @param cv cross-validation fold assignments, if relevant. Changes estimate of standard error.
##' @param res resolution. Returns this many equally-spaced points along the curve. Set res to null to return all points.
##' @return list containing: spec, specificity for res points in every cv fold; sens, sensitivity for res points in every cv fold; auc, areas under the curve for each fold and average (note length is 1 greater than number of CV folds); se, standard error for AUC in each fold and standard error for average auc (note length is 1 greater than number of CV folds)
getroc=function(y,ypred,cv=NULL,res=100,addauc=FALSE,cols=NULL) {
if (is.null(cv)) cv=rep(1,length(y))
if (!(length(y)==length(ypred))) stop("Parameters y and ypred should have the same length")

sens=c(); spec=c(); auc=c(); se=c(); cutoffs=c(); ncut=c()
for (i in 1:max(cv)) {
y0=y[which(cv==i)]; 
ypred0=ypred[which(cv==i)]

yt=sum(y0); yl=length(y0)
opred=order(ypred0)
#ipred=order(opred) # can use ipred to reorder in the order of original ypred

sy=y0[opred]; sp=ypred0[opred]

# Cutoffs and number of samples
cutoffs0=sp
ncut0=1:length(sp)

# Disclosure control
csy=cumsum(sy); csy[which(csy<5)]=0
csy1=cumsum(1-sy); csy1[which(csy1<5)]=0
ncut0[which(ncut0<5)]=0

sens0=1- (csy/yt)
spec0= csy1/(yl-yt)

auc0=integral(sens0,spec0)
se0=aucse(as.numeric(yt),as.numeric(yl-yt),auc0)

if (!is.null(res)) {
  ds=cumsum(sqrt((spec0[1:(yl-1)]-spec0[2:yl])^2 + (sens0[1:(yl-1)]-sens0[2:yl])^2))
  ds=ds/ds[yl-1]
  lsp=(1:(yl-1))/yl
  sub=round(yl*approx(ds,lsp,n=res)$y)
  sens0=sens0[sub]
  spec0=spec0[sub]
  cutoffs0=cutoffs0[sub]
  ncut0=ncut0[sub]
}

auc=c(auc,auc0)
se=c(se,se0)
spec=rbind(spec,spec0)
sens=rbind(sens,sens0)
cutoffs=rbind(cutoffs,cutoffs0)
ncut=rbind(ncut,ncut0)
}

if (length(auc)>1) {
  auc=c(auc,mean(auc))
  se=c(se,ci.cvAUC(ypred,y,folds=cv)$se)
}

out=list(sens=sens,spec=spec,cutoffs=cutoffs,ncut=ncut,auc=auc,se=se)
class(out)="sparraROC"
return(out)
}

# Internal function to compute SE of AUC
aucse=function(n1,n2,auc) {
  q1=auc/(2-auc); q2=2*(auc^2)/(1+auc)
  num=auc*(1-auc) + (n1-1)*(q1- (auc^2)) + (n2-1)*(q2-(auc^2))
  return(sqrt(num/(n1*n2)))
}

##' Plot function for class above
##' @param out output from getroc()
##' @param addauc set to TRUE to add text to the plot showing the (mean) AUC and SE.
##' @param cols colour to draw lines
plot.sparraROC=function(out,addauc=FALSE,cols=rep("black",dim(out$sens)[1]),...) {
  plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="1-Spec.",ylab="Sens.",...)
  ncv=dim(out$spec)[1]
  for (i in 1:ncv) lines(1-out$spec[i,],out$sens[i,],col=cols[i])
  abline(0,1,col="red",lty=2)
  auc=out$auc[length(out$auc)]
  se=out$se[length(out$se)]
  txx=paste0(signif(auc,digits=2),"+/-",signif(se,digits=2))
  if (addauc) text(0.6,0.4,txx)
}


##' getprc() 
##' Comprehensive plotting function for precision-recall curve. Also calculates AUPRC and standard error. 
##' 
##' Rather than returning points corresponding to every cutoff, only returns a representative sample of equally-spaced points along the curve.
##'
##' Does not plot anything. Object can be plotted in a default way.
##'
##' @param y class labels, 0/1 or logical
##' @param ypred predictions Pr(Y=1), numeric vector
##' @param cv cross-validation fold assignments, if relevant. Changes estimate of standard error.
##' @param res resolution. Returns this many equally-spaced points along the curve. Set res to null to return all points.
##' @return list containing: ppv, ppv for res points in every cv fold; sens, sensitivity for res points in every cv fold; auc, areas under the curve for each fold and average (note length is 1 greater than number of CV folds); se, standard error for AUC in each fold and standard error for average auc (note length is 1 greater than number of CV folds)
getprc=function(y,ypred,cv=NULL,res=100,addauc=FALSE,cols=NULL) {
  if (is.null(cv)) cv=rep(1,length(y))
  if (!(length(y)==length(ypred))) stop("Parameters y and ypred should have the same length")
  
  sens=c(); ppv=c(); auc=c(); se=c(); cutoffs=c(); ncut=c()
  for (i in 1:max(cv)) {
    y0=y[which(cv==i)]; 
    ypred0=ypred[which(cv==i)]
    
    yt=sum(y0); yl=length(y0)
    opred=order(ypred0)
    #ipred=order(opred) # can use ipred to reorder in the order of original ypred
    sy=y0[opred]; sp=ypred0[opred]

    cutoffs0=sp
    ncut0=1:length(sp)
    
    # Disclosure control
    crsy=cumsum(rev(sy)); crsy[which(crsy<5)] = 0
    csy=cumsum(sy); csy[which(csy<5)] = 0
    ncut0[which(ncut0<5)] = 0
    
    ppv0=rev(crsy/(1:length(sy)))
    sens0=1- (csy/yt)

    auc0=integral(sens0,ppv0)
    
    
    se0=sqrt(auc0*(1-auc0)/sum(y0))

    if (!is.null(res)) {
      ds=cumsum(sqrt((ppv0[1:(yl-1)]-ppv0[2:yl])^2 + (sens0[1:(yl-1)]-sens0[2:yl])^2))
      ds=ds/ds[yl-1]
      lsp=(1:(yl-1))/yl
      sub=suppressWarnings(round(yl*approx(ds,lsp,n=res)$y))
      sens0=sens0[sub]
      ppv0=ppv0[sub]
      cutoffs0=cutoffs0[sub]
      ncut0=ncut0[sub]
    }
    
    auc=c(auc,auc0)
    se=c(se,se0)
    ppv=rbind(ppv,ppv0)
    sens=rbind(sens,sens0)
    cutoffs=rbind(cutoffs,cutoffs0)
    ncut=rbind(ncut,ncut0)
  }

  if (length(auc)>1) {
    auc=c(auc,mean(auc))
    se=c(se,mean(se)/sqrt(3))
  }
  
    
  out=list(sens=sens,ppv=ppv,cutoffs=cutoffs,ncut=ncut,auc=auc,se=se)
  class(out)="sparraPRC"
  return(out)
}


##' Plot function for class above
##' @param out output from getprc()
##' @param addauc set to TRUE to add text to the plot showing the (mean) AUC and SE.
##' @param cols colour to draw lines
plot.sparraPRC=function(out,addauc=FALSE,cols=rep("black",dim(out$sens)[1]),...) {
  plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="Recall",ylab="Precision",...)
  ncv=dim(out$sens)[1]
  for (i in 1:ncv) lines(out$sens[i,],out$ppv[i,],col=cols[i])
  auc=mean(out$auc)
  se=mean(out$se)/sqrt(3)
  txx=paste0(signif(auc,digits=2),"+/-",signif(se,digits=2))
  if (addauc) text(0.6,0.4,txx)
}









##' integral() 
##' Quick form for trapezoidal integration over range of x
##'
##' @param x x co-ordinates, or nx2 matrix of points 
##' @param y y co-ordinates
##' @return trapezoidal estimate of integral of y[x] over range of x.
integral=function(x,y=NULL) {
	if (is.null(y)) {
		y=x[,2]; x=x[,1]
	}
	ox=order(x); xs=x[ox]; ys=y[ox]
	sum((xs[-1]-xs[-length(xs)])*(ys[-1]+ys[-length(ys)]))/2
}

##' ab() 
##' Shorthand to draw a red x-y line
ab=function(...) abline(0,1,col="red",...)

##' logit() 
##' Logistic function
##' @param x parameter
logit=function(x) 1/(1+exp(-x))

##' ilogit() 
##' Inverse ogistic function
##' @param x parameter between 0 and 1
ilogit=function(x) -log((1/x)-1)





##' roc_2panel
##' Draws a ROC curve (with legend) with a second panel underneath showing sensitivity difference.
##' 
##' @param rocs list of sparraROC objects (if one object, plots folds separately)
##' @param labels labels to use in legend
##' @param col line colours
##' @param lty line type, defaults to 1
##' @param xy_lty line type for x-y line, defaults to 2 (dashed)
##' @param xy_col line colour for x-y line, defaults to red
##' @param ... other parameters passed to legend()
roc_2panel=function(rocs,labels,col,lty=rep(1,length(col)),xy_lty=2,xy_col="red",...) {
  
  if ("sens" %in% names(rocs)) {
    nfold=length(rocs$auc)-1; np=length(rocs$sens)/nfold
    r0=list()
    for (i in 1:nfold) 
      r0[[i]]=list(sens=rocs$sens[i,],
         spec=rocs$spec[i,])
    rocs=r0
  }
  
  # Set up plot parameters
  par(mar=c(1,4,0.1,0.1))
  layout(mat=rbind(matrix(1,4,4),matrix(2,2,4)))
  
  # Initialise
  plot(0,xlim=c(0,1),ylim=c(0,1),xaxt="n",ylab="Sensitivity",type="n")
  abline(0,1,col=xy_col,lty=xy_lty)
  
  # x-values to draw sensitivity at
  xspec=seq(0,1,length=100)[2:99]; xr=c()
  
  # Draw ROC curves on top panel
  for (i in 1:length(rocs)) {
    xy=rocs[[i]]
    lines(1-xy$spec,xy$sens,col=col[i],lty=lty[i])
    xsens=suppressWarnings(approx(1-xy$spec,xy$sens,xspec)$y)
    xr=rbind(xr,xsens)
  }
  
  # Add legend
  legend("bottomright",legend=labels,col=col,lty=lty,...)
  
  # Bottom panel setup
  par(mar=c(4,4,0.1,0.1))
  yrr=range(t(xr)-xr[1,],na.rm=T); if (!is.finite(sum(yrr))) yrr=c(-1,1)
    plot(0,xlim=c(0,1),ylim=yrr,type="n",
    xlab="1-Specificity",ylab=expression(paste(Delta,"(sens.)")),yaxt="n")
  axis(2,at=pretty(yrr,n=2))
  
  # Draw lines on bottom panel
  for (i in 1:length(rocs)) lines(xspec,xr[i,]-xr[1,],col=col[i],lty=lty[i])
  #abline(h=0,col=xy_col,lty=xy_lty)
  
}  




##' prc_2panel
##' Draws a PRC curve (with legend) with a second panel underneath showing precision difference.
##' 
##' @param prcs list of sparraPRC objects. If of length 1, splits into folds
##' @param labels labels to use in legend
##' @param col line colours
##' @param lty line type, defaults to 1
##' @param ... other parameters passed to legend()
prc_2panel=function(prcs,labels,col,lty=rep(1,length(col)),...) {
  
  if ("ppv" %in% names(prcs)) {
    nfold=length(prcs$auc)-1; np=length(prcs$sens)/nfold
    r0=list()
    for (i in 1:nfold) 
      r0[[i]]=list(sens=prcs$sens[i,],
        ppv=prcs$ppv[i,])
    prcs=r0
  }
  
  
  
  # Set up plot 
  par(mar=c(1,4,0.1,0.1))
  layout(mat=rbind(matrix(1,4,4),matrix(2,2,4)))
  
  # Initialise plot
  plot(0,xlim=c(0,1),ylim=c(0,1),xaxt="n",ylab="Precision",type="n")
  xr=c(); 
  
  # X-values at which difference in precision will be plotted
  xsens=seq(0,1,length=100)[2:99]
  
  # Draw curves on top panel
  for (i in 1:length(prcs)) {
    px=prcs[[i]]
    lines(px$sens,px$ppv,col=col[i],lty=lty[i])
    xppv=suppressWarnings(approx(px$sens,px$ppv,xsens)$y)
    xr=rbind(xr,xppv)
  }
  
  # Add legend
  legend("topright",legend=labels,col=col,lty=lty,...)
  
  # Initialise bottom panel
  par(mar=c(4,4,0.1,0.1))
  yrr=range(t(xr)-xr[1,],na.rm=TRUE); if (!is.finite(sum(yrr))) yrr=c(-1,1)
  plot(0,xlim=c(0,1),ylim=yrr,type="n",
    xlab="Recall",ylab=expression(paste(Delta,"(Prec.)")),
    yaxt="n")
  axis(2,at=pretty(range(t(xr)-xr[1,]),n=2))
  
  # Draw lines on bottom panel
  for (i in 1:length(prcs)) lines(xsens,xr[i,]-xr[1,],col=col[i],lty=lty[i])
}



##' cal_2panel
##' Draws calibration curves (with legend) with a second panel underneath showing predicted differences.
##' 
##' @param cals list of calibration objects, output from plotcal(). 
##' @param labels labels to use in legend
##' @param col line colours
##' @param lty line type, defaults to 1
##' @param xy_lty line type for x-y line, defaults to 2 (dashed)
##' @param xy_col line colour for x-y line, defaults to red
##' @param ci_col colours to draw confidence intervals on lower panel; NA to not draw. 
##' @param 
##' @param ... other parameters passed to legend()
cal_2panel=function(cals,labels,col,lty=rep(1,length(col)),xy_lty=2,xy_col="red",ci_col=rep(NA,length(col)),...) {

    
  # Setup plot
  par(mar=c(1,4,0.1,0.1))
  layout(mat=rbind(matrix(1,4,4),matrix(2,2,4)))
  
  # Initialise plot
  plot(0,xlim=c(0,1),ylim=c(0,1),xaxt="n",ylab="Observed",type="n")
  
  # Draw lines on top panel
  xr=c()
  for (i in 1:length(cals)) {
    cx=cals[[i]]
    lines(cx$x,cx$y,col=col[i],lty=lty[i])
    xr=rbind(xr,cx$y-cx$x)
  }
  abline(0,1,col=xy_col,lty=xy_lty,lwd=1)
  
  # Draw legend
  legend("topleft",legend=labels,col=col,lty=lty,...)
  
  
  # Initialise bottom panel
  par(mar=c(4,4,0.1,0.1))
  xpred=seq(0,1,length=100)[2:99]
  yrr=range(xr,na.rm=T); if (!is.finite(sum(yrr))) yrr=c(-1,1)
  plot(0,xlim=c(0,1),ylim=yrr,type="n",
    xlab="Predicted",ylab=expression(paste(Delta,"(cal.)")),
    yaxt="n")
  axis(2,at=pretty(range(xr),n=3))
  
  # Draw lines on bottom panel
  for (i in 1:length(cals)) {
    cx=cals[[i]]
    if (!is.na(ci_col[i])) {
      xx=c(cx$x,rev(cx$x))
      yy=c(cx$lower,rev(cx$upper))-c(cx$x,rev(cx$x))
      wn=which(is.finite(xx+yy)); xx=xx[wn]; yy=yy[wn];
      polygon(xx,yy,
        col=ci_col[i],border=NA)
    }
    lines(cx$x,cx$y-cx$x,col=col[i],lty=lty[i])
  }
  abline(h=0,col=xy_col,lty=xy_lty)
  
}  




