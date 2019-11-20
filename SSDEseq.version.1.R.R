library(MASS) 
library(pastecs)
#library(mixtools)
library(mclust)
library(gplots)

FUN.cross=function(xx, n, mus, sds, sp)
  {n.out=1000
  sx=seq(min(xx), max(xx), length.out = n.out)
  fa=function(ix) sp[1]*dnorm(ix, mus[1], sds[1])
  fb=function(ix) sp[2]*dnorm(ix, mus[2], sds[2])
  fita=fa(sx)
  fitb=fb(sx)
  
  fitts=abs(fita-fitb)
  tpf=turnpoints(fitts)$tppos
  cross=sx[tpf][which(fitts[tpf]==min(fitts[tpf]))]
  
  if (length(tpf)>5) 
    {simx1=rnorm(10000, mus[1], sds[1])
    simx2=rnorm(10000, mus[2], sds[2])
    simx=sort(c(simx1, simx2))
    c1=ecdf(simx1)
    c2=ecdf(simx2)
    d.cdf=abs(c1(simx)-c2(simx))
    md.cdf=max(d.cdf)
    KScross=simx[which(d.cdf==md.cdf)][1]
    diff.2cross=abs(sx[tpf]-KScross)
    cross=sx[tpf][which(diff.2cross==min(diff.2cross))]
    cross=NA
  }
  
  KS.p=d.cdf=NA
  if (!is.na(cross))
    {cdf1=pnorm(cross, mus[1], sds[1])
    cdf2=pnorm(cross, mus[2], sds[2])
    d.cdf=abs(cdf1-cdf2)
    KS.p=exp(-2*d.cdf^2*round(n*sp[1])*round(n*sp[2])/n)
    KS.p[KS.p>1]=1
    }
  return(list(cross=cross, KS.p=KS.p, d.cdf=d.cdf))
  }

FUN.cross1=function(xx, n, mus, sds, sp)
  {n.out=1000
  sx=seq(min(xx), max(xx), length.out = n.out)
  fa=function(ix) sp[1]*dnorm(ix, mus[1], sds[1])
  fb=function(ix) sp[2]*dnorm(ix, mus[2], sds[2])
  fita=fa(sx)
  fitb=fb(sx)

  fitts=abs(fita-fitb)
  tpf=turnpoints(fitts)$tppos
  cross=sx[tpf][which(fitts[tpf]==min(fitts[tpf]))]

  if (length(tpf)>5) 
    {simx1=rnorm(10000, mus[1], sds[1])
    simx2=rnorm(10000, mus[2], sds[2])
    simx=sort(c(simx1, simx2))
    c1=ecdf(simx1)
    c2=ecdf(simx2)
    d.cdf=abs(c1(simx)-c2(simx))
    md.cdf=max(d.cdf)
    KScross=simx[which(d.cdf==md.cdf)][1]
    diff.2cross=abs(sx[tpf]-KScross)
    cross=sx[tpf][which(diff.2cross==min(diff.2cross))]
    cross=NA
    }

  return(list(cross=cross))
  }
FUN.boxcox=function(x, x0)
{box= boxcox(exp(x)~1, lambda=seq(-8, 8, .1), plotit =FALSE)
  cox=data.frame(box$x, box$y)
  cox2=cox[with(cox, order(-cox$box.y)),]
  blambda=cox2[1, "box.x"]
  if (blambda>0)
  {x=(exp(x)^blambda-1)/blambda
  x0=(exp(x0)^blambda-1)/blambda
  } else
    if (blambda<0)
    {x=-1 * exp(x) ^ blambda
    x0=-1 * exp(x0) ^ blambda
    }
  return(list(x, x0))
}

FUN.RM=function(x)
  {x=x[!is.na(x)]
  ox=sort(x)
  if (abs(ox[1]-ox[2])>10*sd(x)) {ox[1]=ox[2]}
  ox=rev(sort(ox))
  if (abs(ox[1]-ox[2])>10*sd(x)) {ox[1]=ox[2]}
  x=ox
  return (x)
  }

FUN.LRT=function (x, x0, boxcox, tit, include.normal, iteration, eps, verbal, plotTF, useKS)
  {K=2
  x=FUN.RM(x)
  x0=x0[!is.na(x0)]  
  n=length(x)

  if (missing(useKS)) {useKS=FALSE}
  if (missing(tit)) {tit=""}
  if (missing(plotTF)) {plotTF=FALSE}
  if (missing(verbal)) {verbal=TRUE}
  if (missing(boxcox)) {boxcox=TRUE}
  if (missing(include.normal)) {include.normal=TRUE}
  if (sum(x==min(x))>min(3, n*.03))
     {stx=sort(unique(round(x,3)))
     range=stx[2]-stx[1]
     x[x==min(x)]=x[x==min(x)]+runif(sum(x==min(x)), -range, range)
     }
  
  if (boxcox)
     {bcx=FUN.boxcox(x,x0)
     x=bcx[[1]]
     x0=bcx[[2]]
     }
  xx=c(x,x0)
  
  if (missing(iteration)) {iteration=1000}
     mx0=mean(x0)
     n0=length(x0)
     se0=sd(x0)/sqrt(n0)
  

  if (missing(eps)) {eps=.0000000001}
  
  CCC=matrix(rep(NA, n*K), ncol=K)
  
  denx=density(x)
  dd=ts(denx$y)
  tp=turnpoints(dd)$tppos
  
  if (length(tp)>1) 
     {xy=data.frame(denx$x[tp], denx$y[tp])
     xy=xy[rev(order(xy[,2])),]
     mus=xy[1:2,1]
     } else
     {dis=tp-c(1,length(dd))
     wdis=which(abs(dis)==min(abs(dis)))
     tp2=c(tp, tp-round(dis/2)[3-wdis])
     mus=denx$x[tp2]
     }
  
  sp=rep(.5, K)
  sds=rep(sd(x), K)*c(.5,.5)
  cross0=loglik0=0
   
  for (B in 1:iteration)
  { 
     for (i in 1:2)
        {CCC[,i]=sp[i]*dnorm(x, mus[i], sds[i])}
        CCC=CCC/rowSums(CCC) 
    
     for (i in 1:2)
        {sp[i]=mean(CCC[,i]) }
    
     sp=sp/sum(sp)
    
     for (i in 1:2)
        {if (min(sp)*n>3)
            {mus[i]=mean(      (CCC[,i])*x)/sp[i]
            sds[i]=sqrt(mean(      (CCC[,i])*(x-mus[i])^2)/sp[i])
            if (sds[i]==0) {break}
            if (abs(mus[i]-mx0)==min(abs(mus-mx0)))
               {if (include.normal) 
                   {mus[i]=(sds[i]^2*mean(x0)+sum( round(CCC[,])* (CCC[,i])*x)*se0^2)/(sds[i]^2+n*sp[i]*se0^2)}
               }
            } 
         }
    
     if (min(sds)>0)
     {loglik=   sum(log(dnorm(x, mus[1], sds[1])*sp[1]
                                    +dnorm(x, mus[2], sds[2])*sp[2]))  
        KS.p=d.cdf=NA
        if (useKS)
        {crosss=FUN.cross1(xx, n, mus, sds, sp)
        cross=crosss[[1]]
        }
     }

      if (useKS) {if (abs(cross-cross0)< eps) {break};cross0=cross} else
        {if (abs(loglik-loglik0)< eps) {break};loglik0=loglik}
  }
  
     muxx=mean(x)
     sdxx=sd(x)
     loglik.null=sum( log(dnorm(x, muxx, sdxx)))
     crosss=FUN.cross(xx, n, mus, sds, sp)
     cross=crosss[[1]]
     KS.p=crosss[[2]]
     d.cdf=crosss[[3]]
     
  if (min(sp)*n<=3|min(sds)==0)
     {sp[sp==min(sp)]=0
     sp[sp==max(sp)]=1
     mus=rep(mean(x, na.rm=TRUE), K)
     sds=rep(sd(x, na.rm=TRUE), K)
     loglik=loglik.null
     KS.p=1
     d.cdf=NA
     } else {if (loglik<loglik.null) {loglik=loglik.null}}
  
  if (verbal) {print(paste("Converge at iteration ", B))}
  
  fa=function(ix) sp[1]*dnorm(ix, mus[1], sds[1])
  fb=function(ix) sp[2]*dnorm(ix, mus[2], sds[2])

  
chisq.p=NA
  if (loglik!=loglik.null)
     {#cross=sx[tpf][which(fitts[tpf]==min(fitts[tpf]))]
      if (include.normal&(!is.na(cross)))
         {gp=c(rep(0,length(x0)), rep(1,length(x)))
         chisq.p=chisq.test(table(data.frame(gp,c(x0<cross,x<cross))))$p.value
         }
     } 

  if (plotTF)
  {xx=c(x0,x)     
  sx=seq(min(xx), max(xx), length.out = 100)
  fita=fa(sx)
  fitb=fb(sx)
  fit=fita+fitb
  
  hx=hist(x, max(30,round(n/5)), plot=FALSE); 
  hx0=hist(x0, max(30,round(n/5)), plot=FALSE);  

  plot( hx, col=4, border= 4, xlim=c(min(xx), max(xx)), xlab="Normalized gene expression",
        main=tit)  # first histogram

  plot( hx0, col=3, border=rgb(0, 1, 0, 1/4), add=TRUE)
  
  points(sx, fita*n/sum(hx[[3]]), lwd=2, col=2, type="l")
  points(sx, fitb*n/sum(hx[[3]]), lwd=2, col=2, type="l")
  points(sx, fit *n/sum(hx[[3]]), lwd=5, col=2, type="l")
  points(sx, fit *n/sum(hx[[3]]), lwd=2, col=7, type="l")
  if (exists(as.character("normalmixEM")))
     {mout=normalmixEM(x, k=2, fast = TRUE, arbmean = TRUE, epsilon = 1e-11)
      mfita=mout$lambda[1]*dnorm(sx, mout$mu[1], mout$sigma[1])*1
      mfitb=mout$lambda[2]*dnorm(sx, mout$mu[2], mout$sigma[2])*1
      mfit=mfita+mfitb
      points(sx, mfita*n/sum(hx[[3]]), lwd=2, col=2, type="l")
      points(sx, mfitb*n/sum(hx[[3]]), lwd=2, col=2, type="l")
      points(sx, mfit *n/sum(hx[[3]]), lwd=4, col=2, type="l")
  }
  points(sx, dnorm(sx, mean(x), sd(x)) *n/sum(hx[[3]]), lwd=1, col=8, type="l")
  points(sx, dnorm(sx, mean(x), sd(x)) *n/sum(hx[[3]]), lwd=1, col=8, type="l")
  abline(v=cross, col=7)
  }
  out=list("Proportion"=sp, "Means"=mus, "SDs"=sds, "loglik.null"=loglik.null, "loglik"=loglik, "sd.sd.ratio"=min(sds)/max(sds), "cutoff"=cross, "x"=x, "x0"=x0, "-2logLRT"=2*(loglik-loglik.null), "chisq.p"=chisq.p, "KS.p"=KS.p)
  
  return(out)
}


LRT.null=function(n, n0, num.Sim, verbal)
{ if (missing(num.Sim)) {num.Sim=100000}
  
  BBLRT=rep(NA, num.Sim)
  for (s in 1:num.Sim)
  {x=rnorm(n, 20)
  x0=rnorm(n0, 20)
  oott= FUN.LRT(x, x0, verbal=verbal)
  BBLRT[s]=oott$`-2logLRT`
  }
  return(BBLRT)
}


SSDEseq=function(Mx, Mx0, plotTF, boxcox, verbal, useKS)
{ if (missing(useKS)) {useKS=FALSE}
 if (missing(verbal)) {verbal=TRUE}
  if (missing(plotTF)) {plotTF=FALSE}
  if (missing(boxcox)) {boxcox=TRUE}
  n=dim(Mx)[[2]]
  n0=dim(Mx0)[[2]]
  
  ngene=nrow(Mx)
  if (table(row.names(Mx)==row.names(Mx0))!=ngene) {print("The two groups have to match!"); break}
  
  gene.loglik <<- data.frame("gene"=row.names(Mx),"-2logLRT"=rep(NA, ngene), "mean expression"=rep(NA, ngene), "minor proportion"=rep(NA, ngene), "cutoff"=rep(NA, ngene), "t.test"=rep(NA, ngene), "t.p"=rep(NA, ngene), "df"=rep(NA, ngene), "chisq.p"=rep(NA, ngene), "sd.sd.ratio"=rep(NA, ngene), "KS.p"=rep(NA, ngene))
  Mx=as.matrix(Mx) 
  Mx0=as.matrix(Mx0) 
  
  for (ng in 1:ngene)
  {if (ng/100==round(ng/100)) {print (paste("Just completed the ", ng, "th gene", sep=""))}
   x0=Mx0[ng,]
   x=Mx[ng,]
   m.tit=row.names(Mx)[ng]

   if (sum(x0==min(x0))<n0*0.05 & sum(x)>n) 
     {oott= FUN.LRT(x, x0, tit=m.tit, boxcox=boxcox, plotTF=plotTF, verbal=verbal, useKS=useKS)

     gene.loglik[ng,2] <<- oott$`-2logLRT`
     gene.loglik[ng,3] <<- mean(x)
     gene.loglik[ng,4] <<- min(oott[[1]])
     gene.loglik[ng,5] <<- oott[[7]]
     nn1=round(n*oott[[1]][1])
     nn2=round(n*oott[[1]][2])
     den=oott[[3]][1]^2/nn1+oott[[3]][2]^2/nn2
     ts=abs(diff(oott[[2]]))/sqrt(den)
     vs=(den)^2/(oott[[3]][1]^4/nn1^2/(nn1-1)+oott[[3]][2]^4/nn2^2/(nn2-1))
     gene.loglik[ng,6] <<- ts
     gene.loglik[ng,7] <<- 2*(1-pt(ts, vs))
     gene.loglik[ng,8] <<- vs
     gene.loglik[ng,9] <<- oott$chisq.p
     gene.loglik[ng,10] <<- round(oott$sd.sd.ratio, 3)
     gene.loglik[ng,11] <<- oott$KS.p
   }
  }
  out = gene.loglik[(!is.na(gene.loglik$X.2logLRT))&(!is.na(gene.loglik$t.p)),]
  
  return(out)
}
   

subgroup.selection=function(Mx, outp, boot.LRT, n.group, plotTF)
{if (missing(plotTF)) {plotTF=FALSE}
  sig=outp[outp$chisq.p<0.05/nrow(outp),]
  sig=sig[(!is.na(sig$`X.2logLRT`))&(!is.na(sig$`X.2logLRT`))&(sig$`X.2logLRT`>boot.LRT),]
  sig=sig[sig$minor.proportion >.03,]
  sig=sig[sig$t.p< (0.05/nrow(outp)),]
  sig=sig[!is.na(sig$cutoff),]
  sig=sig[sig$KS.p<0.05/nrow(outp),]
  ogroup=FUN.LRT(-log10(sig$KS.p),-log10(sig$KS.p), boxcox=FALSE, include.normal = FALSE, useKS=TRUE)
     if (ogroup$`-2logLRT`>30)
     {sig=sig[-log10(sig$KS.p)>ogroup$cutoff,]}
 
  
  Msig=Mx[row.names(Mx)%in%sig$gene, ]

  Mcor=(cor((Msig)))
  c.dist=dist(Mcor)
  c.fit=hclust(c.dist)            
  c.groups <- cutree(c.fit, k=n.group)  
   
  if (plotTF) {plot(c.fit)
              rect.hclust(c.fit, k=n.group, border=6, which=3)
              heatmap.2(Mcor, hclustfun = hclust, col=redgreen(75), trace="none")
              heatmap.2(cor(t(Msig)), hclustfun = hclust, col=redgreen(75), trace="none")
              }
  
  dist = dist(t(Msig))
  fit = hclust(dist)
  groups <- cutree(fit, k=n.group)  
 
  if (plotTF) {fit$labels=substr(fit$labels, 12, 17)
               plot(fit)
               rect.hclust(fit, k=n.group, border=c(4,2), which=c(5,6)) 
               heatmap.2(as.matrix((Msig)), hclustfun = hclust, col=redgreen(75), trace="none")
              }

 
  gid=data.frame(ID=names(groups), group=groups)

  return(gid)
}



