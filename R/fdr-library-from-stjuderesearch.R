##########################################################################
# Function: compute.fdr
# Purpose: Apply the requested fdr method to the set(s) of p-values provided
# Arguments: p - a vector or matrix of p-values
#            opts - a list with the following components
#                   LC = loess.control options for SPLOSH
#                   dplaces = pre-rounding option for SPLOSH
#                   lam = lambda for Storey's q-value
#                   robust = robust option for Storey's q-value
#                   adj.eps = epsilon adjustment for BUM to avoid applying MLE to p-values with 0
#                   n.adj = adjustment factor for BUM
#            notes - a vector of strings giving some notes about the analysis
# Returns: a list with the following components:
#          p - the set of p-values supplied to the function
#          fdr - an estimate of the ratio of the number of false positives to the number of results with lesser or equal p-value
#          q - the q-value or fdr-adjusted p-value
#          cdf - the estimate of the expected proportion of tests with lesser or equal p-value
#          ord - gives the value of indices to order results by ascending p-value
#          pi - estimate of the proportion of tests with a true null hypothesis (n.a. for BH95 & BY01)
#          fp - estimate of the total number of false positives among results with equal or lesser p-value
#          fn - estimate of the total number of false negatives among results with larger p-value
#          toterr - sum of fp and fn
#          fdr.method - indicates which method was used
#          notes - appends string description of method used to notes argument supplied
#########################################################################

compute.fdr<-function(p,
                      fdr.method="PC04",
                      opts=default.opts(fdr.method),
                      notes=NULL
                     )

{
	res<-"Choose a valid fdr.method"
	if (fdr.method=="PC04")    res<-splosh(p,opts$LC,opts$dplaces,notes)
	if (fdr.method=="St02")    res<-storey(p,opts$lam,opts$robust,notes)
	if (fdr.method=="PM03")    res<-bum(p,opts$adj.eps,opts$n.adj,notes)
	if (fdr.method=="BH95")    res<-BH95(p,notes)
	if (fdr.method=="Ch04")    res<-Ch04(p,opts$qpos,notes)
	if (fdr.method=="BY01")    res<-BY01(p,notes)
	return(res)
}

#############################################################
# Function: fdr.plots
# Purpose: Produce diagnostic and illustrative plots for fdr analysis
# Arguments: fdrobj - object returned from compute.fdr
# Returns: NULL
#############################################################

fdr.plots<-function(fdrobj,alpha=1,sub.labels="",main="")

{
	if(is.vector(fdrobj$p)) fdrobj$p<-matrix(fdrobj$p,length(fdrobj$p),1)
	h<-dim(fdrobj$p)[2]
	g<-dim(fdrobj$p)[1]
	
	if (length(main)==1) main<-rep(main,h)
	if (length(sub.labels)==1) 
	{
		if ((sub.labels=="")&&(h>1)) sub.labels<-paste("Null Hyp.",1:h)
		sub.labels<-rep(sub.labels,h)
	}
	nscreen<-3+2*(!is.element(fdrobj$fdr.method,c("BH95","BY01","St02")))+
                 any(!c(is.null(fdrobj$fp),is.null(fdrobj$fn),is.null(fdrobj$toterr)))
        split.screen(c(h,nscreen))
	k<-0
	for (i in 1:h)
	{
		# p-value histogram
		k<-k+1
                screen(k)
		hist(fdrobj$p[,i],prob=T,xlab="p-value",ylab="density",main=main[i],sub=sub.labels[i],col=16)
		if (!is.null(fdrobj$pdf)) 
		{
			H<-hist(fdrobj$p[,i],prob=T,plot=F)
			x<-fdrobj$p[fdrobj$ord[,i],i]
			y<-fdrobj$pdf[fdrobj$ord[,i],i]
			inc<-y<max(H$counts)
			lines(x[inc],y[inc],lwd=5)
		}		
		if (!is.null(fdrobj$pi)) lines(c(0,1),rep(fdrobj$pi[i],2),lwd=5)

		inc<-(fdrobj$p[,i]<=alpha)
		g.alpha<-sum(inc)
		ord.inc<-1:g.alpha
		# P-P plots
		if (!is.element(fdrobj$fdr.method,c("BH95","BY01","St02")))
		{
			k<-k+1
			screen(k)
			plot(ord.inc/g,fdrobj$cdf[fdrobj$ord[ord.inc,i],i],
			     xlab="Empirical p-value CDF",
			     ylab="Estimated p-value CDF",
			     main=main[i],sub=sub.labels[i],type="l")
			lines(c(0,max(g.alpha/g,fdrobj$cdf[fdrobj$ord[ord.inc,i],i])),c(0,max(g.alpha/g,fdrobj$cdf[fdrobj$ord[ord.inc,i],i])))			
			k<-k+1
			screen(k)
			plot(fdrobj$p[fdrobj$ord[ord.inc,i],i],
			     fdrobj$cdf[fdrobj$ord[ord.inc,i],i]-ord.inc/g,
			     xlab="p-value",ylab="p-value CDF Estimate - p-value EDF",
			     main=main[i],sub=sub.labels[i])
			lines(c(0,1),c(0,0))
		}
		
		k<-k+1
		screen(k)
		plot((1:g.alpha)/g,fdrobj$fdr[fdrobj$ord[inc,i],i],
		     xlab="Proportion of Tests Declared Significant",
		     ylab="Raw Empirical FDR",type="p",main=main[i],sub=sub.labels[i])

		k<-k+1
		screen(k)
		plot(fdrobj$p[fdrobj$ord[ord.inc,i],i],fdrobj$q[fdrobj$ord[ord.inc,i],i],
		     type="l",xlab="p-value",ylab="q",sub=sub.labels[i],main=main[i])

		if (any(!c(is.null(fdrobj$fp),is.null(fdrobj$fn),is.null(fdrobj$toterr))))
		{
			r<-range(fdrobj$fp[inc,i],fdrobj$fn[inc,i],fdrobj$toterr[inc,i])
			k<-k+1
			screen(k)
			plot(c(0,alpha),r/g,xlab="p-value threshold",ylab="Error Estimates",
			     type="n",sub=sub.labels[i],main=main[i])
			if (!is.null(fdrobj$fp))     lines(fdrobj$p[fdrobj$ord[ord.inc,i],i],fdrobj$fp[fdrobj$ord[ord.inc,i],i]/g,col=3)
			if (!is.null(fdrobj$fn))     lines(fdrobj$p[fdrobj$ord[ord.inc,i],i],fdrobj$fn[fdrobj$ord[ord.inc,i],i]/g,col=2)
			if (!is.null(fdrobj$toterr)) lines(fdrobj$p[fdrobj$ord[ord.inc,i],i],fdrobj$toterr[fdrobj$ord[ord.inc,i],i]/g,col=1)
			legend(x="topleft",legend=c("FPs","FNs","TEs"),lty=1,col=3:1)
		}
	}
}

#########################################################
# Function: fdr.table
# Purpose: Write results to a tab-delimited text file (Excel will read if the filename has .xls extension)
# Arguments: fdrobj - an object returned by fdr.compute
#            outfile - string or string-vector giving output file name(s)
#            stat - vector or matrix of test statistics corresponding to p-values in fdr.object (default=NULL)
#            top - no. of features to include in file, i.e. include the top 100 (default=100)
#            sep - delimiter for output file(s) (default = "\t", i.e. tab)
#            type - type of table: "abb" for abbreviated or "full" for complete
#            features - vector giving the feature names (default=NULL)
#            annotations - vector, matrix, or data.frame giving the annotations of the features (default=NULL)
# Returns: NULL
#############################################################

fdr.table<-function(fdrobj,outfile,stat=NULL,top=100,sep="\t",type="abb",features=NULL,annotations=NULL)

{
	if(is.vector(fdrobj$p)) fdrobj$p<-matrix(fdrobj$p,length(fdrobj$p),1)
	h<-dim(fdrobj$p)[2]
	g<-dim(fdrobj$p)[1]
	top<-min(top,g)
	if (length(outfile)!=h) return(paste("Your p-value matrix has",h,
		                                  "columns but you have specified",length(outfiles),"output file names.",
		                                  "Please specify",h,"output file names."))
	if (!is.null(stat))
	{
		if(is.vector(stat)) stat<-matrix(stat,length(stat),1)
		if(any(dim(stat)!=dim(fdrobj$p))) return(paste("The p-value matrix has dimension ",dim(p)[1]," by ",dim(p)[2],
			                                    ".  The stat matrix has dimension ",dim(stat)[1]," by ",dim(stat)[2],
			                                    ".  Please give matrices with equal dimension."))
		
	}
	if (is.null(features)) features<-paste("Feature",1:g)
	if (type=="abb")
	{
		for (i in 1:h)
		{
			if (!is.null(stat)) this.stat<-stat[,i]
			else this.stat<-NULL
			ord<-fdrobj$ord[,i]
			X<-cbind(features,stat=this.stat,p=fdrobj$p[,i],q=fdrobj$q[,i],annotations)
			write.table(X[ord[1:top],],outfile[i],sep=sep,row.names=FALSE)
		}
	}
	if (type=="full")
	{
		for (i in 1:h)
		{
			if (!is.null(stat)) this.stat<-stat[,i]
			else this.stat<-NULL	
			
			if (!is.null(fdrobj$fp)) this.fp<-fdrobj$fp[,i]
			else this.fp<-NULL
			
			if (!is.null(fdrobj$fn)) this.fn<-fdrobj$fn[,i]
			else this.fn<-NULL
			
			if (!is.null(fdrobj$toterr)) this.te<-fdrobj$toterr[,i]
			else this.te<-NULL
				
			if (!is.null(fdrobj$ebp)) this.ebp<-fdrobj$ebp[,i]
			else this.ebp<-NULL
				
			if (!is.null(fdrobj$pdf)) this.pdf<-fdrobj$pdf[,i]
			else this.pdf<-NULL	
			
			ord<-fdrobj$ord[,i]
			X<-cbind(features,stat=this.stat,p=fdrobj$p[,i],q=fdrobj$q[,i],fdr=fdrobj$fdr[,i],ebp=this.ebp,fp=this.fp,fn=this.fn,
			         toterr=this.te,cdf=fdrobj$cdf[,i],pdf=this.pdf,annotations)
			write.table(X[ord[1:top],],outfile[i],sep=sep,row.names=FALSE)			
		}
	}
}

default.opts<-function(fdr.method)

{
	if (fdr.method=="PC04") return(list(LC=loess.control(),dplaces=4))
	if (fdr.method=="St02") return(list(lam=NULL,robust=F))
	if (fdr.method=="PM03")    return(list(adj.eps=0.05,n.adj=10000))
	if (fdr.method=="Ch04")    return(list(qpos=c(0.00625,0.01,0.0125,0.025,0.05,0.1,0.25,0.5,0.75)))
}

splosh<-function(p,LC=loess.control(),dplaces=6,notes=NULL)

{
	if (is.vector(p)) p<-matrix(p,length(p),1)
	h<-dim(p)[2]
	g<-dim(p)[1]
	fdr<-q<-cdf<-fp<-fn<-toterr<-pdf<-ebp<-ord<-matrix(0,g,h)
	pi<-p.pi<-rep(0,h)
	for (i in 1:h)
	{
		res<-splosh.1(p[,i],LC,dplaces,notes)
		fdr[,i]<-res$fdr
		q[,i]<-res$q
		cdf[,i]<-res$cdf
		fp[,i]<-res$fp
		fn[,i]<-res$fn
		toterr[,i]<-res$toterr
		pdf[,i]<-res$pdf
		ebp[,i]<-res$ebp
		ord[,i]<-res$ord
		pi[i]<-res$pi
		p.pi[i]<-res$p.pi
	}
	return(list(p=p,fdr=fdr,q=q,cdf=cdf,fp=fp,fn=fn,toterr=toterr,pdf=pdf,ebp=ebp,ord=ord,pi=pi,p.pi=p.pi,
	            dplaces=dplaces,                       # rounded to this no. of decimal places
	            LC=LC,                                 # loess.control options used
	            fdr.method="PC04",
	            notes=c(notes,res$notes)))
	            
}

storey<-function(p,lam=NULL,robust=F,notes=NULL)

{
	if (is.vector(p)) p<-matrix(p,length(p),1)
	h<-dim(p)[2]
	g<-dim(p)[1]
	fdr<-q<-cdf<-fp<-fn<-toterr<-ord<-matrix(0,g,h)
	pi<-rep(0,h)
	for (i in 1:h)
	{
		res<-storey.1(p[,i],lam,robust,notes)
		fdr[,i]<-res$fdr
		q[,i]<-res$q
		cdf[,i]<-res$cdf
		fp[,i]<-res$fp
		fn[,i]<-res$fn
		toterr[,i]<-res$toterr
		ord[,i]<-res$ord
		pi[i]<-res$pi
	}
	return(list(p=p,fdr=fdr,q=q,cdf=cdf,ord=ord,pi=pi,
	            lam=lam,                       
	            robust=robust,                                 
	            fdr.method="St02",
	            notes=c(notes,res$notes)))
	            
}

BH95<-function(p,notes=NULL)

{
	if (is.vector(p)) p<-matrix(p,length(p),1)
	h<-dim(p)[2]
	g<-dim(p)[1]
	fdr<-q<-cdf<-fp<-fn<-toterr<-ord<-matrix(0,g,h)
	pi<-rep(0,h)
	for (i in 1:h)
	{
		res<-BH95.1(p[,i])
		fdr[,i]<-res$fdr
		q[,i]<-res$q
		cdf[,i]<-res$cdf
		fp[,i]<-res$fp
		fn[,i]<-res$fn
		toterr[,i]<-res$toterr
		ord[,i]<-res$ord
		pi[i]<-res$pi
	}
	return(list(p=p,fdr=fdr,q=q,cdf=cdf,ord=ord,
	            fdr.method="BH95",
	            notes=c(notes,res$notes)))
}

BY01<-function(p,notes=NULL)

{
	if (is.vector(p)) p<-matrix(p,length(p),1)
	g<-dim(p)[1]
	h<-dim(p)[2]
	const<-sum(1/(1:g))
	res95<-BH95(p,notes)
	fdr<-const*res95$fdr
	q<-const*res95$q
	fdr[fdr>1]<-1
	q[q>1]<-1
	pi<-rep(const,h)
	cdf<-res95$cdf
	ord<-res95$ord
	return(list(p=p,fdr=fdr,q=q,cdf=cdf,ord=ord,
	            fdr.method="BY01",
	            notes=c("Benjamini and Yekutieli (2001) Used to Control FDR.",notes)))
}

bum<-function(p,adj.eps=0.5,n.adj=length(p),notes=NULL)

{
	if (is.vector(p)) p<-matrix(p,length(p),1)
	h<-dim(p)[2]
	g<-dim(p)[1]
	fdr<-q<-cdf<-fp<-fn<-toterr<-pdf<-ebp<-ord<-matrix(0,g,h)
	a<-lambda<-pi<-rep(0,h)
	for (i in 1:h)
	{
		res<-bum.1(p[,i],adj.eps,n.adj,notes)
		fdr[,i]<-res$fdr
		q[,i]<-res$q
		cdf[,i]<-res$cdf
		fp[,i]<-res$fp
		fn[,i]<-res$fn
		toterr[,i]<-res$toterr
		pdf[,i]<-res$pdf
		ebp[,i]<-res$ebp
		ord[,i]<-res$ord
		pi[i]<-res$pi
		a[i]<-res$a
		lambda[i]<-res$lambda
	}
	return(list(p=p,fdr=fdr,q=q,cdf=cdf,fp=fp,fn=fn,toterr=toterr,pdf=pdf,ebp=ebp,ord=ord,pi=pi,a=a,lambda=lambda,
	            adj.eps=adj.eps,                       # rounded to this no. of decimal places
	            n.adj=n.adj,                                 # loess.control options used
	            fdr.method="PM02",
	            notes=c(notes,res$notes)))
	            
}

Ch04<-function(p,qpos=c(0.00625,0.01,0.0125,0.025,0.05,0.1,0.25,0.5,0.75),notes=NULL)

{
	if (is.vector(p)) p<-matrix(p,length(p),1)
	h<-dim(p)[2]
	g<-dim(p)[1]
	fdr<-q<-cdf<-pdf<-ord<-matrix(0,g,h)
	pi<-rep(NA,h)
	for (i in 1:h)
	{
		res<-Ch04.1(p[,i],qpos)
		fdr[,i]<-res$fdr
		q[,i]<-res$q
		cdf[,i]<-res$cdf
		pdf[,i]<-res$pdf
		ord[,i]<-res$ord
		pi[i]<-res$pi
	}
	return(list(p=p,fdr=fdr,q=q,cdf=cdf,pdf=pdf,pi=pi,
	            fdr.method="Ch04",notes=c(notes,res$notes)))
}



#####################################
# Function: BH95
# Purpose: Compute quantity that Benjamini and Hochberg use to control the FDR
#####################################

BH95.1<-function(p,                  # p-value vector
                 notes=NULL)         # vector of notes gathered to date

{
	m<-length(p)                   # number of genes examined in analysis
	u<-order(p)                    # vector to order results by ascending p-values
	v<-rank(p)                     # rank of p-values
	cdf<-v/m                       # empirical distribution function used to estimate cdf
	fdr<-p/cdf                     # quantity used to control fdr
	fdr[fdr>1]<-1                  # any larger than 1, assign 1

	q<-fdr                         # initialize q-value vector
	q[u[m]]<-min(fdr[u[m]],1)      # get the ball rolling
	q[rev(u)]<-cummin(q[rev(u)])   # find q-value
	fp<-m*cdf*fdr                  # estimated false positive count, very conservative
	fn<-(1-cdf-(1-p))*m            # estimated false negative count, so conservative often near 0
	fn[fn<0]<-0                    # set false negative count to 0 if less than 0
	toterr<-fp+fn                  # compute total error criterion
	# return results
	return(list(fdr=matrix(fdr,m,1),           # fdr estimate
	            q=matrix(q,m,1),               # q-value estimate               
	            cdf=matrix(cdf,m,1),           # cdf estimate
	            fp=matrix(fp,m,1),             # estimated no. of false positives incurred
	            fn=matrix(fn,m,1),             # estimated no. of false negatives incurred
	            toterr=matrix(toterr,m,1),     # total estimated no. of errors
	            ord=u,                         # ordering vector
	            pi=1,                          # implied estimate of null proportion
	            fdr.method="Benjamini and Hochberg (1995)",
	            notes=c(notes,"Benjamini and Hochberg 1995 Used to Estimate FDR.")
	           ))
}

###########################################################
# storey
# computes FDR and q-values according to Storey's method
###########################################################

storey.1<-function(
                 p,            # vector of p-values
                 lam=NULL,     # smoothing parameter for estimating pi
                 robust=F,     # indicates whether to compute robust q-values
                 notes=NULL    # notes or comments on analysis
                )

{
    #This is just some pre-processing
    m <- length(p)
    notes<-c(notes,"Storey's method used to estimate FDR.")
    #These next few functions are the various ways to automatically choose lam
    #and estimate pi0
    if(!is.null(lam)) {
        pi0 <- mean(p>lam)/(1-lam)
        pi0 <- min(pi0,1)
        notes <- c(notes,"The user prespecified lam in the calculation of pi0.")
    }
    else{
        notes <- c(notes,"A smoothing method was used in the calculation of pi0.")
        lam <- seq(0,0.95,0.01)
        pi0 <- rep(0,length(lam))
        for(i in 1:length(lam)) {
        pi0[i] <- mean(p>lam[i])/(1-lam[i])
        }
        spi0 <- smooth.spline(lam,pi0,df=3,w=(1-lam))
        pi0 <- predict(spi0,x=0.95)$y
        pi0 <- min(pi0,1)
    }

    #The q-values are actually calculated here
    u <- order(p)
    v <- rank(p)
    fdr <- pi0*m*p/v
    fdr[fdr>1]<-1
    q<-fdr
    if(robust) 
    {
        q <- pi0*m*p/(v*(1-(1-p)^m))
        notes <- c(notes, "The robust version of the q-value was calculated. See Storey JD (2002) JRSS-B 64: 479-498.")
    }
    q[u[m]] <- min(q[u[m]],1)
    q[rev(u)]<-cummin(q[rev(u)])
    cdf<-v/m
    fp<-m*cdf*fdr
    fn<-(1-cdf-pi0*(1-p))*m
    fn[fn<0]<-0

   return(list(
               fdr=matrix(fdr,m,1),         # FDR estimate
               q=matrix(q,m,1),             # q-values
               cdf=matrix(cdf,m,1),         # CDF estimate (EDF)
               fp=matrix(fp,m,1),           # Est. No. False Positives
               fn=matrix(fn,m,1),           # Est. No. False Negatives
               toterr=matrix(fp+fn,m,1),    # Total of FP and FN
               ord=u,                       # ordering vector, places in order of increasing p-values
               pi=pi0,                      # pi - null proportion

               lam=lam,                     # lam - smoothing parameter for estimating pi
               robust=robust,               # robust q-values or not
               fdr.method="Storey (2002), Storey and Tibshirani (2003)",

               notes=notes      # notes regarding analysis
              ))
}
########################################################################################################
# This S-plus code implements the SPLOSH method for the analysis of microarray data
# described by Pounds and Cheng (Improving False Discovery Rate Estimation - Bioinformatics 2004).
#
# Last Update: July 14, 2004
# 
# Function name: splosh
# Arguments: pvals - a vector of p-values obtained by applying a statistical hypothesis
#                    test to the expression values for each gene individually
#            LC - (optional) a list containing tuning parameters for the loess algorithm, 
#                 see help(loess.control) for more details, defaults to S-plus defaults
#            dplaces - (optional) round p-values to this many decimal places to avoid
#                      numerical difficulties that can arise by division by small numbers, default = 6
#            notes - (optional) may pass any string describing notes regarding the analysis, default = NULL
#
# Outputs: a list with the following elements
#         fdr - a vector of the cfdr estimates at the values specified by the user in pvals
#         q - a vector of corresponding monotone fdr estimates, similar to Storey's q-value
#         cdf - vector of the SPLOSH cdf estimate evaluated at pvals
#         fp - vector of estimated no. of false positives incurred by setting threshold 
#              at corresponding value of pvals
#         fn - vector of estimated no. of false negatives incurred by setting threshold 
#              at corresponding value of pvals
#         toterr - vector giving the sum of fp and fn for each entry
#         ebp - the empirical Bayes posterior that the null hypothesis is true, see Efron et al. (JASA 2001)
#         pi - estimate of the mixing parameter giving the probability of the null hypothesis
#         ord - a vector giving index numbers to order results according to ascending p-values
#         dplaces - echoes value of dplaces used in algorithm
#         LC - echoes value of LC used in algorithm
#         fdr.method - gives bibliographic reference for the method
#         notes - returns notes regarding the analysis
###########################################################################################################

splosh.1<-function(pvals,LC=loess.control(),dplaces=3,notes=NULL)
     
{
	notes<-c(notes,"SPLOSH used to estimate FDR.")
	x<-round(pvals,dplaces)                                      # Prevent numerical difficulties
	n<-length(pvals)                                             # Obtain no. of points
	r<-rank(x)                                                   # rank observations
	ux<-sort(unique(x))                                          # ordered unique p-values
	ur<-(sort(unique(r)))-.5                                     # ordered unique ranks
	#plot(ux,ur)	
	if (max(ux)<1)                                               # edge adjustments
	{
		ux<-c(ux,1)
		ur<-c(ur,n)
	}
	if (min(ux)>0)
	{
		ux<-c(0,ux)
		ur<-c(0,ur)
	}  
	nux<-length(ux)                                              # No. of unique points
	ur<-ur/n                                                     # Divide ranks by n
	dx<-diff(ux)                                                 # spacings of unique p-values
	dr<-diff(ur)                                                 # spacings of unique ranks
	dFdx<-exp(log(dr)-log(dx))                                   # spacing-interval slopes
	mr<-(ur[1:(nux-1)]+ur[2:nux])/(2)                            # Unitized Medians of Unique Rank Intervals
	mx<-(ux[1:(nux-1)]+ux[2:nux])/2                              # Unitized Medians of Unique p-value Intervals
	tr<-asin(2*(mr-.5))                                          # Arc-Sine Transform Rank Midpoints
	fit<-loess(log(dFdx)~tr,loess.control=LC)                    # Lowess of log-transformed EDQF and transformed rank-midpoints
	tr2<-asin(2*(ur-.5))                                         # Arc-Sine Transform Unique
	yhat<-(predict(fit,tr2))                               # Obtain Estimated Derivative of CDF (up to constant multiplier)
	yhat[c(1,nux)]<-approx(tr2[2:(nux-1)],
	                       yhat[2:(nux-1)],
	                       xout=tr2[c(1,nux)],
	                       rule=3)$y                             # Extrapolate to get endpoint estimates
	pdf1<-exp(yhat)                                              # Obtain PDF up to a constant
	trap<-0.5*(pdf1[1:(nux-1)]+pdf1[2:nux])*dx                   # Trapezoid rule terms
	const<-sum(trap)                                             # Constant via trapezoid rule
	pdf2<-pdf1/const                                             # Adjust pdf by constant
	pdf<-approx(ux,pdf2,xout=pvals,rule=3)$y                     # Obtain pdf estimate for p-values
	cdf<-approx(ux,c(0,cumsum(trap)),
	            xout=pvals,rule=2)$y/const                       # Obtain cdf
	pi<-min(pdf)                                                 # Take pi as min of PDF
	p.pi<-max(pvals[pdf==pi])                                    # largest p-value at which pi occurs
	cdf.pi<-max(cdf[pvals==p.pi])                                # CDF evaluated at p.pi
	cfdr<-pi*pvals/cdf                                           # Estimate cFDR for p-values
	ebp<-pi/pdf                                                  # Estimate EBP of null for p-values
	cfdr[pvals==0]<-ebp[pvals==0]                                # L'Hospital's Rule
	cfdr[cfdr>1]<-1                                              # Set max(fdr)=1
	qval<-cfdr                                                   # Initialize q-value vector
	ordering<-order(pvals)                                       # Report Results
	qval[rev(ordering)]<-cummin(qval[rev(ordering)])             # Compute q-value
	fp<-n*cfdr*cdf                                               # Compute Estimated No. False Positives
	fn<-n*(cdf.pi-cdf-pi*(p.pi-pvals))*(pvals<=p.pi)             # Compute Estimated No. False Negatives
	toterr<-fp+fn                                                # Compute Estimated Total No. Errors
	
	return(list(fdr=matrix(cfdr,n,1),                  # cfdr estimate
	            q=matrix(qval,n,1),                    # q-value estimate
	            cdf=matrix(cdf,n,1),                   # cdf estimate
	            fp=matrix(fp,n,1),                     # estimate of no. of false postives
	            fn=matrix(fn,n,1),                     # estimate of no. of false negatives
	            toterr=matrix(toterr,n,1),             # estimated total no. of errors
	            pdf=matrix(pdf,n,1),                   # estimate of pdf
	            ebp=matrix(ebp,n,1),                   # estimate of ebp
	            pi=pi,                                 # estimate of null proportion
	            p.pi=p.pi,                             # p-value at which estimated pdf achieves minimum
	            cdf.pi=cdf.pi,                         # cdf evaluated at p.pi
	            ord=ordering,                          # ordering vector
	            dplaces=dplaces,                       # rounded to this no. of decimal places
	            LC=LC,                                 # loess.control options used
	            fdr.method="Pounds and Cheng (2004)",  # Bibliographic Information
	            notes=notes))                          # Notes regarding analysis
}

bum.1<-function(p,adj.eps=0.5,n.adj=length(p),notes=NULL)

#################################################################
# Function: bum
# Purpose: Compute FDR estimates using beta-uniform mixture model
#          proposed by Pounds and Morris (2003).
#################################################################

{
	m<-length(p)   # No. of genes
	notes<-c(notes,"Pounds and Morris (2003) BUM model used to estimate FDR.")
	# adjust p-values if necessary
	if (any(p==0)) 
	{
		p<-(n.adj*p+adj.eps)/(n.adj+2*adj.eps)
		notes<-c(notes,"The p-values were adjusted because some p-values equal zero.")
	}
	res<-bum.mle(p)                                                 # find MLE's of parameters
	pi<-ext.pi(res$a,res$lambda)                                    # determine pi from MLE's
	fdr<-bum.FDR(p,res$a,res$lambda)                                # determine fdr estimates
	q<-fdr                                                          # q-value = fdr due to implied monotonicity
	cdf<-pbum(p,res$a,res$lambda)                                   # cdf estimate
	fp<-m*bum.false.positive(p,res$a,res$lambda)                    # estimated no. of false positives
	fn<-m*bum.false.negative(p,res$a,res$lambda)                    # estimated no. of false negatives
	toterr<-fp+fn                                                   # estimated total no. of errors
	pdf<-dbum(p,res$a,res$lambda)                                   # estimated pdf
	ebp<-1-bum.emp.post(p,res$a,res$lambda)                         # estimated ebp
	# return results
	return(list(fdr=matrix(fdr,m,1),                                # fdr estimate
	            q=matrix(q,m,1),                                    # q-value estimate
	            cdf=matrix(cdf,m,1),                                # cdf estimate
	            fp=matrix(fp,m,1),                                  # estimated no. of false positives
	            fn=matrix(fn,m,1),                                  # estimated no. of false negatives
	            toterr=matrix(toterr,m,1),                          # total error
	            pdf=matrix(pdf,m,1),                                # pdf estimate
	            ebp=matrix(ebp,m,1),                                # ebp estimate
	            pi=pi,                                              # pi estimate (null proportion)
	            ord=order(p),                                       # ordering vector
	            a=res$a,                                            # MLE of BUM parameter a
	            lambda=res$lambda,                                  # MLE of BUM parameter lambda
	            fdr.method="Pounds and Morris (2003)",              # method used to estimate FDR
	            notes=notes))
}


##########################################################################
# Function: dbum
# Purpose: Compute the pdf of the bum distribution
# Arguments: x - point or vector of points at which to compute pdf
#            a - shape parameter of beta component
#            lambda - mixture parameter, proportion of uniform component
# Returns: value of the pdf of the bum distribution for x
##########################################################################

dbum<-function(x,a,lambda)

{
	return(lambda+(1-lambda)*a*x^(a-1))
}

##########################################################################
# Function: pbum
# Purpose: Compute the cdf of the bum distribution
# Arguments: x - point or vector of points at which to compute the pdf
#            a - shape parameter of the beta component
#            lambda - mixing parameter, weight of uniform component
# Returns: value of the cdf of the bum distribution for x
##########################################################################

pbum<-function(x,a,lambda)
{
	return(lambda*x+(1-lambda)*x^a)
}

##########################################################################
# Function: qbum
# Purpose: Compute the quantile of the bum distribution
# Arguments: p - percentile or vector of percentiles
#            a - shape parameter of the beta component
#            lambda - mixing parameter, weight of uniform component
#            nbisect - the number of bisections to perform
# Returns: the values x such that the pdf of x equals the percentiles
# Notes: Uses the bisection method to find quantiles, a larger value of
#        nbisect will result in more accurate results.  The results will
#        be accurate to within 2^(-nbisect).
##########################################################################

qbum<-function(p,a,lambda,nbisect=20)

{

	n<-length(p)
	mid<-rep(0,n)
	top<-rep(1,n)
	bot<-rep(0,n)
	gohigher<-rep(F,n)
	for (j in 1:nbisect)
	{
		mid<-(top+bot)/2
		gohigher<-(pbum(mid,a,lambda)<p)
		bot[gohigher]<-mid[gohigher]
		top[!gohigher]<-mid[!gohigher]
	}
	return(mid)
}

#########################################################################
# Function: rbum
# Purpose: generate random bum observations
# Arguments: n - number of observations
#            a - parameter a
#            lambda - parameter lambda
# Returns: a vector of random bum observations
##########################################################################

rbum<-function(n,a,lambda)

{
	u<-runif(n)
	return(qbum(u,a,lambda))
}

#########################################################################
# Function: dalt
# Purpose: compute the density of the alternative distribution
# Arguments: x - p-value of interest
#            a - shape parameter of beta component
#            lambda - mixing parameter, proportion of uniform component
# Returns: vector of density at x
########################################################################

dalt<-function(x,a,lambda)

{
	pi<-ext.pi(a,lambda)
	return((dbum(x,a,lambda)-pi)/(1-pi))
}

#########################################################################
# Function: palt
# Purpose: compute the cdf of the extracted alternative component of 
#          the BUM distribution
# Arguments: x - pvalue of interest
#            a - shape parameter of beta component
#            lambda - mixing parameter, proportion of uniform component
# Returns: vector of cdf at x
#########################################################################

palt<-function(x,a,lambda)

{
	pi<-ext.pi(a,lambda)
	return((pbum(x,a,lambda)-pi*x)/(1-pi))
}

##########################################################################
# Function: qalt
# Purpose: compute the quantile of the alternative component of BUM
# Arguments: p - the percentile
#            a - shape parameter of beta component
#            lambda - mixing parameter, proportion of uniform component
# Returns: a vector the p quantiles
###########################################################################

qalt<-function(p,a,lambda,nbisect=20)

{
	n<-length(p)
	mid<-rep(0,n)
	for (i in 1:n)
	{
		top<-1
		bot<-0
		for (j in 1:nbisect)
		{
			mid[i]<-mean(c(top,bot))
			if(palt(mid[i],a,lambda)<p[i]) bot<-mid[i]
			else top<-mid[i]	
		}
	}
	return(mid)
}

#########################################################################
# Function: bum.logL
# Purpose: for a set of p-values x, compute the log-likelihood of 
#          parameters a and lambda
# Arguments: a - the shape parameter of the beta component
#            lambda - mixing parameter, proportion of uniform component
#            x - vector of p-values
# Returns: the log-likelihood of the parameters
# Notes: Permutation techniques can result in p-values of zero.  This
#        routine will provide non-numeric results if a p-value of zero
#        is included.  A recommendation is to adjust permutation p-values
#        before calling the function by letting the new p-values
#        equal (nperms*old p-values + .5)/(nperms+1).
#########################################################################

bum.logL<-function(a,lambda,x)

{
  return(sum(log(dbum(x,a,lambda))))
}

#########################################################################
# Function: logit
# Purpose: Compute the logit of an x in (0,1)
# Arguments: x - point or vector of points in (0,1)
# Returns: logit of x
#########################################################################

logit<-function(x)

{
	return(log(x)-log(1-x))
}

#########################################################################
# Function: inv.logit
# Purpose: Compute the inverse logit of x
# Arguments: x - point or vector of points
# Returns: inverse logit of x
#########################################################################

inv.logit<-function(x)

{
	return(exp(x)/(1+exp(x)))
}

#########################################################################
# Function: neg.bum.logL
# Purpose: Compute the negative log-likelihood for use by nlminb in
#          bum.mle.
# Arguments: x - a vector with two components 
#                x[1] corresponds to logit of a, beta shape parameter
#                x[2] corresponds to logit of lambda, the prop. uniform
#            pvals - a vector of non-zero p-values
# Returns: the negative log-likelihood 
# Notes: used in bum.mle
#########################################################################

neg.bum.logL<-function(x,pvals)

{
	return(-1*sum(log(dbum(pvals,inv.logit(x[1]),inv.logit(x[2])))))
}

#########################################################################
# Function: bum.mle
# Purpose: For a set of p-values, compute MLE's for a and lambda
# Arguments: p-vals - the set of p-values
#            nstartpts - number of starting points to randomly generate
#            starta - vector of logit of starting a's for optimization routine, default = 0
#            startlambda - vector of logit of starting lambdas for optimization, default = 0
#            nadj - number for adjusting p-values if necessary, default = 100000
#            adjeps - epsilon for adjusting pvalues, if necessary, default = 0.5
# Returns: optimal values of a and lambda for given starting values
# Notes: May not provide global mle.  Repeating with multiple starting
#        points may yield different results.
#########################################################################

bum.mle<-function(pvals,nstartpts=0,starta=0,startlambda=0,nadj=100000,adjeps=0.5)
{
	if (nstartpts>0)
	{
		starta<-logit(runif(nstartpts))
		startlambda<-logit(runif(nstartpts))
	}
	adjusted<-min(pvals)<=0
	if (adjusted) pvals<-(nadj*pvals+adjeps)/(nadj+2*adjeps)
	bestlogL<- -Inf
	nstartpts<-length(starta)	
	for (i in 1:nstartpts)
	{
		results <- nlm( neg.bum.logL, c(starta[i],startlambda[i]),pvals=pvals )		#	OK in R
		#    name of value	in R results$minimum as results$objective in S+
		if ( -results$minimum > bestlogL)
		{
			a<-inv.logit(results$estimate[1])		# 	parameters[1])  # S+
			lambda<-inv.logit(results$estimate[2])	# 	parameters[2])	# S+
			logL<- -results$minimum
			nits<- results$iterations
			termination <- results$code				#	message		#	S+
			bestlogL<-logL
		}
	}
	return(list(a=a,lambda=lambda,logL=logL,pvals.adjusted=adjusted,nstartpts=nstartpts,nits=nits,termination=termination))
}

#########################################################################
# Function: qqbum
# Purpose: Produce a quantile-quantile plot for a set of p-values
# Arguments: pvals - a vector of p-values
#            a - shape parameter for beta component of bum distribution
#            lambda - mixing parameter for BUM distribution
#            main - primary plot tile, default = "BUM QQ Plot"
#            xlab - label of x-axis, default = "BUM Expected p-value"
#            ylab - label of y-axis, default = "Observed p-value"
# Returns: a quantile-quantile plot
# Notes: If values of a and lambda are not provided, bum.mle
#        will be used to estimate a and lambda.  May take a few minutes.
#########################################################################


qqbum<-function(pvals,a=NULL,lambda=NULL,main="BUM QQ Plot",xlab="BUM Expected p-value",ylab="Observed p-value",nstartpts=0,starta=0,startlambda=0)

{
	n<-length(pvals)
	pvals<-sort(pvals)	
	if(is.null(a)||is.null(lambda))
	{
		al<-bum.mle(pvals,nstartpts=nstartpts,starta=starta,startlambda=startlambda)
		plot(c(0,1),c(0,1),main=main,xlab=xlab,ylab=ylab,type="n")
		lines(qbum((rank(pvals)-.5)/n,al$a,al$lambda),pvals,lty=2)
		lines(c(0,1),c(0,1))
	}
	else
	{
		plot(c(0,1),c(0,1),main=main,xlab=xlab,ylab=ylab,type="n")
		lines(qbum((rank(pvals)-.5)/n,a,lambda),pvals,lty=2)
		lines(c(0,1),c(0,1))		
	}
}

#########################################################################
# Function: bum.histogram
# Purpose: Compare the fitted bum curve to the histogram
# Arguments: pvalues - vector of p-values
#            a - a for bum curve, default = estimated MLE
#            lambda - lambda for bum curve, default = estimated MLE
#             main - primary title of plot, default = "Histogram"
#             xlab - label of x-axis, default = "p-value
#             ylab - label of y-axis, default = "Density"
#########################################################################

bum.histogram<-function(pvalues,a=NA,lambda=NA,main="Histogram",xlab="p-value",ylab="Density",nstartpts=0,starta=0,startlambda=0)

{
	if(is.na(a)||is.na(lambda))
	{
		MLE<-bum.mle(pvalues,nstartpts=nstartpts,starta=starta,startlambda=startlambda)
		a<-MLE$a
		lambda<-MLE$lambda
	}
	hist(pvalues,probability=T,main=main,xlab=xlab,ylab=ylab)
	x<-1:100/100
	lines(x,dbum(x,a,lambda),lwd=3)
}

#########################################################################
# Function: ext.pi
# Purpose: Extract the maximal uniform componet from a bum density
# Arguments: a - shape parameter of beta component of bum distribution
#            lambda - mixing parameter, component of bum that is uniform
# Returns: the proportion of the density that can be extracted as a uniform
#########################################################################

ext.pi<-function(a,lambda)

{
	lambda[(lambda>1)+(lambda<0)>0]<-NA
	a[(a>1)+(a<0)]<-NA
	return((lambda+(1-lambda)*a))
}

#########################################################################
# Function: bum.emp.post
# Purpose: Compute the upper bound of BUM based empirical Bayes posterior
#          probability of the alternative hypothesis
# Arguments: x - the point or vector of points at which to compute EB post
#            a - shape parameter of beta component of bum distribution
#            lambda - mixing parameter, proportion uniform in bum dist
# Returns: the upper bound of the BUM based empirical Bayes posterior at x
# Notes: Computes an upper bound because maximal extraction of uniform
#        is performed.
#########################################################################

bum.emp.post<-function(x,a,lambda)

{
   pi<-ext.pi(a,lambda)
   return((dbum(x,a,lambda)-pi)/dbum(x,a,lambda))
}

#########################################################################
# Function: bum.FDR
# Purpose: Compute an estimated upper bound for the false discovery rate
#          when significance is determined by comparing a p-value to
#          a threshold tau
# Arguments: tau - the threshold of comparison (point or vector)
#            a - shape parameter of beta component of bum distribution
#            lambda - mixing parameter, proportion uniform in bum dist
# Returns: Estimated upper bound of FDR
# Notes: FDR is the expected proportion of rejections that are false
#        positives, or Type I errors
#########################################################################

bum.FDR<-function(tau,a,lambda)

{
   pi<-ext.pi(a,lambda)
   return(tau*pi/pbum(tau,a,lambda))
}

#########################################################################
# Function: bum.false.positive
# Purpose: Compute an estimated upper bound for the proportion of all tests
#          that will be false positives when significance is determined
#          by comparing the p-value to a threshold tau
# Arguments: tau - the threshold of comparison (point or vector)
#            a - shape parameter of beta component of bum distribution
#            lambda - mixing parameter, proportion uniform in bum dist
# Returns: Estimated upper bound of proportion of all tests that are
#          false positives
#########################################################################

bum.false.positive<-function(tau,a,lambda)

{
   pi<-ext.pi(a,lambda)
   return(tau*pi)
}

#########################################################################
# Function: bum.false.negative
# Purpose: Compute an estimated lower bound for the proportion of all
#          tests resulting in false negatives when significance is
#          determined by comparing the p-value to a threshold tau
# Arguments: tau - threshold
#            a - shape parameter of beta component of bum distribution
#            lambda - mixing parameter, proportion uniform in bum dist
# Returns: Estimated lower bound for proportion of all tests resulting
#          false negatives (Type II errors)
#########################################################################

bum.false.negative<-function(tau,a,lambda)

{
   pi<-ext.pi(a,lambda)
   return(1-pbum(tau,a,lambda)-pi*(1-tau))
}

#########################################################################
# Function: bum.true.positive
# Purpose: Compute an estimated lower bound for the proportion of all
#          tests resulting in true positives when significance is
#          determined by comparing the p-value to a threshold tau
# Arguments: tau - threshold
#            a - shape parameter of beta component of bum distribution
#            lambda - mixing parameter, proportion uniform in bum dist
# Returns: Estimated lower bound for proportion of all tests resulting
#          true positives
#########################################################################

bum.true.positive<-function(tau,a,lambda)

{
   pi<-ext.pi(a,lambda)
   return(pbum(tau,a,lambda)-pi*tau)
}

#########################################################################
# Function: bum.true.negative
# Purpose: Compute an estimated lower bound for the proportion of all
#          tests resulting in true negatives when significance is
#          determined by comparing the p-value to a threshold tau
# Arguments: tau - threshold
#            a - shape parameter of beta component of bum distribution
#            lambda - mixing parameter, proportion uniform in bum dist
# Returns: Estimated lower bound for proportion of all tests resulting
#          true negatives
#########################################################################

bum.true.negative<-function(tau,a,lambda)


{
   pi<-ext.pi(a,lambda)
   return(pi*(1-tau))
}


#########################################################################
# Function: bum.weighted.error
# Purpose: Compute a linear combination of the upper bound for the 
#          proportion of false positives and the lower bound for the
#          proportion of false negatives resulting when significance
#          is determined by comparing p-values to a threshold 
# Arguments: tau - threshold
#            a - shape parameter of beta component of bum distribution
#            lambda - mixing parameter, proportion uniform in bum dist
#            wfp - weight of the false positives
#            wfn - weight of the false negatives
# Returns: the weighted combination of error rates
#########################################################################

bum.weighted.error<-function(tau,a,lambda,wfp=1,wfn=1)

{
   return(wfp*bum.false.positive(tau,a,lambda)+wfn*bum.false.negative(tau,a,lambda))
}

#########################################################################
# Function: find.FDR.threshold
# Purpose: Find a p-value threshold so that all p-values less than the
#          threshold have an FDR lower than a desired level
# Arguments: fdr - desired fdr
#            a - shape parameter of beta component of bum distribution
#            lambda - mixing parameter, proportion uniform in bum dist
# Returns: p-value threshold that maintains the desired fdr
#########################################################################

find.FDR.threshold<-function(fdr,a,lambda)

{
   pi<-ext.pi(a,lambda)
   return(((pi-fdr*lambda)/(fdr*(1-lambda)))^(1/(a-1)))
}



#########################################################################
# Function: find.EB.threshold
# Purpose: Find a p-value threshold so that all p-values less than the
#          threshold have an empirical Bayes (EB) posterior probability
#          greater than a desired level
# Arguments: EB - desired empirical Bayes' posterior
#            a - shape parameter of beta component of bum distribution
#            lambda - mixing parameter, proportion uniform in bum dist
# Returns: p-value threshold that maintains the desired EB posterior
#########################################################################

find.EB.threshold<-function(EB,a,lambda)

{
   pi<-ext.pi(a,lambda)
   return(((EB*lambda+a*(1-lambda))/(a*(1-EB)*(1-lambda)))^(1/(a-1)))
}

##########################################################################
# Function: inv.dbum
# Purpose: Compute the inverse of the pdf  of the bum distribution
# Arguments: y - value of the pdf of the bum of interest
#            a - shape parameter of the beta component
#            lambda - mixing parameter, weight of uniform component
#            nbisect - the number of bisections to perform
# Returns: x so that pbum(x,a,lambda) = y
# Notes: Uses the bisection method.  The results will
#        be accurate to within 2^(-nbisect).  Used by find.WE.threshold.
##########################################################################

inv.dbum<-function(y,a,lambda,nbisect=20)

{

	n<-length(y)
	mid<-rep(0,n)
	for (i in 1:n)
	{
		top<-1
		bot<-0
		for (j in 1:nbisect)
		{
			mid[i]<-mean(c(top,bot))
			if(dbum(mid[i],a,lambda)>y[i]) bot<-mid[i]
			else top<-mid[i]	
		}
	}
	return(mid)
}


#=======================================================================
# Estimate pFDR and FDR using P values and V-D spline cdf/pdf estimators
# Return: list of 
#   PvalFDR -- the evaluation of pFDR and FDR on the alpha.seq grid
#        mx3 matrix, alpha.seq in column 1, pFDR in col 2, FDR in col 3
#   pi0 -- estimate of pi_0
#   cdf -- Pval distr cdf evaluatd on alph.grid
#   pdf -- Pval distr pdf evaluatd on alph..grid
#   alph.grid
#   Pvalues -- Pvals
#   side, error
#=======================================================================
Ch04.1 <- function(Pvals,qpos=c(0.00625,0.01,0.0125,0.025,0.05,0.1,0.25,0.5,0.75),notes=NULL)
{  
  G <- length(Pvals)
  Pgrid <- alph.grid <- alpha.seq <- Pvals
  nPg <- length(Pgrid)
  Z <- VDS.unif.est.Qknots(Pgrid,Pvals,qpos=qpos,degree=4,print=F)
  #pi0 <- min(Z[,2])
  pi0 <- min(1,Z[nPg,2],na.rm=T)
  cdf <- Z[,1]
  pdf <- Z[,2]

# Estimate pFDR and FDR over a grid of alpha values:
  Z <- VDS.unif.est.Qknots(alpha.seq,Pvals,qpos=qpos,degree=4,print=F)
  pFDR <- (pi0*alpha.seq)/(Z[,1]*(1-(1-alpha.seq)^G))
  FDR <- (pi0*alpha.seq)/Z[,1]
  
  ord<-order(Pvals)
  q<-rep(1,G)
  q[ord[G]]<-min(1,FDR[ord[G]])
  q[rev(ord)]<-cummin(FDR[rev(ord)])
  return(list(p=matrix(Pvals,G,1),
              fdr=matrix(FDR,G,1),
              q=matrix(q,G,1),
              cdf=matrix(cdf,G,1),
              pdf=matrix(pdf,G,1),
              ord=matrix(ord,G,1),
              pi=pi0,
              notes=c(notes,"Cheng et al (2004) used to obtain FDR estimates.")))
}

#========================================================
# The eimprical distribution function
# Input: X -- Data (missing values are deleted)
#        x1 -- grid over which to evaluate the EDF
# Output: vector of EDF values over x1
#========================================================
EDF <- function(X,x1) {  
  Xodr <- sort(X)
  n <- length(Xodr)
  if(n==0) return(NA)
  if(n==1) {  
    y <- approx(c(0.5*Xodr[1],Xodr),c(0,1),x1,
           method="constant",rule=2,f=0)
    return(y[[2]])
  }
  LB <- Xodr[1]-(Xodr[2]-Xodr[1])
  uniqueXodr<-unique(Xodr)
  z <- cumsum(table(Xodr))/n
  y <- approx(c(LB,uniqueXodr),c(0,z),x1,method="constant",
         rule=2,f=0)
  return(y[[2]])
}

#========================================================

#==================================================================
# Compute a B-spline approximation of a continuous function over a
# closed interval [a,b]. 
# Input:
#   x -- a grid over which to evaluate the approximation 
#   fv -- a vector of sampled function values
#   knots.int -- the interior knot sequence. 
#            m=length(knots).
#   degree -- spline degree (default 3). degree=order-1
#   a -- left boundary of interval (default=0)
#   b -- right boundary of interval (default=1)
#   The length of fv must be equal to m+degree+1, i.e., the number of 
#   interior knots plus order of the spline; otherwise NULL is returned.
# Output: a vector of the spline values over x.
#====================================================================
Bspline <- function(x,fv,knots.int,degree=3,a=0,b=1) {  
  m <- length(knots.int)
  nf <- length(fv)
  if(nf!=m+degree+1) return(NULL)
  Smat <- bsmat(x,knots.int=knots.int,degree=degree,
            a=a,b=b)
  func.val <- Smat%*%fv
  return(func.val)
}

#=====================================================================
#  Compute B-spline matrix evaluated overa grid of values on a closed 
#  interval [a,b].
# Input:
#   x -- grid
#   knots.int -- the interior knot sequance
#   degree -- spline degree (default 3). degree=order-1
#   a -- left boundary of interval (default=0)
#   b -- right boundary of interval (default=1)
#======================================================================
bsmat <- function(x,knots.int,degree=3,a=0,b=1) {  
# To be comformal with de Boor's notation:
  k <- degree+1  # order of spline
# Number of interior knots:
  m <- length(knots.int)
# Number of Splines:
  n <- m+k
  d <- length(x)
  smat1 <- smat2 <- matrix(0,d,n)
# The initial (i.e., k=1) B spline matrix
  kt <- unique(sort(c(a,knots.int,b))) # initial extended knot seq
  el <- 1; nel <- m+el
  indx <- numeric(0)
  for(j in 1:nel) {  
    indx <- (1:d)[x>=kt[j]&x<kt[j+1]]
    smat1[indx,j] <- 1
  }
# Make the last piece (order-1 B spline) left-continuous at b:
  indx <- (1:d)[x==b]
  smat1[indx,nel] <- 1
  if(k==1) return(smat1)
# Start the iteration, done if el=k:
  while(el<k) {  
    el <- el+1
    nel <- m+el
    kt <- c(a,kt,b)
    for(j in 2:(nel-1)) {  
      smat2[,j] <- ((x-kt[j])/(kt[j+el-1]-kt[j]))*smat1[,j-1]+
                   ((kt[j+el]-x)/(kt[j+el]-kt[j+1]))*smat1[,j]
    }
    smat2[,1] <- ((kt[1+el]-x)/(kt[1+el]-kt[2]))*smat1[,1]
    smat2[,nel] <- ((x-kt[nel])/(kt[nel+el-1]-kt[nel]))*smat1[,nel-1]
    smat1 <- smat2
  }
  return(smat1)
}

#=====================================================================

#==================================================================
# Compute variation-dimishing spline estimators of cdf and pdf on a 
# closed interval [0,1] over a grid of values, under the constraint 
# that the population distr. is stochastically smaller than the U(0,1). 
# The knot sequence is determined by linearly interpolated sample 
# quantiles at positions specified by user.
# Input:
#   x -- grid over which to evaluate the estimators
#   X.data -- vector of data based on which to compute the estimators
#   NOTE: The range of X.data must be in [0,1] in order to get 
#         meaningful estimates.
#   qpos -- quantile pos for knot generation, default (1:16)/17
#   a -- left boundary of interval, default 0
#   b -- right boundary of interval, default 1
#   degree -- degree of B spline for the cdf estimator, default 3
# Output: a two-column matrix; first column contains the cdf estimate
#         over x, 2nd column contains the pdf estimate over x 
#====================================================================
VDS.unif.est.Qknots <- function(x,X.data,qpos=(1:16)/17,a=0,b=1,degree=3,
                          print=T) 
{  
# Get unique interior knots:
# i.knots <- unique(pmin(quantile(X.data,qpos),qpos))
  i.knots <- unique(quantile(X.data,qpos))
  m <- length(i.knots)
# make sure that the interior knots don't touch the boundary:
  indx <- max((1:m)[i.knots<=a])
  if((!is.na(indx))&&(is.finite(indx))) {i.knos <- i.knots[(indx+1):m]; m <- m-indx}
  indx <- min((1:m)[i.knots>=b])
  if((!is.na(indx))&&(is.finite(indx))) {i.knots <- i.knots[1:(indx-1)]; m <- indx-1}
  m <- length(i.knots)
  knots <- unique(c(a,i.knots,b))
  nknots <- length(knots)
# To be conformal with de Boor's notation:
  k <- degree+1 
  n <- nknots+k-2  # so n = number of unique interior knots + order
  knots.ext <- c(rep(a,k-1),knots,rep(b,k-1))
  ts <- apply(matrix(1:n,n,1),1,
          function(i,tk,k) {  
            return(sum(tk[(i+1):(i+k-1)])/(k-1))
          },knots.ext,k
        )
# set the empirical cdf above the unif(0,1) cdf for the estimator
# to obey the stochastic order constraint:
  Ecdf <- pmax(EDF(X.data,ts),ts)
# Ecdf <- EDF(X.data,ts)
# lambn <- (ts[n-1]-ts[n-2])/(ts[n]-ts[n-2])
# Ecdf[n-1] <- lambn*Ecdf[n]+(1-lambn)*Ecdf[n-2]
  Epdf <- (k-1)*(Ecdf[2:n]-Ecdf[1:(n-1)])/
                (knots.ext[(2:n)+k-1]-knots.ext[2:n])
  if(print) {  
    cat("knots =",knots,"\n") 
    cat("knots.ext =",knots.ext,"\n")
    cat("ts =",ts,"\n")
    cat("Ecdf =",Ecdf,"\n")
    cat("Epdf =",Epdf,"\n")
  }
  VDcdf <- Bspline(x,Ecdf,knots[2:(nknots-1)],degree)
  VDpdf <- Bspline(x,Epdf,knots[2:(nknots-1)],degree-1)
  return(cbind(VDcdf,VDpdf))
}
