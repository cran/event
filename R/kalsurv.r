#
#  event : A Library of Special Functions for Event Histories
#  Copyright (C) 1998 J.K. Lindsey
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#  SYNOPSIS
#
#  kalsurv(response, intensity="exponential", distribution="Pareto",
#	depend="independence", update="Markov", mu=NULL, shape=NULL,
#	renewal=T, density=F, censor=NULL, delta=NULL, ccov=NULL, tvcov=NULL,
#	preg=NULL, ptvc=NULL, pbirth=NULL, pintercept=NULL, pshape=NULL,
#	pinitial=1, pdepend=NULL, envir=sys.frame(sys.parent()),
#	print.level=0, ndigit=10, gradtol=0.00001, steptol=0.00001,
#	iterlim=100, fscale=1, typsiz=abs(p), stepmax=10*sqrt(p%*%p))
#
#  DESCRIPTION
#
#    Function to fit various distributions inserted in a Pareto, gamma, or
# Weibull distribution with serial dependence or gamma frailties using
# Kalman-type update for event histories.

kalsurv <- function(response, intensity="exponential", distribution="Pareto",
	depend="independence", update="Markov", mu=NULL, shape=NULL,
	renewal=T, density=F, censor=NULL, delta=NULL, ccov=NULL, tvcov=NULL,
	preg=NULL, ptvc=NULL, pbirth=NULL, pintercept=NULL, pshape=NULL,
	pinitial=1, pdepend=NULL, envir=sys.frame(sys.parent()),
	print.level=0, ndigit=10, gradtol=0.00001, steptol=0.00001,
	iterlim=100, fscale=1, typsiz=abs(p), stepmax=10*sqrt(p%*%p)){
ksurvb <- function(p){
	if(rf)b <- mu1(p)
	if(sf)v <- sh1(p[nps1:np])
	z <- .C("ksurvb",
		p=as.double(p),
		y=as.double(resp$response$y),
		x=as.double(resp$ccov$ccov),
		cens=as.integer(resp$response$censor),
		nind=as.integer(nind),
		nobs=as.integer(resp$response$nobs),
		nbs=as.integer(length(resp$response$y)),
		nccov=as.integer(nccov),
		model=as.integer(mdl),
		density=as.integer(density),
		dep=as.integer(dep),
		birth=as.integer(birth),
		tvc=as.integer(tvc),
		tvcov=as.double(resp$tvcov$tvcov),
		fit=as.integer(0),
		pred=double(length(resp$response$y)),
		rpred=double(length(resp$response$y)),
		renewal=as.integer(renewal),
		rf=as.integer(rf),
		bb=as.double(b),
		sf=as.integer(sf),
		vv=as.double(v),
		like=double(1),
		DUP=F)
	z$like}
ksurvg <- function(p){
	if(rf)b <- mu1(p)
	if(sf)v <- sh1(p[nps1:np])
	z <- .C("ksurvg",
		p=as.double(p),
		y=as.double(resp$response$y),
		x=as.double(resp$ccov$ccov),
		cens=as.integer(resp$response$censor),
		nind=as.integer(nind),
		nobs=as.integer(resp$response$nobs),
		nbs=as.integer(length(resp$response$y)),
		nccov=as.integer(nccov),
		model=as.integer(mdl),
		distribution=as.integer(dst),
		density=as.integer(density),
		dep=as.integer(dep),
		birth=as.integer(birth),
		tvc=as.integer(tvc),
		tvcov=as.double(resp$tvcov$tvcov),
		fit=as.integer(0),
		pred=double(length(resp$response$y)),
		rpred=double(length(resp$response$y)),
		renewal=as.integer(renewal),
		rf=as.integer(rf),
		bb=as.double(b),
		sf=as.integer(sf),
		vv=as.double(v),
		like=double(1),
		DUP=F)
	z$like}
frailb <- function(p){
	if(rf)b <- mu1(p)
	if(sf)v <- sh1(p[nps1:np])
	z <- .C("frailb",
		p=as.double(p),
		y=as.double(resp$response$y),
		x=as.double(resp$ccov$ccov),
		cens=as.integer(resp$response$censor),
		nind=as.integer(nind),
		nobs=as.integer(resp$response$nobs),
		nbs=as.integer(length(resp$response$y)),
		nccov=as.integer(nccov),
		model=as.integer(mdl),
		density=as.integer(density),
		dep=as.integer(dep),
		birth=as.integer(birth),
		tvc=as.integer(tvc),
		tvcov=as.double(resp$tvcov$tvcov),
		fit=as.integer(0),
		pred=double(length(resp$response$y)),
		rpred=double(length(resp$response$y)),
		rf=as.integer(rf),
		bb=as.double(b),
		sf=as.integer(sf),
		vv=as.double(v),
		like=double(1),
		DUP=F)
	z$like}
call <- sys.call()
tmp <- c("exponential","Weibull","gamma","gen logistic",
	"log normal","log logistic","log Cauchy","log Laplace")
mdl <- match(intensity <- match.arg(intensity,tmp),tmp)
tmp <- c("Pareto","gamma","Weibull")
dst <- match(distribution <- match.arg(distribution,tmp),tmp)
depend <- match.arg(depend,c("independence","serial","frailty"))
tmp <- c("elapsed Markov","serial","Markov","event","cumulated","count","kalman","time")
dep <- match(update <- match.arg(update,tmp),tmp)
rf <- !is.null(mu)
sf <- !is.null(shape)
respenv <- inherits(response,"repeated")
envname <- if(respenv)paste(deparse(substitute(response)))
	else NULL
if(!respenv){
	if(!inherits(response,"response")){
		if(is.null(censor)){
			if(is.matrix(response)||is.data.frame(response))
				censor <- rep(1,nrow(response))
			else if(is.list(response))censor <- rep(1,length(response))}
		resp <- restovec(response,censor=censor,delta=delta)}
	else resp <- response
	if(is.null(ccov))nccov <- 0
	else {
		if(!inherits(ccov,"tccov")){
			ccname <- paste(deparse(substitute(ccov)))
			if((is.matrix(ccov)&&is.null(colnames(ccov)))){
				ccname <- paste(deparse(substitute(ccov)))
				if(ncol(ccov)>1){
					tmp <- NULL
					for(i in 1:ncol(ccov))tmp <- c(tmp,paste(ccname,i,sep=""))
					ccname <- tmp}}
			ccov <- tcctomat(ccov,names=ccname)}
		nccov <- if(rf) 0 else ncol(ccov$ccov)}
	if(is.null(tvcov))ttvc <- 0
	else {
		if(!inherits(tvcov,"tvcov")){
			tvcname <- paste(deparse(substitute(tvcov)))
			if(is.list(tvcov)&&ncol(tvcov[[1]])>1){
				if(is.null(colnames(tvcov[[1]]))){
					tvcname <- paste(deparse(substitute(tvcov)))
					tmp <- NULL
					for(i in 1:ncol(tvcov[[1]]))tmp <- c(tmp,paste(tvcname,i,sep=""))
					tvcname <- tmp}
				else tvcname <- colnames(tvcov[[1]])}
			tvcov <- tvctomat(tvcov, names=tvcname)}
		ttvc <- if(rf) 0 else ncol(tvcov$tvcov)}
	resp <- rmna(response=resp, tvcov=tvcov, ccov=ccov)
	if(!is.null(ccov))rm(ccov)
	if(!is.null(tvcov))rm(tvcov)}
else {
	if(!rf){
		resp <- response
		if(is.null(ccov))resp$ccov <- NULL
		else if(inherits(ccov,"formula"))
			resp$ccov$ccov <- attr(finterp(ccov,envir=response,expand=F,name=paste(deparse(substitute(response)))),"model")[,-1,drop=F]
		else stop("ccov must be a W&R formula")
		if(is.null(tvcov))resp$tvcov <- NULL
		else if(inherits(tvcov,"formula"))
			resp$tvcov$tvcov <- attr(finterp(tvcov,envir=response,name=paste(deparse(substitute(response)))),"model")[,-1,drop=F]
		else stop("tvcov must be a W&R formula")}
	else resp <- rmna(response$response)
	nccov <- if(rf||is.null(resp$ccov$ccov)) 0
		 else  ncol(resp$ccov$ccov)
	ttvc <- if(rf||is.null(resp$tvcov$tvcov)) 0
		 else  ncol(resp$tvcov$tvcov)}
if((inherits(envir,"repeated")&&
	(length(resp$response$nobs)!=length(envir$response$nobs)||
	any(resp$response$nobs!=envir$response$nobs)))||
	(inherits(envir,"tvcov")&&
	(length(resp$response$nobs)!=length(envir$tvcov$nobs)||
	any(resp$response$nobs!=envir$tvcov$nobs))))
	stop("response and envir objects are incompatible")
mu3 <- sh3 <- NULL
if(respenv||inherits(envir,"repeated")||inherits(envir,"tccov")){
	type <- if(respenv||inherits(envir,"repeated"))"repeated"
		else if(inherits(envir,"tccov"))"tccov"
		else "tvcov"
	if(is.null(envname))envname <- paste(deparse(substitute(envir)))
	if(inherits(mu,"formula")){
		mu3 <- if(respenv)finterp(mu,envir=response,name=envname)
			else finterp(mu,envir=envir,name=envname)
		class(mu) <- c(class(mu),type)}
	else if(is.function(mu)){
		tmp <- parse(text=paste(deparse(mu))[-1])
		class(mu) <- type
		mu <- if(respenv)fnenvir(mu,envir=response,name=envname)
			else fnenvir(mu,envir=envir,name=envname)
		mu3 <- mu
		if(respenv)attr(mu3,"model") <- tmp}
	if(inherits(shape,"formula")){
		sh3 <- if(respenv)finterp(shape,envir=response,name=envname)
			else finterp(shape,envir=envir,name=envname)
		class(shape) <- c(class(shape),type)}
	else if(is.function(shape)){
		tmp <- parse(text=paste(deparse(shape))[-1])
		class(shape) <- type
		shape <- if(respenv)fnenvir(shape,envir=response,name=envname)
			else fnenvir(shape,envir=envir,name=envname)
		sh3 <- shape
		if(respenv)attr(sh3,"model") <- tmp}}
npreg <- length(preg)
mu1 <- sh1 <- v <- b <- NULL
if(inherits(mu,"formula")){
	pr <- if(npreg>0)preg else ptvc
	npr <- length(pr)
	mu2 <- if(respenv)finterp(mu,envir=response,name=envname,expand=is.null(preg))
		else finterp(mu,envir=envir,name=envname,expand=is.null(preg))
	npt1 <- length(attr(mu2,"parameters"))
	if(is.matrix(attr(mu2,"model"))){
		if(all(dim(attr(mu2,"model"))==1)){
			mu1 <- function(p) p[1]*rep(1,n)
			attributes(mu1) <- attributes(mu2)
			mu2 <- NULL}}
	else {
		if(npr!=npt1&&length(ptvc)!=npt1){
			cat("\nParameters are ")
			cat(attr(mu2,"parameters"),"\n")
			stop(paste("preg or ptvc should have",npt1,"estimates"))}
		if(is.list(pr)){
			if(!is.null(names(pr))){
				o <- match(attr(mu2,"parameters"),names(pr))
				pr <- unlist(pr)[o]
				if(sum(!is.na(o))!=length(pr))stop("invalid estimates for mu - probably wrong names")}
			else pr <- unlist(pr)
			if(npreg>0)preg <- pr else ptvc <- pr}}
	if(!is.null(mu2)){
		if(inherits(envir,"tccov")){
			cv <- covind(response)
			mu1 <- function(p) mu2(p)[cv]
			attributes(mu1) <- attributes(mu2)}
		else {
			mu1 <- mu2
			rm(mu2)}}}
else if(is.function(mu))mu1 <- mu
if(!is.null(mu1)&&is.null(attributes(mu1))){
	attributes(mu1) <- if(is.function(mu)){
		if(!inherits(mu,"formulafn")){
			if(respenv)attributes(fnenvir(mu,envir=response))
			else attributes(fnenvir(mu,envir=envir))}
		else attributes(mu)}
		else {
			if(respenv)attributes(fnenvir(mu1,envir=response))
			else attributes(fnenvir(mu1,envir=envir))}}
nlp <- if(is.function(mu1))length(attr(mu1,"parameters"))
	else if(is.null(mu1))NULL
	else npt1
if(!is.null(nlp)){
	if(is.null(ptvc)&&nlp!=npreg)
		stop(paste("preg should have",nlp,"initial estimates"))
	else if(!is.null(ptvc)&&length(ptvc)!=nlp)
		stop(paste("ptvc should have",nlp,"initial estimates"))}
nps <- length(pshape)
if(inherits(shape,"formula")){
	sh2 <- if(respenv)finterp(shape,envir=response,name=envname)
		else finterp(shape,envir=envir,name=envname)
	npt2 <- length(attr(sh2,"parameters"))
	if(is.matrix(attr(sh2,"model"))){
		if(all(dim(attr(sh2,"model"))==1)){
			sh1 <- function(p) p[npl1]*rep(1,n)
			attributes(sh1) <- attributes(sh2)
			sh2 <- NULL}}
	else {
		if(nps!=npt2){
			cat("\nParameters are ")
			cat(attr(sh2,"parameters"),"\n")
			stop(paste("pshape should have",npt2,"estimates"))}
		if(is.list(pshape)){
			if(!is.null(names(pshape))){
				o <- match(attr(sh2,"parameters"),names(pshape))
				pshape <- unlist(pshape)[o]
				if(sum(!is.na(o))!=length(pshape))stop("invalid estimates for shape - probably wrong names")}
			else pshape <- unlist(pshape)}}
	if(!is.null(sh2)){
		if(inherits(envir,"tccov")){
			cv <- covind(response)
			sh1 <- function(p) sh2(p)[cv]
			attributes(sh1) <- attributes(sh2)}
		else {
			sh1 <- sh2
			rm(sh2)}}}
else if(is.function(shape))sh1 <- shape
if(!is.null(sh1)&&is.null(attributes(sh1)))
	attributes(sh1) <- if(is.function(shape)){
		if(!inherits(shape,"formulafn")){
			if(respenv)attributes(fnenvir(shape,envir=response))
			else attributes(fnenvir(shape,envir=envir))}
		else attributes(shape)}
		else {
			if(respenv)attributes(fnenvir(sh1,envir=response))
			else attributes(fnenvir(sh1,envir=envir))}
nlp <- if(is.function(shape))length(attr(sh1,"parameters"))
	else if(is.null(shape))NULL
	else npt2
if(!is.null(nlp)&&nlp!=nps)
	stop(paste("pshape should have",nlp,"initial estimates"))
if(rf&&!is.function(mu1))stop("mu must be a formula or function")
if(sf&&!is.function(sh1))stop("shape must be a formula or function")
birth <- !is.null(pbirth)
tvc <- length(ptvc)
if(rf&&birth)stop("Birth models cannot be fitted with a mean function")
if(intensity=="exponential"){
	sf <- F
	pshape <- NULL}
else {
	if(is.null(pshape))
		stop("Initial value of the shape parameter must be supplied")
	if(!sf){
		if(pshape<=0)stop("Shape must be positive")
		else pshape <- log(pshape)}}
if(intensity=="gen logistic"){
	if(is.null(pintercept))stop("Initial value of the intercept parameter must be supplied")}
else pintercept <- NULL
np <- 1+npreg+(depend=="serial")+birth+nps+!is.null(pintercept)
if(pinitial<=0)stop("Estimate of initial parameter must be > 0")
else pinitial <- log(pinitial)
if(depend=="independence"){
	pdepend <- NULL
	dep <- 0}
else if(depend=="serial"){
	if(update=="time")stop("time update can only be used with frailty")
	if(is.null(pdepend))
		stop("An estimate of the dependence parameter must be supplied")
	else if(pdepend<=0|pdepend>=1)
		stop("Dependence parameter must be between 0 and 1")
	else pdepend <- log(pdepend/(1-pdepend))}
else if(depend=="frailty"){
	if(update=="time")dep <- 1
	else {
		dep <- 0
		update <- "no"}
	if(!is.null(pdepend))pdepend <- NULL}
if(is.null(resp$response$censor))
	resp$response$censor <- rep(1,length(resp$response$y))
if(rf&&npreg>0)nccov <- npreg-1
if(!rf&&nccov+1!=npreg)
	stop(paste(nccov+1,"regression estimates must be supplied"))
nind <- length(resp$response$nobs)
if(any(resp$response$y<0))stop("All times must be non-negative")
if(ttvc>0)np <- np+tvc
if(!rf&&(ttvc>0&&tvc!=ttvc||ttvc==0&&tvc>0))stop(paste(ttvc,"initial estimates of coefficients for time-varying covariates must be supplied"))
if(rf){
	if(tvc>0&&nccov>0)stop("With a mean function, initial estimates must be supplied either in preg or in ptvc")
	if(tvc>0){
		if(length(mu1(ptvc))!=length(resp$response$y))stop("The mu function must provide an estimate for each observation")
		tvc <- tvc-1
		np <- np+tvc}
	else if(length(mu1(preg))==1){
		if(nccov==0)mu1 <- function(p) rep(p[1],length(resp$response$y))
		else stop("Number of estimates does not correspond to mu function")}
	else if(length(mu1(preg))!=nind)stop("The mu function must provide an estimate for each individual")}
if(sf&&length(sh1(pshape))!=length(resp$response$y))stop("The shape function must provide an estimate for each observation")
nps1 <- np-nps-(!is.null(pintercept))+1
p <- c(preg,pbirth,ptvc,pinitial,pdepend,pshape,pintercept)
if(distribution=="Pareto"){
	if(depend=="frailty")surv <- frailb
	else surv <- ksurvb}
else if(distribution=="gamma"||distribution=="Weibull"){
	if(depend=="frailty")surv <- frailg
	else surv <- ksurvg}
if(fscale==1)fscale <- surv(p)
if(is.na(surv(p)))
	stop("Likelihood returns NAs: probably invalid initial values")
z0 <- nlm(surv, p=p, hessian=T, print.level=print.level,
	typsiz=typsiz, ndigit=ndigit, gradtol=gradtol, stepmax=stepmax,
	steptol=steptol, iterlim=iterlim, fscale=fscale)
p <- z0$estimate
like <- z0$minimum
if(!is.null(resp$response$delta))like <- like-
	if(length(resp$response$delta)==1)
		length(resp$response$y)*log(resp$response$delta)
	else sum(log(resp$response$delta))
if(any(is.na(z0$hessian)))a <- 0
else a <- qr(z0$hessian)$rank
if(a==np)cov <- solve(z0$hessian)
else cov <- matrix(NA,ncol=np,nrow=np)
se <- sqrt(diag(cov))
corr <- cov/(se%o%se)
dimnames(corr) <- list(1:np,1:np)
if(mdl==4)z <- list()
else {
	z <- if(depend=="frailty"){
		if(rf)b <- mu1(p)
		if(sf)v <- sh1(p[nps1:np])
		z <- .C("frailb",
			p=as.double(p),
			y=as.double(resp$response$y),
			x=as.double(resp$ccov$ccov),
			cens=as.integer(resp$response$censor),
			nind=as.integer(nind),
			nobs=as.integer(resp$response$nobs),
			nbs=as.integer(length(resp$response$y)),
			nccov=as.integer(nccov),
			model=as.integer(mdl),
			density=as.integer(density),
			dep=as.integer(dep),
			birth=as.integer(birth),
			tvc=as.integer(tvc),
			tvcov=as.double(resp$tvcov$tvcov),
			fit=as.integer(1),
			pred=double(length(resp$response$y)),
			rpred=double(length(resp$response$y)),
			rf=as.integer(rf),
			bb=as.double(b),
			sf=as.integer(sf),
			vv=as.double(v),
			like=double(1),
			DUP=F)}
	else if(distribution=="Pareto"){
		if(rf)b <- mu1(p)
		if(sf)v <- sh1(p[nps1:np])
		z <- .C("ksurvb",
			p=as.double(p),
			y=as.double(resp$response$y),
			x=as.double(resp$ccov$ccov),
			cens=as.integer(resp$response$censor),
			nind=as.integer(nind),
			nobs=as.integer(resp$response$nobs),
			nbs=as.integer(length(resp$response$y)),
			nccov=as.integer(nccov),
			model=as.integer(mdl),
			density=as.integer(density),
			dep=as.integer(dep),
			birth=as.integer(birth),
			tvc=as.integer(tvc),
			tvcov=as.double(resp$tvcov$tvcov),
			fit=as.integer(1),
			pred=double(length(resp$response$y)),
			rpred=double(length(resp$response$y)),
			renewal=as.integer(renewal),
			rf=as.integer(rf),
			bb=as.double(b),
			sf=as.integer(sf),
			vv=as.double(v),
			like=double(1),
			DUP=F)}
	else {
		if(rf)b <- mu1(p)
		if(sf)v <- sh1(p[nps1:np])
		z <- .C("ksurvg",
			p=as.double(p),
			y=as.double(resp$response$y),
			x=as.double(resp$ccov$ccov),
			cens=as.integer(resp$response$censor),
			nind=as.integer(nind),
			nobs=as.integer(resp$response$nobs),
			nbs=as.integer(length(resp$response$y)),
			nccov=as.integer(nccov),
			model=as.integer(mdl),
			distribution=as.integer(dst),
			density=as.integer(density),
			dep=as.integer(dep),
			birth=as.integer(birth),
			tvc=as.integer(tvc),
			tvcov=as.double(resp$tvcov$tvcov),
			fit=as.integer(1),
			pred=double(length(resp$response$y)),
			rpred=double(length(resp$response$y)),
			renewal=as.integer(renewal),
			rf=as.integer(rf),
			bb=as.double(b),
			sf=as.integer(sf),
			vv=as.double(v),
			like=double(1),
			DUP=F)}
	if(mdl>4){
		z$pred <- exp(z$pred)
		z$rpred <- exp(z$rpred)}
	for(i in 1:length(resp$response$y))if(resp$response$y[i]==0){
		z$pred[i] <- z$pred[i-1]
		z$rpred[i] <- z$rpred[i-1]}}
if(!is.null(mu3))mu1 <- mu3
if(!is.null(sh3))sh1 <- sh3
z <- list(
	call=call,
	intensity=intensity,
	distribution=distribution,
	mu=mu1,
	npr=1+nccov+tvc+birth,
	shape=sh1,
	nps=np-nps,
	density=density,
	depend=depend,
	update=update,
	birth=birth,
	renewal=renewal,
	response=resp$response,
	pred=z$pred,
	rpred=z$rpred,
	ccov=resp$ccov,
	tvcov=resp$tvcov,
	maxlike=like,
	aic=like+np,
	df=length(resp$response$y)-np,
	npt=np,
	npv=npreg,
	coefficients=p,
	se=se,
	cov=cov,
	corr=corr,
	grad=z0$gradient,
	iterations=z0$iterations,
	code=z0$code)
class(z) <- if(mdl==4)"kalsurv" else c("kalsurv","recursive")
return(z)}

coefficients.kalsurv <- function(z) z$coefficients
deviance.kalsurv <- function(z) 2*z$maxlike
fitted.kalsurv <- function(z, recursive=TRUE)
	if(recursive) z$rpred else z$pred
residuals.kalsurv <- function(z, recursive=TRUE)
	if(recursive) z$response$y-z$rpred else z$response$y-z$pred

print.kalsurv <- function(z, digits = max(3, .Options$digits - 3)) {
	tvc <- !is.null(z$tvcov)
	expm <- z$intensity!="exponential"&&!is.function(z$shape)
	glm <- z$intensity=="gen logistic"
	nps <- if(is.function(z$shape)) z$nps else z$npt
	deppar <- (z$depend!="independence")&&(z$depend!="frailty")
	cat("\nCall:",deparse(z$call),sep="\n")
	cat("\n")
	if(z$code>2)cat("Warning: no convergence - error",z$code,"\n\n")
	cat("Number of subjects    ",length(z$resp$nobs),"\n")
	cat("Number of observations",length(z$response$y),"\n")
	cat(z$distribution,"distribution ")
	if(z$renewal){
		if(!z$birth)cat("with renewal process\n")
		else cat("with birth process\n")}
	else {
		cat("with zero origin\n")
		if(z$birth)cat("and birth process\n")}
	if(z$density)cat(z$intensity," density",sep="")
	else cat(z$intensity," intensity",sep="")
	if(z$depend=="independence")cat(" with independence\n")
	else if(z$depend=="frailty")
		cat(" with",z$depend,"dependence and",z$update,"weight\n")
	else cat(" with ",z$update," update\n",sep="")
	if(is.function(z$mu)){
		t <- deparse(z$mu)
		cat("Location function:",t[2:length(t)],sep="\n")}
	cat("\n-Log likelihood   ",z$maxlike,"\n")
	cat("Degrees of freedom",z$df,"\n")
	cat("AIC               ",z$aic,"\n")
	cat("Iterations        ",z$iterations,"\n\n")
	cat("Location parameters\n")
	if(!is.null(attr(z$mu,"formula")))
		cat(deparse(attr(z$mu,"formula")),sep="\n")
	else if(!is.null(attr(z$mu,"model"))){
		t <- deparse(attr(z$mu,"model"))
		t[1] <- sub("expression\\(","",t[1])
		t[length(t)] <- sub("\\)$","",t[length(t)])
		cat(t,sep="\n")}
	coef.table <- cbind(z$coef[1:z$npr],z$se[1:z$npr])
	if(inherits(z$mu,"formulafn"))
		cname <- if(is.matrix(attr(z$mu,"model")))
				colnames(attr(z$mu,"model"))
			else attr(z$mu,"parameters")
	else {
		cname <- "(Intercept)"
		if(z$npv>0)cname <- c(cname,colnames(z$ccov$ccov))
		if(z$birth)cname <- c(cname,"birth")
		if(tvc)cname <- c(cname,colnames(z$tvcov$tvcov))}
	dimnames(coef.table) <- list(cname, c("estimate","se"))
	print.default(coef.table, digits=digits, print.gap=2)
	if(is.function(z$shape))cat("\nDependence parameters\n")
	else cat("\nNonlinear parameters\n")
	coef <- exp(z$coef[(nps-deppar-expm-glm):nps])
	cname <- "initial"
	if(deppar){
		coef[2] <- coef[2]/(1+coef[2])
		cname <- c(cname,"depend")}
	if(glm){
		cname <- c(cname,"asymptote","intercept")
		coef[length(coef)-1] <- 1/coef[length(coef)-1]
		coef[length(coef)] <- NA}
	else if(expm)cname <- c(cname,"shape")
	coef.table <- cbind(z$coef[(nps-deppar-expm-glm):nps],z$se[(nps-deppar-expm-glm):nps],coef)
	dimnames(coef.table) <- list(cname, c("estimate","se","parameter"))
	print.default(coef.table, digits=digits, print.gap=2)
	if(z$depend=="frailty"){
		tmp <- trigamma(exp(-z$coef[nps-deppar-expm]))
		cat("Correlation =",tmp/(tmp+trigamma(1)),"\n")}
	if(is.function(z$sh1)||is.function(z$shape)){
		cat("\nShape parameters\n")
		if(!is.null(attr(z$sh1,"formula")))
			cat(deparse(attr(z$sh1,"formula")),sep="\n")
		else if(!is.null(attr(z$sh1,"model"))){
			t <- deparse(attr(z$sh1,"model"))
			t[1] <- sub("expression\\(","",t[1])
			t[length(t)] <- sub("\\)$","",t[length(t)])
			cat(t,sep="\n")}
		cname <- if(is.matrix(attr(z$sh1,"model")))
				colnames(attr(z$sh1,"model"))
			else attr(z$sh1,"parameters")
		coef.table <- cbind(z$coef[(z$nps+1):z$npt],
			z$se[(z$nps+1):z$npt])
		dimnames(coef.table) <- list(cname, c("estimate","se"))
		print.default(coef.table, digits=digits, print.gap=2)}
	cat("\nCorrelation matrix\n")
	print.default(z$corr, digits=digits)}
