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
#	pp(y, censor=1)
#	ident(y, id)
#	tpast(y)
#	ttime(y, id)
#	bp(y, id, censor=1)
#	tccov(y, x, id)
#	tvcov(y, x, tx)
#	vdm(y, x, id=NULL, tx=NULL, factor=F, time=F)
#	ehr(point, lambda=NULL, linear=NULL, plambda=NULL, wt=1,
#		envir=sys.frame(sys.parent()), print.level=0,
#		typsiz=rep(1,length(plambda)), ndigit=10,
#		gradtol=0.00001, stepmax=max(10*sqrt(plambda%*%plambda),10),
#		steptol=0.0004, iterlim=100, fscale=1)
#
#  DESCRIPTION
#
#    Functions for setting up and fitting counting process models

# point process created from times (y) between events
# y must contain integers
pp <- function(y, censor=1) {
	if(min(y)<=0)stop("All times must be positive")
	if(any(round(y)!=y))stop("Times must be integers")
	if(any(censor!=0&&censor!=1))
		stop("Censor indicator must be zeros and ones")
	if(length(censor)!=1&&length(censor)!=length(y))
		stop("Time and censor vectors must be the same length")
	point <- rep(0, sum(y))
	point[cumsum(y)] <- censor
	point}

# individual identification vector
ident <- function(y, id) {
	if(min(y)<=0)stop("All times must be positive")
	if(length(y)!=length(id))
		stop("Time and id vectors must be the same length")
	rep(id, y)}

# time past since previous event
tpast <- function(y) {
	if(min(y)<=0)stop("All times must be positive")
	unlist(lapply(as.list(y), seq))}

#	sequence(y)}
#sequence <- function(y) unlist(lapply(as.list(y), seq))

# total time elapsed for each individual
ttime <- function(y, id) {
	if(length(idd <- ident(y,id))==1)return(idd)
	z <- collapse(rep(1,length(idd)),idd,cumsum)
	names(z) <- NULL
	z}

# number of previous events for each individual, for birth processes
# add one if process starts at an event
bp <- function(y, id, censor=1) {
	bp1 <- function(i) c(0,cumsum(i)[1:(length(i)-1)])
	if(length(point <- pp(y, censor=censor))==1)return(point)
	if(length(idd <- ident(y, id))==1)return(idd)
	z <- collapse(point, idd, bp1)
	names(z) <- NULL
	z}

# time-constant covariate - id must be numbered consecutively
# x has one value for each distinct id
tccov <- function(y, x, id) {
	if(length(y)!=length(id))stop("Time and id must be the same length")
	if(length(x)!=length(unique(id)))
		stop("There must be one covariate value per individual")
	if(length(idd <- ident(y, id))==1)return(idd)
	x[idd]}

# time-varying covariate - tx gives the times at which x changes
# may also be used to create weight vector
tvcov <- function(y, x, tx) {
	if(min(y)<=0|min(tx)<0)stop("All times must be positive")
	if(length(x)!=length(tx))
		stop("Covariate and time vectors must be the same length")
	if(sum(y)!=sum(tx))
		stop("Total response time must equal total covariate time")
	rep(x, tx)}

# design matrix
vdm <- function(y, x, id=NULL, tx=NULL, factor=F, time=F) {
	if(time) {if(length(xx <- tvcov(y, x, tx))==1)return(xx)}
	else if(length(xx <- tccov(y, x, id))==1)return(xx)
	if(factor)xx <- factor(xx)
	wr(~xx)$design}

# fit an intensity function to event histories, where point is
# produced by point <- pp(y) and lambda is the log intensity function
ehr <- function(point, lambda=NULL, linear=NULL, plambda=NULL, delta=1,
	envir=sys.frame(sys.parent()), print.level=0,
	typsiz=rep(1,length(plambda)), ndigit=10, gradtol=0.00001,
	stepmax=max(10*sqrt(plambda%*%plambda),10), steptol=0.0004,
	iterlim=100, fscale=1){
call <- sys.call()
if(any(point<0))stop("Response vector must be non-negative integers")
n <- length(point)
dt <- any(delta>1)
if(dt){
	if(length(point)!=length(delta))stop("point and delta must be the same length")
	delta <- log(delta)}
if(inherits(lambda,"formula"))lin <- lambda
else if(inherits(linear,"formula"))lin <- linear
else lin <- NULL
lin1a <- lambda3 <- name <- NULL
if(inherits(envir,"repeated")||inherits(envir,"tccov")){
	type <- if(inherits(envir,"repeated"))"repeated"
		else if(inherits(envir,"tccov"))"tccov"
		else "tvcov"
	name <- paste(deparse(substitute(envir)))
	if(inherits(lin,"formula")){
		lin1a <- finterp(lin)
		class(lin) <- c(class(lin),type)}
	if(is.function(lambda)){
		lambda3 <- lambda
		attributes(lambda3) <- attributes(fnenvir(lambda))
		class(lambda) <- type
		lambda <- fnenvir(lambda,envir=envir,name=name)}}
npl <- length(plambda)
if(!is.null(lin)){
	lambda2 <- finterp(lin,envir=envir,name=name)
	npt1 <- length(attr(lambda2,"parameters"))
	if(is.matrix(attr(lambda2,"model"))){
		if(all(dim(attr(lambda2,"model"))==1)){
			if(is.function(lambda)){
				lin1a <- lambda2
				lambda1 <- if(dt)
					function(p) lambda(p,p[1]*rep(1,n))+delta
				else function(p) lambda(p,p[1]*rep(1,n))}
			else {
				lambda1 <- function(p) p[1]*rep(1,n)
				attributes(lambda1) <- attributes(lambda2)}}
		else {
			if(nrow(attr(lambda2,"model"))!=n)stop("lambda model matrix does not match number of response observations")
			if(is.function(lambda)){
				dm1 <- attr(lambda2,"model")
				lin1a <- lambda2
				lambda1 <- if(dt)function(p)
					function(p) lambda(p,dm1%*%p[1:npt1])+delta
				else function(p) lambda(p,dm1%*%p[1:npt1])
				attributes(lambda1) <- attributes(fnenvir(function(p) lambda(p,dm1%*%p[1:npt1])))}
			else {
				if(dt){
					dm1 <- attr(lambda2,"model")
					lambda1 <- function(p) dm1%*%p[1:npt1]+delta
					attributes(lambda1) <- attributes(lambda2)}
				else lambda1 <- lambda2}}}
	else {
		if(npl!=npt1){
			cat("\nParameters are ")
			cat(attr(lambda2,"parameters"),"\n")
			stop(paste("plambda should have",npt1,"estimates"))}
		if(dt){
			tmp <- attributes(lambda2)
			lambda1 <- function(p) lambda2(p)+delta
			attributes(lambda1) <- tmp}
		else lambda1 <- lambda2
		if(is.list(plambda)){
			if(!is.null(names(plambda))){
				o <- match(attr(lambda2,"parameters"),names(plambda))
				plambda <- unlist(plambda)[o]
				if(sum(!is.na(o))!=length(plambda))stop("invalid estimates for lambda - probably wrong names")}
			else plambda <- unlist(plambda)}}
	if(npl<npt1)stop("Not enough initial estimates for lambda")}
else if(!is.function(lambda)){
	lambda1 <- if(dt)function(p) p[1]*rep(1,n)+delta
		else function(p) p[1]*rep(1,n)
	npt1 <- 1}
else {
	if(dt){
		lambda3 <- fnenvir(lambda)
		lambda1 <- function(p) lambda(p)+delta}
	else lambda1 <- lambda}
if(is.null(attributes(lambda1))){
	attributes(lambda1) <- if(is.function(lambda)){
		if(!inherits(lambda,"formulafn"))attributes(fnenvir(lambda))
		else attributes(lambda)}
		else attributes(fnenvir(lambda1))}
nlp <- if(is.function(lambda)){
		if(is.null(lin))length(attr(lambda1,"parameters"))
		else length(attr(lambda1,"parameters"))-1+npt1}
       else npt1
if(nlp!=npl)stop(paste("plambda should have",nlp,"initial estimates"))
fn <- function(p) {
	l <- lambda1(p)
	sum(exp(l)-point*l)}
if(fscale==1)fscale <- fn(plambda)
if(is.na(fn(plambda)))
	stop("Likelihood returns NAs: probably invalid initial values")
z0 <- nlm(fn, p=plambda, hessian=T, print.level=print.level, typsiz=typsiz,
	ndigit=ndigit, gradtol=gradtol, stepmax=stepmax,
	steptol=steptol, iterlim=iterlim, fscale=fscale)
if(any(point>1))z0$minimum <- z0$minimum+sum(lgamma(point+1))
if(length(plambda)==1)cov <- 1/z0$hessian
else {
	a <- qr(z0$hessian)
	if(a$rank==length(plambda))cov <- solve(z0$hessian)
	else cov <- matrix(NA,ncol=length(plambda),nrow=length(plambda))}
se <- sqrt(diag(cov))
if(!is.null(lambda3))lambda1 <- lambda3
if(!is.null(lin1a))lin <- lin1a
z1 <- list(
	call=call,
	intensity=lambda1,
	linear=lin,
	maxlike=z0$minimum,
	aic=z0$minimum+length(plambda),
	coefficients=z0$estimate,
	se=se,
	cov=cov,
	corr=cov/(se%o%se),
	gradient=z0$gradient,
	iterations=z0$iter,
	error=z0$error,
	code=z0$code)
class(z1) <- "intensity"
return(z1)}

coefficients.intensity <- function(z) z$coefficients
deviance.intensity <- function(z) 2*z$maxlike

print.intensity <- function(z) {
	np <- length(z$coefficients)
	cat("\nCall:",deparse(z$call),sep="\n")
	cat("\n")
	if(z$code>2)cat("Warning: no convergence - error",z$code,"\n\n")
	cat("Log intensity function:\n")
	if(!is.null(attr(z$intensity,"formula")))
		cat(deparse(attr(z$intensity,"formula")),sep="\n")
	else if(!is.null(attr(z$intensity,"model"))){
		t <- deparse(attr(z$intensity,"model"))
		t[1] <- sub("expression\\(","",t[1])
		t[length(t)] <- sub("\\)$","",t[length(t)])
		cat(t,sep="\n")}
	if(inherits(z$linear,"formulafn"))
		cat("Linear part: ",deparse(attr(z$linear,"formula")),sep="\n")
	cat("\n-Log likelihood   ",z$maxlike,"\n")
	cat("AIC               ",z$aic,"\n")
	cat("Iterations        ",z$iterations,"\n\n")
	cat("Coefficients:\n")
	cname <- if(is.matrix(attr(z$intensity,"model")))colnames(attr(z$intensity,"model"))
		else if(length(grep("linear",attr(z$intensity,"parameters")))>0)
		attr(z$intensity,"parameters")[grep("\\[",attr(z$intensity,"parameters"))]
		else attr(z$intensity,"parameters")
	if(!is.null(z$linear)&&!is.null(attr(z$linear,"parameters")))
		cname <- c(colnames(attr(z$linear,"model")),cname)
	coef.table <- cbind(z$coefficients, z$se)
	dimnames(coef.table) <- list(cname, c("estimate", "se"))
	print.default(coef.table, digits=4, print.gap=2)
	if(np>1){
		cat("\nCorrelations:\n")
		dimnames(z$corr) <- list(seq(1,np),seq(1,np))
		print.default(z$corr, digits=4)}
	invisible(z)}

# examples of linear log intensity functions
#exponential <- ~1
#Weibull <- ~log(time(y))
#extreme.value <- ~time(y)
#birth1 <- ~bp(y,id)
#birth2 <- ~log(1+bp(y,id))

# examples of nonlinear log intensity functions
#negative.binomial <- function(p) p[1]+log(p[2]+bp(y,id))
#gen.negative.binomial <- function(p) p[1]+p[3]*log(p[2]+bp(y,id))
