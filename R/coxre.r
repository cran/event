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
#	coxre(response, censor, nest=NULL, cov=NULL, stratified=F, cumul=F,
#		estimate=1, iter=10, print.level=0, ndigit=10, gradtol=0.00001,
#		steptol=0.00001, iterlim=100, fscale=1, typsiz=abs(estimate),
#		stepmax=estimate)
#
#  DESCRIPTION
#
#    Cox model with frailites.
# glim code for coxre translated from Clayton, D. (1987) The analysis of
# event history data: a review of progress and outstanding problems.
# Statistics in Medicine 7: 819-841

coxre <- function(response, censor, nest=NULL, cov=NULL, stratified=F, cumul=F,
	estimate=1, iter=10, print.level=0, ndigit=10, gradtol=0.00001,
	steptol=0.00001, iterlim=100, fscale=1, typsiz=abs(estimate),
	stepmax=estimate){
like <- function(g){
	g <- exp(g)
	l <- 0
	for(i in 1:iter){
		fexp <- NULL
		for(ij in split(fit,nnest))fexp <- c(fexp,sum(ij))
		fexp <- fexp/rand
		ffex <- NULL
		for(ij in split(fit,interval))ffex <- c(ffex,sum(ij))
		ffex <- ffex/fixd
		tmp <- l
		l <- sum(counts*log(fit))-sum((fobs+1/g)*log(1+g*fexp)+fobs*log(rand))
		rand <- (1+g*fobs)/(1+g*fexp)
		fixd <- ffob/ffex
		oset <- log(rand[nnest]*fixd[interval])
		z <- glm(lmod,family=poisson,offset=oset)
		fit <- z$fit
		coef <- z$coef
		rm(z)
		if((l-tmp)^2<0.00002)break}
	l <- l+sum(lgamma(1/g+fobs)-lgamma(1/g)+fobs*log(g))
	names(fixd) <- 1:length(fixd)
	names(rand) <- 1:length(rand)
	list(like=-l,
		fixed=fixd,
		random=rand,
		inthaz=ffex,
		coefficients=coef)}
likel <- function(g){
	g <- exp(g)
	l <- 0
	for(i in 1:iter){
		fexp <- NULL
		for(ij in split(fit,nnest))fexp <- c(fexp,sum(ij))
		fexp <- fexp/rand
		ffex <- NULL
		for(ij in split(fit,interval))ffex <- c(ffex,sum(ij))
		ffex <- ffex/fixd
		tmp <- l
		l <- sum(counts*log(fit))-sum((fobs+1/g)*log(1+g*fexp)+fobs*log(rand))
		rand <- (1+g*fobs)/(1+g*fexp)
		fixd <- ffob/ffex
		oset <- log(rand[nnest]*fixd[interval])
		z <- glm(lmod,family=poisson,offset=oset)
		fit <- z$fit
		rm(z)
		if((l-tmp)^2<0.00002)break}
	l <- l+sum(lgamma(1/g+fobs)-lgamma(1/g)+fobs*log(g))
	-l}
call <- sys.call()
if(is.vector(response,mode="numeric")){
	nind1 <- length(response)
	nc <- 1
	if(is.vector(censor,mode="numeric")&&length(censor)==length(response))
		censor <- matrix(censor,ncol=1)
	else stop("censor must be a vector the same length as response")
	response <- matrix(response,ncol=1)
	stratified <- T}
else if(is.matrix(response)){
	nind1 <- nrow(response)
	if(!is.matrix(censor)||ncol(censor)!=ncol(response))
		stop("response and censor must have the same number of columns")
	if(!stratified){
		resp1 <- gl(ncol(response),nrow(response),length(response))
		nest <- as.integer(rep(nest,ncol(response)))
		if(!missing(cov))cov <- rep(cov,ncol(response))
		response <- matrix(response,ncol=1)
		censor <- matrix(censor,ncol=1)
		nc <- 1}
	else nc <- ncol(response)}
else stop("response must be a vector or matrix")
nind <- nrow(response)
if(nrow(censor)!=nind)stop("response and censor must have the same number of rows")
if(!is.null(nest)&&length(nest)!=nind)stop(paste("nest must have length",nind1))
if(!is.null(cov)&&length(cov)!=nind)stop(paste("cov must have length",nind1))
tot <- 0
interval <- nnest <- counts <- event <- ncov <- resp <- NULL
for(k in 1:ncol(response)){
	o <- order(response[,k])
	cc <- censor[o,k]
	i <- nest[o]
	tt <- sort(response[,k])
	tt <- c(cc[1],(cc[2:nind]*(tt[2:nind]>tt[1:(nind-1)])))
	nint <- 1
	id <- i[1]
	nevent <- cc[1]
	nz <- 0
	if(!missing(cov)){
		tcov <- cov[o]
		ncov1 <- tcov[1]}
	if(!stratified){
		tresp <- resp1[o]
		nresp <- tresp[1]}
	for(j in 2:nind){
		nz <- nz+(tt[j]==0)
		ev <- rep(0,j-nz)
		ev[j-nz] <- cc[j]
		nevent <- c(nevent,ev)
		id <- c(id,rep(i[j],j-nz))
		nint <- c(nint,1:(j-nz))
		if(!missing(cov))ncov1 <- c(ncov1,rep(tcov[j],j-nz))
		if(!stratified)nresp <- c(nresp,rep(tresp[j],j-nz))}
	interval <- c(interval,nint+tot)
	nnest <- c(nnest,id)
	counts <- c(counts,nevent)
	event <- c(event,rep(k,length(nint)))
	if(!missing(cov))ncov <- c(ncov,ncov1)
	if(!stratified)resp <- c(resp,nresp)
	tot <- max(interval)}
if(is.vector(response,mode="numeric"))stratified <- T
if(estimate<=0)stop("estimate must be positive")
fcov <- !missing(cov)
if(!fcov){
	if(stratified)lmod <- counts~1
	else lmod <- counts~resp}
else {
	if(stratified)lmod <- counts~ncov
	else lmod <- counts~ncov+resp}
rand <- rep(1,max(nnest))
fixd <- rep(1,max(interval))
fobs <- collapse(counts,nnest)
ffob <- collapse(counts,interval)
z <- glm(lmod,family=poisson)
fit <- z$fit
rm(z)
z1 <- nlm(likel,log(estimate), hessian=T, print.level=print.level,
	typsiz=typsiz, ndigit=ndigit, gradtol=gradtol, stepmax=stepmax,
	steptol=steptol, iterlim=iterlim, fscale=fscale)
z2 <- like(z1$est)
z3 <- list(maxlike=z2$like,
	aic=z2$like+length(z2$coef)+length(z2$fixed)+1,
	df=length(response)-length(z2$coef)-length(z2$fixed)-1,
	iterations=z1$iter,
	code=z1$code,
	call=call,
	gamma=exp(z1$est),
	fixed=z2$fixed,
	random=z2$random,
	inthaz=z2$inthaz,
	coefficients=z2$coef,
	stratified=stratified)
class(z3) <- "llrf"
z3}

print.llrf <- function(z){
	if(z$stratified)cat("Stratified ")
	cat("Cox proportional hazards model with gamma frailty\n")
	cat("\nCall:",deparse(z$call),sep="\n")
	cat("\n")
	if(z$code>2)cat("Warning: no convergence - error",z$code,"\n\n")
	cat("-Log likelihood   ",z$maxlike,"\n")
	cat("Degrees of freedom",z$df,"\n")
	cat("AIC               ",z$aic,"\n")
	cat("Iterations        ",z$iterations,"\n\n")
	tmp <- trigamma(1/z$gamma)
	cat("gamma =      ",z$gamma,"\n")
	cat("correlation =",tmp/(tmp+trigamma(1)),"\n")
	cat("\nRegression coefficients:\n")
	print(z$coef)
	cat("\nFixed effects:\n")
	print(z$fixed)
	cat("\nRandom effects:\n")
	print(z$random)}
