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
#	hboxcox(y,m,s,f)
#	hburr(y,m,s,f)
#	hcauchy(y,m,s)
#	hexp(y,m)
#	hgextval(y,s,m,f)
#	hgamma(y,s,m)
#	hggamma(y,s,m,f)
#	hhjorth(y,m,s,f)
#	hinvgauss(y,m,s)
#	hlaplace(y,m,s)
#	hlnorm(y,m,s)
#	hlogis(y,m,s)
#	hglogis(y,m,s,f)
#	hnorm(y,m,s)
#	hpareto(y,m,s)
#	hstudent(y,m,s,f)
#	hweibull(y,s,m)
#	hgweibull(y,s,m,f)
#
#  DESCRIPTION
#
#    Functions for various log hazards or intensities

# for log distributions, subtract log(y) from the intensity function

# f=1 gives truncated normal
hboxcox <- function(y,m,s,f) {
	y1 <- y^f/f
	-(y1-m)^2/s/2+(f-1)*log(y)-log(2*pi*s)/2-log(1-pnorm(y1,m,sqrt(s))+(f<0)*(1-pnorm(0,m,sqrt(s))))}

hburr <- function(y,m,s,f) {
	y1 <- y/m
	y2 <- y1^s
	log(f*s/m)+(s-1)*log(y1)-log(1+y2)}

hcauchy <- function(y,m,s) log(dcauchy(y,m,s))-log(1-pcauchy(y,m,s))

hexp <- function(y,m) -log(m)*rep(1,length(y))

# f=1 gives truncated extreme value
hgextval <- function(y,s,m,f) {
	y1 <- y^f/f
	ey <-exp(y1)
	log(s)+s*(y1-log(m))-(ey/m)^s+(f-1)*log(y)-log(1-pweibull(ey,s,m)-(f<0)*exp(-m^-s))}

hgamma <- function(y,s,m) log(dgamma(y,s,m))-log(1-pgamma(y,s,m))

hggamma <- function(y,s,m,f) {
	t <- m/s
	u <- t^f
	y1 <- y^f
	v <- s*f
	-v*log(t)-y1/u+log(f)+(v-1)*log(y)-lgamma(s)-log(1-pgamma(y1,s,u))}

hhjorth <- function(y,m,s,f) log(y/m^2+f/(1+s*y))

hinvgauss <- function(y,m,s) {
	t <- y/m
	v <- sqrt(y*s)
	-((t-1)^2/(y*s)+log(2*s*pi*y^3))/2-log(1-pnorm((t-1)/v)
		-exp(2/(m*s))*pnorm(-(t+1)/v))}

hlaplace <- function(y,m,s){
	plp <- function(u){
		t <- exp(-abs(u))/2
		ifelse(u<0,t,1-t)}
	-abs(y-m)/s-log(2*s)-log(1-plp((y-m)/s))}

hlnorm <- function(y,m,s) log(dlnorm(y,m,s))-log(1-plnorm(y,m,s))

hlogis <- function(y,m,s) log(dlogis(y,m,s))-log(1-plogis(y,m,s))

# f=1 gives hlogis
hglogis <- function(y,m,s,f) {
	y1 <- (y-m)/s
	ey <- exp(-y1)
	-log(s/f)-y1-(f+1)*log(1+ey)-log(1-(1+ey)^-f)}

hnorm <- function(y,m,s) log(dnorm(y,m,s))-log(1-pnorm(y,m,s))

hpareto <- function(y,m,s) (s+1)/(m*s+y)

hstudent <- function(y,m,s,f){
	pst <- function(u,f){
		t <- 0.5*pbeta(f/(f+u^2),f/2,0.5)
		ifelse(u<0,t,1-t)}
	t <- (f+1)/2
	u <- (y-m)/s
	lgamma(t)-lgamma(f/2)-log(f)/2-(t)*log(1+u^2/f)
		-log(pi)/2-log(1-pst(u,f))}

hweibull <- function(y,s,m) log(s)+(s-1)*log(y)-s*log(m)

# Mudholkar, Srivastava, & Freimer (1995) Technometrics 37: 436-445
hgweibull <- function(y,s,m,f) {
	y1 <- y/m
	y2 <- y1^s
	y3 <- exp(-y2)
	log(s*f/m)+(s-1)*log(y1)+(f-1)*log(1-y3)-y2-log(1-(1-y3)^f)}
