/*
 *  event : A Library of Special Functions for Event Histories
 *  Copyright (C) 1998 J.K. Lindsey
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  SYNOPSIS
 *
 * void ksurvb(double p[],double y[],double x[],int cens[],int *nind,
 *	    int nobs[],int *nbs,int *nccov,int *model,int *density,int *dep,
 *	    int *birth,int *tvc,double tvcov[],int *fit,double pred[],
 *	    double rpred[],int *renewal,int *rf,double bb[],int *sf,
 *	    double vv[],double *like)
 * void frailb(double p[],double y[],double x[],int cens[],int *nind,
 *	    int nobs[],int *nbs,int *nccov,int *model,int *density,int *dep,
 *	    int *birth,int *tvc,double tvcov[],int *fit,double pred[],
 *	    double rpred[],int *rf,double bb[],int *sf,double vv[],
 *	    double *like)
 *
 *  DESCRIPTION
 *
 *    Function to compute the likelihood function for various distributions
 * inserted in a beta distribution with serial dependence or gamma frailties
 * using Kalman-type update for event histories.
 *
 */

#include <math.h>
#include <stddef.h>

extern double lgamma(double x);
extern double dexp(double x, double scale);
extern double pexp(double x, double scale);
extern double qexp(double x, double scale);
extern double dweibull(double x, double shape, double scale);
extern double pweibull(double x, double shape, double scale);
extern double qweibull(double x, double shape, double scale);
extern double dgamma(double x, double shape, double scale);
extern double pgamma(double x, double shape, double scale);
extern double qgamma(double p, double shape, double scale);
extern double dnorm(double x, double mean, double sd);
extern double pnorm(double x, double mean, double sd);
extern double qnorm(double p, double mu, double sigma);
extern double dlogis(double x, double location, double scale);
extern double plogis(double x, double location, double scale);
extern double qlogis(double x, double location, double scale);
extern double dcauchy(double x, double location, double scale);
extern double pcauchy(double x, double location, double scale);
extern double qcauchy(double x, double location, double scale);
extern double ihgamma(double x, double shape, double scale);
extern double ihlogis(double x, double location, double scale);

void ksurvb(double p[],double y[],double x[],int cens[],int *nind,int nobs[],
	    int *nbs,int *nccov,int *model,int *density,int *dep,int *birth,
	    int *tvc,double tvcov[],int *fit,double pred[],double rpred[],
	    int *renewal,int *rf,double bb[],int *sf,
	    double vv[],double *like){
  int i,j,k,nm,c;
  double a,a1,b,b1,delta,lambda,omega,om,beta,bt,h,yy,kk,tmp,ly,plap,intercept;
  
  *like=0;
  nm=0;
  delta=exp(-p[*nccov+*birth+*tvc+1]);
  if(*dep>0)omega=exp(p[*nccov+*birth+*tvc+2])/(1+exp(p[*nccov+*birth+*tvc+2]));
  if(*model>1&&!*sf){
    if(*model<5)lambda=exp(p[*nccov+*birth+*tvc+2+(*dep>0)]);
    else lambda=exp(p[*nccov+*birth+*tvc+2+(*dep>0)]/2);}
  if(*model==4)intercept=exp(p[*nccov+*birth+*tvc+3+(*dep>0)]);
  for(i=0;i<*nind;i++){
    a=b=delta;
    if(!*rf){
      beta=p[0];
      for(k=0;k<*nccov;k++)beta+=p[k+1]*x[i+k**nind];
      if(*model<4){
	if(beta>40) beta=40;
	if(beta<-40)beta=-40;
	beta=exp(beta);}}
    else if(!*tvc)bt=bb[i];
    c=0;
    yy=0;
    for(j=0;j<nobs[i];j++){
      if(*model>1&&*sf)lambda=vv[nm];
      /* store value and check if ties to follow */
      if(y[nm]>0){
	if(*renewal)yy=y[nm];
	else yy+=y[nm];}
      c+=cens[nm];
      if(j>=nobs[i]-1||y[nm+1]!=0){
	/* if no ties follow, update the likelihood */
	if(*model>=5)ly=log(yy);
	a1=a+c;
	b1=b;
	/* add in birth and time-varying covariates */
	if(!*rf){
	  if(*tvc){
	    bt=0;
	    for(k=0;k<*tvc;k++)bt+=p[*nccov+*birth+k+1]*tvcov[nm+*nbs*k];
	    if(*model<4){
	      if(bt>40) bt=40;
	      if(bt<-40)bt=-40;
	      bt=exp(bt)*beta;}
	    else bt+=beta;}
	  else bt=beta;}
	  if(*birth){
	    if(*model<4)bt*=pow(j+1.,p[*nccov+1]);
	    else bt+=p[*nccov+1]*log(j+1);}
	else if(*tvc)bt=bb[nm];
	if(!*density){
	  /* intensity models */
	  switch(*model){
	  case 1:
	    b1+=yy/bt;
	    h=1/bt;
	    break;
	  case 2:
	    b1+=pow(yy/bt,lambda);
	    h=lambda*pow(yy/bt,lambda-1)/bt;
	    break;
	  case 3:
	    b1+=ihgamma(yy,lambda,bt);
	    h=dgamma(yy,lambda,bt)/(1-pgamma(yy,lambda,bt));
	    break;
	  case 4:
	    b1+=(yy+log(lambda+intercept*exp(-bt*yy))/bt)/lambda;
	    h=1/(lambda+intercept*exp(-bt*yy));
	    break;
	  case 5:
	    b1-=log(1-pnorm(ly,bt,lambda));
	    h=dnorm(ly,bt,lambda)/yy/(1-pnorm(ly,bt,lambda));
	    break;
	  case 6:
	    b1+=ihlogis(ly,bt,lambda);
	    h=dlogis(ly,bt,lambda)/yy/(1-plogis(ly,bt,lambda));
	    break;
	  case 7:
	    b1-=log(1-pcauchy(ly,bt,lambda));
	    h=dcauchy(ly,bt,lambda)/yy/(1-pcauchy(ly,bt,lambda));
	    break;
	  case 8:
	    tmp=exp(-fabs(ly-bt)/lambda)/2;
	    plap=ly<bt?tmp:1-tmp;
	    b1-=log(1-plap);
	    h=tmp/(lambda*yy*(1-plap));
	    break;}}
	else{
	  /* density models */
	  switch(*model){
	  case 1:
	    b1+=pexp(yy,bt);
	    h=dexp(yy,bt);
	    break;
	  case 2:
	    b1+=pweibull(yy,lambda,bt);
	    h=dweibull(yy,lambda,bt);
	    break;
	  case 3:
	    b1+=pgamma(yy,lambda,bt);
	    h=dgamma(yy,lambda,bt);
	    break;
	  case 4:
	    b1+=exp(-yy/lambda)*pow((lambda+intercept)/(lambda+intercept*exp(-bt*yy)),1/(lambda*bt));
	    h=exp(-yy/lambda)*pow((lambda+intercept)/(lambda+intercept*exp(-bt*yy)),1/(lambda*bt)+1);
	    break;
	  case 5:
	    b1+=pnorm(ly,bt,lambda);
	    h=dnorm(ly,bt,lambda)/yy;
	    break;
	  case 6:
	    b1+=plogis(ly,bt,lambda);
	    h=dlogis(ly,bt,lambda)/yy;
	    break;
	  case 7:
	    b1+=pcauchy(ly,bt,lambda);
	    h=dcauchy(ly,bt,lambda)/yy;
	    break;
	  case 8:
	    tmp=exp(-fabs(ly-bt)/lambda)/2;
	    b1+=ly<bt?tmp:1-tmp;
	    h=tmp/lambda/yy;
	    break;}}
	/* calculate likelihood */
	*like-=(c>1?lgamma(a1)-lgamma(a)-lgamma(c+1):c*log(a))+a*log(b)-a1*log(b1)+c*log(h);
	if(*fit){
	  pred[nm-c+cens[nm]]=bt;
	  tmp=b/a;
	  if(!*density){
	    switch(*model){
	    case 1: rpred[nm-c+cens[nm]]=bt*tmp; break;
	    case 2: rpred[nm-c+cens[nm]]=bt*pow(tmp,1/lambda); break;
	    case 3: rpred[nm-c+cens[nm]]=qgamma(1-exp(-tmp),lambda,bt); break;
	    case 5: rpred[nm-c+cens[nm]]=exp(qnorm(1-exp(-tmp),bt,lambda)); break;
	    case 6: rpred[nm-c+cens[nm]]=exp(qlogis(1-exp(-tmp),bt,lambda)); break;
	    case 7: rpred[nm-c+cens[nm]]=exp(qcauchy(1-exp(-tmp),bt,lambda)); break;
	    case 8: rpred[nm-c+cens[nm]]=exp(bt+lambda*log(2*(ly<bt?exp(-tmp):1-exp(-tmp)))); break;}}
	  else{
	    switch(*model){
	    case 1: rpred[nm-c+cens[nm]]=qexp(tmp,bt); break;
	    case 2: rpred[nm-c+cens[nm]]=qweibull(tmp,lambda,bt); break;
	    case 3: rpred[nm-c+cens[nm]]=qgamma(tmp,lambda,bt); break;
	    case 5: rpred[nm-c+cens[nm]]=exp(qnorm(tmp,bt,lambda)); break;
	    case 6: rpred[nm-c+cens[nm]]=exp(qlogis(tmp,bt,lambda)); break;
	    case 7: rpred[nm-c+cens[nm]]=exp(qcauchy(tmp,bt,lambda)); break;
	    case 8: rpred[nm-c+cens[nm]]=exp(bt+lambda*log(2*(ly<bt?tmp:1-tmp))); break;}}}
	      /* update parameters */
        switch(*dep){
	case 1:
	case 2: om=pow(omega,yy); a=om*a1+(1-om)*delta; break;
	case 3:
	case 4: a=omega*a1+(1-omega)*delta; break;
	case 5: a=a1; break;
	case 6:
	case 7: a=omega*a1; break;
	default:}
	switch(*dep){
	case 1: b=om*(b1-b)+delta; break;
	case 2: b=om*b1+(1-om)*delta; break;
	case 3: b=omega*(b1-b)+delta; break;
	case 4: b=omega*b1+(1-omega)*delta; break;
	case 5:
	case 7: b=omega*b1; break;
	default:}
      c=0;}
      nm++;}}
  return;}

void frailb(double p[],double y[],double x[],int cens[],int *nind,int nobs[],
	    int *nbs,int *nccov,int *model,int *density,int *dep,int *birth,
	    int *tvc,double tvcov[],int *fit,double pred[],double rpred[],
	    int *rf,double bb[],int *sf,double vv[],
	    double *like){
  int i,j,k,nm,ns,c,nn;
  double b1,delta,lambda,beta,bt,l1,yy,kk,nb,ly,plap,tmp,intercept,H;
  
  *like=0;
  nm=0;
  delta=exp(p[*nccov+*birth+*tvc+1]);
  if(*model>1&&!*sf){
    if(*model<5)lambda=exp(p[*nccov+*birth+*tvc+2]);
    else lambda=exp(p[*nccov+*birth+*tvc+2]/2);}
  if(*model==4)intercept=exp(p[*nccov+*birth+*tvc+3]);
  for(nn=i=0;i<*nind;i++)nn+=nobs[i];
  for(i=0;i<*nind;i++){
    if(!*rf){
      beta=p[0];
      for(k=0;k<*nccov;k++)beta+=p[k+1]*x[i+k**nind];
      if(*model<4){
	if(beta>40) beta=40;
	if(beta<-40)beta=-40;
	beta=exp(beta);}}
    else if(!*tvc)bt=bb[i];
    ns=b1=0;
    if(!*dep)nb=1;
    else {
      nb=0;
      for(j=0;j<nobs[i];j++)nb+=y[nm+j];
      nb/=beta;}
    for(c=0,j=0;j<nobs[i];j++){
      if(*model>1&&*sf)lambda=vv[nm];
      ns+=cens[nm];
      /* store value and check if ties to follow */
      if(y[nm]>0)yy=y[nm];
      c+=cens[nm];
      if(j>=nobs[i]-1||y[nm+1]!=0){
	if(*model>=5)ly=log(yy);
	l1=log(1+delta*j/nb);
	/* add in birth and time-varying covariates */
	if(!*rf){
	  if(*tvc){
	    bt=0;
	    for(k=0;k<*tvc;k++)bt+=p[*nccov+*birth+k+1]*tvcov[nm+*nbs*k];
	    if(*model<4)bt=exp(bt)*beta;
	    else bt+=beta;}
	  else bt=beta;
	  if(*birth){
	    if(*model<4)bt*=pow(j+1.,p[*nccov+1]);
	    else bt+=p[*nccov+1]*log(j+1);}}
	else if(*tvc)bt=bb[nm];
	if(!*density){
	  /* intensity models */
	  switch(*model){
	  case 1:
	    H=yy/bt;
	    l1+=-log(bt);
	    break;
	  case 2:
	    H=pow(yy/bt,lambda);
	    l1+=log(lambda/bt)+(lambda-1)*log(yy/bt);
	    break;
	  case 3:
	    H=ihgamma(yy,lambda,bt);
	    /*b1-=log(1-pgamma(yy,lambda,bt));*/
	    l1+=log(dgamma(yy,lambda,bt)/(1-pgamma(yy,lambda,bt)));
	    break;
	  case 4:
	    H=(yy+log(lambda+intercept*exp(-bt*yy))/bt)/lambda;
	    l1+=-log(lambda+intercept*exp(-bt*yy));
	    break;
	  case 5:
	    H=-log(1-pnorm(ly,bt,lambda));
	    l1+=log(dnorm(ly,bt,lambda)/yy/(1-pnorm(ly,bt,lambda)));
	    break;
	  case 6:
	    H=ihlogis(ly,bt,lambda);
	    l1+=log(dlogis(ly,bt,lambda)/yy/(1-plogis(ly,bt,lambda)));
	    break;
	  case 7:
	    H=-log(1-pcauchy(ly,bt,lambda));
	    l1+=log(dcauchy(ly,bt,lambda)/yy/(1-pcauchy(ly,bt,lambda)));
	    break;
	  case 8:
	    tmp=exp(-fabs(ly-bt)/lambda)/2;
	    plap=ly<bt?tmp:1-tmp;
	    H=-log(1-plap);
	    l1+=tmp/(lambda*yy*(1-plap));
	    break;}}
	else{
	  /* density models */
	  switch(*model){
	  case 1:
	    H=pexp(yy,bt);
	    l1+=log(dexp(yy,bt));
	    break;
	  case 2:
	    H=pweibull(yy,lambda,bt);
	    l1+=log(dweibull(yy,lambda,bt));
	    break;
	  case 3:
	    H=pgamma(yy,lambda,bt);
	    l1+=log(dgamma(yy,lambda,bt));
	    break;
	  case 4:
	    H=exp(-yy/lambda)*pow((lambda+intercept)/(lambda+intercept*exp(-bt*yy)),1/(lambda*bt));
	    l1+=-yy/lambda+log((lambda+intercept)/(lambda+intercept*exp(-bt*yy)))/((lambda*bt)+1);
	    break;
	  case 5:
	    H=pnorm(ly,bt,lambda);
	    l1+=log(dnorm(ly,bt,lambda)/yy);
	    break;
	  case 6:
	    H=plogis(ly,bt,lambda);
	    l1+=log(dlogis(ly,bt,lambda)/yy);
	    break;
	  case 7:
	    H=pcauchy(ly,bt,lambda);
	    l1+=log(dcauchy(ly,bt,lambda)/yy);
	    break;
	  case 8:
	    tmp=exp(-fabs(ly-bt)/lambda)/2;
	    H=ly<bt?tmp:1-tmp;
	    l1+=log(tmp/lambda/yy);
	    break;}}
	/* calculate likelihood */
	*like-=c*l1;
	if(c>1)*like+=lgamma(c+1);
	if(*fit){
	  pred[nm-c+cens[nm]]=bt;
	  tmp=(b1+nb/(nn*delta))/(nb/(nn*delta)+ns);
	  if(!*density){
	    switch(*model){
	    case 1: rpred[nm-c+cens[nm]]=bt*tmp; break;
	    case 2: rpred[nm-c+cens[nm]]=bt*pow(tmp,1/lambda); break;
	    case 3: rpred[nm]=qgamma(1-exp(-tmp),lambda,bt); break;
	    case 5: rpred[nm-c+cens[nm]]=exp(qnorm(1-exp(-tmp),bt,lambda)); break;
	    case 6: rpred[nm-c+cens[nm]]=exp(qlogis(1-exp(-tmp),bt,lambda)); break;
	    case 7: rpred[nm-c+cens[nm]]=exp(qcauchy(1-exp(-tmp),bt,lambda)); break;
	    case 8: rpred[nm-c+cens[nm]]=exp(bt+lambda*log(2*(ly<bt?exp(-tmp):1-exp(-tmp)))); break;}}
	  else{
	    switch(*model){
	    case 1: rpred[nm-c+cens[nm]]=qexp(tmp,bt); break;
	    case 2: rpred[nm-c+cens[nm]]=qweibull(tmp,lambda,bt); break;
	    case 3: rpred[nm-c+cens[nm]]=qgamma(tmp,lambda,bt); break;
	    case 5: rpred[nm-c+cens[nm]]=exp(qnorm(tmp,bt,lambda)); break;
	    case 6: rpred[nm-c+cens[nm]]=exp(qlogis(tmp,bt,lambda)); break;
	    case 7: rpred[nm-c+cens[nm]]=exp(qcauchy(tmp,bt,lambda)); break;
	    case 8: rpred[nm-c+cens[nm]]=exp(bt+lambda*log(2*(ly<bt?tmp:1-tmp))); break;}}}
	b1+=H;
	c=0;}
      nm++;}
    *like+=(nb/delta+ns)*log(1+delta*b1/nb);}
  return;}
