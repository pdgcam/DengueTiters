/*
 * DS.h
 *
 *  Created on: Jul 20, 2014
 *      Author: hsalje
 */



//#include <omp.h>
#define DO_OMP_PARALLEL 0
#ifdef DO_OMP_PARALLEL
#define NUM_THREADS 1 //function LogLik no longer parallelized
#define CACHE_LINE_SIZE 128
#else
#define NUM_THREADS 1
#define CACHE_LINE_SIZE 1
#endif

//#define isnan(x) ((x) != (x))

//double ran1(long *);
//double runif();
//double rnorm();
//double rexp(double lambda);
//double dexp(double x,double lambda);
//double ldexp(double x,double lambda);
//int rbin(int n, double p);
//int poidev(double x);
//double ldbin(int n, int k, double p);
//double ldmultinomial(int n, double *proba,int size,int* deviate);
//double gammln(double x);
//double gammq(double a, double x);
//void rmultinomial(int n, double *proba,int size,int* deviate);
//int rmultinomialOneValue(double *proba,int size);
//int rnegbin(double n,double p);
//double digammal(double x);
//double trigamma(double x);
//double pgamma(double x,double a, double b);
//double dgamma(double x,double a,double b);
//double rgamma(double a,double b);
//double pnorm(double x,double mean,double var);
//double plnorm(double x,double mean,double var);
//double dlnorm(double x,double mean,double var);
//double logdlnorm(double x,double mean,double var);
//double pweibull(double x,double lambda,double k);
//double dweibull(double x,double lambda,double k);
//double pchisqNCnull(double X, double nonCentralPar);
//double pnorm(double x);
//double logPLnormDiscrete(double x,double mean,double sd);
//double bessi1(double x);
//double dlchisqNCnull(double x,double theta);
//double dlnegbin(int x,double n,double p);
//double dlpois(int x,double r);
//double pnegbin(double k, double a, double b);
//double dnegbin(int x,double n,double p);
//double erff(double x);
//double rTruncNormalRight(double mu,double var,double muMax);
//double choose(int n,int k);
//double pbeta(double x,double a,double b);
//double dbeta(double x,double a, double b);
////============
//
//long SEED; // root random number
//


#ifndef DS_H_
#define DS_H_

//RANDOM NUMBERS
extern double ran1(long *);
extern double runif();
extern double rnorm();
extern double rexp(double lambda);
extern double dexp(double x,double lambda);
extern double pexp(double x,double lambda);
extern double ldexp(double x,double lambda);
extern int rbin(int n, double p);
extern int poidev(double x);
extern double ldbin(int n, int k, double p);
extern double ldmultinomial(int n, double *proba,int size,int* deviate);
extern double gammln(double x);
extern double gammq(double a, double x);
extern void rmultinomial(int n, double *proba,int size,int* deviate);
extern int rmultinomialOneValue(double *proba,int size);
extern int rnegbin(double n,double p);
extern double digammal(double x);
extern double trigamma(double x);
extern double pgamma(double x,double a, double b);
extern double dgamma(double x,double a,double b);
extern double rgamma(double a,double b);
extern double  dnorm(double x, double m, double s);
extern double dnorm(double x,double mean,double var);
extern double plnorm(double x,double mean,double var);
extern double dlnorm(double x,double mean,double var);
extern double logdlnorm(double x,double mean,double var);
extern double pweibull(double x,double lambda,double k);
extern double dweibull(double x,double lambda,double k);
extern double pchisqNCnull(double X, double nonCentralPar);
extern double pnorm(double x);
extern double logPLnormDiscrete(double x,double mean,double sd);
extern double bessi1(double x);
extern double dlchisqNCnull(double x,double theta);
extern double dlnegbin(int x,double n,double p);
extern double dlpois(int x,double r);
extern double pnegbin(double k, double a, double b);
extern double dnegbin(int x,double n,double p);
extern double erff(double x);
extern double rTruncNormalRight(double mu,double var,double muMax);
extern double choose(int n,int k);
extern double pbeta(double x,double a,double b);
extern double dbeta(double x,double a, double b);
//============
#endif /* DS_H_ */


//#ifdef SEED_H_
//#define SEED_H_
//; // root random number
//#endif
