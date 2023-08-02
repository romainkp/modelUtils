
#ifndef __HAS_MODEL_FITTERS__
#define __HAS_MODEL_FITTERS__

#include <math.h>

#include <R.h>
#include <R_ext/Utils.h>
#include <R_ext/Applic.h>
#include <Rmath.h>

#define MODEL_NO_ERROR                  0
#define MODEL_ERROR_AT_BOUNDARY         1
#define MODEL_ERROR_RANK_CHANGE         2
#define MODEL_ERROR_NO_CONVERGENCE      3
#define MODEL_ERROR_UNKNOWN_MODEL_TYPE  4
#define MODEL_ERROR_UNKNOWN             10

/** global constants specifying the debug level. **/
int R_MODEL_DEBUG_LEVEL;

#define R_MODEL_DEBUG_TRACE   1
#define R_MODEL_DEBUG_LOW     2 
#define R_MODEL_DEBUG_WARN    3 

#define MODEL_LINEAR 1
#define MODEL_BINOMIAL 2
#define MODEL_MULTINOMIAL 3

/** starting with R 4.3.0 , DBL_EPSILON is no loger defined in R.h file **/
#ifndef DOUBLE_EPS
# define DOUBLE_EPS 2.2204460492503131e-16
#endif

/** snagged from R. **/
static const double THRESH = 30;
static const double MTHRESH = -30;
static const double INVEPS = 1/DOUBLE_EPS;
static const double MY_ZERO = (DOUBLE_EPS*100);


static R_INLINE double x_d_opx(double x) {
  return x/(1 + x);
}


void __do_copy(double* x,double* y,double *weights, int n, int p, double* xt, int M);


int* model_allocate_int(int model_type, int N, int P, int M);
double* model_allocate_double(int model_type, int N, int P, int M);

void fit_model_with_allocation(int* model_type, double* X, double* Y, int *N, int *P, int *M, 
			       double* weights, double* estimates, double* loss, double* residuals, 
			       int* error, double* tolerance, int* max_iter);

void fit_model(int* model_type, double* X, double* Y, int* N, int* P, int* M, 
	       double* weights, double* estimates, double* loss, double* residuals, 
	       int* error, double* tolerance, int* max_iter, int* work_i, 
	       double* work_d);

void fit_multinomial(double* X, double* Y, int N, int P, int M, double* weights, 
		     double* estimates, double* loss, double* residuals, int* error, 
		     double tolerance, int max_iter, int* work_i, 
		     double* work_d);

void fit_binomial(double* X, double* Y, int N, int P, int M, double* weights, 
		  double* estimates, double* loss, double* residuals, int* error, 
		  double tolerance, int max_iter, int* work_i, 
		  double* work_d);

void fit_linear(double* X, double* Y, int N, int P, int M, double* weights, 
		double* estimates, double* loss, double* residuals, int* error, 
		double tolerance, int max_iter, int* work_i, 
		double* work_d); 


void __do_multinomial(double* X, int N, int P, int M, double* Y, 
		      double tol, double* coefficients, 
		      double* residuals, int* rank, int* pivot,
		      double* weights, double* mu, double* XTWX, double* I,
		      double* grad_theta, int max_iter, double* new_coefficients, 
		      double* user_weights, int* error, double* loss);

void __do_logistic(double* X, int n, int p, double* Y, 
		   double tol, double* coefficients, 
		   double* residuals, double* effects, int* rank, 
		   int* pivot, double* qraux, double* work, double* dev,
		   double* weights, int max_iter, double* new_coefficients, 
		   double* WX, double* y_shifted, int* error);

void __do_linear(double* X, int n, int p, double* Y, 
		 double* weights, double tol, double* coefficients, 
		 double* residuals, double* effects, int* rank, 
		 int* pivot, double* qraux, double* work, double* rss,
		 int* fail);


#endif /** __HAS_MODEL_FITTERS__ **/
