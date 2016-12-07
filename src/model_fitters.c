
#include <math.h>

#include <Rinternals.h>
#include <R.h>
#include <R_ext/Utils.h>
#include <R_ext/Applic.h>
#include <Rmath.h>

#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#include "model_fitters.h"



void get_epsilon_compare(double* eps) {
  *eps = MY_ZERO;
}

int* model_allocate_int(int model_type, int n, int p, int m)
{
  return (int*) R_alloc(m*p, sizeof(int));
}

double* model_allocate_double(int model_type, int n, int p, int m)
{
  double* work_d;
  int dcopy = n*(p+m+1);

  if (model_type == MODEL_LINEAR) {
    work_d = (double*) R_alloc(n + /* effects */
			      p + /* qraux */
			      2*p + /* work */
			      dcopy, /* X/Y copy and weights*/
			      sizeof(double));
  }
  else if (model_type == MODEL_BINOMIAL) {
    work_d = (double*) R_alloc(n + /* effects */
			      p + /* qraux */
			      2*p + /* work */
			      n*p + /* WX */
			      n + /* Y shifted */
			      p + /* new coefficients */
			      dcopy, /* X/Y copy  and weights*/
			      sizeof(double));
  }
  else if (model_type == MODEL_MULTINOMIAL) {
    work_d = (double*) R_alloc(n*m*m + /* W */
			      n*m + /* mu */
			      m*p*m*p + /* Hessian */
			      m*p + /* Gradient */
			      m*p + /* next estimates */ 
			      m*p*m*p + /* I */
			      dcopy, /* X/Y copy and weights*/
			      sizeof(double));
  }
  else {
    Rprintf("error: unknown model type!\n");
  }
  return work_d;
}

void fit_model_with_allocation(int* model_type, double* X, double* Y, int* N, int* P, int* M, 
			       double* weights, double* estimates, double* loss, 
			       double* residuals, int* error, double* tolerance, 
			       int* max_iter)
{
  int* work_i = model_allocate_int(*model_type, *N, *P, *M);
  double* work_d = model_allocate_double(*model_type, *N, *P, *M);
  fit_model(model_type, X, Y, N, P, M, weights, estimates, loss, 
	    residuals, error, tolerance, max_iter, work_i, work_d);

}

void fit_model(int* model_type, double* X, double* Y, int* N, int* P, int* M, 
	       double* weights, double* estimates, double* loss, double* residuals, 
	       int* error, double* tolerance, int* max_iter, int* work_i, 
	       double* work_d)
{
  int n, m, p;
  n = *N; m = *M; p = *P;
  
  if (*model_type == MODEL_LINEAR) {
    fit_linear(X, Y, n, p, m, weights, estimates, loss, residuals, error,
	       *tolerance, *max_iter, work_i, work_d);
  }
  else if (*model_type == MODEL_BINOMIAL) {
    fit_binomial(X, Y, n, p, m, weights, estimates, loss, residuals, error,
		 *tolerance, *max_iter, work_i, work_d);
  }
  else if (*model_type == MODEL_MULTINOMIAL) {
    fit_multinomial(X, Y, n, p, m, weights, estimates, loss, residuals, error,
		    *tolerance, *max_iter, work_i, work_d); 
  }
  else {
    Rprintf("Error, unknown model type!\n");
    *error = MODEL_ERROR_UNKNOWN_MODEL_TYPE;
  }
}

void fit_multinomial(double* X, double* Y, int N, int P, int M, double* weights, 
		     double* estimates, double* loss, double* residuals, int* error, 
		     double tolerance, int max_iter, int* work_i, 
		     double* work_d) 
{
  /** unchunk memory. **/
  double* W = work_d;
  double* mu = (W + (N*M*M));
  double* H = (mu + (N*M));
  double* G = (H + (M*P)*(M*P));
  double* next_estimates = (G + (M*P));
  double* I = (next_estimates + (M*P));
  double* t_space = (I + (M*P)*(M*P));
  double *tmp,*tmp2;

  int* pivots = work_i;
  int* rank; 
  int j = 0, k = 0; 

  for (j = 0; j < (M*P); j++) {
    for (k = 0; k < (M*P); k++) {
      if (k == j) {
	I[j*(M*P) + k] = 1.0;
      }
      else {
	I[j*(M*P) + k] = 0.0;
      }
    }
    pivots[j] = (j+1);
  }

  __do_copy(X,Y,weights, P, N, t_space,M);
  tmp = X;
  X = t_space;
  t_space = tmp;

  tmp = Y;
  Y = X+N*P;
  tmp2 = weights;
  weights = X+N*(P+M);

  __do_multinomial(X, N, P, M, Y, tolerance, estimates, residuals, rank,
		   pivots, W, mu, H, I, G, max_iter, next_estimates, weights, 
		   error, loss);

  X = t_space;
  Y = tmp;
  weights = tmp2;
}

void fit_binomial(double* X, double* Y, int N, int P, int M, double* weights, 
		  double* estimates, double* loss, double* residuals, int* error, 
		  double tolerance, int max_iter, int* work_i, 
		  double* work_d) 
{
  double* effects = work_d;
  double* qraux = (work_d + N);
  double* work = (qraux + P);
  double* WX = (work + 2*P);
  double* Y_shifted = (WX + N*P);
  double* next_estimates = (Y_shifted + N);
  double* t_space = (next_estimates + P); 
  double *tmp,*tmp2;

  int* pivots = work_i;
  int rank = P;
  int j;

  for (j = 0; j < P; j++) { 
    pivots[j] = j + 1;
  }
    
  __do_copy(X,Y, weights,P, N, t_space,M);
  tmp = X;
  X = t_space;
  t_space = tmp;

  tmp = Y;
  Y = X+N*P;
  tmp2 = weights;
  weights = X+N*(P+M);

  __do_logistic(X, N, P, Y, tolerance, estimates, residuals, effects,
		&rank, pivots, qraux, work, loss, weights, max_iter, 
		next_estimates, WX, Y_shifted, error);
    
  X = t_space;
  Y = tmp;
  weights = tmp2;
}

void fit_linear(double* X, double* Y, int N, int P, int M, double* weights, 
		double* estimates, double* loss, double* residuals, int* error, 
		double tolerance, int max_iter, int* work_i, 
		double* work_d)
{
  double* effects = work_d;
  double* qraux = (effects + N);
  double* work = (qraux + P);
  double* t_space = (work + (2*P)); 
  double *tmp,*tmp2;

  int* pivots = work_i;
  int rank = P; 
  int j,i;

  for (j = 0; j < P; j++) { 
    pivots[j] = j + 1;
  }

  __do_copy(X,Y,weights, P, N, t_space,M);
  tmp = X;
  X = t_space;
  t_space = tmp;

  tmp = Y;
  Y = X+N*P;
  tmp2 = weights;
  weights = X+N*(P+M);

  __do_linear(X, N, P, Y, weights, tolerance, estimates, residuals, effects, 
	      &rank, pivots, qraux, work, loss, error); 

  X = t_space;
  Y = tmp;
  weights = tmp2;
}

double __get_multinomial_likelihood(double* betas, int N, int P, int M, double* X, double* Y,
				    double* weights, double* mus) 
{
  int n, m, p;
  double denom, sm, eta, yy;
  double ll = 0;

  for (n = 0; n < N; n++) {
    denom = sm = yy = 0;
    
    for (m = 0; m < M; m++) {
      eta = 0;
      for (p = 0; p < P; p++) {
	eta += (X[p*N + n] * betas[p + m*P]);
      }
      mus[m] = ((eta < MTHRESH) ? DOUBLE_EPS :
		((eta > THRESH) ? INVEPS : exp(eta)));
      denom += mus[m];
      yy += Y[m*N + n];
    }
  
    for (m = 0; m < M; m++) {
      mus[m] = mus[m]/(1 + denom);
    }

    for (m = 0; m < M; m++) {
      ll += Y[m*N + n] * log(mus[m])*weights[n];
    }
    ll += (1 - yy) * log(1/(1+denom))*weights[n];
    
  }
  return ll;
}


void __do_multinomial(double* X, int N, int P, int M, double* Y, 
		      double tol, double* coefficients, 
		      double* residuals, int* rank, int* pivot,
		      double* W, double* mu, double* XTWX, double* I,
		      double* grad_theta, int max_iter, 
		      double* new_coefficients, double* user_weights, 
		      int* error, double* likelihood)
{
  int M_2 = M*M;
  int P_2 = P*P;
  int MP  = M*P;
  int n,m,mm,p,pp,iter = 0;
  double denom, eta, tmp, sm, abs_diff, ll_current, ll_previous, ll_trial;
  int info = 0,NnonNA=N;
  double alpha = 1.0;

  *error = MODEL_NO_ERROR;

  /** for safety. **/
  for (m = 0; m < MP; m++) {
    coefficients[m] = 0;
  }

  /*Remove Missing Values*/
  for(n = 0; n < N; n++)
    {
      if(ISNA(Y[n]))
	{
	  for(m = 0; m < M; m++)Y[m*N + n]=0.0;
	  user_weights[n]=0.0; /*if one Ycat missing then the others also are missing*/
	}
      for (p = 0; p < P; p++) {
	if(ISNA(X[n + p*N])){
	  user_weights[n]=X[n + p*N]=0.0;
	}
      }
      if(fabs(user_weights[n])<MY_ZERO)
	NnonNA--;
    }


  ll_current = ll_previous = __get_multinomial_likelihood(coefficients, N, P, M, X, Y,
							  user_weights, mu);

  while (1) {

    /** compute eta or eta = XB. and save the exponentiated result
	for the softmax denominator. **/
    for(n = 0; n < N; n++) {
      denom = 0.0;
      for(m = 0; m < M; m++) {
	eta = 0.0; 
	for(p = 0; p < P; p++) {
	  eta = eta + (X[p*N + n] * coefficients[p + m*P]);
	}
	/*mu[m*N + n] = exp(eta);*/
	mu[m*N + n] = ((eta < MTHRESH) ? DOUBLE_EPS :
		       ((eta > THRESH) ? INVEPS : exp(eta)));

	denom += mu[m*N + n];
      }
      
      /** Compute the softmax denominator for each unit. */
      for(m = 0; m < M; m++) {
	/** XXX: need to test whether this is 0 or 1; esp denom and maybe warn ?? **/
	mu[m*N + n] = ((mu[m*N + n])/(1.0 + denom));
      }
      
      /** now compute the block diagonal matrix W. **/
      for(m = 0; m < M; m++) {
	for(mm = 0; mm < M; mm++) {
	  if (mm == m) {
	    tmp = (1 - mu[m*N + n]);
	  }
	  else {
	    tmp = -mu[mm*N + n];
	  }
	  W[n*M_2 + m*M + mm] = tmp*mu[m*N + n]*user_weights[n];
	}
      }
    }
    
    /** Get XTWX */
    for(m = 0; m < M; m++) {
      for(mm = 0; mm < M; mm++) {
	for(p = 0; p < P; p++) {
	  for(pp = 0; pp < P; pp++) {
	    sm = 0;
	    for (n = 0; n < N; n++) {
	      sm += (X[p*N + n] * W[n*M_2 + m*M + mm] * X[pp*N + n]);
	    }
	    XTWX[M*P*p + P*mm + m*P_2*M + pp] = sm;
	  }
	}
      }
    }
    
    /** Now get the gradient. */
    for(m = 0; m < M; m++) {
      for(p = 0; p < P; p++) {
	sm = 0;
	for(n = 0; n < N; n++) {
	  sm += (X[p*N + n] * user_weights[n] * (Y[m*N + n] - mu[m*N + n]));
	}
	grad_theta[m*P + p] = sm;
      }
    }
    
    /** invert the matrix. **/
    F77_CALL(dgesv)(&MP, &MP, XTWX, &MP, pivot, I, &MP, &info); 

    if (info != 0) {
      if (R_MODEL_DEBUG_LEVEL <= R_MODEL_DEBUG_WARN)
	Rprintf("Error: Unable to invert matrix in \n");
      
      *likelihood = R_PosInf;
    }

    alpha = 1.0;

    /** now we step-halve until we get something better. **/
    while (1) {
      /** Multiply the gradient by the inverse hessian and get your coefficient update. */
      for(m = 0; m < MP; m++) {
	sm = 0;
	for(p = 0; p < MP; p++) {
	  sm += (I[m*MP + p] * grad_theta[p]);
	
	  /** set the diagonal matrix back. **/
	  if (m == p) {
	    I[m*MP + p] = 1;
	  }
	  else {
	    I[m*MP + p] = 0;
	  }
	}
	new_coefficients[m] = coefficients[m] + (alpha*sm);
      }
      ll_trial = __get_multinomial_likelihood(new_coefficients, N, P, M, X, Y,
					      user_weights, mu);

      if (ll_trial >= ll_current) {
	ll_current = ll_trial;
	break;
      }
      else if (alpha < tol) {
	break; /** probably should do something more drastic here.**/
      }
      else {
	alpha = alpha/2.0;
      }
    }

    if (iter++ >= 25 ||  fabs(ll_current - ll_previous) < tol) {
      break;
    }
    else {
      ll_previous = ll_current;
      for (m = 0; m < MP; m++) {
	coefficients[m] = new_coefficients[m];
      }
    }
  }

  /**
   * error checking on the way out.
   */
  if (iter > max_iter) {
    *error = MODEL_ERROR_NO_CONVERGENCE;
    *likelihood = R_PosInf;
  }
  else {
    *likelihood = -ll_current; /** we use -log p(.) loss. */
  }

  *pivot=NnonNA;
}


/* taken right from family.c in the R code. */
double y_log_y(double y, double mu)
{
  return (y) ? (y * log(y/mu)) : 0;
}

void __do_copy(double* x,double* y,double* weights, int n, int p, double* xt,int M) {
  int i;
  
/*   if(transpose==0) */
/*     { */
/*       for (i = 0; i < n*p; i++) */
/* 	xt[i] = x[(i / p) + (i % p) * n]; */
/*     } */
/*   else  */
/*     { */
  memcpy(xt,x,sizeof(double)*n*p);
/*     } */
  memcpy(xt+n*p,y,sizeof(double)*p*M); /*p is here the number of observations and n the number of variables*/
  memcpy(xt+(n+M)*p,weights,sizeof(double)*p); /*p is here the number of observations and n the number of variables*/
}

void __do_linear (double* X, int n, int p, double* Y, 
		  double* weights, double tol, double* coefficients, 
		  double* residuals, double* effects, int* rank, 
		  int* pivot, double* qraux, double* work, double* rss,
		  int* fail)
{
  if (R_MODEL_DEBUG_LEVEL <= R_MODEL_DEBUG_LOW)
    Rprintf("entering: __do_linear with: n: %d, p: %d\n", n, p);

  int k = 0, j = 0, i = 0, one = 1;
  double RSS = 0, sroot = 0;
  int two = 2,NnonNA=n;
  
  *fail = MODEL_NO_ERROR;

  /** Really probably don't need to do this, but it is cleaner. **/
  for (i = 0; i < p; i++) {
    pivot[i] = i + 1;
  }
  memset(qraux, 0, sizeof(double)*p);
  memset(effects, 0, sizeof(double)*n);
  memset(work, 0, sizeof(double)*2*p);
  memset(coefficients, 0, sizeof(double)*p);
  
  /* apply weights to the X/Y matrix. */
  for (j = 0; j < n; j++) {
    if(ISNA(Y[j]))weights[j]=0.0;
    else
      {
	for (i = 0; i < p; i++) {
	  if(ISNA(X[j + i*n])){
	    weights[j]=0.0;
	    break;
	  }
	}
    }
    if(fabs(weights[j])<MY_ZERO)
      {
	weights[j]=Y[j]=0.0;
	for (i = 0; i < p; i++)X[j + i*n]=0.0;
	NnonNA--;
      }
    else 
      {
	sroot = sqrt(weights[j]);
	Y[j] = Y[j] * sroot;
	for(i = 0; i < p; i++)
	  X[j + i*n] = X[j + i*n] * sroot;
      }
  }

  F77_CALL(dqrls)(X, &n, &p, Y, &one, &tol, coefficients, residuals, 
		  effects, rank, pivot, qraux, work);

  if (*rank < p) { 
    if (R_MODEL_DEBUG_LEVEL <= R_MODEL_DEBUG_WARN)
      Rprintf("Rank deficiency; setting RSS to infinity.");
    
    *rss = R_PosInf;
    *fail = MODEL_ERROR_RANK_CHANGE;
  }
  else {
    for(k = 0; k < n; k++) {
      RSS += R_pow_di(residuals[k], two);
    }
    *rss = RSS; 
  } 
  
  /* pivot the coefficients. */
  for (j = 0; j < p; j++) {
    work[j] = coefficients[pivot[j] - 1]; 
  }
  for (j = 0; j < p; j++) {
    coefficients[j] = work[j];
  }

  /* unapply weights on the way out. (is my unweighting correct?) NO - X is destroyed by dqrls*/
/*   for (j = 0; j < n; j++) { */
/*     sroot = sqrt(weights[j]); */
/*     for (i = 0; i < p; i++) { */
/*       X[j + i*n] = (sroot != 0) ? (X[j + i*n] / sroot) : 0; */
/*     } */
/*     Y[j] = (sroot != 0) ? (Y[j] / sroot) : 0; */
/*   } */
  
  if (R_MODEL_DEBUG_LEVEL <= R_MODEL_DEBUG_LOW)
    Rprintf("Calculated rss of: %f\n", RSS);

  *pivot=NnonNA;
}

void __do_logistic(double* X, int n, int p, double* Y, 
		   double tol, double* coefficients, 
		   double* residuals, double* effects, int* rank, 
		   int* pivot, double* qraux, double* work, double* dev,
		   double* weights, int max_iter, double* new_coefficients, 
		   double* WX, double* y_shifted, int* error)
{
  if (R_MODEL_DEBUG_LEVEL <= R_MODEL_DEBUG_LOW)
    Rprintf("Entering __do_logistic with: n: %d, p: %d\n", n, p);
  
  int iter = 0, k = 0, j = 0, i = 0; 
  double sum = 0, theta_j = 0, pi_j = 0, var_pi_j = 0;
  double abs_diff = 0.0;
  double DEV = 0.0, DEV_OLD = 0.0;
  int one = 1,NnonNA=n;
  double sqrt_var_pi_j,sqrt_weight_j;

  *error = MODEL_NO_ERROR;
  
  if (tol == 0) { 
    if (R_MODEL_DEBUG_LEVEL <= R_MODEL_DEBUG_WARN)
      Rprintf("tolerance is precisely 0 resetting to something sensible.\n");
    tol = 1e-8;
  }

  /* zero the coefficients in case they haven't been zeroed. */
  memset(coefficients, 0.0, sizeof(double)*p); 

  while (1) {
    if (R_MODEL_DEBUG_LEVEL <= R_MODEL_DEBUG_LOW)
      Rprintf("Iteration: %d\n", iter);
    
    if (iter++ >= max_iter) {
      *error = MODEL_ERROR_NO_CONVERGENCE;
      DEV = R_PosInf;
      break;
    }
    
    DEV = 0.0;
    
    for (j = 0; j < n; j++) {
      if(iter==1)
	{
	  if(ISNA(Y[j]))weights[j]=Y[j]=0.0;
	  for (i = 0; i < p; i++) {
	    if(ISNA(X[j + i*n])){
	      weights[j]=X[j + i*n]=0.0;
	    }
	  }
	  if(fabs(weights[j])<MY_ZERO)
	    NnonNA--;
	}

      theta_j = 0;
      for (i = 0; i < p; i++) {
	if (R_MODEL_DEBUG_LEVEL <= R_MODEL_DEBUG_TRACE)
	  Rprintf("%f * %f\n", X[j + i*n], coefficients[i]);
	
 	theta_j += X[j + i*n] * coefficients[i];
      }
      
      pi_j = x_d_opx((theta_j < MTHRESH) ? DOUBLE_EPS :
		     ((theta_j > THRESH) ? INVEPS : exp(theta_j)));

      var_pi_j = pi_j * (1 - pi_j); 
      sqrt_var_pi_j = sqrt(var_pi_j);
      sqrt_weight_j = sqrt(weights[j]);

      /* avoid dividing by 0, this effectively ignores the row for this fit. */
      y_shifted[j] = (sqrt_var_pi_j > MY_ZERO) ?
	sqrt_weight_j * (Y[j] - pi_j)/sqrt_var_pi_j :
	0;

      DEV += 2*weights[j]*(y_log_y(Y[j], pi_j) + y_log_y(1 - Y[j], 1 - pi_j));
      
      for (i = 0; i < p; i++) {
	WX[j + i*n] = X[j + i*n] * sqrt_var_pi_j * sqrt_weight_j; 
      }
      
      if (R_MODEL_DEBUG_LEVEL <= R_MODEL_DEBUG_TRACE) {
	Rprintf("theta_j: %f\n", theta_j);
	Rprintf("pi_j: %f\n", pi_j);
	Rprintf("weight[%d] = %f\n", j, weights[j]);
	Rprintf("y_shifted[%d] = %f\n", j, y_shifted[j]);
      }
    }

    F77_CALL(dqrls)(WX, &n, &p, y_shifted, &one, &tol, new_coefficients, residuals, 
		    effects, rank, pivot, qraux, work);

    if (*rank < p) {
      *error = MODEL_ERROR_RANK_CHANGE;
      DEV = R_PosInf;
      break;
    }

    abs_diff = 0.0;

    for (j = 0; j < p; j++) {
      /* converge based on the absolute difference in the coefficients. */
      /* abs_diff += fabs(new_coefficients[j]); */
    
      /* pivot the coefficients. */
      coefficients[j] += new_coefficients[pivot[j] - 1]; 

      /* zero the coefficients. */
      new_coefficients[pivot[j] - 1] = 0;

      /* reset the pivots. */
      pivot[j] = j + 1; 
    }
    
    if (R_MODEL_DEBUG_LEVEL <= R_MODEL_DEBUG_LOW)
      Rprintf("deviance: %g\n", DEV);

    if (fabs(DEV_OLD - DEV)/(0.1 + fabs(DEV)) < tol)
      break;

    /** R uses the difference in deviance to decide when to exit. */
    DEV_OLD = DEV;
  }
  
  /* set some things on the way out. */
  *dev = DEV;

  *pivot=NnonNA;
}


