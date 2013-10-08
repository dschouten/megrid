
// #include <stdio.h>

#include <cmath>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_poly.h>

#include "logger.hh"

// supernewton method from libAGF ( (C) Peter Mills)
//////////////////////////////////////////////////////////////////////////////
// note: for smooth ("well-behaved") functions with local non-zero third
// moments, the algorithm has a tendency to land repeatedly on one side
// of the root, slowing down convergence.  This issue needs to be addressed.
//
// minimization function which brackets the root
// and then approximates root by fitting a third-
// order polynomial
//////////////////////////////////////////////////////////////////////////////

template <class FLOAT_T>
FLOAT_T supernewton(void (*funcd) (FLOAT_T, void *, FLOAT_T *, FLOAT_T *),
		    void* params, // function parameters
		    FLOAT_T xa,	  // first bracket
		    FLOAT_T xb,	  // second bracket
		    FLOAT_T xtol, // desired tolerance in x direction
		    FLOAT_T ytol, // desired tolerance in y direction
		    long maxiter, // maximum number of iterations
		    long& err,	  // error code, positive for success
		    FLOAT_T &ya,  // to avoid re-calculating these items
		    FLOAT_T &dydxa, 
		    FLOAT_T yb,
		    FLOAT_T dydxb)
{
  //function evaluations at the brackets:
  /*
    FLOAT_T ya;
    FLOAT_T dydxa;
    FLOAT_T yb;
    FLOAT_T dydxb;
  */
  
  //polynomial coefficients:
  double a(0), b(0), c(0), d(0);	//might as well carry over from double to double
  
  FLOAT_T x0((xa+xb)/2); //solution
  
  FLOAT_T y0;		//solution value for y
  FLOAT_T dydx0;	//solution value for dydx
  
  //intermediates in calculation:
  FLOAT_T xa_2, xa_3, xb_2, xb_3;
  long nroot;
  
  //roots of polynomial:
  double x0_1, x0_2, x0_3;	//these need to be double because they're 
  				//returned from a GSL routine
  FLOAT_T x0_4;
  
  //number of iterations:
  long i;
  
  //x and y error:
  FLOAT_T xerr, yerr;

  //FLOAT_T xerr1, xerr2;
  
  //approximate slope of the function:
  //FLOAT_T m;
  
  gsl_vector *yp;
  gsl_matrix *A;
  
  int gsl_status;	//error state from GSL calls
  
  yp=gsl_vector_alloc(4);
  A=gsl_matrix_alloc(4, 4);
  
  err=0;
  i=0;

  while( true )
  {
    i++;
  
    /*
      printf("Brackets: [%g, %g]\n", xa, xb);
      printf("          [%g, %g]\n", ya, yb);
      printf("          [%g, %g]\n", dydxa, dydxb);
    */
      
    //populate the matrix and solution vector:
    xa_2=xa*xa;
    xa_3=xa_2*xa;
    xb_2=xb*xb;
    xb_3=xb_2*xb;
    
    gsl_matrix_set(A, 0, 0, 1);
    gsl_matrix_set(A, 0, 1, xa);
    gsl_matrix_set(A, 0, 2, xa_2);
    gsl_matrix_set(A, 0, 3, xa_3);
    
    gsl_vector_set(yp, 0, ya);
    
    gsl_matrix_set(A, 1, 0, 1);
    gsl_matrix_set(A, 1, 1, xb);
    gsl_matrix_set(A, 1, 2, xb_2);
    gsl_matrix_set(A, 1, 3, xb_3);
    
    gsl_vector_set(yp, 1, yb);
    
    gsl_matrix_set(A, 2, 0, 0);
    gsl_matrix_set(A, 2, 1, 1);
    gsl_matrix_set(A, 2, 2, 2*xa);
    gsl_matrix_set(A, 2, 3, 3*xa_2);
    
    gsl_vector_set(yp, 2, dydxa);
    
    gsl_matrix_set(A, 3, 0, 0);
    gsl_matrix_set(A, 3, 1, 1);
    gsl_matrix_set(A, 3, 2, 2*xb);
    gsl_matrix_set(A, 3, 3, 3*xb_2);
    
    gsl_vector_set(yp, 3, dydxb);
    
    //solve using Householder transformations:
    gsl_status=gsl_linalg_HH_svx(A, yp);
    
    if (gsl_status != 0)
    {
      logger::log() << msg::ERROR << "Supernewton: gsl_linalg_HH_svx returned the following error:";
      logger::log() << msg::ERROR << gsl_strerror(gsl_status) << " : " << gsl_status;
      //use a simple bisection step:
      x0=(xa+xb)/2;
      // err = -i;
      // break;
    }
    else 
    {      
      a=gsl_vector_get(yp, 0);
      b=gsl_vector_get(yp, 1);
      c=gsl_vector_get(yp, 2);
      d=gsl_vector_get(yp, 3);
      
      //printf("Polynomial coeffs.:%f, %f, %f, %f\n", a, b, c, d);
      
      //solve the cubic:
      nroot=gsl_poly_solve_cubic(c/d, b/d, a/d, &x0_1, &x0_2, &x0_3);
      
      if (nroot==1)
      {
        //if there is only one root AND it's finite AND it advances the solution
        //we use that one:
        x0=(xa+xb)/2;
        if (std::isfinite(x0_1))
	{
          if (std::fabs(x0_1-x0) < std::fabs((xa+xb)/2)) 
	    x0=x0_1;
        }
        //printf("Root: %f\n", x0);
      }
      else if (nroot==3) 
      {
        FLOAT_T xdiff, xdiffmin;
	
        x0_4=(xa+xb)/2;
        //printf("Roots: %f, %f, %f\n", x0_1, x0_2, x0_3);
	
        x0=x0_4;
        xdiffmin=std::fabs((xa-xb)/2);
        if (std::isfinite(x0_1)) 
	{
          xdiff=std::fabs(x0_1-x0_4);
          if (xdiff < xdiffmin)
	  {
            xdiffmin=xdiff;
            x0=x0_1;
          }
        }
        if (std::isfinite(x0_2))
	{
          xdiff=std::fabs(x0_2-x0_4);
          if (xdiff < xdiffmin) 
	  {
            xdiffmin=xdiff;
            x0=x0_2;
          }
        }
        if (std::isfinite(x0_3)) 
	{
          xdiff=std::fabs(x0_3-x0_4);
          if (xdiff < xdiffmin) 
	  {
            xdiffmin=xdiff;
            x0=x0_3;
          }
        }
      } 
      else
      {
	logger::log() << msg::ERROR << "Supernewton: gsl_poly_solve_cubic returned the following error:";
	logger::log() << msg::ERROR << gsl_strerror(gsl_status) << " : " << gsl_status;
        err=-i;
        break;
      }
    }
    
    //evaluate the function at the approximated root
    (*funcd) (x0, params, &y0, &dydx0);
    
    //m=(yb-ya)/(xb-xa);
    yerr=std::fabs(y0);
    
    //approximate the x error from the y error and approx. slope of func.:
    //xerr=std::fabs(y0/m);
    
    //xerr1=std::fabs(x0-xa);
    //xerr2=std::fabs(x0-xb);
    
    xerr=std::fabs(xa-xb);
    
    //if (xerr1 < xerr2) xerr=xerr1; else xerr=xerr2;
    
    //test for convergence:
    // ( *note* this convergence test is designed for speed, NOT accuracy )
    if (yerr < ytol || xerr < xtol) {
      //if (yerr < ytol) {
      err=i;
      break;
    }
    
    if (i > maxiter)
    {
      logger::log() << msg::ERROR << "Maximum number of iterations exceeded in supernewton:";
      err = -maxiter;
      // printf("\tBrackets: [%f, %f]\n", xa, xb);
      // printf("\tPolynomial coeffs.:%f, %f, %f, %f\n", a, b, c, d);
      // printf("\tRoots: %f, %f, %f\n", x0_1, x0_2, x0_3);
      break;
    }
    
    //rebracket the "true" root:
    if (y0*ya > 0) 
    {
      xa=x0;
      ya=y0;
      dydxa=dydx0;
    }
    else 
    {
      xb=x0;
      yb=y0;
      dydxb=dydx0;
    }    
  }
  
  //place the function values in the fourth and third last parameters:
  ya=y0;
  dydxa=dydx0;
  
  //printf("Supernewton: %d iterations required to reach convergence\n", i);
  
  //clean up:
  gsl_matrix_free(A);
  gsl_vector_free(yp);
  
  return x0;
}

template float supernewton<float>(void (*funcd) (float, void *, float *, float *),
				  void *params,
				  float xa,
				  float xb,
				  float xtol,
				  float ytol,
				  long maxiter,
				  long &err,
				  float &ya,
				  float &dydxa,
				  float yb,
				  float dydxb);

template double supernewton<double>(void (*funcd) (double, void *, double *, double *),
				    void *params,
				    double xa,
				    double xb,
				    double xtol,
				    double ytol,
				    long maxiter,
				    long &err,
				    double &ya,
				    double &dydxa,
				    double yb,
				    double dydxb);
