/*
    GSL Interpolation Wrapper by Pak Ki Henry Tsang 
    
    Jan-2021
    
    this header provides a simple C++ class interface to the GNU Scientific Library (GSL) interpolation routines
    
    Example: to do cubic spline interpolation with Num_points long data x_array, y_array:
    
    	#include "interpolate.h"
    	...
    	int Num_points=6;
    	double[6] x_array={0.0,0.2,0.4,0.6,0.8,1.0};
    	double[6] y_array={0.0,0.5,0.25,0.0.4,0.6,1.0};
    	interpolate interpolating_function(x_array, y_array, Num_points , interpolate::splines::cubic ); //initialing interpolation object doing cubic spline
    	double x_interpl = interpolating_function.eval(0.5); //evaluate the interpolated function at x=0.5
    
    The possible choices of splines are 
    
		  1. linear
		  2. polynomial
		  3. cubic
		  4. cubic_periodic
		  5. akima 
		  6. akima_periodic
		  7. steffen
    
    and one specifies them as interpolate::splines::linear, interpolate::splines::cubic , ... 
*/

//GSL headers
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

class interpolate{

  public:
  
		enum class splines{ linear,	polynomial, cubic, cubic_periodic, akima, akima_periodic, steffen };
		
  	interpolate(double* x_array,double* y_array,size_t size,const interpolate::splines spline_type){ //constructor

			//Local variables
			this->size=size;
			
			//possible splines - 1. gsl_interp_linear 2. gsl_interp_polynomial 3.gsl_interp_cspline 4. gsl_interp_cspline_periodic 5.gsl_interp_akima 6. gsl_interp_akima_periodic 7.gsl_interp_steffen
			
			switch (spline_type) {
				case splines::linear:
					spline = gsl_spline_alloc (gsl_interp_linear, size);
					break;
				case splines::polynomial:
					spline = gsl_spline_alloc (gsl_interp_polynomial, size);
					break;
				case splines::cubic:
					spline = gsl_spline_alloc (gsl_interp_cspline, size);
					break;
				case splines::cubic_periodic:
					spline = gsl_spline_alloc (gsl_interp_cspline_periodic, size);
					break;
				case splines::akima:
					spline = gsl_spline_alloc (gsl_interp_akima, size);
					break;
				case splines::akima_periodic:
					spline = gsl_spline_alloc (gsl_interp_akima_periodic, size);
					break;
				case splines::steffen:
					spline = gsl_spline_alloc (gsl_interp_steffen, size);
					break;
				default:
					spline = gsl_spline_alloc (gsl_interp_linear, size);
					break;
			}
			
			acc = gsl_interp_accel_alloc();
			gsl_spline_init (spline, x_array, y_array, size);
		};
		
  	~interpolate(){ gsl_spline_free (spline);	gsl_interp_accel_free (acc); }; //Destructor
  	
  	//Standard spline evaluation, will abort if x is outside of interpolation range
  	double eval(double x){ return gsl_interp_eval(spline->interp, spline->x, spline->y, x, acc); };

  	double eval_ignore_bounds(double x){ //Custom implementation that returns zero on out-of-bound x
  		if (x < spline->interp->xmin || x > spline->interp->xmax) { return 0.0; }
  		else { return gsl_interp_eval(spline->interp, spline->x, spline->y, x, acc); }
  	};
  	
  	//derivatives, only works inside interpolation range
  	double eval_deriv(double x){ gsl_interp_eval_deriv(spline->interp, spline->x, spline->y, x, acc); };
  	double eval_deriv2(double x) { gsl_interp_eval_deriv2(spline->interp, spline->x, spline->y, x, acc); };
  	
  	//integral (trapezoidal rule), only works inside interpolation range
		double eval_integ(double a, double b){ return gsl_interp_eval_integ(spline->interp, spline->x, spline->y, a, b, acc); };

	private:
	
		size_t size;
		gsl_interp_accel *acc;
		gsl_spline *spline;

};


