#ifndef ATTACH_H
#include "Attach.h"
#endif

// Declaration of the Runge-Kutta method, order four, for first and second order ODEs
// Code for second order ODEs will also work for two coupled ODEs
// R. Sheehan 2 - 3 - 2014

void ode_solver::RK4_First_Order_ODE(double start, double end, int N_steps, double alpha, double (*f)(double x, double y),  double *x_vals, double *sol_vals)
{
	// Runge-Kutta Method, Order Four, for single ODE
	// the ode to be solved is y'(x) = f(x , y)
	// calculation takes place at Nsteps in the domain start <= x <= end
	// the initial condition for y is alpha
	// the positions at which the solution is computed are stored in x_vals
	// the solution values are stored in sol_vals

	// 1. Compute the step size that will be used
	double h = (end - start) / (static_cast<double>(N_steps - 1)); 
	double half_f = 0.5*h; 

	// 2. Declare some dummy variables that you can use in the calculation
	double pos, y_old, y_new;
	double k1, k2, k3, k4, rk_cons_sum; // constants to be used in computing ynew

	// 3. Loop over the position in [start, end] and compute the solution using RK4
	
	// initialise the calculation
	pos = start; 
	y_new = y_old = alpha; 

	cout<<"x = "<<pos<<" , "<<"y = "<<y_new<<endl;

	// store the solution for the first step
	x_vals[1] = pos; 	sol_vals[1] = y_new; 

	for(int j=2; j<=N_steps; j++){
		// update the position variable
		pos = start + (j-1)*h;

		// compute the Runge-Kutta constants
		k1 = h*(*f)(pos, y_old); 
		
		k2 = h*(*f)(pos+half_f, y_old + 0.5*k1); 

		k3 = h*(*f)(pos+half_f, y_old + 0.5*k2); 

		k4 = h*(*f)(pos+h, y_old + k3); 

		// compute the sum over the Runge-Kutta constants
		rk_cons_sum = ( k1 + 2.0*k2 + 2.0*k3 + k4 ); 

		// compute the solution at the next step
		y_new = y_old + (rk_cons_sum/6.0); 

		if(j%10 == 0){
			cout<<"x = "<<pos<<" , "<<"y = "<<y_new<<endl;
		}

		// store the compute solution
		x_vals[j] = pos; 	sol_vals[j] = y_new; 

		// update y_old for the next calculation
		y_old = y_new;
	}

	cout<<"Solution computed using RK4\n";
}


void ode_solver::RK4_Coupled_Pair(double start, double end, int N_steps, double alpha_1, double alpha_2, double (*f1)(double x, double u1, double u2), double (*f2)(double x, double u1, double u2),  double *x_vals, double *sol_1_vals, double *sol_2_vals)
{
	// Runge-Kutta Method, Order Four, for a pair of first order ODEs	
	// the ode to be solved is either the second order ode y''(x) = F(x, y)
	// or the pair of coupled first order ODEs u1'(x) = f1(x, u1, u2) and u2'(x) = f2(x, u1, u2)
	// In the case of a second order equation it is assumed that the user has re-written the second order ODE
	// as a pair of coupled first order equations. 
	// calculation takes place at Nsteps in the domain start <= x <= end
	// the initial condition for u1 is alpha_1
	// the initial condition for u2 is alpha_2
	// the positions at which the solution is computed are stored in x_vals
	// the solution values for u1 are stored in sol_vals_1
	// the solution values for u2 are stored in sol_vals_2
	
	// 1. Compute the step size that will be used
	double h = (end - start) / (static_cast<double>(N_steps - 1)); 
	double half_f = 0.5*h; 

	// 2. Declare some dummy variables that you can use in the calculation
	double pos, u1_old, u1_new, u2_old, u2_new;
	double k11, k21, k31, k41, rk_cons_sum_1; // constants to be used in computing u1_new
	double k12, k22, k32, k42, rk_cons_sum_2; // constants to be used in computing u2_new

	// 3. Loop over the position in [start, end] and compute the solution using RK4
	
	// initialise the calculation
	pos = start; 
	u1_new = u1_old = alpha_1; 
	u2_new = u2_old = alpha_2; 

	cout<<"x = "<<pos<<" , "<<"u1 = "<<u1_new<<" , "<<"u2 = "<<u2_new<<endl;

	// store the solution for the first step
	x_vals[1] = pos; sol_1_vals[1] = u1_new; sol_2_vals[1] = u2_new; 

	for(int j=2; j<=N_steps; j++){
		// update the position variable
		pos = start + (j-1)*h;

		// compute the Runge-Kutta constants
		k11 = h*(*f1)(pos, u1_old, u2_old); 
		k12 = h*(*f2)(pos, u1_old, u2_old); 
		
		k21 = h*(*f1)(pos+half_f, u1_old + 0.5*k11, u2_old + 0.5*k12); 
		k22 = h*(*f2)(pos+half_f, u1_old + 0.5*k11, u2_old + 0.5*k12); 

		k31 = h*(*f1)(pos+half_f, u1_old + 0.5*k21, u2_old + 0.5*k22); 
		k32 = h*(*f2)(pos+half_f, u1_old + 0.5*k21, u2_old + 0.5*k22); 

		k41 = h*(*f1)(pos+h, u1_old + k31, u2_old + k32); 
		k42 = h*(*f2)(pos+h, u1_old + k31, u2_old + k32); 

		// compute the sum over the Runge-Kutta constants
		rk_cons_sum_1 = ( k11 + 2.0*k21 + 2.0*k31 + k41 ); 
		rk_cons_sum_2 = ( k12 + 2.0*k22 + 2.0*k32 + k42 ); 

		// compute the solution at the next step
		u1_new = u1_old + (rk_cons_sum_1/6.0); 

		u2_new = u2_old + (rk_cons_sum_2/6.0); 

		if(j%10 == 0){
			cout<<"x = "<<pos<<" , "<<"u1 = "<<u1_new<<" , "<<"u2 = "<<u2_new<<endl;
		}

		// store the compute solution
		x_vals[j] = pos; sol_1_vals[j] = u1_new; sol_2_vals[j] = u2_new;

		// update y_old for the next calculation
		u1_old = u1_new; 
		u2_old = u2_new; 
	}
	
	cout<<"Solution computed using RK4\n"; 

}