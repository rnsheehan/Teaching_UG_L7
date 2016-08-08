#ifndef ODE_SOLVERS_H
#define ODE_SOLVERS_H

// Declaration of the Runge-Kutta method, order four, for first and second order ODEs
// Code for second order ODEs will also work for two coupled ODEs
// R. Sheehan 2 - 3 - 2014

namespace ode_solver{

	// Runge-Kutta Method, Order Four, for single ODE
	void RK4_First_Order_ODE(double start, double end, int N_steps, double alpha, double (*f)(double x, double y), double *x_vals, double *sol_vals); 

	// Runge-Kutta Method, Order Four, for a pair of first order ODEs
	void RK4_Coupled_Pair(double start, double end, int N_steps, double alpha_1, double alpha_2, double (*f1)(double x, double u1, double u2), double (*f2)(double x, double u1, double u2),  double *x_vals, double *sol_1_vals, double *sol_2_vals); 

}

#endif