/* ============================================================================
 * D Y N I B E X - Example for a classical DAE : the pendulum
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Julien Alexandre dit Sandretto and Alexandre Chapoutot
 * Created     : Oct, 2015
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, 
 * THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */
#include "ibex.h"

using namespace ibex;

int main(){

const int n= 4;

IntervalVector yinit(n);
yinit[0]= Interval(1.0);
yinit[1]= 0.0;
yinit[2]= 0.0;
yinit[3]= 0.0;


Variable y(n);

IntervalVector xinit(1);
xinit[0]= Interval(-0.01,0.01);

Variable x(1);

Interval m = Interval(1.0);
Interval g = Interval(9.81);
Interval l = Interval(1.0);

Function ydot = Function(y,x,Return(y[2],
				    y[3],
				    -y[0]*x[0]/m,
				    -(y[1]*x[0] + g)/m));

Function fg = Function(y,x,m*((y[2]*y[2]) + (y[3]*y[3]))-g*y[1] -l*l*x[0]);

//for this version, additive constraints are required
NumConstraint csp1(y,sqr(y[0])+sqr(y[1])-sqr(l) = 0);
NumConstraint csp2(y,y[0]*y[2] + y[1]*y[3]=0);
NumConstraint csp3(y,x,m*(sqr(y[2]) + sqr(y[3]))-g*y[1] -l*l*x[0]=0);
NumConstraint csp4(y,l - sqrt(sqr(y[0])+sqr(y[1])) = 0);

Array<NumConstraint> csp(csp1,csp2,csp3,csp4);

ivp_dae_h1 prob_dae = ivp_dae_h1(ydot,fg,0.0,yinit,xinit,csp);
simulation simu = simulation(&prob_dae,1.6,RADAU3_DAE,1e-15);

simu.run_simulation();


return 0;

}