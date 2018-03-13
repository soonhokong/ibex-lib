/* ============================================================================
 * D Y N I B E X - Example for a second validated simulation of a DAE with additive constraints
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

const int n= 3;

IntervalVector yinit(3);
yinit[0]= Interval(5.0);
yinit[1]= Interval(1.0);
yinit[2]= Interval(0.0);

Variable y(3);

IntervalVector xinit(2);
xinit[0]= Interval(-1.0);
xinit[1]= Interval(-1e-14,1e-14);


Variable x(2);

Function ydot = Function(y,x,Return(-y[2]*y[1]-(1+y[2])*x[0],
				    y[2]*y[0]-(1+y[2])*x[1],Interval(1)));

Function g = Function(y,x,Return((y[0]-x[1])/5.0 - cos(sqr(y[2])/2.0),
				  (y[1]+x[0])/5.0 - sin(sqr(y[2])/2.0)));

NumConstraint csp1(y,x,sqr(x[0]) + sqr(x[1]) - 1.0 = 0);

Array<NumConstraint> csp(csp1);


ivp_dae_h1 problem = ivp_dae_h1(ydot,g,0.0,yinit,xinit,csp);

simulation simu = simulation(&problem,0.5,RADAU3_DAE, 1e-14);

simu.run_simulation();

return 0;

}