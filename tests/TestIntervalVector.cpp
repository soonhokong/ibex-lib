/* ============================================================================
 * I B E X - Interval Vector Tests
 * ============================================================================
 * Copyright   : Ecole des Mines de Nantes (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Gilles Chabert
 * Created     : Mar 2, 2012
 * ---------------------------------------------------------------------------- */

#include "TestIntervalVector.h"
#include "ibex_Interval.h"
#include "utils.h"

using namespace std;

void TestIntervalVector::cons01() {
	IntervalVector x(2);
	x[0]=Interval::ALL_REALS;
	x[1]=Interval::ALL_REALS;
	check(x,IntervalVector(2));
	check(x,IntervalVector(x));
	check(x,IntervalVector(2)=x);
}

void TestIntervalVector::cons02() {
	IntervalVector x(2);
	x[0]=Interval(0,1);
	x[1]=Interval(0,1);
	check(x,IntervalVector(2,Interval(0,1)));
	check(x,IntervalVector(x));
	check(x,IntervalVector(2)=x);
}

void TestIntervalVector::cons03() {
	IntervalVector x(2);
	x[0]=Interval(0,1);
	x[1]=Interval(2,3);
	check(x,IntervalVector(x));
	check(x,IntervalVector(2)=x);
}

void TestIntervalVector::cons04() {
	double bounds[][2] = {{0,1},{2,3}};
	IntervalVector x(2);
	x[0]=Interval(0,1);
	x[1]=Interval(2,3);
	check(x,IntervalVector(2,bounds));
	check(x,IntervalVector(2)=x);
}

void TestIntervalVector::cons05() {
	IntervalVector x(2);
	x[0].set_empty();
	x[1].set_empty();
	//check(x,IntervalVector::empty(2));
	CPPUNIT_ASSERT(x.is_empty());
	//check(x,IntervalVector(2)=x);
	CPPUNIT_ASSERT((IntervalVector(2)=x).is_empty());
}

void TestIntervalVector::set_empty01() {
	IntervalVector x(2);
	CPPUNIT_ASSERT(!x.is_empty());
	x.set_empty();
	CPPUNIT_ASSERT(x.is_empty());
}


void TestIntervalVector::is_empty01() {
	CPPUNIT_ASSERT(IntervalVector::empty(2).is_empty());
}

void TestIntervalVector::is_empty02() {
	CPPUNIT_ASSERT(!IntervalVector(2).is_empty());
}

void TestIntervalVector::resize01() {
	IntervalVector x(1);
	x[0]=Interval(1,2);
	x.resize(3);
	CPPUNIT_ASSERT(x.size()==3);
	check(x[0],Interval(1,2));
	check(x[1],Interval::ALL_REALS);
	check(x[2],Interval::ALL_REALS);
}

void TestIntervalVector::resize02() {
	IntervalVector x(1);
	x[0]=Interval(1,2);
	x.resize(1);
	CPPUNIT_ASSERT(x.size()==1);
	check(x[0],Interval(1,2));
}

void TestIntervalVector::resize03() {
	IntervalVector x(2);
	x[0]=Interval(1,2);
	x.set_empty();
	x.resize(3);
	CPPUNIT_ASSERT(x.size()==3);
	CPPUNIT_ASSERT(x.is_empty());
	CPPUNIT_ASSERT(x[2]==Interval::ALL_REALS);
}

void TestIntervalVector::resize04() {
	IntervalVector x(5);
	x[0]=Interval(1,2);
	x[1]=Interval(3,4);
	x.resize(2);
	CPPUNIT_ASSERT(x.size()==2);
	check(x[0],Interval(1,2));
	check(x[1],Interval(3,4));
}

static double _x[][2]={{0,1},{2,3},{4,5}};

void TestIntervalVector::subvector01() {
	double _x01[][2]={{0,1},{2,3}};
	CPPUNIT_ASSERT(IntervalVector(3,_x).subvector(0,1)==IntervalVector(2,_x01));
}

void TestIntervalVector::subvector02() {
	double _x12[][2]={{2,3},{4,5}};
	CPPUNIT_ASSERT(IntervalVector(3,_x).subvector(1,2)==IntervalVector(2,_x12));
}

void TestIntervalVector::subvector03() {
	double _x11[][2]={{2,3}};
	CPPUNIT_ASSERT(IntervalVector(3,_x).subvector(1,1)==IntervalVector(1,_x11));
}

void TestIntervalVector::subvector04() {
	double _x22[][2]={{4,5}};
	CPPUNIT_ASSERT(IntervalVector(3,_x).subvector(2,2)==IntervalVector(1,_x22));
}

void TestIntervalVector::subvector05() {
	CPPUNIT_ASSERT(IntervalVector(3,_x).subvector(0,2)==IntervalVector(3,_x));
}

void TestIntervalVector::subvector06() {
	CPPUNIT_ASSERT(IntervalVector::empty(3).subvector(1,2).is_empty());
}

void TestIntervalVector::cart_prod01() {
	CPPUNIT_ASSERT(cart_prod(IntervalVector(3,_x),IntervalVector::empty(3)).is_empty());
	CPPUNIT_ASSERT(cart_prod(IntervalVector::empty(3),IntervalVector(3,_x)).is_empty());
}

void TestIntervalVector::inter01() {
	double _x1[][2]={{0,2},{4,6}};
	double _x2[][2]={{1,3},{5,7}};
	double _res[][2]={{1,2},{5,6}};
	CPPUNIT_ASSERT(((IntervalVector(2,_x1)) &=IntervalVector(2,_x2))==IntervalVector(2,_res));
	CPPUNIT_ASSERT(((IntervalVector(2,_x1)) & IntervalVector(2,_x2))==IntervalVector(2,_res));
}

void TestIntervalVector::inter02() {
	double _x1[][2]={{0,2},{4,6}};
	double _x2[][2]={{1,3},{7,8}};
	CPPUNIT_ASSERT(((IntervalVector(2,_x1)) &=IntervalVector(2,_x2)).is_empty());
	CPPUNIT_ASSERT(((IntervalVector(2,_x1)) & IntervalVector(2,_x2)).is_empty());
}

void TestIntervalVector::inter03() {
	double _x1[][2]={{0,2},{4,6}};
	CPPUNIT_ASSERT(((IntervalVector(2,_x1)) &=IntervalVector::empty(2)).is_empty());
	CPPUNIT_ASSERT(((IntervalVector(2,_x1)) & IntervalVector::empty(2)).is_empty());
}

void TestIntervalVector::inter04() {
	double _x1[][2]={{0,2},{4,6}};
	CPPUNIT_ASSERT(((IntervalVector::empty(2)) &=IntervalVector(2,_x1)).is_empty());
	CPPUNIT_ASSERT(((IntervalVector::empty(2)) & IntervalVector(2,_x1)).is_empty());
}

void TestIntervalVector::hull01() {
	double _x1[][2]={{0,1},{4,5}};
	double _x2[][2]={{2,3},{6,7}};
	double _res[][2]={{0,3},{4,7}};
	CPPUNIT_ASSERT(((IntervalVector(2,_x1)) |=IntervalVector(2,_x2))==IntervalVector(2,_res));
	CPPUNIT_ASSERT(((IntervalVector(2,_x1)) | IntervalVector(2,_x2))==IntervalVector(2,_res));
}

void TestIntervalVector::hull02() {
	double _x1[][2]={{0,1},{4,5}};
	IntervalVector x1(2,_x1);
	CPPUNIT_ASSERT((x1 |= x1)==x1);
	CPPUNIT_ASSERT((x1 | x1)==x1);
}

void TestIntervalVector::hull03() {
	double _x1[][2]={{0,2},{4,6}};
	IntervalVector x1(2,_x1);
	CPPUNIT_ASSERT((x1 |=IntervalVector::empty(2))==x1);
	CPPUNIT_ASSERT((x1 | IntervalVector::empty(2))==x1);
}

void TestIntervalVector::hull04() {
	double _x1[][2]={{0,2},{4,6}};
	IntervalVector x1(2,_x1);
	CPPUNIT_ASSERT(((IntervalVector::empty(2)) |= x1)==x1);
	CPPUNIT_ASSERT(((IntervalVector::empty(2)) | x1)==x1);
}

void TestIntervalVector::eq01() {
	IntervalVector x(3,_x);
	CPPUNIT_ASSERT(x==x);
	CPPUNIT_ASSERT(!(x!=x));
}

void TestIntervalVector::eq02() {
	IntervalVector x(3,_x);
	double _x01[][2]={{0,1},{2,3}};
	IntervalVector x1(2,_x01);
	CPPUNIT_ASSERT(!(x==x1));
	CPPUNIT_ASSERT(x!=x1);
}

void TestIntervalVector::eq03() {
	double _x1[][2]={{0,1},{4,5}};
	double _x2[][2]={{2,3},{6,7}};
	IntervalVector x1(2,_x1);
	IntervalVector x2(2,_x2);
	x1.set_empty();
	x2.set_empty();
	CPPUNIT_ASSERT(x1==x2);
	CPPUNIT_ASSERT(!(x1!=x2));
}

void TestIntervalVector::eq04() {
	CPPUNIT_ASSERT(IntervalVector::empty(2)==IntervalVector::empty(2));
	CPPUNIT_ASSERT(IntervalVector::empty(2)!=IntervalVector::empty(3));
	IntervalVector x(2);
	x.set_empty();
	CPPUNIT_ASSERT(IntervalVector::empty(2)==x);
}

void TestIntervalVector::mid01() {
	IntervalVector x(3,_x);
	IntervalVector m=x.mid();
	check(m[0],0.5);
	check(m[1],2.5);
	check(m[2],4.5);
}

void TestIntervalVector::is_flat01() {
	CPPUNIT_ASSERT(!IntervalVector(3,_x).is_flat());
}

void TestIntervalVector::is_flat02() {
	CPPUNIT_ASSERT(IntervalVector::empty(3).is_flat());
}

void TestIntervalVector::is_flat03() {
	CPPUNIT_ASSERT(IntervalVector(1,Interval(0,0)).is_flat());
	CPPUNIT_ASSERT(!IntervalVector(1,Interval(0,1)).is_flat());
}

void TestIntervalVector::is_flat04() {
	double _x1[][2]={{0,1},{2,2},{3,4}};
	CPPUNIT_ASSERT(IntervalVector(3,_x1).is_flat());
}

void TestIntervalVector::is_flat05() {
	double _x1[][2]={{0,1},{2,3},{4,4}};
	CPPUNIT_ASSERT(IntervalVector(3,_x1).is_flat());
}

void TestIntervalVector::is_unbounded01() {
	CPPUNIT_ASSERT(!IntervalVector::empty(3).is_unbounded());
}

void TestIntervalVector::is_unbounded02() {
	double _x1[][2]={{0,1},{0,2},{NEG_INFINITY,0}};
	CPPUNIT_ASSERT(IntervalVector(3,_x1).is_unbounded());
}

void TestIntervalVector::is_unbounded03() {
	double _x1[][2]={{0,1},{0,2}};
	CPPUNIT_ASSERT(!IntervalVector(2,_x1).is_unbounded());
}

void TestIntervalVector::is_unbounded04() {
	CPPUNIT_ASSERT(IntervalVector(1).is_unbounded());
}

void TestIntervalVector::is_subset01() {
	double _x1[][2]={{0,2},{2,4}};
	double _x2[][2]={{0,1},{3,4}};
	IntervalVector x1(2,_x1);
	IntervalVector x2(2,_x2);

	CPPUNIT_ASSERT(x1.is_superset(x2));
	CPPUNIT_ASSERT(x2.is_subset(x1));
	CPPUNIT_ASSERT(x1.is_strict_superset(x2));
	CPPUNIT_ASSERT(!x2.is_strict_interior_subset(x1));
}

void TestIntervalVector::is_subset02() {
	double _x1[][2]={{0,2},{2,4}};
	double _x2[][2]={{1,1},{3,4}};
	IntervalVector x1(2,_x1);
	IntervalVector x2(2,_x2);

	CPPUNIT_ASSERT(x1.is_superset(x2));
	CPPUNIT_ASSERT(x2.is_subset(x1));
	CPPUNIT_ASSERT(x1.is_strict_superset(x2));
	CPPUNIT_ASSERT(!x2.is_strict_interior_subset(x1));
}

void TestIntervalVector::is_subset03() {
	double _x1[][2]={{0,2},{2,4}};
	double _x2[][2]={{0,1},{3,3}};
	IntervalVector x1(2,_x1);
	IntervalVector x2(2,_x2);

	CPPUNIT_ASSERT(x1.is_superset(x2));
	CPPUNIT_ASSERT(x2.is_subset(x1));
	CPPUNIT_ASSERT(x1.is_strict_superset(x2));
	CPPUNIT_ASSERT(!x2.is_strict_interior_subset(x1));
}

void TestIntervalVector::is_subset04() {
	double _x1[][2]={{0,2},{2,4}};
	double _x2[][2]={{1,1},{3,3}};
	IntervalVector x1(2,_x1);
	IntervalVector x2(2,_x2);

	CPPUNIT_ASSERT(x1.is_superset(x2));
	CPPUNIT_ASSERT(x2.is_subset(x1));
	CPPUNIT_ASSERT(x1.is_strict_superset(x2));
	CPPUNIT_ASSERT(x2.is_strict_interior_subset(x1));
}

void TestIntervalVector::is_subset05() {
	double _x1[][2]={{0,2},{2,4}};
	IntervalVector x1(2,_x1);
	IntervalVector x2(IntervalVector::empty(2));

	CPPUNIT_ASSERT(x1.is_superset(x2));
	CPPUNIT_ASSERT(x2.is_subset(x1));
	CPPUNIT_ASSERT(x1.is_strict_superset(x2));
	CPPUNIT_ASSERT(x2.is_strict_interior_subset(x1));
}

void TestIntervalVector::is_subset06() {

	double _x2[][2]={{1,1},{3,3}};

	IntervalVector x1(IntervalVector::empty(2));
	IntervalVector x2(2,_x2);

	CPPUNIT_ASSERT(!x1.is_superset(x2));
	CPPUNIT_ASSERT(!x2.is_subset(x1));
	CPPUNIT_ASSERT(!x1.is_strict_superset(x2));
	CPPUNIT_ASSERT(!x2.is_strict_interior_subset(x1));
}

void TestIntervalVector::is_subset07() {
	double _x1[][2]={{0,2},{2,4}};
	double _x2[][2]={{1,1},{3,5}};

	IntervalVector x1(2,_x1);
	IntervalVector x2(2,_x2);

	CPPUNIT_ASSERT(!x1.is_superset(x2));
	CPPUNIT_ASSERT(!x2.is_subset(x1));
	CPPUNIT_ASSERT(!x1.is_strict_superset(x2));
	CPPUNIT_ASSERT(!x2.is_strict_interior_subset(x1));
}

void TestIntervalVector::extr_diam_index01() {
	 double _x1[][2]={{0,2},{0,1},{0,3}};
	IntervalVector x1(3,_x1);
	CPPUNIT_ASSERT(x1.extr_diam_index(true)==1);
	CPPUNIT_ASSERT(x1.extr_diam_index(false)==2);
	CPPUNIT_ASSERT(x1.min_diam()==1);
	CPPUNIT_ASSERT(x1.max_diam()==3);
}

void TestIntervalVector::extr_diam_index02() {
	double _x1[][2]={{0,1},{0,3},{0,2}};
	IntervalVector x1(3,_x1);
	CPPUNIT_ASSERT(x1.extr_diam_index(true)==0);
	CPPUNIT_ASSERT(x1.extr_diam_index(false)==1);
	CPPUNIT_ASSERT(x1.min_diam()==1);
	CPPUNIT_ASSERT(x1.max_diam()==3);
}

void TestIntervalVector::extr_diam_index03() {
	double _x1[][2]={{0,3},{0,2},{0,1}};
	IntervalVector x1(3,_x1);
	CPPUNIT_ASSERT(x1.extr_diam_index(true)==2);
	CPPUNIT_ASSERT(x1.extr_diam_index(false)==0);
	CPPUNIT_ASSERT(x1.min_diam()==1);
	CPPUNIT_ASSERT(x1.max_diam()==3);
}

void TestIntervalVector::extr_diam_index04() {
	double _x1[][2]={{0,1},{0,2},{NEG_INFINITY,0}};
	IntervalVector x1(3,_x1);
	CPPUNIT_ASSERT(x1.extr_diam_index(true)==0);
	CPPUNIT_ASSERT(x1.extr_diam_index(false)==2);
	CPPUNIT_ASSERT(x1.min_diam()==1);
	CPPUNIT_ASSERT(x1.max_diam()==POS_INFINITY);
}

void TestIntervalVector::extr_diam_index05() {
	double _x1[][2]={{NEG_INFINITY,0}};
	IntervalVector x1(1,_x1);
	CPPUNIT_ASSERT(x1.extr_diam_index(true)==0);
	CPPUNIT_ASSERT(x1.extr_diam_index(false)==0);
	CPPUNIT_ASSERT(x1.min_diam()==POS_INFINITY);
	CPPUNIT_ASSERT(x1.max_diam()==POS_INFINITY);
}

void TestIntervalVector::extr_diam_index06() {
	double _x1[][2]={{NEG_INFINITY,0},{0,1},{NEG_INFINITY,1},{1,3}};
	IntervalVector x1(4,_x1);
	CPPUNIT_ASSERT(x1.extr_diam_index(true)==1);
	CPPUNIT_ASSERT(x1.extr_diam_index(false)==2);
	CPPUNIT_ASSERT(x1.min_diam()==1);
	CPPUNIT_ASSERT(x1.max_diam()==POS_INFINITY);
}

void TestIntervalVector::extr_diam_index07() {
	double _x1[][2]={{NEG_INFINITY,0},{-2,POS_INFINITY},{NEG_INFINITY,1}};
	IntervalVector x1(3,_x1);
	CPPUNIT_ASSERT(x1.extr_diam_index(true)==0);
	CPPUNIT_ASSERT(x1.extr_diam_index(false)==1);
	CPPUNIT_ASSERT(x1.min_diam()==POS_INFINITY);
	CPPUNIT_ASSERT(x1.max_diam()==POS_INFINITY);
}

void TestIntervalVector::extr_diam_index08() {
	double _x1[][2]={{NEG_INFINITY,0},{NEG_INFINITY,1},{-2,POS_INFINITY}};
	IntervalVector x1(3,_x1);
	CPPUNIT_ASSERT(x1.extr_diam_index(true)==0);
	CPPUNIT_ASSERT(x1.extr_diam_index(false)==2);
	CPPUNIT_ASSERT(x1.min_diam()==POS_INFINITY);
	CPPUNIT_ASSERT(x1.max_diam()==POS_INFINITY);
}

void TestIntervalVector::extr_diam_index09() {
	double _x1[][2]={{-2,POS_INFINITY},{NEG_INFINITY,0},{NEG_INFINITY,1}};
	IntervalVector x1(3,_x1);
	CPPUNIT_ASSERT(x1.extr_diam_index(true)==1);
	CPPUNIT_ASSERT(x1.extr_diam_index(false)==0);
	CPPUNIT_ASSERT(x1.min_diam()==POS_INFINITY);
	CPPUNIT_ASSERT(x1.max_diam()==POS_INFINITY);
}

void TestIntervalVector::extr_diam_index10() {
	double _x1[][2]={{-2,POS_INFINITY},{NEG_INFINITY,1},{NEG_INFINITY,0}};
	IntervalVector x1(3,_x1);
	CPPUNIT_ASSERT(x1.extr_diam_index(true)==2);
	CPPUNIT_ASSERT(x1.extr_diam_index(false)==0);
	CPPUNIT_ASSERT(x1.min_diam()==POS_INFINITY);
	CPPUNIT_ASSERT(x1.max_diam()==POS_INFINITY);
}

void TestIntervalVector::sort_indices01() {
	double _x[][2]={{0,2},{NEG_INFINITY,0},{0,1},{3,3},{-10,10}};
	IntervalVector x(5,_x);
	int tab[5];
	x.sort_indices(true,tab);
	CPPUNIT_ASSERT(tab[0]==3);
	CPPUNIT_ASSERT(tab[1]==2);
	CPPUNIT_ASSERT(tab[2]==0);
	CPPUNIT_ASSERT(tab[3]==4);
	CPPUNIT_ASSERT(tab[4]==1);
}

void TestIntervalVector::sort_indices02() {
	double _x[][2]={{0,2},{NEG_INFINITY,0},{0,1},{3,3},{-10,10}};
	IntervalVector x(5,_x);
	int tab[5];
	x.sort_indices(false,tab);
	CPPUNIT_ASSERT(tab[0]==1);
	CPPUNIT_ASSERT(tab[1]==4);
	CPPUNIT_ASSERT(tab[2]==0);
	CPPUNIT_ASSERT(tab[3]==2);
	CPPUNIT_ASSERT(tab[4]==3);
}


void TestIntervalVector::rel_distance01() {
  IntervalVector box1(3);
  IntervalVector box2(3);
  box1[0]=Interval(0,0);
  box2[0]=Interval(0,0);
  box1[1]=Interval(-1,0);
  box2[1]=Interval(-1,0);
  box1[2]=Interval(1,4);
  box2[2]=Interval(1.5,3);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/3.0,box1.rel_distance(box2),ERROR);
}

void TestIntervalVector::perimeter01() {
  IntervalVector box1(3);
  IntervalVector box2(3);
  box1[0]=Interval(0,0);
  box2[0]=Interval(0,0);
  box1[1]=Interval(-1,0);
  box2[1]=Interval(-1,0);
  box1[2]=Interval(1,4);
  box2[2]=Interval(1.5,3);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0,box1.perimeter(),ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.5,box2.perimeter(),ERROR);
}

void TestIntervalVector::perimeter02() {
	double _x1[][2]={{0,1},{0,POS_INFINITY}};
	CPPUNIT_ASSERT(IntervalVector(2,_x1).perimeter()==POS_INFINITY);
}

void TestIntervalVector::volume01() {
	double _x1[][2]={{0,1},{0,POS_INFINITY}};
	CPPUNIT_ASSERT(IntervalVector(2,_x1).volume()==POS_INFINITY);
}

void TestIntervalVector::volume02() {
	double _x1[][2]={{0,1},{1,1}};
	CPPUNIT_ASSERT(IntervalVector(2,_x1).volume()==0);
}

void TestIntervalVector::volume03() {
	double _x1[][2]={{0,2},{2,5},{4,8}};
	CPPUNIT_ASSERT_DOUBLES_EQUAL(24.0,IntervalVector(3,_x1).volume(),ERROR);
}

void TestIntervalVector::minus01() {
	double _x1[][2]={{0,3},{0,2},{0,1}};
	double _x2[][2]={{-3,0},{-2,0},{-1,0}};
	CPPUNIT_ASSERT((-IntervalVector(3,_x1))==IntervalVector(3,_x2));
}

void TestIntervalVector::minus02() {
	double _x1[][2]={{0,1},{0,POS_INFINITY}};
	double _x2[][2]={{-1,0},{NEG_INFINITY,0}};
	check(-IntervalVector(2,_x1),IntervalVector(2,_x2));
}

void TestIntervalVector::minus03() {
	check(-IntervalVector::empty(2),IntervalVector::empty(2));
}

void TestIntervalVector::add01() {
	double _x1[][2]={{0,3},{0,2},{0,1}};
	double _x2[][2]={{0,1},{0,1},{0,1}};
	double _x3[][2]={{0,4},{0,3},{0,2}};

	IntervalVector x1(3,_x1);
	IntervalVector x2(3,_x2);
	IntervalVector x3(3,_x3);
	IntervalVector e(IntervalVector::empty(3));

	check(x1+x2,x3);
	//check(x1+e,e);
	CPPUNIT_ASSERT((x1+e).is_empty());
	//check(IntervalVector(x1)+=e,e);
	CPPUNIT_ASSERT((IntervalVector(x1)+=e).is_empty());

	//check(e+x1,e);
	CPPUNIT_ASSERT((e+x1).is_empty());
	//check(e+=x1,e);
	CPPUNIT_ASSERT((e+=x1).is_empty());
	//check(e+e,e);
	CPPUNIT_ASSERT((e+e).is_empty());
	//check(e+=e,e);
	CPPUNIT_ASSERT((e+=e).is_empty());

	check(IntervalVector(x1)+=x2,x3);
	//check(IntervalVector(x1)+=e,e);
	CPPUNIT_ASSERT((IntervalVector(x1)+=e).is_empty());

	check(IntervalVector(x2)+=x1,x3);
}

void TestIntervalVector::sub01() {
	double _x1[][2]={{0,3},{0,2},{0,1}};
	double _x2[][2]={{0,1},{0,1},{0,1}};
	double _x3[][2]={{-1,3},{-1,2},{-1,1}};
	IntervalVector x1(3,_x1);
	IntervalVector x2(3,_x2);
	IntervalVector x3(3,_x3);
	IntervalVector e(IntervalVector::empty(3));

	check(x1-x2,x3);
	check(x2-x1,-x3);
	//check(x1-e,e);
	CPPUNIT_ASSERT((x1-e).is_empty());
	//check(IntervalVector(x1)-=e,e);
	CPPUNIT_ASSERT((IntervalVector(x1)-=e).is_empty());

	//check(e-x1,e);
	CPPUNIT_ASSERT((e-x1).is_empty());
	//check(e-=x1,e);
	CPPUNIT_ASSERT((e-=x1).is_empty());
	//check(e-e,e);
	CPPUNIT_ASSERT((e-e).is_empty());
	//check(e-=e,e);
	CPPUNIT_ASSERT((e-=e).is_empty());

	check(IntervalVector(x1)-=x2,x3);
	check(IntervalVector(x2)-=x1,-x3);
}

void TestIntervalVector::compl01() {
	double _b[][2]={{0,1},{0,1}};
	IntervalVector b(2,_b);
	IntervalVector* c;
	int n=b.complementary(c);

	CPPUNIT_ASSERT(n==4);
	CPPUNIT_ASSERT(c[0].size()==2);

	CPPUNIT_ASSERT(c[0][0]==Interval::NEG_REALS);
	CPPUNIT_ASSERT(c[0][1]==Interval::ALL_REALS);

	CPPUNIT_ASSERT(c[1][0]==Interval(1,POS_INFINITY));
	CPPUNIT_ASSERT(c[1][1]==Interval::ALL_REALS);

	CPPUNIT_ASSERT(c[2][0]==Interval(0,1));
	CPPUNIT_ASSERT(c[2][1]==Interval::NEG_REALS);

	CPPUNIT_ASSERT(c[3][0]==Interval(0,1));
	CPPUNIT_ASSERT(c[3][1]==Interval(1,POS_INFINITY));

	delete[] c;
}

/**
 * complementary of an empty box = (-oo,oo)x...(-oo,oo)
 */
void TestIntervalVector::compl02() {

	IntervalVector* c;
	int n=IntervalVector::empty(2).complementary(c);
	CPPUNIT_ASSERT(n==1);
	CPPUNIT_ASSERT(c[0].size()==2);

	CPPUNIT_ASSERT(c[0][0]==Interval::ALL_REALS);
	CPPUNIT_ASSERT(c[0][1]==Interval::ALL_REALS);

	delete[] c;
}

/** b\b = emptyset */
void TestIntervalVector::diff01() {
	double _b[][2]={{0,1},{0,1}};
	IntervalVector b(2,_b);
	IntervalVector* c;
	int n=b.diff(b,c);

	CPPUNIT_ASSERT(n==0);

	CPPUNIT_ASSERT(c[0].size()==2);
	CPPUNIT_ASSERT(c[0].is_empty());

	delete[] c;
}

/** b\emptyset = b */
void TestIntervalVector::diff02() {
	double _b[][2]={{0,1},{0,1}};
	IntervalVector b(2,_b);
	IntervalVector* c;

	int n=b.diff(IntervalVector::empty(2),c);

	CPPUNIT_ASSERT(n==1);
	CPPUNIT_ASSERT(c[0].size()==2);
	CPPUNIT_ASSERT(c[0]==b);

	delete[] c;
}

/**
 * [-7,7]x[-7,7]  \ [-5,5]x[-5,5] =
 *  {  [-7,-5[x[-7,7] ; ]5,7]x[-7,7]  ; [-5,5]x[-7,-5[ ; [-5,5]x]5,7] }
 */
void TestIntervalVector::diff03() {
	double _b1[][2]={{-7,7},{-7,7}};
	double _b2[][2]={{-5,5},{-5,5}};
	IntervalVector b1(2,_b1);
	IntervalVector b2(2,_b2);
	IntervalVector* c;

	int n=b1.diff(b2,c);

	CPPUNIT_ASSERT (n==4);

	double _b3[][2]={{-7,-5}, {-7, 7}};
	double _b4[][2]={{ 5, 7}, {-7, 7}};
	double _b5[][2]={{-5, 5}, {-7,-5}};
	double _b6[][2]={{-5, 5}, { 5, 7}};
	IntervalVector b3(2,_b3);
	IntervalVector b4(2,_b4);
	IntervalVector b5(2,_b5);
	IntervalVector b6(2,_b6);

	CPPUNIT_ASSERT(c[0]==b3);
	CPPUNIT_ASSERT(c[1]==b4);
	CPPUNIT_ASSERT(c[2]==b5);
	CPPUNIT_ASSERT(c[3]==b6);

	delete[] c;
}

/**
 * [-7,7]x[-7,7]  \ [-7,7]x[-7,5] =
 *  {  [-7,7]x]5,7] }
 */
void TestIntervalVector::diff04() {
	double _b1[][2]={{-7,7},{-7,7}};
	double _b2[][2]={{-7,7},{-7,5}};
	IntervalVector b1(2,_b1);
	IntervalVector b2(2,_b2);
	IntervalVector* c;

	int n=b1.diff(b2,c);

	CPPUNIT_ASSERT (n==1);

	double _b3[][2]={{-7,7}, {5, 7}};
	IntervalVector b3(2,_b3);

	CPPUNIT_ASSERT(c[0]==b3);

	delete[] c;

}

void TestIntervalVector::diff05() {
	double _b[][2]={{0,1},{0,1}};
	IntervalVector b(2,_b);
	IntervalVector* c;

	int n=IntervalVector::empty(2).diff(b,c);

	CPPUNIT_ASSERT(n==0);

	delete[] c;

}

void TestIntervalVector::random01() {
	double _b[][2]={{0,1},{0,1}};
	IntervalVector b(2,_b);
	IntervalVector r=b.random();

	CPPUNIT_ASSERT(r.size()==2);
	CPPUNIT_ASSERT(r.max_diam()==0);
	CPPUNIT_ASSERT(b.is_superset(r));
}

void TestIntervalVector::random02() {
	double _b[][2]={{1,1},{2,2}};
	IntervalVector b(2,_b);
	IntervalVector r=b.random();

	CPPUNIT_ASSERT(b==r);
}
