//============================================================================
//                                  I B E X                                   
// File        : ibex_Set.cpp
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : 13 juil. 2014
//============================================================================

#include "ibex_SetInterval.h"
#include "ibex_SetLeaf.h"
#include "ibex_SetBisect.h"
#include "ibex_CellHeap.h"

using namespace std;

namespace ibex {

SetInterval::SetInterval(const IntervalVector& bounding_box, double eps, bool inner) : root(new SetLeaf(inner? __IBEX_IN__: __IBEX_UNK__)), eps(eps), bounding_box(bounding_box) {

}

bool SetInterval::is_empty() const {
	return root==NULL;
}

void SetInterval::sync(Separator& sep) {
	try {
		root = root->sync(bounding_box, sep, eps);
	} catch(NoSet& e) {
		delete root;
		root = NULL;
		throw e;
	}
}

void SetInterval::contract(Separator& sep) {
	root->set_in_tmp();
	root = root->inter(bounding_box, sep, eps);
	root->unset_in_tmp();
}

SetInterval& SetInterval::operator&=(const SetInterval& set) {
	root->set_in_tmp();
	root = root->inter(bounding_box, set.root, set.bounding_box, eps);
	root->unset_in_tmp();
	return *this;
}

SetInterval& SetInterval::operator|=(const SetInterval& set) {
	root = root->union_(bounding_box, set.root, set.bounding_box, eps);
	return *this;
}


void SetInterval::visit_leaves(SetNode::leaf_func func) const {
	root->visit_leaves(func, bounding_box);
}

std::ostream& operator<<(std::ostream& os, const SetInterval& set) {
	set.root->print(os,set.bounding_box, 0);
	return os;
}

namespace {


class NodeAndDist : public Backtrackable {
public:
	NodeAndDist() : node(NULL), dist(-1) { }

	NodeAndDist(SetNode* node) : node(node), dist(-1) { }

	/**
	 * calculate the square of the distance to pt
	 * for the box of the current cell (box given in argument)
	 */
	void set_dist(const IntervalVector& box, const Vector pt) {
		assert(box.size()==pt.size());

		Interval d=Interval::ZERO;
		for (int i=0; i<pt.size(); i++) {
			d += sqr(box[i]-pt[i]);
		}
		dist=d.lb();
	}

	virtual std::pair<Backtrackable*,Backtrackable*> down() {
		assert(!node->is_leaf());

		SetBisect& b=*((SetBisect*) node);
		return std::pair<NodeAndDist*,NodeAndDist*>(new NodeAndDist(b.left),
													  new NodeAndDist(b.right));
	}

	SetNode* node;
	double dist;
};

/**
 * Cell heap where the criterion is the distance to "pt"
 */
class CellHeapDist : public CellHeap {
public:

	virtual double cost(const Cell& c) const {
		return c.get<NodeAndDist>().dist;
	}
};
}

double SetInterval::dist(const Vector& pt, bool inside) const {
	CellHeapDist heap;

	Cell* root_cell =new Cell(bounding_box);
	root_cell->add<NodeAndDist>();
	root_cell->get<NodeAndDist>().node = root;
	root_cell->get<NodeAndDist>().set_dist(bounding_box,pt);

	heap.push(root_cell);

	double lb = POS_INFINITY;

	while (!heap.empty()) {
		Cell* c = heap.pop();
		SetNode* node = c->get<NodeAndDist>().node;
		if (node->status==__IBEX_OUT__) {
			double d=root_cell->get<NodeAndDist>().dist;
			if (d<lb) {
				lb=d;
				heap.contract_heap(lb);
			}
		} else if (!node->is_leaf() && possibly_contains_out(node->status)) {
			SetBisect b= *((SetBisect*) node);
			IntervalVector left=b.left_box(c->box);
			IntervalVector right=b.right_box(c->box);

			std::pair<Cell*,Cell*> p=c->bisect(left,right);

			p.first->get<NodeAndDist>().set_dist(left,pt);
			if (p.first->get<NodeAndDist>().dist<=lb) heap.push(p.first);

			p.second->get<NodeAndDist>().set_dist(right,pt);
			if (p.second->get<NodeAndDist>().dist<=lb) heap.push(p.second);
		}
	}
	return lb;
}

SetInterval::~SetInterval() {
	delete root;
}

} // namespace ibex
