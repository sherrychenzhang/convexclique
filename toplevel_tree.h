/* $Id$
 * 
 * Copyright (C) 2007 Jose Alberto Cisneros Perez. All rights reserved.
 *
 * This file may be used under the terms of the GNU General Public License
 * versions 3.0 as published by the Free Software Foundation and
 * appearing in the COPYING file included in the packaging of this file.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 *  WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
*/

#ifndef __TOPLEVEL_TREE_H__
#define __TOPLEVEL_TREE_H__

#include <iostream>
#include <cmath>
#include "cqueue.h"
#include <vector>
#include<math.h>
#include<cmath>
#include<algorithm>
#include<iomanip>
using namespace std;

// forward declaration
class ConcatenableQueue;

struct point
{
	point() : x_coord(0.0), y_coord(0.0)
	{
	}
	point(const  double &x, const double &y) : x_coord(x), y_coord(y)
	{
	}
	point(const point &other) : x_coord(other.x_coord), y_coord(other.y_coord)
	{
	}

	bool operator == (const point &rhs) const
	{
		return (x_coord == rhs.x_coord && y_coord == rhs.y_coord);
	}
	bool operator != (const point &rhs) const
	{
		return !(*this == rhs);
	}
  point& operator = (const point &rhs)
	{
		this->x_coord = rhs.x_coord;
		this->y_coord = rhs.y_coord;
		return *this;
	}

  bool operator <(const point &p) const {
  		return x_coord < p.x_coord || (x_coord == p.x_coord && y_coord < p.y_coord);
  	}
	double x_coord, y_coord;
};

struct tree_node
{
	tree_node() : p(), bridge1_lc(), bridge2_lc(), bridge1_rc(), bridge2_rc(), label(0.0),
			Ql(), Qr(), leftchild(NULL), rightchild(NULL), parent(NULL),
			balance_factor(0)
	{
	}
	point p, bridge1_lc, bridge2_lc, bridge1_rc, bridge2_rc; //only leaf nodes have a point stored in them
	double label;//it is the y-coord of the point stored rightmost leaf of the subtree rooted on the left child
	ConcatenableQueue Ql, Qr;
	tree_node *leftchild, *rightchild, *parent;
	int balance_factor;
};

class Tree
{
public:
	Tree();
	void addNode(const point &new_point);
	void deleteNode(const point &point_to_delete);
	tree_node* searchNode(const point &point_to_search);
	void print();
	tree_node* root() const;
	void set_show(const bool &show);
	tree_node *root_;

	vector<point> convex_hull(vector<point> P);
	double cross(const point &O, const point &A, const point &B);
	void addQueuePointsToSet(vector<point> &v, node *current_node);

	double xmulti(point p1, point p2, point p0);
	bool isInConvex(vector<point> &vec, point p0);

	void test();

private:
	tree_node *createLabel(const double &new_label, tree_node *left_child, tree_node *right_child);
	tree_node *createObject(const point &new_point);
	void clockwiseAnticlockDoubleRotation(tree_node *z_node);
	void anticlockwiseClockDoubleRotation(tree_node *z_node);
	void clockwiseRotation(tree_node *z_node);
	void anticlockwiseRotation(tree_node *z_node);
	void updateBalanceFactor(tree_node *current, bool right);
	void updateBalanceFactorDelete(tree_node*, bool right, bool update_label, bool update_label_done, const double &label);
	void printValues(tree_node *current_node, const int &indent);
	void updateBridgeHull(tree_node *v_node, const int &hull);
	void updateHull(tree_node *current_node);
	void buildChildrensHulls(tree_node *parent_node);
	void getPointsRightHalf(node *v_node, point &p, point &p0, point &p1, int &flag); 
	void getPointsLeftHalf(node *v_node, point &q, point &q0, point &q1, int &flag);
	int rightTurn(const point &p0, const point &p1, const point &p2) const;
	double case9(const point &p, const point &p0, const point &p1,
					const point &q, const point &q0, const point &q1, const int &flag_lh, const int &flag_uh) const ;
	int bridgeCasesLcHull(const point &p, const point &p0, const point &p1,
					const point &q, const point &q0, const point &q1,
					const double &half, const int &flag_lh, const int &flag_uh) const;
	int bridgeCasesRcHull(const point &p, const point &p0, const point &p1,
					const point &q, const point &q0, const point &q1,
					const double &half, const int &flag_lh, const int &flag_uh) const;
	node *newSSLhLeft(node *current, point &q, point &q0, point &q1, bool &onepoint);
	node *newSSUhRight(node *current, point &p, point &p0, point &p1, bool &onepoint);
	node *newSSLhRight(node *current, point &q, point &q0, point &q1, bool &onepoint, bool &out);
	node *newSSUhLeft(node *current, point &p, point &p0, point &p1, bool &onepoint, bool &out);


	bool show_;
};


#endif /*__TOPLEVEL_TREE_H__*/
