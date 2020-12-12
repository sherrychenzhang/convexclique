#ifndef STRUCTURE_H_
#define STRUCTURE_H_

#include<map>
#include<cstring>
#include<string>
#include<vector>
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<algorithm>
#include<set>
#include<math.h>
#include "toplevel_tree.h"
#include "cqueue.h"
#include <queue>
#include<stack>
#include<time.h>

using namespace std;


//vector<point> VCoord_backup;
#define eps 1e-10

class Graph{
public:
	//vector<point> VCoord;
	vector<vector<int>*> con;
	bool **adjacencyMatrix;
	int n;
	int m;

	vector<vector<int> > intVV;
	bool *aa;

	int *vertexSets;
	int* vertexLookup;

	long cliqueCount;

	Tree tree1;

	vector<int> partialClique;
	vector<int> pAdd;

	vector<vector<point> > dynamicConvex;


	double duration1;


	int CountVerAdd;





public:
	Graph();

public:
	void readGraph(string input_path,double dist_r);

	int tomitaAlgo();
	void listAllMaximalCliquesMatrixRecursive(int beginX, int beginP, int beginR);
	int findBestPivotNonNeighborsMatrix( int** pivotNonNeighbors, int* numNonNeighbors, int beginX, int beginP, int beginR);
	void moveToRMatrix(int vertex, int* pBeginX, int *pBeginP, int *pBeginR, int* pNewBeginX, int* pNewBeginP, int *pNewBeginR);
	void moveFromRToXMatrix(int vertex,int* pBeginX, int *pBeginP, int *pBeginR);

	int tomitaConvex();
	void CliqueWithConvex(int beginX, int beginP, int beginR);
	void addQueuePointsToSetLeft(vector<point> &v , node *current_node);
	void addQueuePointsToSetRight(vector<point> &v , node *current_node);

	double xmulti(point p1, point p2, point p0);
	bool isInConvex(vector<point> &vec, point p0);


	void CliqueWithIncrementalConvex(int beginX, int beginP, int beginR);
	void addPoint(vector<point> &a, point p);
	int sqDist(point p1, point p2);
	int orientation(point a,point b,point c);
	int tomitaWithIncrementalConvex();

    bool pointCmp( point &a, point &b, point &center);
	void clockwiseSortPoints(vector<point> &vPoints);

	double cross(const point &O, const point &A, const point &B);
	vector<point> convex_hull(vector<point> P);


	void test();
	bool compare(int a,int b);



	void Merge(int *A,int *L,int leftCount,int *R,int rightCount,point center);
	void MergeSort(int *A,int n, point center);

};



#endif /* STRUCTURE_H_ */
