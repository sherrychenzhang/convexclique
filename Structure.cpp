#include "Structure.h"
#include <iostream>
#include <fstream>

using namespace std;

Graph::Graph() {
	n = 0;
	m = 0;
	tree1=Tree();
}

//vector<point> VCoord_backup;

	vector<point> VCoord;

void Graph::readGraph(string input_path, double dist_r) {

	con.clear();

	FILE *fin = fopen(input_path.c_str(), "r");
	fscanf(fin, "%d", &n);
	printf("number of vertices : %d\n", n);

	adjacencyMatrix = new bool*[n];

	for (int i = 0; i < n; i++) {
		adjacencyMatrix[i] = new bool[n];
		memset(adjacencyMatrix[i], false, n * sizeof(bool));
	}

	for (int i = 0; i < n; i++) {
		double x, y;
		fscanf(fin, "%lf", &x);
		fscanf(fin, "%lf", &y);
		con.push_back(new vector<int>());
		point p = point();
		p.x_coord = x;
		p.y_coord = y;
		VCoord.push_back(p);

	}
	double dist_square = dist_r * dist_r;

	printf("dist_square:%lf\n", dist_square);

	double max=0;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {

//			if (VCoord[j].x_coord - VCoord[i].x_coord > 2 * dist_r) {
//				break;
//			}

			double a1 = (VCoord[i].x_coord - VCoord[j].x_coord)
					* (VCoord[i].x_coord - VCoord[j].x_coord);
			double a2 = (VCoord[i].y_coord - VCoord[j].y_coord)
					* (VCoord[i].y_coord - VCoord[j].y_coord);

			double dist = a1 + a2;

			if (dist <= dist_square) {
				con[i]->push_back(j);
				con[j]->push_back(i);
				adjacencyMatrix[i][j] = true;
				adjacencyMatrix[j][i] = true;
				m++;
			}



			if(max<dist){
				max=dist;
			}

		}

	}


	cout<<"size::"<<VCoord.size()<<endl;



	cout<<"graph distance:"<<sqrt(max)<<endl;




	m = 2 * m;


//        for(int i=0;i<n;i++){
//        	for(int j=0;j<con[i]->size();j++){
//        		printf("%d\t", (*con[i])[j]);
//        	}
//        	printf("\n");
//        }
//
//        for(int i=0;i<n;i++){
//        	for(int j=0;j<n;j++){
//        		printf("%d\t",adjacencyMatrix[i][j] );
//        	}
//        	printf("\n");
//        }

//        for(int i=0;i<VCoord.size();i++){
//        	printf("%lf\t",VCoord[i].x_coord);
//        	printf("%lf\n",VCoord[i].y_coord);
//        }

	fclose(fin);

}

int Graph::tomitaAlgo() {

	partialClique.clear();
	// vertex sets are stored in an array like this:
	// |--X--|--P--|
	vertexSets = new int[n];
	// vertex i is stored in vertexSets[vertexLookup[i]]
	vertexLookup = new int[n];

	partialClique.clear();

	for (int i = 0; i < n; i++) {
		vertexSets[i] = i;
		vertexLookup[i] = i;
	}

	int beginX = 0;
	int beginP = 0;
	int beginR = n;
	cliqueCount = 0;

	listAllMaximalCliquesMatrixRecursive(beginX, beginP, beginR);

	printf("clique Number::%d\n", cliqueCount);


	return cliqueCount;
}

void Graph::listAllMaximalCliquesMatrixRecursive(int beginX, int beginP,
		int beginR) {


	if (beginX >= beginP && beginP >= beginR) {
		cliqueCount++;

//		vector<int> test(partialClique);
//			std::sort(test.begin(), test.end());
//
//		    intVV.push_back(test);


//		if(partialClique[0]==261){
//			for(int i=0;i<partialClique.size();i++){
//
//
//				cout<<partialClique[i]<<"\t";
//
//				//cout<<VCoord[partialClique[i]].x_coord<<"\t"<<VCoord[partialClique[i]].y_coord<<endl;
//			}

//
//		}



//		vector<point> convex;
//		for(int i=0;i<partialClique.size();i++){
//			convex.push_back(VCoord[partialClique[i]]);
//
//
//		}
//
//		convex=convex_hull(convex);



		return;
	}

	if (beginP >= beginR)
		return;

	int* myCandidatesToIterateThrough;
	int numCandidatesToIterateThrough;

	findBestPivotNonNeighborsMatrix(&myCandidatesToIterateThrough,
			&numCandidatesToIterateThrough, beginX, beginP, beginR);

	if (numCandidatesToIterateThrough != 0) {
		int iterator = 0;
		while (iterator < numCandidatesToIterateThrough) {
			int vertex = myCandidatesToIterateThrough[iterator];

			int newBeginX, newBeginP, newBeginR;

			partialClique.push_back(vertex);

			moveToRMatrix(vertex, &beginX, &beginP, &beginR, &newBeginX,
					&newBeginP, &newBeginR);

			listAllMaximalCliquesMatrixRecursive(newBeginX, newBeginP,
					newBeginR);

			partialClique.pop_back();

			moveFromRToXMatrix(vertex, &beginX, &beginP, &beginR);

			iterator++;

		}

		iterator = 0;
		while (iterator < numCandidatesToIterateThrough) {
			int vertex = myCandidatesToIterateThrough[iterator];
			int vertexLocation = vertexLookup[vertex];

			beginP--;
			vertexSets[vertexLocation] = vertexSets[beginP];
			vertexSets[beginP] = vertex;
			vertexLookup[vertex] = beginP;
			vertexLookup[vertexSets[vertexLocation]] = vertexLocation;

			iterator++;
		}

			delete(myCandidatesToIterateThrough);

	}

}

int Graph::findBestPivotNonNeighborsMatrix(int** pivotNonNeighbors,
		int* numNonNeighbors, int beginX, int beginP, int beginR) {
	int pivot = -1;
	int maxIntersectionSize = -1;

	int i = beginX;

	// loop through all vertices in P union X, find u
	while (i < beginR) {
		int vertex = vertexSets[i];
		int neighborCount = 0;

		// count the number of neighbors vertex has in P.
		int j = beginP;
		while (j < beginR) {
			if (adjacencyMatrix[vertexSets[j]][vertex])
				neighborCount++;
			j++;
		}

		// if vertex has more neighbors in P, then update the pivot
		if (neighborCount > maxIntersectionSize) {
			maxIntersectionSize = neighborCount;
			pivot = vertex;
		}
		i++;

	}

	*numNonNeighbors = 0;

	int tmp = beginR - beginP - maxIntersectionSize;
	if (tmp > 0) {
		*pivotNonNeighbors = new int[tmp];

		int j = beginP;

		while (j < beginR) {
			if (!adjacencyMatrix[pivot][vertexSets[j]]) {
				(*pivotNonNeighbors)[*numNonNeighbors] = vertexSets[j];
				(*numNonNeighbors)++;
			}

			j++;
		}
	}

	return pivot;
}

void Graph::moveToRMatrix(int vertex, int* pBeginX, int *pBeginP, int *pBeginR,
		int* pNewBeginX, int* pNewBeginP, int *pNewBeginR) {

	int vertexLocation = vertexLookup[vertex];

	// swap vertex into R and update beginR
	(*pBeginR)--;
	vertexSets[vertexLocation] = vertexSets[*pBeginR];
	vertexLookup[vertexSets[*pBeginR]] = vertexLocation;
	vertexSets[*pBeginR] = vertex;
	vertexLookup[vertex] = *pBeginR;

	*pNewBeginX = *pBeginP;
	*pNewBeginP = *pBeginP;
	*pNewBeginR = *pBeginP;

	// for each vertex in X, ask if it has vertex as a neighbor,
	// if it does, keep it in X. Otherwise, swap it out.
	int j = *pBeginX;
	while (j < *pNewBeginX) {
		int neighbor = vertexSets[j];
		int neighborLocation = j;

		if (adjacencyMatrix[vertex][neighbor]) {
			// swap into new X territory
			(*pNewBeginX)--;
			vertexSets[neighborLocation] = vertexSets[*pNewBeginX];
			vertexLookup[vertexSets[*pNewBeginX]] = neighborLocation;
			vertexSets[*pNewBeginX] = neighbor;
			vertexLookup[neighbor] = *pNewBeginX;
		} else {
			j++;
		}

	}

	// for each vertex in P, ask if it has vertex as a neighbor,
	// if it does, keep it in P. Otherwise, swap it out.
	j = *pBeginP;
	while (j < *pBeginR) {
		int neighbor = vertexSets[j];
		int neighborLocation = j;

		if (adjacencyMatrix[vertex][neighbor]) {
			// swap into new P territory
			vertexSets[neighborLocation] = vertexSets[*pNewBeginR];
			vertexLookup[vertexSets[*pNewBeginR]] = neighborLocation;
			vertexSets[*pNewBeginR] = neighbor;
			vertexLookup[neighbor] = *pNewBeginR;

			(*pNewBeginR)++;
		}

		j++;
	}

}

void Graph::moveFromRToXMatrix(int vertex, int* pBeginX, int *pBeginP,
		int *pBeginR) {

	int vertexLocation = vertexLookup[vertex];

	//swap vertex into X and increment beginP and beginR
	vertexSets[vertexLocation] = vertexSets[*pBeginP];
	vertexLookup[vertexSets[*pBeginP]] = vertexLocation;
	vertexSets[*pBeginP] = vertex;
	vertexLookup[vertex] = *pBeginP;

	*pBeginP = *pBeginP + 1;
	*pBeginR = *pBeginR + 1;

}

void Graph::addQueuePointsToSetLeft(vector<point> &v, node *current_node) {
	if (current_node != NULL) {
		point p(current_node->x_coord, current_node->data2);
		if (current_node->type == LEAF_NODE
				&& find(v.begin(), v.end(), p) == v.end()) {

			v.push_back(p);
		}

		addQueuePointsToSetLeft(v, current_node->left);

		addQueuePointsToSetLeft(v, current_node->middle);

		addQueuePointsToSetLeft(v, current_node->right);

	}

}

void Graph::addQueuePointsToSetRight(vector<point> &v, node *current_node) {
	if (current_node != NULL) {
		point p(current_node->x_coord, current_node->data2);
		if (current_node->type == LEAF_NODE
				&& find(v.begin(), v.end(), p) == v.end()) {

			v.push_back(p);
		}

		addQueuePointsToSetRight(v, current_node->right);

		addQueuePointsToSetRight(v, current_node->middle);

		addQueuePointsToSetRight(v, current_node->left);

	}

}


int Graph::tomitaConvex() {

   partialClique.clear();
	// vertex sets are stored in an array like this:
	// |--X--|--P--|
	vertexSets = new int[n];
	// vertex i is stored in vertexSets[vertexLookup[i]]
	vertexLookup = new int[n];

	for (int i = 0; i < n; i++) {
		vertexSets[i] = i;
		vertexLookup[i] = i;
	}

	int beginX = 0;
	int beginP = 0;
	int beginR = n;
	cliqueCount=0;


	CliqueWithConvex(beginX, beginP, beginR);

	printf("clique Number::%d\n", cliqueCount);

	return cliqueCount;

}

//求p1p0和p2p0的叉积,如果大于0,则p1在p2的顺时针方向
double Graph::xmulti(point p1, point p2, point p0) {



	return (p1.x_coord - p0.x_coord) * (p2.y_coord - p0.y_coord)- (p2.x_coord - p0.x_coord) * (p1.y_coord - p0.y_coord);
}


double Graph::cross(const point &O, const point &A, const point &B) {
	return (A.x_coord - O.x_coord) * (B.y_coord - O.y_coord)- (A.y_coord - O.y_coord) * (B.x_coord - O.x_coord);
}

//vec follow clockwise
bool Graph::isInConvex(vector<point> &vec, point p0) {

	int sz = vec.size();

	if (xmulti(p0, vec[1], vec[0]) == 0
			&& (min(vec[1].x_coord, vec[0].x_coord) <= p0.x_coord)
			&& (p0.x_coord <= max(vec[1].x_coord, vec[0].x_coord))
			&& (min(vec[1].y_coord, vec[0].y_coord) <= p0.y_coord)
			&& (p0.y_coord <= max(vec[1].y_coord, vec[0].y_coord)))
		return true;

	if (xmulti(p0, vec[sz - 1], vec[0]) == 0
			&& min(vec[sz - 1].x_coord, vec[0].x_coord) <= p0.x_coord
			&& p0.x_coord <= max(vec[sz - 1].x_coord, vec[0].x_coord)
			&& min(vec[sz - 1].y_coord, vec[0].y_coord) <= p0.y_coord
			&& p0.y_coord <= max(vec[sz - 1].y_coord, vec[0].y_coord))
		return true;


	if (xmulti(p0, vec[1], vec[0]) < eps|| xmulti(p0, vec[sz - 1], vec[0]) > -eps)
		return false;


	int left = 1, right = sz - 1;
	while (right - left != 1) {
		int mid = (left + right) / 2;
		if (xmulti(p0, vec[mid], vec[0]) > eps)
			left = mid;
		else
			right = mid;
	}


	  if(xmulti( p0, vec[right], vec[left]) ==0 && min(vec[left].x_coord , vec[right].x_coord) <= p0.x_coord && p0.x_coord <= max(vec[left].x_coord , vec[right].x_coord) && min(vec[left].y_coord , vec[right].y_coord) <= p0.y_coord && p0.y_coord <= max(vec[left].y_coord , vec[right].y_coord) )
	     return true;


	if (xmulti(p0, vec[right], vec[left]) < eps)
		return false;

	return true;
}


//return convex hull in clockwise order
vector<point> Graph::convex_hull(vector<point> P) {
	int n = P.size(), k = 0;
	if (n == 1|| n==0)
		return P;
	vector<point> H(2 * n);

	// Sort points lexicographically
	sort(P.begin(), P.end());

	// Build lower hull
	for (int i = 0; i < n; ++i) {
		while (k >= 2 && cross(H[k - 2], H[k - 1], P[i]) <= 0)
			k--;
		H[k++] = P[i];
	}

	// Build upper hull
	for (int i = n - 2, t = k + 1; i >= 0; i--) {
		while (k >= t && cross(H[k - 2], H[k - 1], P[i]) <= 0)
			k--;
		H[k++] = P[i];
	}

	H.resize(k - 1);
	reverse(H.begin(),H.end());

	return H;
}


bool sortby(int a, int b){
	return ( VCoord[a].x_coord > VCoord[b].x_coord);
}



void Graph::CliqueWithConvex(int beginX, int beginP, int beginR) {

	if (beginX >= beginP && beginP >= beginR) {
		cliqueCount++;


//		vector<int> test(partialClique);
//		std::sort(test.begin(), test.end());


//
//		if(find(intVV.begin(),intVV.end(), test)==intVV.end())
//		{
//
//			cout<<"false"<<endl;
//		}


//
//        // cout<<partialClique.size()<<endl;
//		for(int i=0;i<partialClique.size();i++){
//			cout<<partialClique[i]<<" ";
//		}
//
//
//////
//		cout<<"\n"<<endl;

		//cout<<"clique Size:"<<partialClique.size()<<endl;

//		vector<point> convex1;
//		for(int i=0;i<partialClique.size();i++){
//			convex1.push_back(VCoord[partialClique[i]]);
//		}
//
//		convex1=convex_hull(convex1);
//
//
//		cout<<pAdd.size()<<"------------"<<partialClique.size()<< "------------"<<convex1.size()<<endl;
		return;
	}

	if (beginP >= beginR){


		cout<<"beigidajlj"<<endl;
		return;


	}


	int* myCandidatesToIterateThrough;
	int numCandidatesToIterateThrough;

    findBestPivotNonNeighborsMatrix(&myCandidatesToIterateThrough,
			&numCandidatesToIterateThrough, beginX, beginP, beginR);


	if (numCandidatesToIterateThrough != 0) {
		int iterator = 0;
		while (iterator < numCandidatesToIterateThrough) {
			int vertex = myCandidatesToIterateThrough[iterator];
			int newBeginX, newBeginP, newBeginR;


			pAdd.push_back(vertex);

			vector<point> vp;

			for (int i = 0; i < pAdd.size(); i++) {
				vp.push_back(VCoord[pAdd[i]]);
			}

			vector<point> convex = convex_hull(vp);

			bool tag=true;

			for(int j=beginX;j<beginP;j++){
				 int v=vertexSets[j];
				 if (isInConvex(convex, VCoord[vertexSets[j]]))
				 {
					// cout<<"false";
                     tag=false;
					 break;
				 }

			}

			int vertexLocation = vertexLookup[vertex];



			if(tag==false){

				cout<<"false"<<endl;

				int vertexLocation = vertexLookup[vertex];

				vertexSets[vertexLocation] = vertexSets[beginP];
				vertexLookup[vertexSets[beginP]] = vertexLocation;

				vertexSets[beginP] = vertex;
				vertexLookup[vertex] = beginP;

				beginP++;


				pAdd.pop_back();
                iterator++;
                continue;
			}


			partialClique.push_back(vertex);

			moveToRMatrix(vertex, &beginX, &beginP, &beginR, &newBeginX, &newBeginP, &newBeginR);

			int newnewX, newnewP, newnewR;
			newnewX = newBeginX;
			newnewP = newBeginP;
			newnewR = newBeginR;

			int end = newBeginR;

			int addsize = 1;




			if (convex.size() >= 3) {
				for (int j = beginP; j < end; j++) {

					int v = vertexSets[j];

				    if (isInConvex(convex, VCoord[vertexSets[j]])) {


						partialClique.push_back(vertexSets[j]);

						newBeginX = newnewX;
						newBeginP = newnewP;
						newBeginR = newnewR;

						moveToRMatrix(vertexSets[j], &newBeginX, &newBeginP,
								&newBeginR, &newnewX, &newnewP, &newnewR);

						end = newnewR;
						addsize++;
						j--;
					}
				}
			}




			CliqueWithConvex( newnewX, newnewP, newnewR);

			moveFromRToXMatrix(vertex, &beginX, &beginP, &beginR);

			for (int i = 0; i < addsize; i++)
				partialClique.pop_back();

			pAdd.pop_back();

			iterator++;
		}

		iterator = 0;
		while (iterator < numCandidatesToIterateThrough) {
			int vertex = myCandidatesToIterateThrough[iterator];
			int vertexLocation = vertexLookup[vertex];

			beginP--;
			vertexSets[vertexLocation] = vertexSets[beginP];
			vertexSets[beginP] = vertex;
			vertexLookup[vertex] = beginP;
			vertexLookup[vertexSets[vertexLocation]] = vertexLocation;

			iterator++;
		}

	}

}

int Graph::orientation(point a,point b,point c)
{
    int res = (b.y_coord-a.y_coord)*(c.x_coord-b.x_coord) -
              (c.y_coord-b.y_coord)*(b.x_coord-a.x_coord);

    if (res == 0)
        return 0;
    if (res > 0)
        return 1;
    return -1;
}

// Returns the square of distance between two input points
int Graph::sqDist(point p1, point p2)
{
    return (p1.x_coord-p2.x_coord)*(p1.x_coord-p2.x_coord) +
           (p1.y_coord-p2.y_coord)*(p1.y_coord-p2.y_coord);
}

// Adds a point p to given convex hull a[]
void Graph::addPoint(vector<point> &a, point p)
{

    // point having minimum distance from the point p
    int ind = 0;
    int n = a.size();
    for (int i=1; i<n; i++)
        if (sqDist(p, a[i]) < sqDist(p, a[ind]))
            ind = i;

    // Find the upper tangent
    int up = ind;
    while (orientation(p, a[up], a[(up+1)%n])>=0)
        up = (up + 1) % n;

    // Find the lower tangent
    int low = ind;
    while (orientation(p, a[low], a[(n+low-1)%n])<=0)
        low = (n+low - 1) % n;

    // Initialize result
    vector<point>ret;

    // making the final hull by traversing points
    // from up to low of given convex hull.
    int curr = up;
    ret.push_back(a[curr]);
    while (curr != low)
    {
        curr = (curr+1)%n;
        ret.push_back(a[curr]);
    }

    // Modify the original vector
    ret.push_back(p);
    a.clear();
    for (int i=0; i<ret.size(); i++)
        a.push_back(ret[i]);
}



bool Graph::pointCmp( point &a, point &b, point &center)
{
	if (a.x_coord - center.x_coord >= 0 && b.x_coord - center.x_coord < 0)
	        return true;
    if (a.x_coord - center.x_coord < 0 && b.x_coord - center.x_coord >= 0)
	        return false;

	if (a.x_coord - center.x_coord == 0 && b.x_coord - center.x_coord == 0) {
	        if (a.y_coord - center.y_coord >= 0 || b.y_coord - center.y_coord >= 0)
	            return a.y_coord > b.y_coord;
	        return b.y_coord > a.y_coord;
	    }

	    // compute the cross product of vectors (center -> a) x (center -> b)
	    //从A到B的角度
	    double det = (a.x_coord - center.x_coord) * (b.y_coord - center.y_coord) - (b.x_coord - center.x_coord) * (a.y_coord - center.y_coord);


	    if (det < 0)
	        return true;
	    if (det > 0)
	        return false;

	   // points a and b are on the same line from the center
	  //  check which point is closer to the center
	    double d1 = (a.x_coord - center.x_coord) * (a.x_coord - center.x_coord) + (a.y_coord - center.y_coord) * (a.y_coord - center.y_coord);
	    double d2 = (b.x_coord - center.x_coord) * (b.x_coord - center.x_coord) + (b.y_coord - center.y_coord) * (b.y_coord - center.y_coord);
	    return d1 > d2;
 }


void Graph::clockwiseSortPoints(vector<point> &vPoints)
{

    point center;

    double x = 0,y = 0;
    for (int i = 0;i < vPoints.size();i++)

     {
        x += vPoints[i].x_coord;
        y += vPoints[i].y_coord;
     }
     center.x_coord = x/vPoints.size();
     center.y_coord = y/vPoints.size();



     for(int i=0;i<vPoints.size()-1;i++){
    	 for(int j=0;j<vPoints.size()-i-1;j++){
    		 if(!pointCmp(vPoints[j],vPoints[j+1],center)){
    			 point tmp=vPoints[j];
    			 vPoints[j]=vPoints[j+1];
    			 vPoints[j+1]=tmp;
    		 }
    	 }

     }


 }

int Graph::tomitaWithIncrementalConvex() {

    partialClique.clear();
	// vertex sets are stored in an array like this:
	// |--X--|--P--|
	vertexSets = new int[n];
	// vertex i is stored in vertexSets[vertexLookup[i]]
	vertexLookup = new int[n];

	for (int i = 0; i < n; i++) {
		vertexSets[i] = i;
		vertexLookup[i] = i;
	}

	int beginX = 0;
	int beginP = 0;
	int beginR = n;
	cliqueCount=0;

	CountVerAdd=0;


	CliqueWithIncrementalConvex(beginX, beginP, beginR);


	cout<<cliqueCount<<endl;
	return cliqueCount;

}



void Graph::CliqueWithIncrementalConvex(int beginX, int beginP, int beginR) {

	if (beginX >= beginP && beginP >= beginR) {


		//cout << "find a clipque" << endl;

		//cout<<cliqueCount<<endl;

		cliqueCount++;


		//cout<<partialClique.size()<<endl;

		for(int i=0;i<partialClique.size();i++){
			cout<<partialClique[i]<<" ";
		}

		cout<<"\n"<<endl;

//		vector<point> convex1;
//		for(int i=0;i<partialClique.size();i++){
//			convex1.push_back(VCoord[partialClique[i]]);
//		}
//
//		convex1=convex_hull(convex1);
//
//
//		cout<<pAdd.size()<<"------------"<<partialClique.size()<< "------------"<<  convex1.size()<<endl;


		return;
	}

	if (beginP >= beginR)
		return;


	int* myCandidatesToIterateThrough;
	int numCandidatesToIterateThrough;

    findBestPivotNonNeighborsMatrix(&myCandidatesToIterateThrough,
			&numCandidatesToIterateThrough, beginX, beginP, beginR);


	vector<point> convex1;
//
//	if(dynamicConvex.size()>0){
//		convex1= dynamicConvex[dynamicConvex.size()-1];
//	}
//
//
//		    point center;
//		    double x = 0,y = 0;
//		    for (int i = 0;i < convex1.size();i++)
//
//		     {
//		        x += convex1[i].x_coord;
//		        y += convex1[i].y_coord;
//		     }
//		     center.x_coord = x/convex1.size();
//		     center.y_coord = y/convex1.size();
//
//		    MergeSort(myCandidatesToIterateThrough, numCandidatesToIterateThrough, center);



	if (numCandidatesToIterateThrough != 0) {
		int iterator = 0;
		while (iterator < numCandidatesToIterateThrough) {



			int vertex = myCandidatesToIterateThrough[iterator];

			int newBeginX, newBeginP, newBeginR;

			partialClique.push_back(vertex);
			pAdd.push_back(vertex);


			moveToRMatrix(vertex, &beginX, &beginP, &beginR, &newBeginX,&newBeginP, &newBeginR);


		    vector<point> convex;

			if(dynamicConvex.size()>0){
				convex= dynamicConvex[dynamicConvex.size()-1];
			}


			if(convex.size()<4){
				convex.push_back(VCoord[vertex]);
				convex = convex_hull(convex);
			}
			else{
				addPoint(convex, VCoord[vertex]);
				clockwiseSortPoints(convex);
			}

			dynamicConvex.push_back(convex);

			int newnewX, newnewP, newnewR;
			newnewX = newBeginX;
			newnewP = newBeginP;
			newnewR = newBeginR;


			int end = newBeginR;

			int addsize = 1;




			//int countConvexIn=0;

			if (convex.size() >= 3) {
				for (int j = beginP; j < end; j++) {

					int v = vertexSets[j];

				    if (isInConvex(convex, VCoord[vertexSets[j]])) {

				    	//countConvexIn++;
						partialClique.push_back(vertexSets[j]);

						newBeginX = newnewX;
						newBeginP = newnewP;
						newBeginR = newnewR;



						moveToRMatrix(vertexSets[j], &newBeginX, &newBeginP,
								&newBeginR, &newnewX, &newnewP, &newnewR);

						end = newnewR;
						addsize++;
						j--;
					}
				}
			}



//			cout<<"convex size"<<convex.size();
//			cout<<"convex in"<<countConvexIn<<endl;

			CliqueWithIncrementalConvex( newnewX, newnewP, newnewR);

			moveFromRToXMatrix(vertex, &beginX, &beginP, &beginR);

			for (int i = 0; i < addsize; i++)
				partialClique.pop_back();

			pAdd.pop_back();
//
//			cout<<"back tracking"<<endl;
			dynamicConvex.pop_back();
			CountVerAdd--;

			iterator++;
		}

		iterator = 0;
		while (iterator < numCandidatesToIterateThrough) {
			int vertex = myCandidatesToIterateThrough[iterator];
			int vertexLocation = vertexLookup[vertex];

			beginP--;
			vertexSets[vertexLocation] = vertexSets[beginP];
			vertexSets[beginP] = vertex;
			vertexLookup[vertex] = beginP;
			vertexLookup[vertexSets[vertexLocation]] = vertexLocation;

			iterator++;
		}

	}

}

void Graph::Merge(int *A,int *L,int leftCount,int *R,int rightCount,point center){
    int i,j,k;

    i = 0; j = 0; k =0;

    while(i<leftCount && j< rightCount) {

		double da=(VCoord[L[i]].x_coord-center.x_coord)*(VCoord[L[i]].x_coord-center.x_coord)+(VCoord[L[i]].y_coord-center.y_coord)*(VCoord[L[i]].y_coord-center.y_coord);
		double db=(VCoord[R[j]].x_coord-center.x_coord)*(VCoord[R[j]].x_coord-center.x_coord)+(VCoord[R[j]].y_coord-center.y_coord)*(VCoord[R[j]].y_coord-center.y_coord);

      if(da>db)
        	A[k++] = L[i++];
        else A[k++] = R[j++];
    }
    while(i < leftCount) A[k++] = L[i++];
    while(j < rightCount) A[k++] = R[j++];

}



void Graph::MergeSort(int *A,int n, point center){
	int mid,i, *L, *R;
	    if(n < 2) return;

	    mid = n/2;

	    L = new int[mid];
	    R = new int [n - mid];

	    for(i = 0;i<mid;i++) L[i] = A[i];
	    for(i = mid;i<n;i++) R[i-mid] = A[i];

	    MergeSort(L,mid,center);
	    MergeSort(R,n-mid,center);
	    Merge(A,L,mid,R,n-mid,center);

	    delete [] R;
	    delete [] L;

}

void Graph::test(){
    point center;
    center.x_coord=0.0;
    center.y_coord=0.0;



    point a0;
    a0.x_coord=29.22105923;
    a0.y_coord=-94.67852218;


    point a1;
    a1.x_coord=29.2648275;
    a1.y_coord=-94.83292222;


    point a2;
    a2.x_coord=29.27238034;
    a2.y_coord=-94.81815934;

    point a3;
    a3.x_coord=29.27355438;
   a3.y_coord=-94.81559563;

    point a4;
    a4.x_coord=29.2735575;
    a4.y_coord=-94.81557418;

    point a5;
     a5.x_coord=29.27574525;
     a5.y_coord=-94.83065808;


     VCoord.push_back(a0);
     VCoord.push_back(a1);
     VCoord.push_back(a2);
     VCoord.push_back(a3);
     VCoord.push_back(a4);
     VCoord.push_back(a5);
    int *A=new int[6]{0,1,2,3,4,5};

    MergeSort(A,6,center);


    for(int i=0;i<6;i++){
    	cout<<A[i]<<" "<<endl;
    }


}
