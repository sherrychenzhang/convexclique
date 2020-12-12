#include "Structure.h"

#include<time.h>

using namespace std;



int main(int argc, char *argv[]) {


	double dist_r=1;

// string path = "./data/example";
 string path = "./data/Gowalla_1166";
//string path = "./data/Gowalla";
	//string path = "./data/tst.txt";
	//string path = "./data/aaa";
//	string path = "./data/OLout.txt";
//string path = "./data/CRout.txt";
	//string path = "./data/TGout.txt";
	//string path = "./data/OLout.txt";
	printf(" %lf\n", dist_r);

	Graph* g=new Graph();

//	/Tree *t= new Tree();

	clock_t start;
	double duration;


	int convexTime=0;

    start = clock();
	g->readGraph(path,dist_r);



//
// g->tomitaAlgo();

////
//   cout<<"_________________"<<endl;
  g->tomitaConvex();
//////

//g->tomitaWithIncrementalConvex();


// g->test();

    duration = ( clock() - start ) / (double)1000000;
	printf("runing time : %f (s)\n", duration);

}
