#include<iostream>
#include<fstream>
#include<stdio.h>
#include<assert.h>
#include<vector>
#include<time.h>
#include <algorithm>
#include "queue.h"
//#include "nonUniformArray.h"
//#include "nonUniformArraySmaller.h"
using namespace std;

static clock_t timeAlg = 0;

static double zeta; // the large shortest distance of two unreachable nodes

// for DFS
static bool visited[MAX_NODES]; // whether it has been visited before

// for shortest distance computation
static  int distances[MAX_NODES]; // distances to a source node
static int numShortestPathsBefore[MAX_NODES];
static int numShortestPathsAfter[MAX_NODES];

static bool isFoundforGBC[MAX_NODES];
static double betweenCentralities[MAX_NODES];
static double currrentBetweenCentralities[MAX_NODES];

const size_t numSamples = 500000;
//const size_t maxNodesEachSample = MAX_ROWS;
//static int sampleShortestPathsArray[numSamples][maxNodesEachSample];

static int nodesInEachSampleSetS[numSamples];
int* sampleShortestPathsArraySetS[numSamples];
int totalSampledShortestPathsSetS = 0;

static int nodesInEachSampleSetT[numSamples];
int* sampleShortestPathsArraySetT[numSamples];
int totalSampledShortestPathsSetT = 0;

//static Queue  bfsQueue; // 

// for dominator algorithm: see paper `A fast algorithm for finding domnators in a flow graph'
static int reachableSet[MAX_NODES]; // the set of reachable nodes from some source node
static int numReachNodes; // number of Reachable nodes

static nonUniformArray<int> propagationGraph;
static nonUniformArray<int> predecessorList;


const size_t maxSampledSources = 12800;
nonUniformArraySmaller<int>* sourcesSuccessorsAdjList[maxSampledSources];
nonUniformArraySmaller<int>* shortestPathsAdjLists;
int* sourcesNumShortestPathsBefore[maxSampledSources];
int* oneSourceNumShortestPathsBefore;
int* sourcesNumShortestPathsAfter[maxSampledSources];
int* oneSourceNumShortestPathsAfter;
int* sourcesIn2OutIndex[maxSampledSources];

nonUniformArraySmaller<int>* sourcesPredecessorsAdjList[maxSampledSources];
nonUniformArraySmaller<int>* shortestPathsPredecessorAdjLists;
double* sourcesGBC[maxSampledSources];
double* oneSourceGBC;
int* oneSourceDecreasedNumShortestPaths;

static bool visitedBackwards[MAX_NODES];
static  int distancesBackwards[MAX_NODES];
static int numShortestPathsBackward[MAX_NODES]; 

struct PAIR{
	int s;
	int t;
	PAIR(int src, int sink) { s = src; t = sink; }
	PAIR() { s = t = 0; }
};

struct EdgeWithWeight {
	int u;
	int v;
	double weight;
	int direction; // 1: forward, 2: backward
	EdgeWithWeight(int startNode, int endNode, int di, double edgeWeight)
	{
		u = startNode;
		v = endNode; 
		direction = di;
		weight = edgeWeight;
	}
};

// another, fast way to represent the buckets
struct LinkNode{ 
	int node;
	int nextNodeIndex; // the index of next node in the list. -1 indicates that it is the last element
};

static int visitNumber = 0;

static bool isReachable[MAX_NODES];

static nonUniformArray<int> parentSCCNeighbors;
static nonUniformArray<int> sccNeighborGroups;
static nonUniformArray<int> sccNeighborSizes;


static int reachSCCsSet[MAX_NODES];
static int numReachSCCs;

static int globalDegrees[MAX_NODES];

static double sumEachSCC[MAX_NODES];

struct Node;
static vector<Node> sortSNs;

static nonUniformArray<int> oneSCC;
static nonUniformArray<int> SCCs;
static nonUniformArray<int> sccGraph;
static nonUniformArray<int> sccTransposeGraph;
static nonUniformArray<int> sccUnderlyingGraph;


// for relabelling the propagation graph 
static int outIndex2InIndex[MAX_NODES];
static int inIndex2OutIndex[MAX_NODES];

class GroupBetweenesssCentrality;

struct Node
{
	int id;
	double key;
	Node(int a, double b){ id = a; key = b; }
	Node(){}
	bool operator>(Node &e)
	{
		if (key > e.key) return true;
		else					return false;
	}
	bool operator>=(Node &e)
	{
		if (key > e.key) return true;
		else					return false;
	}
	bool operator<(Node &e)
	{
		if (key < e.key) return true;
		else					return false;
	}
	bool operator<=(Node &e)
	{
		if (key < e.key) return true;
		else					return false;
	}
};

inline bool compLess(Node& a, Node& b) // operator <
{
	return a.key < b.key;
}

inline bool comp(Node & a, Node & b) // operator >
{
	return a.key > b.key;
}

//static int Mode = 9973;
//static int Mode = 5499979;
inline double Uniform()
{
	//return (rand() % Mode) / (double)Mode;
	return double(rand()) / RAND_MAX;
}

int Uniform(int a, int b)
{
	return  floor(double(rand()) / RAND_MAX * (b - a) + a);
}

inline void swap(int &u, int &v)
{
	int tmp = u;
	u = v;
	v = tmp;
}

struct Edge{
	int u, v;
	Edge(int a, int b){ u = a; v = b; }
	Edge(){}
};

class Graph{
public:
	friend class GroupBetweenesssCentrality;
	friend class HAM;
public:
	void readGraph(const char *graphFile, bool isDirected);
	void generateEdgeNodeWeights(double minEdgeWeight, double maxEdgeWeight, double minNodeWeight, double maxNodeWeight);
	void generateEdgeNodeWeightsCaseStudy(double minNodeWeight, double maxNodeWeight);
	void printGraph();
	void printGraphTopology();
	int findAllCCsDFS(nonUniformArray<int> &allCCs); // find all connected components
	int findAllCCsBFS(nonUniformArray<int> &allCCs); // find all connected components
	void findAPsAndLowerBounds();
	void cleanGraph(const char *graphFile);
	Graph();
	Graph(const Graph &g);
	Graph & operator=(const Graph &g);
	~Graph();

public:
// for shortest distance
public:
	int numOfEdges();
	int numOfNodes(){ return n; }

protected:
	
	double rankByGroupBCWithSampleShortestPaths(vector<PAIR>& randomPairs, vector<int>& foundNodes, int K, 
		vector<double>& foundTimes);
	double rankByGroupBCWithSampleShortestPathsBiBFS(vector<PAIR>& randomPairs, vector<int>& foundNodes, int K,
		vector<double>& foundTimes);

	double rankByGroupBCWithSampleShortestPathsBiBFSunbiasedKDD23(vector<PAIR>& randomPairs, vector<int>& foundNodes, int K,
		vector<double>& foundTimes);
	double evaluateGBCbySampleShoretstPaths(vector<int>& foundNodes, vector<PAIR>& SamplePairs);

	void generateRandomPermulation(int n, vector<int>& randomizedSeq);

	double findTopKGBCbyKDD16(int K, double epsilon, double errorProb, vector<int>& foundNodes, int& sampledPaths);

	double findTopKGBCbyKDD23(int K, double epsilon, double errorProb, vector<int>& foundNodes, int& sampledPaths);

	double findTopKGBCAdaptiveSampling(int K, double epsilon, double errorProb, vector<int>& foundNodes, int& sampledPaths);

	void generateRandomPairs(int n, int L, vector<PAIR>& randomizedPairs);
protected:
	// single-source shortest distance
	void shortestDistances(int src, bool markSPT);
protected:
	void modifiedDFS(int u, nonUniformArray<int> & allCCs, int numNodesValid);
	inline void DFS(int u, int CCID, int & sizeOfCC);
	void clearGraph();
protected:
	nonUniformArray<int> adjLists;
	nonUniformArray<double> adjWeights;
	vector<double> nodesWeights;
	int n, m; // number of nodes and edges in the graph
	// whether a node is valid in the graph, as it may be identified as a structural hole
	// and removed
	bool *isValid;

public: //for DFS
	void initilizeVariables();
	void clearVariables();
	void allocateVariablesForDFS();
private: // for DFS
	int *numChildren; // number of chidren in the DFS tree
	int *numNodesSubtree; // number of nodes in the subtree rooted at each node in the DFS tree
	int *parent; // the parent of each node in the DFS tree
	bool *isAP; // whether the node is an Articulation point
	int *discoveredTime; // the discovered time of each vertex
	int *lowestTime; // the smallest discovered time of any neighbor of u’s descendants (through a back-edge)
	double *lowerBounds; // the lower bound for each node
	int time; // global time counter

protected:
	vector<int> sizeofEachCC; // the number of nodes in each connected component
	int *ccBelongto; // the cc id that each node belongs to

private:
	bool debug; //
};



int Graph::numOfEdges()
{
	int edges = 0;
	for (int i = 0; i < adjLists.size(); ++i)
		edges += adjLists.size(i);
	return edges / 2;
}




void Graph::generateRandomPermulation(int n, vector<int>& randomizedSeq)
{ // generate a random permulation of n numbers 0, 1, ..., n-1
	assert(n > 0);
	randomizedSeq.clear();
	randomizedSeq.reserve(n);
	int i;
	for (i = 0; i < n; ++i)
		randomizedSeq.push_back(i);

	int tmp;
	int index;
	double p;
	for (i = n - 1; i > 0; --i)
	{
		p = Uniform();
		index = int(floor(p * i)); // a random number between 0 and i
		// swap randomizedSeq[index] and randomizedSeq[i]
		tmp = randomizedSeq[index];
		randomizedSeq[index] = randomizedSeq[i];
		randomizedSeq[i] = tmp;
		//printf("swap %d th node and %d th node\n", index, i);
	}
	/*for (i = 0; i < n; ++i)
	{
		printf("%d\t", randomizedSeq[i]);
		if ((i + 1) % 10 == 0) printf("\n");
	}*/
}


double Graph::rankByGroupBCWithSampleShortestPathsBiBFSunbiasedKDD23(vector<PAIR>& randomPairs, vector<int>& foundNodes, int K, vector<double>& foundTimes)
{ // find the sample with bidirectional BFS
	clock_t  st = clock();

	foundNodes.clear();
	foundNodes.reserve(K);
	foundTimes.clear();
	assert(K > 0 && K <= n);
	int L = randomPairs.size();
	assert(L <= numSamples);
	PAIR onePair;

	int i, j, k;
	int destID;
	bool debug = false;

	int sz;
	int* adj;
	int u, v;

	vector<int> degrees; // degree of each node 
	degrees.reserve(n);
	for (i = 0; i < n; ++i) degrees.push_back(adjLists.size(i));
 
	int** originalGraph = new int* [n];
	for (i = 0; i < n; ++i) originalGraph[i] = new int[degrees[i] + 1]; // some times degrees[i] is zero
	int* originalGraphSize = new int[n];
	for (i = 0; i < n; ++i) originalGraphSize[i] = 0;
	for (i = 0; i < n; ++i)
	{
		sz = adjLists.size(i);
		adj = &adjLists.access(i, 0);
		for (j = 0; j < sz; ++j) {
			v = adj[j];
			originalGraph[i][j] = v;
		}
		originalGraphSize[i] = sz;
	}

	int  **predecessorsAdjLists = new int *[n]; // record predecessors in the shortest paths
	for (i = 0; i < n; ++i)
		predecessorsAdjLists[i] = new int[degrees[i] + 1];
	int* predecessorsAdjListsSize = new int[n];

	int src, dst;
	// the number of nodes that their shortest distances to src have been found
	//int numDistancesFound = 0;
	

	double totalGBC = 0;
	double largestIncreasedBC;
	int largestNodeID;

	//printf("num nodes: %d\n", n);
	int newDistance;
	double p;
	double delta;
	double numEdegesExplord = 0;

	int* sampleShortestPathsArray[numSamples];
	const int maxNodesInPath = 100;
	for (i = 0; i < L; ++i) sampleShortestPathsArray[i] = new int[maxNodesInPath];

	degrees.clear();
	for (i = 0; i < n; ++i) degrees.push_back(0); // degrees in the transpose graph
	for (i = 0; i < n; ++i)
	{
		sz = adjLists.size(i);
		adj = &adjLists.access(i, 0);
		for (j = 0; j < sz; ++j)
			degrees[adj[j]]++;
	}
	int** transposeGraph = new int* [n];
	for (i = 0; i < n; ++i) transposeGraph[i] = new int[degrees[i] + 1];
	int* transposeGraphSize = new int[n];
	for (i = 0; i < n; ++i)transposeGraphSize[i] = 0;

	for (i = 0; i < n; ++i)
	{
		sz = adjLists.size(i);
		adj = &adjLists.access(i, 0);
		for (j = 0; j < sz; ++j) {
			v = adj[j];
			transposeGraph[ v ][ transposeGraphSize[v] ] = i;
			transposeGraphSize[v]++;
		}
	}
	int** predecessorAdjListsBackwards = new int* [n];
	for (i = 0; i < n; ++i) predecessorAdjListsBackwards[i] = new int[degrees[i]+1];
	int* predecessorAdjListsBackwardsSize = new int[n];

	//printf("original graph:\n"); adjLists.printArray();
	//printf("\n transpose graph:\n"); transposeGraph.printArray();

	
	clock_t ft = st;
	int middleNode;
	int tempNode;
	vector<int> randomSeq;
	randomSeq.reserve(100);

	NodeInTwoQuques nodeU, nodeV;
	int stopBFSdistanceForward, stopBFSdistanceBackward;
	vector<EdgeWithWeight> touchEdges;
	double totalWeight;
	double sumWeight;

	for (i = 0; i < L; ++i) // randomly generate L shortest paths
	{
		//if (i == 40) debug = true;
		//else debug = false;
		//debug = true;

		onePair = randomPairs[i];
		src = onePair.s;
		dst = onePair.t;

		/*if (i % 1000 == 0) {
			printf("L: %d, find the %d th sample, time used: %.2lf, delta time: %.2lf\n", L, i + 1, (clock() - st) / 1000.0,
				(clock() - ft) / 1000.0);
			ft = clock();
		}*/
		for (j = 0; j < n; ++j) predecessorsAdjListsSize[j] = 0; 
		for (j = 0; j < n; ++j) predecessorAdjListsBackwardsSize[j] = 0;

		// find the shortest paths from src 
		for (j = 0; j < n; ++j) visited[j] = false;
		visited[src] = true;
		for (j = 0; j < n; ++j) visitedBackwards[j] = false;
		visitedBackwards[dst] = true;
		for (j = 0; j < n; ++j) numShortestPathsBefore[j] = 0;
		numShortestPathsBefore[src] = 1;
		for (j = 0; j < n; ++j) numShortestPathsBackward[j] = 0;
		numShortestPathsBackward[dst] = 1;

		distances[src] = 0;
		distancesBackwards[dst] = 0;
		biQueue_clear(); // clear the queue
		biQueue_push(NodeInTwoQuques(src, 1));  //first queue
		biQueue_push(NodeInTwoQuques(dst, 2));  // second queue
		stopBFSdistanceForward = stopBFSdistanceBackward = n - 1;
		touchEdges.clear();
		if (debug) printf("\npair(%d, %d)******************\n", src, dst);

		while (biQueue_size() > 0) // bfs
		{
			nodeU = biQueue_pop();
			u = nodeU.v;
			if (debug) printf("exploring node %d ", u);

			if (1 == nodeU.qID) // forward search
			{
				if (debug) printf("forward \n");
				if (distances[u] >= stopBFSdistanceForward) continue;

				sz = originalGraphSize[u]; 
				adj = &originalGraph[u][0]; 
				newDistance = distances[u] + 1;

				generateRandomPermulation(sz, randomSeq);
				for (j = 0; j < sz; ++j)
				{
					v = adj[randomSeq[j]]; 
					v = adj[j];
					if (false == visited[v])
					{
						visited[v] = true;
						distances[v] = newDistance;
						biQueue_push(NodeInTwoQuques(v, 1));
						if (debug) printf("node %d is added to queue\n", v);
					}
					if (distances[v] == newDistance) {
						predecessorsAdjLists[v][predecessorsAdjListsSize[v]] = u;
						predecessorsAdjListsSize[v]++;
						numShortestPathsBefore[v] += numShortestPathsBefore[u];
					}

					if (true == visitedBackwards[v])
					{
						stopBFSdistanceForward = newDistance;
						stopBFSdistanceBackward = distancesBackwards[v];
						touchEdges.push_back(EdgeWithWeight(u, v, 1, 0));  // forward edge
						if (debug) printf("at node %d, found a middle node: %d\n", u, v);
					}
				}
				numEdegesExplord += sz;
			}
			else { // backward search
				if (debug) printf("backward \n");
				if (distancesBackwards[u] >= stopBFSdistanceBackward) continue;
				sz = transposeGraphSize[u]; 
				adj = &transposeGraph[u][0];
				newDistance = distancesBackwards[u] + 1;
				generateRandomPermulation(sz, randomSeq);
				for (j = 0; j < sz; ++j)
				{
					v = adj[randomSeq[j]];
					//v = adj[j];
					if (false == visitedBackwards[v])
					{
						visitedBackwards[v] = true;
						distancesBackwards[v] = newDistance;
						biQueue_push(NodeInTwoQuques(v, 2));
						if (debug) printf("node %d is added to queue\n", v);
					}
					if (distancesBackwards[v] == newDistance) {
						predecessorAdjListsBackwards[v][predecessorAdjListsBackwardsSize[v]] = u;
						predecessorAdjListsBackwardsSize[v]++;
						numShortestPathsBackward[v] += numShortestPathsBackward[u];
					}

					if (true == visited[v])
					{
						stopBFSdistanceForward = distances[v];
						stopBFSdistanceBackward = newDistance;
						touchEdges.push_back(EdgeWithWeight(u, v, 2, 0)); // backward edge
						if (debug) printf("at node %d, found a middle node: %d\n", u, v);
					}
				}
				numEdegesExplord += sz;
			}
		}
		destID = onePair.t;
		nodesInEachSampleSetS[i] = 0;
		if (touchEdges.size() == 0) {
			printf("pair: %d, %d, not reachable\n", src, destID);
			sampleShortestPathsArray[i][0] = src;
			nodesInEachSampleSetS[i] = 1;
			continue;
		}

		totalWeight = 0;
		for (j = 0; j < touchEdges.size(); ++j)
		{
			if (1 == touchEdges[j].direction) // forward edge
			{
				u = touchEdges[j].u;
				v = touchEdges[j].v;
			}
			else { // backward edge
				u = touchEdges[j].v;
				v = touchEdges[j].u;
			}
			touchEdges[j].weight = double(numShortestPathsBefore[u]) * numShortestPathsBackward[v];
			totalWeight += touchEdges[j].weight;
			if (debug) printf("touch edge: (%d, %d), weight: %.1lf\n", u, v, touchEdges[j].weight);
		}
		p = Uniform();
		if (debug) printf("p: %.3lf\n", p);
		p *= totalWeight;
		sumWeight = 0;

		for (j = 0; j < touchEdges.size(); ++j)
		{
			sumWeight += touchEdges[j].weight;
			if (sumWeight >= p) break;
		}
		assert(j != touchEdges.size());
		middleNode = touchEdges[j].v;
		if (debug) printf("chosen edge: (% d, % d), weight: % .1lf\n", touchEdges[j].u, touchEdges[j].v, touchEdges[j].weight);

		if (debug)
		{
			printf("pair: (%d, %d)\n predecessorsAdjLists:\n", src, onePair.t);
			//predecessorsAdjLists.printArray();
			printf("predecessorAdjListsBackwords:\n");
			//predecessorAdjListsBackwards.printArray();
			printf("middle node: %d\n", middleNode);
		}

		tempNode = middleNode;
		while (tempNode != src)
		{
			assert(tempNode >= 0 && tempNode < n);
			if (nodesInEachSampleSetS[i] >= maxNodesInPath) abort();
			sampleShortestPathsArraySetS[i][nodesInEachSampleSetS[i]] = tempNode;
			nodesInEachSampleSetS[i]++;
			p = Uniform();
			sz = predecessorsAdjListsSize[tempNode];
			if (sz == 0) {
				printf("error: pair: %d, %d\npredecessorsAdjLists:\n", onePair.s, onePair.t);
				//predecessorsAdjLists.printArray();
				printf("%d th sample\n", i + 1);
				abort();
			}
			delta = 1.0 / sz;
			j = int(floor(p / delta)); // round bin gaming
			if (j == sz) j--;
			assert(j >= 0 && j < sz);
			//if (debug) printf("---INSERT node %d\n", tempNode);

			tempNode = predecessorsAdjLists[tempNode][j];
		}
		if (nodesInEachSampleSetS[i] >= maxNodesInPath) abort();
		sampleShortestPathsArraySetS[i][nodesInEachSampleSetS[i]] = src;
		nodesInEachSampleSetS[i]++;
		//if (debug) printf("---INSERT node %d\n", src);

		tempNode = middleNode;
		if (tempNode == dst) continue;
		p = Uniform();
		sz = predecessorAdjListsBackwardsSize[tempNode];
		delta = 1.0 / sz;
		j = int(floor(p / delta)); // round bin gaming
		if (j == sz) j--;
		assert(j >= 0 && j < sz);
		tempNode = predecessorAdjListsBackwards[tempNode][j];

		while (tempNode != dst)
		{
			assert(tempNode >= 0 && tempNode < n);
			if (nodesInEachSampleSetS[i] >= maxNodesInPath) abort();
			sampleShortestPathsArraySetS[i][nodesInEachSampleSetS[i]] = tempNode;
			nodesInEachSampleSetS[i]++;
			p = Uniform();
			sz = predecessorAdjListsBackwardsSize[tempNode];
			if (sz == 0) {
				printf("error: pair: %d, %d \n predecessorAdjListsBackwards:\n", onePair.s, onePair.t);
				//predecessorAdjListsBackwards.printArray();
				printf("%d th sample\n", i + 1);
				abort();
			}
			delta = 1.0 / sz;
			j = int(floor(p / delta)); // round bin gaming
			if (j == sz) j--;
			//if (j < 0) j = 0;
			assert(j >= 0 && j < sz);
			//if (debug) printf("---INSERT node %d\n", tempNode);

			tempNode = predecessorAdjListsBackwards[tempNode][j];
		}
		if (nodesInEachSampleSetS[i] >= maxNodesInPath) abort();
		sampleShortestPathsArraySetS[i][nodesInEachSampleSetS[i]] = dst;
		nodesInEachSampleSetS[i]++;
		//if (debug) printf("---INSERT node %d\n", dst);
		if (debug)
		{
			printf("%d th sample: ", i + 1);
			for (j = 0; j < nodesInEachSampleSetS[i]; ++j)
				printf("%d,", sampleShortestPathsArraySetS[i][j]);
			printf(", middle Point: %d\n", middleNode);
		}
	}
	foundTimes.push_back((clock() - st) / 1000.0);
	numEdegesExplord /= L;
	//printf("bi BFS numEdegesExplord: %.0lf\n", numEdegesExplord);

	for (i = 0; i < n; ++i) isFoundforGBC[i] = false;
	for (k = 0; k < K; ++k)
	{
		//if (k % 10 == 0) printf("finding the %d node...\n", k + 1);
		for (i = 0; i < n; ++i) betweenCentralities[i] = 0;
		for (i = 0; i < L; ++i)
		{
			for (j = 0; j < nodesInEachSampleSetS[i]; ++j)
				if (isFoundforGBC[sampleShortestPathsArray[i][j]] == true) break;
			if (j != nodesInEachSampleSetS[i]) continue; // already covered
			for (j = 0; j < nodesInEachSampleSetS[i]; ++j)
				betweenCentralities[sampleShortestPathsArray[i][j]]++;
		}
		largestIncreasedBC = 0;
		for (i = 0; i < n; ++i)
		{
			if (betweenCentralities[i] > 0)
				assert(isFoundforGBC[i] == false);
			if (betweenCentralities[i] > largestIncreasedBC)
			{
				largestIncreasedBC = betweenCentralities[i];
				largestNodeID = i;
			}
		}
		if (largestIncreasedBC < 0.5) break;// no nodes are found
		isFoundforGBC[largestNodeID] = true;
		foundNodes.push_back(largestNodeID);
		totalGBC += largestIncreasedBC;

		if ((k + 1) % 10 == 0) foundTimes.push_back((clock() - st) / 1000.0);
		//printf("%d th node: %d, covered paths: %.3lf\n", k + 1, 	largestNodeID, largestIncreasedBC);
	}// end for (int k = 0; k < K; ++k)
	for (i = 0; i < n; ++i) // may be less than K nodes are chosen
	{
		if (foundNodes.size() >= K) break;
		if (isFoundforGBC[i] == true) continue;
		foundNodes.push_back(i); // arbitraries 
	}

	totalGBC = (totalGBC / L);
	for (i = 0; i < L; ++i) delete[]sampleShortestPathsArray[i];

	for (i = 0; i < n; ++i) {
		delete[]originalGraph[i];
		delete[] transposeGraph[i];
		delete[] predecessorsAdjLists[i];
		delete[] predecessorAdjListsBackwards[i];
	}
	delete[]originalGraphSize;
	delete[]transposeGraphSize;
	delete[]predecessorsAdjListsSize;
	delete[]predecessorAdjListsBackwardsSize;

	return totalGBC; // normalized GBC
}

double Graph::evaluateGBCbySampleShoretstPaths(vector<int>& foundNodes, vector<PAIR>& SamplePairs)
{
	clock_t  st = clock();
	int K = foundNodes.size();
	
	assert(K > 0 );
	int L = SamplePairs.size();
	assert(L <= numSamples);
	PAIR onePair;
	//printf("num of sampled paths: %d for evaluation\n", L);

	int i, j, k;
	int destID;
	bool debug = false;
	nonUniformArray<int>  predecessorsAdjLists; // record predecessors in the shortest paths
	vector<int> degrees; // degree of each node 
	degrees.reserve(n);
	for (i = 0; i < n; ++i) degrees.push_back(adjLists.size(i));
	predecessorsAdjLists.allocateMemory(degrees);

	int src, dst;
	int u, v;
	// the number of nodes that their shortest distances to src have been found
	//int numDistancesFound = 0;
	int sz;
	int* adj;

	double totalGBC = 0;

	//printf("num nodes: %d\n", n);
	int newDistance;
	double p;
	double delta;

	//int* sampleShortestPathsArray[numSamples];
	const int maxNodesInPath = 100;
	for (i = totalSampledShortestPathsSetT; i < L; ++i) sampleShortestPathsArraySetT[i] = new int[maxNodesInPath];

	degrees.clear();
	for (i = 0; i < n; ++i) degrees.push_back(0); // degrees in the transpose graph
	for (i = 0; i < n; ++i)
	{
		sz = adjLists.size(i);
		adj = &adjLists.access(i, 0);
		for (j = 0; j < sz; ++j)
			degrees[adj[j]]++;
	}
	nonUniformArray<int> transposeGraph;
	transposeGraph.allocateMemory(degrees);
	for (i = 0; i < n; ++i)
	{
		sz = adjLists.size(i);
		adj = &adjLists.access(i, 0);
		for (j = 0; j < sz; ++j)
			transposeGraph.push_back(adj[j], i);
	}
	nonUniformArray<int>  predecessorAdjListsBackwards;
	predecessorAdjListsBackwards.allocateMemory(degrees);
	//printf("original graph:\n"); adjLists.printArray();
	//printf("\n transpose graph:\n"); transposeGraph.printArray();


	clock_t ft = st;
	int middleNode;
	int tempNode;

	NodeInTwoQuques nodeU, nodeV;
	int stopBFSdistanceForward, stopBFSdistanceBackward;
	vector<EdgeWithWeight> touchEdges;
	double totalWeight;
	double sumWeight;

	vector<int> nodesInForwardBFS, nodesInBackwardBFS;
	nodesInForwardBFS.reserve(n);
	nodesInBackwardBFS.reserve(n);

	for (j = 0; j < n; ++j) predecessorsAdjLists.clear(j);
	for (j = 0; j < n; ++j) predecessorAdjListsBackwards.clear(j);
	for (j = 0; j < n; ++j) visited[j] = false;
	for (j = 0; j < n; ++j) visitedBackwards[j] = false;
	for (j = 0; j < n; ++j) numShortestPathsBefore[j] = 0;
	for (j = 0; j < n; ++j) numShortestPathsBackward[j] = 0;

	for (i = totalSampledShortestPathsSetT; i <  L; ++i) // randomly generate L shortest paths
	{
		//if (i == 40) debug = true;
		//else debug = false;
		//debug = true;

		onePair = SamplePairs[i];
		src = onePair.s;
		dst = onePair.t;

		/*if ((i + 1) % 100000 == 0) {
			printf("L: %d, find the %d th sample, time used: %.2lf, delta time: %.2lf\n", L, i + 1, (clock() - st) / 1000.0,
				(clock() - ft) / 1000.0);
			ft = clock();
		}*/
		

		// find the shortest paths from src 
		sz = nodesInForwardBFS.size();
		for (j = 0; j < sz; ++j) {
			v = nodesInForwardBFS[j];
			predecessorsAdjLists.clear(v);
			visited[v] = false;
			numShortestPathsBefore[v] = 0;
		}
		visited[src] = true;
		numShortestPathsBefore[src] = 1;
		nodesInForwardBFS.clear();
		nodesInForwardBFS.push_back(src);

		// find the shortest paths from dst 
		sz = nodesInBackwardBFS.size();
		for (j = 0; j < sz; ++j) {
			v = nodesInBackwardBFS[j];
			predecessorAdjListsBackwards.clear(v);
			visitedBackwards[v] = false;
			numShortestPathsBackward[v] = 0;
		}
		visitedBackwards[dst] = true;
		numShortestPathsBackward[dst] = 1;
		nodesInBackwardBFS.clear();
		nodesInBackwardBFS.push_back(dst);

		distances[src] = 0;
		distancesBackwards[dst] = 0;
		biQueue_clear(); // clear the queue
		biQueue_push(NodeInTwoQuques(src, 1));  //first queue
		biQueue_push(NodeInTwoQuques(dst, 2));  // second queue
		stopBFSdistanceForward = stopBFSdistanceBackward = n - 1;
		touchEdges.clear();
		//if (debug) printf("\npair(%d, %d)******************\n", src, dst);

		while (biQueue_size() > 0) // bfs
		{
			nodeU = biQueue_pop();
			u = nodeU.v;
			//if (debug) printf("exploring node %d ", u);

			if (1 == nodeU.qID) // forward search
			{
				//if (debug) printf("forward \n");
				if (distances[u] >= stopBFSdistanceForward) continue;

				sz = adjLists.size(u);
				adj = &adjLists.access(u, 0);
				newDistance = distances[u] + 1;

				for (j = 0; j < sz; ++j)
				{
					v = adj[j];
					if (false == visited[v])
					{
						visited[v] = true;
						distances[v] = newDistance;
						biQueue_push(NodeInTwoQuques(v, 1));
						nodesInForwardBFS.push_back(v);
						//if (debug) printf("node %d is added to queue\n", v);
					}
					if (distances[v] == newDistance) {
						predecessorsAdjLists.push_back(v, u);
						numShortestPathsBefore[v] += numShortestPathsBefore[u];
					}

					if (true == visitedBackwards[v])
					{
						stopBFSdistanceForward = newDistance;
						stopBFSdistanceBackward = distancesBackwards[v];
						touchEdges.push_back(EdgeWithWeight(u, v, 1, 0));  // forward edge
						//if (debug) printf("at node %d, found a middle node: %d\n", u, v);
					}
				}
			}
			else { // backward search
				//if (debug) printf("backward \n");
				if (distancesBackwards[u] >= stopBFSdistanceBackward) continue;
				sz = transposeGraph.size(u);
				adj = &transposeGraph.access(u, 0);
				newDistance = distancesBackwards[u] + 1;
				for (j = 0; j < sz; ++j)
				{
					v = adj[j];
					if (false == visitedBackwards[v])
					{
						visitedBackwards[v] = true;
						distancesBackwards[v] = newDistance;
						biQueue_push(NodeInTwoQuques(v, 2));
						nodesInBackwardBFS.push_back(v);
						//if (debug) printf("node %d is added to queue\n", v);
					}
					if (distancesBackwards[v] == newDistance) {
						predecessorAdjListsBackwards.push_back(v, u);
						numShortestPathsBackward[v] += numShortestPathsBackward[u];
					}

					if (true == visited[v])
					{
						stopBFSdistanceForward = distances[v];
						stopBFSdistanceBackward = newDistance;
						touchEdges.push_back(EdgeWithWeight(u, v, 2, 0)); // backward edge
						//if (debug) printf("at node %d, found a middle node: %d\n", u, v);
					}
				}
			}
		}
		destID = onePair.t;
		nodesInEachSampleSetT[i] = 0;
		if (touchEdges.size() == 0) {
			//printf("pair: %d, %d, not reachable\n", src, destID);
			//abort();
			//nodesInEachSample[i] = 0;
			continue;
		}

		totalWeight = 0;
		for (j = 0; j < touchEdges.size(); ++j)
		{
			if (1 == touchEdges[j].direction) // forward edge
			{
				u = touchEdges[j].u;
				v = touchEdges[j].v;
			}
			else { // backward edge
				u = touchEdges[j].v;
				v = touchEdges[j].u;
			}
			touchEdges[j].weight = double(numShortestPathsBefore[u]) * numShortestPathsBackward[v];
			totalWeight += touchEdges[j].weight;
			//if (debug) printf("touch edge: (%d, %d), weight: %.1lf\n", u, v, touchEdges[j].weight);
		}
		p = Uniform();
		//if (debug) printf("p: %.3lf\n", p);
		p *= totalWeight;
		sumWeight = 0;

		for (j = 0; j < touchEdges.size(); ++j)
		{
			sumWeight += touchEdges[j].weight;
			if (sumWeight >= p) break;
		}
		assert(j != touchEdges.size());
		middleNode = touchEdges[j].v;
		//if (debug) printf("chosen edge: (% d, % d), weight: % .1lf\n", touchEdges[j].u, touchEdges[j].v, touchEdges[j].weight);

		/*if (debug)
		{
			printf("pair: (%d, %d)\n predecessorsAdjLists:\n", src, onePair.t);
			predecessorsAdjLists.printArray();
			printf("predecessorAdjListsBackwords:\n");
			predecessorAdjListsBackwards.printArray();
			printf("middle node: %d\n", middleNode);
		}*/

		tempNode = middleNode;
		while (tempNode != src)
		{
			//assert(tempNode >= 0 && tempNode < n);
			//if (nodesInEachSample[i] >= maxNodesInPath) abort();
			sampleShortestPathsArraySetT[i][nodesInEachSampleSetT[i]] = tempNode;
			nodesInEachSampleSetT[i]++;
			p = Uniform();
			sz = predecessorsAdjLists.size(tempNode);
			/*if (sz == 0) {
				printf("error: pair: %d, %d\npredecessorsAdjLists:\n", onePair.s, onePair.t);
				//predecessorsAdjLists.printArray();
				printf("%d th sample\n", i + 1);
				abort();
			}*/
			delta = 1.0 / sz;
			j = int(floor(p / delta)); // round bin gaming
			if (j == sz) j--;
			//assert(j >= 0 && j < sz);
			//if (debug) printf("---INSERT node %d\n", tempNode);

			tempNode = predecessorsAdjLists.access(tempNode, j);
		}
		//if (nodesInEachSample[i] >= maxNodesInPath) abort();
		sampleShortestPathsArraySetT[i][nodesInEachSampleSetT[i]] = src;
		nodesInEachSampleSetT[i]++;
		//if (debug) printf("---INSERT node %d\n", src);

		tempNode = middleNode;
		if (tempNode == dst) continue;
		p = Uniform();
		sz = predecessorAdjListsBackwards.size(tempNode);
		delta = 1.0 / sz;
		j = int(floor(p / delta)); // round bin gaming
		if (j == sz) j--;
		//assert(j >= 0 && j < sz);
		tempNode = predecessorAdjListsBackwards.access(tempNode, j);

		while (tempNode != dst)
		{
			//assert(tempNode >= 0 && tempNode < n);
			//if (nodesInEachSample[i] >= maxNodesInPath) abort();
			sampleShortestPathsArraySetT[i][nodesInEachSampleSetT[i]] = tempNode;
			nodesInEachSampleSetT[i]++;
			p = Uniform();
			sz = predecessorAdjListsBackwards.size(tempNode);
			/*if (sz == 0) {
				printf("error: pair: %d, %d \n predecessorAdjListsBackwards:\n", onePair.s, onePair.t);
				//predecessorAdjListsBackwards.printArray();
				printf("%d th sample\n", i + 1);
				abort();
			}*/
			delta = 1.0 / sz;
			j = int(floor(p / delta)); // round bin gaming
			if (j == sz) j--;
			//if (j < 0) j = 0;
			//assert(j >= 0 && j < sz);
			//if (debug) printf("---INSERT node %d\n", tempNode);

			tempNode = predecessorAdjListsBackwards.access(tempNode, j);
		}
		//if (nodesInEachSample[i] >= maxNodesInPath) abort();
		sampleShortestPathsArraySetT[i][nodesInEachSampleSetT[i]] = dst;
		nodesInEachSampleSetT[i]++;
		//if (debug) printf("---INSERT node %d\n", dst);
		/*if (debug)
		{
			printf("%d th sample: ", i + 1);
			for (j = 0; j < nodesInEachSample[i]; ++j)
				printf("%d,", sampleShortestPathsArray[i][j]);
			printf(", middle Point: %d\n", middleNode);
		}*/
	}
	totalSampledShortestPathsSetT = L;
	
	//printf("bi BFS numEdegesExplord: %.0lf\n", numEdegesExplord);

	for (i = 0; i < n; ++i) isFoundforGBC[i] = false;
	for (i = 0; i < K; ++i) isFoundforGBC[foundNodes[i]] = true;

	bool isCovered;
	double numCoveredPaths = 0;
	for (i = 0; i < L; ++i)
	{
		isCovered = false;
		for (j = 0; j < nodesInEachSampleSetT[i]; ++j)
		{
			v = sampleShortestPathsArraySetT[i][j];
			if (true == isFoundforGBC[v]) {
				isCovered = true;
				break;
			}
		}
		if (true == isCovered)
			numCoveredPaths++;
	}
	totalGBC = numCoveredPaths / L; 
	
	//for (i = 0; i < L; ++i) delete[]sampleShortestPathsArray[i];
	return totalGBC; // normalized GBC
}

void Graph::generateRandomPairs(int n, int L, vector<PAIR>& randomizedPairs)
{// Post: generate L pairs, with indices between 0 and n-1
	assert(n > 0);
	randomizedPairs.clear();
	randomizedPairs.reserve(L);

	int i, j;
	int s, t;
	double p;
	PAIR pair;
	for (i = 0; i < L; ++i)
	{
		p = Uniform();
		s = int(floor(p * n));
		if (s == n) s--;
		pair.s = s;
		while (true) {
			p = Uniform();
			t = int(floor(p * n));
			if (t == n) t--;
			if (t != s) break;
		}
		pair.t = t;
		randomizedPairs.push_back(pair);
	}
	/*for (i = 0; i < L; ++i)
	{
		s = randomizedPairs[i].s;
		t = randomizedPairs[i].t;
		assert(s >= 0 && s < n);
		assert(t >= 0 && t < n && s!= t);
	}*/
}


double Graph::findTopKGBCAdaptiveSampling(int K, double epsilon, double errorProb, vector<int>& foundNodes, int& sampledPaths)
{
	assert(K > 0 && K < n);
	assert(epsilon > 0 && epsilon < 0.632);
	assert(errorProb > 0 && errorProb = 1);
	foundNodes.clear();
	foundNodes.reserve(K);

	double e = 2.71828;
	double alpha = epsilon / (2.0 - 1 / e);
	double c2 = (2 + alpha) / (alpha * alpha);
	double b = 3 * c2 + 2 + sqrt(18*c2+4);
	b /= 3 * c2 - 2;
	if (b < 1.1) b = 1.1;
	//b = 2; // for testing

	double Qmax = log(double(n) * (n - 1.0));
	Qmax = ceil( Qmax / log(b) );
	//printf("\nn: %d, base b: %.3lf, Qmax: %.0lf\n", n, b, Qmax);

	double d1;
	 d1 = log(2.0 / errorProb) + log(Qmax);

	double cnt = 0;
	int q;
	
	double Lq;
	double gq = 1.0;
	double theta = d1 * c2; // (ln 2/gamma+ln Qmax ) (2+alpha)/alpha^2
	Lq = theta;
	vector<PAIR> randomPairs;
	vector<int> curFoundNodes;
	vector<double> foundTime;
	double curBiasedGBC, curUnbiasedGBC;
	double totalSamples = 0;

	double epsilon1,  epsilonSum;
	double beta;
	double betaMax = epsilon / (1.0 - 1.0 / e);
	double timesLessThanOPT;
	double c1;

	totalSampledShortestPathsSetS = 0;
	totalSampledShortestPathsSetT = 0;
	for (q = 1; q <= Qmax; ++q)
	{
		gq /= b;
		Lq *= b;

		generateRandomPairs(n, int(Lq), randomPairs);
		curBiasedGBC = rankByGroupBCWithSampleShortestPathsBiBFS(randomPairs, curFoundNodes, K, foundTime);

		generateRandomPairs(n, int(Lq), randomPairs);
		curUnbiasedGBC = evaluateGBCbySampleShoretstPaths(curFoundNodes, randomPairs);
		beta = 1.0 - curUnbiasedGBC / curBiasedGBC;

		//printf("%dth try, biased gbc: %.4lf, unbiased gbc: %.4lf, gq: %.4lf, L: %d, totalSampledShortestPaths: %d\n", q, curBiasedGBC, curUnbiasedGBC, gq, (int)Lq, totalSampledShortestPathsSetS);
		if (curUnbiasedGBC >= gq) cnt++;
		if (cnt >= 2) {
			timesLessThanOPT = pow(b, cnt - 2);
			c1 = log(4.0 / errorProb) / (theta * timesLessThanOPT);
			epsilon1 = (2*c1/3.0 + sqrt(4*c1 * c1/9.0 + 8 * c1)) / 2.0;
			//epsilon2 = sqrt(2*log(4.0/errorProb)/(theta * timesLessThanOPT) );
			epsilonSum = beta * (1 - 1.0 / e) * (1 - epsilon1) + (2 - 1.0/e) * epsilon1;
			betaMax = 1.0 - (1 - 1 / e - epsilon + epsilon1) / ((1 - 1 / e) * (1 - epsilon1));
			//printf("cnt: %.0lf, epsilon1: %.4lf,  beta: %.4lf, betaMax: %4lf, epsilonSum: %.4lf\n", cnt,	epsilon1,  beta, betaMax, epsilonSum);
			if (epsilonSum > epsilon) continue;

			for (int k = 0; k < totalSampledShortestPathsSetS; ++k)
				delete[]sampleShortestPathsArraySetS[k];
			for (int k = 0; k < totalSampledShortestPathsSetT; ++k)
				delete[]sampleShortestPathsArraySetT[k];
			sampledPaths = totalSampledShortestPathsSetS + totalSampledShortestPathsSetT;
			for (int k = 0; k < K; ++k) foundNodes.push_back(curFoundNodes[k]);
			//sampledPaths = totalSamples;
			
			return curUnbiasedGBC;
		}
	}
}


/*
double Graph::findTopKGBCbyKDD23(int K, double epsilon, double errorProb, vector<int>& foundNodes, int& sampledPaths)
{
	assert(K > 0 && K < n);
	assert(epsilon > 0 && epsilon < 0.632);
	assert(errorProb > 0 && errorProb = 1);
	foundNodes.clear();
	foundNodes.reserve(K);

	double e = 2.71828;
	double alpha = epsilon / (2.0 - 1 / e);
	double c1 = (2 + alpha) / (alpha * alpha);

	double Qmax = ceil(log(double(n) * (n - 1.0)));
	double c2 = log(2.0 / errorProb) + log(Qmax);

	double cnt = 0;
	int q;
	double Lq;
	double gq = 1.0;
	Lq = (c2 + K * log( 3.0 * K)) * c1; // (ln 1/gamma+ln Qmax + K ln n) (2+alpha)/alpha^2
	vector<PAIR> randomPairs;
	vector<int> curFoundNodes;
	vector<double> foundTime;
	double curGBC;
	double totalSamples = 0;
	for (q = 1; q <= Qmax; ++q)
	{
		gq /= 2.0;
		Lq *= 2.0;
		totalSamples += Lq;

		generateRandomPairs(n, int(Lq), randomPairs);
		curGBC = rankByGroupBCWithSampleShortestPathsBiBFSunbiased(randomPairs, curFoundNodes, K, foundTime);

		//printf("%dth try, gbc: %.4lf, gq: %.4lf, Lq: %d\n", q, curGBC, gq, (int)Lq);
		if (curGBC >= gq) cnt++;
		if (cnt >= 2) {
			for (int k = 0; k < K; ++k) foundNodes.push_back(curFoundNodes[k]);
			sampledPaths = totalSamples;
			return curGBC;
		}
	}
}
/**/

double Graph::findTopKGBCbyKDD16(int K, double epsilon, double errorProb, vector<int>& foundNodes, int &sampledPaths)
{
	assert(K > 0 && K < n);
	assert(epsilon > 0 && epsilon < 0.632);
	assert(errorProb > 0 && errorProb = 1);
	foundNodes.clear();
	foundNodes.reserve(K);

	double e = 2.71828;
	double alpha = epsilon / (2.0 - 1 / e);
	double c1 = (2 + alpha) / (alpha * alpha);

	double Qmax = ceil(log(double(n) * (n - 1.0)));
	double c2 = log(2.0 / errorProb) + log(Qmax);

	double cnt = 0;
	int q;
	double Lq;
	double gq = 1.0;
	Lq = (c2 + K * log(double(n))) * c1; // (ln 1/gamma+ln Qmax + K ln n) (2+alpha)/alpha^2
	vector<PAIR> randomPairs;
	vector<int> curFoundNodes;
	vector<double> foundTime;
	double curGBC;
	double totalSamples = 0;

	for (q = 1; q <= Qmax; ++q)
	{
		gq /= 2.0;
		Lq *= 2.0;
		totalSamples += Lq;

		generateRandomPairs(n, int(Lq), randomPairs);
		curGBC = rankByGroupBCWithSampleShortestPathsBiBFS(randomPairs, curFoundNodes, K, foundTime);
		
		//printf("%dth try, gbc: %.4lf, gq: %.4lf, Lq: %d\n", q, curGBC, gq, (int)Lq);
		if (curGBC >= gq) cnt++;
		if (cnt >= 2){
			for (int k = 0; k < K; ++k) foundNodes.push_back(curFoundNodes[k]);
			sampledPaths = totalSamples;
			return curGBC;
		}
	}
}



double Graph::rankByGroupBCWithSampleShortestPathsBiBFS(vector<PAIR>& randomPairs, vector<int>& foundNodes, int K, vector<double>& foundTimes)
{ // find the sample with bidirectional BFS
	clock_t  st = clock();

	foundNodes.clear();
	foundNodes.reserve(K);
	foundTimes.clear();
	assert(K > 0 && K <= n);
	int L = randomPairs.size();
	assert(L <= numSamples);
	PAIR onePair;

	int i, j, k;
	int destID;
	bool debug = false;
	nonUniformArray<int>  predecessorsAdjLists; // record predecessors in the shortest paths
	vector<int> degrees; // degree of each node 
	degrees.reserve(n);
	for (i = 0; i < n; ++i) degrees.push_back(adjLists.size(i));
	predecessorsAdjLists.allocateMemory(degrees);

	int src, dst;
	int u, v;
	// the number of nodes that their shortest distances to src have been found
	//int numDistancesFound = 0;
	int sz;
	int* adj;

	double totalGBC = 0;
	double largestIncreasedBC;
	int largestNodeID;

	//printf("num nodes: %d\n", n);
	int newDistance;
	double p;
	double delta;

	
	const int maxNodesInPath = 100;
	for (i = totalSampledShortestPathsSetS; i < L; ++i) sampleShortestPathsArraySetS[i] = new int[maxNodesInPath];

	degrees.clear();
	for (i = 0; i < n; ++i) degrees.push_back(0); // degrees in the transpose graph
	for (i = 0; i < n; ++i)
	{
		sz = adjLists.size(i);
		adj = &adjLists.access(i, 0);
		for (j = 0; j < sz; ++j)
			degrees[adj[j]]++;
	}
	nonUniformArray<int> transposeGraph;
	transposeGraph.allocateMemory(degrees);
	for (i = 0; i < n; ++i)
	{
		sz = adjLists.size(i);
		adj = &adjLists.access(i, 0);
		for (j = 0; j < sz; ++j)
			transposeGraph.push_back(adj[j], i);
	}
	nonUniformArray<int>  predecessorAdjListsBackwards;
	predecessorAdjListsBackwards.allocateMemory(degrees);
	//printf("original graph:\n"); adjLists.printArray();
	//printf("\n transpose graph:\n"); transposeGraph.printArray();

	int middleNode;
	int tempNode;

	NodeInTwoQuques nodeU, nodeV;
	int stopBFSdistanceForward, stopBFSdistanceBackward;
	vector<EdgeWithWeight> touchEdges;
	touchEdges.reserve(100);
	double totalWeight;
	double sumWeight;

	clock_t ft = st;
	//printf("time used for initialization: %.3lf\n", (ft - st) / 1000.0);

	vector<int> nodesInForwardBFS, nodesInBackwardBFS;
	nodesInForwardBFS.reserve(n); 
	nodesInBackwardBFS.reserve(n);

	for (j = 0; j < n; ++j) predecessorsAdjLists.clear(j);
	for (j = 0; j < n; ++j) predecessorAdjListsBackwards.clear(j);
	for (j = 0; j < n; ++j) visited[j] = false;
	for (j = 0; j < n; ++j) visitedBackwards[j] = false;
	for (j = 0; j < n; ++j) numShortestPathsBefore[j] = 0;
	for (j = 0; j < n; ++j) numShortestPathsBackward[j] = 0;

	for (i = totalSampledShortestPathsSetS; i < L; ++i) // randomly generate L shortest paths
	{
		//if (i == 40) debug = true;
		//else debug = false;
		//debug = true;

		onePair = randomPairs[i];
		src = onePair.s;
		dst = onePair.t;

		if ((i+1) % 1000000 == 0)printf("L: %d, find the %d th sample, time used: %.2lf\n", L, i + 1, (clock() - st) / 1000.0);


		// find the shortest paths from src 
		sz = nodesInForwardBFS.size();
		for (j = 0; j < sz; ++j) {
			v = nodesInForwardBFS[j];
			predecessorsAdjLists.clear(v);
			visited[v] = false;
			numShortestPathsBefore[ v ] = 0;
		}
		visited[src] = true;
		numShortestPathsBefore[src] = 1;
		nodesInForwardBFS.clear();
		nodesInForwardBFS.push_back(src);
		
		// find the shortest paths from dst 
		sz = nodesInBackwardBFS.size();
		for (j = 0; j < sz; ++j) {
			v = nodesInBackwardBFS[j];
			predecessorAdjListsBackwards.clear(v);
			visitedBackwards[v] = false;
			numShortestPathsBackward[v] = 0;
		}
		visitedBackwards[dst] = true;
		numShortestPathsBackward[dst] = 1;
		nodesInBackwardBFS.clear();
		nodesInBackwardBFS.push_back(dst);

		distances[src] = 0;
		distancesBackwards[dst] = 0;
		biQueue_clear(); // clear the queue
		biQueue_push(NodeInTwoQuques(src, 1));  //first queue
		biQueue_push(NodeInTwoQuques(dst, 2));  // second queue
		stopBFSdistanceForward = stopBFSdistanceBackward = n - 1;
		touchEdges.clear();
		if( debug ) printf("\npair(%d, %d)******************\n", src, dst);

		while (biQueue_size() > 0) // bfs
		{
			nodeU = biQueue_pop();
			u = nodeU.v;
			//if (debug) printf("exploring node %d ", u);

			if (1 == nodeU.qID ) // forward search
			{
				//if (debug) printf("forward \n");
				if (distances[u] >= stopBFSdistanceForward) continue;

				sz = adjLists.size(u);
				adj = adjLists.startAddress(u);
				newDistance = distances[u] + 1;
				for (j = 0; j < sz; ++j)
				{
					v = adj[j];
					if (false == visited[v])
					{
						visited[v] = true;
						distances[v] = newDistance;
						biQueue_push(NodeInTwoQuques(v,1));
						nodesInForwardBFS.push_back(v);
						//if (debug) printf("node %d is added to queue\n", v);
					}
					if (distances[v] == newDistance) {
						predecessorsAdjLists.push_back(v, u);
						numShortestPathsBefore[v] += numShortestPathsBefore[u];
					}
						
					if (true == visitedBackwards[v])
					{
						stopBFSdistanceForward = newDistance;
						stopBFSdistanceBackward = distancesBackwards[v];
						touchEdges.push_back(EdgeWithWeight(u, v, 1, 0));  // forward edge
						//if (debug) printf("at node %d, found a middle node: %d\n", u, v);
					}
				}
			}
			else { // backward search
				//if (debug) printf("backward \n");
				if (distancesBackwards[ u ] >= stopBFSdistanceBackward) continue;
				sz = transposeGraph.size(u);
				adj = transposeGraph.startAddress(u);
				newDistance = distancesBackwards[u] + 1;
				for (j = 0; j < sz; ++j)
				{
					v = adj[j];
					if (false == visitedBackwards[v])
					{
						visitedBackwards[v] = true;
						distancesBackwards[v] = newDistance;
						biQueue_push(NodeInTwoQuques(v,2));
						nodesInBackwardBFS.push_back(v);
						//if (debug) printf("node %d is added to queue\n", v);
					}
					if (distancesBackwards[v] == newDistance) {
						predecessorAdjListsBackwards.push_back(v, u);
						numShortestPathsBackward[v] += numShortestPathsBackward[u];
					}
						
					if (true == visited[v])
					{
						stopBFSdistanceForward = distances[v];
						stopBFSdistanceBackward = newDistance;
						touchEdges.push_back(EdgeWithWeight(u, v, 2, 0)); // backward edge
						//if (debug) printf("at node %d, found a middle node: %d\n", u, v);	
					}
				}
			}
		}
		destID = onePair.t;
		nodesInEachSampleSetS[i] = 0;
		if (touchEdges.size() == 0) {
			//printf("pair: %d, %d, not reachable\n", src, destID);
			//abort();
			continue;
		}

		totalWeight = 0;
		for (j = 0; j < touchEdges.size(); ++j)
		{
			if (1 == touchEdges[j].direction) // forward edge
			{
				u = touchEdges[j].u;
				v = touchEdges[j].v;
			}else { // backward edge
				u = touchEdges[j].v;
				v = touchEdges[j].u;
			}
			touchEdges[j].weight = double(numShortestPathsBefore[u]) * numShortestPathsBackward[v];
			totalWeight += touchEdges[j].weight;
			//if (debug) printf("touch edge: (%d, %d), weight: %.1lf\n", u, v, touchEdges[j].weight);
		}
		p = Uniform();
		//if (debug) printf("p: %.3lf\n", p);
		p *= totalWeight;
		sumWeight = 0;
		for (j = 0; j < touchEdges.size(); ++j)
		{
			sumWeight += touchEdges[j].weight;
			if (sumWeight >= p) break;
		}
		//assert(j != touchEdges.size());
		middleNode = touchEdges[j].v;
		//if (debug) printf("chosen edge: (% d, % d), weight: % .1lf\n", touchEdges[j].u, touchEdges[j].v, touchEdges[j].weight);

		/*if (debug)
		{
			printf("pair: (%d, %d)\n predecessorsAdjLists:\n", src, onePair.t);
			predecessorsAdjLists.printArray();
			printf("predecessorAdjListsBackwords:\n");
			predecessorAdjListsBackwards.printArray();
			printf("middle node: %d\n", middleNode);
		}*/

		tempNode = middleNode;
		while (tempNode != src)
		{
			//assert(tempNode >= 0 && tempNode < n);
			//if (nodesInEachSample[i] >= maxNodesInPath) abort();
			sampleShortestPathsArraySetS[i][nodesInEachSampleSetS[i]] = tempNode;
			nodesInEachSampleSetS[i]++;
			p = Uniform();
			sz = predecessorsAdjLists.size(tempNode);
			/*if (sz == 0) {
				printf("error: pair: %d, %d\npredecessorsAdjLists:\n", onePair.s, onePair.t);
				//predecessorsAdjLists.printArray();
				printf("%d th sample\n", i + 1);
				abort();
			}*/
			delta = 1.0 / sz;
			j = int(floor(p / delta)); // round bin gaming
			if (j == sz) j--;
			//assert(j >= 0 && j < sz);
			//if (debug) printf("---INSERT node %d\n", tempNode);

			tempNode = predecessorsAdjLists.access(tempNode, j);
		}
		//if (nodesInEachSample[i] >= maxNodesInPath) abort();
		sampleShortestPathsArraySetS[i][nodesInEachSampleSetS[i]] = src;
		nodesInEachSampleSetS[i]++;
		//if (debug) printf("---INSERT node %d\n", src);

		tempNode = middleNode;
		if (tempNode == dst) continue;
		p = Uniform();
		sz = predecessorAdjListsBackwards.size(tempNode);
		delta = 1.0 / sz;
		j = int(floor(p / delta)); // round bin gaming
		if (j == sz) j--;
		//assert(j >= 0 && j < sz);
		tempNode = predecessorAdjListsBackwards.access(tempNode, j);

		while (tempNode != dst)
		{
			//assert(tempNode >= 0 && tempNode < n);
			//if (nodesInEachSample[i] >= maxNodesInPath) abort();
			sampleShortestPathsArraySetS[i][nodesInEachSampleSetS[i]] = tempNode;
			nodesInEachSampleSetS[i]++;
			p = Uniform();
			sz = predecessorAdjListsBackwards.size(tempNode);
			/*if (sz == 0) {
				printf("error: pair: %d, %d \n predecessorAdjListsBackwards:\n", onePair.s, onePair.t);
				//predecessorAdjListsBackwards.printArray();
				printf("%d th sample\n", i + 1);
				abort();
			}*/
			delta = 1.0 / sz;
			j = int(floor(p / delta)); // round bin gaming
			if (j == sz) j--;
			//if (j < 0) j = 0;
			//assert(j >= 0 && j < sz);
			//if (debug) printf("---INSERT node %d\n", tempNode);

			tempNode = predecessorAdjListsBackwards.access(tempNode, j);
		}
		//if (nodesInEachSample[i] >= maxNodesInPath) abort();
		sampleShortestPathsArraySetS[i][nodesInEachSampleSetS[i]] = dst;
		nodesInEachSampleSetS[i]++;
		//if (debug) printf("---INSERT node %d\n", dst);
		/*if (debug)
		{
			printf("%d th sample: ", i + 1);
			for (j = 0; j < nodesInEachSample[i]; ++j)
				printf("%d,", sampleShortestPathsArray[i][j]);
			printf(", middle Point: %d\n", middleNode);
		}*/
	}
	foundTimes.push_back((clock() - st) / 1000.0);
	totalSampledShortestPathsSetS = L;

	for (i = 0; i < n; ++i) isFoundforGBC[i] = false;
	int* pointerOnePath;
	for (k = 0; k < K; ++k)
	{
		//if ( (k+1) % 10 == 0) printf("finding the %d node...\n", k + 1);
		for (i = 0; i < n; ++i) betweenCentralities[i] = 0;
		for (i = 0; i < L; ++i)
		{
			pointerOnePath = sampleShortestPathsArraySetS[i];
			for (j = 0; j < nodesInEachSampleSetS[i]; ++j)
				if (isFoundforGBC[pointerOnePath[j]] == true) break;
			if (j != nodesInEachSampleSetS[i]) continue; // already covered
			for (j = 0; j < nodesInEachSampleSetS[i]; ++j)
				betweenCentralities[pointerOnePath[j]]++;
		}
		largestIncreasedBC = 0;
		for (i = 0; i < n; ++i)
		{
			//if (betweenCentralities[i] > 0)
				//assert(isFoundforGBC[i] == false);
			if (betweenCentralities[i] > largestIncreasedBC)
			{
				largestIncreasedBC = betweenCentralities[i];
				largestNodeID = i;
			}
		}
		if (largestIncreasedBC < 0.5) break;// no nodes are found
		isFoundforGBC[largestNodeID] = true;
		foundNodes.push_back(largestNodeID);
		totalGBC += largestIncreasedBC;

		if ((k + 1) % 10 == 0) foundTimes.push_back((clock() - st) / 1000.0);
		//printf("%d th node: %d, covered paths: %.3lf\n", k + 1, 	largestNodeID, largestIncreasedBC);
	}// end for (int k = 0; k < K; ++k)
	for (i = 0; i < n; ++i) // may be less than K nodes are chosen
	{
		if (foundNodes.size() >= K) break;
		if (isFoundforGBC[i] == true) continue;
		foundNodes.push_back(i); // arbitrarily
	}
	//printf("time used for finding: %.3lf\n", (clock() - ft) / 1000.0);

	//for (i = 0; i < L; ++i) delete[]sampleShortestPathsArray[i];

	totalGBC /= L;
	return totalGBC; // normalized GBC
}


/*
double Graph::rankByGroupBCWithSampleShortestPaths(vector<PAIR>& randomPairs, vector<int>& foundNodes, int K, vector<double>& foundTimes)
{
	foundNodes.clear();
	foundNodes.reserve(K);
	foundTimes.clear();
	assert(K > 0 && K <= n);
	int L = randomPairs.size();
	assert(L <= numSamples  );
	PAIR onePair;

	int i, j, k;
	int destID;
	bool debug = false;
	nonUniformArray<int>  predecessorsAdjLists; // record predecessors in the shortest paths
	vector<int> degrees; // degree of each node 
	degrees.reserve(n);
	for (i = 0; i < n; ++i) degrees.push_back(adjLists.size(i));
	predecessorsAdjLists.allocateMemory(degrees);
	int src;
	int u, v;
	// the number of nodes that their shortest distances to src have been found
	//int numDistancesFound = 0;
	int sz;
	int* adj;

	double totalGBC = 0;
	double largestIncreasedBC;
	int largestNodeID;

	//printf("num nodes: %d\n", n);
	int newDistance;
	double p;
	double delta;

	int* sampleShortestPathsArray[numSamples];
	const int maxNodesInPath = 100;
	for (i = 0; i < L; ++i) sampleShortestPathsArray[i] = new int[maxNodesInPath];

	clock_t  st = clock();
	clock_t ft = st;
	double numEdegesExplord = 0;
	for (i = 0; i < L; ++i) // randomly generate L shortest paths
	{
		onePair = randomPairs[i];
		src = onePair.s;
		
		if (i % 1000 == 0) {
			printf("L: %d, find the %d th sample, time used: %.2lf, delta time: %.2lf\n", L, i + 1, (clock()-st)/1000.0, 
				(clock() - ft) / 1000.0);
			ft = clock();
		}
		for (j = 0; j < n; ++j) predecessorsAdjLists.clear(j);
	
		// find the shortest paths from src 
		for (j = 0; j < n; ++j) visited[j] = false;
		visited[src] = true;
		distances[src] = 0;
		queue_clear(); // clear the queue
		queue_push(src);
		while (queue_size() > 0) // bfs
		{
			u = queue_pop();
			if (u == onePair.t) break;

			sz = adjLists.size(u);
			adj = &adjLists.access(u, 0);
			newDistance = distances[u] + 1;
			for (j = 0; j < sz; ++j)
			{
				v = adj[j];
				if (false == visited[v])
				{
					visited[v] = true;
					distances[v] = newDistance;
					queue_push(v);
				}
				if (distances[v] == newDistance) 
					predecessorsAdjLists.push_back(v, u);
			}
			numEdegesExplord += sz;
		}
		if (debug)
		{
			printf("pair: (%d, %d)\n predecessorsAdjLists:\n", src, onePair.t);
			predecessorsAdjLists.printArray();
		}
		destID = onePair.t;
		nodesInEachSample[i] = 0;
		if (false == visited[destID]){
			//printf("pair: %d, %d, not reachable\n", src, destID);
			//abort();
			sampleShortestPathsArray[i][0] = src;
			nodesInEachSample[i] = 1;
			continue;
 		}
		
		while (destID != src)
		{
			assert(destID >= 0 && destID < n);
			if (nodesInEachSample[i] >= maxNodesInPath) abort();
			sampleShortestPathsArray[i][ nodesInEachSample[i] ] = destID;
			nodesInEachSample[i]++;
			p = Uniform();
			sz = predecessorsAdjLists.size(destID);
			if (sz ==0) {
				//printf("error: pair: %d, %d\npredecessorsAdjLists:\n", onePair.s, onePair.t);
				//predecessorsAdjLists.printArray();
				abort();
			}
			delta = 1.0 / sz ;
			j = int (floor(p / delta)); // round bin gaming
			if (j == sz ) j--;
			//if (j < 0) j = 0;
			assert(j >= 0 && j < sz);

			destID = predecessorsAdjLists.access(destID, j);
		}
		if (nodesInEachSample[i] >= maxNodesInPath) abort();
		sampleShortestPathsArray[i][nodesInEachSample[i]] = src;
		nodesInEachSample[i]++;
		if (debug)
		{
			printf("%d th sample: ", i + 1);
			for (j = 0; j < nodesInEachSample[i]; ++j)
				printf("%d,", sampleShortestPathsArray[i][j]);
			printf("\n");
		}
	}
	foundTimes.push_back( (clock()-st)/1000.0);
	numEdegesExplord /= L;
	//printf("numEdegesExplord: %.0lf\n", numEdegesExplord);

	for (i = 0; i < n; ++i) isFoundforGBC[i] = false;
	for (k = 0; k < K; ++k)
	{
		//if( k%10 == 0) printf("finding the %d node...\n", k + 1);
		for (i = 0; i < n; ++i) betweenCentralities[i] = 0;
		for (i = 0; i < L; ++i)
		{
			for (j = 0; j < nodesInEachSample[i]; ++j)
				if (isFoundforGBC[ sampleShortestPathsArray[i][j] ] == true) break;
			if (j != nodesInEachSample[i]) continue; // already covered
			for (j = 0; j < nodesInEachSample[i]; ++j)
				betweenCentralities[ sampleShortestPathsArray[i][j] ]++;
		}
		largestIncreasedBC = 0;
		for (i = 0; i < n; ++i)
		{
			if (betweenCentralities[i] > 0)
				assert(isFoundforGBC[i] == false);
			if (betweenCentralities[i] > largestIncreasedBC)
			{
				largestIncreasedBC = betweenCentralities[i];
				largestNodeID = i;
			}
		}
		if (largestIncreasedBC < 0.5) break;// no nodes are found
		isFoundforGBC[largestNodeID] = true;
		foundNodes.push_back(largestNodeID);
		totalGBC += largestIncreasedBC;

		if( (k+1)%10 == 0 ) foundTimes.push_back((clock() - st) / 1000.0);
		//printf("%d th node: %d, covered paths: %.3lf\n", k + 1, 	largestNodeID, largestIncreasedBC);
	}// end for (int k = 0; k < K; ++k)
	for (i = 0; i < n; ++i) // may be less than K nodes are chosen
	{
		if (foundNodes.size() >= K) break;
		if (isFoundforGBC[i] == true) continue;
		foundNodes.push_back(i); // arbitraries 
	}

	totalGBC = (totalGBC / L);
	for (i = 0; i < L; ++i) delete []sampleShortestPathsArray[i];
	return totalGBC; // normalized GBC
}
/**/


void Graph::shortestDistances(int src, bool markSPT)
{
	int i, j;
	for (i = 0; i < n; ++i) visited[i] = false;

	visited[src] = true;
	distances[src] = 0;
	queue_clear(); // clear the queue
	queue_push(src);

	int u, v;
	// the number of nodes that their shortest distances to src have been found
	//int numDistancesFound = 0;
	int sz;
	int *adj; 
	while (queue_size() > 0 )
	{
		u = queue_pop();

		//++numDistancesFound;
		//if (debug) printf("node %d, dis: %d\n", u, (int)distances[u]);

		sz = adjLists.size(u);
		adj = & adjLists.access(u, 0);
		for (j = 0; j < sz; ++j)
		{
			//v = adjLists.access(u, j);
			if ( visited[ adj[j] ] ) continue;
			v = adj[ j ];
			visited[v] = true;
			distances[v] = distances[u] + 1;
			queue_push(v);
		}
	}
}

void Graph::allocateVariablesForDFS()
{
	assert(n > 0);
	if (NULL == ccBelongto)
	{
		numChildren = new int[n];
		numNodesSubtree = new int[n];
		parent = new int[n];
		isAP = new bool[n];
		discoveredTime = new int[n];
		lowestTime = new int[n];
		ccBelongto = new int[n];
		lowerBounds = new double[n];
	}

	for (int i = 0; i < n; ++i) // initialize
	{
		visited[i] = false;
		numChildren[i] = 0;
		numNodesSubtree[i] = 0;
		parent[i] = -1; // -1 indicates that the node is the root
		isAP[i] = false;
		lowerBounds[i] = 0;
	}
	time = 0;
}


void Graph::initilizeVariables()
{
	numChildren = NULL;
	numNodesSubtree = NULL;
	parent = NULL;
	isAP = NULL;
	discoveredTime = NULL;
	lowestTime = NULL;
	ccBelongto = NULL;
	lowerBounds = NULL;
	time = 0;
}

void Graph::clearVariables()
{
	if (NULL != ccBelongto)
	{
		delete []numChildren;
		delete []numNodesSubtree;
		delete []parent;
		delete []isAP;
		delete []discoveredTime;
		delete []lowestTime;
		delete []ccBelongto;
		delete []lowerBounds;
	}
	initilizeVariables();
}

void Graph::findAPsAndLowerBounds()
{
	nonUniformArray<int> allCCs;
	findAllCCsDFS(allCCs);
	
	int numNodesValid = 0; 
	int i;
	for (i = 0; i < allCCs.size(); ++i)
		numNodesValid += allCCs.size(i);

	allocateVariablesForDFS();

	for (i = 0; i < n; ++i)
	{
		if (false == isValid[i]) continue;
		if (true == visited[i]) continue;
		
		modifiedDFS(i, allCCs, numNodesValid);
		//if ( sizeOfCC > 10 ) printf("CC %d, size: %d\n", numofCCs, sizeOfCC);
	}
	
	double base = 0;
	for (i = 0; i < allCCs.size(); ++i)
		base += allCCs.size(i) * (numNodesValid - allCCs.size(i));
	base -= 2 * numNodesValid;

	//printf("the set of APs:\n");
	int numOfAPs = 0;
	for (i = 0; i < n; ++i)
	{
		if (false == isValid[i])continue;
		//if (false == isAP[i]) continue;
		++numOfAPs;
		lowerBounds[i] += base; 

		lowerBounds[i] *= zeta;
		//printf("%d ", i);
	}
	//printf("\n");
	//printf("num of APs: %d, percentage: %.2lf \n",  numOfAPs, 100.0*numOfAPs / n);

	/*printf("\nDFS tree: \n");
	for (i = 0; i < n; ++i)
		printf("parent of %d: %d\n", i, parent[i]);

	printf("\n");
	for (i = 0; i < n; ++i)
		printf("node %d, numofchildren: %d, numofdesc: %d\n", i, numChildren[i], numNodesSubtree[i]);

	printf("\n");
	for (i = 0; i < n; ++i)
		printf("node %d, discoverTime: %d, lowestTime: %d\n", i, discoveredTime[i], lowestTime[i]);
    */
}

int Graph::findAllCCsDFS(nonUniformArray<int> &allCCs)
{
	int i; 
	int numOfCCs = 0;
	sizeofEachCC.clear();
	sizeofEachCC.reserve( 1000 );

	allocateVariablesForDFS();

	int sizeOfCC;
	int maxCCsize = 0;
	for (i = 0; i < n; ++i)
	{
		if (false == isValid[i]) continue;
		if (true == visited[i]) continue;

		sizeOfCC = 0;
		DFS(i, numOfCCs, sizeOfCC); // dfs search
		//assert(sizeOfCC > 0);
		sizeofEachCC.push_back( sizeOfCC );
		++numOfCCs;
		//if ( sizeOfCC >= 30 )   printf("CC %d, size: %d\n", numOfCCs, sizeOfCC);
		if ( sizeOfCC > maxCCsize) maxCCsize = sizeOfCC;
	}
	if (debug) 
		printf("numofCCs: %d, maxCCsize: %d (%.2lf)\n", numOfCCs, maxCCsize, 100.0 * maxCCsize / n);
	
	if (debug)
		for (i = 0; i < n; ++i) printf("node: %d, in CC %d\n", i, ccBelongto[i]);

	allCCs.allocateMemory( sizeofEachCC );
	for (i = 0; i < n; ++i)
	{
		if (false == isValid[i]) continue; 
		allCCs.push_back(ccBelongto[i], i);
	}

	if (debug)
		for (i = 0; i < allCCs.size(); ++i)
		{
			printf("CC %d: ", i);
			for (int j = 0; j < allCCs.size(i); ++j)
				printf("%d ", allCCs.access(i, j ));
			printf("\n");
		}
	return numOfCCs;
}

int Graph::findAllCCsBFS(nonUniformArray<int> &allCCs)
// find all connected components
{
	int i, j;
	int numofCCs = 0;
	sizeofEachCC.clear();
	sizeofEachCC.reserve(1000);

	allocateVariablesForDFS();

	queue_clear();

	int sizeOfCC;
	int maxCCsize = 0;
	int u, v;
	int sz;
	int *adj;
	for (i = 0; i < n; ++i)
	{
		if (true == visited[i]) continue;
		if (false == isValid[i]) continue;

		// search the cc that contains vertex i
		sizeOfCC = 1;
		ccBelongto[ i ] = numofCCs;
		queue_push(i);
		while (queue_size() > 0 )
		{
			u = queue_pop();

			sz = adjLists.size(u);
			adj = & adjLists.access(u, 0);
			for (j = 0; j < sz; ++j)
			{
				v = adj[ j ];
				if (true == visited[v]) continue;
				if (false == isValid[v]) continue;

				queue_push(v);

				visited[ v ] = true;
				ccBelongto[ v ] = numofCCs;
				++sizeOfCC;
			}
		}
		assert(sizeOfCC > 0);

		sizeofEachCC.push_back(sizeOfCC);
		++numofCCs;
		//if ( sizeOfCC > 10 )   printf("CC %d, size: %d\n", numofCCs, sizeOfCC);
		if (sizeOfCC > maxCCsize) maxCCsize = sizeOfCC;
	}
	if (debug)
		printf("numofCCs: %d, maxCCsize: %d (%.2lf)\n", numofCCs, maxCCsize, 100.0 * maxCCsize / n);

	if (debug)
		for (i = 0; i < n; ++i) printf("node: %d, in CC %d\n", i, ccBelongto[i]);

	allCCs.allocateMemory( sizeofEachCC );
	for (i = 0; i < n; ++i)
	{
		if (false == isValid[i]) continue;
		allCCs.push_back(ccBelongto[i], i);
	}

	if (debug)
		for (i = 0; i < allCCs.size(); ++i)
		{
		printf("CC %d: ", i);
		for (int j = 0; j < allCCs.size(i); ++j)
			printf("%d ", allCCs.access(i, j));
		printf("\n");
		}
	return numofCCs;
}


inline void Graph::DFS(int u, int CCID, int & sizeOfCC)
{
	visited[u] = true;
	ccBelongto[u] = CCID;
	++sizeOfCC;

	int v;
	int sz = adjLists.size(u);
	int *adj = & adjLists.access(u, 0);
	for (int j = 0; j < sz ; ++j)
	{
		v = adj[ j ];
		//assert(v >= 0 && v < n);

		if (isValid[v] && !visited[v] )
		{
			parent[v] = u;
			DFS(v, CCID, sizeOfCC);
			numNodesSubtree[u] += numNodesSubtree[v];
		}
	}
	++numNodesSubtree[u]; // including itself
}

void Graph::modifiedDFS(int u, nonUniformArray<int> & allCCs, int numNodesValid)
{
	visited[ u ] = true;
	++time;
	discoveredTime[ u ] = lowestTime[ u ] = time;
	//printf("time %d, node %d\n", time, u);
	int ccSize = allCCs.size(ccBelongto[u]);
	int cc0 = ccSize - 1;

	int v;
	int j;
	int sz = adjLists.size(u);
	int *adj = &adjLists.access(u, 0);
	int nv; 
	for (j = 0; j < sz; ++j)
	{
		v = adj[j];
		if (false == isValid[v]) continue;

		if (false == visited[v] )
		{
			parent[ v ] = u;
			++numChildren[ u ];

			modifiedDFS(v, allCCs, numNodesValid);
			numNodesSubtree[u] += numNodesSubtree[v];
			
			if (lowestTime[v] < lowestTime[u])
				lowestTime[u] = lowestTime[v];
			// u is an articulation point in following cases

			// (1) u is root of DFS tree and has two or more chilren.
			if (parent[u] == -1 && numChildren[u] >= 2 ||
				// (2) If u is not root and low value of one of its child
				// is more than discovery value of u.
				parent[u] != -1 && lowestTime[v] >= discoveredTime[u])
			{
				isAP[u] = true;

				nv = numNodesSubtree[v];
				assert(nv < ccSize);
				lowerBounds[u] += nv * ( ccSize - 1 - nv);
				cc0 -= nv;
				assert(cc0 >= 0);
			}
		}
		else if (v != parent[u]) // a back edge is found
		{
			if (discoveredTime[v] < lowestTime[u])
				lowestTime[u] = discoveredTime[v];
		}
	}
	++numNodesSubtree[u]; // including itself

	if (true == isAP[u])
	{
		lowerBounds[u] += cc0 * (ccSize - 1 - cc0); // *zeta;
		lowerBounds[u] += 2 * ccSize;
	}else  // that is the upper bound for non-APs
		lowerBounds[u] = 2 * ccSize + 1; // *zeta;
	
}

void Graph::printGraph()
{
	int i, j;
	for (i = 0; i < n; ++i)
	{
		printf("node %d, weight:%f,  adjList: \n", i, nodesWeights[i]);
		for (j = 0; j < adjLists.size(i); ++j)
			printf("(%d, %d, %.4lf)\t", i, adjLists.access(i,j), adjWeights.access(i,j) );
		printf("\n\n");
	}
}

void Graph::printGraphTopology()
{
	int i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < adjLists.size(i); ++j)
			printf("(%d, %d)\t", i, adjLists.access(i, j));
		printf("\n\n");
	}
}

void Graph::cleanGraph(const char *graphFile)
{
	int i, j;
	ifstream in(graphFile);
	if (!in) {
		printf("Graph file `%s' was not found.\n", graphFile);
		return;
	}
	clearGraph();

	in >> n >> m;
	assert(n > 0 && m >= 0);
	assert(n <= MAX_NODES);

	vector<int> degrees; // the degree of each node 
	degrees.resize(n, 0);

	vector<Edge> allEdges;
	allEdges.reserve(m);

	int u, v;
	int numEdges = 0;
	for (i = 0; i < m; ++i)
	{
		if (i % 1000000 == 0) printf("%d (1000,000) th  edge \n", (int)(i / 1000000));

		in >> u >> v;

		if (u == v) {
			printf("A self loop occurs.\n");  
			continue;
		}
		if (u > v) swap(u, v); // make sure that u < v

		if (u < 0 || u >= n){ printf("u exceeds range: %d\n", u); continue; }
		if (v < 0 || v >= n){ printf("v exceeds range: %d\n", v); continue; }

		++degrees[u];
		allEdges.push_back(Edge(u, v));
	}
	in.close();

	adjLists.allocateMemory( degrees );
	for (i = 0; i < allEdges.size(); ++i)
	{
		u = allEdges[i].u;
		v = allEdges[i].v;
		assert(u < v);
		if (true == adjLists.find(u, v)) continue;
		//printf("Parallel edge (%d, %d) occurs\n", u, v); 
		adjLists.push_back(u, v); // insert vertex v to the adjlist of vertex u
		++numEdges;
	}
	printf("number of edges (before): %d, (after): %d\n", m, numEdges);
	m = numEdges;

	char dstFile[100];
	sprintf_s(dstFile, "datasets/%s", graphFile);
	ofstream out(dstFile);
	if (!out) abort();

	out << n << '\t' << m << '\n';
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < adjLists.size(i); ++j)
			out << i << '\t' << adjLists.access(i,j) << '\n';
	}
	out.close();

	printf("clean finished\n");
}


void Graph::generateEdgeNodeWeightsCaseStudy(double minNodeWeight, double maxNodeWeight)
{
	assert(0 <= minNodeWeight && minNodeWeight <= maxNodeWeight);
	vector<int> degrees; // degree of each node
	degrees.reserve(numOfNodes());
	int i, j;
	for (i = 0; i < adjLists.size(); ++i)
		degrees.push_back(adjLists.size(i));
	adjWeights.allocateMemory(degrees);
	double w;
	int u, v;
	for (i = 0; i < adjLists.size(); ++i)
	{
		u = i;
		for (j = 0; j < adjLists.size(i); ++j)
		{
			v = adjLists.access(i, j);
			if (u == 9 || v == 9)	w = 0.1;
			else w = 0.2;
			
			adjWeights.push_back(i, w);
			//if (i <= 1) printf("(%d,%d, %f)\n", i, adjLists.access(i, j), w);
		}
	}


	nodesWeights.clear();
	nodesWeights.reserve(numOfNodes());
	for (j = 0; j < numOfNodes(); ++j)
	{
		w = Uniform() *(maxNodeWeight - minNodeWeight) + minNodeWeight;
		nodesWeights.push_back(w);
		//if (j < 10) printf("(%d, %f)\n", j, w);
	}
	//printGraph();
}

void Graph::generateEdgeNodeWeights(double minEdgeWeight, double maxEdgeWeight, double minNodeWeight, double maxNodeWeight)
{
	assert(0 <= minEdgeWeight&& minEdgeWeight <= maxEdgeWeight && maxEdgeWeight <= 1.0);
	assert(0 <= minNodeWeight && minNodeWeight <= maxNodeWeight);
	vector<int> degrees; // degree of each node
	degrees.reserve(numOfNodes());
	int i, j;
	for (i = 0; i < adjLists.size(); ++i)
		degrees.push_back(adjLists.size(i));
	adjWeights.allocateMemory(degrees);
	double w;
	int degree;
	for (i = 0; i < adjLists.size(); ++i)
	{
		degree = adjLists.size(i);
		for (j = 0; j < adjLists.size(i); ++j)
		{
			w = Uniform() * (maxEdgeWeight - minEdgeWeight) + minEdgeWeight;
			adjWeights.push_back(i, w);
			//if (i <= 1) printf("(%d,%d, %f)\n", i, adjLists.access(i, j), w);
		}
	}
	
	nodesWeights.clear();
	nodesWeights.reserve(numOfNodes());
	for (j = 0; j < numOfNodes(); ++j)
	{
		w = Uniform() *(maxNodeWeight - minNodeWeight) + minNodeWeight;
		nodesWeights.push_back(w);
		//if (j < 10) printf("(%d, %f)\n", j, w);
	}
}

void Graph::readGraph(const char *graphFile, bool isDirected)
{
	FILE *f; 
	errno_t err = fopen_s(&f, graphFile, "r");
	if (err != 0)
	{
		printf("file: %s was not found!\n", graphFile);
		abort();
	}

	clearGraph();

	fscanf_s(f, "%d%d", &n, &m);
	
	assert(n > 0 && m >= 0);
	assert(n <= MAX_NODES);

	isValid = new bool[n];
	int i;
	for (i = 0; i < n; ++i) isValid[i] = true;

	vector<Edge> allEdges;
	allEdges.reserve(m);

	vector<int> degrees; // degree of each node
	degrees.resize(n, 0);

	int u, v;
	int numEdges = 0;
	for (i = 0; i < m; ++i)
	{
		if (i % 1000000 == 0 && i > 0 ) printf("%d (1000,000) th  edge \n", (int)(i / 1000000));

		fscanf_s(f, "%d%d", &u, &v);
		assert(u >= 0 && u < n);
		assert(v >= 0 && v < n);

		
		++numEdges;

		++degrees[ u ];
		if (false == isDirected) ++degrees[v];
		allEdges.push_back( Edge(u, v) );
	}
	fclose(f);
	assert(m == numEdges);

	adjLists.allocateMemory(degrees);
	assert( n == adjLists.size() );

	for (i = 0; i < allEdges.size(); ++i)
	{
		u = allEdges[i].u; 
		v = allEdges[i].v;
		adjLists.push_back(u, v);
		if (false == isDirected)
			adjLists.push_back(v, u);
	}
	vector<int> adjListSorting;
	adjListSorting.reserve(n);
	int j;
	for (i = 0; i < n; ++i)
	{
		adjListSorting.clear();
		for (j = 0; j < adjLists.size(i); ++j)
			adjListSorting.push_back(adjLists.access(i, j));
		sort(adjListSorting.begin(), adjListSorting.end());// sort in ascending order
		adjLists.clear(i);
		for (j = 0; j < adjListSorting.size(); ++j)
			adjLists.push_back(i, adjListSorting[j]);
	}
	printf("graph has been created.\n");
}

void Graph::clearGraph()
{
	if (NULL != isValid)  delete[]isValid;
	isValid = NULL;
	adjLists.clear();
	n = m = 0;
	clearVariables();
}

Graph::Graph()
{
	isValid = NULL;
	n = m = 0;
	debug = false;

	initilizeVariables();
}

Graph::~Graph()
{
	clearGraph();
	clearVariables();
}

Graph::Graph(const Graph &g)
{
	assert(g.n > 0);
	n = g.n; 
	m = g.m;
	adjLists = g.adjLists;
	isValid = new bool[n];
	for (int i = 0; i < n; ++i) isValid[i] = g.isValid[i];
	debug = g.debug;
	
	initilizeVariables();
}

Graph & Graph::operator=(const Graph &g)
{
	if (this == &g) return *this;

	clearGraph();
	clearVariables();

	assert(g.n > 0);
	n = g.n;
	m = g.m;
	adjLists = g.adjLists;

	isValid = new bool[n];
	for (int i = 0; i < n; ++i) isValid[i] = g.isValid[i];
	debug = g.debug;

	initilizeVariables();
	
	return *this;
}



