
#include "graph.h"
using namespace std;

static int lambda = 2; //lambda core



class GroupBetweenesssCentrality{ 
public:
	void readGraph(const char *graphFile, bool isCaseStudy, double maxWeight);

	void readGBC(char *holeFileName, vector<int> & foundHoles, double *costs, double *maxCCsize, double *foundTimes);

	void testErrorProbByRandomGroupBetweenCentrality(int K, vector<int>& foundHoles, 
		vector<double>& foundTimes, const char*graphName);
	void generateRandomPairs(int n, int L, vector<PAIR>& randomizedPairs);

	void testErrorRatioBeta(const char* graphName); // the relative error between the biased and unbiased estimations

	void calcGBCupperBound(const char* graphName); // find a (1-1/e)-approximate solution

	void compareDifferentSamplePathsMethods(const char* graphName);
public: 
	void allocateMemoriesForGraphs();

	void createSCCGraph(int *sccID, bool debug);
	void createSCCTransposeGraph(bool debug);
	void createSCCUnderlyingGraph(int numSCCs);

	void calSumEachSCC();
	
	void findReachableNodes(nonUniformArray<int> &V_sources, bool debug);
	void findReachableNodes(int src);
	
	void removeAdjacentEdges(int chosenSCCid);

public: // for maxBlock
	
	void generateRandomPermulation(int n, vector<int> & randomizedSeq);
protected: // for maximum influence
	void ObtainTransformGraph(Graph &gT, Graph &g);

public:
	void findLargestCC(vector<int> &largestCC);

	int maxCCsize(Graph &g);
public:
	Graph originalGraph;
};

void GroupBetweenesssCentrality::createSCCGraph(int *sccID, bool debug)
{
	int n = SCCs.size();
	int i, j;
	
	for (i = 0; i < n; ++i) globalDegrees[i] = 0;

	int sccID1, sccID2;
	int sz; 
	int *adj; 
	int numOriginalGraph = propagationGraph.size();
	for (i = 0; i < numOriginalGraph; ++i)
	{
		sccID1 = sccID[i];
		sz = propagationGraph.size(i);
		adj = &propagationGraph.access(i, 0);
		for (j = 0; j < sz; ++j)
		{
			sccID2 = sccID[ adj[j] ];
			if (sccID1 == sccID2) continue;
			++globalDegrees[sccID1];
		}
	}
	sccGraph.allocateMemory(globalDegrees, n);

	for (i = 0; i < numOriginalGraph; ++i)
	{
		sccID1 = sccID[i];
		sz = propagationGraph.size(i);
		adj = &propagationGraph.access(i, 0);
		for (j = 0; j < sz; ++j)
		{
			sccID2 = sccID[adj[j]];
			if (sccID1 == sccID2) continue;
			if (false == sccGraph.find(sccID1, sccID2))
				sccGraph.push_back(sccID1, sccID2);
		}
	}
	int m = 0;
	for (i = 0; i < n; ++i)
		m += sccGraph.size(i);
	printf("size of scc graph: %d nodes, %d edges\n", sccGraph.size(), m);

	if (debug)
	{
		printf("------------------------\nscc graph: \n");
		for (i = 0; i < sccGraph.size(); ++i)
		{
			printf("supernode %d: ", i);
			for (j = 0; j < sccGraph.size(i); ++j)
				printf("%d, ", sccGraph.access(i, j));
			printf("\n");
		}
	}
}

void GroupBetweenesssCentrality::createSCCTransposeGraph(bool debug)
{
	int n = SCCs.size();
	int i, j;
	for (i = 0; i < n; ++i) globalDegrees[i] = 0;
	int sz;
	int *adj;
	for (i = 0; i < n; ++i)
	{
		sz = sccGraph.size(i);
		adj = &sccGraph.access(i, 0);
		for (j = 0; j < sz; ++j)
			++globalDegrees[adj[j]];
	}
	sccTransposeGraph.allocateMemory(globalDegrees, n);

	for (i = 0; i < n; ++i)
	{
		sz = sccGraph.size(i);
		adj = &sccGraph.access(i, 0);
		for (j = 0; j < sz; ++j)
			sccTransposeGraph.push_back(adj[j], i);
	}

	if (debug)
	{
		printf("transpose scc graph: \n");
		for (i = 0; i < n; ++i)
		{
			printf("node %d: ", i);
			sz = sccTransposeGraph.size(i);
			adj = &sccTransposeGraph.access(i, 0);
			for (j = 0; j < sz; ++j)
				printf("(%d, %d), ", i, adj[j]);
			printf("\n");
		}
	}
}


void GroupBetweenesssCentrality::createSCCUnderlyingGraph(int numSCCs)
{
	int i, j;
	int sz;

	for (i = 0; i < numSCCs; ++i)
		globalDegrees[i] = sccGraph.size(i); // the sum of outgoing neighbors;
	for (i = 0; i < numSCCs; ++i)
	{
		sz = sccGraph.size(i);
		for (j = 0; j < sz; ++j)
			++globalDegrees[sccGraph.access(i, j)]; // plus incoming neighbors
	}
		
	sccUnderlyingGraph.allocateMemory(globalDegrees, numSCCs);
	int *adj;
	for (i = 0; i < numSCCs; ++i)
	{
		sz = sccGraph.size(i);
		adj = &sccGraph.access(i, 0);
		for (j = 0; j < sz; ++j) // undirected graph
		{
			sccUnderlyingGraph.push_back(i, adj[j]); 
			sccUnderlyingGraph.push_back(adj[j], i);
		}
	}

	for (i = 0; i < numSCCs; ++i)
		globalDegrees[i] = sccUnderlyingGraph.size(i);
	parentSCCNeighbors.allocateMemory(globalDegrees, numSCCs);
	sccNeighborGroups.allocateMemory(globalDegrees, numSCCs);
	sccNeighborSizes.allocateMemory(globalDegrees, numSCCs);

	/*int j, sz;
	int *adj;
	printf("-----------------------\n scc underlying graph:\n");
	for (i = 0; i < numSCCs; ++i)
	{
		sz = sccUnderlyingGraph.size(i);
		adj = &sccUnderlyingGraph.access(i, 0);
		printf("super node %d: ", i);
		for (j = 0; j < sz; ++j)
			printf("%d, ", sccUnderlyingGraph.access(i, j));
		printf("\n");
	}*/
}

void GroupBetweenesssCentrality::allocateMemoriesForGraphs()
{
	int n = originalGraph.n;
	int i, j;
	vector<int> degrees;
	vector<int> inDegrees;
	degrees.reserve(n);
	inDegrees.resize(n, 0);
	for (i = 0; i < n; ++i)
	{
		degrees.push_back(originalGraph.adjLists.size(i));
		for (j = 0; j < originalGraph.adjLists.size(i); ++j)
			++inDegrees[originalGraph.adjLists.access(i,j)];
	}
		
	propagationGraph.allocateMemory(degrees);
	predecessorList.allocateMemory(inDegrees);

	int trivalVector = n;
	oneSCC.allocateMemory(&trivalVector, 1);
	SCCs.reserve(n);
	int numOfEdges = originalGraph.numOfEdges();
	sccGraph.reserve(2 * numOfEdges);
	sccTransposeGraph.reserve(2 * numOfEdges);
	sccUnderlyingGraph.reserve(2 * numOfEdges);
	parentSCCNeighbors.reserve(2 * numOfEdges);
	sccNeighborGroups.reserve(2 * numOfEdges);
	sccNeighborSizes.reserve(2 * numOfEdges);

	sortSNs.reserve(n);
}

void GroupBetweenesssCentrality::calSumEachSCC()
{
	int nSCCs = SCCs.size();
	int  i, j;
	double sum;
	int sz;
	int *adj;
	int outIndex;
	for (i = 0; i < nSCCs; ++i)
	{
		sz = SCCs.size(i);
		adj = &SCCs.access(i, 0);
		sum = 0;
		for (j = 0; j < sz; ++j)
		{
			outIndex = inIndex2OutIndex[ adj[j] ];
			sum += originalGraph.nodesWeights[outIndex];
		}
			
		sumEachSCC[i] = sum;
	}
}

 void GroupBetweenesssCentrality::findReachableNodes(int src)
{
	 int v, w;
	 int i, j;
	 int sz;
	 int *adj;

	 numReachSCCs = 1;
	 reachSCCsSet[0] = src;
	 queue_clear();
	 queue_push(src);
	 isReachable[src] = true;

	 while (queue_size() >0)
	 {
		 v = queue_pop();
		 sz = sccGraph.size(v);
		 adj = &sccGraph.access(v, 0);
		// numOpersDom += 6 * sz + 5;
		 for (i = 0; i < sz; ++i, ++adj)
		 {
			 w = *adj;
			 if (true == isReachable[w]) continue;
			 isReachable[w] = true;
			 queue_push(w);
			 reachSCCsSet[numReachSCCs] = w;
			 ++numReachSCCs;
			 //if (debug) printf("%d,", w);
		 }
	 }
	 //if (debug) printf("\n");

	 numReachNodes = 0;
	 for (i = 0; i < numReachSCCs; ++i)
	 {
		 v = reachSCCsSet[i];
		 isReachable[v] = false;// reset for later usage

		 sz = SCCs.size(v);
		 adj = &SCCs.access(v, 0);
		 for (j = 0; j < sz; ++j, ++numReachNodes)
			 reachableSet[numReachNodes] = adj[j];
	 }
}

void GroupBetweenesssCentrality::findReachableNodes(nonUniformArray<int> &V_sources, bool debug)
{
	int v, w;
	int i;
	int sz;
	int *adj;
	numReachNodes = V_sources.size(0);
	queue_clear();
	//if (debug) printf("reachable nodes: ");
	for (i = 0; i < numReachNodes; ++i)
	{
		v = V_sources.access(0, i);
		reachableSet[i] = v;
		isReachable[v] = true;
		queue_push(v);
		//if (debug) printf("%d,", v);
	}
	while (queue_size() >0 )
	{
		v = queue_pop();
		sz = propagationGraph.size(v);
		adj = & propagationGraph.access(v, 0);
		for (i = 0; i < sz; ++i)
		{
			w = adj[i];
			if (true == isReachable[w]) continue;
			isReachable[w] = true;
			queue_push(w);
			reachableSet[numReachNodes] = w;
			++numReachNodes;
			//if (debug) printf("%d,", w);
		}
	}
	//if (debug) printf("\n");
	for (i = 0; i < numReachNodes; ++i)
		isReachable[reachableSet[i]] = false;
}



void GroupBetweenesssCentrality::removeAdjacentEdges(int chosenSCCid)
{
	int i;
	int v;
	int sz = sccGraph.size(chosenSCCid);
	if (sccTransposeGraph.size(chosenSCCid) > 0) // there are incoming neighbors
	{
		for (i = 0; i< sz; ++i)
		{
			v = sccGraph.access(chosenSCCid, i);
			sccTransposeGraph.erase(v, chosenSCCid);
			sccUnderlyingGraph.erase(v, chosenSCCid);
			sccUnderlyingGraph.erase(chosenSCCid, v);
		}
		sccGraph.clear(chosenSCCid);
	}else{ // there are no incoming neighbors
		for (i = 0; i< sz; ++i)
		{
			v = sccGraph.access(chosenSCCid, i);
			sccTransposeGraph.erase(v, chosenSCCid);
			sccUnderlyingGraph.erase(v, chosenSCCid);
		}
		sccGraph.clear(chosenSCCid);
		sccUnderlyingGraph.clear(chosenSCCid);
	}
}


void GroupBetweenesssCentrality::generateRandomPairs(int n, int L, vector<PAIR>& randomizedPairs)
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
		while(true){
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

void GroupBetweenesssCentrality::generateRandomPermulation(int n, vector<int> & randomizedSeq)
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


void GroupBetweenesssCentrality::ObtainTransformGraph(Graph &gT, Graph &g)
{
	gT.n = g.n;
	gT.m = g.m;

	vector<int> degrees;
	degrees.reserve(g.numOfNodes());
	int i, j;
	for (i = 0; i < g.n; ++i)
		degrees.push_back(g.adjLists.size(i)); // it is bidirected graph, i.e., <u,v> implies that <v,u> in graph g
	gT.adjLists.allocateMemory(degrees);
	gT.adjWeights.allocateMemory(degrees);

	int u, v;
	double w;
	assert(g.n == g.adjLists.size());
	for (i = 0; i < g.adjLists.size(); ++i)
	{
		u = i;
		for (j = 0; j < g.adjLists.size(i); ++j)
		{
			v = g.adjLists.access(i, j);
			w = g.adjWeights.access(i, j);
			gT.adjLists.push_back(v, u);
			gT.adjWeights.push_back(v, w);
		}
	}
	gT.nodesWeights = g.nodesWeights;
}

void GroupBetweenesssCentrality::readGBC(char *holeFileName, vector<int> & foundHoles, double *costs, double *maxCCsize, double *foundTimes)
{
	int K = 50;
	ifstream in(holeFileName);
	if (!in) abort();
	int index;
	int hole;
	double c, ccSize, duration;
	foundHoles.clear();  foundHoles.reserve(K);
	int n = originalGraph.n;
	for (int i = 0; i < K; ++i)
	{
		in >> index >> hole >> c >> ccSize >> duration;
		foundHoles.push_back(hole);
		costs[i] = c;
		maxCCsize[i] = ccSize * (n - (i + 1));
		foundTimes[i] = duration;
	}
	in.close();
}

void GroupBetweenesssCentrality::readGraph(const char *graphFile, bool isCaseStudy, double maxWeight)
{
	originalGraph.clearGraph();
	originalGraph.clearVariables();

	originalGraph.readGraph(graphFile, false);
	if (false == isCaseStudy) originalGraph.generateEdgeNodeWeights(0, maxWeight, 0, 1);
	//else //originalGraph.generateEdgeNodeWeightsCaseStudy(0, 1);
	else originalGraph.generateEdgeNodeWeights(0, maxWeight, 1, 1);
	//originalGraph.printGraph();
	double n = originalGraph.n;
}



int GroupBetweenesssCentrality::maxCCsize(Graph &g)
{
	nonUniformArray<int>  allCCs;
	g.findAllCCsDFS(allCCs);
	//g.findAllCCsBFS(allCCs);
	int numOfCCs = allCCs.size();
	int maxCC = 0;
	for (int i = 0; i < numOfCCs; ++i)
		if (allCCs.size(i) > maxCC) maxCC = allCCs.size(i);

	return maxCC;
}



void GroupBetweenesssCentrality::compareDifferentSamplePathsMethods(const char* graphName)
{ 
	vector<int> largestCC;
	findLargestCC(largestCC); // find the largest CC in the original graph
	assert(largestCC.size() > 0);

	int i, j;
	int outIndex;
	int n = originalGraph.n;
	int* out2In = new int[n]; // index
	for (i = 0; i < n; ++i) out2In[i] = -1; //
	for (i = 0; i < largestCC.size(); ++i)
		out2In[largestCC[i]] = i;

	vector<int> degrees; //  the degree of each node 
	degrees.reserve(n);
	Graph curCC;
	curCC.n = largestCC.size();

	// allocate memory
	int deg;
	for (i = 0; i < largestCC.size(); ++i)
	{
		outIndex = largestCC[i];
		deg = originalGraph.adjLists.size(outIndex);
		degrees.push_back(deg); // degree of each node
	}
	curCC.adjLists.allocateMemory(degrees);

	// create the graph formed by the largest CC
	int v;
	int neighbor;
	int* adj;
	for (i = 0; i < largestCC.size(); ++i)
	{
		outIndex = largestCC[i];
		deg = originalGraph.adjLists.size(outIndex);
		adj = &originalGraph.adjLists.access(outIndex, 0);
		for (j = 0; j < deg; ++j)
		{
			//neighbor = originalGraph.adjLists.access(outIndex, j);
			neighbor = adj[j];
			if (false == originalGraph.isValid[neighbor]) continue;
			v = out2In[neighbor];
			//assert(v != -1);
			curCC.adjLists.push_back(i, v);
		}
	}

	vector<int> foundNodes;
	foundNodes.reserve(200);
	// find the centrality of each node in the largest CC

	double scaleDownFactor = double(curCC.n) * (curCC.n - 1.0) / originalGraph.n;
	scaleDownFactor /= (originalGraph.n - 1.0);
	printf("scaleDownFactor: %.4lf\n", scaleDownFactor);

	vector<PAIR> randomPairs;
	vector<double> foundTimes;

	clock_t st = clock();
	clock_t ft;

	double epsilon = 0.3;
	int K = 100;
	double epsilonKDD23, epsilonKDD16;

	double gamma = 0.01;
	
	double nodesFractionForEvaluation =0.005;
	int testTimes = 1;

	double curGBC=0;

	char comparisonOutGBCFileName[200];
	char comparisonOutSampleNumberFileName[200];

	int comparisonType = 1;  // 1: compare different values of K, 2: compare different values of epsilon
	
	if (1 == comparisonType) {
		sprintf_s(comparisonOutGBCFileName, "%s-comparisonK-epsion=%.1lf-GBC.txt", graphName, epsilon);
		sprintf_s(comparisonOutSampleNumberFileName, "%s-comparisonK-epsion=%.1lf-numPaths.txt", graphName, epsilon);
	}
	else {
		sprintf_s(comparisonOutGBCFileName, "%s-comparisonEpsilon-K=%d-GBC.txt", graphName, K);
		sprintf_s(comparisonOutSampleNumberFileName, "%s-comparisonEpsilon-K=%d-numPaths.txt", graphName, K);
	}

	/*vector<double> foundTime;
	int L = 10000;
	char comparisonOutRunTimeFileName[200];
	
	srand(11 + K);
	generateRandomPairs(curCC.n, L, randomPairs);

	st = clock();
	totalSampledShortestPathsSetS = 0;
	curGBC = curCC.rankByGroupBCWithSampleShortestPathsBiBFS(randomPairs, foundNodes, K, foundTime);
	for (j = 0; j < totalSampledShortestPathsSetS; ++j) delete[]sampleShortestPathsArraySetS[j];

	curGBC *= scaleDownFactor;
	ft = clock();
	printf("graph: %s, K: %d, curGBC: %.4lf, time: %.3lf\n", graphName, K, curGBC, (ft - st) / 1000.0);
	sprintf_s(comparisonOutRunTimeFileName, "%s-runningTime-L=%d.txt", graphName, L);
	ofstream outRunTime(comparisonOutRunTimeFileName, ios::app);
	outRunTime << L << '\t' << (ft - st) / 1000.0;
	outRunTime.close();
	/**/


	/*double alpha;
	double e = 2.71828;
	double gbc[101];
	gbc[20] = 0.185301;
	gbc[40] = 0.161855;
	gbc[50] = 0.18778;
	gbc[60]	= 0.192596;
	gbc[80]	= 0.213626;
	gbc[100] = 0.230036;


	epsilon = 0.3;
	for (K = 20; K <= 100; K += 20) gbc[K] /= scaleDownFactor;
	gbc[50] /= scaleDownFactor;
	double Qmax = ceil(log(double(n) * (n - 1.0)));
	double L;
	printf("n: %d, Qmax: %.0lf\n", curCC.n, Qmax);
	
	if (1 == comparisonType) { //1: compare different values of K
		epsilon = 0.3;
		printf("\n----epsilon: %.1lf\n", epsilon);
		ofstream outNumPaths(comparisonOutSampleNumberFileName, ios::app);
		for (K = 20; K <= 20; K += 20)
		{
			outNumPaths << K << '\t';
			alpha = epsilon / (2 - 1.0 / e);
			L = (log(2.0 / gamma) + K * log(double(curCC.n))) * (2 + alpha) / (alpha * alpha * gbc[K]);
			printf("KDD 16, K: %d,\t L: %.0lf\n", K, L);
			outNumPaths << L << '\t';

			L = (log(2.0 / gamma) + 0.8 * K * log(3.0 * K)) * (2 + alpha) / (alpha * alpha * gbc[K]);
			printf("KDD 23, K: %d,\t L: %.0lf\n", K, L);
			outNumPaths << L << '\n';

			L = (log(2.0 / gamma) + log (Qmax)  ) * (2 + alpha) / (alpha * alpha * gbc[K]);
			printf("Ada, K: %d,\t L: %.0lf\n\n", K, L);
		}
		outNumPaths.close();
	}
	else {
		K = 50;
		printf("\n----K: %d\n", K);
		ofstream outNumPaths(comparisonOutSampleNumberFileName, ios::app);
		for (epsilon=0.05; epsilon<= 0.5; epsilon+=0.05)
		{
			outNumPaths << epsilon << '\t';
			alpha = epsilon / (2 - 1.0 / e);
			L = (log(2.0 / gamma) + K * log(double(curCC.n))) * (2 + alpha) / (alpha * alpha * gbc[K]);
			printf("KDD 16, epsilon: %.2lf, L: %.0lf\n", epsilon, L);
			outNumPaths << L << '\t';

			L = (log(2.0 / gamma) + 0.8 * K * log( 3.0*K)) * (2 + alpha) / (alpha * alpha * gbc[K]);
			printf("KDD 23, epsilon: %.2lf, L: %.0lf\n\n", epsilon, L);
			outNumPaths << L << '\n';
		}
		outNumPaths.close();
	}
	*/
	
	int sampledPaths = 0;
	vector<int> foundSets;
	foundSets.reserve(K * testTimes);
	double minGBCKDD16 = 0, minGBCKDD23=0, minGBCada;
	double maxGBCKDD16=0, maxGBCKDD23 = 0, maxGBCada;
	double avgGBCKDD16=0, avgGBCKDD23 = 0, avgGBCada;
	double samplePathsKDD16=0, samplePathsKDD23 = 0, samplePathsAda;

	int printDelta = 10;

	//double opt = 0.185301 / scaleDownFactor;
	//printf("\nopt: %.4lf\n", opt);

	if (1 == comparisonType) { //1: compare different values of K
		epsilon = 0.3;
		printf("\nepsilon: %.1lf\n", epsilon);
		for (K = 20; K<= 100; K+=20)
		{
			printf("\nK: %d\n", K);
			
			foundSets.clear();

			samplePathsAda = 0;
			avgGBCada = 0;
			st = clock();
			for (i = 0; i < testTimes; ++i)
			{
				srand(11 + K +  i);
				curGBC = curCC.findTopKGBCAdaptiveSampling(K, epsilon, gamma, foundNodes, sampledPaths);
				
				avgGBCada += curGBC;
				for (j = 0; j < K; ++j) foundSets.push_back(foundNodes[j]);
				samplePathsAda += sampledPaths;
				if ((i + 1) % printDelta == 0)
					printf("Ada, %d th try, sampledPaths: %d, time used: %.0lf\n\n", i + 1, sampledPaths, (clock() - st) / 1000.0);
			}
			samplePathsAda /= testTimes;
			ft = clock();

			avgGBCada /= testTimes;
			avgGBCada *= scaleDownFactor;
 			
			printf("samplePathsAda: %.0lf, avgGBCada: %.4lf, time: %.4lf\n", samplePathsAda, avgGBCada, (ft-st)/1000.0);
			
			ofstream outNumPaths(comparisonOutSampleNumberFileName, ios::app);
			outNumPaths << K << '\t' << samplePathsAda << '\n';
			outNumPaths.close();
			
			/*samplePathsKDD23 = 0;
			epsilonKDD23 = epsilon * 0.9;
			for (i = 0; i < testTimes; ++i)
			{
				srand(11 + i);
				curGBC = curCC.findTopKGBCAdaptiveSampling(K, epsilonKDD23, gamma, foundNodes, sampledPaths);
				for (j = 0; j < K; ++j) foundSets.push_back(foundNodes[j]);
				samplePathsKDD23 += sampledPaths;
				if ((i + 1) % printDelta == 0)
					printf("KDD23, %d th try, sampledPaths: %d, time used: %.0lf\n", i + 1, sampledPaths, (clock() - st) / 1000.0);
			}
			samplePathsKDD23 /= testTimes;
			
			samplePathsKDD16 = 0;
			epsilonKDD16 = epsilonKDD23 * 0.9;
			for (i = 0; i < testTimes; ++i)
			{
				srand(11 + i);
				curGBC = curCC.findTopKGBCAdaptiveSampling(K, epsilonKDD16, gamma, foundNodes, sampledPaths);
				for (j = 0; j < K; ++j) foundSets.push_back(foundNodes[j]);
				samplePathsKDD16 += sampledPaths;
				if ((i + 1) % printDelta == 0)
					printf("KDD16, %d th try, sampledPaths: %d, time used: %.0lf\n", i + 1, sampledPaths, (clock() - st) / 1000.0);
			}
			samplePathsKDD16 /= testTimes;
			
			// max: KDD16, avg: KDD23, min: Ada
			avgGBCKDD23 = curCC.groupBC_BatchEvaluation(foundSets, K, avgGBCada,
				avgGBCKDD16, scaleDownFactor, NULL, 0, testTimes, nodesFractionForEvaluation, false);
			

			ofstream outGBC(comparisonOutGBCFileName, ios::app);
			ofstream outNumPaths(comparisonOutSampleNumberFileName, ios::app);
			outGBC << K << '\t'; // compare different values of K
			outGBC << avgGBCKDD16 << '\t' << avgGBCKDD23 << '\t' << avgGBCada << '\n';
			outNumPaths << K << '\t';
			outNumPaths << samplePathsKDD16 << '\t' << samplePathsKDD23 << '\t' << samplePathsAda << '\n';
			outGBC.close();
			outNumPaths.close();
	        /**/
			
		}
	}else if (2 == comparisonType) { // 2: compare different values of epsilon
		K = 100;
		for (epsilon = 0.1; epsilon <= 0.5; epsilon += 0.1) {
			printf("\nepsilon: %.1lf\n", epsilon);

			foundSets.clear();

			samplePathsAda = 0;
			avgGBCada = 0; 
			for (i = 0; i < testTimes; ++i)
			{
				srand(11 + K + i);
				curGBC = curCC.findTopKGBCAdaptiveSampling(K, epsilon, gamma, foundNodes, sampledPaths);
				avgGBCada += curGBC;
				samplePathsAda += sampledPaths;
				for (j = 0; j < K; ++j) foundSets.push_back(foundNodes[j]);
				if ((i + 1) % printDelta == 0)
					printf("Ada, %d th try, sampledPaths: %d, time used: %.0lf\n", i + 1, sampledPaths, (clock() - st) / 1000.0);
			}
			samplePathsAda /= testTimes;
			avgGBCada /= testTimes;
			avgGBCada *= scaleDownFactor;

			printf("samplePathsAda: %.0lf, avgGBCada: %.4lf\n", samplePathsAda, avgGBCada);
			
			ofstream outNumPaths(comparisonOutSampleNumberFileName, ios::app);
			outNumPaths << epsilon << '\t' << samplePathsAda << '\n';
			outNumPaths.close();
			printf("samplePathsAda: %.0lf\n", samplePathsAda);
			

			/*samplePathsKDD23 = 0;
			epsilonKDD23 = epsilon * 0.9;
			for (i = 0; i < testTimes; ++i)
			{
				srand(11 + i);
				curGBC = curCC.findTopKGBCAdaptiveSampling(K, epsilonKDD23, gamma, foundNodes, sampledPaths);
				for (j = 0; j < K; ++j) foundSets.push_back(foundNodes[j]);
				samplePathsKDD23 += sampledPaths;
				if ((i + 1) % printDelta == 0)
					printf("KDD23, %d th try, sampledPaths: %d, time used: %.0lf\n", i + 1, sampledPaths, (clock() - st) / 1000.0);
			}
			samplePathsKDD23 /= testTimes;
			
			samplePathsKDD16 = 0;
			epsilonKDD16 = epsilonKDD23 * 0.9;
			for (i = 0; i < testTimes; ++i)
			{
				srand(11 + i);
				curGBC = curCC.findTopKGBCAdaptiveSampling(K, epsilonKDD16, gamma, foundNodes, sampledPaths);
				for (j = 0; j < K; ++j) foundSets.push_back(foundNodes[j]);
				samplePathsKDD16 += sampledPaths;
				if ((i + 1) % printDelta == 0)
					printf("KDD16, %d th try, sampledPaths: %d, time used: %.0lf\n", i + 1, sampledPaths, (clock() - st) / 1000.0);
			}
			samplePathsKDD16 /= testTimes;
			
			// max: KDD16, avg: KDD23, min: Ada
			avgGBCKDD23 = curCC.groupBC_BatchEvaluation(foundSets, K, avgGBCada, avgGBCKDD16, scaleDownFactor, NULL, 0, testTimes, nodesFractionForEvaluation, false);

			ofstream outGBC(comparisonOutGBCFileName, ios::app);
			ofstream outNumPaths(comparisonOutSampleNumberFileName, ios::app);
			outGBC << epsilon << '\t'; // compare different values of K
			outGBC << avgGBCKDD16 << '\t' << avgGBCKDD23 << '\t' << avgGBCada << '\n';
			outNumPaths << epsilon << '\t';
			outNumPaths << samplePathsKDD16 << '\t' << samplePathsKDD23 << '\t' << samplePathsAda << '\n';
			outGBC.close();
			outNumPaths.close();
			/**/
			
		}
	}/**/

	delete[]out2In;
}


void GroupBetweenesssCentrality::calcGBCupperBound(const char* graphName)
{ // find a (1-1/e)-approximate solution
	vector<int> largestCC;
	findLargestCC(largestCC); // find the largest CC in the original graph
	assert(largestCC.size() > 0);

	int i, j;
	int outIndex;
	int n = originalGraph.n;
	int* out2In = new int[n]; // index
	for (i = 0; i < n; ++i) out2In[i] = -1; //
	for (i = 0; i < largestCC.size(); ++i)
		out2In[largestCC[i]] = i;

	vector<int> degrees; //  the degree of each node 
	degrees.reserve(n);
	Graph curCC;
	curCC.n = largestCC.size();

	// allocate memory
	int deg;
	for (i = 0; i < largestCC.size(); ++i)
	{
		outIndex = largestCC[i];
		deg = originalGraph.adjLists.size(outIndex);
		degrees.push_back(deg); // degree of each node
	}
	curCC.adjLists.allocateMemory(degrees);

	// create the graph formed by the largest CC
	int v;
	int neighbor;
	int* adj;
	for (i = 0; i < largestCC.size(); ++i)
	{
		outIndex = largestCC[i];
		deg = originalGraph.adjLists.size(outIndex);
		adj = &originalGraph.adjLists.access(outIndex, 0);
		for (j = 0; j < deg; ++j)
		{
			//neighbor = originalGraph.adjLists.access(outIndex, j);
			neighbor = adj[j];
			if (false == originalGraph.isValid[neighbor]) continue;
			v = out2In[neighbor];
			//assert(v != -1);
			curCC.adjLists.push_back(i, v);
		}
	}

	vector<int> foundNodes;
	foundNodes.reserve(200);
	// find the centrality of each node in the largest CC

	double scaleDownFactor = double(curCC.n) * (curCC.n - 1.0) / originalGraph.n;
	scaleDownFactor /= (originalGraph.n - 1.0);
	printf("scaleDownFactor: %.4lf\n", scaleDownFactor);

	vector<PAIR> randomPairs;

	clock_t st = clock();
	clock_t ft;

	double epsilon = 0.2;
	double gamma = 0.0001;

	int maxK = 100;
	int K;

	char upperBoundOutFileName[200];
	sprintf_s(upperBoundOutFileName, "%s-upperBound.txt", graphName);

	double numberChosenPaths;
	double delta = 100;
	int  sampledPaths;
	const int numDifferentKvalues = 5;
	int valuesOfK[numDifferentKvalues] = {20, 40, 60, 80, 100 };

	int testTimes = 10;
	int t;
	double maxGBC, curGBC;

	for (i =0; i< numDifferentKvalues; ++i)
	{
		K = valuesOfK[i];
		printf("\nK: %d, epsilon: %.2lf\n", K, epsilon);
		maxGBC = 0;
		for (t = 1; t <= testTimes; ++t) {
			srand(11 + K + t);
			curGBC = curCC.findTopKGBCAdaptiveSampling(K, epsilon, gamma, foundNodes, sampledPaths);
			curGBC *= scaleDownFactor;
			if (curGBC > maxGBC) maxGBC = curGBC;
			printf("K: %d, curGBC: %.6lf,  time used: %.2lf\n", K, curGBC, (clock() - st) / 1000.0);
		}
		printf("K: %d, maxGBC: %.6lf,  time used: %.2lf\n", K, maxGBC, (clock() - st) / 1000.0);
		ofstream out(upperBoundOutFileName, ios::app);
		out << K << '\t' << maxGBC << '\n';
		out.close();
	}

	delete[]out2In;
}



void GroupBetweenesssCentrality::testErrorRatioBeta(const char* graphName)
{ //  the relative error between the biased and unbiased estimations
	vector<int> largestCC;
	findLargestCC(largestCC); // find the largest CC in the original graph
	assert(largestCC.size() > 0);

	int i, j;
	int outIndex;
	int n = originalGraph.n;
	int* out2In = new int[n]; // index
	for (i = 0; i < n; ++i) out2In[i] = -1; //
	for (i = 0; i < largestCC.size(); ++i)
		out2In[largestCC[i]] = i;

	vector<int> degrees; //  the degree of each node 
	degrees.reserve(n);
	Graph curCC;
	curCC.n = largestCC.size();

	// allocate memory
	int deg;
	for (i = 0; i < largestCC.size(); ++i)
	{
		outIndex = largestCC[i];
		deg = originalGraph.adjLists.size(outIndex);
		degrees.push_back(deg); // degree of each node
	}
	curCC.adjLists.allocateMemory(degrees);

	// create the graph formed by the largest CC
	int v;
	int neighbor;
	int* adj;
	for (i = 0; i < largestCC.size(); ++i)
	{
		outIndex = largestCC[i];
		deg = originalGraph.adjLists.size(outIndex);
		adj = &originalGraph.adjLists.access(outIndex, 0);
		for (j = 0; j < deg; ++j)
		{
			//neighbor = originalGraph.adjLists.access(outIndex, j);
			neighbor = adj[j];
			if (false == originalGraph.isValid[neighbor]) continue;
			v = out2In[neighbor];
			//assert(v != -1);
			curCC.adjLists.push_back(i, v);
		}
	}

	vector<int> foundNodes;
	foundNodes.reserve(200);
	// find the centrality of each node in the largest CC

	double scaleDownFactor = double(curCC.n) * (curCC.n - 1.0) / originalGraph.n;
	scaleDownFactor /= (originalGraph.n - 1.0);
	printf("scaleDownFactor: %.4lf\n", scaleDownFactor);

	vector<PAIR> randomPairs;
	vector<double> foundTimes;

	clock_t st = clock();  
	clock_t ft;

	double curBiasedGBC;
	double curUnbiasedGBC;

	double  maxErrorRatioBeta;
	double avgErrorRatioBeta;
	double curErrorRatioBeta;

	int testTimes = 100;
	int K = 100;

	char errorRatioBetaOutFileName[200];
	sprintf_s(errorRatioBetaOutFileName, "%s-errorRatioBeta-K=%d-testTimes=%d.txt", graphName, K, testTimes);
	//ofstream out(errorRatioBetaOutFileName, ios::app);

	int numberChosenPaths;
	int minChosenPaths = 500, maxChosenPaths = 16000;
	for (numberChosenPaths = minChosenPaths; numberChosenPaths <= maxChosenPaths; numberChosenPaths *= 2)
	{
		printf("numberChosenPaths L: %d \n", numberChosenPaths);
		maxErrorRatioBeta = -10000;
		avgErrorRatioBeta = 0;
		for (i = 0; i < testTimes; ++i)
		{
			srand(i + 11);
			generateRandomPairs(curCC.n, numberChosenPaths, randomPairs);
			curBiasedGBC = curCC.rankByGroupBCWithSampleShortestPathsBiBFS(randomPairs, foundNodes, K, foundTimes);
			curBiasedGBC *= scaleDownFactor;
			
			srand(1997 + i);
			generateRandomPairs(curCC.n, numberChosenPaths, randomPairs);
			curUnbiasedGBC = curCC.evaluateGBCbySampleShoretstPaths(foundNodes, randomPairs);
			curUnbiasedGBC *= scaleDownFactor;
			if (curUnbiasedGBC > curBiasedGBC)  curUnbiasedGBC = curBiasedGBC;

			curErrorRatioBeta = 1.0 - curUnbiasedGBC / curBiasedGBC;
			avgErrorRatioBeta += curErrorRatioBeta;
			if (curErrorRatioBeta > maxErrorRatioBeta) maxErrorRatioBeta = curErrorRatioBeta;

			if ((i + 1) % 10 == 0)
				printf("maxL: %d, cur L: %d, testTimes: %d, cur Times: %d, time used: %.2lf\n",
					maxChosenPaths, numberChosenPaths, testTimes, i + 1, (clock() - st) / 1000.0);

			//printf("curBiasedGBC: %.6lf, curUnbiasedGBC: %.6lf, curErrorRatioBeta: %.4lf%%, maxErrorRatioBeta: %.4lf%%\n",curBiasedGBC, curUnbiasedGBC, 100 * curErrorRatioBeta, 100 * maxErrorRatioBeta);
		}
		avgErrorRatioBeta /= testTimes;

		//printf("avg ActualGBC: %.6lf, avg difference: %.6lf, ratio: %.3lf%%\n", sumActualGBC, difference, 100 * difference / sumActualGBC);
		printf("K: %d, L: %d, avgErrorRatioBeta: %.4lf%%, maxErrorRatioBeta: %.4lf%%\n", K, numberChosenPaths,
			100 * avgErrorRatioBeta, 100 * maxErrorRatioBeta);

		ofstream out(errorRatioBetaOutFileName, ios::app);
		out << numberChosenPaths << '\t' << avgErrorRatioBeta << '\t' << maxErrorRatioBeta << '\n';
		out.close();
	}

	delete[]out2In;
}


void GroupBetweenesssCentrality::testErrorProbByRandomGroupBetweenCentrality(int K, vector<int>& foundHoles, 
	vector<double>& foundTimes, const char* graphName)
{
	assert(K > 0 && K <= originalGraph.n);
	foundHoles.clear();
	foundHoles.reserve(K);
	foundTimes.clear();
	foundTimes.reserve(K);
	clock_t st = clock();

	vector<int> largestCC;
	findLargestCC(largestCC); // find the largest CC in the original graph
	assert(largestCC.size() > 0);

	int i, j;
	int outIndex;
	int n = originalGraph.n;
	int* out2In = new int[n]; // index
	for (i = 0; i < n; ++i) out2In[i] = -1; //
	for (i = 0; i < largestCC.size(); ++i)
		out2In[largestCC[i]] = i;

	vector<int> degrees; //  the degree of each node 
	degrees.reserve(n);
	Graph curCC;
	curCC.n = largestCC.size();
	curCC.isValid = new bool[curCC.n];
	for (i = 0; i < curCC.n; ++i) curCC.isValid[i] = true;

	// allocate memory
	int deg;
	for (i = 0; i < largestCC.size(); ++i)
	{
		outIndex = largestCC[i];
		deg = originalGraph.adjLists.size(outIndex);
		degrees.push_back(deg); // degree of each node
	}
	curCC.adjLists.allocateMemory(degrees);

	// create the graph formed by the largest CC
	int v;
	int neighbor;
	int* adj;
	for (i = 0; i < largestCC.size(); ++i)
	{
		outIndex = largestCC[i];
		deg = originalGraph.adjLists.size(outIndex);
		adj = &originalGraph.adjLists.access(outIndex, 0);
		for (j = 0; j < deg; ++j)
		{
			//neighbor = originalGraph.adjLists.access(outIndex, j);
			neighbor = adj[j];
			if (false == originalGraph.isValid[neighbor]) continue;
			v = out2In[neighbor];
			//assert(v != -1);
			curCC.adjLists.push_back(i, v);
		}
	}

	//originalGraph.printGraphTopology();
	//printf("\n");
	//curCC.printGraphTopology();

	vector<int> foundNodes;
	foundNodes.reserve(K);
	// find the centrality of each node in the largest CC
	vector<int> sourceNodes;
	sourceNodes.reserve(curCC.n);
	for (i = 0; i < curCC.n; ++i) sourceNodes.push_back(i);
	double totalGBC;
	double scaleDownFactor = double(curCC.n) * (curCC.n - 1.0) / originalGraph.n;
	scaleDownFactor /= (originalGraph.n - 1.0);
	printf("scaleDownFactor: %.4lf\n", scaleDownFactor);
	
	//totalGBC = curCC.rankByGroupBCWithSampleSources(sourceNodes, foundNodes, K);
	//totalGBC *= scaleDownFactor;

	
	//totalGBC = 0.2177; // top-10
	totalGBC = 0.496; // top-50
	printf("1-1/e ratio: gbc: %.4lf\n", totalGBC);

	double ratio = 1.0 - pow(1 - 1.0 / K, K); // no less than 1-1/e
	double opt_up = totalGBC / ratio;
	if (opt_up > 1) opt_up = 1;
	double errorRatio = 0.5;
	double alpha = errorRatio / (1 + ratio);
	double errorProb = 0.05;
	int numChosenNodes = (int) ceil((log(1.0/errorProb)+log( 2.0 )) *(2+alpha)/(alpha*alpha*opt_up));
	double errorThreshold = opt_up * (ratio - errorRatio);
	
	vector<int> randomizedSeq;
	randomizedSeq.reserve(curCC.n);
	int count = 0;
	double curGBC = 0;
	double estimatedGBC;
	
	int testTimes = 10;

	vector<PAIR> randomPairs;
	double minGBC = 0;
	vector<int> foundSets;
	foundSets.reserve(K* testTimes);
	double avergageGBC = 0;


	double nodesFractionForEvaluation = 1;

	st = clock();
	double difference = 0;
	double sumActualGBC = 0;

	char pathSampleFile[100] = "pathSampling";
	printf("pathSampling\n");
	int numberChosenPaths = numChosenNodes;
	numberChosenPaths = 200;
	
	//for (numberChosenPaths = 100; numberChosenPaths <= 1600.1; numberChosenPaths *= 2)
	 {
		printf("numberChosenPaths L: %d, threshold: %.3lf, \n", numberChosenPaths, errorThreshold);
		//foundSets.clear();
		for (i = 0; i < testTimes; ++i)
		{
			foundSets.clear();
			if (i % 10 == 0) printf("%dth simulation\n", i + 1);
			srand(i+11);
			generateRandomPairs(curCC.n, numberChosenPaths, randomPairs);
			//estimatedGBC = curCC.rankByGroupBCWithSampleShortestPaths(randomPairs, foundNodes, K, foundTimes);
			estimatedGBC = curCC.rankByGroupBCWithSampleShortestPathsBiBFS(randomPairs, foundNodes, K, foundTimes);
			//estimatedGBC = curCC.rankByGroupBCWithSampleShortestPathsBiBFSunbiasedKDD23(randomPairs, foundNodes, K, foundTimes);
			estimatedGBC *= scaleDownFactor;
			for (j = 0; j < K; ++j) foundSets.push_back(foundNodes[j]);
			printf("%d estimated GBC: %.6lf\n", i+1, estimatedGBC);
			
			srand(999 + numberChosenPaths);
			generateRandomPairs(curCC.n, numberChosenPaths, randomPairs);
			avergageGBC = curCC.evaluateGBCbySampleShoretstPaths(foundNodes, randomPairs);
			avergageGBC *= scaleDownFactor;
			//avergageGBC = curCC.groupBC_BatchEvaluation(foundSets, K, minGBC, scaleDownFactor, pathSampleFile, numberChosenPaths, testTimes, nodesFractionForEvaluation, false);
			printf("actual GBC:   \t %.6lf, ratio: %.2lf%%, minGBC: %.6lf, ratio: %.3lf %%\n\n", avergageGBC, 100 * avergageGBC / totalGBC, minGBC, 100 * minGBC / totalGBC);
			difference += abs(estimatedGBC - avergageGBC);
			sumActualGBC += avergageGBC;
		}
		
	}
	 difference /= testTimes;
	 sumActualGBC /= testTimes;
	 printf("avg ActualGBC: %.6lf, avg difference: %.6lf, ratio: %.3lf%%\n", sumActualGBC, difference, 100 * difference / sumActualGBC);
	 /**/



	clock_t ft = clock();
	printf("----------\ntime used: %.3lf,  normalized gbc: %.6lf\n----------\n", (ft - st) / 1000.0,  totalGBC);

	

	
	for (i = 0; i < K; ++i)
	{
		assert(foundNodes[i] >= 0 && foundNodes[i] < largestCC.size());
		outIndex = largestCC[foundNodes[i]];
		foundHoles.push_back(outIndex);
	}/**/

	/*double e = 2.71828;
	double L1, L2;
	double Lmax;
	for (double x = 1.2; x > 0.8; x -= 0.005)
	{
		L1 = (2 * x + 1.0) / (2 * x - 1.0);
		L1 /= (x - 0.5);
		Lmax = L1;

		L2 = 1.0 - e * x / (4 * e - 4.0);
		L2 = 2.0 / (L2 * L2);
		if (L2 > Lmax) Lmax = L2;
		printf("x: %.3lf, L1: %.6lf, L2: %.6lf, Lmax: %.6lf\n", x, L1, L2, Lmax);
	}*/

	delete[]out2In;

	double duration = double(clock() - st) / CLOCKS_PER_SEC;
	foundTimes.resize(K, duration);
}



void GroupBetweenesssCentrality::findLargestCC(vector<int> &largestCC)
{
	nonUniformArray<int> allCCs;
	originalGraph.findAllCCsBFS( allCCs );
	int maxCCid;
	int maxCC = 0;
	int i;
	for (i = 0; i < allCCs.size(); ++i)
	{
		if (allCCs.size(i) > maxCC)
		{
			maxCC = allCCs.size(i);
			maxCCid = i;
		}
	}
	largestCC.clear();
	largestCC.reserve( maxCC );
	for (i = 0; i < maxCC; ++i)
		largestCC.push_back( allCCs.access(maxCCid, i) );
}

