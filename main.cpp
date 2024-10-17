#include "GBC.h"
#pragma comment(linker, "/STACK:1500000000,1500000000") 

const int MAX_EDGES = 27440444;

void testVector();

void findPrime();

void testUniform();

int main(const int argc, const char *argv[])
{
	//  GR-QC Epinions Twitter  DBLP-2011 LiveJournal dblp testGraph    
	int dataset = 13;  // 0-10
	char graphNames[14][30] = { "GR-QC", "Epinions", "Twitter", "DBLP-2011",
		"LiveJournal", "dblp", "Email-euAll", "Facebook", "testGraph", "exampleGraph", "dominatorTestGraph", 
		"testGraphForBC", "syntheticNetwork-BA", "syntheticNetwork-WS"};
	// GR-QC: 0, 
	//facebook: 7, 
	//Twitter: 2, 
	//Email-euAll: 6, 
	//DBLP-11: 3, 
	//LiveJournal: 4, 
	//dblp: 5
	// Epinions: 1
	//syntheticNetwork-BA: 12
	// syntheticNetwork-WS: 13
	char testFileName[50];
	sprintf_s(testFileName, "datasets/%s.txt", graphNames[dataset]);

	const int maxMethods = 14;
	char methodNames[maxMethods][30] = { "maxBlock", "maxBlockFast", "PathCount", "2Step", "PageRank", "Constraint",
		"APGreedy", "maxInfluence", "HAM", "HIS", "MaxD", "shortestPath", 
		"BetweennessCentrality", "GroupBetweenCentrality"};
	

	srand(1);
	//srand(time(NULL));

	const int K = 50; // number of to-be-found nodes
	vector<int> foundGBC;
	vector<double> foundTimes; // record the found time for each hole node


	// 0: evaluation error ratio beta
	// 1: find a (1-1/e)-approximate value
	// 2: compare different shortest-path sampling methods
	// 12: betweeenness centrality
	// 13: group betweenness centrality
	int gbcMethod = 2;

	double maxWeight = 0.1;
	bool isOutputAll = false;

	printf("dataset: %s, maxweight: %.2f\n", graphNames[dataset], maxWeight);

	char suffix[100];
	if (gbcMethod == 0 || gbcMethod == 1 || gbcMethod == 7 || gbcMethod == 11)
		sprintf_s(suffix, "-%.2f", maxWeight );
	else sprintf_s(suffix, "");
	char outFileName[100];
	char blockedInfoEachNodeFileName[100];
	if (gbcMethod != -1)
	{
		sprintf_s(outFileName, "results/%s-%s%s.txt", graphNames[dataset], methodNames[gbcMethod], suffix);
		sprintf_s(blockedInfoEachNodeFileName, "results/%s-maxBlock%s-allNodes.txt", graphNames[dataset], suffix);
	}

	GroupBetweenesssCentrality GBCInstance;
	bool isCaseStudy = false;
	if (9 == dataset) isCaseStudy = true;
	
	
	GBCInstance.readGraph(testFileName, isCaseStudy, maxWeight);

	vector<Node> blockedInforAllNodes;

	clock_t st = clock();

	int n = GBCInstance.originalGraph.numOfNodes();

	
	float L = 2 * log2(1.0*n);
	if (true == isOutputAll)
	{
		L *= 10;
	}
	L = floor(L);
	L =1;

	// the 7 th dataset is 'dblp', notice that HIS and MaxD are supervised methods
	switch (gbcMethod)
	{
	case 0:
		GBCInstance.testErrorRatioBeta(graphNames[dataset]);
		break;
	case 1:
		GBCInstance.calcGBCupperBound(graphNames[dataset]);
		break;
	case 2:
		GBCInstance.compareDifferentSamplePathsMethods(graphNames[dataset]);
		break;
	case 3:
		break;
	case 4:
		break;
	case 5:
		break;
	case 6:
		break;
	case 7:
		break;
	case 8:
		break;
	case 9: // HIS
	case 10: // MAX
		break;
	case 11:
		break;
	case 12:
	
	case 13:
		GBCInstance.testErrorProbByRandomGroupBetweenCentrality(K, foundGBC, foundTimes, graphNames[dataset]);
	case -1:
		break;
	default:
		abort();
	}

	if (gbcMethod != -1)
	{
		//if (1 != holeMethod || false == isOutputAll)
		{
			ofstream out(outFileName, ios::app);
			for (int i = 0; i < foundGBC.size(); ++i)
				out << i + 1 << '\t' << foundGBC[i] << '\t' << foundTimes[i] << "\n";
			out.close();
		}
	}else{ // evaluate 
		int n = GBCInstance.originalGraph.numOfNodes();
		vector<int> allGBCMethods;
		allGBCMethods.reserve(11);
		int i, j;
		// we have results for the first eight algorithms for each dataset
		for (i = 0; i < 8; ++i)allGBCMethods.push_back(i); 
		if (dataset == 0 || dataset == 9) // dataset "QR-QC"
			allGBCMethods.push_back(8); // algorithm HAM
		if (dataset == 5) // dataset "dblp"
		{
			allGBCMethods.push_back(9); // algorithm HIS
			allGBCMethods.push_back(10); // algorithm MaxD
		}
		int numMethods =allGBCMethods.size();
		printf("methods: ");
		for (i = 0; i < numMethods; ++i) printf("%d,",allGBCMethods[i]);
		printf("\n");

		vector<int>allGBCNodes[11];
		vector<double>allGBCBlockedInfo[11];
		vector<float>allGBCInfluence[11];
		vector<float>allGBCFoundTime[11];
		for (i = 0; i < numMethods; ++i)allGBCNodes[i].reserve(K);
		for (i = 0; i < numMethods; ++i)allGBCBlockedInfo[i].reserve(K);
		for (i = 0; i < numMethods; ++i)allGBCInfluence[i].reserve(K);
		for (i = 0; i < numMethods; ++i)allGBCFoundTime[i].reserve(K);
		
		/**********--START--read candadiate hole nodes by each alg**********************/
		int curMethod;
		char resultFileName[100];
		int index, nodeId;
		float time;
		for (i = 0; i < numMethods; ++i)
		{
			curMethod =allGBCMethods[i];
			if (curMethod == 0 || curMethod == 1 || curMethod == 7)
				sprintf_s(suffix, "-%.2f.txt", maxWeight);
			else sprintf_s(suffix, ".txt");
			sprintf_s(resultFileName, "results/%s-%s%s", graphNames[dataset], methodNames[curMethod], suffix);
			//printf("%s\n", resultFileName);
			for (j = 0; j < n; ++j) visited[j] = false;
			ifstream in(resultFileName);
			for (j = 0; j < K; ++j)
			{
				in >> index >> nodeId >> time;
				if (false == visited[nodeId]) visited[nodeId] = true;
				else
				{
					printf("method %s, node: %d\n", methodNames[curMethod], nodeId);
					abort();
				}
				//if(j < 2) printf("%d, %d, %f\n", index, nodeId, time);
				allGBCNodes[i].push_back(nodeId);
				allGBCFoundTime[i].push_back(time);
			}
			in.close();
		}
		/**********--END--read candadiate hole nodes by each algorithm**********************/

		/**********--START--read blocked information propagations by each node****************/
		float* blockedInfo = new float[n];
		sprintf_s(blockedInfoEachNodeFileName, "results/%s-maxBlock-%.2f-allNodes.txt", graphNames[dataset], maxWeight);
		ifstream in(blockedInfoEachNodeFileName);
		int node;
		float w;
		for (i = 0; i < n; ++i)
		{
			in >> node >> w;
			blockedInfo[node] = w;
		}
		double blockedInfoSum;
		for (i = 0; i < numMethods; ++i)
		{
			printf("%s\n", methodNames[allGBCMethods[i]]);  // each algorithm
			blockedInfoSum = 0;
			for (j = 0; j < K; ++j)
			{
				node =allGBCNodes[i][j]; // the found node by the algorithm
				assert (blockedInfo[node] >= 0);
				
				blockedInfoSum += blockedInfo[node]; // the blocked info by the top-j nodes
				allGBCBlockedInfo[i].push_back(blockedInfoSum);
				if ((j + 1) % 10 == 0)
					printf("top %d nodes, block: %.3lf info\n", j + 1, blockedInfoSum);
			}
		}
		in.close();
		delete[]blockedInfo;
		/**********--END--read blocked information propagations by each node**********************/
	

		/******--START--output comparison results**********************************/
		
		ofstream compBlockedInfo, compInfluence, compTime;
		char fileNames[100];
		sprintf_s(fileNames, "results/comp-%s-%.2f-blockedInfo.txt", graphNames[dataset], maxWeight);
		compBlockedInfo.open(fileNames);
		sprintf_s(fileNames, "results/comp-%s-%.2f-influence.txt", graphNames[dataset], maxWeight);
		compInfluence.open(fileNames);
		sprintf_s(fileNames, "results/comp-%s-%.2f-time.txt", graphNames[dataset], maxWeight);
		compTime.open(fileNames);

		int k;
		for (k = 0; k < K; ++k)
		{
			if ((k + 1) % 5 != 0) continue; 
			compBlockedInfo << k + 1 << '\t';
			for (i = 0; i < numMethods; ++i)
				compBlockedInfo <<allGBCBlockedInfo[i][k] << '\t';
			compBlockedInfo << '\n';
		}

		for (k = 0; k < K; ++k)
		{
			if ((k + 1) % 5 != 0) continue;
			compInfluence << k + 1 << '\t';
			for (i = 0; i < numMethods; ++i)
				compInfluence <<allGBCInfluence[i][k] << '\t';
			compInfluence << '\n';
		}

		for (k = 0; k < K; ++k)
		{
			if ((k + 1) % 5 != 0) continue;
			compTime << k + 1 << '\t';
			for (i = 0; i < numMethods; ++i)
				compTime <<allGBCFoundTime[i][k] << '\t';
			compTime << '\n';
		}
		compBlockedInfo.close();
		compInfluence.close();
		compTime.close();
		/******--END--output comparison results**********************************/
	}
	                                
	printf("time: %.2lf (s)\n", (clock() - st) / 1000.0);
	printf("timeALG: %.3f (s)\n", timeAlg / 1000.0);
	
	
	int a;
	cin >> a;
	return 0;
}

void testUniform()
{
	const int maxNum = 100;
	 int arrayNum[maxNum];
	int i;
	for (i = 0; i < maxNum; ++i) arrayNum[i] = 0;
	int times = 50000;
	double p;
	int index;
	srand(time(NULL));
	for (i = 0; i < times * maxNum; ++i)
	{
		p = Uniform();
		index = floor(p * maxNum);
		if (index >= maxNum) index = maxNum - 1;
		++arrayNum[index];
	}
	for (i = 0; i < maxNum; ++i)
	{
		printf("%d\t times: %d\n", i, arrayNum[i]);
	}
}




void findPrime()
{
	const int maxNum = 32899;
	 bool arrayNum[maxNum];
	int i, j;
	for (i = 0; i < maxNum; ++i) arrayNum[i] = true;
	int step;
	for (i = 2; i < maxNum; ++i)
	{
		if (false == arrayNum[i]) continue;
		step = i;
		for (j = 2 * i; j < maxNum; j += step)
			arrayNum[j] = false;
	}
	for (i = maxNum-200; i < maxNum; ++i)
		if (true == arrayNum[i]) printf("%d\n", i);
}

void testEigenValues()
{
	/*const int n = 2;
	const int arraySize = n * n;
	double A[arraySize] = { 1, 1, 1, 0 };
	double eigenVectors[arraySize];
	double eigenValues[n];

	HAM ham;
	ham.JacbiCor(A, n, eigenVectors, eigenValues, 0.1, 2);
	int i, j;
	printf("eigen values: ");
	for (i = 0; i < n; ++i) printf("%.4lf\t", eigenValues[i]);
	printf("\n");
	printf("eigen vectors: \n");
	for (i = 0; i < n; ++i)
	{
	for (j = 0; j < n; ++j)
	printf("%.4lf\t", eigenVectors[i*n + j]);
	printf("\n");
	}*/
}

void testVector()
{
	
	const int SIZE = 1000000;
	int arr[SIZE];
	int *newPt = new int[SIZE];
	vector<int> testV;
	testV.reserve(SIZE);
	vector<int> v;
	v.push_back(SIZE);
	nonUniformArray<int> adj;
	adj.allocateMemory(v);
	int i, j;

	const int loopTimes = 1000;

	clock_t st;
	int re;

	st = clock();
	for (j = 0; j < SIZE; ++j)
		arr[j] = 0;
	for (i = 0; i < loopTimes; ++i)
	{
		for (j = 0; j < SIZE; ++j)
		{
			re = i * 30 + j;
			arr[j] = re;
		}
	}
	printf("time1: %d (ms), %d \n", (int)(clock() - st), arr[1]);

	st = clock();
	for (j = 0; j < SIZE; ++j)
		newPt[j] = 0;
	for (i = 0; i < loopTimes; ++i)
	{
		for (j = 0; j < SIZE; ++j)
		{
			re = i * 30 + j;
			newPt[j] = re;
		}
	}
	printf("time2: %d (ms), %d\n", (int)(clock() - st), newPt[1]);

	st = clock();
	for (j = 0; j < SIZE; ++j)
		testV.push_back(0);
	for (i = 0; i < loopTimes; ++i)
	{
		for (j = 0; j < SIZE; ++j)
		{
			re = i * 30 + j;
			testV[j] = re;
		}
	}
	printf("time3: %d (ms), %d\n", (int)(clock() - st), testV[1]);

	st = clock();
	for (j = 0; j < SIZE; ++j)
		adj.push_back(0, 0);
	for (i = 0; i < loopTimes; ++i)
	{
		for (j = 0; j < SIZE; ++j)
		{
			re = i * 30 + j;
			adj.access(0, j) = re;
		}
	}
	printf("time4: %d (ms), %d\n", (int)(clock() - st), adj.access(0, 1));

	delete[]newPt;
}