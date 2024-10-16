#include<assert.h>
#include<iostream> 
#include "nonUniformArraySmaller.h"
// define a traditional queue and a bidirectional queue

static const int MAX_NODES = MAX_ROWS;


static int queueArray[MAX_NODES];

static int q_size = 0;
static int q_first = 0; // point to the first and last elements, resp

inline int queue_size()
{ 
	return q_size; 
}

inline int queue_pop()
{
	//if (0 == q_size) abort();
	int element = queueArray[q_first];
	++q_first;
	--q_size;
	//if (0 == q_size) q_first = 0;
	//return queueArray[q_first-1];
	return element;
}

inline void queue_push(int element)
{
	//assert(q_first + q_size < MAX_NODES);

	queueArray[q_first + q_size] = element;

	++q_size;
}

inline void queue_clear()
{
	q_size = 0;
	q_first = 0;
}

struct NodeInTwoQuques{
	int v;
	int qID;
	NodeInTwoQuques(int Node,  int ququeID)
	{
		v = Node;
		qID = ququeID;
	}
	NodeInTwoQuques() {}
};

static NodeInTwoQuques biQueueArray[2*MAX_NODES];

static int biQ_size = 0;
static int biQ_first = 0; // point to the first and last elements, resp

inline int biQueue_size()
{
	return biQ_size;
}

inline NodeInTwoQuques biQueue_pop()
{
	NodeInTwoQuques element = biQueueArray[biQ_first];
	++biQ_first;
	--biQ_size;
	return element;
}

inline void biQueue_push(NodeInTwoQuques element)
{
	biQueueArray[biQ_first + biQ_size] = element;

	++biQ_size;
}

inline void biQueue_clear()
{
	biQ_size = 0;
	biQ_first = 0;
}




