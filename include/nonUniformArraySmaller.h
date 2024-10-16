#include<stdio.h>
#include<assert.h>
#include<vector>
#include<unordered_set>
#include "nonUniformArray.h"
using namespace std;

const size_t MAX_ROWS_OneCC = 800000;

//0: GrQc:             6000
//7: Facebook:        70000
//2: Twitter:        100000
//6: Email - euAll:  300000
//3: DBLP - 2011:   1000000
//4: LiveJournal:   5400000
//5: dblp:            60000
// 2: Epinions      80000
//12: syntheticNetwork-BA: 800000
// 13: syntheticNetwork-WS: 800000

// this class defines a data structure, called non-uniform array,
// which is a two-dimensional array, but the sizes of different rows
// may vary significantly.
// The reason for defining this data structure is that it is very time-comsuming
// to do the two operations of `new' and 'delete', since a 'new' and and a 'delete'
// operaton are need for the adjacent list of each vertex. 
// Therefore, there are O(n) such operations.
// The interal data structure of the non-uniform array in fact is a one-dementional array.
// Only O(1) 'new' and 'delete' operations are needed.
// Three basic methods: 
// (1) void allocateMemory(const vector<int> &sizesOfRows)
//      the vector tells the size of row, which are the maximum numbers of 
//        possible elements that can be stored at each row 
// (2) void push_back(int i, int element)
//     store the element at the end of row i
// (3) int & access(int i, int j)
//      access the (j+1) th element at row (i+1)
template <class Type>
class nonUniformArraySmaller{
public:
	void allocateMemory(const vector<int> &sizesOfRows);
	void allocateMemory(int *sizesOfRows, int rows);
	int size(); // how many rows
	int size(int i); // how many elements in row i
	// store the element at the end of row i
	inline void push_back(int i, Type element);
	void pop_back(int i, int sz); // pop the last sz elements
	Type & access(int i, int j); // access the element at row i and column j
	inline Type* startAddress(int i); // return the starting address of the first the element at row i
	void erase(int i, Type element);  // erase the element from row i
	bool find(int i, Type element); // whether the element is contained in row i
	void printArray();
	void clear(); // clear all  elements
	void clear(int i); // clear row i
	void reserve(int s);
public:
	nonUniformArraySmaller();
	~nonUniformArraySmaller();
	nonUniformArraySmaller(const nonUniformArraySmaller<Type> & other);
	nonUniformArraySmaller & operator=(const nonUniformArraySmaller<Type> & other);
protected:
	Type *storeArray; // store the two dimensional array into one sequential array
	int arraySize;
	//int maxSize[MAX_ROWS_OneCC]; // the max size of each row
	int sizeEachRow[MAX_ROWS_OneCC]; // how many elements stored at each row
	int baseEachRow[MAX_ROWS_OneCC]; // the base of each row
	int numRows; // total number of rows
};

template <class Type>
void nonUniformArraySmaller<Type>::reserve(int s)
{
	assert(s > 0);
	if (s <= arraySize) return;

	Type *newArray = new int[s];
	for (int i = 0; i < arraySize; ++i) // copy
		newArray[i] = storeArray[i]; 
	if (NULL != storeArray) delete[]storeArray;
	storeArray = newArray;
	arraySize = s;
}

template <class Type>
void nonUniformArraySmaller<Type>::clear(int i) // clear row i
{
	//assert(i >= 0 && i < numRows);
	sizeEachRow[i] = 0;
}

template <class Type>
void nonUniformArraySmaller<Type>::pop_back(int i, int sz)
{
	assert(i >= 0 && i < numRows);
	assert(sz >=0 && sz <= sizeEachRow[i]);
	sizeEachRow[i] -= sz;
}

template <class Type>
bool nonUniformArraySmaller<Type>::find(int i, Type element)
{
	//assert(i >= 0 && i < numRows);
	int sz = sizeEachRow[i];
	int base = baseEachRow[i];
	for (int j = 0; j < sz ; ++j)
		if (storeArray[ base + j] == element) return true;
	return false;
}

template <class Type>
void nonUniformArraySmaller<Type>::erase(int i, Type element)
{  // erase the element from row i
	assert(i >= 0 && i < numRows);
	int j;
	int base = baseEachRow[i];
	int sz = sizeEachRow[i];
	for (j = 0; j < sz; ++j)
		if (storeArray[ base + j ] == element) break;
	assert(j != sz); // the element must be found
	for (j; j <= sz - 2; ++j)
		storeArray[base + j] = storeArray[base + (j + 1)];
	--sizeEachRow[i];
}

template <class Type>
 Type & nonUniformArraySmaller<Type>::access(int i, int j)
{ // access the element at row i and column j
	//assert(i >= 0 && i < numRows);
	//assert(j >= 0 && j < sizeEachRow[i]);
	return storeArray[ baseEachRow[i] + j];
}

template <class Type>
inline Type* nonUniformArraySmaller<Type>::startAddress(int i)
{// return the starting address of the first the element at row i
	//assert(i >= 0 && i < numRows);
	return &storeArray[baseEachRow[i]];
}



template <class Type>
inline void nonUniformArraySmaller<Type>::push_back(int i, Type element)
{ // store the element at the end of row i
	//assert(i >= 0 && i < numRows);
	//assert(sizeEachRow[i] < maxSize[i]); // there are some storage left
	storeArray[ baseEachRow[i] + sizeEachRow[i]++ ] = element; // store the element
	//++sizeEachRow[i];
}

template <class Type>
int nonUniformArraySmaller<Type>::size(int i)
{ // how many elements in row i
	//assert(i >= 0 && i < numRows);
	return sizeEachRow[i];
}

template <class Type>
int nonUniformArraySmaller<Type>::size() // how many rows
{
	return numRows;
}

template <class Type>
void nonUniformArraySmaller<Type>::allocateMemory(const vector<int> &sizesOfRows)
{
	assert(sizesOfRows.size() > 0 &&
		sizesOfRows.size() <= MAX_ROWS_OneCC);
	
	int i;
	int numElements = 0;
	for (i = 0; i < sizesOfRows.size(); ++i)
	{
		//assert(sizesOfRows[i] >= 0);
		//maxSize[i] = sizesOfRows[i];
		numElements += sizesOfRows[i]; // count how many elements are needed
		sizeEachRow[i] = 0;
	}

	//assert(numElements <= 55000000);
	
	if (numElements > arraySize)
	{
		// clear allocated memory
		if ( NULL != storeArray) delete[] storeArray; 
		storeArray = new Type[numElements];
		arraySize = numElements;
	} // else 
	numRows = sizesOfRows.size();
	
	int base = 0;
	for (i = 0; i < numRows; ++i)
	{
		baseEachRow[i] = base;
		base += sizesOfRows[i];
	}
}

template <class Type>
void nonUniformArraySmaller<Type>::allocateMemory(int *sizesOfRows, int rows)
{
	//assert(rows > 0 && rows <= MAX_ROWS_OneCC);
	int i;
	int numElements = 0;
	for (i = 0; i < rows; ++i)
	{
		//assert(sizesOfRows[i] >= 0);
		//maxSize[i] = sizesOfRows[i];
		numElements += sizesOfRows[i]; // count how many elements are needed
		sizeEachRow[i] = 0;
	}
	//assert(numElements <= 55000000);

	if (numElements > arraySize)
	{
		// clear allocated memory
		if (NULL != storeArray) delete[] storeArray;
		storeArray = new Type[numElements];
		arraySize = numElements;
	} // else 
	numRows = rows;

	int base = 0;
	for (i = 0; i < numRows; ++i)
	{
		baseEachRow[i] = base;
		base += sizesOfRows[i];
	}
}

template <class Type>
nonUniformArraySmaller<Type> & nonUniformArraySmaller<Type>::operator = (const nonUniformArraySmaller<Type> & other)
{
	if (this == &other) return *this;

	// clear allocated memory
	if (NULL != storeArray) delete[]storeArray;

	if (other.arraySize == 0)
	{
		numRows = 0;
		storeArray = NULL;
		arraySize = 0;
		return *this;
	}

	numRows = other.numRows;
	arraySize = other.arraySize;
	int i;
	for (i = 0; i < numRows; ++i)
	{
		//maxSize[i] = other.maxSize[i];
		sizeEachRow[i] = other.sizeEachRow[i];
		baseEachRow[i] = other.baseEachRow[i];
	}
	storeArray = new Type[arraySize];
	for (i = 0; i < arraySize; ++i)
		storeArray[i] = other.storeArray[i];

	return *this;
}

template <class Type>
nonUniformArraySmaller<Type>::nonUniformArraySmaller(const nonUniformArraySmaller<Type> & other)
{
	if (other.arraySize == 0)
	{
		numRows = 0;
		arraySize = 0;
		storeArray = NULL;
		return;
	}

	numRows = other.numRows;
	arraySize = other.arraySize;
	int i;
	for (i = 0; i < numRows; ++i)
	{
		//maxSize[i] = other.maxSize[i];
		sizeEachRow[i] = other.sizeEachRow[i];
		baseEachRow[i] = other.baseEachRow[i];
	}
	
	storeArray = new Type[arraySize];
	for (i = 0; i < arraySize; ++i)
		storeArray[i] = other.storeArray[i];
}

template <class Type>
void nonUniformArraySmaller<Type>::printArray()
{
	assert(numRows > 0);
	int i, j;
	for (i = 0; i < numRows; ++i)
	{
		printf("row %d:\t", i);
		for (j = 0; j < sizeEachRow[i]; ++j)
			printf("%d ", storeArray[ baseEachRow[i]+j ] );
		printf("\n");
	}
}

template <class Type>
nonUniformArraySmaller<Type>::nonUniformArraySmaller()
{
	storeArray = NULL;
	arraySize = 0;
	numRows = 0;
}

template <class Type>
nonUniformArraySmaller<Type>::~nonUniformArraySmaller()
{
	if (NULL != storeArray) delete[] storeArray;
}

template <class Type>
void nonUniformArraySmaller<Type>::clear()
{
	numRows = 0;
}

