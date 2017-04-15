#pragma once
#ifndef _BIN_HEAP_
#define _BIN_HEAP_

#include <vector>
#include <iostream>
#include <algorithm>
#include <iterator>
using namespace std;

template <class T>
class BinHeap {
public: 
	BinHeap(int iSize): heapElements(iSize + 1){
		maxSize = iSize;
		heapSize = 0;
	}

	int getHeapSize() { return heapSize; }
	bool isEmpty();
	void insert(T& elem);
	T extractMax();
	void maxHeapify(int iPos);
	void buildHeap(vector<T>& iArray);
private:
	vector<T> heapElements;
	int heapSize;
	int maxSize;
};



template<class T>
bool BinHeap<T>::isEmpty() {
	return (heapSize == 0);
}

template<class T>
void BinHeap<T>::insert( T& elem)
{
	heapSize++;
	int N = heapSize;
	heapElements[heapSize] = elem;
	while ((N > 1) && (heapElements[N / 2] < heapElements[N])) {
		T temp = heapElements[N];
		heapElements[N] = heapElements[N / 2];
		heapElements[N / 2] = temp;
		N = (N / 2);
	}
}

template<class T>
T BinHeap<T>::extractMax()
{
	int root = 1;
	T max = heapElements[root];
	heapElements[root] = heapElements[heapSize--];
	maxHeapify(root);
	return max;
}


template<class T>
void BinHeap<T>::maxHeapify(int iPos)
{
	int l = 2 * iPos;
	int r = 2 * iPos + 1;
	int largest = 0;

	if ((l <= heapSize) && (heapElements[l] > heapElements[iPos])) {
		largest = l;
	}
	else {
		largest = iPos;
	}
	if ((r <= heapSize) && (heapElements[r] > heapElements[largest])) {
		largest = r;
	}
	if (largest != iPos) {
		T temp = heapElements[iPos];
		heapElements[iPos] = heapElements[largest];
		heapElements[largest] = temp;
		maxHeapify(largest);
	}
}

template<class T>
void BinHeap<T>::buildHeap(vector<T>& iArray)
{
	for (unsigned int j = 0; j < iArray.size(); j++) {
		heapElements[j + 1] = iArray[j];
		heapSize++;
	}
	for (int i = (iArray.size() / 2); i >= 1; i--)
	{
		maxHeapify(i);
	}
}
#endif
