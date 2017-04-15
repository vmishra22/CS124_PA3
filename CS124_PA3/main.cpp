#include <iostream>
#include <sstream>
#include <time.h>
#include <limits.h>
#include <fstream>
#include <cstdlib>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <vector>
#include <random>
#include "BinHeap.h"

#define MAX_ITER 25000

using namespace std;

void generate_random_ints(string filename) {
	ofstream output_file(filename, std::ofstream::trunc);
	if (!output_file.is_open()) {
		cout << "Output file could not be opened! Terminating!" << endl;
		return;
	}

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<unsigned long long> dis(1, (unsigned long long)pow(10, 12));

	for (int n = 0; n<100; ++n)
		output_file << dis(gen) << endl;
	output_file.close();
}

unsigned long long KarmarkarKarpImplementation(vector<unsigned long long>& iElements) {
	unsigned long long retVal = 0;
	BinHeap<unsigned long long>* heap = new BinHeap<unsigned long long>(iElements.size());
	heap->buildHeap(iElements);

	while (heap->getHeapSize() > 1) {
		unsigned long long val1 = heap->extractMax();
		unsigned long long val2 = heap->extractMax();
		unsigned long long diff = val1 - val2;
		heap->insert(diff);
	}

	//After the diff loop only one element should remain in the heap. 
	retVal = heap->extractMax();
	return retVal;
}

int main(int argc, char** argv) {

	using namespace chrono;
	chrono::steady_clock::time_point tStartKK;
	chrono::steady_clock::time_point tEndKK;

	//STEP1-Run Karmarkar-Karp for one set of 100 random integers.
	//Read the input file for matrix elements.
	ifstream file;
	file.open(argv[1]);
	if (!file.is_open()) return -1;
	vector<unsigned long long> elements;
	unsigned long long elementInput;
	while (file >> elementInput){
		elements.push_back(elementInput);
	}

	unsigned long long residue = KarmarkarKarpImplementation(elements);
	cout << "Karmarkar Karp Residue: " << residue << endl;
	file.close();

	//Result.txt file for storing the timing and result of algorithm run.
	ofstream result_file("Result.txt", std::ofstream::trunc);
	if (!result_file.is_open()) {
		cout << "Output file could not be opened! Terminating!" << endl;
		return -1;
	}

	//Generate 100 random instances of the problem
	for (int i = 1; i <= 100; i++) {
		string fileName = "input" + to_string(i);
		generate_random_ints(fileName);
	}
	for (int i = 1; i <= 100; i++) {
		string fileName = "input" + to_string(i);
		file.open(fileName);
		if (!file.is_open()) return -1;
		vector<unsigned long long> numberSet;
		unsigned long long num;
		while (file >> num) {
			numberSet.push_back(num);
		}
		tStartKK = steady_clock::now();
		unsigned long long residueKK = KarmarkarKarpImplementation(numberSet);
		tEndKK = steady_clock::now();
		duration<double, std::milli> diffKK = tEndKK - tStartKK;

		result_file << "Karmarkar Karp Residue: " << residueKK << " Time: " << diffKK.count() << endl;
		file.close();
	}
	//Representation 1
	//Generate a random solution S pf +1 and -1 values
	vector<int> rand_Solution1;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> dis(0, 1);

	for (int n = 0; n < 100; ++n) {
		int val = dis(gen);
		if(val == 0)
			rand_Solution1.push_back(-1);
		else
			rand_Solution1.push_back(val);
	}

	//Repeated Random
	for (int i = 0; i < MAX_ITER; i++) {

	}

	result_file.close();
	return 0;
}
