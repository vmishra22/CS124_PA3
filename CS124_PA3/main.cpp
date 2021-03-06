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
using namespace chrono;

void generate_random_ints(string filename) {
	ofstream output_file(filename, std::ofstream::trunc);
	if (!output_file.is_open()) {
		cout << "Output file could not be opened! Terminating!" << endl;
		return;
	}

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<long long> dis(1, (long long)pow(10, 12));

	for (int n = 0; n<100; ++n)
		output_file << dis(gen) << endl;
	output_file.close();
}

long long KarmarkarKarpImplementation(vector<long long>& iElements) {
	long long retVal = 0;
	BinHeap<long long>* heap = new BinHeap<long long>(iElements.size());
	heap->buildHeap(iElements);

	while (heap->getHeapSize() > 1) {
		long long val1 = heap->extractMax();
		long long val2 = heap->extractMax();
		long long diff = val1 - val2;
		heap->insert(diff);
	}

	//After the diff loop only one element should remain in the heap. 
	retVal = heap->extractMax();

	delete heap;
	return retVal;
}

double T(int iter) {
	return (pow(10, 10)*(pow(0.8, iter/300)));
}

//Residue Calculate For Rep 1
long long CalculateResidue(vector<long long>& iNumberSet, vector<int>& iSolutionSet) {
	long long retVal = 0;
	int setSize = iNumberSet.size();
	for (int n = 0; n < setSize; ++n) {
		retVal += (iNumberSet[n] * iSolutionSet[n]);
	}
	return abs(retVal);
}

//Residue Calculate For Rep 2
long long CalculateResidueForPrePartition(vector<long long>& iNumberSet, vector<int>& iSolutionSet) {
	long long retVal = 0;
	int setSize = iNumberSet.size();
	vector<long long> iNumberSet2(setSize, 0);
	for (int j = 0; j < setSize; j++) {
		int solutionIndex = iSolutionSet[j];
		iNumberSet2[solutionIndex] += iNumberSet[j];
	}
	retVal = KarmarkarKarpImplementation(iNumberSet2);
	return retVal;
}

void RunAndAnalyzeKarmarkarKarpRep1(vector<long long>& iNumberSet, ofstream& oResultfile, vector<long long>& oResidues, vector<double> oTotalTime) {
	chrono::steady_clock::time_point tStartKK;
	chrono::steady_clock::time_point tEndKK;
	tStartKK = steady_clock::now();
	long long residueKK = KarmarkarKarpImplementation(iNumberSet);
	tEndKK = steady_clock::now();
	duration<double, std::milli> diffKK = tEndKK - tStartKK;
	oResidues[0] += residueKK;
	oTotalTime[0] += diffKK.count();
	oResultfile << "Karmarkar Karp Rep1: " << residueKK << " Time: " << diffKK.count() << endl;
}

//void RunAndAnalyzeKarmarkarKarpRep2(vector<long long>& iNumberSet, vector<int>& iInitialSolution, ofstream& oResultfile) {
//	chrono::steady_clock::time_point tStartKK;
//	chrono::steady_clock::time_point tEndKK;
//	tStartKK = steady_clock::now();
//	long long retVal = 0;
//	int setSize = iNumberSet.size();
//
//	vector<long long> iNumberSet2(setSize, 0);
//	for (int j = 0; j < setSize; j++) {
//		int solutionIndex = iInitialSolution[j];
//		iNumberSet2[solutionIndex] += iNumberSet[j];
//	}
//	retVal = KarmarkarKarpImplementation(iNumberSet2);
//	tEndKK = steady_clock::now();
//	duration<double, std::milli> diffKK = tEndKK - tStartKK;
//	oResultfile << "Karmarkar Karp Rep2: " << retVal << " Time: " << diffKK.count() << endl;
//}

//Repeated Random Rep 1
void RunAndAnalyzeRepeatedRandomRep1(vector<long long>& iNumberSet, vector<int>& iInitialSolution, ofstream& oResultfile, 
							  vector<long long>& oResidues, vector<double> oTotalTime) {
	//Generate a random solution S of +1 and -1 values
	chrono::steady_clock::time_point tStartRepReandom1;
	chrono::steady_clock::time_point tEndRepReandom1;
	tStartRepReandom1 = steady_clock::now();
	int setSize = iNumberSet.size(), val = 0;
	vector<int> rand_Solution1(iInitialSolution);
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<int> dis(0, 1);
	
	long long residue1 = CalculateResidue(iNumberSet, rand_Solution1);
	vector<int> rand_Solution2;

	for (int i = 0; i < MAX_ITER; i++) {
		//Create a new random solution S'.
		for (int n = 0; n < setSize; ++n) {
			val = dis(gen);
			if (val == 0)
				rand_Solution2.push_back(-1);
			else
				rand_Solution2.push_back(val);
		}
		
		long long residue2 = CalculateResidue(iNumberSet, rand_Solution2);
		if (residue2 < residue1) {
			rand_Solution1 = rand_Solution2;
			residue1 = residue2;
		}
		rand_Solution2.clear();
	}
	tEndRepReandom1 = steady_clock::now();
	duration<double, std::milli> diffRepReandom1 = tEndRepReandom1 - tStartRepReandom1;
	oResidues[1] += residue1;
	oTotalTime[1] += diffRepReandom1.count();
	oResultfile << "Repeated Random Rep1: " << residue1 << " Time: " << diffRepReandom1.count() << endl;
}

//Repeated Random Rep 2
void RunAndAnalyzeRepeatedRandomRep2(vector<long long>& iNumberSet, vector<int>& iInitialSolution, ofstream& oResultfile,
							  vector<long long>& oResidues, vector<double> oTotalTime) {
	//Generate a random solution P of [1,n]
	chrono::steady_clock::time_point tStartRepReandom1;
	chrono::steady_clock::time_point tEndRepReandom1;
	tStartRepReandom1 = steady_clock::now();
	int setSize = iNumberSet.size(), val = 0;
	vector<int> rand_Solution1(iInitialSolution);
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<int> dis(0, setSize-1);
	long long residue1 = CalculateResidueForPrePartition(iNumberSet, rand_Solution1);
	vector<int> rand_Solution2;
	for (int i = 0; i < MAX_ITER; i++) {
		//Create a new random solution P'.
		for (int n = 0; n < setSize; ++n) {
			val = dis(gen);
			rand_Solution2.push_back(val);
		}
		long long residue2 = CalculateResidueForPrePartition(iNumberSet, rand_Solution2);
		if (residue2 < residue1) {
			rand_Solution1 = rand_Solution2;
			residue1 = residue2;
		}
		rand_Solution2.clear();
	}
	tEndRepReandom1 = steady_clock::now();
	duration<double, std::milli> diffRepReandom1 = tEndRepReandom1 - tStartRepReandom1;
	oResidues[4] += residue1;
	oTotalTime[4] += diffRepReandom1.count();
	oResultfile << "Repeated Random Rep2: " << residue1 << " Time: " << diffRepReandom1.count() << endl;
}

//Hill Climbing Rep 1
void RunAndAnalyzeHillClimbingRep1(vector<long long>& iNumberSet, vector<int>& iInitialSolution, ofstream& oResultfile,
						  vector<long long>& oResidues, vector<double> oTotalTime) {
	chrono::steady_clock::time_point tStartHillClimb1;
	chrono::steady_clock::time_point tEndHillClimb1;
	tStartHillClimb1 = steady_clock::now();
	//Generate a random solution S of +1 and -1 values
	int setSize = iNumberSet.size();
	vector<int> rand_Solution1(iInitialSolution);
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<int> dis(0, 1);
	
	long long residue1 = CalculateResidue(iNumberSet, rand_Solution1);

	int idx1 = 0, idx2 = 0, swapToss = 0;
	vector<int> rand_Solution2;
	rand_Solution2 = rand_Solution1;
	//TO CHECK: if the distribution needs to be inside the loop
	uniform_int_distribution<int> disRep(0, setSize-1);
	for (int i = 0; i < MAX_ITER; i++){
		idx1 = disRep(gen);
		idx2 = disRep(gen);
		while (idx1 == idx2)
			idx2 = disRep(gen);
		//With probability 1/2 one element is moved, with probability 1/2 two elements are swapped/moved.
		rand_Solution2[idx1] *= -1;
		swapToss = dis(gen);
		if (swapToss == 0)
			rand_Solution2[idx2] *= -1;
		long long residue2 = CalculateResidue(iNumberSet, rand_Solution2);
		if (residue2 < residue1) {
			rand_Solution1 = rand_Solution2;
			residue1 = residue2;
		}
		else {
			rand_Solution2 = rand_Solution1;
		}
	}
	tEndHillClimb1 = steady_clock::now();
	duration<double, std::milli> diffHillClimb1 = tEndHillClimb1 - tStartHillClimb1;
	oResidues[2] += residue1;
	oTotalTime[2] += diffHillClimb1.count();
	oResultfile << "Hill Climb Rep1: " << residue1 << " Time: " << diffHillClimb1.count() << endl;
}

//Hill Climbing Rep 2
void RunAndAnalyzeHillClimbingRep2(vector<long long>& iNumberSet, vector<int>& iInitialSolution, ofstream& oResultfile, 
						  vector<long long>& oResidues, vector<double> oTotalTime) {
	chrono::steady_clock::time_point tStartHillClimb1;
	chrono::steady_clock::time_point tEndHillClimb1;
	tStartHillClimb1 = steady_clock::now();
	//Generate a random solution S of +1 and -1 values
	int setSize = iNumberSet.size();
	vector<int> rand_Solution1(iInitialSolution);
	random_device rd;
	mt19937 gen(rd());
	
	long long residue1 = CalculateResidueForPrePartition(iNumberSet, rand_Solution1);
	int idx1 = 0, idx2 = 0;
	vector<int> rand_Solution2;
	rand_Solution2 = rand_Solution1;
	//TO CHECK: if the distribution needs to be inside the loop
	uniform_int_distribution<int> disRep(0, setSize - 1);
	for (int i = 0; i < MAX_ITER; i++) {
		idx1 = disRep(gen);
		idx2 = disRep(gen);
		while (rand_Solution1[idx1] == idx2)
			idx2 = disRep(gen);
		
		rand_Solution2[idx1] = idx2;
		long long residue2 = CalculateResidueForPrePartition(iNumberSet, rand_Solution2);
		if (residue2 < residue1) {
			rand_Solution1 = rand_Solution2;
			residue1 = residue2;
		}
		else {
			rand_Solution2 = rand_Solution1;
		}
	}
	tEndHillClimb1 = steady_clock::now();
	duration<double, std::milli> diffHillClimb1 = tEndHillClimb1 - tStartHillClimb1;
	oResidues[5] += residue1;
	oTotalTime[5] += diffHillClimb1.count();
	oResultfile << "Hill Climb Rep2: " << residue1 << " Time: " << diffHillClimb1.count() << endl;
}

//Simulated Annealing Rep 1
void RunAndAnalyzeSimulatedAnnealingRep1(vector<long long>& iNumberSet, vector<int>& iInitialSolution, ofstream& oResultfile, 
							   vector<long long>& oResidues, vector<double> oTotalTime) {
	chrono::steady_clock::time_point tStartSimulatedAnnealing1;
	chrono::steady_clock::time_point tEndSimulatedAnnealing1;
	tStartSimulatedAnnealing1 = steady_clock::now();
	//Generate a random solution S of +1 and -1 values
	int setSize = iNumberSet.size();
	vector<int> rand_Solution1(iInitialSolution);
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<int> dis(0, 1);
	uniform_real_distribution<> disReal(0, 1);
	
	long long residue1 = CalculateResidue(iNumberSet, rand_Solution1);
	//S" initialized with S
	vector<int> rand_Solution3(rand_Solution1);
	long long residue3 = residue1;

	int idx1 = 0, idx2 = 0, swapToss = 0;
	//S' - a random neighbor of S
	vector<int> rand_Solution2;
	rand_Solution2 = rand_Solution1;
	//TO CHECK: if the distribution needs to be inside the loop
	uniform_int_distribution<int> disRep(0, setSize - 1);
	for (int i = 0; i < MAX_ITER; i++) {
		idx1 = disRep(gen);
		idx2 = disRep(gen);
		while (idx1 == idx2)
			idx2 = disRep(gen);
		
		rand_Solution2[idx1] *= -1;
		swapToss = dis(gen);
		if (swapToss == 0)
			rand_Solution2[idx2] *= -1;
		long long residue2 = CalculateResidue(iNumberSet, rand_Solution2);
		if (residue2 < residue1) {
			rand_Solution1 = rand_Solution2;
			residue1 = residue2;
		}
		else {
			double prob = disReal(gen);
			long long diffResidue = residue2 - residue1;
			double coolingSchedule = T(i+1);
			double coolingProb = exp(-(diffResidue)/coolingSchedule);
			if ((coolingProb - prob) > 1.0e-12) {
				rand_Solution1 = rand_Solution2;
				residue1 = residue2;
			}
		}
		if (residue1 < residue3) {
			rand_Solution3 = rand_Solution1;
			residue3 = residue1;
		}
	}

	tEndSimulatedAnnealing1 = steady_clock::now();
	duration<double, std::milli> diffSimulatedAnnealing1 = tEndSimulatedAnnealing1 - tStartSimulatedAnnealing1;
	oResidues[3] += residue3;
	oTotalTime[3] += diffSimulatedAnnealing1.count();
	oResultfile << "Simulated Annealing Rep1: " << residue3 << " Time: " << diffSimulatedAnnealing1.count() << endl;
}

//Simulated Annealing Rep 2
void RunAndAnalyzeSimulatedAnnealingRep2(vector<long long>& iNumberSet, vector<int>& iInitialSolution, ofstream& oResultfile,
							   vector<long long>& oResidues, vector<double> oTotalTime) {
	chrono::steady_clock::time_point tStartSimulatedAnnealing1;
	chrono::steady_clock::time_point tEndSimulatedAnnealing1;
	tStartSimulatedAnnealing1 = steady_clock::now();
	//Generate a random solution S of +1 and -1 values
	int setSize = iNumberSet.size();
	vector<int> rand_Solution1(iInitialSolution);
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> disReal(0, 1);
	
	long long residue1 = CalculateResidueForPrePartition(iNumberSet, rand_Solution1);
	//S" initialized with S
	vector<int> rand_Solution3(rand_Solution1);
	long long residue3 = residue1;

	int idx1 = 0, idx2 = 0;
	//S' - a random neighbor of S
	vector<int> rand_Solution2;
	rand_Solution2 = rand_Solution1;
	//TO CHECK: if the distribution needs to be inside the loop
	uniform_int_distribution<int> disRep(0, setSize - 1);
	for (int i = 0; i < MAX_ITER; i++) {
		idx1 = disRep(gen);
		idx2 = disRep(gen);
		while (rand_Solution1[idx1] == idx2)
			idx2 = disRep(gen);

		rand_Solution2[idx1] = idx2;
		long long residue2 = CalculateResidueForPrePartition(iNumberSet, rand_Solution2);

		if (residue2 < residue1) {
			rand_Solution1 = rand_Solution2;
			residue1 = residue2;
		}
		else {
			double prob = disReal(gen);
			long long diffResidue = residue2 - residue1;
			double coolingSchedule = T(i + 1);
			double coolingProb = exp(-(diffResidue) / coolingSchedule);
			if ((coolingProb - prob) > 1.0e-12) {
				rand_Solution1 = rand_Solution2;
				residue1 = residue2;
			}
		}
		if (residue1 < residue3) {
			rand_Solution3 = rand_Solution1;
			residue3 = residue1;
		}
	}

	tEndSimulatedAnnealing1 = steady_clock::now();
	duration<double, std::milli> diffSimulatedAnnealing1 = tEndSimulatedAnnealing1 - tStartSimulatedAnnealing1;
	oResidues[6] += residue3;
	oTotalTime[6] += diffSimulatedAnnealing1.count();
	oResultfile << "Simulated Annealing Rep2: " << residue3 << " Time: " << diffSimulatedAnnealing1.count() << endl;
}

int main(int argc, char** argv) {
	(void)argc;
	//STEP1-Run Karmarkar-Karp for one set of 100 random integers.
	//Read the input file for matrix elements.
	ifstream file;
	file.open(argv[1]);
	if (!file.is_open()) return -1;
	vector<long long> elements;
	long long elementInput;
	while (file >> elementInput){
		elements.push_back(elementInput);
	}
	int setSize = elements.size();
	long long residue = KarmarkarKarpImplementation(elements);
	cout << "Karmarkar Karp Residue: " << residue << endl;
	file.close();

	//Result.txt file for storing the timing and result of algorithm run.
	ofstream result_file("Result.txt", std::ofstream::trunc);
	if (!result_file.is_open()) {
		cout << "Output file could not be opened! Terminating!" << endl;
		return -1;
	}

	//STEP2- Running heuristics on 100 data sets.
	//Generate 100 random instances of the problem
	for (int i = 1; i <= setSize; i++) {
		string fileName = "input" + to_string(i);
		generate_random_ints(fileName);
	}
	vector<long long> Residues(7, 0);
	vector<double> Total_Time(7, 0);
	//Running different heuristics on 100 datasets.
	for (int i = 1; i <= setSize; i++) {
		string fileName = "input" + to_string(i);
		file.open(fileName);
		if (!file.is_open()) return -1;
		vector<long long> numberSet;
		long long num;
		while (file >> num) {
			numberSet.push_back(num);
		}

		RunAndAnalyzeKarmarkarKarpRep1(numberSet, result_file, Residues, Total_Time);

		vector<int> rand_SolutionRep1;
		random_device rd;
		mt19937 gen(rd());
		uniform_int_distribution<int> dis(0, 1);
		for (int n = 0; n < setSize; ++n) {
			int val = dis(gen);
			if (val == 0)
				rand_SolutionRep1.push_back(-1);
			else
				rand_SolutionRep1.push_back(val);
		}
		RunAndAnalyzeRepeatedRandomRep1(numberSet, rand_SolutionRep1, result_file, Residues, Total_Time);
		RunAndAnalyzeHillClimbingRep1(numberSet, rand_SolutionRep1, result_file, Residues, Total_Time);
		RunAndAnalyzeSimulatedAnnealingRep1(numberSet, rand_SolutionRep1, result_file, Residues, Total_Time);

		vector<int> rand_SolutionRep2;
		uniform_int_distribution<int> disRep2(0, setSize - 1);
		for (int n = 0; n < setSize; ++n) {
			int val = disRep2(gen);
			rand_SolutionRep2.push_back(val);
		}
		//RunAndAnalyzeKarmarkarKarpRep2(numberSet, rand_SolutionRep2, result_file);
		RunAndAnalyzeRepeatedRandomRep2(numberSet, rand_SolutionRep2, result_file, Residues, Total_Time);
		RunAndAnalyzeHillClimbingRep2(numberSet, rand_SolutionRep2, result_file, Residues, Total_Time);
		RunAndAnalyzeSimulatedAnnealingRep2(numberSet, rand_SolutionRep2, result_file, Residues, Total_Time);

		result_file << "**********************" << endl;

		file.close();
	}
	//Calculate Averages
	for (int k = 0; k < (int)Residues.size(); k++) {
		Residues[k] = Residues[k] / 7;
		Total_Time[k] = Total_Time[k] / 7;
	}

	result_file << "Karmarkar-Karp Average Residue: " << Residues[0] << "Time Average: " << Total_Time[0] << endl;
	result_file << "Repeated Random Rep1 Residue: " << Residues[1] << "Time Average: " << Total_Time[1] << endl;
	result_file << "Hill Climbing Rep1 Average Residue: " << Residues[2] << "Time Average: " << Total_Time[2] << endl;
	result_file << "Simulated Annealing Rep1 Average Residue: " << Residues[3] << "Time Average: " << Total_Time[3] << endl;
	result_file << "Repeated Random Rep2 Average Residue: " << Residues[4] << "Time Average: " << Total_Time[4] << endl;
	result_file << "Hill Climbing Rep2 Average Residue: " << Residues[5] << "Time Average: " << Total_Time[5] << endl;
	result_file << "Simulated Annealing Rep2 Average Residue: " << Residues[6] << "Time Average: " << Total_Time[6] << endl;

	result_file.close();
	return 0;
}
