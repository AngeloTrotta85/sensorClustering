//============================================================================
// Name        : lambdaSensorsClustering.cpp
// Author      : Angelo Trotta
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <stdlib.h>
#include <stdio.h>
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <algorithm>    // std::find
#include <vector>       // std::vector
#include <list>       // std::list
#include <cstdlib>
#include <ctime>

#include <boost/range/irange.hpp>
#include "MyCoord.h"

#define ALGOVERSION 0

using namespace std;
using namespace boost;

const char* COLOR_LIST_10[] = {
		"green",
		"blue",
		"red",
		"gold",
		"magenta",
		"brown",
		"black",
		"darkorange",
		"salmon",
		"greenyellow"
};


double weightFunciont (MyCoord p1, MyCoord p2) {
	double ris;
	//cout << "weightFunciont - p1: " << p1 << "; p2: " << p2 << endl; fflush(stdout);

	ris = ( 1.0 / pow((p1.distance(p2)), 2.0) );
	//ris = ( 1.0 / p1.distance(p2) );
	//ris = ( 1.0 / pow((p1.distance(p2)), 0.5) );
	//ris = ( 1.0 / pow((p1.distance(p2)), 0.2) );

	if (ris > 1) ris = 1;
	//ris = ( p1.distance(p2) );

	return ris;
}

class InputParser{
public:
	InputParser (int &argc, char **argv){
		for (int i=1; i < argc; ++i)
			this->tokens.push_back(std::string(argv[i]));
	}
	const std::string& getCmdOption(const std::string &option) const{
		std::vector<std::string>::const_iterator itr;
		itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
		if (itr != this->tokens.end() && ++itr != this->tokens.end()){
			return *itr;
		}
		static const std::string empty_string("");
		return empty_string;
	}
	bool cmdOptionExists(const std::string &option) const{
		return std::find(this->tokens.begin(), this->tokens.end(), option)
		!= this->tokens.end();
	}
private:
	std::vector <std::string> tokens;
};

class CoordCluster{
public:
	CoordCluster (){ }
	CoordCluster (std::list<MyCoord> &pl){
		//for (std::list<MyCoord>::iterator it = pl.begin(); it != pl.end(); it++) {
		for (auto& cc : pl) {
			MyCoord newPos = MyCoord(cc.x, cc.y);
			pointsList.push_back(newPos);
		}
	}

	void printCluter(void) {
		for (auto it1 = pointsList.begin(); it1 != pointsList.end(); it1++){
			bool ch = (*it1 == clusterHead);
			if (ch) {
				cout << "*";
			}
			cout << *it1;
			if (ch) {
				cout << "*";
			}
			cout << " ";
		}
	}

	void chooseFirstClusterHead(void){
		if (pointsList.size() > 0) {
			int ris = std::rand() % pointsList.size();
			int idx = 0;
			for (auto& vt : pointsList){
				if (idx == ris) {
					clusterHead.x = vt.x;
					clusterHead.y = vt.y;
				}
				idx++;
			}
		}
	}

	double getMaxCorrelation(void) {
		double ris = 0;
		for (auto it1 = pointsList.begin(); it1 != pointsList.end(); it1++){
			for (auto it2 = it1; it2 != pointsList.end(); it2++){
				if (it1 != it2) {
					double corr = weightFunciont(*it1, *it2);
					if (corr > ris) ris = corr;
				}
			}
		}
		return ris;
	}

	double getAvgCorrelation(void) {
		double ris = 0;
		double count = 0;
		for (auto it1 = pointsList.begin(); it1 != pointsList.end(); it1++){
			for (auto it2 = it1; it2 != pointsList.end(); it2++){
				if (it1 != it2) {
					ris += weightFunciont(*it1, *it2);
					count += 1;
				}
			}
		}
		if (count == 0) {
			return 0;
		}
		else {
			return (ris/count);
		}
	}

	double getMoransI(MyCoord exclude = MyCoord::NIL) {
		MyCoord meanP = MyCoord::ZERO;
		double pCount = 0;
		for (auto& pp : pointsList) {
			if (pp != exclude) {
				meanP += pp;
				pCount += 1;
			}
		}
		meanP = meanP / pCount;

		double sumUp = 0;
		double sumDown = 0;
		double sumW = 0;
		for (auto& pp1 : pointsList) {
			for (auto& pp2 : pointsList) {
				if ((pp1 != exclude) && (pp2 != exclude)) {
					double small_w = 0;
					if (pp1 != pp2) {
						small_w = weightFunciont(pp1, pp2);
					}
					sumW += small_w;

					sumUp += small_w * (pp1.distance(meanP)) * (pp2.distance(meanP));
				}
			}
			sumDown += pow(pp1.distance(meanP), 2.0);
		}

		return ( (pCount / sumW) * (sumUp / sumDown) );
	}

public:
	std::list<MyCoord> pointsList;
	MyCoord clusterHead;
};

double getClusterMaxCorrelation(CoordCluster &cc, MyCoord coord, bool inclusive) {
	double ris = 0;
	for (auto it1 = cc.pointsList.begin(); it1 != cc.pointsList.end(); it1++){
		for (auto it2 = it1; it2 != cc.pointsList.end(); it2++){
			if (it1 != it2) {
				if ((*it1 != coord) && (*it2 != coord)) {
					double corr = weightFunciont(*it1, *it2);
					if (corr > ris) ris = corr;
				}
			}
		}
	}
	if (inclusive) {
		for (auto it1 = cc.pointsList.begin(); it1 != cc.pointsList.end(); it1++){
			if (*it1 != coord) {
				double corr = weightFunciont(*it1, coord);
				if (corr > ris) ris = corr;
			}
		}
	}
	return ris;
}

double weightCoordClusterV1(MyCoord coord, CoordCluster &cc) {
	double noIncl = getClusterMaxCorrelation(cc, coord, false);
	double incl = getClusterMaxCorrelation(cc, coord, true);

	//cout << "Including me; " << incl << "; excluding me: " << noIncl << endl;

	//return noIncl - incl;
	return incl - noIncl;
}

double weightCoordClusterV0(MyCoord coord, CoordCluster &cc) {
	double ris = 0;

	for (auto& pc : cc.pointsList) {
		if (coord != pc) {
			double actW = weightFunciont(coord, pc);
			if (actW > ris) ris = actW;
		}
	}

	return ris;
}

double weightCoordClusterV2(MyCoord coord, std::vector<CoordCluster> &cv, int index) {
	double ris = 0;

	for (int j = 0; j < (int)cv.size(); j++){
	//for (auto& cc : cv) {
		for (auto& pc : cv [j].pointsList) {
			if (coord != pc) {
				double actW = weightFunciont(coord, pc);
				if (actW > ris) ris = actW;
			}
		}
	}

	return ris;
}

//double weightCoordCluster(MyCoord coord, CoordCluster &cc, int version) {
double weightCoordCluster(MyCoord coord, std::vector<CoordCluster> &cv, int index, int version) {
	double ris = 0;

	switch (version) {
	default:
	case 0:
		ris = weightCoordClusterV0(coord, cv[index]);
		break;
	case 1:
		ris = weightCoordClusterV1(coord, cv[index]);
		break;
	case 2:
		ris = weightCoordClusterV2(coord, cv, index);
		break;

	}

	return ris;
}

void generateOutputRandomPoints(std::string fn, int ni, int ss) {
	std::ofstream fout(fn, std::ofstream::out);

	if (fout.is_open()) {

		for (int i : boost::irange(0,ni)) { // i goes from 0 to (ni-1) inclusive
			if (i >= 0) {
				fout << (std::rand() % ss) << ";" << (std::rand() % ss) << endl;
			}
		}

		fout.close();
	}

}


void generateRandomPoints(std::list<MyCoord> &pl, int ss) {
	for (int i : boost::irange(0,100)) { // i goes from 0 to 99 inclusive
		if (i >= 0) {
			MyCoord newCoord = MyCoord(std::rand() % ss, std::rand() % ss);
			newCoord += MyCoord(1,1);
			pl.push_back(newCoord);
			//cout << "Coord: " << i << " --> [" << newCoord.x << ":" << newCoord.y << "]" << endl;
			cout << newCoord.x << ";" << newCoord.y << endl;
		}
	}
}

void importPoints(std::string inputFileName, std::list<MyCoord> &pl) {
	std::ifstream fileInput(inputFileName, std::ifstream::in);
	std::string str;
	double x, y;

	if(fileInput.is_open()) {
		while (std::getline(fileInput, str)) {
			sscanf(str.c_str(), "%lf;%lf", &x, &y);

			pl.push_back(MyCoord(x, y));

			//cout << "From file: " << str << ". Parsed: " << x << ";" << y << endl;
		}

		fileInput.close();
	}
}


void getAndRemoveMaxTillL (int l, MyCoord &vi, std::vector<CoordCluster> &cv) {
	double maxW = -1;
	int clusterMaxIndex = -1;

	for (int j=0; j<=l; j++){
		for (auto it = cv[j].pointsList.begin(); it != cv[j].pointsList.end(); it++){
			MyCoord vt = *it;
			if (vt != cv[j].clusterHead) {
				//double actW = weightFunciont(vt, cv[j].clusterHead);
				//double actW = weightCoordCluster(vt, cv[j], ALGOVERSION);  // ANGELO new
				double actW = weightCoordCluster(vt, cv, j, ALGOVERSION);  // ANGELO new
				if (actW > maxW) {
					vi = vt;
					clusterMaxIndex = j;
				}
			}
		}
	}

	if (clusterMaxIndex < 0) {
		cerr << "Error choosing next ch" << endl;
	}
	else {
		cv[clusterMaxIndex].pointsList.remove(vi);
	}
}

void printAllCluters(std::vector<CoordCluster> &cv, const char *msg) {
	cout << msg << "  CLUSTERS: " << endl;
	for (auto& cc : cv){
		cout << "||| ";
		cc.printCluter();
		cout << endl;
	}
}

double getSystemMaxCorrelation(std::vector<CoordCluster> &clustersVec) {
	double maxCorrelation = 0;
	for (auto& cv : clustersVec) {
		double actMaxCorrelation = cv.getMaxCorrelation();

		if (actMaxCorrelation > maxCorrelation) maxCorrelation = actMaxCorrelation;

	}
	return maxCorrelation;
}

double getSystemAvgCorrelation(std::vector<CoordCluster> &clustersVec) {
	double  sumElements = 0;
	double maxAvgCorrelation = 0;
	for (auto& cv : clustersVec) {
		double actAvgCorrelation = cv.getAvgCorrelation();

		maxAvgCorrelation += actAvgCorrelation;

		sumElements += 1;
	}
	return (maxAvgCorrelation / sumElements);
}

void k_2MMalgo(std::vector<CoordCluster> &cv, int k, int iterations) {

	//printAllCluters(cv, "BEGIN");

	for (int l=0; l<(k-1); l++){
		MyCoord vi;
		getAndRemoveMaxTillL(l, vi, cv);

		cv[l+1].pointsList.push_back(vi);
		cv[l+1].clusterHead = vi;

		//printAllCluters(cv, "New CH");

		//cout << "k_2MMalgo - New Cluster head: " << vi << ". round l: " << l << endl; fflush(stdout);

		// add other points
		for (int j=0; j<=l; j++){
			//for (auto& vt : cv[j].pointsList){
			bool erased;
			do {
				erased = false;
				for (auto it = cv[j].pointsList.begin(); it != cv[j].pointsList.end(); it++){
					MyCoord vt = *it;
					if (vt != cv[j].clusterHead) {
					//if (cv[j].pointsList.size() >1) {
						double wHere = weightCoordCluster(vt, cv, j, ALGOVERSION);
						double wThere = weightCoordCluster(vt, cv, l+1, ALGOVERSION);
						//double wHere = weightCoordCluster(vt, cv[j], ALGOVERSION);
						//double wThere = weightCoordCluster(vt, cv[l+1], ALGOVERSION);
						//if ( weightFunciont(vt, cv[j].clusterHead) >= weightFunciont(vt, vi)) {
						if ( wHere >= wThere) {  // ANGELO new
							//move vt to the new cluster
							cv[l+1].pointsList.push_back(MyCoord(vt.x, vt.y));
							//cv[j].pointsList.erase(it);

							cv[j].pointsList.remove(vt);

							erased = true;
							break;
						}
					}
				}
			} while(erased);
		}

		//printAllCluters(cv, "END CICLE");
	}

	//double systemW = getSystemMaxCorrelation(cv);

	for (int round = 0; round < iterations; round++) {

		bool movingSomeone = false;

		//shuffle
		//for (int ss=0; ss<k; ss++){
		//	std::random_shuffle(cv[ss].pointsList.begin(), cv[ss].pointsList.end());
		//}

		//go back
		for (int l=(k-1); l>0; l--){

			for (int j=(k-1); j>=l; j--){
				//for (auto& vt : cv[j].pointsList){
				bool erased;
				do {
					erased = false;
					for (auto it = cv[j].pointsList.begin(); it != cv[j].pointsList.end(); it++){
						MyCoord vt = *it;
						//if (vt != cv[j].clusterHead) {
						if (cv[j].pointsList.size() > 1) {
							//int idxThere = l-1;
							//double wHere = weightCoordCluster(vt, cv, j, ALGOVERSION);
							//double wThere = 1;
							//double wHere = weightCoordCluster(vt, cv[j], ALGOVERSION);
							//double wThere = weightCoordCluster(vt, cv[l-1], ALGOVERSION);
							double wHere = weightCoordCluster(vt, cv, j, ALGOVERSION);
							double wThere = weightCoordCluster(vt, cv, l-1, ALGOVERSION);
							//for (int jj = 0; jj < k; jj++) {
							//	if (jj != j) {
							//		double wJ = weightCoordCluster(vt, cv, jj, ALGOVERSION);
							//		if (wJ < wThere) {
							//			wThere = wJ;
							//			idxThere = jj;
							//		}
							//	}
							//}
							//if ( weightFunciont(vt, cv[j].clusterHead) >= weightFunciont(vt, vi)) {
							if ( wHere >= wThere) {  // ANGELO new
								//move vt to the new cluster
								cv[l-1].pointsList.push_back(MyCoord(vt.x, vt.y));
								//cv[idxThere].pointsList.push_back(MyCoord(vt.x, vt.y));
								//cv[j].pointsList.erase(it);

								//cv[j].pointsList.remove(vt);
								it = cv[j].pointsList.erase(it);

								movingSomeone = true;

								//erased = true;
								//break;
							}
						}
					}
				} while(erased);
			}
		}

		//go ahead
		for (int l=0; l<(k-1); l++){
			for (int j=0; j<=l; j++){
				//for (auto& vt : cv[j].pointsList){
				bool erased;
				do {
					erased = false;
					for (auto it = cv[j].pointsList.begin(); it != cv[j].pointsList.end(); it++){
						MyCoord vt = *it;
						if (cv[j].pointsList.size() > 1) {
						//if (vt != cv[j].clusterHead) {
							double wHere = weightCoordCluster(vt, cv, j, ALGOVERSION);
							double wThere = weightCoordCluster(vt, cv, l+1, ALGOVERSION);
							//int idxThere = l-1;
							//double wHere = weightCoordCluster(vt, cv, j, ALGOVERSION);
							//double wThere = 1;
							//double wHere = weightCoordCluster(vt, cv[j], ALGOVERSION);
							//double wThere = weightCoordCluster(vt, cv[l+1], ALGOVERSION);
							//for (int jj = 0; jj < k; jj++) {
							//	if (jj != j) {
							//		double wJ = weightCoordCluster(vt, cv, jj, ALGOVERSION);
							//		if (wJ < wThere) {
							//			wThere = wJ;
							//			idxThere = jj;
							//		}
							//	}
							//}
							//if ( weightFunciont(vt, cv[j].clusterHead) >= weightFunciont(vt, vi)) {
							if ( wHere >= wThere) {  // ANGELO new
								//move vt to the new cluster
								cv[l+1].pointsList.push_back(MyCoord(vt.x, vt.y));
								//cv[idxThere].pointsList.push_back(MyCoord(vt.x, vt.y));
								//cv[j].pointsList.erase(it);

								//cv[j].pointsList.remove(vt);
								it = cv[j].pointsList.erase(it);

								movingSomeone = true;

								//erased = true;
								//break;
							}
						}
					}
				} while(erased);
			}
		}

		//double actW = getSystemMaxCorrelation(cv);
		//if (actW == systemW) {
		if (!movingSomeone) {
			cout << "Stopping at the " << round << " iteration" << endl;
			break;
		}
		//systemW = actW;
	}
}

void generateDOTfile(std::string outFileName, std::vector<CoordCluster> &clustVec, double pSize){
	std::ofstream fout(outFileName, std::ofstream::out);
	int count = 1;

	if (fout.is_open()) {

		/*for (int i : boost::irange(0,ni)) { // i goes from 0 to (ni-1) inclusive
			if (i >= 0) {
				fout << (std::rand() % 10000) << ";" << (std::rand() % 10000) << endl;
			}
		}*/

		fout << "graph G{" << endl;

		for (int i = 0; i < (int) clustVec.size(); i++) {
			std::string color = std::string(COLOR_LIST_10[i%10]);

			for (auto& p : clustVec[i].pointsList) {
				fout << "S" << count << " [shape=\"point\" color=\"" << color << "\" pos=\""
						<< p.x << "," << p.y << "!\" width=" << pSize << ", height=" << pSize << "]" << endl;

				count++;
			}
		}

		fout << "}" << endl;

		fout.close();
	}
}

int main(int argc, char **argv) {
	std::list<MyCoord> pointsList;
	std::vector<CoordCluster> clustersVec;
	int k = 5;
	int scenarioSize = 100;
	int iterationNum = 10;

	cout << "Begin!!!" << endl;

	std::srand(std::time(0)); // use current time as seed for random generator

	InputParser input(argc, argv);

	const std::string &inputFileName = input.getCmdOption("-f");
	const std::string &outputGenerateFileName = input.getCmdOption("-g");
	const std::string &numberOfItemToGenerate = input.getCmdOption("-n");
	const std::string &dotFileOutput = input.getCmdOption("-o");
	const std::string &kValue = input.getCmdOption("-k");
	const std::string &scenarioMaxVal = input.getCmdOption("-s");
	const std::string &iterationNumber = input.getCmdOption("-i");

	if (!iterationNumber.empty()) {
		iterationNum = atoi(iterationNumber.c_str());
	}
	if (!scenarioMaxVal.empty()) {
		scenarioSize = atoi(scenarioMaxVal.c_str());
	}
	if (!kValue.empty()) {
		k = atoi(kValue.c_str());
	}

	if ( (!outputGenerateFileName.empty()) && (!numberOfItemToGenerate.empty()) ) {
		int nItems = atoi(numberOfItemToGenerate.c_str());
		generateOutputRandomPoints(outputGenerateFileName, nItems, scenarioSize);
	}

	if (inputFileName.empty()) {
		generateRandomPoints(pointsList, scenarioSize);
	}
	else {
		importPoints(inputFileName, pointsList);
	}

	clustersVec.resize(k);

	/*for (auto& cv : clustersVec) {
		cout << cv.pointsList.size() << " " << cv.clusterHead
				<< " - MaxCorrelation: " << cv.getMaxCorrelation()
				<< " - AvgCorrelation: " << cv.getAvgCorrelation()
				<< endl;
	}*/

	clustersVec[0] = CoordCluster(pointsList);
	clustersVec[0].chooseFirstClusterHead();

	//cout << "Start k_2MMalgo" << endl; fflush(stdout);
	k_2MMalgo(clustersVec, k, iterationNum);
	//cout << "End k_2MMalgo" << endl; fflush(stdout);

	int sumElements = 0;
	//double maxCorrelation = 0;
	//double maxAvgCorrelation = 0;
	for (auto& cv : clustersVec) {
		//double actMaxCorrelation = cv.getMaxCorrelation();
		//double actAvgCorrelation = cv.getAvgCorrelation();
		//cout << cv.pointsList.size() << " " << cv.clusterHead << " - MaxCorrelation: " << actMaxCorrelation << " - AvgCorrelation: " << actAvgCorrelation << " - Moran's I: " << cv.getMoransI() << endl;
		cout << cv.pointsList.size() << " " << cv.clusterHead << " - Moran's I: " << cv.getMoransI() << " - MAX corr: " << cv.getMaxCorrelation() << endl;

		//if (actMaxCorrelation > maxCorrelation) maxCorrelation = actMaxCorrelation;
		//if (actAvgCorrelation > maxAvgCorrelation) maxAvgCorrelation = actAvgCorrelation;

		sumElements += cv.pointsList.size();
	}
	cout << "Total elements: " << sumElements << endl;
	cout << "Maximum MAX correlation: " << getSystemMaxCorrelation(clustersVec) << endl;
	cout << "Maximum AVG correlation: " << getSystemAvgCorrelation(clustersVec) << endl;

	if (!dotFileOutput.empty()) {
		generateDOTfile(dotFileOutput, clustersVec, ((double) scenarioSize)/50.0);
	}

	cout << "End!!!" << endl;
	return EXIT_SUCCESS;
}

