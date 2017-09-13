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
#include <boost/math/special_functions/factorials.hpp>
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

	//ris = ( 1.0 / pow((p1.distance(p2)), 2.0) );
	ris = ( 1.0 / p1.distance(p2) );
	//ris = ( 1.0 / pow((p1.distance(p2)), 0.5) );
	//ris = ( 1.0 / pow((p1.distance(p2)), 0.2) );

	if (ris > 1) ris = 1;
	//ris = ( p1.distance(p2) );

	return ris;
}



class ExtCoord : public MyCoord {
public:
	int clustBelonging;
	std::list<ExtCoord *> toAvoid;
};

class NodeArc {
public:
	static bool compare(const NodeArc& first, const NodeArc& second) {
		return (first.w > second.w);
	}

public:
	ExtCoord *p1;
	ExtCoord *p2;
	double w;
};

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
	CoordCluster (){ clusterID = -1; algo3freeze = false;}
	CoordCluster (std::list<MyCoord> &pl){
		//for (std::list<MyCoord>::iterator it = pl.begin(); it != pl.end(); it++) {
		for (auto& cc : pl) {
			MyCoord newPos = MyCoord(cc.x, cc.y);
			pointsList.push_back(newPos);
		}

		clusterID = -1;
		algo3freeze = false;
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
		for (auto& it1 : pointsList) {
		//for (auto it1 = pointsList.begin(); it1 != pointsList.end(); it1++){
			for (auto& it2 : pointsList) {
			//for (auto it2 = it1; it2 != pointsList.end(); it2++){
				if (it1 != it2) {
					//double corr = weightFunciont(*it1, *it2);
					double corr = weightFunciont(it1, it2);
					if (corr > ris) ris = corr;
				}
			}
		}
		return ris;
	}

	/*double getMaxCorrelation(void) {
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
	}*/

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
	//std::list<MyCoord> pointsList;
	std::vector<MyCoord> pointsList;
	MyCoord clusterHead;
	int clusterID;
	bool algo3freeze;
};

double calculateFullMoranIndex(std::vector<CoordCluster> &cv){
	double xMM = (((double) cv.size()) - 1.0) / 2.0;

	double sumUp = 0;
	double sumDown = 0;
	double sumW = 0;
	double pCount = 0;
	//for (auto& cc1 : cv) {
	for (unsigned int i1 = 0; i1 < cv.size(); i1++) {
		CoordCluster *cc1 = &(cv[i1]);
		for (auto& pp1 : cc1->pointsList){

			//for (auto& cc2 : cv) {
			for (unsigned int i2 = 0; i2 < cv.size(); i2++) {
				CoordCluster *cc2 = &(cv[i2]);
				for (auto& pp2 : cc2->pointsList){

					double small_w = 0;
					if (pp1 != pp2) {
						small_w = weightFunciont(pp1, pp2);
					}
					sumW += small_w;

					//sumUp += small_w * (pp1.distance(cc1->clusterHead)) * (pp2.distance(cc2->clusterHead));
					sumUp += small_w * (i1 - xMM) * (i2 - xMM);
				}
			}

			//sumDown += pow(pp1.distance(cc1->clusterHead), 2.0);
			sumDown += pow(i1 - xMM, 2.0);
			pCount += 1.0;
		}
	}

	return ( (pCount / sumW) * (sumUp / sumDown) );
}

double calculateFullMoranIndex_const(std::vector<CoordCluster> &cv){
	double z = 1;

	double sumUp = 0;
	double sumDown = 0;
	double sumW = 0;
	double pCount = 0;
	//for (auto& cc1 : cv) {
	for (unsigned int i1 = 0; i1 < cv.size(); i1++) {
		CoordCluster *cc1 = &(cv[i1]);
		for (auto& pp1 : cc1->pointsList){

			//for (auto& cc2 : cv) {
			for (unsigned int i2 = 0; i2 < cv.size(); i2++) {
				CoordCluster *cc2 = &(cv[i2]);
				for (auto& pp2 : cc2->pointsList){

					double small_w = 0;
					if (pp1 != pp2) {
						small_w = weightFunciont(pp1, pp2);
					}
					sumW += small_w;

					//sumUp += small_w * (pp1.distance(cc1->clusterHead)) * (pp2.distance(cc2->clusterHead));
					sumUp += small_w * (z) * (z);
				}
			}

			//sumDown += pow(pp1.distance(cc1->clusterHead), 2.0);
			sumDown += pow(z, 2.0);
			pCount += 1.0;
		}
	}

	return ( (pCount / sumW) * (sumUp / sumDown) );
}

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

void generateOutputRandomGridPoints(std::string fn, int ni, int ss) {
	std::ofstream fout(fn, std::ofstream::out);

	if (fout.is_open()) {

		int remainingP = ni;

		int sq = sqrt((double) ni);
		if ((sq*sq) < ni) ++sq;

		double iterDist = ((double) ss) / ((double) sq);

		for (int i = 0; i < sq; ++i) {
			for (int j = 0; j < sq; ++j) {
				if (remainingP > 0) {

					fout << (iterDist * i) << ";" << (iterDist * j) << endl;

					--remainingP;
				}
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
	std::vector<MyCoord>::iterator maxIT;

	for (int j=0; j<=l; j++){
		if (cv[j].algo3freeze) continue;
		for (auto it = cv[j].pointsList.begin(); it != cv[j].pointsList.end(); it++){
			MyCoord vt = *it;
			if (vt != cv[j].clusterHead) {
				//double actW = weightFunciont(vt, cv[j].clusterHead);
				//double actW = weightCoordCluster(vt, cv[j], ALGOVERSION);  // ANGELO new
				double actW = weightCoordCluster(vt, cv, j, ALGOVERSION);  // ANGELO new
				if (actW > maxW) {
					vi = vt;
					clusterMaxIndex = j;
					maxIT = it;
				}
			}
		}
	}

	if (clusterMaxIndex < 0) {
		cerr << "Error choosing next ch" << endl;
	}
	else {
		//cv[clusterMaxIndex].pointsList.remove(vi);
		cv[clusterMaxIndex].pointsList.erase(maxIT);
	}
}

void printAllCluters(std::vector<CoordCluster> &cv, const char *msg, bool printPoints) {
	cout << msg << "  CLUSTERS: " << endl;
	for (auto& cc : cv){
		cout << "|| " << cc.clusterID <<
				" | " << cc.algo3freeze <<
				" | " << cc.pointsList.size() <<
				" | " << cc.getMaxCorrelation();

		if (printPoints) {
			cout << " | ";
			cc.printCluter();
		}
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

							//cv[j].pointsList.remove(vt);
							cv[j].pointsList.erase(it);

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
		for (auto& ss : cv){
			//std::random_shuffle(cv[ss].pointsList.begin(), cv[ss].pointsList.end());
			std::random_shuffle(ss.pointsList.begin(), ss.pointsList.end());
		}

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
								cv[j].pointsList.erase(it);
								//it = cv[j].pointsList.erase(it);

								movingSomeone = true;

								erased = true;
								break;
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
								cv[j].pointsList.erase(it);
								//it = cv[j].pointsList.erase(it);

								movingSomeone = true;

								erased = true;
								break;
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
	//printAllCluters(cv, "FINAL CLUSTERS");
}

void k_2MMalgo3(std::vector<CoordCluster> &cv, int k, int iterations) {
	for (int l=0; l<(k-1); l++){
		MyCoord vi;

		if (cv[l].algo3freeze) continue;

		int nextIdx = -1;
		for (int kk = (l+1); kk < k; kk++) {
			if (!cv[kk].algo3freeze) {
				nextIdx = kk;
				break;
			}
		}
		if (nextIdx < 0) continue;


		getAndRemoveMaxTillL(l, vi, cv);

		cv[nextIdx].pointsList.push_back(vi);
		cv[nextIdx].clusterHead = vi;

		//printAllCluters(cv, "New CH");

		//cout << "k_2MMalgo - New Cluster head: " << vi << ". round l: " << l << endl; fflush(stdout);

		// add other points
		for (int j=0; j<=l; j++){
			//for (auto& vt : cv[j].pointsList){
			if (cv[j].algo3freeze) continue;

			bool erased;
			do {
				erased = false;
				for (auto it = cv[j].pointsList.begin(); it != cv[j].pointsList.end(); it++){
					MyCoord vt = *it;
					if (vt != cv[j].clusterHead) {
					//if (cv[j].pointsList.size() >1) {
						double wHere = weightCoordCluster(vt, cv, j, ALGOVERSION);
						double wThere = weightCoordCluster(vt, cv, nextIdx, ALGOVERSION);
						//double wHere = weightCoordCluster(vt, cv[j], ALGOVERSION);
						//double wThere = weightCoordCluster(vt, cv[l+1], ALGOVERSION);
						//if ( weightFunciont(vt, cv[j].clusterHead) >= weightFunciont(vt, vi)) {
						if ( wHere >= wThere) {  // ANGELO new
							//move vt to the new cluster
							cv[nextIdx].pointsList.push_back(MyCoord(vt.x, vt.y));
							//cv[j].pointsList.erase(it);

							//cv[j].pointsList.remove(vt);
							cv[j].pointsList.erase(it);

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
		for (auto& ss : cv){
			//std::random_shuffle(cv[ss].pointsList.begin(), cv[ss].pointsList.end());
			std::random_shuffle(ss.pointsList.begin(), ss.pointsList.end());
		}

		//go back
		for (int l=(k-1); l>0; l--){

			if (cv[l].algo3freeze) continue;

			int nextIdx = -1;
			for (int kk = (l-1); kk >= 0; kk--) {
				if (!cv[kk].algo3freeze) {
					nextIdx = kk;
					break;
				}
			}
			if (nextIdx < 0) continue;

			for (int j=(k-1); j>=l; j--){
				//for (auto& vt : cv[j].pointsList){
				if (cv[j].algo3freeze) continue;

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
							double wThere = weightCoordCluster(vt, cv, nextIdx, ALGOVERSION);
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
								cv[nextIdx].pointsList.push_back(MyCoord(vt.x, vt.y));
								//cv[idxThere].pointsList.push_back(MyCoord(vt.x, vt.y));
								//cv[j].pointsList.erase(it);

								//cv[j].pointsList.remove(vt);
								cv[j].pointsList.erase(it);
								//it = cv[j].pointsList.erase(it);

								movingSomeone = true;

								erased = true;
								break;
							}
						}
					}
				} while(erased);
			}
		}

		//go ahead
		for (int l=0; l<(k-1); l++){
			if (cv[l].algo3freeze) continue;

			int nextIdx = -1;
			for (int kk = (l+1); kk < k; kk++) {
				if (!cv[kk].algo3freeze) {
					nextIdx = kk;
					break;
				}
			}
			if (nextIdx < 0) continue;

			for (int j=0; j<=l; j++){
				if (cv[j].algo3freeze) continue;
				//for (auto& vt : cv[j].pointsList){
				bool erased;
				do {
					erased = false;
					for (auto it = cv[j].pointsList.begin(); it != cv[j].pointsList.end(); it++){
						MyCoord vt = *it;
						if (cv[j].pointsList.size() > 1) {
						//if (vt != cv[j].clusterHead) {
							double wHere = weightCoordCluster(vt, cv, j, ALGOVERSION);
							double wThere = weightCoordCluster(vt, cv, nextIdx, ALGOVERSION);
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
								cv[nextIdx].pointsList.push_back(MyCoord(vt.x, vt.y));
								//cv[idxThere].pointsList.push_back(MyCoord(vt.x, vt.y));
								//cv[j].pointsList.erase(it);

								//cv[j].pointsList.remove(vt);
								cv[j].pointsList.erase(it);
								//it = cv[j].pointsList.erase(it);

								movingSomeone = true;

								erased = true;
								break;
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
	//printAllCluters(cv, "FINAL CLUSTERS");
}

void generateDOTfile(std::string outFileName, std::vector<CoordCluster> &clustVec, double pSize){
	std::ofstream fout(outFileName, std::ofstream::out);
	int count = 1;

	if (fout.is_open()) {

		double maxCorrelation = 0;
		int maxIdx = 0;
		for (int i = 0; i < (int) clustVec.size(); i++) {
		//for (auto& cv : clustVec) {
			double actMaxCorrelation = clustVec[i].getMaxCorrelation();

			if (actMaxCorrelation > maxCorrelation) {
				maxCorrelation = actMaxCorrelation;
				maxIdx = i;
			}

		}
		//return maxCorrelation;

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

				if (i == maxIdx) {
					fout << "S" << count << "_rad [shape=\"circle\" color=\"" << "black" << "\" style=\"dotted\" label=\"\" pos=\""
							<< p.x << "," << p.y << "!\" width=" << (2.0/maxCorrelation) << ", height=" << (2.0/maxCorrelation) << "]" << endl;
				}

				count++;
			}
		}

		fout << "}" << endl;

		fout.close();
	}
}

//void k_2MMalgo(std::vector<CoordCluster> &cv, int k, int iterations) {
void equalize1(std::vector<CoordCluster> &cv, int pointsNumber, int k) {
	int n4cluster = pointsNumber / k;
	//int remainingP = pointsNumber % k;

	//cout << "Equalizing1 the " << k << " clusters; " << n4cluster << " each; remaining points: " << remainingP << endl;

	bool removed;
	do {
		removed = false;
		for (unsigned int idx1 = 0; idx1 < cv.size(); idx1++) {
			if (((int) cv[idx1].pointsList.size()) < n4cluster) {
				double corrMinGlob = std::numeric_limits<double>::max();
				unsigned int clustMin = 0;
				std::vector<MyCoord>::iterator itMinGlob;

				for (unsigned int idx2 = 0; idx2 < cv.size(); idx2++) {
					if (idx1 != idx2) {
						if (((int) cv[idx2].pointsList.size()) > n4cluster) {
							std::vector<MyCoord>::iterator itMinLoc;
							double corrMinLoc = std::numeric_limits<double>::max();

							for (auto itP2 = cv[idx2].pointsList.begin(); itP2 != cv[idx2].pointsList.end(); itP2++) {
								//for (auto& pp : cv[idx2].pointsList) {
								for (auto itP1 = cv[idx1].pointsList.begin(); itP1 != cv[idx1].pointsList.end(); itP1++) {
									double corr = weightFunciont(*itP1, *itP2);
									if (corr < corrMinLoc) {
										corrMinLoc = corr;
										itMinLoc = itP2;
									}
								}
							}

							if (corrMinLoc < corrMinGlob) {
								corrMinGlob = corrMinLoc;
								itMinGlob = itMinLoc;
								clustMin = idx2;
							}
						}
					}
				}

				if (corrMinGlob < std::numeric_limits<double>::max()) {
					cv[idx1].pointsList.push_back(MyCoord(itMinGlob->x, itMinGlob->y));

					cv[clustMin].pointsList.erase(itMinGlob);

					removed = true;
					break;
				}
				else {
					cerr << "Error in equalizing the clusters" << endl;
					exit(EXIT_FAILURE);
				}
			}
		}
	} while (removed);
}

void equalize2(std::vector<CoordCluster> &cv, int pointsNumber, int k) {
	int n4cluster = pointsNumber / k;
	//int remainingP = pointsNumber % k;

	//cout << "Equalizing2 the " << k << " clusters; " << n4cluster << " each; remaining points: " << remainingP << endl;

	bool removed;
	do {

		double corrMinGlob = std::numeric_limits<double>::max();
		unsigned int clustMinSrc = 0;
		unsigned int clustMinDst = 0;
		std::vector<MyCoord>::iterator itMinGlobSrc;

		removed = false;
		for (unsigned int idx1 = 0; idx1 < cv.size(); idx1++) {
			if (((int) cv[idx1].pointsList.size()) < n4cluster) {
				for (unsigned int idx2 = 0; idx2 < cv.size(); idx2++) {
					if (idx1 != idx2) {
						if (((int) cv[idx2].pointsList.size()) > n4cluster) {
							for (auto itP2 = cv[idx2].pointsList.begin(); itP2 != cv[idx2].pointsList.end(); itP2++) {
								//for (auto& pp : cv[idx2].pointsList) {
								for (auto itP1 = cv[idx1].pointsList.begin(); itP1 != cv[idx1].pointsList.end(); itP1++) {
									double corr = weightFunciont(*itP1, *itP2);

									if (corr < corrMinGlob) {
										corrMinGlob = corr;

										clustMinDst = idx1;
										clustMinSrc = idx2;

										itMinGlobSrc = itP2;
									}
								}
							}
						}
					}
				}
			}
		}


		if (corrMinGlob < std::numeric_limits<double>::max()) {
			cv[clustMinDst].pointsList.push_back(MyCoord(itMinGlobSrc->x, itMinGlobSrc->y));

			cv[clustMinSrc].pointsList.erase(itMinGlobSrc);

			removed = true;
		}
		/*else {
			cerr << "Error in equalizing the clusters" << endl;
			exit(EXIT_FAILURE);
		}*/

	} while (removed);
}

void equalize3(std::vector<CoordCluster> &cv, int pointsNumber, int k, unsigned int n4cluster, int remainingP) {
	//int n4cluster = pointsNumber / k;
	//int remainingP = pointsNumber % k;
	unsigned int minClusterSize = n4cluster;
	if (remainingP > 0) {
		++minClusterSize;
	}

	//cout << "Equalizing3 the " << k << " clusters; " << n4cluster << " each; remaining points: " << remainingP << endl;

	unsigned int minSize = std::numeric_limits<unsigned int>::max();
	unsigned int idxActCluster = 0;
	for (unsigned int idxC = 0; idxC < cv.size(); ++idxC) {
		if (!cv[idxC].algo3freeze) {
			if ((cv[idxC].pointsList.size() >= minClusterSize) && (cv[idxC].pointsList.size() < minSize)) {
				idxActCluster = idxC;
				minSize = cv[idxC].pointsList.size();
			}
		}
	}

	if (minSize != std::numeric_limits<int>::max()){
		cv[idxActCluster].algo3freeze = true;

		if (minSize > minClusterSize) {
			unsigned int destClust = 0;

			for (unsigned int idxC = 0; idxC < cv.size(); ++idxC) {
				if (!cv[idxC].algo3freeze) {
					destClust = idxC;
					break;
				}
			}

			do {

				std::vector<MyCoord>::iterator itRemove;
				double minCorrToRemove = std::numeric_limits<double>::max();

				for (auto pp = cv[idxActCluster].pointsList.begin(); pp != cv[idxActCluster].pointsList.end(); pp++) {
					//for (auto& pp : cv[idxActCluster].pointsList) {

					for (auto& cc : cv) {
						if (!cc.algo3freeze) {
							for (auto ppCheck = cc.pointsList.begin(); ppCheck != cc.pointsList.end(); ppCheck++) {
								//for (auto& ppCheck : cc.pointsList) {
								double corr = weightFunciont(*pp, *ppCheck);
								if (corr < minCorrToRemove) {
									minCorrToRemove = corr;
									itRemove = pp;
								}
							}
						}
					}
				}

				if (minCorrToRemove != std::numeric_limits<double>::max()) {
					cv[destClust].pointsList.push_back(MyCoord(itRemove->x, itRemove->y));

					cv[idxActCluster].pointsList.erase(itRemove);
				}
				else {
					cerr << "Error in Equalizing3" << endl;
					exit(EXIT_FAILURE);
				}

			} while(cv[idxActCluster].pointsList.size() > minClusterSize);

		}
	}
}

bool allEqualized3(std::vector<CoordCluster> &cv) {
	bool ris = true;

	for (auto& c : cv) {
		if (!c.algo3freeze) {
			ris = false;
		}
	}

	return ris;
}

int numNotEqualized3(std::vector<CoordCluster> &cv) {
	int ris = 0;

	for (auto& c : cv) {
		if (!c.algo3freeze) {
			ris++;
		}
	}

	return ris;
}

void randomizeClusters(std::vector<CoordCluster> &cv, unsigned int nCluster, unsigned int n4cluster, int remainingP) {
	for (unsigned int i = 1; i < nCluster; ++i) {
		for (unsigned int j = 0; j < n4cluster; ++j) {

			unsigned int r = rand() % cv[0].pointsList.size();

			cv[i].pointsList.push_back(MyCoord(cv[0].pointsList[r].x, cv[0].pointsList[r].y));

			cv[0].pointsList.erase (cv[0].pointsList.begin()+r);
		}
	}
}

bool isBetterChange(double corr2beat, std::vector<MyCoord>::iterator &itP1, std::vector<MyCoord> &pl1, std::vector<MyCoord>::iterator &itP2, std::vector<MyCoord> &pl2) {
	bool ris = false;

	double corr1Here, corr1There, corr2Here, corr2There;

	corr1Here = corr1There = corr2Here = corr2There = 0;

	for (auto it = pl1.begin(); it != pl1.end(); it++) {
		if (it != itP1) {
			double corr1 = weightFunciont(*itP1, *it);
			if (corr1 > corr1Here) {
				corr1Here = corr1;
			}

			double corr2 = weightFunciont(*itP2, *it);
			if (corr2 > corr2There) {
				corr2There = corr2;
			}
		}
	}

	for (auto it = pl2.begin(); it != pl2.end(); it++) {
		if (it != itP2) {
			double corr2 = weightFunciont(*itP2, *it);
			if (corr2 > corr2Here) {
				corr2Here = corr2;
			}

			double corr1 = weightFunciont(*itP1, *it);
			if (corr1 > corr1There) {
				corr1There = corr1;
			}
		}
	}

	//if (corr1There < corr2beat) {
	//if (corr2There < corr2beat) {
	if ( (corr1There < corr2beat) && (corr2There < corr2beat) ) {
		ris = true;
	}

	return ris;
}

bool isBetterChange3CW(double corr2beat, std::vector<MyCoord>::iterator &itP1, std::vector<MyCoord> &pl1, std::vector<MyCoord>::iterator &itP2, std::vector<MyCoord> &pl2, std::vector<MyCoord>::iterator &itP3, std::vector<MyCoord> &pl3) {
	bool ris = false;
	double corr1Here, corr1There, corr2Here, corr2There, corr3Here, corr3There;

	corr1Here = corr1There = corr2Here = corr2There = corr3Here = corr3There = 0;

	for (auto it = pl1.begin(); it != pl1.end(); it++) {
		if (it != itP1) {
			double corr1 = weightFunciont(*itP1, *it);
			if (corr1 > corr1Here) {
				corr1Here = corr1;
			}

			double corr3 = weightFunciont(*itP3, *it);
			if (corr3 > corr3There) {
				corr3There = corr3;
			}
		}
	}

	for (auto it = pl2.begin(); it != pl2.end(); it++) {
		if (it != itP2) {
			double corr2 = weightFunciont(*itP2, *it);
			if (corr2 > corr2Here) {
				corr2Here = corr2;
			}

			double corr1 = weightFunciont(*itP1, *it);
			if (corr1 > corr1There) {
				corr1There = corr1;
			}
		}
	}

	for (auto it = pl3.begin(); it != pl3.end(); it++) {
		if (it != itP3) {
			double corr3 = weightFunciont(*itP3, *it);
			if (corr3 > corr3Here) {
				corr3Here = corr3;
			}

			double corr2 = weightFunciont(*itP2, *it);
			if (corr2 > corr2There) {
				corr2There = corr2;
			}
		}
	}


	if ( (corr1There < corr2beat) && (corr2There < corr2beat) && (corr3There < corr2beat) ) {
		ris = true;
	}

	return ris;
}

bool isBetterChange3CCW(double corr2beat, std::vector<MyCoord>::iterator &itP1, std::vector<MyCoord> &pl1, std::vector<MyCoord>::iterator &itP2, std::vector<MyCoord> &pl2, std::vector<MyCoord>::iterator &itP3, std::vector<MyCoord> &pl3) {
	bool ris = false;
	double corr1Here, corr1There, corr2Here, corr2There, corr3Here, corr3There;

	corr1Here = corr1There = corr2Here = corr2There = corr3Here = corr3There = 0;

	for (auto it = pl1.begin(); it != pl1.end(); it++) {
		if (it != itP1) {
			double corr1 = weightFunciont(*itP1, *it);
			if (corr1 > corr1Here) {
				corr1Here = corr1;
			}

			double corr2 = weightFunciont(*itP2, *it);
			if (corr2 > corr2There) {
				corr2There = corr2;
			}
		}
	}

	for (auto it = pl2.begin(); it != pl2.end(); it++) {
		if (it != itP2) {
			double corr2 = weightFunciont(*itP2, *it);
			if (corr2 > corr2Here) {
				corr2Here = corr2;
			}

			double corr3 = weightFunciont(*itP3, *it);
			if (corr3 > corr3There) {
				corr3There = corr3;
			}
		}
	}

	for (auto it = pl3.begin(); it != pl3.end(); it++) {
		if (it != itP3) {
			double corr3 = weightFunciont(*itP3, *it);
			if (corr3 > corr3Here) {
				corr3Here = corr3;
			}

			double corr1 = weightFunciont(*itP1, *it);
			if (corr1 > corr1There) {
				corr1There = corr1;
			}
		}
	}


	if ( (corr1There < corr2beat) && (corr2There < corr2beat) && (corr3There < corr2beat) ) {
		ris = true;
	}

	return ris;
}

void makeSwaps(std::vector<CoordCluster> &cv, unsigned int k) {

	int swapNumber = 0;
	bool madeSwap;

	do {
		unsigned int worstCluster = 0;
		double corrMax = 0;
		std::vector<MyCoord>::iterator itMinSrc1, itMinSrc2;

		madeSwap = false;

		/*
		 * for (unsigned int idx1 = 0; idx1 < cv.size(); idx1++) {
			if (((int) cv[idx1].pointsList.size()) < n4cluster) {
				for (unsigned int idx2 = 0; idx2 < cv.size(); idx2++) {
					if (idx1 != idx2) {
						if (((int) cv[idx2].pointsList.size()) > n4cluster) {
							for (auto itP2 = cv[idx2].pointsList.begin(); itP2 != cv[idx2].pointsList.end(); itP2++) {
								//for (auto& pp : cv[idx2].pointsList) {
								for (auto itP1 = cv[idx1].pointsList.begin(); itP1 != cv[idx1].pointsList.end(); itP1++) {
									double corr = weightFunciont(*itP1, *itP2);
		 */

		for (unsigned int idx1 = 0; idx1 < cv.size(); idx1++) {
			for (auto itP1 = cv[idx1].pointsList.begin(); itP1 != cv[idx1].pointsList.end(); itP1++) {
				for (auto itP2 = cv[idx1].pointsList.begin(); itP2 != cv[idx1].pointsList.end(); itP2++) {
					if (itP1 != itP2) {
						double corr = weightFunciont(*itP1, *itP2);

						//cout << "Check correlation: " << corr << endl;

						if (corr > corrMax) {
							corrMax = corr;

							itMinSrc1 = itP1;
							itMinSrc2 = itP2;
							worstCluster = idx1;
						}
					}
				}
			}
		}

		//cout << "Worst correlation: [" << itMinSrc1->info() << "] - ["
		//		<< itMinSrc2->info() << "] in cluster "  << worstCluster << endl;

		for (unsigned int idx1 = 0; idx1 < cv.size(); idx1++) {
			if (idx1 != worstCluster) {

				for (auto itP1 = cv[idx1].pointsList.begin(); itP1 != cv[idx1].pointsList.end(); itP1++) {
					if (isBetterChange(corrMax, itP1, cv[idx1].pointsList, itMinSrc1, cv[worstCluster].pointsList)) {

						MyCoord newCoordInHere = MyCoord(itMinSrc1->x, itMinSrc1->y);
						MyCoord newCoordInWorse = MyCoord(itP1->x, itP1->y);

						cv[worstCluster].pointsList.erase (itMinSrc1);
						cv[idx1].pointsList.erase (itP1);

						cv[idx1].pointsList.push_back(newCoordInHere);
						cv[worstCluster].pointsList.push_back(newCoordInWorse);

						//cout << "Making swap" << endl;

						madeSwap = true;
						++swapNumber;
						break;

					} else if (isBetterChange(corrMax, itP1, cv[idx1].pointsList, itMinSrc2, cv[worstCluster].pointsList)) {

						MyCoord newCoordInHere = MyCoord(itMinSrc2->x, itMinSrc2->y);
						MyCoord newCoordInWorse = MyCoord(itP1->x, itP1->y);

						cv[worstCluster].pointsList.erase (itMinSrc2);
						cv[idx1].pointsList.erase (itP1);

						cv[idx1].pointsList.push_back(newCoordInHere);
						cv[worstCluster].pointsList.push_back(newCoordInWorse);

						//cout << "Making swap" << endl;

						madeSwap = true;
						++swapNumber;
						break;
					}
				}
			}

			if (madeSwap) break;
		}

		if (!madeSwap) {
			//try a triangular swap
			for (unsigned int idx1 = 0; idx1 < cv.size(); idx1++) {
				for (unsigned int idx2 = 0; idx2 < cv.size(); idx2++) {
					if ((idx1 != worstCluster) && (idx2 != worstCluster) && (idx1 != idx2)) {

						for (auto itP1 = cv[idx1].pointsList.begin(); itP1 != cv[idx1].pointsList.end(); itP1++) {
							for (auto itP2 = cv[idx2].pointsList.begin(); itP2 != cv[idx2].pointsList.end(); itP2++) {
								// ClockWise p1-p2-Min1 --> Min1-p1-p2
								// CounterClockWise p1-p2-Min1 --> p2-Min1-p1
								if (isBetterChange3CW(corrMax, itP1, cv[idx1].pointsList, itP2, cv[idx2].pointsList, itMinSrc1, cv[worstCluster].pointsList)) {

									MyCoord newCoordP1 = MyCoord(itP1->x, itP1->y);
									MyCoord newCoordP2 = MyCoord(itP2->x, itP2->y);
									MyCoord newCoordWorse = MyCoord(itMinSrc1->x, itMinSrc1->y);

									cv[worstCluster].pointsList.erase (itMinSrc1);
									cv[idx1].pointsList.erase (itP1);
									cv[idx2].pointsList.erase (itP2);

									cv[idx1].pointsList.push_back(newCoordWorse);
									cv[idx2].pointsList.push_back(newCoordP1);
									cv[worstCluster].pointsList.push_back(newCoordP2);

									madeSwap = true;
									cout << "Making CW1 swap" << endl;
									++swapNumber;
									break;
								} else	if (isBetterChange3CCW(corrMax, itP1, cv[idx1].pointsList, itP2, cv[idx2].pointsList, itMinSrc1, cv[worstCluster].pointsList)) {

									MyCoord newCoordP1 = MyCoord(itP1->x, itP1->y);
									MyCoord newCoordP2 = MyCoord(itP2->x, itP2->y);
									MyCoord newCoordWorse = MyCoord(itMinSrc1->x, itMinSrc1->y);

									cv[worstCluster].pointsList.erase (itMinSrc1);
									cv[idx1].pointsList.erase (itP1);
									cv[idx2].pointsList.erase (itP2);

									cv[idx1].pointsList.push_back(newCoordP2);
									cv[idx2].pointsList.push_back(newCoordWorse);
									cv[worstCluster].pointsList.push_back(newCoordP1);

									madeSwap = true;
									cout << "Making CCW1 swap" << endl;
									++swapNumber;
									break;
								} else if (isBetterChange3CW(corrMax, itP1, cv[idx1].pointsList, itP2, cv[idx2].pointsList, itMinSrc2, cv[worstCluster].pointsList)) {

									MyCoord newCoordP1 = MyCoord(itP1->x, itP1->y);
									MyCoord newCoordP2 = MyCoord(itP2->x, itP2->y);
									MyCoord newCoordWorse = MyCoord(itMinSrc2->x, itMinSrc2->y);

									cv[worstCluster].pointsList.erase (itMinSrc2);
									cv[idx1].pointsList.erase (itP1);
									cv[idx2].pointsList.erase (itP2);

									cv[idx1].pointsList.push_back(newCoordWorse);
									cv[idx2].pointsList.push_back(newCoordP1);
									cv[worstCluster].pointsList.push_back(newCoordP2);

									madeSwap = true;
									cout << "Making CW2 swap" << endl;
									++swapNumber;
									break;
								} else if (isBetterChange3CCW(corrMax, itP1, cv[idx1].pointsList, itP2, cv[idx2].pointsList, itMinSrc2, cv[worstCluster].pointsList)) {

									MyCoord newCoordP1 = MyCoord(itP1->x, itP1->y);
									MyCoord newCoordP2 = MyCoord(itP2->x, itP2->y);
									MyCoord newCoordWorse = MyCoord(itMinSrc2->x, itMinSrc2->y);

									cv[worstCluster].pointsList.erase (itMinSrc2);
									cv[idx1].pointsList.erase (itP1);
									cv[idx2].pointsList.erase (itP2);

									cv[idx1].pointsList.push_back(newCoordP2);
									cv[idx2].pointsList.push_back(newCoordWorse);
									cv[worstCluster].pointsList.push_back(newCoordP1);

									madeSwap = true;
									cout << "Making CCW2 swap" << endl;
									++swapNumber;
									break;
								}
							}
							if (madeSwap) break;
						}
						if (madeSwap) break;
					}
				}
				if (madeSwap) break;
			}
		}

	} while (madeSwap);

	cout << "Total swaps: " << swapNumber << endl;
}

void optAlgo(std::list<MyCoord> &pl, std::vector<CoordCluster> &cv, unsigned int n4cPlus, unsigned int n4c, unsigned int k) {
	std::vector<unsigned int> idxs(pl.size(), 0);
	std::vector<unsigned int> checkVec(k, 0);
	//bool finish = false;
	//unsigned int nRounds = pow((double)pl.size(), (double)k);
	unsigned long long int nRounds = pow((double)k, (double)pl.size());
	unsigned long long int nRealTest = 0;
	double bestCorr = std::numeric_limits<double>::max();
	std::vector<CoordCluster> bestCombination;
	bestCombination.resize(k);

	checkVec[0] = pl.size();

	cout << "Start opt algorithm with " << k << "^" << pl.size() << " = " << nRounds << " operations" << endl;
	cout << "Opt algorithm with factorial " << boost::math::factorial<long double>((long double)pl.size())  << endl;

	for (unsigned long long int nr = 0; nr < nRounds; ++nr) {
	//do {
		//fprintf(stdout, "Optimal algorithm %lf%%\r", (((double) (nr + 1)) / ((double) nRounds)) * 100.0);// fflush(stdout);

		bool formationOK = true;
		for (unsigned int i_check = 0; i_check < k; ++i_check){
			if (	((i_check == 0) && (checkVec[i_check] != n4cPlus)) ||
					((i_check > 0) && (checkVec[i_check] != n4c)) ) {
				formationOK = false;
				break;
			}
		}

		if(formationOK) {
			++nRealTest;
			for (auto& cc : cv) {
				cc.pointsList.clear();
			}
			unsigned int idxPL = 0;
			for (auto& p : pl) {
				MyCoord newPos = MyCoord(p.x, p.y);
				cv[idxs[idxPL]].pointsList.push_back(newPos);
				++idxPL;
			}

			double worstCorr = 0;
			for (auto& cc : cv) {
				double actCorr = cc.getMaxCorrelation();
				if (actCorr > worstCorr) {
					worstCorr = actCorr;
				}
			}
			if (worstCorr < bestCorr) {
				bestCorr = worstCorr;

				for (auto& cc : bestCombination) {
					cc.pointsList.clear();
				}
				std::vector<CoordCluster>::iterator itcBest = bestCombination.begin();
				for (std::vector<CoordCluster>::iterator itc = cv.begin(); itc != cv.end(); itc++) {

					for (auto& pp : itc->pointsList) {
						MyCoord newPos = MyCoord(pp.x, pp.y);
						itcBest->pointsList.push_back(newPos);
					}
					itcBest->clusterID = itc->clusterID;

					++itcBest;
				}
			}
		}


		/*
		for (auto& cc : cv) {
			cc.pointsList.clear();
		}
		unsigned int idxPL = 0;
		for (auto& p : pl) {
			MyCoord newPos = MyCoord(p.x, p.y);
			cv[idxs[idxPL]].pointsList.push_back(newPos);
			++idxPL;
		}



		// check formation
		//bool formationOK = true;
		for (unsigned iform = 0; iform < cv.size(); ++iform){

			if (	((iform == 0) && (cv[iform].pointsList.size() != n4cPlus)) ||
					((iform > 0) && (cv[iform].pointsList.size() != n4c)) ) {
				formationOK = false;
				break;
			}
		}
		if(formationOK) {
			double worstCorr = 0;
			for (auto& cc : cv) {
				double actCorr = cc.getMaxCorrelation();
				if (actCorr > worstCorr) {
					worstCorr = actCorr;
				}
			}
			if (worstCorr < bestCorr) {
				bestCorr = worstCorr;

				for (auto& cc : bestCombination) {
					cc.pointsList.clear();
				}
				std::vector<CoordCluster>::iterator itcBest = bestCombination.begin();
				for (std::vector<CoordCluster>::iterator itc = cv.begin(); itc != cv.end(); itc++) {

					for (auto& pp : itc->pointsList) {
						MyCoord newPos = MyCoord(pp.x, pp.y);
						itcBest->pointsList.push_back(newPos);
					}
					itcBest->clusterID = itc->clusterID;

					++itcBest;
				}
			}
		}
		*/



		//increase the idxs counter
		for (unsigned int idxCheck = 0; idxCheck < pl.size(); ++idxCheck) {

			--checkVec[idxs[idxCheck]];
			++idxs[idxCheck];
			if (idxs[idxCheck] < k) {
				++checkVec[idxs[idxCheck]];
				break;
			}
			else {
				idxs[idxCheck] = 0;
				++checkVec[idxs[idxCheck]];
			}
		}

		//cout << "Vector clusters: ";
		//for (auto& ii : checkVec) {
		//	cout << ii << " ";
		//}
		//cout << endl;
		//cout << "Vector points: ";
		//for (auto& ii : idxs) {
		//	cout << ii << " ";
		//}
		//cout << endl << endl;

		// check if every combination is done
		//finish = true;
		//for (auto& idd : idxs) {
		//	if (idd < (k-1)) {
		//		finish = false;
		//		break;
		//	}
		//}
	//} while (finish);
	}

	for (auto& cc : cv) {
		cc.pointsList.clear();
	}
	std::vector<CoordCluster>::iterator itc = cv.begin();
	for (std::vector<CoordCluster>::iterator itcBest = bestCombination.begin(); itcBest != bestCombination.end(); itcBest++) {

		for (auto& pp : itcBest->pointsList) {
			MyCoord newPos = MyCoord(pp.x, pp.y);
			itc->pointsList.push_back(newPos);
		}
		itc->clusterID = itcBest->clusterID;

		++itc;
	}

	cout << endl;
	cout << "Opt algo end. Total tests: " << nRounds << ". Total real tests: " << nRealTest << " (" <<
			(((long double) nRealTest)/((long double)nRounds)) * 100.0 << "%)" << endl;
}

void optGoWorst(std::list<MyCoord> &pl, std::vector<CoordCluster> &cv, unsigned int n4cPlus, unsigned int n4c, unsigned int k) {
	std::list<NodeArc> arcs;
	std::list<ExtCoord> ecl;

	for (auto& cc : cv) {
		cc.pointsList.clear();
	}
	for (std::list<MyCoord>::iterator itP = pl.begin(); itP != pl.end(); ++itP) {
		ExtCoord ec;
		ec.x = itP->x;
		ec.y = itP->y;
		ec.clustBelonging = -1;
		ecl.push_back(ec);
	}

	for (std::list<ExtCoord>::iterator itP1 = ecl.begin(); itP1 != ecl.end(); ++itP1) {
		for (std::list<ExtCoord>::iterator itP2 = itP1; itP2 != ecl.end(); ++itP2) {
			if (itP1 != itP2) {
				NodeArc na;
				na.p1 = &(*itP1);
				na.p2 = &(*itP2);
				na.w = weightFunciont(*itP1, *itP2);
				arcs.push_back(na);
			}
		}
	}

	arcs.sort(NodeArc::compare);

	unsigned int cPrint = 0;
	for (std::list<NodeArc>::iterator itNA = arcs.begin(); itNA != arcs.end(); ++itNA) {
		cout << "Arco " << ++cPrint << " " << *(itNA->p1) << " - " <<  *(itNA->p2) << " --> " << itNA->w << " (" << 1.0/itNA->w  << ")m" << endl;
	}

	unsigned int clIdx = 0;
	for (std::list<NodeArc>::iterator itNA = arcs.begin(); itNA != arcs.end(); ++itNA) {

		if ((itNA->p1->clustBelonging < 0) && (itNA->p2->clustBelonging < 0)) {
			MyCoord nc1 = MyCoord(itNA->p1->x, itNA->p1->y);
			itNA->p1->clustBelonging = clIdx;
			cv[clIdx].pointsList.push_back(nc1);
			itNA->p1->toAvoid.push_back(itNA->p2);
			clIdx = (clIdx + 1) % k;

			MyCoord nc2 = MyCoord(itNA->p2->x, itNA->p2->y);
			itNA->p2->clustBelonging = clIdx;
			cv[clIdx].pointsList.push_back(nc2);
			itNA->p2->toAvoid.push_back(itNA->p1);
			clIdx = (clIdx + 1) % k;
		}
		else if ((itNA->p1->clustBelonging >= 0) && (itNA->p2->clustBelonging >= 0)) {
			// try to move someone
			itNA->p1->toAvoid.push_back(itNA->p2);
			itNA->p2->toAvoid.push_back(itNA->p1);

			if (itNA->p1->clustBelonging == itNA->p2->clustBelonging) {
				bool movedLink, movedOthers;

				do {
					movedLink = movedOthers = false;

					int idx1 = (itNA->p1->clustBelonging + 1) % k;
					while (idx1 != itNA->p1->clustBelonging) {
						bool hereOk = true;
						for (auto& pp : ecl) {
							if (pp.clustBelonging == idx1) {
								for (auto& toa : itNA->p1->toAvoid) {
									if ((*toa) == (pp)) {
										hereOk = false;
										break;
									}
								}
							}
							if (!hereOk) break;
						}
						if (hereOk) {
							for (std::vector<MyCoord>::iterator itt1 = cv[itNA->p1->clustBelonging].pointsList.begin(); itt1 != cv[itNA->p1->clustBelonging].pointsList.end(); ++itt1) {
								if ((*itt1) == (*(itNA->p1))) {
									cv[itNA->p1->clustBelonging].pointsList.erase(itt1);
									break;
								}
							}
							MyCoord nc = MyCoord(itNA->p1->x, itNA->p1->y);
							itNA->p1->clustBelonging = idx1;
							cv[idx1].pointsList.push_back(nc);
							movedLink = true;
							//cout << "Breaking someone 1" << endl;
							break;
						}
						idx1 = (idx1 + 1) % k;
					}

					if (idx1 == itNA->p1->clustBelonging) {
						int idx2 = (itNA->p2->clustBelonging + 1) % k;
						while (idx2 != itNA->p2->clustBelonging) {
							bool hereOk = true;
							for (auto& pp : ecl) {
								if (pp.clustBelonging == idx2) {
									for (auto& toa : itNA->p2->toAvoid) {
										if ((*toa) == (pp)) {
											hereOk = false;
											break;
										}
									}
								}
								if (!hereOk) break;
							}
							if (hereOk) {
								for (std::vector<MyCoord>::iterator itt2 = cv[itNA->p2->clustBelonging].pointsList.begin(); itt2 != cv[itNA->p2->clustBelonging].pointsList.end(); ++itt2) {
									if ((*itt2) == (*(itNA->p2))) {
										cv[itNA->p2->clustBelonging].pointsList.erase(itt2);
										break;
									}
								}
								MyCoord nc = MyCoord(itNA->p2->x, itNA->p2->y);
								itNA->p2->clustBelonging = idx2;
								cv[idx2].pointsList.push_back(nc);
								movedLink = true;
								//cout << "Breaking someone 2" << endl;
								break;
							}
							idx2 = (idx2 + 1) % k;
						}
					}

					//if (false) {
					if (!movedLink) {
						// try to move someone other

						int idx1 = (itNA->p1->clustBelonging + 1) % k;
						while (idx1 != itNA->p1->clustBelonging) {
							int p2avoid, p2move;
							p2avoid = p2move = 0;
							for (auto& pp : ecl) {
								if (pp.clustBelonging == idx1) {
									for (auto& toa : itNA->p1->toAvoid) {
										if ((*toa) == (pp)) {
											// move the point 'pp'
											++p2avoid;
											int idxPP = (idx1 + 1) % k;
											while (idxPP != idx1) {
												bool hereOk = true;
												for (auto& ppPP : ecl) {
													if (ppPP.clustBelonging == idxPP) {
														for (auto& toaPP : pp.toAvoid) {
															if ((*toaPP) == (ppPP)) {
																hereOk = false;
																break;
															}
														}
													}
													if (!hereOk) break;
												}
												if (hereOk) {
													//cout << "1 - I can move " << pp << " from " << idx1 << " to " << idxPP << endl;
													++p2move;
													break;
												}
												idxPP = (idxPP + 1) % k;
											}
										}
									}
								}
							}

							if (p2move == p2avoid) {
								//cout << "To move " << p2move << " points 1" << endl;
								for (auto& pp : ecl) {
									if (pp.clustBelonging == idx1) {
										for (auto& toa : itNA->p1->toAvoid) {
											if ((*toa) == (pp)) {
												// move the point 'pp'
												int idxPP = (idx1 + 1) % k;
												while (idxPP != idx1) {
													bool hereOk = true;
													for (auto& ppPP : ecl) {
														if (ppPP.clustBelonging == idxPP) {
															for (auto& toaPP : pp.toAvoid) {
																if ((*toaPP) == (ppPP)) {
																	hereOk = false;
																	break;
																}
															}
														}
														if (!hereOk) break;
													}
													if (hereOk) {
														MyCoord nc = MyCoord(pp.x, pp.y);
														for (std::vector<MyCoord>::iterator itt1 = cv[pp.clustBelonging].pointsList.begin(); itt1 != cv[pp.clustBelonging].pointsList.end(); ++itt1) {
															if ((*itt1) == (pp)) {
																cv[pp.clustBelonging].pointsList.erase(itt1);
																break;
															}
														}
														pp.clustBelonging = idxPP;
														cv[idxPP].pointsList.push_back(nc);

														//cout << "Moving someone 1" << endl;
														break;
													}
													idxPP = (idxPP + 1) % k;
												}
											}
										}
									}
								}
								movedOthers = true;
								break;
							}

							idx1 = (idx1 + 1) % k;
						}

						//cout << "Finish check 1. Moved? " << movedOthers << endl;

						if (!movedOthers) {
							int idx2 = (itNA->p2->clustBelonging + 1) % k;
							while (idx2 != itNA->p2->clustBelonging) {
								int p2avoid, p2move;
								p2avoid = p2move = 0;
								for (auto& pp : ecl) {
									if (pp.clustBelonging == idx2) {
										for (auto& toa : itNA->p2->toAvoid) {
											if ((*toa) == (pp)) {
												// move the point 'pp'
												++p2avoid;
												int idxPP = (idx2 + 1) % k;
												while (idxPP != idx2) {
													bool hereOk = true;
													for (auto& ppPP : ecl) {
														if (ppPP.clustBelonging == idxPP) {
															for (auto& toaPP : pp.toAvoid) {
																if ((*toaPP) == (ppPP)) {
																	hereOk = false;
																	break;
																}
															}
														}
														if (!hereOk) break;
													}
													if (hereOk) {
														//cout << "2 - I can move " << pp << " from " << idx2 << " to " << idxPP << endl;
														++p2move;
														break;
													}
													idxPP = (idxPP + 1) % k;
												}
											}
										}
									}
								}

								if (p2move == p2avoid) {
									//cout << "To move " << p2move << " points 2" << endl;
									for (auto& pp : ecl) {
										if (pp.clustBelonging == idx2) {
											for (auto& toa : itNA->p2->toAvoid) {
												if ((*toa) == (pp)) {
													// move the point 'pp'
													int idxPP = (idx2 + 1) % k;
													while (idxPP != idx2) {
														bool hereOk = true;
														for (auto& ppPP : ecl) {
															if (ppPP.clustBelonging == idxPP) {
																for (auto& toaPP : pp.toAvoid) {
																	if ((*toaPP) == (ppPP)) {
																		hereOk = false;
																		break;
																	}
																}
															}
															if (!hereOk) break;
														}
														if (hereOk) {
															MyCoord nc = MyCoord(pp.x, pp.y);
															for (std::vector<MyCoord>::iterator itt2 = cv[pp.clustBelonging].pointsList.begin(); itt2 != cv[pp.clustBelonging].pointsList.end(); ++itt2) {
																if ((*itt2) == (pp)) {
																	cv[pp.clustBelonging].pointsList.erase(itt2);
																	break;
																}
															}
															pp.clustBelonging = idxPP;
															cv[idxPP].pointsList.push_back(nc);

															//cout << "Moving someone 2" << endl;
															break;
														}
														idxPP = (idxPP + 1) % k;
													}
												}
											}
										}
									}
									movedOthers = true;
									break;
								}

								idx2 = (idx2 + 1) % k;
							}
							//cout << "Finish check 2. Moved? " << movedOthers << endl;
						}
					}
				} while ((!movedLink) && (movedOthers));
			}
		}
		else {
			if (itNA->p1->clustBelonging < 0) {
				MyCoord nc1 = MyCoord(itNA->p1->x, itNA->p1->y);

				if ((int)clIdx == itNA->p2->clustBelonging) {
					clIdx = (clIdx + 1) % k;
				}

				cv[clIdx].pointsList.push_back(nc1);
				itNA->p1->clustBelonging = clIdx;

				itNA->p1->toAvoid.push_back(itNA->p2);
				itNA->p2->toAvoid.push_back(itNA->p1);

				clIdx = (clIdx + 1) % k;
			}
			else {
				MyCoord nc2 = MyCoord(itNA->p2->x, itNA->p2->y);

				if ((int)clIdx == itNA->p1->clustBelonging) {
					clIdx = (clIdx + 1) % k;
				}

				cv[clIdx].pointsList.push_back(nc2);
				itNA->p2->clustBelonging = clIdx;

				itNA->p2->toAvoid.push_back(itNA->p1);
				itNA->p1->toAvoid.push_back(itNA->p2);

				clIdx = (clIdx + 1) % k;
			}
		}


		/*if (itNA->p1->clustBelonging < 0) {
			MyCoord nc = MyCoord(itNA->p1->x, itNA->p1->y);
			itNA->p1->clustBelonging = clIdx;
			cv[clIdx].pointsList.push_back(nc);
			clIdx = (clIdx + 1) % k;
		}

		if (itNA->p2->clustBelonging < 0) {
			MyCoord nc = MyCoord(itNA->p2->x, itNA->p2->y);
			itNA->p2->clustBelonging = clIdx;
			cv[clIdx].pointsList.push_back(nc);
			clIdx = (clIdx + 1) % k;
		}*/
	}

	cout << "AFTER" << endl;
	cPrint = 0;
	for (std::list<NodeArc>::iterator itNA = arcs.begin(); itNA != arcs.end(); ++itNA) {
		++cPrint;
		if (itNA->p1->clustBelonging == itNA->p2->clustBelonging) {
			cout << "Arco " << cPrint << " " << *(itNA->p1) << " - " <<  *(itNA->p2) << " --> " << itNA->w << " (" << 1.0/itNA->w  << ")m" << endl;
		}
	}
}

int main(int argc, char **argv) {
	std::list<MyCoord> pointsList;
	std::vector<CoordCluster> clustersVec;
	int k = 5;
	int scenarioSize = 100;
	int iterationNum = 10;
	int equalize_t = 1;
	int lam = 1;

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
	const std::string &equalizeType = input.getCmdOption("-e");
	const std::string &randomGenerationType = input.getCmdOption("-r");
	const std::string &lambdaInput = input.getCmdOption("-l");

	if (!iterationNumber.empty()) {
		iterationNum = atoi(iterationNumber.c_str());
	}
	if (!scenarioMaxVal.empty()) {
		scenarioSize = atoi(scenarioMaxVal.c_str());
	}
	if (!kValue.empty()) {
		k = atoi(kValue.c_str());
	}

	if (!lambdaInput.empty()) {
		lam = atoi(lambdaInput.c_str());
	}
	else {
		cerr << "Insert the lambda value -l" << endl;
		exit(EXIT_FAILURE);
	}

	if ( (!outputGenerateFileName.empty()) && (!numberOfItemToGenerate.empty()) ) {
		int nItems = atoi(numberOfItemToGenerate.c_str());
		int randomType = 0; //0: random; 1: grid

		if (!randomGenerationType.empty()) {
			int rt = atoi(randomGenerationType.c_str());
			switch (rt) {
			default:
			case 0:
				randomType = 0;
				break;
			case 1:
				randomType = 1;
				break;
			}
		}
		switch (randomType) {
		default:
		case 0:
			generateOutputRandomPoints(outputGenerateFileName, nItems, scenarioSize);
			break;

		case 1:
			generateOutputRandomGridPoints(outputGenerateFileName, nItems, scenarioSize);
			break;
		}
		//generateOutputRandomPoints(outputGenerateFileName, nItems, scenarioSize);
	}

	if (inputFileName.empty()) {
		generateRandomPoints(pointsList, scenarioSize);
	}
	else {
		importPoints(inputFileName, pointsList);
	}

	if (!equalizeType.empty()) {
		equalize_t = atoi(equalizeType.c_str());
		if ((equalize_t < 0) || (equalize_t > 7)) {
			cerr << "Wrong equalyze type [0..7]" << endl;
			exit(EXIT_FAILURE);
		}
	}

	k = pointsList.size() / lam;

	clustersVec.resize(k);

	cout << "Scenario with " << pointsList.size() << " sensors.  Making " << k << " cluster having lambda = " << lam << endl;

	/*for (auto& cv : clustersVec) {
		cout << cv.pointsList.size() << " " << cv.clusterHead
				<< " - MaxCorrelation: " << cv.getMaxCorrelation()
				<< " - AvgCorrelation: " << cv.getAvgCorrelation()
				<< endl;
	}*/

	clustersVec[0] = CoordCluster(pointsList);
	clustersVec[0].chooseFirstClusterHead();

	int idd = 0;
	for (auto& ccc : clustersVec) {
		ccc.clusterID = idd++;
	}

	if (equalize_t == 0) {
		exit(EXIT_SUCCESS);
	}
	else if (equalize_t < 3) {
		//cout << "Start k_2MMalgo" << endl; fflush(stdout);
		k_2MMalgo(clustersVec, k, iterationNum);
		//cout << "End k_2MMalgo" << endl; fflush(stdout);

		//cout << "Start equalizer" << endl; fflush(stdout);
		if (equalize_t == 1) {
			equalize1(clustersVec, pointsList.size(), k);
		}
		else {
			equalize2(clustersVec, pointsList.size(), k);
		}
		//cout << "End equalizer" << endl; fflush(stdout);
	}
	else if (equalize_t == 3) {
		//while (!allEqualized3(clustersVec)) {
		unsigned int n4cluster = ((unsigned int) pointsList.size()) / k;
		int remainingP = ((int) pointsList.size()) % k;

		do {
			//cout << "Start k_2MMalgo" << endl; fflush(stdout);
			k_2MMalgo3(clustersVec, k, iterationNum);
			//cout << "End k_2MMalgo" << endl; fflush(stdout);

			//printAllCluters(clustersVec, "BEFORE EQUALIZER", false);

			//cout << "Start equalizer" << endl; fflush(stdout);
			equalize3(clustersVec, pointsList.size(), k, n4cluster, remainingP);
			if (remainingP > 0) {
				--remainingP;
			}
			//for (auto& c : clustersVec) { c.algo3freeze = true; }
			//cout << "End equalizer" << endl; fflush(stdout);

			//printAllCluters(clustersVec, "AFTER EQUALIZER", false);

			if (numNotEqualized3(clustersVec) <= 1) {
				break;
			}
		} while (true);
	}
	else if (equalize_t == 4) {
		// new algorithm
		unsigned int n4cluster = lam;
		int remainingP = ((int) pointsList.size()) % lam;

		// randomize the input in the clusters
		randomizeClusters(clustersVec, k, n4cluster, remainingP);

		makeSwaps(clustersVec, k);
	}
	else if (equalize_t == 5) {
		// opt algorithm
		unsigned int n4cluster = lam;
		unsigned int remainingP = ((int) pointsList.size()) % lam;
		unsigned int n4clusterPlus = lam + remainingP;

		optAlgo(pointsList, clustersVec, n4clusterPlus, n4cluster, k);

	}
	else if (equalize_t == 6) {
		// iterate over worst link
		unsigned int n4cluster = lam;
		unsigned int remainingP = ((int) pointsList.size()) % lam;
		unsigned int n4clusterPlus = lam + remainingP;

		optGoWorst(pointsList, clustersVec, n4clusterPlus, n4cluster, k);

	}
	else if (equalize_t == 7) {
		// just random
		unsigned int n4cluster = lam;
		int remainingP = ((int) pointsList.size()) % lam;

		// randomize the input in the clusters
		randomizeClusters(clustersVec, k, n4cluster, remainingP);
	}

	int sumElements = 0;
	//double maxCorrelation = 0;
	//double maxAvgCorrelation = 0;
	for (auto& cv : clustersVec) {
		//double actMaxCorrelation = cv.getMaxCorrelation();
		//double actAvgCorrelation = cv.getAvgCorrelation();
		//cout << cv.pointsList.size() << " " << cv.clusterHead << " - MaxCorrelation: " << actMaxCorrelation << " - AvgCorrelation: " << actAvgCorrelation << " - Moran's I: " << cv.getMoransI() << endl;
		//cout << "[" << cv.clusterID << "] "<< cv.pointsList.size() << " " << cv.clusterHead << " - Moran's I: " << cv.getMoransI() << " - MAX corr: " << cv.getMaxCorrelation() << endl;

		//if (actMaxCorrelation > maxCorrelation) maxCorrelation = actMaxCorrelation;
		//if (actAvgCorrelation > maxAvgCorrelation) maxAvgCorrelation = actAvgCorrelation;

		sumElements += cv.pointsList.size();
	}
	double maxTotalCorr = getSystemMaxCorrelation(clustersVec);
	cout << "Total elements: " << sumElements << endl;
	cout << "Maximum MAX correlation in " << inputFileName << " is: " << maxTotalCorr << " cioè " << 1.0 / maxTotalCorr << " metri" << endl;
	cout << "StatMaxCorr " << maxTotalCorr << " " << 1.0 / maxTotalCorr << endl;
	cout << "Maximum AVG correlation: " << getSystemAvgCorrelation(clustersVec) << endl;
	cout << "END Moran Full: " << calculateFullMoranIndex(clustersVec) << endl;
	cout << "END Moran Full Const: " << calculateFullMoranIndex_const(clustersVec) << endl;

	if (!dotFileOutput.empty()) {
		generateDOTfile(dotFileOutput, clustersVec, ((double) scenarioSize)/50.0);
	}

	cout << "End!!!" << endl;
	return EXIT_SUCCESS;
}

