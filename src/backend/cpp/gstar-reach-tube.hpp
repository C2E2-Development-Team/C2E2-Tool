/*
 * generalized-reach-tube.hpp
 *
 *  Created on: October 21, 2018
 *      Author: 
 */

#ifndef GSTARREACHTUBE_H_
#define GSTARREACHTUBE_H_

#include <iostream>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <cstdio>
#include <time.h>
#include <stdlib.h>

#include "generalized-star.hpp"
#include "reach-tube.hpp"
#include "simulator.hpp"
#include "Eigen/Dense"

using namespace std;

class gReachTube
{
private:
  // the dimension of the reach tube
  int dimension;

  // the mode corresponding to the reach tube
  int mode;

  // the tube object it self
  // stored in a vector of generalized star
  vector<gstar> tube;

  // a pointer to the simulator of the model
  // the simulator is inherited from other part of the program
  Simulator modelSim;

  // the initial set of the reachtube
  // represented by a generalized star
  gstar initialSet;

  // whether the reachtube is already bloated
  int bloated;

  // the epsilon value used to bloat the reachtube
  // mentioned in theorem 4 of the paper
  // the simulation trace between time steps is bounded by this value
  double epsilon;

  vector<double> computeBasis(vector<double> ptn, vector<double> center);
  double computeEps(vector<double> ptn1, vector<double> ptn2);
  void bloatBasis(vector<double> &basis);

public:
  // Constructor
  gReachTube();
  gReachTube(gstar star);
  // gReachTube(ReachTube rtube);

  // Data access functions
  int getDimension();
  int getMode();
  vector<gstar> getTube();
  gstar getStar(int index);
  double getEpsilon();
  int isBloated();

  // Data modification functions
  void setDimension(int newDimension);
  void setMode(int newMode);
  void addToTube(gstar star);
  void setSimulator(Simulator simu);

  // other member functions
  void computeSimulationTube();
  void computeSimulationTube_noSim();
  void doRandomSimulation();
  vector<Point> parseSimRes(const char *filename, int hasMode);
  void printTubeFile(string fn, int flag);
  void printBloatTubeFile(string fn, int flag);
  void tubeToReachTube(ReachTube &rt);
  void bloatTube();
};

#endif /* INITIALSET_H_ */
