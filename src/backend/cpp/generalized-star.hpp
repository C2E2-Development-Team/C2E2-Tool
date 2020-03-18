/*
 * generalized-star.hpp
 *
 *  Created on: October 19, 2018
 *      Author: 
 */

#ifndef GENERALIZEDSTAR_H_
#define GENERALIZEDSTAR_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>

#include "rep-point.hpp"
#include "point.hpp"

using namespace std;

class gstar
{
  private:
    // the basis vector for the generalized star (time is not included)
    // 2-d vector with size dimension*dimension
    // first dimension is the index for each basis vector (b1-bn)
    // second dimension is the index for each entry for each basis vector (b1[1]-b1[n])
    vector<vector<double>> basis;

    // the center for the generalized star (time is not included)
    // the length of the vector is equal to dimension
    vector<double> center;

    // the dimension of the generalized star. Time is not counted as one of the dimension
    int dimension;

    // the mode the generalized star belong to
    int mode;

    // the time for the generalilzed star
    double time_loc;

    // the predicate of the generalized star
    // In the paper, the generalized star is defined in page 4
    // the alpha in the definition is given by predicate
    // P(\alpha)=C\alpha<=d
    // Since the predicate is constant all the time
    // we alpha is the same all the time.
    // Since all the initial set for C2E2 is linear (or even rectangular)
    // it is possible to easily compute alpha for once and use it later
    vector<double> alpha;

  public:
    // Constructor
    gstar();
    gstar(int dimension);
    gstar(RepPoint ptn);
    gstar(RepPoint *ptn);

    // Data access function
    vector<vector<double>> getBasis();
    vector<double> getBasis(int i);
    double getBasis(int i, int j);
    vector<double> getCenter();
    double getCenter(int idx);
    int getDimension();
    int getMode();
    double getTime();
    vector<double> getAlpha();
    double getAlpha(int idx);

    // Data modification function
    void setBasis(vector<vector<double>> newBasis);
    void setBasis(vector<double> baseVec, int varIdx);
    void setBasis(double baseVal, int varIdx, int vecIdx);
    void setCenter(vector<double> newCenter);
    void setCenter(double centerVal, int varIdx);
    void setDimension(int newDimension);
    void setMode(int newMode);
    void setTime(double newTime);
    void setAlpha(vector<double> newAlpha);
    void setAlpha(double newAlpha, int idx);

    // Other member functions
    void normalize();
    void print();
    vector<Point> getSimPoint();
};

class bloated_gstar
{
  public:
    double time_step;
    vector<double> upper_bound;
    vector<double> lower_bound;
};

#endif /* INITIALSET_H_ */
