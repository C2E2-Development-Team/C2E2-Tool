/**
 * \file   simulator.hpp
 * \class  Simulator
 *
 * \author parasara
 * \author Lucas Brown
 * \date   Feb 1, 2013
 * \date   April 1, 2019
 * 
 * \brief  LMBTODO
 */


#ifndef SIMULATOR_H_
#define SIMULATOR_H_

#include <ppl.hh>
#include <stack>
#include <string>
#include <vector>

#include "annotation.hpp"
#include "checker.hpp"
#include "linear-set.hpp"
#include "point.hpp"
#include "rep-point.hpp"

class Simulator
{
  public:
    Simulator();
    ~Simulator();

    int simulate();
    int hybridSimulation();
    int hybridSimulation(Point origin);
    void simulatePoint(Point point, int mode);
    void simulatePoint(Point point, int mode, int idx);

    std::vector<Point> getRepresentativeCover();  // Unsure / Simulator
    Point getPointFromPoly(std::NNC_Polyhedron Poly);  // Unsure
    void print();

    int getDimensions();
    void setDimensions(int val);
    int getInitialModeIdx();
    void setInitialModeIdx(int val);

    double getAbsoluteError();
    void setAbsoluteError(double val);
    double getRelativeError();
    void setRelativeError(double val);
    double getTimeStep();
    void setTimeStep(double val);
    double getTimeHorizon();
    void setTimeHorizon(double val);

    std::string getVisualizeFilename();
    void setVisualizeFilename(std::string str);
    std::string getExecutable();
    void setExecutable(std::string str);

    LinearSet getInitialSet();
    void setInitialSet(LinearSet obj);
    LinearSet getUnsafeSet();
    void setUnsafeSet(LinearSet obj);

    std::vector<int> getModeLinear();  // LMBTODO: Band-aid
    void setModeLinear(std::vector<int> vec);  // LMBTODO: Band-aid
    std::vector<double> getKConsts();  // LMBTODO: Band-aid
    void setKConsts(std::vector<double> vec);  // LMBTODO: Band-aid

    Checker getChecker();

  protected:

    int dimensions;
    int initial_mode_idx;

    double absolute_error;
    double relative_error;
    double time_step;
    double time_horizon;

    std::string visualize_filename;    
    std::string executable;

    LinearSet initial_set;
    LinearSet unsafe_set;

    Checker checker;

    std::vector<int> mode_linear;
    std::vector<double> k_consts;
};

#endif /* SIMULATOR_H_ */