/**
 * \file   rep-point.hpp
 * \class  RepPoint
 *
 * \author Bolun
 * \author Lucas Brown
 * \date   Feb 15, 2016
 * \date   April 1, 2019
 * 
 * \brief  LMBTODO
 */

#ifndef REPPOINT_H_
#define REPPOINT_H_

#include <vector>

#include "point.hpp"

class RepPoint
{
  public:
    RepPoint();
    ~RepPoint();

    bool hasParentState();
    void print();

    int getDimensions();
    void setDimensions(int val);
    int getRefineTime();
    void setRefineTime(int val);
    int getMode();
    void setMode(int val);
    int getParentMode();
    void setParentMode(int val);

    std::vector<double> getDeltaArray();
    void setDeltaArray(std::vector<double> vec);
    std::vector<double> getParentDeltaArray();
    void setParentDeltaArray(std::vector<double> vec);
    
    Point getState();
    void setState(Point obj);
    Point getParentState();
    void setParentState(Point obj);

  private:
    
    bool has_parent_state;
    int dimensions;
    int refine_time;
    int mode;
    int parent_mode;
    
    std::vector<double> delta_array;
    std::vector<double> parent_delta_array;

    Point state;
    Point parent_state;
};

#endif /* REPPOINT_H_ */