/**
 * \file   point.hpp
 * \class  Point
 * 
 * \author parasara
 * \author Lucas Brown
 * \date   Feb 12, 2012
 * \date   April 1, 2019
 * 
 * \brief   LMBTODO
 */

#ifndef POINT_H_
#define POINT_H_

#include <vector>

class Point 
{
  public:
    Point();
    Point(int dimensions);
    Point(int dimensions, std::vector<double> coordinates);
    ~Point();

    void print();

    int getDimensions();
    void setDimensions(int val);
    std::vector<double> getCoordinates();
    void setCoordinates(std::vector<double> vec);
    void setCoordinate(int index, double val);
    double getCoordinate(int index);    

  private:
  
    int dimensions;
    std::vector<double> coordinates;
};

#endif /* POINT_H_ */