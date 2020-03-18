/*
 * \file   linear-set.hpp
 * \class  LinearSet
 * 
 * \author parasara
 * \author Lucas Brown
 * \date   July 28, 2014
 * \date   April 1, 2019
 * 
 * \brief LMBTODO 
 */

#ifndef LINEARSET_H_
#define LINEARSET_H_

#include <vector>

#include "point.hpp"

class LinearSet 
{
  public:
	LinearSet();
	~LinearSet();
	
    int isInternal(Point point);
    double getMax(int dimension);
	double getMin(int dimension);
    int hasIntersection(Point lower_point, Point upper_point);
    void print();
    
    int getDimensions();
    void setDimensions(int val);
    int getNumEqns();
    void setNumEqns(int val);

    std::vector<double> getMatrix();
    void setMatrix(std::vector<double> mat);
    double getMatrixElement(int row, int col);
    void setMatrixElement(int row, int col, double value);
	
    std::vector<double> getB();
    void setB(std::vector<double> vec);
	double getBElement(int index);
	void setBElement(int index, double value);

  private:

	int dimensions;
	int number_of_equations;
    std::vector<double> matrix;
    std::vector<double> b;

    double findMaxMin(int dimension_ID, bool max);
};

#endif /* LINEARSET_H_ */
