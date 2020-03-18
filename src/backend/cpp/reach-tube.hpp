/**
 * \file   reach-tube.hpp
 * \class  ReachTube
 * 
 * \author parasara
 * \author Lucas Brown
 * \date   July, 2014
 * \date   April 2, 2019
 * 
 * \brief LMBTODO
 */

#ifndef REACHTUBE_H_
#define REACHTUBE_H_

#include <ppl.hh>
#include <stack>
#include <string>
#include <vector>

#include "annotation.hpp"
#include "point.hpp"
#include "rep-point.hpp"

class ReachTube
{
  public:
	ReachTube();
	~ReachTube();

    int getSize();
    Point getUpperBound(int index);
	Point getLowerBound(int index);
    void parseInvariantTube(char const* filename, int hasMode);
    void printReachTube(const std::string, int flag);
    void clear(int from);
    ReachTube bloatReachTube(std::vector<double> delta_array, Annotation annotation);

    int getNextSetStack(std::stack<RepPoint>& itr_stack, RepPoint parent_rep_point);
    void addGuards(std::vector<std::pair<std::NNC_Polyhedron, int> > guards);
    double getMinCoordinate(int dim, int cur_mode);
    double getMaxCoordinate(int dim, int cur_mode);
    int checkIntersection(int cur_mode, Point cur_point, 
        std::vector<double> delta_array);
    double getMinTime(int cur_mode, Point cur_point, 
        std::vector<double> delta_array);

    int getDimensions();
    void setDimensions(int val);
	int getMode();
	void setMode(int val);  // Setting mode also sets the isReachTube bool
    std::vector<int> getModeVec();
    void setModeVec(std::vector<int> vec);
    void addLowerBoundState(Point obj);
    void addUpperBoundState(Point obj);
	
  private:
  
	int dimensions;
	int isReachTube;
	int reachTubeMode;
	std::vector<int> color;
	std::vector<int> mode;
	std::vector<Point> upper_bound;
	std::vector<Point> lower_bound;
  
};

#endif /* REACHTUBE_H_ */
