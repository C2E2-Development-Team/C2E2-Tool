/**
 * \file   rep-point.cpp
 * \class  RepPoint
 *
 * \author Bolun
 * \author Lucas Brown
 * \date   Feb 15, 2016
 * \date   April 1, 2019
 * 
 * \brief  LMBTODO
 */

#include <iostream>

#include "rep-point.hpp"

RepPoint::RepPoint()
{ 
    dimensions = -1;
    mode = -1;
    refine_time = 0;  // This must be intitialized for our verification loop
    has_parent_state = false;
}

RepPoint::~RepPoint()
{ }

bool RepPoint::hasParentState()
{
    return has_parent_state;
}
void RepPoint::print()
{
    std::cout << "********** Rep Point **********" << std::endl;
    std::cout << "    dimensions: " << dimensions << std::endl;
    std::cout << "    refine_time: " << refine_time << std::endl;
    std::cout << "    mode: " << mode << std::endl;
    std::cout << "    delta_arraay: [";
    for (auto v: delta_array)
    {
        std::cout << v << " ";
    }
    std::cout << "]" << std::endl;
    std::cout << "----- State -----" << std::endl;
    state.print();
    std::cout << "    has_parent_state: " << has_parent_state << std::endl;
    if (has_parent_state)
    {
        std::cout << "    parent_mode: " << parent_mode << std::endl;
        std::cout << "    parent_delta_arraay: [";
        for (auto v: parent_delta_array)
        {
            std::cout << v << " ";
        }
        std::cout << "]" << std::endl;
        std::cout << "----- Parent State -----" << std::endl;
        parent_state.print();
    }    
    std::cout << "********** End Rep Point **********" << std::endl;
}

int RepPoint::getDimensions()
{
    return dimensions;
}
void RepPoint::setDimensions(int val)
{
    dimensions = val;
}
int RepPoint::getRefineTime()
{
    return refine_time;
}
void RepPoint::setRefineTime(int val)
{
    refine_time = val;
}
int RepPoint::getMode()
{
    return mode;
}
void RepPoint::setMode(int val)
{
    mode = val;
}
int RepPoint::getParentMode()
{
    return parent_mode;
}
void RepPoint::setParentMode(int val)
{
    parent_mode = val;
}
std::vector<double> RepPoint::getDeltaArray()
{
    return delta_array;
}
void RepPoint::setDeltaArray(std::vector<double> vec)
{
    delta_array = vec;
}
std::vector<double> RepPoint::getParentDeltaArray()
{
    return parent_delta_array;
}
void RepPoint::setParentDeltaArray(std::vector<double> vec)
{
    parent_delta_array = vec;
}
Point RepPoint::getState()
{
    return state;
}
void RepPoint::setState(Point obj)
{
    state = obj;
}
Point RepPoint::getParentState()
{
    return parent_state;
}
void RepPoint::setParentState(Point obj)
{
    has_parent_state = true;
    parent_state = obj;
}