/**
 * \file   point.cpp
 * \class  Point
 *
 * \author parasara
 * \author Lucas Brown
 * \date   Feb 12, 2012
 * \date   April 1, 2019
 * 
 * \brief LMBTODO
 */

#include <iostream>
#include <vector>

#include "point.hpp"

Point::Point()
{ }

Point::Point(int dimensions)
{ 
    this->dimensions = dimensions;
    this->coordinates = std::vector<double>(dimensions, 0);
}

Point::Point(int dimensions, std::vector<double> coordinates)
{ 
    this->dimensions = dimensions;
    this->coordinates = coordinates;
}

Point::~Point()
{ }

void Point::print()
{ 
    std::cout << "********** Point **********" << std::endl;
    std::cout << "    Dimensions: " << dimensions << std::endl;
    std::cout << "    Coordinates: " << " { ";
    for (auto p: coordinates)
    {
        std::cout << p << " ";
    }
    std::cout << "}" << std::endl;
    std::cout << "********** End Point **********" << std::endl;
}

int Point::getDimensions()
{
    return dimensions;
}
void Point::setDimensions(int val)
{
    dimensions = val;
}
std::vector<double> Point::getCoordinates()
{
    return coordinates;
}
void Point::setCoordinates(std::vector<double> vec)
{
    coordinates = vec;
}

double Point::getCoordinate(int index)
{
    // LMBTODO: Add error checking for index
    return coordinates[index];
}
void Point::setCoordinate(int index, double val)
{
    // LMBTODO: Addd error checking for index
    coordinates[index] = val;
}