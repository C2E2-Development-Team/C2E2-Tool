/*
 * generalized-star.cpp
 *
 *  Created on: October 19, 2018
 *      Author: 
 */

#include "generalized-star.hpp"

using namespace std;

/*  
gstar()
    Input:  
        None
    Return value:
        None
    Description:
        Create an enpty generalized star
*/
gstar::gstar()
{
    basis = vector<vector<double>>();
    center = vector<double>();
    alpha = vector<double>();
    mode = -1;
    dimension = 0;
    time_loc = -1;
}

/*
gstar(RepPoint ptn)
    Input:
        RepPoint ptn: a reppoint
    Return value:
        None
    Description
        Create a generalized star from a reppoint
*/
gstar::gstar(RepPoint ptn)
{
    dimension = ptn.getDimensions();
    mode = ptn.getMode();
    vector<double> temp = ptn.getState().getCoordinates();
    time_loc = temp[0];
    for (int i = 0; i < dimension; i++)
    {
        center.push_back(temp[i + 1]);
    }
    temp = ptn.getDeltaArray();
    for (int i = 0; i < dimension; i++)
    {
        vector<double> temp_vec = vector<double>();
        for (int j = 0; j < dimension; j++)
        {
            temp_vec.push_back(0);
        }
        temp_vec[i] = temp[i];
        basis.push_back(temp_vec);
    }
    for (int i = 0; i < dimension; i++)
    {
        alpha.push_back(1.0);
    }
}

/*
gstar(RepPoint *ptn)
    Input:
        RepPoint *ptn: A pointer to reppoint
    Return value:
        None
    Description
        Create a generalized star from a reppoint pointer
*/
gstar::gstar(RepPoint *ptn)
{
    dimension = ptn->getDimensions();
    mode = ptn->getMode();
    vector<double> temp = ptn->getState().getCoordinates();
    time_loc = temp[0];
    for (int i = 0; i < dimension; i++)
    {
        center.push_back(temp[i + 1]);
    }
    temp = ptn->getDeltaArray();
    for (int i = 0; i < dimension; i++)
    {
        vector<double> temp_vec = vector<double>();
        for (int j = 0; j < dimension; j++)
        {
            temp_vec.push_back(0);
        }
        temp_vec[i] = temp[i];
        basis.push_back(temp_vec);
    }
    for (int i = 0; i < dimension; i++)
    {
        alpha.push_back(1.0);
    }
}

/*
vector<vector<double>> getBasis()
    Input:
        None
    Return value:
        vector<vector<double>>: A 2-d vector contain all basis vector for the star
    Description:
        return a 2-d vector contain all basis vector of the star
*/
vector<vector<double>> gstar::getBasis()
{
    return basis;
}

vector<double> gstar::getBasis(int i)
{
    if (i >= basis.size() || i < 0)
    {
        throw std::invalid_argument("index out of range");
        return vector<double>();
    }
    return basis[i];
}

/*
double getBasis(int i, int j)
    Input:
        int i: the index for the basis vector
        int j: the index for the specific value in a basis vector
    Return value:
        double: the specific value for basis vector i and variable j
    Description:
        return a basis value for specific basis vector i for variable j
*/
double gstar::getBasis(int i, int j)
{
    if (i >= basis.size() || j >= basis[0].size() || i < 0 || j < 0)
    {
        throw std::invalid_argument("index out of range");
        return 0;
    }
    return basis[i][j];
}

/*
vector<double> getCenter()
    Input:
        None
    Return value:
        vector<double>: the center of star
    Description:
        return the center of the star
*/
vector<double> gstar::getCenter()
{
    return center;
}

/*
double getCenter(int idx)
    Input:
        int idx: a variable index
    Return value:
        double: the center value for variable idx
    Description:
        return the center value for variable idx
*/
double gstar::getCenter(int idx)
{
    if (idx >= center.size() || idx < 0)
    {
        throw std::invalid_argument("index out of range");
        return 0;
    }
    return center[idx];
}

/*
int getDimension()
    Input:
        None
    Return value:
        int: dimension of the star
    Description:
        return the dimension of the star
*/
int gstar::getDimension()
{
    return dimension;
}

/*
int getMode()
    Input:
        None
    Return value:
        int: mode of the star
    Description:
        return the mode of the star
*/
int gstar::getMode()
{
    return mode;
}

/*
int getTime()
    Input:
        None
    Return value:
        int: time of the star
    Description:
        return the time of the star
*/
double gstar::getTime()
{
    return time_loc;
}

/*
vector<double> getAlpha()
    Input:
        None
    Return value:
        vector<double>: a vector of alpha value 
    Description:
        return all alpha value (predicate) of the star
*/
vector<double> gstar::getAlpha()
{
    return alpha;
}

/*
double getAlpha(int idx)
    Input: 
        int idx: variable idx
    Return value:
        double: alpha value for variable idx
    Description:
        return the alpha value value for variable idx
*/
double gstar::getAlpha(int idx)
{
    if (idx >= alpha.size())
    {
        throw std::invalid_argument("index out of range");
        return 0;
    }
    return alpha[idx];
}

/*
void setBasis(vector<vector<double>> newBasis)
    Input:
        vector<vector<double>> newBasis: new basis vectors
    Return value:
        None
    Description:
        set the basis value of the star
*/
void gstar::setBasis(vector<vector<double>> newBasis)
{
    basis = newBasis;
}

/*
void setBasis(vector<double> baseVec, int varIdx)
    Input:
        vector<double> baseVec: basis vector for specific dimension
        int varIdx: index to specify the dimension
    Return value:
        None
    Description:
        set the basis vector for dimension specified by varIdx
*/
void gstar::setBasis(vector<double> baseVec, int varIdx)
{
    basis[varIdx] = baseVec;
}

/*
void setBasis(double baseVal, int varIdx, int vecIdx)
    Input:
        vector<double> baseVec: basis vector for specific dimension
        int varIdx: index to specify the dimension
    Return value:
        None
    Description:
        set the basis vector for dimension specified by varIdx
*/
void gstar::setBasis(double baseVal, int varIdx, int vecIdx)
{
    basis[varIdx][vecIdx] = baseVal;
}

/* 
void setCenter(vector<double> newCenter)
    Input: 
        vector<double> newCenter: vector contain new center
    Return value:
        None
    Description:
        set new center value
*/
void gstar::setCenter(vector<double> newCenter)
{
    center = newCenter;
}

/* 
void setCenter(double centerVal, int varIdx)
    Input: 
        double centerVal: value of the center
        int varIdx: index of variable
    Return value:
        None
    Description:
        set the center value for variable given by varIdx
*/
void gstar::setCenter(double centerVal, int varIdx)
{
    center[varIdx] = centerVal;
}

/*
void setDimension(int newDimension)
    Description:
        set dimension value of the star
*/
void gstar::setDimension(int newDimension)
{
    dimension = newDimension;
}

/*
void setMode(int newMode)
    Description:
        Set mode value information of the star
*/
void gstar::setMode(int newMode)
{
    mode = newMode;
}

/*
void setTime(double newTime)
    Description:
        Set time information of the star
*/
void gstar::setTime(double newTime)
{
    time_loc = newTime;
}

/*
void setAlpha(vector<double> newAlpha)
    Description:
        Set all alpha value (predicate) of the model
*/
void gstar::setAlpha(vector<double> newAlpha)
{
    alpha = newAlpha;
}

/*
void setAlpha(double newAlpha, int idx)
    Input:
        double newAlpha: the new alpha value
        int idx: index of variable
    Return value:
        None
    Description:
        Set the alpha value for variable idx to value newAlpha
*/
void gstar::setAlpha(double newAlpha, int idx)
{
    alpha[idx] = newAlpha;
}

/*
void normalize()
    Input:
        None
    Output:
        None
    Description:
        In part 4 of the paper, a method to bloat the reachable set is mentioned
        The requirement is max(||alpha||)<=1/n for initial set
        This is the purpose of this function
        The function will increase the length of basis vector while decrease the alpha value 
        by the same factor to achieve the requirement above
        The function only need to be called for the initial set of the model
*/
void gstar::normalize()
{
    double max_alpha = *(max_element(alpha.begin(), alpha.end()));
    double factor = max_alpha * ((double)dimension);
    for (int i = 0; i < dimension; i++)
    {
        alpha[i] = alpha[i] / factor;
        for (int j = 0; j < dimension; j++)
        {
            basis[i][j] = basis[i][j] * factor;
        }
    }
}

/*
void print()
    Input:
        None
    Output:
        None
    Description:
        Print the star to termainl. Mainly for debugging purpose
*/
void gstar::print()
{
    cout << "*************" << endl;
    cout << "Time: " << time_loc << endl;
    cout << "Center: ( ";
    for (int i = 0; i < dimension; i++)
    {
        cout << center[i] << " ";
    }
    cout << ")" << endl;
    cout << "Basis: " << endl;
    for (int i = 0; i < dimension; i++)
    {
        cout << "( ";
        for (int j = 0; j < dimension; j++)
        {
            cout << basis[i][j] << " ";
        }
        cout << ")" << endl;
    }
    cout << "Alpha: ( ";
    for (int i = 0; i < dimension; i++)
    {
        cout << alpha[i] << " ";
    }
    cout << ")" << endl;
    cout << "*************" << endl;
}

/*
vector<Point> getSimPoint()
    Input:
        None
    Output:
        vector<Point> getSimPoint(): a vector of point, will do simulation starting from these points
    Description:
        This function will return a vector of dimension+1 points. 
        The returned points will be used as starting points of simulation.
        The first point in the vector is the center of the star.
        The rest are given by center+alpha_i*basis_i where center and basis_i are vector and alpha_i scalar, is integer and 0<=i<dimension
*/
vector<Point> gstar::getSimPoint()
{
    vector<Point> simPoint = vector<Point>();

    // push the center of the star into the list of simulation point
    vector<double> cenArray;
    cenArray.push_back(time_loc);
    for (int i = 0; i < dimension; i++)
    {
        cenArray.push_back(center[i]);
    }
    simPoint.push_back(Point(dimension + 1, cenArray));

    // other simulatioin points are given by x_i=x_0+v_i
    for (int i = 0; i < dimension; i++)
    {
        vector<double> coArray;
        coArray.push_back(time_loc);
        for (int j = 1; j < dimension + 1; j++)
        {
            coArray.push_back(center[j - 1] + basis[i][j - 1]);
        }
        Point tempPoint = Point(dimension + 1, coArray);
        simPoint.push_back(tempPoint);
    }
    return simPoint;
}