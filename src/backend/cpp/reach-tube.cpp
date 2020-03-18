/**
 * \file   reach-tube.cpp
 * \class  ReachTube
 * 
 * \author parasara
 * \author Lucas Brown
 * \date   July, 2014
 * \date   April 2, 2019
 * 
 * \brief LMBTODO
 */


#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "annotation.hpp" 
#include "linear-set.hpp"
#include "reach-tube.hpp"
#include "simulator.hpp"
#include "verifier.hpp"

ReachTube::ReachTube()
{ }

ReachTube::~ReachTube()
{ }

int ReachTube::getSize()
{
	if (upper_bound.size() == lower_bound.size())
    {
		return upper_bound.size();
    }
	else
    {
		return 0;
    }
}

Point ReachTube::getUpperBound(int index)
{
	return upper_bound.at(index);
}

Point ReachTube::getLowerBound(int index)
{
	return lower_bound.at(index);
}

void ReachTube::parseInvariantTube(char const* filename, int has_mode)
{
    FILE* tRFile;  // LMBTODO Rename
	tRFile = fopen(filename, "r");

	if (tRFile == NULL)
    {
		std::cout << "ERROR: Invariant tube file is NULL" << std::endl;
	}
	
	if (has_mode == 1)
    {
		int mode;
		fscanf(tRFile, "%d", &mode);
	}

    int breakloop = 0;
    int number = 0;
    double buffer_reader;
    int flag = 0;
	while (fscanf(tRFile, "%lf", &buffer_reader) != EOF && breakloop == 0)
    {
		number++;
		Point trace_point(dimensions+1);
		trace_point.setCoordinate(0, buffer_reader);
		for (int j=0; j < dimensions; j++)
        {
			if (fscanf(tRFile, "%lf", &buffer_reader) != EOF)
            {
				trace_point.setCoordinate(j+1, buffer_reader);
			}
			else
            {
				breakloop = 1;
			}
		}

		if (breakloop == 0 && flag == 0)
        {
			lower_bound.push_back(trace_point);
		}
		if (breakloop == 0 && flag == 1)
        {
			upper_bound.push_back(trace_point);
		}
		flag = 1-flag;
	}
	fclose(tRFile);
}

void ReachTube::printReachTube(const std::string file_name, int flag)
{
	std::ofstream invFile;
	if (flag == 0)
    {
		invFile.open(file_name);
		invFile << " " << reachTubeMode << std::endl;
	}
	if (flag == 1)
    {
		if (lower_bound.size() >= 1)
        {
			invFile.open(file_name, std::ios::app);
			invFile << " % mode " << reachTubeMode << std::endl;
		}
	}
	if (flag == 2)
    {
		if (lower_bound.size() >= 1)
        {
			invFile.open(file_name, std::ios::app);
			invFile << " % mode " << reachTubeMode << " unsafe"<< std::endl;
		}
	}

    for (int i=0; i < lower_bound.size(); i++)
    {
        Point temp_point = lower_bound.at(i);
        for (int j=0; j < temp_point.getDimensions(); j++)
        {
            invFile << " " << temp_point.getCoordinate(j); 
        }
        invFile << std::endl;

        temp_point = upper_bound.at(i);
        for (int j=0; j < temp_point.getDimensions(); j++)
        {
            invFile << " " << temp_point.getCoordinate(j);
        }
        invFile << std::endl;
    }

    invFile.close();
}

void ReachTube::clear(int from)
{
	upper_bound.erase(upper_bound.begin()+from, upper_bound.end());
	lower_bound.erase(lower_bound.begin()+from, lower_bound.end());
}

ReachTube ReachTube::bloatReachTube(std::vector<double> delta_array, 
    Annotation annotation)
{
    ReachTube bloated_tube;
    bloated_tube.setDimensions(dimensions);
    bloated_tube.setMode(reachTubeMode);

    for (int i=0; i < lower_bound.size(); i++)
    {
        Point lower_point_i = lower_bound.at(i);
        Point upper_point_i = upper_bound.at(i);

        Point lower_point(lower_point_i.getDimensions());
        Point upper_point(upper_point_i.getDimensions());
        lower_point.setCoordinate(0, lower_point_i.getCoordinate(0));
        upper_point.setCoordinate(0, upper_point_i.getCoordinate(0));

        for (int j=1; j < lower_point.getDimensions(); j++)
        {
            double delta = delta_array[j-1];
            double epsilon = annotation.getED(
                delta, lower_point_i.getCoordinate(0),
                upper_point_i.getCoordinate(0), reachTubeMode);
            
            double a = lower_point_i.getCoordinate(j);
            double b = upper_point_i.getCoordinate(j);
            double temp;
            if (a > b)
            {
                temp = b;
                b = a;
                a = temp;
            }
            
            lower_point.setCoordinate(j, a - epsilon);
            upper_point.setCoordinate(j, b + epsilon);
        }
        bloated_tube.addLowerBoundState(lower_point);
        bloated_tube.addUpperBoundState(upper_point);
    }
    return bloated_tube;
}

int ReachTube::getNextSetStack(std::stack<RepPoint>& itr_stack, 
    RepPoint parent_rep_point)
{
    std::vector<int> modes_in_set;
	Point parent_point;
	std::vector<double> parent_delta_array;
	std::stack<RepPoint> temp_stack;
	RepPoint temp_rep_point;
	Point temp_point;
	int refinetime = parent_rep_point.getRefineTime();
	int init_mode;

	if (parent_rep_point.hasParentState())
    {
		parent_point = parent_rep_point.getParentState();
		parent_delta_array = parent_rep_point.getParentDeltaArray();
		init_mode = parent_rep_point.getParentMode();
	}
	else
    {
		parent_point = parent_rep_point.getState();
		parent_delta_array = parent_rep_point.getDeltaArray();
		init_mode = parent_rep_point.getMode();
	}

    std::vector<double> deltas(dimensions);
    std::vector<double> b(2 * dimensions);
    std::vector<double> matrix(2 * dimensions * dimensions, 0);

	for (int modeIndex = mode.size() - 1; modeIndex >= 0; modeIndex--)
    {
        int tempMode;
		tempMode = mode.at(modeIndex);
		int foundFlag = 0;
		for (int j=0; j < modes_in_set.size(); j++)
        {
			if (tempMode == modes_in_set.at(j))
            {
				foundFlag = 1;
			}
		}
		if (foundFlag == 0)
        {
			modes_in_set.push_back(tempMode);
		}
	}

	if (modes_in_set.size() == 0)
    {
		return 0;
    }

	std::vector<int> thinvector;
	std::string thinprop;
	std::ifstream thinpropfile("../work-dir/ThinVarProp");
	thinvector.push_back(1);

	if (thinpropfile.is_open())
    {
		while (getline(thinpropfile, thinprop))
        {
			thinvector.push_back(thinprop == "1");
		}
		thinpropfile.close();
	}
	else
    { 
        std::cout << "Uable to open thin property file" << std::endl;
    }
	std::cout << "dimensions = " << dimensions << std::endl;

	for (int modeIndex = 0; modeIndex < modes_in_set.size(); modeIndex++)
    {
        int tempMode;
		tempMode = modes_in_set.at(modeIndex);

		for (int indexDimension=1; indexDimension <= dimensions; 
            indexDimension++)
        {
			double minValue = getMinCoordinate(indexDimension, tempMode);
			double maxValue = getMaxCoordinate(indexDimension, tempMode);
			if (thinvector[indexDimension] == 1)
            {
				minValue = getMaxCoordinate(indexDimension,tempMode);
            }
			deltas[indexDimension-1] = fabs(maxValue-minValue)/2;
	
            b[2 * (indexDimension - 1)] = -1 * minValue;
            b[2 * (indexDimension - 1) + 1] = maxValue;

            matrix[2 * (indexDimension - 1) * dimensions + 
                (indexDimension - 1)] = -1;
            matrix[(2 * (indexDimension - 1) + 1) * dimensions +
                (indexDimension - 1)] = 1;
		}

		LinearSet modeLinear;
		modeLinear.setDimensions(dimensions);
		modeLinear.setNumEqns(2*dimensions);
		modeLinear.setB(b);
		modeLinear.setMatrix(matrix);
		
		// temp_stack = modeLinear.getCoverStack(deltas, tempMode, refinetime);
        Simulator modeSimulator;
        modeSimulator.setInitialModeIdx(tempMode);
        modeSimulator.setInitialSet(modeLinear);
        modeSimulator.setDimensions(dimensions);
        Annotation temp_annotation;
        Verifier temp_verifier(modeSimulator, temp_annotation);
        temp_verifier.generateCoverStack();
        temp_stack = temp_verifier.getCoverStack();

		while (!temp_stack.empty())
        {
			temp_rep_point = temp_stack.top();
			temp_stack.pop();
            Point temp_point(temp_rep_point.getState());
			int curMode = tempMode;
            // FIXME 
			if (checkIntersection(curMode, temp_point, deltas)){
				double tmin = getMinTime(curMode, temp_point, deltas);
				temp_point.setCoordinate(0,tmin);
				temp_rep_point.setState(temp_point);
				temp_rep_point.setDimensions(dimensions);
				temp_rep_point.setDeltaArray(deltas);
				temp_rep_point.setParentState(parent_point);
				temp_rep_point.setParentDeltaArray(parent_delta_array);
				temp_rep_point.setParentMode(init_mode);
				temp_rep_point.setRefineTime(refinetime);
                RepPoint obj = temp_rep_point;
				itr_stack.push(obj);

				std::cout<<"====Find next region Information, recalculating deltaArray, generate one RepPoint, Push to stack===="<<std::endl;
			}
		}
	}
	return 1;
}

void ReachTube::addGuards(std::vector<std::pair<std::NNC_Polyhedron, int> > guards)
{
    isReachTube = 0;
	reachTubeMode = -1;

	int tempMode;
	std::NNC_Polyhedron poly;
	Point ptU;
    Point ptL;

	int guardthinprop;
	std::vector<int> thinvector;
	std::string thinprop;
	std::ifstream thinpropfile("../work-dir/ThinVarProp");
	thinvector.push_back(1);
	if(thinpropfile.is_open()){
		while(getline(thinpropfile,thinprop)){
			thinvector.push_back(thinprop == "1");
		}
		thinpropfile.close();
	}
	else std::cout<< "unable to open thin property file" << std::endl;

	for(int i=0; i<guards.size(); i++){
		poly = guards[i].first;
		tempMode = guards[i].second;

		Point ptL(dimensions+1);
		Point ptU(dimensions+1);
		
		for(int j=0; j<dimensions+1; j++){
			ptL.setCoordinate(j, DBL_MAX);
			ptU.setCoordinate(j, -DBL_MAX);
		}

		std::Generator_System gs=poly.minimized_generators();
		std::Generator_System::const_iterator k;
		for(k=gs.begin();k!=gs.end();++k)
		{
			if(k->is_point())
			{
			 	double divisor=mpz_get_d(k->divisor().get_mpz_t());
			  	int dim=int(k->space_dimension());
			  	for(int j=0;j<dim;j++)
			  	{
			    	double dividend=mpz_get_d(k->coefficient(std::Variable(j)).get_mpz_t());
			    	double num = dividend/divisor;

			    	if(thinvector[j] == 1){
			    		if (num > ptU.getCoordinate(j)){
			    			ptL.setCoordinate(j,num);
			    			ptU.setCoordinate(j,num);
			    		}
			    	}
			    	else{
			    		if(num < ptL.getCoordinate(j))
			    			ptL.setCoordinate(j, num);
				    	if(num > ptU.getCoordinate(j))
				    		ptU.setCoordinate(j, num);

			    	}
			  	}
			}
		}
		mode.push_back(tempMode);
		upper_bound.push_back(ptU);
		lower_bound.push_back(ptL);
	}
}

double ReachTube::getMinCoordinate(int dim, int cur_mode)
{
    double min_value = 100000;
    Point temp_point;
    
    for (int pnt_ndx = lower_bound.size()-1; pnt_ndx >= 0; pnt_ndx--)
    {
        if (mode.at(pnt_ndx) == cur_mode)
        {
            temp_point = lower_bound.at(pnt_ndx);
            if (min_value >= temp_point.getCoordinate(dim))
            {
                min_value = temp_point.getCoordinate(dim);
            }
        }
    }
    return  min_value;
}

double ReachTube::getMaxCoordinate(int dim, int cur_mode)
{
    double max_value = -100000;
    Point temp_point;

    for (int pnt_ndx = upper_bound.size()-1; pnt_ndx >= 0; pnt_ndx--)
    {
        if (mode.at(pnt_ndx) == cur_mode)
        {
            temp_point = upper_bound.at(pnt_ndx);
            if (max_value <= temp_point.getCoordinate(dim))
            {
                max_value = temp_point.getCoordinate(dim);
            }
        }
    }
    return max_value;
}

int ReachTube::checkIntersection(int cur_mode, Point cur_point, 
    std::vector<double> delta_array)
{
    int has_intersection;
    Point upper_bound_point;
    Point lower_bound_point;
    double max_dim, min_dim;
    double max_point_dim, min_point_dim;
    double delta;

    for (int i=0; i < mode.size(); i++)
    {
        if (mode.at(i) == cur_mode)
        {
            has_intersection = 1;
            upper_bound_point = upper_bound.at(i);
            lower_bound_point = lower_bound.at(i);
            for (int j=0; j < dimensions; j++)
            {
                delta = delta_array[j];
                max_dim = upper_bound_point.getCoordinate(j+1);
                min_dim = lower_bound_point.getCoordinate(j+1);
                max_point_dim = cur_point.getCoordinate(j+1) + delta;
                min_point_dim = cur_point.getCoordinate(j+1) - delta;
                if (min_dim > max_point_dim || max_dim < min_point_dim)
                {
                    has_intersection = 0;
                }
            }
            if (has_intersection == 1)
            {
                return 1;
            }
        }
    }
    return 0;
}

double ReachTube::getMinTime(int cur_mode, Point cur_point, 
    std::vector<double> delta_array)
{
    int has_intersection;
    Point upper_bound_point;
    Point lower_bound_point;
    double max_dim, min_dim;
    double max_point_dim, min_point_dim;
    double delta;
    double min_time = 1000000;

    for (int i=0; i < mode.size(); i++)
    {
        if (mode.at(i) == cur_mode)
        {
            has_intersection = 1;
            upper_bound_point = upper_bound.at(i);
            lower_bound_point = lower_bound.at(i);
            for (int j=0; j < dimensions; j++)
            {
                delta = delta_array[j];
                max_dim = upper_bound_point.getCoordinate(j+1);
                min_dim = lower_bound_point.getCoordinate(j+1);
                max_point_dim = cur_point.getCoordinate(j+1) + delta;
                min_point_dim = cur_point.getCoordinate(j+1) - delta;
                if (min_dim > max_point_dim || max_dim < min_point_dim)
                {
                    has_intersection = 0;
                }
            }
            if (has_intersection == 1)
            {
                if (min_time > lower_bound_point.getCoordinate(0))
                {
                    min_time = lower_bound_point.getCoordinate(0);
                }
            }
        }
    }
    return min_time;
}
        

//////////////////// Getters and Setters ////////////////////

void ReachTube::setDimensions(int val)
{
	dimensions = val;
}

int ReachTube::getDimensions()
{
	return dimensions;
}

void ReachTube::setMode(int val)
{
	isReachTube = 1;
	reachTubeMode = val;
}

int ReachTube::getMode()
{
	return reachTubeMode;
}

std::vector<int> ReachTube::getModeVec()
{
    return mode;
}

void ReachTube::setModeVec(std::vector<int> vec)
{
    mode = vec;
}

void ReachTube::addLowerBoundState(Point obj)
{
    lower_bound.push_back(obj);
}
void ReachTube::addUpperBoundState(Point obj)
{
    upper_bound.push_back(obj);
}