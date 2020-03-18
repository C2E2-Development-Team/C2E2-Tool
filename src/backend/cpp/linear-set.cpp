/*
 * \file   linear-set.cpp
 * \class  LinearSet
 * 
 * \author parasara
 * \author Lucas Brown
 * \date   July 28, 2014
 * \date   April 1, 2019
 * 
 * \brief LMBTODO 
 */

#include <iostream>
#include <cstdlib>
#include <glpk.h>
#include <vector>

#include "linear-set.hpp"

LinearSet::LinearSet() 
{ }

LinearSet::~LinearSet() 
{ }

int LinearSet::isInternal(Point point)
{
	glp_prob* feas;
	feas = glp_create_prob();

	int* irow;
    int* icol;
	double* icoeffs;
	int row_col_size = 1 + ((number_of_equations+1) * (dimensions+1));
    
    irow = (int*)malloc(row_col_size * sizeof(int));
	icol = (int*)malloc(row_col_size * sizeof(int));
	icoeffs = (double*)malloc(row_col_size * sizeof(double));
	for (int i=0; i < row_col_size; i++)
    {
        irow[i] = 0; 
        icol[i] = 0; 
        icoeffs[i] = 0;
	}

	glp_set_prob_name(feas, "Feasibility checker");
	glp_set_obj_dir(feas, GLP_MAX);
	glp_add_rows(feas, number_of_equations+1);
	glp_add_cols(feas, dimensions+1);

	for (int i=1; i <= dimensions; i++)
    {
        // LMBTODO We shouldn't have to do this, or at least do it without 
        //         having to evaluate an if statement every iteration. What
        //         do we expect to happen if point.getDimensions isn't either 
        //         equal to dimensions or dimensions + 1?
		if (point.getDimensions() == dimensions)
        {
			glp_set_col_bnds(feas, i, GLP_FX, point.getCoordinate(i-1), 
                point.getCoordinate(i-1));
		}
		else if (point.getDimensions() == dimensions+1)
        {
			glp_set_col_bnds(feas, i, GLP_FX, point.getCoordinate(i),
                point.getCoordinate(i));
		}
	}
	for (int i=1; i <= number_of_equations; i++)
    {
		glp_set_row_bnds(feas, i, GLP_UP, -10000, getBElement(i-1));
	}

	int count = 0;
	for (int i=0; i < number_of_equations; i++)
    {
		for (int j=0; j < dimensions; j++)
        {
			irow[count+1] = i+1; 
            icol[count+1] = j+1;
			icoeffs[count+1] = getMatrixElement(i, j);
			count++;
		}
	}

	glp_load_matrix(feas, count, irow, icol, icoeffs);
	glp_set_obj_coef(feas, 1, 1.0);
	glp_term_out(GLP_OFF);
	glp_simplex(feas, NULL);

	int status = glp_get_status(feas);
    glp_delete_prob(feas);
    free(irow);
    free(icol);
    free(icoeffs);
	if (status == GLP_INFEAS || status == GLP_NOFEAS || status == GLP_UNDEF)
    {
		return 0;
	}
	else
    {
		return 1;
	}
}

double LinearSet::getMax(int dimension)
{
    return findMaxMin(dimension, true);	
}

double LinearSet::getMin(int dimension)
{
	return findMaxMin(dimension, false);
}

int LinearSet::hasIntersection(Point lower_point, Point upper_point)
{
	glp_prob* feas;
	feas = glp_create_prob();

	int* irow;
    int* icol;
	double* icoeffs;
    int row_col_size = 1 + (number_of_equations + 1) * (dimensions + 1);

	irow = (int*)malloc(row_col_size * sizeof(int));
	icol = (int*)malloc(row_col_size * sizeof(int));
	icoeffs = (double*)malloc(row_col_size * sizeof(double));
	for(int i=0; i < row_col_size; i++)
    {
		irow[i] = 0;
        icol[i] = 0;
        icoeffs[i] = 0;
	}

	glp_set_prob_name(feas, "Feasibility checker");
	glp_set_obj_dir(feas, GLP_MAX);
	glp_add_rows(feas, number_of_equations+1);
	glp_add_cols(feas, dimensions+1);

	for(int i=1; i <= dimensions; i++)
    {
        // LMBTODO We shouldn't have to do this, or at least do it without 
        //         having to evaluate an if statement every iteration. What
        //         do we expect to happen if point.getDimensions isn't either 
        //         equal to dimensions or dimensions + 1?
		if(lower_point.getDimensions() == dimensions)
        {
			glp_set_col_bnds(feas, i, GLP_DB, lower_point.getCoordinate(i-1),
                upper_point.getCoordinate(i-1));
		}
		else if(lower_point.getDimensions() == dimensions+1)
        {
			glp_set_col_bnds(feas, i, GLP_DB, lower_point.getCoordinate(i-1),
                upper_point.getCoordinate(i));
		}
	}
	for(int i=1; i <= number_of_equations; i++)
    {
		glp_set_row_bnds(feas, i, GLP_UP, -10000, getBElement(i-1));
	}

	int count = 0;
	for (int i=0; i < number_of_equations; i++)
    {
		for (int j=0; j < dimensions; j++)
        {
			irow[count+1] = i+1; 
            icol[count+1] = j+1;
			icoeffs[count+1] = getMatrixElement(i, j);
			count++;
		}
	}

	glp_load_matrix(feas, count, irow, icol, icoeffs);
	glp_set_obj_coef(feas, 1, 1.0);
	glp_term_out(GLP_OFF);
	glp_simplex(feas, NULL);

    int status = glp_get_status(feas);
    glp_delete_prob(feas);
	free(irow);
    free(icol);
    free(icoeffs);
	if (status == GLP_INFEAS || status == GLP_NOFEAS || status == GLP_UNDEF)
    {
		return 0;
	}
	else
    {
		return 1;
	}
}

void LinearSet::print()
{
    std::cout << "Printing LinearSet...\n";
    std::cout << "  Dimensions: " << dimensions << "\n";
    std::cout << "  number_of_equations: " << number_of_equations << "\n";
    std::cout << "  Matrix: \n    ";
    for (auto v: matrix)
    {
        std::cout << v << " ";
    }
    std::cout << "\n  B: \n    ";
    for (auto v: b)
    {
        std::cout << v << " ";
    }
    std::cout << std::endl;
}

//////////////////// Getters and Setters ////////////////////

int LinearSet::getDimensions()
{
	return dimensions;
}
void LinearSet::setDimensions(int val)
{
	dimensions = val;
}

int LinearSet::getNumEqns()
{
	return number_of_equations;
}
void LinearSet::setNumEqns(int val)
{
	number_of_equations = val;
}

std::vector<double> LinearSet::getMatrix()
{
    return matrix;
}
void LinearSet::setMatrix(std::vector<double> vec)
{
    matrix = vec;
}
double LinearSet::getMatrixElement(int row, int col)
{
	return matrix[row*dimensions + col];
}
void LinearSet::setMatrixElement(int row, int col, double value)
{
    matrix[row*dimensions + col] = value;
}

std::vector<double> LinearSet::getB()
{
    return b;
}
void LinearSet::setB(std::vector<double> vec)
{
    b = vec;
}
double LinearSet::getBElement(int index)
{
	return b[index];
}
void LinearSet::setBElement(int index, double value)
{
    b[index] = value;
}

//////////////////// Private Functions ////////////////////

double LinearSet::findMaxMin(int dim, bool max)
{
    glp_prob* feas;
	feas = glp_create_prob();

	int* irow;
    int* icol;
	double* icoeffs;
    int row_col_size = 1 + ((number_of_equations+1) * (dimensions+1));

	irow = (int*)malloc(row_col_size * sizeof(int));
	icol = (int*)malloc(row_col_size * sizeof(int));
	icoeffs = (double*)malloc(row_col_size * sizeof(double));
	for (int i=0; i < row_col_size; i++)
    {
        irow[i] = 0; 
        icol[i] = 0; 
        icoeffs[i] = 0;
	}

	glp_set_prob_name(feas, "Feasibility checker");
    if (max)
    {
        glp_set_obj_dir(feas, GLP_MAX);
    }
    else
    {
        glp_set_obj_dir(feas, GLP_MIN);
    }
	glp_add_rows(feas, number_of_equations+1);
	glp_add_cols(feas, dimensions+1);

	for (int i=1; i <= dimensions; i++)
    {
		glp_set_col_bnds(feas, i, GLP_FR, -1000000, 1000000);
	}
	for (int i=1; i <= number_of_equations; i++)
    {
		glp_set_row_bnds(feas, i, GLP_UP, -10000, getBElement(i-1));
	}

	int count = 0;
	for (int i=1; i <= number_of_equations; i++)
    {
		for (int j=1; j <= dimensions; j++)
        {
			irow[count+1] = i; 
            icol[count+1] = j;
			icoeffs[count+1] = getMatrixElement(i-1, j-1);
			count++;
		}
	}

	glp_load_matrix(feas, count, irow, icol, icoeffs);
	glp_set_obj_coef(feas, dim+1, 1.0);
	glp_term_out(GLP_OFF);
	glp_simplex(feas, NULL);

	int status = glp_get_status(feas);
    glp_delete_prob(feas);
    free(irow);
    free(icol);
    free(icoeffs);
	if (status == GLP_INFEAS || status == GLP_NOFEAS || status == GLP_UNDEF)
    {
		//std::cout << "Says that its infeasible! " <<  std::endl;
		return 0;
	}
	else
    {
		return glp_get_obj_val(feas);
	}
}