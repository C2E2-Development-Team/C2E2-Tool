/*
 * \file   annotation.cpp
 * \class  Annotation
 * 
 * \author parasara
 * \author Lucas Brown
 * \date   Feb 1, 2013
 * \date   March 29, 2019
 * 
 * \brief LMBTODO 
 */

#include <iostream>
#include <math.h>

#include "annotation.hpp"

Annotation::Annotation() 
{ 
    type = -1;
}

Annotation::~Annotation() 
{ }

double Annotation::getED(double delta, double t1, double t2, int mode)
{
	double delta1;
	double delta2;
	double finRet;  // LMBTODO Rename
	
    double k1 = k_consts[mode-1];  // LMBTODO Rename
    double gamma1 = gammas[mode-1];  // LMBTODO Rename

	if (type == 1)
    {
		delta1 = 1.1*k1*gamma1*t1*delta;
		delta2 = 1.1*k1*gamma1*t2*delta;
		if (delta1 > delta2)
        {
			finRet = delta1;
		}
		else
        {
			finRet = delta2;
		}

		if (finRet < 0)
        {
			finRet = -1*finRet;
		}
		else if (finRet == 0)
        {
			finRet = delta;
		}

		return finRet;
	}
	if (type ==2 || type == 3)
    {
		delta1 = 1.1*k1*exp(gamma1*t1)*delta;
		delta2 = 1.1*k1*exp(gamma1*t2)*delta;

		if(delta1 > delta2)
        {
			finRet = delta1;
		}
		else
        {
			finRet = delta2;
		}

		if (finRet < 0)
        {
			finRet = -1*finRet;
		}
		else if (finRet == 0)
        {
			finRet = delta;
		}

		return finRet;
	}
	else
    {
		return 1.1*delta;
	}
}

//////////////////// Getters and Setters ////////////////////

int Annotation::getType()
{
	return type;
}
void Annotation::setType(int val)
{
	type = val;
}

std::string Annotation::getAnnotation()
{
	return Annot;
}
void Annotation::setAnnotation(std::string val)
{
	Annot = val;
}
std::string Annotation::getBeta()
{
	return Beta;
}
void Annotation::setBeta(std::string val)
{
	Beta = val;
}

double Annotation::getKVal(int mode)
{
    return k_consts[mode-1];
}
void Annotation::setKVal(int mode, double val)
{
    k_consts[mode-1] = val;
}
std::vector<double> Annotation::getKConsts()
{
    return k_consts;
}
void Annotation::setKConsts(std::vector<double> vec)
{
    k_consts = vec;
}

double Annotation::getGammaVal(int mode)
{
    return gammas[mode-1];
}
void Annotation::setGammaVal(int mode, double val)
{
    gammas[mode-1] = val;
}
std::vector<double> Annotation::getGammas()
{
    return gammas;
}
void Annotation::setGammas(std::vector<double> vec)
{
    gammas = vec;
}