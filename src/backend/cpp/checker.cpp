/*
 * \file   checker.cpp
 * \class  Checker
 * 
 * \author parasara
 * \author Lucas Brown
 * \date   Feb 2, 2013
 * \date   April 1, 2019
 * 
 * \brief LMBTODO 
 */

#include "checker.hpp"
#include "point.hpp"

Checker::Checker()
{ }

Checker::~Checker()
{ }

int Checker::checkHybridSimulation(ReachTube simulation_tube, 
    LinearSet unsafe_set)
{
	Point upper_pt;
	int size = simulation_tube.getSize();
	for (int i=0; i < size; i++)
    {
		upper_pt = simulation_tube.getUpperBound(i);
        // LMBTODO Why do we only check the upper point?
		if (unsafe_set.isInternal(upper_pt))
        {
			return -1;
		}
	}
	return 1;
}

int Checker::check(ReachTube invariant_tube, LinearSet unsafe_set)
{
    // LMBTODO: Need to ask about this 
	int is_safe = 1;

	int size_reach_set = invariant_tube.getSize();
	for (int i=0; i < size_reach_set; i++)
    {
		Point upper_point = invariant_tube.getUpperBound(i);
		Point lower_point = invariant_tube.getLowerBound(i);

		if (unsafe_set.hasIntersection(lower_point, upper_point) == 1)
        {
			// Check for 2^n possibilities!
            // LMBTODO Reveiw the above comment - what is n? Is it accurate?
			if (unsafe_set.isInternal(upper_point) && 
                unsafe_set.isInternal(lower_point))
            {
				is_safe = -1;
                break;
			}
			else
            {
				is_safe = 0;
			}
		}
	}
	return is_safe;
}