/*
 * \file   checker.hpp
 * \class  Checker
 * 
 * \author parasara
 * \author Lucas Brown
 * \date   Feb 2, 2013
 * \date   April 1, 2019
 * 
 * \brief LMBTODO 
 */

#ifndef CHECKER_H_
#define CHECKER_H_

#include "linear-set.hpp"
#include "reach-tube.hpp"

class Checker 
{
  public:
	Checker();
	~Checker();

    // LMBTODO Both of these classes need more descriptive names
    int checkHybridSimulation(ReachTube simulation_tube, LinearSet unsafe_set);
    int check(ReachTube invariant_tube, LinearSet unsafe_set);

  private:

};

#endif /* CHECKER_H_ */