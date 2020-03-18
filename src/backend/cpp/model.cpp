/**
 * \file   model.cpp
 * \class  Model
 *
 * \author Matt Potok
 * \author Lucas Brown
 * \date   April 1, 2019
 * 
 * \brief  The Model data structure contains all information pertinent of a 
 * particular model.
 */

#include <ctime>
#include <iostream>
#include <Python.h>
#include <string>
#include <vector>

#include "model.hpp"
#include "simulator.hpp"
#include "verifier.hpp"

Model::Model() 
{ }

Model::~Model() 
{ } 

int Model::simulateVerify()
{
    std::time_t start_time = std::time(NULL);
    // freopen("../work-dir/log.txt", "w", stdout);

    initialize();
    printModel();

    Simulator simulator;
    simulator.setDimensions(dimensions);
    simulator.setInitialModeIdx(initial_mode_idx);
    simulator.setInitialSet(initial_set);
    simulator.setUnsafeSet(unsafe_set);
    simulator.setVisualizeFilename(visualize_filename);  // Rename
    simulator.setAbsoluteError(abs_error);  
    simulator.setRelativeError(rel_error); 
    simulator.setTimeStep(time_step);
    simulator.setTimeHorizon(time_horizon);
    simulator.setExecutable(executable);

    simulator.print();

    int safety;
    if (is_simulation)
    {
        safety = simulator.simulate();
    }
    else
    {
        Annotation annotation;
        annotation.setType(annot_type);
        annotation.setKConsts(k_consts);
        annotation.setGammas(gammas);

        Verifier verifier(simulator, annotation);
        verifier.setNumUnsafeEqns(num_unsafe_eqns);  // LMBTODO: Band-aid
        verifier.setModeLinear(mode_linear); // LMBTODO: Band-aid
        verifier.setKConsts(k_consts);  // LMBTODO: Band-aid
        verifier.setUnsafeMatrix(unsafe_matrix);  // LMBTODO: Band-aid

        verifier.setRefineStrat(refine_strat);
        verifier.setRefineOrder(refine_order);

        safety = verifier.verify();
    }

    std::cout << "Execution Time: "
        << std::difftime(std::time(NULL), start_time) << " seconds." 
        << std::endl;
    return safety;
}

void Model::initialize()
{
    initial_set.setDimensions(dimensions);
    initial_set.setNumEqns(num_initial_eqns);
    initial_set.setMatrix(initial_matrix);
    initial_set.setB(initial_b);

    unsafe_set.setDimensions(dimensions);
    unsafe_set.setNumEqns(num_unsafe_eqns);
    unsafe_set.setB(unsafe_b);
    unsafe_set.setMatrix(unsafe_matrix);

    Py_Initialize();
}

void Model::printModel() 
{
    std::cout << "\n-----------------------------\n"
              << "Model Parameters\n"
              << "-----------------------------" << std::endl;
    
    std::cout << "is_simulation: " << is_simulation << "\n"
              << "refine_strat: " << refine_strat << "\n" 
              << "dimensions: " << dimensions << "\n"
              << "num_modes: " << num_modes << "\n"
              << "initial_mode_idx" << initial_mode_idx << "\n"
              << "num_initial_eqns: " << num_initial_eqns << "\n"
              << "num_unsafe_eqns: " << num_unsafe_eqns << "\n"
              << "annot_type: " << annot_type << "\n"
              << std::endl;

    std::cout << "abs_error: " << abs_error << "\n"
              << "rel_error: " << rel_error << "\n"
              << "time_step: " << time_step << "\n"
              << "time_horizon: " << time_horizon << "\n"
              << std::endl;

    print_vector<int>("mode_linear: ", mode_linear);
    
    print_vector<double>("gammas: ", gammas);
    print_vector<double>("k_consts: ", k_consts);
    print_vector<double>("initial_matrix: ", initial_matrix);
    print_vector<double>("initial_b: ", initial_b);
    print_vector<double>("unsafe_matrix: ", unsafe_matrix);
    print_vector<double>("unsafe_b: ", unsafe_b);

    std::cout << "annot_str: " << annot_str << "\n"
              << "beta_str: " << beta_str << "\n"
              << "opt_str: " << opt_str << "\n"
              << "visualize_filename: " << visualize_filename << "\n"
              << std::endl;

    std::cout << "-----------------------------" << std::endl;
}

//////////////////// Getters and Setters ////////////////////

////////// Bool and Int //////////
bool Model::isSimulation() 
{
    return is_simulation;
}
void Model::setSimulationBool(bool val)
{
    is_simulation = val;
}
int Model::getRefineStrat()
{
    return refine_strat;
}
void Model::setRefineStrat(int val)
{
    refine_strat = val;
}
int Model::getDimensions()
{
    return dimensions;
}
void Model::setDimensions(int val)
{
    dimensions = val;
}
int Model::getNumModes()
{
    return num_modes;
}
void Model::setNumModes(int val)
{
    num_modes = val;
}
int Model::getInitialModeIdx()
{
    return initial_mode_idx;
}
void Model::setInitialModeIdx(int val)
{
    initial_mode_idx = val;
}
int Model::getNumInitialEqns()
{
    return num_initial_eqns;
}
void Model::setNumInitialEqns(int val)
{
    num_initial_eqns = val;
}
int Model::getNumUnsafeEqns()
{
    return num_unsafe_eqns;
}
void Model::setNumUnsafeEqns(int val)
{
    num_unsafe_eqns = val;
}
int Model::getAnnotType()
{
    return annot_type;
}
void Model::setAnnotType(int val)
{
    annot_type = val;
}

////////// Double //////////
double Model::getAbsError()
{
    return abs_error;
}
void Model::setAbsError(double val)
{
    abs_error = val;
}
double Model::getRelError()
{
    return rel_error;
}
void Model::setRelError(double val)
{
    rel_error = val;
}
double Model::getTimeStep()
{
    return time_step;
}
void Model::setTimeStep(double val)
{
    time_step = val;
}
double Model::getTimeHorizon()
{
    return time_horizon;
}
void Model::setTimeHorizon(double val)
{
    time_horizon = val;
}

////////// Int Vector //////////
std::vector<int> Model::getModeLinear()
{
    return mode_linear;
}
void Model::setModeLinear(std::vector<int> vec)
{
    mode_linear = vec;
}

void Model::setRefineOrder(std::vector<int> vec)
{
    refine_order = vec;
}

////////// Double Vector /////////////
std::vector<double> Model::getGammas()
{
    return gammas;
}
void Model::setGammas(std::vector<double> vec)
{
    gammas = vec;
}
std::vector<double> Model::getKConsts()
{
    return k_consts;
}
void Model::setKConsts(std::vector<double> vec)
{
    k_consts = vec;
}
std::vector<double> Model::getInitialMatrix()
{
    return initial_matrix;
}
void Model::setInitialMatrix(std::vector<double> vec)
{
    initial_matrix = vec;
}
std::vector<double> Model::getInitialB()
{
    return initial_b;
}
void Model::setInitialB(std::vector<double> vec)
{
    initial_b = vec;
}
std::vector<double> Model::getUnsafeMatrix()
{
    return unsafe_matrix;
}
void Model::setUnsafeMatrix(std::vector<double> vec)
{
    unsafe_matrix = vec;
}
std::vector<double> Model::getUnsafeB()
{
    return unsafe_b;
}
void Model::setUnsafeB(std::vector<double> vec)
{
    unsafe_b = vec;
}

////////// String //////////
std::string Model::getAnnotStr()
{
    return annot_str;
}
void Model::setAnnotStr(std::string val)
{
    annot_str = val;
}
std::string Model::getBetaStr()
{
    return beta_str;
}
void Model::setBetaStr(std::string val)
{
    beta_str = val;
}
std::string Model::getOptStr()
{
    return opt_str;
}
void Model::setOptStr(std::string val)
{
    opt_str = val;
}
std::string Model::getVisualizeFilenmae()
{
    return visualize_filename;
}
void Model::setVisualizeFilename(std::string val)
{
    visualize_filename = val;
}
std::string Model::getExecutable()
{
    return executable;
}
void Model::setExecutable(std::string val)
{
    executable = val;
}