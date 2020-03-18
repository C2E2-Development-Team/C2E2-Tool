/**
 * \file   simulator.cpp
 * \class  Simulator
 *
 * \author parasara
 * \author Lucas Brown
 * \date   Feb 1, 2013
 * \date   April 1, 2019
 * 
 * \brief  LMBTODO
 */

#include <algorithm>
#include <cmath>
#include <dlfcn.h>
#include <fstream>
#include <ppl.hh>
#include <Python.h>
#include <stack>
#include <string>
#include <vector>

#include "rep-point.hpp"
#include "simulator.hpp"

Simulator::Simulator()
{
}

Simulator::~Simulator()
{
}

int Simulator::simulate()
{
    return hybridSimulation();
}

int Simulator::hybridSimulation()
{
    std::ofstream invFile;            // LMBTODO Need a better filename
    invFile.open(visualize_filename); //LMBTODO Need a better name for this too
    invFile.close();

    std::vector<Point> pts = getRepresentativeCover();

    int is_safe = -1;
    for (int i = 0; i < pts.size(); i++)
    {
        std::cout << "Hybrid Simulation " << i + 1 << " -> Point: \n";
        pts[i].print();

        is_safe = hybridSimulation(pts[i]);

        if (is_safe == -1)
        {
            std::cout << "Hybrid Simulation " << i + 1 << " unsafe" << std::endl;
            break;
        }
        else if (is_safe == 1)
        {
            std::cout << "Hybrid Simulation " << i + 1 << " safe" << std::endl;
        }

        // LMBTODO What do we do in the event the simulation is not -1 or 1?
    }

    std::ofstream resultStream;
    resultStream.open("../work-dir/Result.dat");
    resultStream << is_safe << std::endl;
    resultStream.close();

    return is_safe;
}

int Simulator::hybridSimulation(Point origin)
{
    void *lib = dlopen("../work-dir/libhybridsim.so", RTLD_LAZY);
    // LMBTODO handle failed library opening.

    typedef std::vector<std::pair<std::NNC_Polyhedron, int>> (*guard_fn)(int, double *, double *);
    typedef bool (*inv_fn)(int, double *, double *);
    guard_fn guards = (guard_fn)dlsym(lib, "hitsGuard");
    inv_fn invs = (inv_fn)dlsym(lib, "invariantSatisfied");

    int mode = initial_mode_idx;
    int is_safe = -1;
    double *lower_point = NULL;
    double *upper_point = NULL;
    while (true)
    {
        std::cout << "Simulating Mode " << mode << " from Point:" << std::endl;
        origin.print();
        simulatePoint(origin, mode);
        std::cout << "Simualtion complete, performing analysis." << std::endl;

        ReachTube simulation_tube;
        simulation_tube.setDimensions(dimensions);
        simulation_tube.setMode(mode);
        simulation_tube.parseInvariantTube("../work-dir/SimuOutput", 0);

        int size = simulation_tube.getSize();
        if (size == 0)
            break;

        std::vector<std::pair<std::NNC_Polyhedron, int>> guards_hit;
        for (int i = 0; i < size; i++)
        {
            std::vector<double> lower =
                simulation_tube.getLowerBound(i).getCoordinates(); // LMBTODO: Need to make sure this returns a double
            std::vector<double> upper =
                simulation_tube.getUpperBound(i).getCoordinates();
            lower_point = &lower[0]; // Valid for C++11 and on
            upper_point = &upper[0];
            guards_hit = guards(mode, lower_point, upper_point);

            if (!guards_hit.empty())
            {
                std::pair<std::NNC_Polyhedron, int> guard_taken =
                    guards_hit[rand() % guards_hit.size()];
                mode = guard_taken.second;
                origin = getPointFromPoly(guard_taken.first);
                simulation_tube.clear(i + 1);
                break;
            }
            if (!invs(mode, lower_point, upper_point))
            {
                simulation_tube.clear(i);
                break;
            }
        }

        // LMBTODO: Why use trace_safe instead of using is_safe?
        int trace_safe;
        trace_safe = checker.checkHybridSimulation(simulation_tube, unsafe_set);
        if (trace_safe == 1)
        {
            simulation_tube.printReachTube(visualize_filename, 1);
            is_safe = 1;
        }
        else if (trace_safe == -1)
        {
            simulation_tube.printReachTube(visualize_filename, 2);
            is_safe = -1;
            break;
        }
        else
        {
            std::cout << "<ERROR> UNKNOWN TUBE IN HYBRID SYSTEM" << std::endl;
        }

        if (guards_hit.empty())
        {
            break;
        }
    }
    lower_point = NULL;
    upper_point = NULL;
    dlclose(lib);

    std::cout << "is_safe: " << is_safe << std::endl;
    return is_safe;
}

void Simulator::simulatePoint(Point point, int mode)
{
    std::ofstream simulation_config;
    simulation_config.open("../work-dir/Config");

    for (int i = 0; i < point.getDimensions(); i++)
    {
        simulation_config << point.getCoordinate(i) << "\n";
    }
    simulation_config << absolute_error << "\n";
    simulation_config << relative_error << "\n";
    simulation_config << time_step << "\n";
    simulation_config << time_horizon << "\n";
    simulation_config << mode << "\n";
    simulation_config.close();

    std::string exec_command = "./" + executable + " < ../work-dir/Config > " +
                               "../work-dir/SimuOutput";
    std::cout << "@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
    std::cout << exec_command << std::endl;
    system(exec_command.c_str());
}

void Simulator::simulatePoint(Point point, int mode, int idx)
{
    std::ofstream simulation_config;
    simulation_config.open("../work-dir/Config");

    for (int i = 0; i < point.getDimensions(); i++)
    {
        simulation_config << point.getCoordinate(i) << "\n";
    }
    simulation_config << absolute_error << "\n";
    simulation_config << relative_error << "\n";
    simulation_config << time_step << "\n";
    simulation_config << time_horizon << "\n";
    simulation_config << mode << "\n";
    simulation_config.close();

    std::string exec_command = "./" + executable + " < ../work-dir/Config > " +
                               "../work-dir/SimuOutput" + std::to_string(idx);
    std::cout << "@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
    std::cout << exec_command << std::endl;
    system(exec_command.c_str());
}

std::vector<Point> Simulator::getRepresentativeCover() // LMBTODO Consider renaming
{
    int n = 3; // LMBTODO: Move to config
    std::vector<double> minimums(dimensions);
    std::vector<double> maximums(dimensions);
    for (int i = 0; i < dimensions; i++)
    {
        minimums[i] = initial_set.getMin(i);
        maximums[i] = initial_set.getMax(i);
    }

    // Count all dimensions that have thickness
    int thick_dimensions = 0;
    for (int i = 0; i < dimensions; i++)
    {
        if (minimums[i] != maximums[i])
        {
            thick_dimensions++;
        }
    }

    // Initialize cover_points vector
    int num_points = (int)pow(n, thick_dimensions); // pow(base, exponent)
    std::vector<Point> cover(num_points);
    for (int i = 0; i < num_points; i++)
    {
        cover[i] = Point(dimensions + 1); // 0th dimensions reserved for time.
    }

    // Create Cover
    for (int dim = 0, thick_dim_count = 0; dim < dimensions; dim++)
    {
        if (minimums[dim] == maximums[dim])
        {
            for (int i = 0; i < num_points; i++)
            {
                cover[i].setCoordinate(dim + 1, minimums[dim]);
            }
        }
        else
        {
            double start = minimums[dim];
            double step_size = (maximums[dim] - minimums[dim]) / (n - 1);
            int block_size = (int)pow(n, thick_dim_count);
            int i = 0;
            int j = 0;
            while (i < num_points)
            {
                double val = start + (j % n) * step_size;
                for (int k = 0; k < block_size; k++)
                {
                    cover[i].setCoordinate(dim + 1, val);
                    i++;
                }
                j++;
            }
            thick_dim_count++;
        }
    }

    // Print Cover
    std::cout << "\n-----------------------------\n"
              << "Representative Cover\n"
              << "-----------------------------" << std::endl;
    for (auto pt : cover)
    {
        pt.print();
    }
    std::cout << "-----------------------------" << std::endl;

    return cover;
}

// LMBTOD: This function has not been reviewed.
Point Simulator::getPointFromPoly(std::NNC_Polyhedron poly)
{
    Point point(dimensions + 1);

    std::Generator_System gs = poly.minimized_generators();
    std::Generator_System::const_iterator k;
    for (k = gs.begin(); k != gs.end(); ++k)
    {
        if (k->is_point())
        {
            double divisor = mpz_get_d(k->divisor().get_mpz_t());
            int dim = int(k->space_dimension());
            for (int j = 0; j < dim; j++)
            {
                double dividend =
                    mpz_get_d(k->coefficient(std::Variable(j)).get_mpz_t());
                double num = dividend / divisor;
                point.setCoordinate(j, num);
            }
        }
    }

    return point;
}

void Simulator::print()
{
    std::cout << "\n-----------------------------\n"
              << "Simulator Parameters\n"
              << "-----------------------------" << std::endl;
    std::cout << "Dimensions: " << dimensions << "\n"
              << "Initial Mode Index: " << initial_mode_idx << "\n"
              << "Time Step: " << time_step << "\n"
              << "Time Horizon: " << time_horizon << "\n"
              << "Absolute Error: " << absolute_error << "\n"
              << "Relative Error: " << relative_error << "\n"
              << "Visualize Filename: " << visualize_filename << "\n"
              << "Executable: " << executable << std::endl;
}

//////////////////// Getters and Setters ////////////////////

int Simulator::getDimensions()
{
    return dimensions;
}
void Simulator::setDimensions(int val)
{
    dimensions = val;
}
int Simulator::getInitialModeIdx()
{
    return initial_mode_idx;
}
void Simulator::setInitialModeIdx(int val)
{
    initial_mode_idx = val;
}

double Simulator::getAbsoluteError()
{
    return absolute_error;
}
void Simulator::setAbsoluteError(double val)
{
    absolute_error = val;
}
double Simulator::getRelativeError()
{
    return relative_error;
}
void Simulator::setRelativeError(double val)
{
    relative_error = val;
}
double Simulator::getTimeStep()
{
    return time_step;
}
void Simulator::setTimeStep(double val)
{
    time_step = val;
}
double Simulator::getTimeHorizon()
{
    return time_horizon;
}
void Simulator::setTimeHorizon(double val)
{
    time_horizon = val;
}

std::string Simulator::getVisualizeFilename()
{
    return visualize_filename;
}
void Simulator::setVisualizeFilename(std::string str)
{
    visualize_filename = str;
}
std::string Simulator::getExecutable()
{
    return executable;
}
void Simulator::setExecutable(std::string str)
{
    executable = str;
}

LinearSet Simulator::getInitialSet()
{
    return initial_set;
}
void Simulator::setInitialSet(LinearSet obj)
{
    initial_set = obj;
}
LinearSet Simulator::getUnsafeSet()
{
    return unsafe_set;
}
void Simulator::setUnsafeSet(LinearSet obj)
{
    unsafe_set = obj;
}

// LMBTODO: Band-aid functions
std::vector<int> Simulator::getModeLinear()
{
    return mode_linear;
}
void Simulator::setModeLinear(std::vector<int> vec)
{
    mode_linear = vec;
}
std::vector<double> Simulator::getKConsts()
{
    return k_consts;
}
void Simulator::setKConsts(std::vector<double> vec)
{
    k_consts = vec;
}

Checker Simulator::getChecker()
{
    return checker;
}