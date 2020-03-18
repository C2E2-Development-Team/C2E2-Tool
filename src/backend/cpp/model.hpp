/**
 * \file   model.hpp
 * \class  Model
 * 
 * \author Matt Potok
 * \author Lucas Brown
 * \date   April 1, 2019
 * 
 * \brief  The Model data structure contains all information pertinent of a 
 * particular model.
 */

#ifndef MODEL_H
#define MODEL_H

#include <vector>
#include <iostream>
#include <string>

#include "linear-set.hpp"

class Model
{
  public:
    Model();
    ~Model();

    int simulateVerify();
    void initialize();
    void printModel();

    bool isSimulation();
    void setSimulationBool(bool val);
    int getRefineStrat();
    void setRefineStrat(int val);
    int getDimensions();
    void setDimensions(int val);
    int getNumModes();
    void setNumModes(int val);
    int getInitialModeIdx();
    void setInitialModeIdx(int val);
    int getNumInitialEqns();
    void setNumInitialEqns(int val);
    int getNumUnsafeEqns();
    void setNumUnsafeEqns(int val);
    int getAnnotType();  // LMBTODO: Add to config, move to Simulator
    void setAnnotType(int val);  // LMBTODO: Add to config, move to Simulator

    double getAbsError();  // LMBTODO: Add to config, move to Simulator
    void setAbsError(double val);  // LMBTODO: Add to config, move to Simulator
    double getRelError();  // LMBTODO: Add to config, move to Simulator
    void setRelError(double val);  // LMBTODO: Add to config, move to Simulator
    double getTimeStep();
    void setTimeStep(double val);
    double getTimeHorizon();
    void setTimeHorizon(double val);

    std::vector<int> getModeLinear();
    void setModeLinear(std::vector<int> vec);

    void setRefineOrder(std::vector<int> vec);
    
    std::vector<double> getGammas();
    void setGammas(std::vector<double> vec);
    std::vector<double> getKConsts();
    void setKConsts(std::vector<double> vec);
    std::vector<double> getInitialMatrix();
    void setInitialMatrix(std::vector<double> vec);
    std::vector<double> getInitialB();
    void setInitialB(std::vector<double> vec);
    std::vector<double> getUnsafeMatrix();
    void setUnsafeMatrix(std::vector<double> vec);
    std::vector<double> getUnsafeB();
    void setUnsafeB(std::vector<double> vec);

    std::string getAnnotStr(); // LMBTODO: Add to config, move to Simulator???
    void setAnnotStr(std::string val);  // LMBTODO: Add to config, move to Simulator???
    std::string getBetaStr();  // LMBTODO: Add to config, move to Simulator???
    void setBetaStr(std::string val);  // LMBTODO: Add to config, move to Simulator???
    std::string getOptStr();  // LMBTODO: Add to config, move to Simulator???
    void setOptStr(std::string val);  // LMBTODO: Add to config, move to Simulator???
    std::string getVisualizeFilenmae();
    void setVisualizeFilename(std::string val);
    std::string getExecutable();
    void setExecutable(std::string val);
    
  private:

    bool is_simulation;
    int refine_strat; /** TODOLMB Refine strategy? */
    int dimensions;  /** Number of local variables */
    int num_modes;
    int initial_mode_idx;  // Previously initial_mode
    int num_initial_eqns; 
    int num_unsafe_eqns;
    int annot_type;  // LMBTODO: Add to config, move to Simulator

    double abs_error;  // LMBTODO: Add to config, move to Simulator
    double rel_error; // LMBTODO: Add to config, move to Simulator
    double time_step;  // Previously "delta_time" in C++ code
    double time_horizon;  // Previously "end_time" in C++ code

    std::vector<int> mode_linear;
    std::vector<int> refine_order;

    std::vector<double> gammas;
    std::vector<double> k_consts;
    std::vector<double> initial_matrix;
    std::vector<double> initial_b;
    std::vector<double> unsafe_matrix;
    std::vector<double> unsafe_b;

    std::string annot_str;  // LMBTODO: Add to config, move to Simulator???
    std::string beta_str;  // LMBTODO: Add to config, move to Simulator???
    std::string opt_str;  // LMBTODO: Add to config, move to Simulator???
    std::string visualize_filename;
    std::string executable;

    LinearSet initial_set;
    LinearSet unsafe_set;

    template <typename T> void print_vector(std::string str, 
        std::vector<T> vec)
    {
        std::cout << str;
        for (auto v : vec)
        {
            std::cout << v << " ";
        }
        std::cout << std::endl;
    }
};

#endif /* MODEL_H */