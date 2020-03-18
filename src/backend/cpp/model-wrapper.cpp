/**
 * \file   model-wrapper.cpp
 * 
 * \author Matt Potok
 * \author Lucas Brown
 * \date   April 1, 2019
 * 
 * \breif Interface between C++ and Python
 */

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "model.hpp"

namespace bp = boost::python;


BOOST_PYTHON_MODULE(libc2e2)
{
    bp::class_<Model>("Model")
        // Functions
        .def("print_model", &Model::printModel)
        .def("simulate_verify", &Model::simulateVerify)

        // Intenger / Boolean
        .def("set_simulation_bool", &Model::setSimulationBool)
        .def("set_refine_strat", &Model::setRefineStrat)
        .def("set_dimensions", &Model::setDimensions)
        .def("set_num_modes", &Model::setNumModes)
        .def("set_initial_mode_idx", &Model::setInitialModeIdx)
        .def("set_num_initial_eqns", &Model::setNumInitialEqns)
        .def("set_num_unsafe_eqns", &Model::setNumUnsafeEqns)
        .def("set_annot_type", &Model::setAnnotType)  // TODOLMB Do we need this? It's hardcoded to 3. Move to config file?
    
        // Double
        .def("set_abs_err", &Model::setAbsError)
        .def("set_rel_err", &Model::setRelError)
        .def("set_time_step", &Model::setTimeStep)
        .def("set_time_horizon", &Model::setTimeHorizon)

        // Integer Vector
        .def("set_mode_linear", &Model::setModeLinear)
        .def("set_refine_order", &Model::setRefineOrder)

        // Double Vector
        .def("set_gammas", &Model::setGammas)
        .def("set_k_consts", &Model::setKConsts)
        .def("set_initial_matrix", &Model::setInitialMatrix)
        .def("set_initial_b", &Model::setInitialB)
        .def("set_unsafe_matrix", &Model::setUnsafeMatrix)
        .def("set_unsafe_b", &Model::setUnsafeB)

        // String
        .def("set_annot_str", &Model::setAnnotStr)
        .def("set_beta_str", &Model::setBetaStr)
        .def("set_opt_str", &Model::setOptStr)
        .def("set_visualize_filename", &Model::setVisualizeFilename)
        .def("set_executable", &Model::setExecutable)
    ;

    bp::class_<std::vector<int>>("IntegerVector")
        .def(bp::vector_indexing_suite<std::vector<int>>())
    ;
    bp::class_<std::vector<double>>("DoubleVector")
        .def(bp::vector_indexing_suite<std::vector<double>>())
    ;
};