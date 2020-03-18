/**
 * \file   verifier.cpp
 * \class  Verifier
 * 
 * \author Lucas Brown
 * \date   April 9, 2019
 *
 * \brief  LMBTODO
 */

#include <iostream>
#include <dlfcn.h>
#include <fstream>
#include <ppl.hh>
#include <Python.h>

#include "verifier.hpp"

using namespace std;

Verifier::Verifier(Simulator simu, Annotation annotation)
{
    dimensions = simu.getDimensions();
    initial_mode_idx = simu.getInitialModeIdx();
    initial_set = simu.getInitialSet();
    unsafe_set = simu.getUnsafeSet();
    checker = simu.getChecker();
    mode_linear = simu.getModeLinear();
    k_consts = simu.getKConsts();
    visualize_filename = simu.getVisualizeFilename();
    simulator = simu;

    this->annotation = annotation;

    // LMBTODO Handle failed library opening.
    lib = dlopen("../work-dir/libbloatedsim.so", RTLD_LAZY);
    guards = (guard_fn)dlsym(lib, "hitsGuard");
    invs = (inv_fn)dlsym(lib, "invariantSatisfied");

    libJaThin = dlopen("../work-dir/libJaThin.so", RTLD_LAZY);
    if (!libJaThin)
    {
        cerr << "Cannot open library: " << dlerror() << '\n';
    }
    thinHandle = (thin_fn)dlsym(libJaThin, "thinHandle");
    jcalc = (jcalc_fn)dlsym(libJaThin, "jcalc");
}

Verifier::~Verifier()
{
    dlclose(lib);
    dlclose(libJaThin);
}

int Verifier::verify()
{
    // int LMBTODO = simulator.hybridSimulation();
    // Open bloating library.
    // void *lib = dlopen("../work-dir/libbloatedsim.so", RTLD_LAZY);
    // LMBTODO: Hanlde failed library opening - need to provide an error
    //          and return an "unknown" result.

    // LMBTODO
    // typedef std::vector< std::pair< std::NNC_Polyhedron, int>>
    //     (*guard_fn)(int, double*, double*);
    // typedef bool (*inv_fn)(int, double*, double*);
    // guard_fn guards = (guard_fn) dlsym(lib, "hitsGuard");
    // inv_fn invs = (inv_fn) dlsym(lib, "invariantSatisfied");

    // Clear output file
    std::ofstream invFile;
    invFile.open(visualize_filename);
    invFile.close();

    // LMBTODO Add custom refine order?
    // int refine_method = 0;

    int num_sample_points = 0;
    int refine_threshold = 100; // LMBTODO Move this to config file.
    int num_refinements = 0;
    std::vector<ReachTube> result_tube;
    std::vector<ReachTube> trace_tube;
    generateCoverStack();

    int refine_unsafe_flag = 0; // LMBTODO
    int index_itr = 0;          // LMBTODO    

    int iteration = 0;
    while (!cover_stack.empty())
    {
        std::cout << "-----------------------------\n"
                  << "  Iteration: " << iteration++
                  << "\n-----------------------------" << std::endl;
        RepPoint itr_rep_point = cover_stack.top();
        cover_stack.pop();
        itr_rep_point.print(); // LMBTODO Redo printing. Inlucde refine time.

        // Return unknown result if refine time exceeds threshold.
        if (itr_rep_point.getRefineTime() > refine_threshold)
        {
            std::cout << "\n***** Refine Threshold Exceeded. *****\n"
                      << "\n      Refine Time: " << itr_rep_point.getRefineTime()
                      << "\n      Threshold: " << refine_threshold << std::endl;
            std::ofstream result_stream;
            result_stream.open("../work-dir/Result.dat");
            result_stream << "0" << std::endl;
            result_stream.close();
            return 0;
        }

        Point simulation_point = itr_rep_point.getState();
        int mode_simulated = itr_rep_point.getMode();
        std::vector<double> ref_delta_array = itr_rep_point.getDeltaArray();

        simulator.simulatePoint(simulation_point, mode_simulated);
        ReachTube simulation_tube;
        simulation_tube.setDimensions(dimensions);
        simulation_tube.setMode(mode_simulated);
        simulation_tube.parseInvariantTube("../work-dir/SimuOutput", 0);
        rewriteSimulationOutput(simulation_tube);

        bloatReachTube(simulation_tube, ref_delta_array, itr_rep_point);
        ReachTube guard_set;
        std::vector<int> index_in_unsafe_set;
        checkInvariantsGuards(simulation_tube, guard_set, index_in_unsafe_set);

        int trace_safe_flag = checker.check(simulation_tube, unsafe_set);
        if (trace_safe_flag == 0)
        {
            // Tube unknown, trace to the origin and refine immediately
            std::cout << "Tube Unknown. Refining." << std::endl;
            trace_tube.clear();

            std::vector<double> origin_delta_array;
            if (itr_rep_point.hasParentState())
            {
                origin_delta_array = itr_rep_point.getParentDeltaArray();
            }
            else
            {
                origin_delta_array = itr_rep_point.getDeltaArray();
            }

            if (refine_strat == 0)
            {
                std::cout << "********** Using Default Refinement Strategy **********" << std::endl;
                if (refine_unsafe_flag < 4)
                {
                    std::cout << index_in_unsafe_set.at(
                                     index_itr % index_in_unsafe_set.size())
                              << std::endl; // DEBUG
                    refine(itr_rep_point, index_in_unsafe_set.at(
                                              index_itr % index_in_unsafe_set.size()));
                    num_refinements++;
                    index_itr++;
                    refine_unsafe_flag++;
                }
                else
                {
                    double max = 0;
                    int refine_idx = 0;
                    for (int i = 0; i < dimensions; i++)
                    {
                        if (origin_delta_array[i] > max)
                        {
                            max = origin_delta_array[i];
                            refine_idx = i;
                        }
                    }
                    std::cout <<"refine index: " << refine_idx << std::endl; // DEBUG
                    refine(itr_rep_point, refine_idx);
                    refine_unsafe_flag = 0;
                }
            }
            else
            {
                // std::cout << "Custom Refinement not yet implemented"
                //           << std::endl;
                // LMBTODO;
                cout << "********** Using Custom Refinement Strategy **********" << endl;
                int refine_idx = refine_order[itr_rep_point.getRefineTime() % refine_order.size()] - 1;
                cout <<"refine index: "<< refine_idx << endl;
                refine(itr_rep_point, refine_idx);
            }
        }
        else if (trace_safe_flag == 1)
        {
            std::cout << "Tube Safe! Checking if there are transitions for the"
                      << " next mode." << std::endl;
            trace_tube.push_back(simulation_tube);
            int if_next_set = guard_set.getNextSetStack(cover_stack,
                                                        itr_rep_point);
            std::cout << "if_next_set: " << if_next_set << std::endl;
            if (!if_next_set)
            {
                std::cout << "No more transitions" << std::endl;
                result_tube.reserve(result_tube.size() + trace_tube.size());
                result_tube.insert(result_tube.end(), trace_tube.begin(),
                                   trace_tube.end());
                trace_tube.clear();
            }
        }
        else if (trace_safe_flag == -1)
        {
            std::cout << "Tube Unsafe, Break" << std::endl;
            trace_tube.push_back(simulation_tube);
            result_tube.reserve(result_tube.size() + trace_tube.size());
            result_tube.insert(result_tube.end(), trace_tube.begin(),
                               trace_tube.end());
            trace_tube.clear();
            ReachTube inv_tube;
            for (int i = 0; i < result_tube.size(); i++)
            {
                inv_tube = result_tube.at(i);
                inv_tube.printReachTube(visualize_filename, 1);
            }
            std::cout << "The system is unsafe." << std::endl;
            std::ofstream result_stream;
            result_stream.open("../work-dir/Result.dat");
            result_stream << "-1" << std::endl;
            result_stream.close();
            return -1;
        }
    }

    ReachTube inv_tube;
    std::cout << "Size: " << result_tube.size() << std::endl;
    for (int i = 0; i < result_tube.size(); i++)
    {
        inv_tube = result_tube.at(i);
        inv_tube.printReachTube(visualize_filename, 1);
    }
    std::cout << "System is safe" << std::endl;

    std::ofstream result_stream;
    result_stream.open("../work-dir/Result.dat");
    result_stream << "1" << std::endl;
    result_stream << num_sample_points << std::endl;
    result_stream << num_refinements << std::endl;
    result_stream.close();
    //dlclose(lib);
    return 1;
}

void Verifier::generateCoverStack()
{
    std::cout << "Generating Initial Cover Stack... ";
    RepPoint rep_pnt;

    // Create Delta Array
    std::vector<double> delta_array(dimensions, 0);
    for (int i = 0; i < dimensions; i++)
    {
        delta_array[i] = (initial_set.getMax(i) - initial_set.getMin(i)) / 2;
    }

    // Create Point
    Point pnt(dimensions + 1);
    pnt.setCoordinate(0, 0);
    for (int i = 1; i < dimensions + 1; i++)
    {
        pnt.setCoordinate(i, initial_set.getMin(i - 1) + delta_array[i - 1]);
    }

    rep_pnt.setState(pnt);
    rep_pnt.setDimensions(dimensions);
    rep_pnt.setMode(initial_mode_idx);
    rep_pnt.setDeltaArray(delta_array);

    cover_stack.push(rep_pnt);

    std::cout << "Complete." << std::endl;
}

void Verifier::rewriteSimulationOutput(ReachTube &simulation_tube)
{
    std::cout << "Rewriting Simulation Output... ";
    int sim_tube_size = simulation_tube.getSize();
    if (sim_tube_size)
    {
        // LMBTODO
        // void *lib = dlopen("../work-dir/libbloatedsim.so", RTLD_LAZY);
        // typedef std::vector< std::pair< std::NNC_Polyhedron, int>>
        // (*guard_fn)(int, double*, double*);
        // guard_fn guards = (guard_fn) dlsym(lib, "hitsGuard");

        double *lower_point = NULL;
        double *upper_point = NULL;
        std::vector<std::pair<std::NNC_Polyhedron, int>> guards_hit;
        for (int i = 0; i < sim_tube_size; i++)
        {
            std::vector<double> lower =
                simulation_tube.getLowerBound(i).getCoordinates();
            std::vector<double> upper =
                simulation_tube.getUpperBound(i).getCoordinates();
            lower_point = &lower[0]; // C++11 and later
            upper_point = &upper[0];
            guards_hit = guards(simulation_tube.getMode(),
                                lower_point, upper_point);

            if (!guards_hit.empty())
            {
                simulation_tube.clear(i + 15); // LMBTODO Config?
                break;
            }
        }
        std::ofstream sim_output_rewrite;
        sim_output_rewrite.open("../work-dir/SimuOutput");
        for (int i = 0; i < simulation_tube.getSize(); i++)
        {
            // Print Box Coordinates
            Point lower_bound = simulation_tube.getLowerBound(i);
            Point upper_bound = simulation_tube.getUpperBound(i);
            for (int j = 0; j < lower_bound.getDimensions(); j++)
            {
                sim_output_rewrite << lower_bound.getCoordinate(j) << " ";
            }
            sim_output_rewrite << std::endl;
            for (int j = 0; j < upper_bound.getDimensions(); j++)
            {
                sim_output_rewrite << upper_bound.getCoordinate(j) << " ";
            }
            sim_output_rewrite << std::endl;
        }
        sim_output_rewrite.close();
        //dlclose(lib);
    }
    std::cout << "Complete." << std::endl;
}

void Verifier::bloatReachTube(ReachTube &simulation_tube,
                              std::vector<double> ref_delta_array, RepPoint curItrRepPoint)
{
    std::cout << "Bloating Reach Tube... " << std::endl;
    int mode_simulated = simulation_tube.getMode();
    if (mode_linear[mode_simulated - 1] == 0)
    {
        // Py_Initialize();
        // char input_buff[128];
        // char input_buff2[32];
        // char input_buff3[32];
        // char input_buff4[32];
        // FILE *fid = NULL;
        // const char *file_name = "../work-dir/ComputeLDF.py";
        // double CT_step = k_consts[mode_simulated-1];
        // int modeforpython = mode_simulated;
        // strcpy (input_buff, "delta = [");
        // for (int i = 0; i< dimensions; i++)
        // {
        //     char temp [8];
        //     sprintf(temp,"%f", ref_delta_array[i]);
        //     strcat(input_buff,temp);
        //     if (i<dimensions-1)
        //         strcat(input_buff,",");
        // }
        // strcat(input_buff,"]");
        // PyRun_SimpleString(input_buff);
        // sprintf(input_buff2,"CT_step = int(%f)", CT_step);
        // PyRun_SimpleString(input_buff2);
        // sprintf(input_buff3,"state = '%d'", modeforpython);
        // PyRun_SimpleString(input_buff3);
        // sprintf(input_buff4,"Is_linear = int(%d)",
        //         mode_linear[mode_simulated-1]);
        // PyRun_SimpleString(input_buff4);
        // fid = fopen(file_name, "r");
        // PyRun_SimpleFile(fid, file_name);
        // fclose(fid);
        // std::cout << input_buff << std::endl;
        // std::cout << input_buff2 << std::endl;
        // std::cout << input_buff3 << std::endl;
        // std::cout << input_buff4 << std::endl;
        int CT_step = (int)k_consts[mode_simulated - 1];
        int modeforpython = mode_simulated;
        ComputeLdf(ref_delta_array, CT_step, modeforpython);

        simulation_tube = ReachTube();
        simulation_tube.setDimensions(dimensions);
        simulation_tube.setMode(mode_simulated);
        simulation_tube.parseInvariantTube("../work-dir/reachtube.dat", 1);
    }
    else
    {
        // ReachTube bloated_tube;
        // bloated_tube = simulation_tube.bloatReachTube(ref_delta_array,
        //                                               annotation);
        // bloated_tube.printReachTube("../work-dir/reachtube.dat", 0);
        simulation_tube = ReachTube();
        gstar star = gstar(curItrRepPoint);
        star.normalize();
        gReachTube tube(star);
        tube.setSimulator(simulator);
        // tube.setSimulatorStr(getExecutable());
        tube.computeSimulationTube();
        tube.bloatTube();
        tube.tubeToReachTube(simulation_tube);
    }
    std::cout << "Bloating Complete." << std::endl;
}

void Verifier::checkInvariantsGuards(ReachTube &simulation_tube,
                                     ReachTube &guard_set, std::vector<int> &index_in_unsafe_set)
{
    std::cout << "Checking Invariants and Guards... " << std::endl;
    //void *lib = dlopen("../work-dir/libbloatedsim.so", RTLD_LAZY);
    // LMBTODO
    // typedef std::vector< std::pair< std::NNC_Polyhedron, int>>
    //     (*guard_fn)(int, double*, double*);
    // typedef bool (*inv_fn)(int, double*, double*);
    // guard_fn guards = (guard_fn) dlsym(lib, "hitsGuard");
    // inv_fn invs = (inv_fn) dlsym(lib, "invariantSatisfied");

    // LMBTODO: This seems a bit redundant - we seem to be creating a reach
    // tube from the bloated simulation_tube, printing it out, then
    // creating a new reach tube by reading in what we just printed out.
    // Need to look into if this is the best way to handle this.

    guard_set.setDimensions(dimensions);

    // Initialize refinement dimensions
    std::vector<int> refine_flag(dimensions, 0);
    for (int i = 0; i < num_unsafe_equations * dimensions; i++)
    {
        if (unsafe_matrix[i] != 0.0)
        {
            refine_flag[i % dimensions] = 1;
        }
    }
    //std::vector<int> index_in_unsafe_set;
    for (int i = 0; i < dimensions; i++)
    {
        if (refine_flag[i])
        {
            index_in_unsafe_set.push_back(i);
        }
    }

    double *lower_point = NULL;
    double *upper_point = NULL;
    std::vector<std::pair<std::NNC_Polyhedron, int>> guards_hit;
    bool hit_guard = false;
    for (int i = 0; i < simulation_tube.getSize(); i++)
    {
        // Print Box Corners
        std::vector<double> lower =
            simulation_tube.getLowerBound(i).getCoordinates();
        std::vector<double> upper =
            simulation_tube.getUpperBound(i).getCoordinates();
        lower_point = &lower[0]; // Valid for C++11 and on
        upper_point = &upper[0];
        bool inv_true = invs(simulation_tube.getMode(), lower_point, upper_point);
        if (!inv_true)
        {
            // Invariant not satisfied.
            simulation_tube.clear(i);
            break;
        }
        guards_hit = guards(simulation_tube.getMode(), lower_point, upper_point);
        if (!guards_hit.empty())
        {
            guard_set.addGuards(guards_hit);
            hit_guard = true;
        } // LMBTODO: Add documentation for the else if logic
        else if (hit_guard == true)
        {
            simulation_tube.clear(i);
            break;
        }
    }
    //dlclose(lib);
    std::cout << "Invarnts and Guard Checking Complete." << std::endl;
}

void Verifier::refine(RepPoint rep_point, int dim_index)
{
    int dim = rep_point.getDimensions();
    if (dim_index >= dim || dim_index < 0)
    {
        std::cout << "Invalid dim_index" << std::endl; // Add error handling
        return;
    }

    double max[dim];
    double min[dim];
    std::vector<double> delta_array;
    Point temp_point;
    int init_mode;
    int refine_time = rep_point.getRefineTime();

    if (rep_point.hasParentState())
    {
        std::cout << "===========Trace back to init mode, refine, Push two rep point to Stack!==============="
                  << std::endl;
        delta_array = rep_point.getParentDeltaArray();
        temp_point = rep_point.getParentState();
        init_mode = rep_point.getParentMode();
    }
    else
    {
        std::cout << "======refine happen in init mode, using it's own state, Push two rep point to Stack!===="
                  << std::endl;
        delta_array = rep_point.getDeltaArray();
        temp_point = rep_point.getState();
        init_mode = rep_point.getMode();
    }

    for (int i = 0; i < dim; i++)
    {
        max[i] = temp_point.getCoordinate(i + 1) + delta_array[i];
        min[i] = temp_point.getCoordinate(i + 1) - delta_array[i];
    }

    delta_array[dim_index] = delta_array[dim_index] / 2;

    Point p1(dim + 1);
    RepPoint rep1;
    p1.setCoordinate(0, 0);
    for (int i = 1; i <= dim; i++)
    {
        p1.setCoordinate(i, min[i - 1] + delta_array[i - 1]);
    }
    rep1.setState(p1);
    rep1.setDimensions(dim);
    rep1.setMode(init_mode);

    // LMBTODO This was marked as a temporary fix in the original code
    std::vector<double> temp_delta_array(delta_array.begin(),
                                         delta_array.end());

    rep1.setDeltaArray(temp_delta_array);
    rep1.setRefineTime(refine_time + 1);
    cover_stack.push(rep1);

    Point p2(dim + 1);
    RepPoint rep2;
    p2.setCoordinate(0, 0);
    for (int i = 1; i <= dim; i++)
    {
        p2.setCoordinate(i, max[i - 1] - delta_array[i - 1]);
    }
    rep2.setState(p2);
    rep2.setDimensions(dim);
    rep2.setMode(init_mode);
    rep2.setDeltaArray(temp_delta_array);
    rep2.setRefineTime(refine_time + 1);
    cover_stack.push(rep2);

    return;
}

void Verifier::setNumUnsafeEqns(int val)
{
    num_unsafe_equations = val;
}
void Verifier::setUnsafeMatrix(std::vector<double> vec)
{
    unsafe_matrix = vec;
}

void Verifier::setModeLinear(std::vector<int> vec)
{
    mode_linear = vec;
}

void Verifier::setRefineStrat(int val)
{
    refine_strat = val;
}

void Verifier::setRefineOrder(std::vector<int> order)
{
    refine_order = order;
}

void Verifier::setKConsts(std::vector<double> vec)
{
    k_consts = vec;
}

std::stack<RepPoint> Verifier::getCoverStack()
{
    return cover_stack;
}

int Verifier::ComputeLdf(vector<double> delta, int CT_step, int state)
{
    vector<double> notbloating;
    vector<double> thin_indicate;
    vector<double> bloating;
    vector<double> not_thin;
    cout << "Begin compute LDF" << endl;
    thinHandle(state, notbloating, thin_indicate, bloating, not_thin);
    vector<double> thin_Jac_ind;
    for (int not_thin_idx = 0; not_thin_idx < not_thin.size(); not_thin_idx++)
    {
        double new_non_ele = not_thin[not_thin_idx];
        for (int not_bloat_idx = 0; not_bloat_idx < notbloating.size(); not_bloat_idx++)
        {
            if (notbloating[not_bloat_idx] < not_thin[not_thin_idx])
            {
                new_non_ele -= 1;
            }
        }
        thin_Jac_ind.push_back(new_non_ele);
    }
    for (int i = 0; i < thin_indicate.size(); i++)
    {
        thin_indicate[i]++;
    }

    for (int i = 0; i < bloating.size(); i++)
    {
        bloating[i]++;
    }

    for (int i = 0; i < not_thin.size(); i++)
    {
        not_thin[i]++;
    }

    for (int i = 0; i < notbloating.size(); i++)
    {
        notbloating[i]++;
    }
    cout << "read simulation output" << endl;
    ifstream fin("../work-dir/SimuOutput");
    string line;
    int curline = 0;
    Eigen::MatrixXd Simulation_data;
    while (getline(fin, line))
    {
        Eigen::RowVectorXd dataLine;
        stringstream lineStream(line);
        double value = 0;
        while (lineStream >> value)
        {
            dataLine.conservativeResize(dataLine.size() + 1);
            dataLine(dataLine.size() - 1) = value;
        }

        Simulation_data.conservativeResize(Simulation_data.rows() + 1, dataLine.size());
        Simulation_data.row(Simulation_data.rows() - 1) = dataLine;
    }
    fin.close();

    cout << "start bloating reach tube" << endl;
    int numofvar = Simulation_data.cols();
    CT_step = min(CT_step, (int)Simulation_data.rows());
    Eigen::MatrixXd blowting = Eigen::MatrixXd::Zero((int)Simulation_data.rows(), (int)not_thin.size());
    vector<double> deltacopy = delta;
    vector<double> tmp_list;
    for (int tmp_list_idx = 0; tmp_list_idx < not_thin.size(); tmp_list_idx++)
    {
        tmp_list.push_back(not_thin[tmp_list_idx]);
    }

    for (int tmp = 0; tmp < tmp_list.size(); tmp++)
    {
        blowting(0, tmp) = delta[tmp];
    }
    numofvar -= (1 + notbloating.size());
    vector<vector<double>> Reach_tube;
    if (numofvar > 0)
    {
        for (int i = 0; i < Simulation_data.rows(); i += CT_step)
        {
            vector<double> mean_array;
            Eigen::MatrixXd Simulation_data_truncate = Simulation_data.block(i, 0, min(CT_step, (int)Simulation_data.rows() - i), Simulation_data.cols());
            for (int j = 0; j < bloating.size(); j++)
            {
                double mean = Simulation_data_truncate.col(bloating[j]).mean();
                mean_array.push_back(mean);
            }
            if (mean_array.size() != numofvar)
            {
                throw runtime_error("Error in calculating the mean matrix for coordinated transformation");
                ;
            }
            vector<double> Jacobian_CT_vec = jcalc(mean_array, state);
            Eigen::MatrixXd Jacobian_CT1 = Eigen::Map<Eigen::MatrixXd>(Jacobian_CT_vec.data(), numofvar, numofvar);
            Jacobian_CT1.transposeInPlace();
            Eigen::MatrixXd Jacobian_CT2 = Eigen::MatrixXd::Zero((int)Jacobian_CT1.rows(), thin_Jac_ind.size());
            for (int j = 0; j < thin_Jac_ind.size(); j++)
            {
                Jacobian_CT2.col(j) = Jacobian_CT1.col(thin_Jac_ind[j]);
            }

            Eigen::MatrixXd Jacobian_CT = Eigen::MatrixXd::Zero(thin_Jac_ind.size(), (int)Jacobian_CT2.cols());
            for (int j = 0; j < thin_Jac_ind.size(); j++)
            {
                Jacobian_CT.row(j) = Jacobian_CT2.row(thin_Jac_ind[j]);
            }
            Eigen::EigenSolver<Eigen::MatrixXd> Jacobian_CT_solver(Jacobian_CT);
            Eigen::MatrixXcd CT_matrix = Jacobian_CT_solver.eigenvectors();
            Eigen::MatrixXcd eigvalues = Jacobian_CT_solver.eigenvalues();
            double Local_Lipschitz = Jacobian_CT.norm();
            for (int j = i; j < min(i + CT_step, (int)Simulation_data.rows()); j += 2)
            {
                if (j + 1 > Simulation_data.rows())
                {
                    throw runtime_error("Error accessing the last rectangle");
                    ;
                }
                double Delta_time = Simulation_data(j + 1, 0) - Simulation_data(j, 0);
                Eigen::MatrixXd Simulation_box = Eigen::MatrixXd::Zero(2, bloating.size());
                for (int k = 0; k < bloating.size(); k++)
                {
                    Simulation_box.col(k) = Simulation_data.block(j, bloating[k], 2, 1);
                }
                Eigen::MatrixXd meanVec = (Simulation_box.row(0) + Simulation_box.row(1)) / 2;
                vector<double> Current_Jacobian_vec = jcalc(vector<double>(meanVec.data(), meanVec.data() + meanVec.size()), state);
                Eigen::MatrixXd Current_Jacobian1 = Eigen::Map<Eigen::MatrixXd>(Current_Jacobian_vec.data(), numofvar, numofvar);
                Current_Jacobian1.transposeInPlace();

                Eigen::MatrixXd Current_Jacobian2 = Eigen::MatrixXd::Zero((int)Current_Jacobian1.rows(), thin_Jac_ind.size());
                for (int k = 0; k < thin_Jac_ind.size(); k++)
                {
                    Current_Jacobian2.col(k) = Current_Jacobian1.col(thin_Jac_ind[k]);
                }

                Eigen::MatrixXd Current_Jacobian = Eigen::MatrixXd::Zero(thin_Jac_ind.size(), (int)Current_Jacobian2.cols());
                for (int k = 0; k < thin_Jac_ind.size(); k++)
                {
                    Current_Jacobian.row(k) = Current_Jacobian2.row(thin_Jac_ind[k]);
                }
                Eigen::MatrixXcd Current_Jordan = CT_matrix.inverse() * Current_Jacobian * CT_matrix;
                Eigen::MatrixXcd Current_Jordan_transpose = Current_Jordan.transpose();
                Eigen::MatrixXcd Current_Jordan_Hermitian = Current_Jordan + Current_Jordan_transpose;
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> Current_Jordan_solver(Current_Jordan_Hermitian);
                Eigen::MatrixXd Current_Jordan_eivalues = Current_Jordan_solver.eigenvalues();
                double Current_lambda = Current_Jordan_eivalues.maxCoeff() / 2;
                double Max_blowting = blowting.row(j).maxCoeff();

                Eigen::MatrixXd Simulation_box_max = Simulation_box.colwise().maxCoeff();
                Simulation_box_max = Simulation_box_max.array() + Max_blowting / 2;
                vector<double> Max_Jacobian_vec = jcalc(vector<double>(Simulation_box_max.data(), Simulation_box_max.data() + Simulation_box_max.size()), state);
                Eigen::MatrixXd Max_Jacobian = Eigen::Map<Eigen::MatrixXd>(Max_Jacobian_vec.data(), numofvar, numofvar);
                Max_Jacobian.transposeInPlace();

                Eigen::MatrixXd Simulation_box_min = Simulation_box.colwise().minCoeff();
                Simulation_box_min = Simulation_box_min.array() - Max_blowting / 2;
                vector<double> Min_Jacobian_vec = jcalc(vector<double>(Simulation_box_min.data(), Simulation_box_min.data() + Simulation_box_min.size()), state);
                Eigen::MatrixXd Min_Jacobian = Eigen::Map<Eigen::MatrixXd>(Min_Jacobian_vec.data(), numofvar, numofvar);
                Min_Jacobian.transposeInPlace();

                Eigen::MatrixXd Disturb_matrix = Max_Jacobian - Min_Jacobian;
                double Disturbance = Disturb_matrix.operatorNorm() * exp(Local_Lipschitz * Delta_time);
                Current_lambda = Current_lambda + Disturbance;
                if (j == 0)
                {
                    blowting.row(1) = blowting.row(0) * exp(Current_lambda * Delta_time);
                }
                else
                {
                    blowting.row(j) = blowting.row(j - 1);
                    blowting.row(j + 1) = blowting.row(j) * exp(Current_lambda * Delta_time);
                }
            }
            for (int cnt = 0; cnt < min(CT_step, (int)Simulation_data.rows() - i); cnt++)
            {
                double cond = CT_matrix.operatorNorm() * CT_matrix.inverse().operatorNorm();
                blowting.row(i + cnt) = blowting.row(i + cnt) * cond;
            }
        }
        for (int i = 0; i < Simulation_data.rows(); i += 2)
        {
            vector<double> reachtube_i(Simulation_data.cols(), 0.0);
            vector<double> reachtube_i1(Simulation_data.cols(), 0.0);
            for (int j = 0; j < not_thin.size(); j++)
            {
                reachtube_i[not_thin[j]] = min(Simulation_data(i, not_thin[j]), Simulation_data(i + 1, not_thin[j])) - blowting(i, j);
                reachtube_i1[not_thin[j]] = max(Simulation_data(i, not_thin[j]), Simulation_data(i + 1, not_thin[j])) + blowting(i, j);
            }

            for (int j = 0; j < notbloating.size(); j++)
            {
                reachtube_i[notbloating[j]] = Simulation_data(i, notbloating[j]) - deltacopy[notbloating[j] - 1];
                reachtube_i1[notbloating[j]] = Simulation_data(i + 1, notbloating[j]) + deltacopy[notbloating[j] - 1];
            }

            for (int j = 0; j < thin_indicate.size(); j++)
            {
                reachtube_i[thin_indicate[j]] = Simulation_data(i, thin_indicate[j]);
                reachtube_i1[thin_indicate[j]] = Simulation_data(i + 1, thin_indicate[j]);
            }

            Reach_tube.push_back(reachtube_i);
            Reach_tube.push_back(reachtube_i1);
        }
    }
    for (int i = 0; i < Simulation_data.rows(); i += 2)
    {
        Reach_tube[i][0] = Simulation_data(i, 0);
        Reach_tube[i + 1][0] = Simulation_data(i + 1, 0);
    }
    ofstream fout("../work-dir/reachtube.dat");
    fout.precision(numeric_limits<double>::digits10 + 1);
    fout << state << endl;
    for (int i = 0; i < Reach_tube.size(); i++)
    {
        for (int j = 0; j < Reach_tube[0].size(); j++)
        {
            fout << Reach_tube[i][j] << " ";
        }
        fout << endl;
    }
    fout.close();
    return 0;
}
