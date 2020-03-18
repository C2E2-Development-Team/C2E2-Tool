/**
 * \file   verifier.hpp
 * \class  Verifier
 * 
 * \author Lucas Brown
 * \date   April 9, 2019
 *
 * \brief  LMBTODO
 */

#include <ppl.hh>
#include <stack>
#include <vector>
#include <Eigen/Dense>

#include "annotation.hpp"
#include "rep-point.hpp"
#include "simulator.hpp"
#include "gstar-reach-tube.hpp"
#include "generalized-star.hpp"

using namespace std;

class Verifier
{
public:
  Verifier(Simulator simu, Annotation annotation);
  ~Verifier();

  int verify();
  void generateCoverStack();
  void rewriteSimulationOutput(ReachTube &simulation_tube);
  void bloatReachTube(ReachTube &simulation_tube,
                      std::vector<double> ref_delta_array, RepPoint curItrRepPoint);
  void checkInvariantsGuards(ReachTube &simulation_tube,
                             ReachTube &guard_set, std::vector<int> &index_in_unsafe_set);
  void refine(RepPoint rep_point, int dim_index);

  void setNumUnsafeEqns(int val);
  void setUnsafeMatrix(std::vector<double> vec);

  std::stack<RepPoint> getCoverStack();

  int ComputeLdf(vector<double> delta, int CT_step, int state);

  void setModeLinear(std::vector<int> vec); // LMBTODO: Band-aid
  void setKConsts(std::vector<double> vec);  // LMBTODO: Band-aid

  void setRefineStrat(int val);
  void setRefineOrder(std::vector<int> vec);

private:
  int num_unsafe_equations; // LMBTODO Replace this with a size function
                            //         Not sure if there's a reason we can't
  std::vector<double> unsafe_matrix;
  std::stack<RepPoint> cover_stack;
  Annotation annotation;

  void *lib;
  typedef std::vector<std::pair<std::NNC_Polyhedron, int>> (*guard_fn)(int, double *, double *);
  typedef bool (*inv_fn)(int, double *, double *);
  guard_fn guards;
  inv_fn invs;

  void *libJaThin;
  typedef std::vector<double> (*jcalc_fn)(vector<double>, int);
  typedef void (*thin_fn)(int, vector<double> &, vector<double> &, vector<double> &, vector<double> &);
  jcalc_fn jcalc;
  thin_fn thinHandle;

  Simulator simulator;
  Checker checker;
  int dimensions;
  LinearSet unsafe_set;
  LinearSet initial_set;
  int initial_mode_idx;
  std::vector<int> mode_linear;
  std::vector<double> k_consts;
  std::string visualize_filename;

  int refine_strat;
  vector<int> refine_order;
};