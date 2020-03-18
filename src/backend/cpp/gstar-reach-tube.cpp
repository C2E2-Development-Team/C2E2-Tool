#include "gstar-reach-tube.hpp"

using namespace std;
using namespace Eigen;

/*
gReachTube()
    Description:
        Default constructor, create an empty reach tube.
*/
gReachTube::gReachTube()
{
    dimension = 0;
    mode = -1;
    tube = vector<gstar>();
    initialSet = gstar();
    bloated = 0;
    epsilon = 0;
}

/*
gReachTube(gstar star)
    Input:
        gstar star: a generalized star
    Description:
        create a reach tube from an initial set represented by gstar star
*/
gReachTube::gReachTube(gstar star)
{
    dimension = star.getDimension();
    mode = star.getMode();
    tube = vector<gstar>();
    initialSet = star;
    bloated = 0;
    epsilon = 0;
}

/*
int getDimension()
    Description:
        return the dimension of the reach tube
*/
int gReachTube::getDimension()
{
    return dimension;
}

/*
int getMode()
    Description:
        return the mode of the reach tube
*/
int gReachTube::getMode()
{
    return mode;
}

/*
vector<gstar> getTube()
    Description:
        return the vector of gstar containing the reach tube
*/
vector<gstar> gReachTube::getTube()
{
    return tube;
}

/*
gstar getStar(int index)
    Input:
        int index: the index of the generalized star to be accessed
    Description:
        return the generalized star in the reach tube given by index
*/
gstar gReachTube::getStar(int index)
{
    if (index >= tube.size() || index < 0)
    {
        throw std::invalid_argument("index out of range");
        return gstar();
    }
    return tube[index];
}

/*
double getEpsilon()
    Description:
        return the epsilon value of the reach tube
*/
double gReachTube::getEpsilon()
{
    return epsilon;
}

/*
double isBloated()
    Description:
        check if the the reach tube is already bloated
*/
int gReachTube::isBloated()
{
    return bloated;
}

/*
void setDimension(int newDimension()
    Input:
        int newDimension: the new dimension of the reach tube
    Description:
        set the dimension of the reach tube
*/
void gReachTube::setDimension(int newDimension)
{
    dimension = newDimension;
}

/*
void setMode(int newMode)
    Input:
        int newMode: the new mode of the reach tube
    Description:
        set the mode of the reach tube
*/
void gReachTube::setMode(int newMode)
{
    mode = newMode;
}

/*
void addToTube(gstar star)
    Input:
        gstar star: a generalized to be added to the tube
    Description:
        push a generazlied star star to the back of the reach tube
*/
void gReachTube::addToTube(gstar star)
{
    tube.push_back(star);
}

/*
void setSimulator(Simulator *simu)
    Input:
        Simulator *simu: a pointer to a c2e2 simulator
    Description:
        set the simulator for the reach tube
*/
void gReachTube::setSimulator(Simulator simu)
{
    // if (simu == NULL)
    // {
    //     throw std::invalid_argument("invalid simulator");
    //     return;
    // }
    modelSim = simu;
}

/*
double computeEps(double *ptn1, double *ptn2)
    Input:
        double *ptn1: an array of double represent a point 
            the size of the array is equal to dimension
        double *ptn2: an array of double represent another point
            the size of the array is equal to dimension
    Output:
        double: the distance between the two points
    Description:
        According to theorem 4 of the paper, a bloating algorithm is proposed.
        The algorithm require an epsilon value that define simulation trace between different time steps.
        Here, since we use rectangular to handle the value between time steps, the diagonal of the rectangular
        polytope, which is the distance between two simulation points can be used as the epsilon value.
*/
double gReachTube::computeEps(vector<double> ptn1, vector<double> ptn2)
{
    double res = 0;
    // Map<VectorXd> vec1(ptn1, dimension + 1);
    // Map<VectorXd> vec2(ptn2, dimension + 1);
    Map<VectorXd> vec1(&ptn1[1], dimension);
    Map<VectorXd> vec2(&ptn2[1], dimension);
    res = (vec1 - vec2).norm() / 2;
    return res;
}

/*
void computeSimulationTube()
    Description:
        This function compute a discrete simulation tube from simulation.
        The function will first perform dimension+1 simulations from the initial set.
        Then compute generalized star from each time step of the simulations.
        The result will be stored in tube class member variable
*/
void gReachTube::computeSimulationTube()
{
    vector<vector<Point>> simuResult;
    // points for simulation from initial set
    vector<Point> simulationPoint = initialSet.getSimPoint();
    cout << "doing simulation" << endl;
    for (int i = 0; i < simulationPoint.size(); i++)
    {
        modelSim.simulatePoint(simulationPoint[i], mode);
        // parse simulation result into vector
        vector<Point> test = parseSimRes("../work-dir/SimuOutput", 0);
        // push result in simuResult
        simuResult.push_back(test);
    }
    cout << "converting simulation result to gstar" << endl;
    //iterate through each time step in simulation
    for (int i = 0; i < simuResult[0].size(); i++)
    {
        gstar star;
        vector<vector<double>> temp_basis;
        Point temp_point = simuResult[0][i];

        // compute what happened between each time step
        if (i > 0)
        {
            Point pre_point = simuResult[0][i - 1];
            double temp_eps = computeEps(pre_point.getCoordinates(), temp_point.getCoordinates());
            if (temp_eps > epsilon)
            {
                epsilon = temp_eps;
            }
        }

        // create generalized star from each time step information
        vector<double> temp_center = temp_point.getCoordinates();
        star.setDimension(dimension);
        star.setMode(mode);
        star.setTime(temp_center[0]);
        star.setAlpha(initialSet.getAlpha());
        star.setCenter(vector<double>(&temp_center[1], &temp_center[1 + dimension]));
        // calculate generalized star and stor result in tube
        for (int j = 0; j < dimension; j++)
        {
            temp_point = simuResult[j + 1][i];
            vector<double> temp_coord = temp_point.getCoordinates();
            vector<double> temp_basis_entry = computeBasis(temp_coord, temp_center);
            temp_basis.push_back(temp_basis_entry);
        }
        star.setBasis(temp_basis);
        tube.push_back(star);
    }

    // // remove created local variable
    // for (int i = 0; i < simuResult[0].size(); i++)
    // {
    //     for (int j = 0; j < simulationPoint.size(); j++)
    //     {
    //         delete simuResult[j][i];
    //     }
    // }
}

/*
void doRandomSimulation()
    Description:
        Perform dimension*dimension simulation.
        Each simulation started from a point randomly chosen from the initial set.
        The simulation result is stored in file SimuOutput0-SimuOutput{dimension*dimension-1} in work-dir.
*/
void gReachTube::doRandomSimulation()
{
    vector<double> alpha = initialSet.getAlpha();
    vector<vector<double>> basis = initialSet.getBasis();
    vector<double> center = initialSet.getCenter();
    for (int i = 0; i < dimension * dimension; i++)
    {
        vector<double> temp;
        temp.push_back(initialSet.getTime());
        for (int j = 0; j < dimension; j++)
        {
            temp.push_back(center[j]);
            for (int k = 0; k < dimension; k++)
            {
                srand(time(NULL));
                temp[j + 1] += (basis[k][j] * (2 * alpha[j] * (double)rand() / (double)RAND_MAX - alpha[j]));
            }
        }
        Point simu_pnt = Point(dimension + 1, temp);
        modelSim.simulatePoint(simu_pnt, mode, i);
        // delete[] temp;
    }
}

/*
vector<double> computeBasis(double *ptn, double *center)
    Input:
        double *ptn: an array of double which contains the coordinate of a point in simulation
        double *center: an array of double which corresponding to the simulation of the center
    Output:
        vector<double>: a vector of double which contains the a basis vector for a generalized star for a variable
    Description:
        This function take in the coordinate of the center of the generalized and also another simulation point.
        It calculate the basis vector corresponding to that simulation point
*/
vector<double> gReachTube::computeBasis(vector<double> ptn, vector<double> center)
{
    // vector<double> res_comp;
    // for (int i = 0; i < dimension; i++)
    // {
    //     res_comp.push_back(ptn[i + 1] - center[i + 1]);
    // }
    Map<VectorXd> vec_ptn(&ptn[1], dimension);
    Map<VectorXd> cen_ptn(&center[1], dimension);
    vec_ptn = vec_ptn - cen_ptn;

    vector<double> res(&ptn[1], &ptn[dimension + 1]);
    return res;
}

/*
vector<Point *> parseSimRes(const char *filename, int hasMode)
    Input:
        const char *filename: a char array storing the file name of the simulation result
        int hasMode: a flag that indicate if the smulation have mode
    Output:
        vector<Point *>: a vector of pointer to Point that represent a simulation trace
    Description:
        This function parse the simulation result generated by the simulator.
        It read in the result file and convert into a vector of point object.
*/
vector<Point> gReachTube::parseSimRes(const char *filename, int hasMode)
{
    vector<Point> simu = vector<Point>();
    if (filename == NULL)
    {
        throw invalid_argument("invalid file name");
        return simu;
    }
    FILE *fd = fopen(filename, "r");
    if (fd == NULL)
    {
        throw invalid_argument("simulation file not found");
        return simu;
    }

    double bufReader = 0;
    double lastTime = -1;
    while (fscanf(fd, "%lf", &bufReader) != EOF)
    {
        if (lastTime != bufReader)
        {
            Point tracePoint(dimension + 1);
            lastTime = bufReader;
            tracePoint.setCoordinate(0, bufReader);
            for (int i = 0; i < dimension; i++)
            {
                fscanf(fd, "%lf", &bufReader);
                tracePoint.setCoordinate(i + 1, bufReader);
            }
            simu.push_back(tracePoint);
        }
        else
        {
            lastTime = bufReader;
            for (int i = 0; i < dimension; i++)
            {
                fscanf(fd, "%lf", &bufReader);
            }
        }
    }
    fclose(fd);
    return simu;
}

/*
void printTubeFile(string fn, int flag)
    Input:
        string fn: file name of the output file
        int flag: flag to determine if additional information should be printed
    Description:
        This function will print the current tube to file. 
        It will print out the exact reach tube including the dimension, alpha value, basis vector and center for each generalized star.
        The first value in the file is the dimension of the model.
        Then the dimension values are the alpha values corresponding to the reach tube.
        Then there will be a line specify the mode of the reach tube in format "% mode #".
        Then there will be data sections. The first column of the data section is the time information. The 2:end columns are the center and basis vector of star
        The 1:(dimension+1):end line of the data section corresponding to the center of each of the generalized star.
        Other lines corresponding to the basis vector of the generalized star. Each row of the data corresponding to a basis vector.
*/
void gReachTube::printTubeFile(string fn, int flag)
{
    if (tube.empty())
    {
        cout << "simulation tube is empty, return" << endl;
    }
    cout << "start print result to file" << endl;
    ofstream opfile;
    if (flag == 0)
    {
        opfile.open(fn);
        opfile << dimension << endl;
        for (int i = 0; i < initialSet.getAlpha().size(); i++)
        {
            opfile << initialSet.getAlpha()[i] << ' ';
        }
        opfile << endl;
    }
    else
    {
        opfile.open(fn, ios::app);
    }
    opfile << " % mode " << mode << endl;
    for (int i = 0; i < tube.size(); i++)
    {
        // print center to file
        opfile << tube[i].getTime();
        for (int j = 0; j < dimension; j++)
        {
            opfile << " " << tube[i].getCenter()[j];
        }
        opfile << endl;

        // print each basis to file
        for (int j = 0; j < dimension; j++)
        {
            opfile << tube[i].getTime();
            for (int k = 0; k < dimension; k++)
            {
                opfile << " " << tube[i].getBasis()[j][k];
            }
            opfile << endl;
        }
    }
}

/*
void printTubeFile(string fn, int flag)
    Input:
        string fn: file name of the output file
        int flag: flag to determine if additional information should be printed
    Description:
        This function will print the current tube to file. 
        It will print out a bloated reach tube.
        The idea is to use a tightly fit rectangular box to over approximate the polytope. 
        The first value in the file is the dimension of the model.
        Then there will be a line specify the mode of the reach tube in format "% mode #".
        Then there will be data sections. The first column of the data section is the time information. 
        The 2:end columns are the maximum and minium value on each dimension for the generalized star.
        The 2:2:end columns are the minimum value on each dimension for the generalized star.
        The 3:2:end columns are the maximum value on each dimension for the generalized star.
*/
void gReachTube::printBloatTubeFile(string fn, int flag)
{
    if (tube.empty())
    {
        cout << "simulation tube is empty, return" << endl;
    }
    cout << "start print result to file" << endl;
    ofstream opfile;
    if (flag == 0)
    {
        opfile.open(fn);
        opfile << dimension << endl;
        for (int i = 0; i < initialSet.getAlpha().size(); i++)
        {
            opfile << initialSet.getAlpha()[i] << ' ';
        }
        opfile << endl;
    }
    else
    {
        opfile.open(fn, ios::app);
    }

    // print the mode information of the tube
    opfile << " % mode " << mode << endl;

    // start printing the data section
    for (int i = 0; i < tube.size(); i++)
    {
        gstar star = tube[i];
        // first column of the data section is the time infromation
        opfile << star.getTime();

        // calculate maximum and minimum value for each dimension
        for (int j = 0; j < dimension; j++)
        {
            double upper_val = star.getCenter(j);
            double lower_val = upper_val;
            for (int k = 0; k < dimension; k++)
            {
                double basis_val = star.getBasis(k, j);
                if (basis_val > 0)
                {
                    upper_val = upper_val + star.getAlpha(k) * basis_val;
                    lower_val = lower_val - star.getAlpha(k) * basis_val;
                }
                else
                {
                    upper_val = upper_val - star.getAlpha(k) * basis_val;
                    lower_val = lower_val + star.getAlpha(k) * basis_val;
                }
            }
            opfile << " " << lower_val;
            opfile << " " << upper_val;
        }
        opfile << endl;
    }
}

/*
void tubeToReachTube(ReachTube *rt)
    Input: 
        ReachTube *rt: a pointer to an existing reachtube object
    Description:
        This function will convert a generalized star reach tube to a c2e2 reachtube using the same
            bloating method used in function printBloatTubeFile.
        The upper bound for each time step is the maximum value in two adjacent generalized star.
        The lower bound for each time step is the minimum value in two adjacent genrealized star.
*/
void gReachTube::tubeToReachTube(ReachTube &rt)
{
    if (tube.empty())
    {
        cout << "simulation tube is empty, return" << endl;
    }

    // if (rt == NULL)
    // {
    //     throw std::invalid_argument("reach tube doesn't exist");
    //     return;
    // }

    rt.setDimensions(dimension);
    rt.setMode(mode);

    for (int i = 0; i < tube.size() - 1; i++)
    {
        gstar star1 = tube[i];
        gstar star2 = tube[i + 1];
        vector<double> upper_arr;
        vector<double> lower_arr;
        upper_arr.push_back(star2.getTime());
        lower_arr.push_back(star1.getTime());

        for (int j = 0; j < dimension; j++)
        {
            double upper_val1 = star1.getCenter(j);
            double lower_val1 = upper_val1;

            double upper_val2 = star2.getCenter(j);
            double lower_val2 = upper_val2;

            for (int k = 0; k < dimension; k++)
            {
                double basis_val1 = star1.getBasis(k, j);
                double basis_val2 = star2.getBasis(k, j);
                if (isnan(basis_val1) || isnan(basis_val2))
                {
                    cout << "stop here" << endl;
                }
                if (basis_val1 > 0)
                {
                    upper_val1 = upper_val1 + star1.getAlpha(k) * basis_val1;
                    lower_val1 = lower_val1 - star1.getAlpha(k) * basis_val1;
                }
                else
                {
                    upper_val1 = upper_val1 - star1.getAlpha(k) * basis_val1;
                    lower_val1 = lower_val1 + star1.getAlpha(k) * basis_val1;
                }
                if (basis_val2 > 0)
                {
                    upper_val2 = upper_val2 + star2.getAlpha(k) * basis_val2;
                    lower_val2 = lower_val2 - star2.getAlpha(k) * basis_val2;
                }
                else
                {
                    upper_val2 = upper_val2 - star2.getAlpha(k) * basis_val2;
                    lower_val2 = lower_val2 + star2.getAlpha(k) * basis_val2;
                }
            }
            upper_arr.push_back(max(upper_val1, upper_val2));
            lower_arr.push_back(min(lower_val1, lower_val2));
        }

        Point upper_pnt = Point(dimension + 1, upper_arr);
        Point lower_pnt = Point(dimension + 1, lower_arr);
        rt.addUpperBoundState(upper_pnt);
        rt.addLowerBoundState(lower_pnt);
    }
}

void gReachTube::bloatTube()
{
    if (tube.empty())
    {
        cout << "tube is empty, return" << endl;
        return;
    }
    if (bloated)
    {
        cout << "already bloated, return" << endl;
        return;
    }
    for (int i = 0; i < tube.size(); i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            // if(j==10)
            // {
            //     cout<<"stop here"<<endl;
            // }
            vector<double> basis = tube[i].getBasis(j);
            bloatBasis(basis);
            tube[i].setBasis(basis, j);
        }
    }
}

void gReachTube::bloatBasis(vector<double> &basis)
{
    Map<VectorXd> basis_vec(basis.data(), dimension);
    double val = basis_vec.norm();
    if (val != 0)
    {
        basis_vec = basis_vec + basis_vec / val * epsilon * 3;
    }
}