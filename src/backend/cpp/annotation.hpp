/*
 * \file   annotation.hpp
 * \class  Annotation
 * 
 * \author parasara
 * \author Lucas Brown
 * \date   Feb 1, 2013
 * \date   March 29, 2019
 * 
 * \brief LMBTODO 
 */

#ifndef ANNOTATION_H_
#define ANNOTATION_H_

#include <string>
#include <vector>

class Annotation 
{
  public:
	Annotation();
	virtual ~Annotation();

    double getED(double delta, double t1, double t2, int mode);

    int getType();
	void setType(int val);
	
    std::string getAnnotation();
	void setAnnotation(std::string val);
	std::string getBeta();
	void setBeta(std::string val);

    double getKVal(int mode);
    void setKVal(int mode, double val);    
    std::vector<double> getKConsts();
    void setKConsts(std::vector<double> vec);

    double getGammaVal(int mode);
    void setGammaVal(int mode, double val);
    std::vector<double> getGammas();
    void setGammas(std::vector<double> vec);
	
  private:
  
	int type;
	std::string Annot;
	std::string Beta;

    std::vector<double> k_consts;
    std::vector<double> gammas;
};

#endif /* ANNOTATION_H_ */
