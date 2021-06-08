#ifndef CUSTOMMODULE_H
#define CUSTOMMODULE_H

#include <vector>
#include "lorentz.h"
#include "JetEnergyLossModule.h"
#include "JetScapeConstants.h"

using namespace Jetscape; 

class CustomModule : public JetEnergyLossModule<CustomModule> //, public std::enable_shared_from_this<CustomModule>
{
 private: 
	double alpha_s;
  	double g;
  	double pcut;        // below this scale, no further Eloss
  	double M = 0.; 
  	double mu; 
	std::normal_distribution<double> white_noise;
	std::default_random_engine generator;

        static RegisterJetScapeModule<CustomModule> reg;

 public:  
  	CustomModule();
  	virtual ~CustomModule();

  	//main//
  	void Init();
  	void DoEnergyLoss(double deltaT, double Time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut);
 
	double qperp(double E, double T, int id);
	double qpara(double E, double T, int id);
	void Ito_update( double dt, double T, std::vector<double> v, const fourvec & pIn, fourvec & pOut, int id);
}; 
#endif

