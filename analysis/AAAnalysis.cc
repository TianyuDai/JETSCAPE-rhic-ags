#include <iostream>
#include <fstream>
#include <memory>
#include <chrono>
#include <thread>

#include "gzstream.h"
#include "PartonShower.h"
#include "JetScapeLogger.h"
#include "JetScapeReader.h"
#include "JetScapeBanner.h"

#include <GTL/dfs.h>

using namespace std;

using namespace Jetscape;


// -------------------------------------

std::string getcwd_string( void ) {
	   char buff[PATH_MAX];
	      getcwd( buff, PATH_MAX );
	         std::string cwd( buff );
		    return cwd;
}

int main(int argc, char** argv)
{
    int nEvents = 0, i_event=0; 

    vector <double> pTBin{10., 20., 30., 40., 50., 60.};
    vector <double> pTHatBin{60., 100.}; 
    // vector <double> pTHatBin{5., 10., 20., 40., 60., 100.}; 
    vector <double> sigmaGen(pTHatBin.size()-1);  
    vector <double> sigmaErr(pTHatBin.size()-1, 100.);  

    std::ofstream jet_output (argv[2]);
    // std::ofstream jet_output (getcwd_string()+"/../../JETSCAPE-output/AuAu200/AuAu200_chargedHadron.txt");
    std::vector <double> hadron_cs(pTBin.size()-1, 0.), hadron_cs_err(pTBin.size()-1, 0.), pTAvg(pTBin.size()-1, 0.), total_sigma(pTBin.size()-1, 0.); 
    std::vector<std::vector <int>> hadron_ct(pTHatBin.size()-1);
    std::vector<std::vector <int>> hadron_ct_sq(pTHatBin.size()-1);
    for (unsigned int iBin = 0; iBin < pTHatBin.size() - 1; iBin++)
    {
	std::cout << argv[1]+std::to_string(pTHatBin[iBin]) << "\n"; 
        // auto reader=make_shared<JetScapeReaderAscii>(getcwd_string()+"/../../JETSCAPE-output/AuAu200/"+std::to_string(pTHatBin[iBin])+".dat"); 
        auto reader=make_shared<JetScapeReaderAscii>(argv[1]+std::to_string(pTHatBin[iBin])+".dat");  
        hadron_ct[iBin] = std::vector<int>(pTBin.size()-1, 0); 
        hadron_ct_sq[iBin] = std::vector<int>(pTBin.size()-1, 0);
	std::vector<double> pTSum(pTBin.size()-1, 0.);
        while (!reader->Finished())
	{
            std::vector<int> hadron_ct_s(pTBin.size()-1, 0); 
	    reader->Next();

	    i_event = reader->GetCurrentEvent(); 
	    auto sigma_gen = reader->GetSigmaGen(); 
	    auto sigma_err = reader->GetSigmaErr(); 
	    if (sigma_err < sigmaErr[iBin])
            {
                sigmaGen[iBin] = sigma_gen; 
		sigmaErr[iBin] = sigma_err; 
	    }

            auto hadrons = reader->GetHadrons(); 
	    auto pdghelper = JetScapeParticleBase::InternalHelperPythia.particleData; 
            for (unsigned int iHadron = 0; iHadron < hadrons.size(); iHadron++)
            {
                if (hadrons[iHadron]->pt() < pTBin[0]) continue; 
                if (fabs(hadrons[iHadron]->eta()) > 1.) continue; 
                auto ID = hadrons[iHadron]->pid(); 
                auto charge = pdghelper.charge( ID );
                if (charge == 0) continue; 
		for (unsigned int ipT = 0; ipT < pTBin.size()-1; ipT++)
		    if (hadrons[iHadron]->pt() > pTBin[ipT] && hadrons[iHadron]->pt() <= pTBin[ipT+1])
		    {
		        hadron_ct_s[ipT]++; 
                        pTSum[ipT] += hadrons[iHadron]->pt(); 
		        break; 
		    }
            }
            for (unsigned int ipT = 0; ipT < pTBin.size() - 1; ipT++)
            {
                hadron_ct[iBin][ipT] += hadron_ct_s[ipT]; 
                hadron_ct_sq[iBin][ipT] += hadron_ct_s[ipT] * hadron_ct_s[ipT]; 
            }
	}
	nEvents = i_event + 1;
	for (unsigned int ipT = 0; ipT < pTBin.size() - 1; ipT++)
	{
            if (hadron_ct[iBin][ipT] > 0)
            {
               pTSum[ipT] /= hadron_ct[iBin][ipT]; 
               total_sigma[ipT] += sigmaGen[iBin]; 
            }
            pTAvg[ipT] += pTSum[ipT] * sigmaGen[iBin]; 
	    hadron_cs[ipT] += double(hadron_ct[iBin][ipT]) / nEvents * sigmaGen[iBin];
            double err_ipT; 
            err_ipT = sqrt((double(hadron_ct_sq[iBin][ipT])/nEvents - pow(double(hadron_ct[iBin][ipT])/nEvents, 2))/nEvents);
            hadron_cs_err[ipT] += sqrt(sigmaErr[iBin]*sigmaErr[iBin] * (double(hadron_ct[iBin][ipT])/nEvents)*(double(hadron_ct[iBin][ipT])/nEvents) + sigmaGen[iBin]*sigmaGen[iBin]*err_ipT*err_ipT); 
	}
	reader->Close(); 
    }
    for (unsigned int ipT = 0; ipT < pTBin.size()-1; ipT++)
	jet_output << pTAvg[ipT] << " " << hadron_cs[ipT] / (pTBin[ipT+1] - pTBin[ipT]) / pTAvg[ipT] << " " << hadron_cs_err[ipT] / (pTBin[ipT+1] - pTBin[ipT]) / pTAvg[ipT] << endl; 

	cout<<"Finished!"<<endl; 
}
