#include <iostream>
#include <fstream>
#include <time.h>

// JetScape Framework includes ...
#include "JetScape.h"
#include "JetEnergyLoss.h"
#include "JetEnergyLossManager.h"
#include "JetScapeWriterStream.h"
#ifdef USE_HEPMC
#include "JetScapeWriterHepMC.h"
#endif

// User modules derived from jetscape framework clasess
#include "TrentoInitial.h"
#include "AdSCFT.h"
#include "Matter.h"
#include "Martini.h"
#include "Tequila.h"
#include "NullPreDynamics.h"
#include "FreestreamMilneWrapper.h"
#include "Brick.h"
#include "GubserHydro.h"
#include "HydroFromFile.h"
#include "MusicWrapper.h"
#include "PGun.h"
#include "PythiaGun.h"
#include "HadronizationManager.h"
#include "Hadronization.h"
#include "ColoredHadronization.h"
#include "ColorlessHadronization.h"
#include "HybridHadronization.h"

#include <chrono>
#include <thread>

using namespace Jetscape;
// Forward declaration
void Show();

// -------------------------------------
void RunEvents(double scale, double alpha_s, double omegacut, int N)
{

    //modify the init.xml file
    JetScapeXML::Instance()->OpenXMLFile("for_hydro_jetscape_init.xml");
    // JetScapeXML::Instance()->OpenXMLFile("AA_init.xml");
    tinyxml2::XMLElement *scalexml=JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" )->FirstChildElement("Tequila" )->FirstChildElement("muqperp_over_T"); 
    scalexml->SetText(scale);
    tinyxml2::XMLElement *alphaxml= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" )->FirstChildElement("Tequila" )->FirstChildElement("alpha_s");
    alphaxml->SetText(alpha_s);
    tinyxml2::XMLElement *omegaxml= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" )->FirstChildElement("Tequila" )->FirstChildElement("muomega_over_T");
    omegaxml->SetText(omegacut); 

    // JetScapeXML::Instance()->CloseXMLFile();
	  
    double deltaT = 0.01; 
    tinyxml2::XMLElement *dtxml= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" )->FirstChildElement("deltaT" );
    dtxml->SetText(deltaT);

    double pT_hat_min, pT_hat_max; 
    tinyxml2::XMLElement *pythiaxml= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Hard" )->FirstChildElement("PythiaGun");
    pythiaxml->FirstChildElement("pTHatMin")->QueryDoubleText(&pT_hat_min);
    pythiaxml->FirstChildElement("pTHatMax")->QueryDoubleText(&pT_hat_max);
    auto jetscape = make_shared<JetScape>("for_hydro_jetscape_init.xml", N); 
    // auto jetscape = make_shared<JetScape>("AA_init.xml", N); 
    jetscape->SetId("primary");
    jetscape->SetReuseHydro (true);
    jetscape->SetNReuseHydro (1000000);

    // Initial conditions and hydro
    auto trento = make_shared<TrentoInitial> ();
    auto preeqdyn = make_shared<NullPreDynamics> ();
    auto freestream = make_shared<FreestreamMilneWrapper> ();
    auto pythiaGun = make_shared<PythiaGun> ();
    auto pGun = make_shared<PGun> ();
    auto hydro = make_shared<MpiMusic> ();
    jetscape->Add(trento);
    jetscape->Add(preeqdyn); 
                // jetscape->Add(freestream); 
		jetscape->Add(pythiaGun);
		// jetscape->Add(pGun);
		jetscape->Add(hydro);

		// Energy loss
		auto jlossmanager = make_shared<JetEnergyLossManager> ();
		auto jloss = make_shared<JetEnergyLoss> ();
		auto matter = make_shared<Matter> (); 
		auto tequila = make_shared<Tequila> ();
		// auto martini = make_shared<Martini> (); 

		jloss->Add(matter); 
		// jloss->Add(tequila);
		// jloss->Add(martini); 
		jlossmanager->Add(jloss);  
		jetscape->Add(jlossmanager);

		auto printer = make_shared<PartonPrinter> ();
		jetscape->Add(printer);

		// Hadronization
		// This helper module currently needs to be added for hadronization.
                auto hadroMgr = make_shared<HadronizationManager> ();
                auto hadro = make_shared<Hadronization> ();
                auto hadroModule = make_shared<ColoredHadronization> ();
                auto colorless = make_shared<ColorlessHadronization> ();
                hadro->Add(colorless);
                //hadro->Add(hybridHadr);
                hadroMgr->Add(hadro);
                jetscape->Add(hadroMgr);


                // Output
		// auto writer= make_shared<JetScapeWriterAscii> (("../build/test/"+std::to_string(pT_hat_min)+".dat").c_str());
		auto writer= make_shared<JetScapeWriterAscii> ("test_out.dat");

		jetscape->Add(writer);

		jetscape->Init();
		jetscape->Exec();
		jetscape->Finish();
}

int main(int argc, char** argv)
{
	clock_t t; t = clock();
	time_t start, end; time(&start);

	cout<<endl;

	// DEBUG=true by default and REMARK=false
	// can be also set also via XML file (at least partially)
	JetScapeLogger::Instance()->SetInfo(true);
	JetScapeLogger::Instance()->SetDebug(false);
	JetScapeLogger::Instance()->SetRemark(false);
	//SetVerboseLevel (9 adds a lot of additional debug output)
	//If you want to suppress it: use SetVerboseLevle(0) or max  SetVerboseLevle(9) or 10
	JetScapeLogger::Instance()->SetVerboseLevel(0);

  	Show();

	double scale_list[1] = {1.}; 
  	double alpha_list[1] = {0.3}; 
  	double omegacut_list[1] = {1.}; 
	for (int i = 0; i < 1; i++)
  	{
  	  
	  double scale = scale_list[i]; 
	  for (int j = 0; j < 1; j++)
	  {
	  	double alpha_s = alpha_list[j]; 
	  	for (int k = 0; k < 1; k++)
	  	{
	  		double omegacut = omegacut_list[k]; 
	  		RunEvents(scale, alpha_s, omegacut, 10000); 
	  	}
	  }
  	}


	INFO_NICE<<"Finished!";
	cout<<endl;

	t = clock() - t;
	time(&end);
	printf ("CPU time: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);
	printf ("Real time: %f seconds.\n",difftime(end,start));
	return 0;
}

// -------------------------------------

void Show()
{
	INFO_NICE<<"------------------------------------";
	INFO_NICE<<"| Brick Test JetScape Framework ... |";
	INFO_NICE<<"------------------------------------";
	INFO_NICE;
}
