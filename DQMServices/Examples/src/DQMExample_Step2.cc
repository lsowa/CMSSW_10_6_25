#include "DQMServices/Examples/interface/DQMExample_Step2.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <sstream>
#include <math.h>

//
// -------------------------------------- Constructor --------------------------------------------
//
DQMExample_Step2::DQMExample_Step2(const edm::ParameterSet& ps)
{
  edm::LogInfo("DQMExample_Step2") <<  "Constructor  DQMExample_Step2::DQMExample_Step2 " << std::endl;

  // Get parameters from configuration file
  numMonitorName_      =  ps.getParameter<std::string>("numMonitorName");
  denMonitorName_      =  ps.getParameter<std::string>("denMonitorName");

}

//
// -- Destructor
//
DQMExample_Step2::~DQMExample_Step2()
{
  edm::LogInfo("DQMExample_Step2") <<  "Destructor DQMExample_Step2::~DQMExample_Step2 " << std::endl;
}

//
// -------------------------------------- beginJob --------------------------------------------
//
void DQMExample_Step2::beginJob()
{
  edm::LogInfo("DQMExample_Step2") <<  "DQMExample_Step2::beginJob " << std::endl;
}
//
// -------------------------------------- bookHistograms --------------------------------------------
//
void DQMExample_Step2::bookHistograms(DQMStore::IBooker & ibooker_)
{
  std::cout << "DQMExample_Step2::bookHistograms" << std::endl;

  // create and cd into new folder
  ibooker_.setCurrentFolder("What_I_do_in_the_client/Ratio");

  //get available histograms
  MonitorElement* numerator = dbe_->get(numMonitorName_);
  MonitorElement* denominator = dbe_->get(denMonitorName_);

  if (!numerator || !denominator)
    {
      edm::LogError("DQMExample_Step2") <<  "MEs not found!" << std::endl;
      return;
    }


  //book new histogram
  h_ptRatio = ibooker_.book1D("ptRatio","pt ratio pf matched objects",50,0.,100.);
  h_ptRatio->setAxisTitle("pt [GeV]");

  for (int iBin=1; iBin<numerator->getNbinsX(); ++iBin)
    {
      if(denominator->getBinContent(iBin) == 0)
	h_ptRatio->setBinContent(iBin, 0.);
      else
	h_ptRatio->setBinContent(iBin, numerator->getBinContent(iBin) / denominator->getBinContent(iBin));
    }
}

//
// -------------------------------------- endJob --------------------------------------------
//
void DQMExample_Step2::dqmEndJob()
{

  std::cout << "DQMExample_Step2::dqmEndJob" << std::endl;
  edm::LogInfo("DQMExample_Step2") <<  "DQMExample_Step2::dqmEndJob" << std::endl;
}


