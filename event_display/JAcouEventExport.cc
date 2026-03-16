#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>

#include "TROOT.h"
#include "TFile.h"

#include "JLang/JPredicate.hh"
#include "JLang/JComparator.hh"
#include "JLang/JComparison.hh"

#include "JTools/JQuantile.hh"

#include "JSupport/JMultipleFileScanner.hh"
#include "JSupport/JFileRecorder.hh"
#include "JSupport/JMeta.hh"

#include "JASPER/JSupport.hh"     
#include "JASPER/JAcouHit.hh"         
#include "JASPER/JAcouEvent.hh"        

#include "Jeep/JContainer.hh"
#include "Jeep/JParser.hh"
#include "Jeep/JMessage.hh"


/**
 * \file
 *
 * Example application to export acoustics event.
 * \author hbaetsen
 */
int main(int argc, char **argv)
{
  using namespace std;
  using namespace JPP;
  using namespace JASPER;

  JMultipleFileScanner<JAcouEvent>  inputFile;
  JLimit_t&                         numberOfEvents = inputFile.getLimit();
  int                               debug;
  string                            outputFile; 

  try {

    JParser<> zap("Example application to print acoustics event.");

    zap['f'] = make_field(inputFile)        ;
    zap['n'] = make_field(numberOfEvents)   = JLimit::max();
    zap['d'] = make_field(debug)            = 2;
    zap['o'] = make_field(outputFile)       = "EventOutFileTest.root";

    zap(argc, argv);
  }
  catch(const exception &error) {
    FATAL(error.what() << endl);
  }

  TFile out(outputFile.c_str(), "recreate");
  TTree *tree = new TTree("tree", "Hits");

  Double_t EventID, HydrophoneID, TemplateID, ToA, SNR;

  tree->Branch("EventID",       &EventID);
  tree->Branch("HydrophoneID",  &HydrophoneID);
  tree->Branch("TemplateID",    &TemplateID);
  tree->Branch("ToA",           &ToA);
  tree->Branch("SNR",           &SNR);

  while (inputFile.hasNext()) {

    const JAcouEvent* evt = inputFile.next();

    // map<int, vector<JAcouTransmission>> buffer; // vector<JAcouTransmission>

    EventID = evt->getCounter();
    STATUS("EVENT ID = " << EventID << endl);

    for (JAcouEvent::const_iterator hit = evt->begin(); hit != evt->end(); ++hit) {

      // buffer[hit->getID()].push_back(*hit);
      HydrophoneID  = hit->getID();
      TemplateID    = hit->getTemplateID();
      ToA           = hit->getToA();
      SNR           = hit->getSNR();

      tree->Fill();
    }
    
    // size_t duplicates = 0;
    
    // for (map<int, vector<JAcouTransmission> >::iterator ps = buffer.begin(); ps != buffer.end(); ++ps) {

    //   sort(ps->second.begin(), ps->second.end(), make_comparator(&JAcouTransmission::getSNR, JComparison::gt()));

	  //   duplicates += 1;
    // } 
  }
  
  out.Write();
  out.Close();
}
