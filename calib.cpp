
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TApplication.h"
#include "TStyle.h"

#include "MiceTofCalib.hh"
#include "TofT0Calib.hh"
#include "TofTriggerCalib.hh"
#include "TofTWCalib.hh"

void doTriggerCalib();
void doTWCalib();
void doStationsCalib();
bool getParams();

TofTWCalib *TWcalib = NULL;
TofTriggerCalib *TofTrigger = NULL;
TofT0Calib *Tof1calib = NULL;
TofT0Calib *Tof2calib = NULL;

const bool twCalib = true;
const bool triggerCalib = true;
const bool t0Calib = true;
const bool verboseMode = false;
unsigned int trStation;
int tof0_sXmax, tof0_sYmax, tof1_sXmax, tof1_sYmax, tof2_sXmax, tof2_sYmax;
int min_entries_tw, min_entries_trig, min_entries_t0;
int tof1_len, tof2_len;

TFile *datafile;
TFile *histofile;

int main(int argc, char **argv) {
  int one = 1;

  // first read the configuration data cards
  // we cannot proceed unless we have the parameters
  if (! getParams()) return 1;
  
  std::cout << "got Params: " << trStation << " " << tof1_sXmax << std::endl;
  TApplication *app;
  if(verboseMode)
    app = new TApplication("appKey",&one,argv);

  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(1);
  
  string fileName(argv[1]);
  datafile = new TFile(fileName.c_str(), "update");
  histofile = new TFile("tofcalib_histos.root", "update");

  datafile->ls();
  TTree* dataTree = (TTree*)datafile->Get("dataTree");
  //dataTree->Print();

  if (twCalib) {
    TWcalib = new TofTWCalib(dataTree);
    doTWCalib();
  } else {
    TWcalib = new TofTWCalib();
    TWcalib->LoadParameters("tofTWcalib.txt", TWcalib->GetNPar(),-1);
    cout<<"Time Walk calibration - loaded"<<endl;
    //TWcalib->PrintPMTParameters();
  }

  if (triggerCalib) {
    TofTrigger = new TofTriggerCalib(trStation, TWcalib, dataTree);
    doTriggerCalib();
  } else {
    TofTrigger = new TofTriggerCalib(trStation,TWcalib);
    TofTrigger->LoadParameters("tofT0calib.txt", 2, trStation);
    TofTrigger->LoadTriggerDelays("tofTriggercalib.txt");
    cout<<"Trigger calibration - loaded"<<endl;
    //TofTrigger->PrintPMTParameters();
    //TofTrigger->PrintTriggerDelays();
  }

  if (t0Calib) {
    // arg 1 = station
    Tof1calib = new TofT0Calib(0, TWcalib, TofTrigger, dataTree);
    Tof2calib = new TofT0Calib(2, TWcalib, TofTrigger, dataTree);
    doStationsCalib();
  }

  histofile->Close();
  datafile->Close();

  return 0;
}

void doTWCalib() {
  ofstream f_tw("tofTWcalib.txt");
  TWcalib->SetTxtFile(&f_tw);
  TWcalib->SetMinEntries(min_entries_tw);

  histofile->cd();
  TDirectory *tw;
  histofile->GetObject("TW",tw);
  if(tw) histofile->rmdir("TW");
  tw = histofile->mkdir("TW","Time Walk");
  tw->cd();
  TWcalib->FullCalib(0,tof0_sXmax,tof0_sYmax,verboseMode,true); // mod
  TWcalib->FullCalib(1,tof1_sXmax,tof1_sYmax,verboseMode,true);
  TWcalib->FullCalib(2,tof2_sXmax,tof2_sYmax,verboseMode,false);
}

void doTriggerCalib() {
  ofstream f_trigger("tofTriggercalib.txt");
  ofstream f_t0("tofT0calib.txt");
  histofile->cd();
  TDirectory *trigger;
  histofile->GetObject("Trigger",trigger);
  if(trigger) histofile->rmdir("Trigger");
  trigger = histofile->mkdir("Trigger","Trigger Calibration");
  trigger->cd();

  TofTrigger->SetTxtFile(&f_t0);
  TofTrigger->SetTriggerTxtFile(&f_trigger);
  //TofTrigger->MarkSlab(0,4);
  TofTrigger->SetMinEntries(min_entries_trig);

  if(trStation==0) TofTrigger->FullCalib(tof0_sXmax,tof0_sYmax,verboseMode);
  if(trStation==1) TofTrigger->FullCalib(tof1_sXmax,tof1_sYmax,verboseMode);
  TofTrigger->Test();  
}

void doStationsCalib() {
  //Tof1calib->SetL(-7693.4*mm);
  Tof1calib->SetL(tof1_len*mm);
  Tof1calib->SetMinEntries(min_entries_t0);

  //Tof2calib->SetL(2500.0*mm);
  Tof2calib->SetL(tof2_len*mm);
  Tof2calib->SetMinEntries(min_entries_t0);

  //Tof1calib->SetNBins(500);
  ofstream f_t0("tofT0calib.txt");

  TofTrigger->SetTxtFile(&f_t0);
  TofTrigger->SaveAllToFile(2);
  Tof1calib->SetTxtFile(&f_t0);
  Tof2calib->SetTxtFile(&f_t0);

  histofile->cd();
  TDirectory *tof;
  histofile->GetObject("TOF",tof);
  if(tof) histofile->rmdir("TOF");
  tof = histofile->mkdir("TOF","TOF Calibration");
  tof->cd();
  Tof1calib->FullCalib(tof0_sXmax,tof0_sYmax,verboseMode);
  Tof2calib->FullCalib(tof1_sXmax,tof1_sYmax,verboseMode); // mod
  Tof1calib->Test();
  Tof2calib->Test();
  //Tof1calib->DevTest();
}

bool getParams() {
    ifstream inputFile("tofcalib_cards.txt");
    string line;
    std::vector<std::string> keystr; 
    std::map<std::string, std::string> cfgmap;
    std::map<string, string>::iterator it;
    while (getline(inputFile, line))
    {
        istringstream ss(line);
        std::string key;
        if (std::getline(ss, key, '='))
        {
	    key.erase( key.find_last_not_of( " \f\t\v" ) + 1 );
	    keystr.push_back(key);
            std::string value;
            if (std::getline(ss, value)) {
	      value.erase( 0, value.find_first_not_of( " \f\t\v" ));
	      cfgmap[key] = value;
            }
        }
    }
    trStation = atoi(cfgmap["trigger_station"].c_str());
    tof0_sXmax = atoi(cfgmap["tof0_sXmax"].c_str());
    tof0_sYmax = atoi(cfgmap["tof0_sYmax"].c_str());
    tof1_sXmax = atoi(cfgmap["tof1_sXmax"].c_str());
    tof1_sYmax = atoi(cfgmap["tof1_sYmax"].c_str());
    tof2_sXmax = atoi(cfgmap["tof2_sXmax"].c_str());
    tof2_sYmax = atoi(cfgmap["tof2_sYmax"].c_str());
    min_entries_tw = atoi(cfgmap["min_entries_tw"].c_str());
    min_entries_trig = atoi(cfgmap["min_entries_trig"].c_str());
    min_entries_t0 = atoi(cfgmap["min_entries_t0"].c_str());
    tof1_len = atoi(cfgmap["tof1_len"].c_str());
    tof2_len = atoi(cfgmap["tof2_len"].c_str());
    return true;
}

