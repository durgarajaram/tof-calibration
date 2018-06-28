//Y.Karadzhov 02.2009

#ifndef TOF_TW_CALIB
#define TOF_TW_CALIB

#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCut.h"
#include "TMinuit.h"
#include "TVirtualFitter.h"
#include "TPaveStats.h"


#include "MiceTofCalib.hh"

class TofTriggerCalib;
class TofT0Calib;
using namespace std;

class TofTWCalib : public MiceTofCalib
{
 public:

  TofTWCalib(string n="");
  TofTWCalib(TTree*, string n="");
  virtual ~TofTWCalib();

  void SetTrigger(TofTriggerCalib *trigger) {triggerCalib = trigger;};
  void SetT0(TofT0Calib *t0) {t0Calib = t0;};
  void SetName(string n) {_name = n;};
  bool ReferenceSlabCalib(int station, int slabA, int slabB, bool verboseMode, bool Minuit);
  bool SlabCalib(int station, int plane, int slab, int refslab, bool verboseMode);
  void FullCalib(int station, int refSlA,int refSlB, bool verboseMode, bool Minuit);
  void ResetHistograms();
  double GetResolution(int st ,int slabA, int slabB, bool verboseMode, int& Entris);
  int GetNPar() { return nPar; };
  int TWCorrection(double adc, double* par);
  int TWCorrection(double adc, int station, int plane, int slab, int pmt);

  TH2F* GetPMTHistogram(int station, int plane, int slab, int pmt);
  TProfile* GetPMTProfile(int station, int plane, int slab, int pmt);
  TH1F* GetRawTHistogram(int station, int plane, int slab, bool RefCalib);
  TH1F* GetResolHistogram(int station, int plane, int slab, bool Minuit);
  void PrintPMTParameters();
  void SetAdcCut(double cut)  { if(cut<1.) adcCut = cut; };

  void FillHistograms(int st, int refSlA, int refSlB);
 private:

  void InitializeMinuit();
  void MakeHistograms();
  void ClearHistograms();

  void FillRefHistograms(int st, int refSlA, int refSlB);

  void RootFit(int Station, int Plane, int Slab, int Pmt, bool verboseMode);
  TF1* RootFit(TProfile* pfh, bool verboseMode);
  void MinuitFit(int Station, int Plane, int Slab, int Pmt, bool verboseMode);
  TF1* MinuitFit(TProfile* pfh, double adc_mean, double adc_rms, bool verboseMode);

  void CompareFits(int Station, int Plane, int Slab, int refSl);
  TProfile* MakeProfile(int Station, int Plane, int Slab, int Pmt);
  TProfile* MakeProfile(TH2F* h2);
  void DrawMaxAdc(TH2F* h2);
  int nPar;
  double adcCut;
  double ReferenceSlabPar[4][5];
  double PmtPar[4][5];
  TF1 *f1;

  TofTriggerCalib *triggerCalib;
  TofT0Calib *t0Calib;
  TMinuit * minuit;
  double arglist[10];
  int ierflg;

  string _name;
  TCanvas* c1_tw;

  vector< vector< vector< TH1F* > > > raw_t;
  vector< vector< vector< TH1F* > > > t1;
  vector< vector< vector< TH1F* > > > t2;
  vector< vector< vector< TH1F* > > > t2m;

  vector< vector< vector<  TH2F*  > > > R_hists;
  vector< vector< vector<TProfile*> > > R_profiles;

  vector< vector< vector< vector<  TH2F*  > > > > hists;
  vector< vector< vector< vector<TProfile*> > > > profiles;
};

#endif

