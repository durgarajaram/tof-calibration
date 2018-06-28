//Y.Karadzhov 02.2009

#ifndef TOF_T0_CALIB
#define TOF_T0_CALIB

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TSpectrum.h"
#include "TMath.h"
#include "MiceTofCalib.hh"

#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Units/PhysicalConstants.h"
using CLHEP::picosecond;
using CLHEP::cm;
using CLHEP::mm;
using CLHEP::c_light;

class TofTWCalib;
class TofTriggerCalib;
using namespace std;

class TofT0Calib : public MiceTofCalib
{
 public:
  TofT0Calib(int station, TofTWCalib*, TofTriggerCalib*);
  TofT0Calib(int station, TofTWCalib*, TofTriggerCalib* , TTree*);
  virtual ~TofT0Calib();

  void   SetTW(TofTWCalib *tw)  {TWCalib = tw;};
  void   SetTrigger(TofTriggerCalib *tr)  {TriggerCalib = tr;};
  void   SetdT(int dt) {dT = dt;};
  void   SetL(int l) {L = l; dT = int( l/(c_light*picosecond) );};
  void   SetNBins(int n)  { NHistoBins = n; };

  double GetL()   { return L; };
  int    GetdT()  { return dT; };
  int    T0Correction(int station, int plane, int slab, int pmt);
  void   FullCalib(int refSlA, int refSlB, bool verboseMode);
  void   PlaneCalib(int plane, int refSl, bool verboseMode);
  void   Test();
  void   DevTest();
  void   PrintPMTParameters();

 private:

  void MakeHistograms();
  void FillHistograms();

  TofTWCalib *TWCalib;
  TofTriggerCalib *TriggerCalib;
  int Station;
  int nSlabs;
  int dT;
  double L;
  int NHistoBins;

  vector< vector< vector<TH1F*> > > hists;
  
  TCanvas* c1_t0;
};

#endif


