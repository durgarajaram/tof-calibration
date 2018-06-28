//Y.Karadzhov 02.2009

#ifndef TOF_TRIGGER_CALIB
#define TOF_TRIGGER_CALIB

#include "TF1.h"
#include "TH1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TApplication.h"

#include "MiceTofCalib.hh"
#include "src/common_cpp/Utils/TOFCalibrationMap.hh"

class TofTWCalib;
using namespace std;

class TofTriggerCalib : public MiceTofCalib
{
 public:

  TofTriggerCalib(int, TofTWCalib*);
  TofTriggerCalib(int, TofTWCalib*, TTree*);
  ~TofTriggerCalib();

  bool TriggerCalib(int Plane, int Slab, double refTime, int refSl, bool verboseMode);
  void FullCalib(int refSlX, int refSlY, bool verboseMode);
  void MakeHistograms();
  void MakeTriggerTree();
  void FillHistograms();
  void SetTW(TofTWCalib *tw)  {TWCalib = tw;};
  void SetTriggerTree(TTree *t) {triggerT0Tree = t;};
  TTree* GetTriggerTree()  {return triggerT0Tree;};
  void SetTriggerTxtFile(ofstream *f)    {Triggerfile = f;};

  void SetStation(int st)   { Station = st; };
  int GetStation() {return Station;};

  TH1F* GetPMTHistogram(int plane, int slab, int refslab, int pmt);
  TH1F* GetDeltaTHistogram(int plane, int slab, int refslab);
  void SetHistoRange(TH1F *hist);

  int GetNCalibrated() const {return nGoodPixels;};

  int GetTriggerT0(int slabX, int slabY);
  int GetTriggerT0(TTree *tree, int slabX, int slabY);
  int GetChannelT0(int station, int plane, int slab, int pmt, int &refslab);
  int T0Correction(int station, int plane, int slab, int pmt, int crossSlab);

  void Test();
  void PrintTriggerDelays();
  void SaveTriggerDelay(int slX, int slY, int tr_t0);
  bool LoadTriggerDelays(string fileName);
  void PrintPMTParameters();
  void MarkSlab(int pl, int sl) {MarkedSlab.plane = pl; MarkedSlab.slab = sl;};

 private:

  double refA_mean_L;
  double refA_mean_R;
  double refA_mean_LmR;
  double refB_mean_L;
  double refB_mean_R;
  double refB_mean_LmR;
  int Station;
  int nSlabs;
  int nGoodPixels;

  vector< vector< vector<TH1F*> > > tr_t0;
  vector< vector< vector<TH1F*> > > t_LmR;
  vector< vector< double > > tr_delays;

  TofTWCalib *TWCalib;
  TriggerT0  _trT0;
  TTree *triggerT0Tree;
  ofstream *Triggerfile;
  
  pmtParam MarkedSlab;
  TCanvas* c1_tr;
};
#endif


