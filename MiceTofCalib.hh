//Y.Karadzhov 02.2009

#ifndef TOF_CALIB
#define TOF_CALIB

#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "TTree.h"
#include "TCanvas.h"
#include "TH2I.h"
#include "TPaveStats.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TColor.h"

#include "src/common_cpp/Utils/TOFChannelMap.hh"

using namespace std;
  typedef struct
  {
    int slabA;
    int slabB;
    double t0;
    double t1;
    double t2;
    double t3;
    int adc0;
    int adc1;
    int adc2;
    int adc3;

    int slabC;
    int slabD;
    double t4;
    double t5;
    double t6;
    double t7;
    int adc4;
    int adc5;
    int adc6;
    int adc7;

    int slabE;
    int slabF;
    double t8;
    double t9;
    double t10;
    double t11;
    int adc8;
    int adc9;
    int adc10;
    int adc11;
  } tofData;

  typedef struct
  {
    int st;
    int plane;
    int slab;
    int pmt;
    double par[5];
  } pmtParam;

  typedef struct
  {
    int st;
    int slabX;
    int slabY;
    int TrDelay;
  } TriggerT0;

class MiceTofCalib
{
  public:

  MiceTofCalib();
  ~MiceTofCalib();

  void SetDataTree(TTree *t);
  void SetParamTree(TTree *t)     {channelParTree = t;};

  void SetAddresses(int Station);
  void SetTxtFile(ofstream *f)    {file = f;};
  void SetMinEntries(int n)       {minEntries = n;};

  TTree* GetDataTree() { return dataTree; };
  TTree* GetParamTree() { return channelParTree; };
  double* GetParameters(int station, int plane, int slab, int pmt);

  void PrintPMTParameters(int);
  void SaveToFile(int st, int pl, int sl, int pmt, double *par, int npar);
  void SaveAllToFile(int nPar);
  bool LoadParameters(string fileName, int npar, int station);
  void AutoSave()  { dataTree->AutoSave(); };

  void InitDataMonitor();
  void DrawDataMonitor();

  protected:

  void MakeDataTree();
  void MakeParamTree();

  double pmtPar[5];
  TTree *dataTree;
  TTree *channelParTree;
  
  tofData _data;
  pmtParam _pmtPar;
  ofstream *file;
  int minEntries;

  int _slX, _slY, _adc0, _adc1, _adc2, _adc3;
  double _t0, _t1, _t2, _t3;
  int buffer0, buffer1,buffer2,buffer3,buffer4,buffer5,buffer6,buffer7,buffer8,buffer9;
  int buffer10, buffer11,buffer12,buffer13,buffer14,buffer15,buffer16,buffer17,buffer18,buffer19;

  TCanvas *c1_mon;
  TH2I tof0;
  TH2I tof1;
  TH2I tof2;
};
#endif
