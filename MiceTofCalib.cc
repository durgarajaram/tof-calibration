//Y.Karadzhov 02.2009

#include "MiceTofCalib.hh"

MiceTofCalib::MiceTofCalib()
{
  MakeDataTree();
  MakeParamTree();

  minEntries = 1000;

  for(int i=0;i<5;i++)
    _pmtPar.par[i] = 0.;

}

MiceTofCalib::~MiceTofCalib()
{
  //delete dataTree;
  //delete channelParTree;
  //delete c1_mon;
  if( tof0.GetEntries() )
  {
    tof0.Write();
    tof1.Write();
    tof2.Write();
  }
}

void MiceTofCalib::MakeDataTree()
{
  dataTree = new TTree("dataTree","tofdata");
  dataTree->Branch("slabA", &_data.slabA, "slabA/I");
  dataTree->Branch("slabB", &_data.slabB, "slabB/I");
  dataTree->Branch("t0", &_data.t0, "t0/I");
  dataTree->Branch("t1", &_data.t1, "t1/I");
  dataTree->Branch("t2", &_data.t2, "t2/I");
  dataTree->Branch("t3", &_data.t3, "t3/I");
  dataTree->Branch("adc0", &_data.adc0, "adc0/I");
  dataTree->Branch("adc1", &_data.adc1, "adc1/I");
  dataTree->Branch("adc2", &_data.adc2, "adc2/I");
  dataTree->Branch("adc3", &_data.adc3, "adc3/I");

  dataTree->Branch("slabC", &_data.slabC, "slabC/I");
  dataTree->Branch("slabD", &_data.slabD, "slabD/I");
  dataTree->Branch("t4", &_data.t4, "t4/I");
  dataTree->Branch("t5", &_data.t5, "t5/I");
  dataTree->Branch("t6", &_data.t6, "t6/I");
  dataTree->Branch("t7", &_data.t7, "t7/I");
  dataTree->Branch("adc4", &_data.adc4, "adc4/I");
  dataTree->Branch("adc5", &_data.adc5, "adc5/I");
  dataTree->Branch("adc6", &_data.adc6, "adc6/I");
  dataTree->Branch("adc7", &_data.adc7, "adc7/I");

  dataTree->Branch("slabE", &_data.slabE, "slabE/I");
  dataTree->Branch("slabF", &_data.slabF, "slabF/I");
  dataTree->Branch("t8", &_data.t8, "t8/I");
  dataTree->Branch("t9", &_data.t9, "t9/I");
  dataTree->Branch("t10", &_data.t10, "t10/I");
  dataTree->Branch("t11", &_data.t11, "t11/I");
  dataTree->Branch("adc8", &_data.adc8, "adc8/I");
  dataTree->Branch("adc9", &_data.adc9, "adc9/I");
  dataTree->Branch("adc10", &_data.adc10, "adc10/I");
  dataTree->Branch("adc11", &_data.adc11, "adc11/I");
}

void MiceTofCalib::MakeParamTree()
{
  channelParTree = new TTree("channelParTree","PmtParameters");
  channelParTree->Branch("station", &_pmtPar.st, "station/I");
  channelParTree->Branch("plane", &_pmtPar.plane, "plane/I");
  channelParTree->Branch("slab", &_pmtPar.slab, "slab/I");
  channelParTree->Branch("pmt", &_pmtPar.pmt, "pmt/I");
  channelParTree->Branch("par", &_pmtPar.par, "par[5]/D");
}

void MiceTofCalib::SetDataTree(TTree *t)
{
  if(dataTree) dataTree->Delete();
  //dataTree = t->CloneTree();
  dataTree = t;
}

double* MiceTofCalib::GetParameters(int station, int plane, int slab, int pmt)
{
  int Station, Plane, Slab, Pmt;
  double Par[5];

  //channelParTree->Scan();
  channelParTree->GetBranch("station")->SetAddress(&Station);
  channelParTree->GetBranch("plane")->SetAddress(&Plane);
  channelParTree->GetBranch("slab")->SetAddress(&Slab);
  channelParTree->GetBranch("pmt")->SetAddress(&Pmt);
  channelParTree->GetBranch("par")->SetAddress(&Par);

  int Par_entries = channelParTree->GetEntries();
  int i = 0;

  while(i<Par_entries)
  {
    channelParTree->GetEntry(i);
    i++;

    if(station==Station && plane == Plane && slab==Slab && pmt==Pmt)
    {
      //cout<<"St"<<station<<"Pl"<<plane<<"Sl"<<slab<<"Pmt"<<pmt<<" ";
      for(int i=0;i<5;i++)
      {
        //cout<<Par[i]<<" ";
        pmtPar[i] = Par[i];
      }
      //cout<<endl;
      return pmtPar;
    }
  }

  return NULL;
}

void MiceTofCalib::PrintPMTParameters(int nPar)
{
  int Station, Plane, Slab, Pmt;
  double Par[5];

  //channelParTree->Scan();
  channelParTree->GetBranch("station")->SetAddress(&Station);
  channelParTree->GetBranch("plane")->SetAddress(&Plane);
  channelParTree->GetBranch("slab")->SetAddress(&Slab);
  channelParTree->GetBranch("pmt")->SetAddress(&Pmt);
  channelParTree->GetBranch("par")->SetAddress(&Par);

  int Par_entries = channelParTree->GetEntries();
  int i = 0;

  while(i<Par_entries)
  {
    channelParTree->GetEntry(i);
    i++;
    cout<<"St"<<Station<<"  Pl"<<Plane<<"  Sl"<<Slab<<"  Pmt"<<Pmt<<"  (";
    for(int p=0;p<(nPar-1);p++)
      cout<<Par[p]<<",";
      cout<<Par[nPar-1]<<")"<<endl;
  }
  cout<<"====================================="<<endl;
}

void MiceTofCalib::SetAddresses(int Station)
{
  if(Station == 0)
  {
    dataTree->GetBranch("slabA")->SetAddress(&_slX);
    dataTree->GetBranch("slabB")->SetAddress(&_slY);
    dataTree->GetBranch("adc0")->SetAddress(&_adc0);
    dataTree->GetBranch("adc1")->SetAddress(&_adc1);
    dataTree->GetBranch("adc2")->SetAddress(&_adc2);
    dataTree->GetBranch("adc3")->SetAddress(&_adc3);
    dataTree->GetBranch("t0")->SetAddress(&_t0);
    dataTree->GetBranch("t1")->SetAddress(&_t1);
    dataTree->GetBranch("t2")->SetAddress(&_t2);
    dataTree->GetBranch("t3")->SetAddress(&_t3);

    dataTree->GetBranch("slabC")->SetAddress(&buffer0);
    dataTree->GetBranch("slabD")->SetAddress(&buffer1);
    dataTree->GetBranch("adc4")->SetAddress(&buffer2);
    dataTree->GetBranch("adc5")->SetAddress(&buffer3);
    dataTree->GetBranch("adc6")->SetAddress(&buffer4);
    dataTree->GetBranch("adc7")->SetAddress(&buffer5);
    dataTree->GetBranch("t4")->SetAddress(&buffer6);
    dataTree->GetBranch("t5")->SetAddress(&buffer7);
    dataTree->GetBranch("t6")->SetAddress(&buffer8);
    dataTree->GetBranch("t7")->SetAddress(&buffer9);

    dataTree->GetBranch("slabE")->SetAddress(&buffer10);
    dataTree->GetBranch("slabF")->SetAddress(&buffer11);
    dataTree->GetBranch("adc8")->SetAddress(&buffer12);
    dataTree->GetBranch("adc9")->SetAddress(&buffer13);
    dataTree->GetBranch("adc10")->SetAddress(&buffer14);
    dataTree->GetBranch("adc11")->SetAddress(&buffer15);
    dataTree->GetBranch("t8")->SetAddress(&buffer16);
    dataTree->GetBranch("t9")->SetAddress(&buffer17);
    dataTree->GetBranch("t10")->SetAddress(&buffer18);
    dataTree->GetBranch("t11")->SetAddress(&buffer19);
  }
  if(Station == 1)
  {
    dataTree->GetBranch("slabC")->SetAddress(&_slX);
    dataTree->GetBranch("slabD")->SetAddress(&_slY);
    dataTree->GetBranch("adc4")->SetAddress(&_adc0);
    dataTree->GetBranch("adc5")->SetAddress(&_adc1);
    dataTree->GetBranch("adc6")->SetAddress(&_adc2);
    dataTree->GetBranch("adc7")->SetAddress(&_adc3);
    dataTree->GetBranch("t4")->SetAddress(&_t0);
    dataTree->GetBranch("t5")->SetAddress(&_t1);
    dataTree->GetBranch("t6")->SetAddress(&_t2);
    dataTree->GetBranch("t7")->SetAddress(&_t3);

    dataTree->GetBranch("slabA")->SetAddress(&buffer0);
    dataTree->GetBranch("slabB")->SetAddress(&buffer1);
    dataTree->GetBranch("adc0")->SetAddress(&buffer2);
    dataTree->GetBranch("adc1")->SetAddress(&buffer3);
    dataTree->GetBranch("adc2")->SetAddress(&buffer4);
    dataTree->GetBranch("adc3")->SetAddress(&buffer5);
    dataTree->GetBranch("t0")->SetAddress(&buffer6);
    dataTree->GetBranch("t1")->SetAddress(&buffer7);
    dataTree->GetBranch("t2")->SetAddress(&buffer8);
    dataTree->GetBranch("t3")->SetAddress(&buffer9);

    dataTree->GetBranch("slabE")->SetAddress(&buffer10);
    dataTree->GetBranch("slabF")->SetAddress(&buffer11);
    dataTree->GetBranch("adc8")->SetAddress(&buffer12);
    dataTree->GetBranch("adc9")->SetAddress(&buffer13);
    dataTree->GetBranch("adc10")->SetAddress(&buffer14);
    dataTree->GetBranch("adc11")->SetAddress(&buffer15);
    dataTree->GetBranch("t8")->SetAddress(&buffer16);
    dataTree->GetBranch("t9")->SetAddress(&buffer17);
    dataTree->GetBranch("t10")->SetAddress(&buffer18);
    dataTree->GetBranch("t11")->SetAddress(&buffer19);
  }

  if(Station == 2)
  {
    dataTree->GetBranch("slabE")->SetAddress(&_slX);
    dataTree->GetBranch("slabF")->SetAddress(&_slY);
    dataTree->GetBranch("adc8")->SetAddress(&_adc0);
    dataTree->GetBranch("adc9")->SetAddress(&_adc1);
    dataTree->GetBranch("adc10")->SetAddress(&_adc2);
    dataTree->GetBranch("adc11")->SetAddress(&_adc3);
    dataTree->GetBranch("t8")->SetAddress(&_t0);
    dataTree->GetBranch("t9")->SetAddress(&_t1);
    dataTree->GetBranch("t10")->SetAddress(&_t2);
    dataTree->GetBranch("t11")->SetAddress(&_t3);

    dataTree->GetBranch("slabA")->SetAddress(&buffer0);
    dataTree->GetBranch("slabB")->SetAddress(&buffer1);
    dataTree->GetBranch("adc0")->SetAddress(&buffer2);
    dataTree->GetBranch("adc1")->SetAddress(&buffer3);
    dataTree->GetBranch("adc2")->SetAddress(&buffer4);
    dataTree->GetBranch("adc3")->SetAddress(&buffer5);
    dataTree->GetBranch("t0")->SetAddress(&buffer6);
    dataTree->GetBranch("t1")->SetAddress(&buffer7);
    dataTree->GetBranch("t2")->SetAddress(&buffer8);
    dataTree->GetBranch("t3")->SetAddress(&buffer9);

    dataTree->GetBranch("slabC")->SetAddress(&buffer10);
    dataTree->GetBranch("slabD")->SetAddress(&buffer11);
    dataTree->GetBranch("adc4")->SetAddress(&buffer12);
    dataTree->GetBranch("adc5")->SetAddress(&buffer13);
    dataTree->GetBranch("adc6")->SetAddress(&buffer14);
    dataTree->GetBranch("adc7")->SetAddress(&buffer15);
    dataTree->GetBranch("t4")->SetAddress(&buffer16);
    dataTree->GetBranch("t5")->SetAddress(&buffer17);
    dataTree->GetBranch("t6")->SetAddress(&buffer18);
    dataTree->GetBranch("t7")->SetAddress(&buffer19);
  }
}
void MiceTofCalib::SaveAllToFile(int nPar)
{
  int Station, Plane, Slab, Pmt;
  double Par[5];

  //channelParTree->Scan();
  channelParTree->GetBranch("station")->SetAddress(&Station);
  channelParTree->GetBranch("plane")->SetAddress(&Plane);
  channelParTree->GetBranch("slab")->SetAddress(&Slab);
  channelParTree->GetBranch("pmt")->SetAddress(&Pmt);
  channelParTree->GetBranch("par")->SetAddress(&Par);

  int Par_entries = channelParTree->GetEntries();
  int i = 0;

  while(i<Par_entries)
  {
    channelParTree->GetEntry(i);
    i++;
    SaveToFile(Station, Plane, Slab, Pmt, Par, nPar);
  }
}

void MiceTofCalib::SaveToFile(int st, int pl, int sl, int pmt, double *par, int npar)
{
  /*
  (*file) << setw(4) <<st;
  (*file) << setw(4) <<pl;
  (*file) << setw(4) <<sl;
  (*file) << setw(4) <<pmt;
  */
  stringstream ss;
  ss << "tof" <<st;
  MAUS::TOFChannelKey key(st, pl, sl, pmt, ss.str());
  (*file) << key << "  ";
  for(int i=0;i<npar;i++)
    (*file) << setw(17) <<par[i];

  (*file)<<std::endl;
}

bool MiceTofCalib::LoadParameters(string fileName, int npar, int station)
{
  ifstream inf( fileName.c_str() );
  if( ! inf )
  {
    std::cerr << "Can't open calibration file " << fileName << std::endl;
    exit(1);
  }
  if( channelParTree->GetEntries() != 0 )
  {
    delete channelParTree;
    MakeParamTree();
  }

  while ( ! inf.eof() )
  {
    /*
    inf >> _pmtPar.st;
    inf >> _pmtPar.plane;
    inf >> _pmtPar.slab;
    inf >> _pmtPar.pmt;
    */
    MAUS::TOFChannelKey key;
    inf >> key;

    _pmtPar.st    = key.station();
    _pmtPar.plane = key.plane();
    _pmtPar.slab  = key.slab();
    _pmtPar.pmt   = key.pmt();
    for(int i=0;i<npar;i++)
      inf >>_pmtPar.par[i];

    if(_pmtPar.par[0] && _pmtPar.par[1] && ( station==-1 || station==_pmtPar.st) )
      channelParTree->Fill();
  }

  return true;
}

void MiceTofCalib::InitDataMonitor()
{
  c1_mon = new TCanvas("DataMonitor","Data Monitor", 700, 700);
  c1_mon->Divide(2,2);
  tof0.SetNameTitle("TOF0","TOF0");
  tof0.SetBins(10,-0.5,9.5,10,-0.5,9.5);
  tof0.GetXaxis()->SetTitle("Slab Num");
  tof0.GetYaxis()->SetTitle("Slab Num");

  tof1.SetNameTitle("TOF1","TOF1");
  tof1.SetBins(7,-0.5,6.5,7,-0.5,6.5);
  tof1.GetXaxis()->SetTitle("Slab Num");
  tof1.GetYaxis()->SetTitle("Slab Num");

  tof2.SetNameTitle("TOF2","TOF2");
  tof2.SetBins(10,-0.5,9.5,10,-0.5,9.5);
  tof2.GetXaxis()->SetTitle("Slab Num");
  tof2.GetYaxis()->SetTitle("Slab Num");

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  double stops[NRGBs] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
  double red[NRGBs]   = { 1.00, 0.00, 0.90, 1.00, 0.80 };
  double green[NRGBs] = { 1.00, 0.80, 1.00, 0.70, 0.00 };
  double blue[NRGBs]  = { 1.00, 1.00, 0.15, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}

void MiceTofCalib::DrawDataMonitor()
{
  c1_mon->cd(1);
  tof0.SetMarkerSize(1.7);
  tof0.Draw();
  tof0.SetDrawOption("text&&colz");
  gPad->Update();
  TPaveStats *st_tof0 = (TPaveStats*)tof0.FindObject("stats");
  st_tof0->SetOptStat(11);
  st_tof0->SetY1NDC(0.89); st_tof0->SetY2NDC(1.);

  c1_mon->cd(2);
  tof1.SetMarkerSize(1.7);
  tof1.Draw();
  tof1.SetDrawOption("text&&colz");
  gPad->Update();
  TPaveStats *st_tof1 = (TPaveStats*)tof1.FindObject("stats");
  st_tof1->SetOptStat(11);
  st_tof1->SetY1NDC(0.89); st_tof1->SetY2NDC(1.);

  c1_mon->cd(3);
  tof2.SetMarkerSize(1.7);
  tof2.Draw();
  tof2.SetDrawOption("text&&colz");
  gPad->Update();
  TPaveStats *st_tof2 = (TPaveStats*)tof2.FindObject("stats");
  st_tof2->SetOptStat(11);
  st_tof2->SetY1NDC(0.89); st_tof2->SetY2NDC(1.);

  c1_mon->Update();
}

