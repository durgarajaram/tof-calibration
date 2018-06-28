//Y.Karadzhov 02.2009

#include "TofTriggerCalib.hh"
#include "TofTWCalib.hh"

TofTriggerCalib::TofTriggerCalib(int station, TofTWCalib* tw)
: TWCalib(tw)
{
  MakeTriggerTree();
  c1_tr = new TCanvas("TriggerCanvas","TriggerCanvas",1100,700);

  tr_delays.resize(10);
  for(int sl=0;sl<10;sl++)
    tr_delays[sl].resize(10);

  Station = station;
  nSlabs = (Station==0 || Station==2)? 10 : 7;
  SetAddresses(Station);
  minEntries = 200;

  MarkedSlab.plane = -99;
  MarkedSlab.slab = -99;
  nGoodPixels = 0;
}

TofTriggerCalib::TofTriggerCalib(int station, TofTWCalib* tw, TTree* tree)
: TWCalib(tw)
{
  MakeTriggerTree();
  SetDataTree(tree);
  c1_tr = new TCanvas("TriggerCanvas","TriggerCanvas",1100,700);

  tr_delays.resize(10);
  for(int sl=0;sl<10;sl++)
	 tr_delays[sl].resize(10);

  Station = station;
  nSlabs = (Station==0 || Station==2)? 10 : 7;
  SetAddresses(Station);
  minEntries = 200;

  MarkedSlab.plane = -99;
  MarkedSlab.slab = -99;
  nGoodPixels = 0;
}

TofTriggerCalib::~TofTriggerCalib()
{
  tr_t0.erase( tr_t0.begin(), tr_t0.end() );
  t_LmR.erase( t_LmR.begin(), t_LmR.end() );

  //delete triggerT0Tree;
}

void TofTriggerCalib::FullCalib(int refSlX,int refSlY, bool verboseMode)
{
  MakeHistograms();
  FillHistograms();

  _trT0.slabX = refSlX;
  _trT0.slabY = refSlY;
  _trT0.st = Station;
  _trT0.TrDelay = 0;
// remove duplicate for ref
//  SaveTriggerDelay(_trT0.slabX, _trT0.slabY, _trT0.TrDelay);
//  triggerT0Tree->Fill();
  nGoodPixels ++;

  TriggerCalib(0,refSlX,0,refSlY,verboseMode);
  TriggerCalib(1,refSlY,0,refSlX,verboseMode);

  //TTree* cloneTree = triggerT0Tree->CloneTree();
  for(int i=0;i<nSlabs;i++)
  {
    int t0_trSlX=0, t0_trSlY=0;
    if( tr_delays[i][refSlY] ) t0_trSlX = tr_delays[i][refSlY];
    if( tr_delays[refSlX][i] ) t0_trSlY = tr_delays[refSlX][i];

    if(i!=refSlX && t0_trSlX) TriggerCalib(0,i,t0_trSlX,refSlY,verboseMode);
    if(i!=refSlY && t0_trSlY) TriggerCalib(1,i,t0_trSlY,refSlX,verboseMode);

  }

  for(int i=0;i<nSlabs;i++)
  for(int j=0;j<nSlabs;j++)
  {
    if( tr_delays[i][j] )
    {
      _trT0.st = Station;
      _trT0.slabX = i;
      _trT0.slabY = j;
      _trT0.TrDelay = int( tr_delays[i][j] );
      SaveTriggerDelay(_trT0.slabX, _trT0.slabY, _trT0.TrDelay);
      triggerT0Tree->Fill();
      nGoodPixels ++;
    }
    else
    {
      //SaveTriggerDelay(i,j,int(tr_delays[i][j]));
    }


    c1_tr->Clear();
    c1_tr->Divide(2,2);
    c1_tr->cd(1);
    GetPMTHistogram(0,i,j,0)->Draw();
    c1_tr->cd(2);
    GetPMTHistogram(0,i,j,1)->Draw();
    c1_tr->cd(3);
    GetPMTHistogram(1,j,i,0)->Draw();
    c1_tr->cd(4);
    GetPMTHistogram(1,j,i,1)->Draw();

    stringstream pixelname;
    pixelname<<"Pixel_h"<<i<<"v"<<j;
    c1_tr->SetName( pixelname.str().c_str() );
    c1_tr->SetTitle( pixelname.str().c_str() );
    c1_tr->Update();
    c1_tr->Write();
  }

  cout<< "Number of calibrated pixels is : "<< nGoodPixels <<endl;
  //Test();
}

bool TofTriggerCalib::TriggerCalib(int Plane, int Slab, double refTime, int refSl, bool verboseMode)
{
  TF1 fgaus("fgaus","gaus");
  fgaus.SetLineWidth(1);
  fgaus.SetLineColor(4);

  double NoCalib[5] = {0.,-99999,0.,0.,0.};
  double pixel0_Sl_mean_0, pixel0_Sl_mean_1;

  TH1F* ref_h0 = this->GetPMTHistogram(Plane,Slab,refSl,0);
  if(ref_h0->GetEntries()>minEntries)
  {
    SetHistoRange(ref_h0);
    fgaus.SetParameters(0.,0.,0.);
    ref_h0->Fit("fgaus","Q");
    pixel0_Sl_mean_0 = fgaus.GetParameter("Mean");
    //pixel0_Sl_mean_0 = ref_h0->GetMean();

    TH1F* ref_h1 = this->GetPMTHistogram(Plane,Slab,refSl,1);
    SetHistoRange(ref_h1);
    fgaus.SetParameters(0.,0.,0.);
    ref_h1->Fit("fgaus","Q");
    pixel0_Sl_mean_1 = fgaus.GetParameter("Mean");
    //pixel0_Sl_mean_1 = ref_h1->GetMean();
  }

  for(int i=0;i<nSlabs;i++)
  {
    double thispixel_Sl_mean_0=0, thispixel_Sl_mean_1=0;
    TH1F* pixel_h0 = this->GetPMTHistogram(Plane,Slab,i,0);
    if(pixel_h0->GetEntries()>minEntries)
    {
      SetHistoRange(pixel_h0);
      fgaus.SetParameters(0.,0.,0.);
      pixel_h0->Fit("fgaus","Q");
      thispixel_Sl_mean_0 = fgaus.GetParameter("Mean");
      //thispixel_Sl_mean_0 = pixel_h0->GetMean();

      TH1F* pixel_h1 = this->GetPMTHistogram(Plane,Slab,i,1);
      SetHistoRange(pixel_h1);
      fgaus.SetParameters(0.,0.,0.);
      pixel_h1->Fit("fgaus","Q");
      thispixel_Sl_mean_1 = fgaus.GetParameter("Mean");
      //thispixel_Sl_mean_1 = pixel_h1->GetMean();
    }

    double Tr_delay=0.;

    _pmtPar.st = Station;
    _pmtPar.plane = !Plane;
    _pmtPar.slab = i;

    if( (!thispixel_Sl_mean_0 || !thispixel_Sl_mean_1) && !refTime)
    {
      SaveToFile(_pmtPar.st, _pmtPar.plane, _pmtPar.slab, 0, NoCalib, 2);
      SaveToFile(_pmtPar.st, _pmtPar.plane, _pmtPar.slab, 1, NoCalib, 2);
        _pmtPar.pmt = 0;
        _pmtPar.par[0] = 0;
        _pmtPar.par[1] = NoCalib[1];
        channelParTree->Fill();
        _pmtPar.pmt = 1;
        channelParTree->Fill();
    }

    if(thispixel_Sl_mean_0 && thispixel_Sl_mean_1)
    {
      Tr_delay = ( (pixel0_Sl_mean_0 - thispixel_Sl_mean_0) + (pixel0_Sl_mean_1 - thispixel_Sl_mean_1) )/2.;

      if(Plane==0)
      {
        if( tr_delays[Slab][i] )
        {
          double t1 = tr_delays[Slab][i];
          double t2 = Tr_delay + refTime;
          double trD = (t1 + t2)/2.;
          if(fabs(t1-t2)>30. && verboseMode)
          {
            cout<<"WARNING in pixel h"<<Slab<<"v"<<i<<" ";
            cout<<"tr_delay1 = "<<t1<<"  tr_delay2 = "<<t2<<"   delta = "<<t1-t2<<endl;
          }
          tr_delays[Slab][i] = trD;
        }
        else tr_delays[Slab][i] = Tr_delay + refTime;
      }
      else
      {
        if(tr_delays[i][Slab])
        {
          double t1 = tr_delays[i][Slab];
          double t2 = Tr_delay + refTime;
          double trD = (t1 + t2)/2.;
          if(fabs(t1-t2)>30. && verboseMode)
          {
            cout<<"WARNING in pixel h"<<i<<"v"<<Slab<<" ";
            cout<<"tr_delay1 = "<<t1<<"  tr_delay2 = "<<t2<<"   delta = "<<t1-t2<<endl;
          }
          tr_delays[i][Slab] = trD;
        }
        else tr_delays[i][Slab] = Tr_delay + refTime;
      }

      if(!refTime)
      {
        TH1F* pixel_h2 = this->GetPMTHistogram(!Plane,i,Slab,0);
        SetHistoRange(pixel_h2);
        pixel_h2->SetFillColor(25);
        fgaus.SetParameters(0.,0.,0.);
        pixel_h2->Fit("fgaus","Q");
        _pmtPar.pmt = 0;
        _pmtPar.par[0] = fgaus.GetParameter("Mean");
        //_pmtPar.par[0] = pixel_h2->GetMean();
        _pmtPar.par[1] = Slab;
        channelParTree->Fill();
        SaveToFile(_pmtPar.st, _pmtPar.plane, _pmtPar.slab, _pmtPar.pmt, _pmtPar.par, 2);

        TH1F* pixel_h3 = this->GetPMTHistogram(!Plane,i,Slab,1);
        SetHistoRange(pixel_h3);
        pixel_h3->SetFillColor(25);
        fgaus.SetParameters(0.,0.,0.);
        pixel_h3->Fit("fgaus","Q");
        _pmtPar.pmt = 1;
        _pmtPar.par[0] = fgaus.GetParameter("Mean");
        //_pmtPar.par[0] = pixel_h3->GetMean();
        _pmtPar.par[1] = Slab;
        channelParTree->Fill();
        SaveToFile(_pmtPar.st, _pmtPar.plane, _pmtPar.slab, _pmtPar.pmt, _pmtPar.par, 2);
      }
    }
  }
  return true;
}

void TofTriggerCalib::MakeTriggerTree()
{
  triggerT0Tree = new TTree("triggerT0Tree","TriggerDelay");
  triggerT0Tree->Branch("station", &_trT0.st, "station/I");
  triggerT0Tree->Branch("slabX", &_trT0.slabX, "slabX/I");
  triggerT0Tree->Branch("slabY", &_trT0.slabY, "slabY/I");
  triggerT0Tree->Branch("TrDelay", &_trT0.TrDelay, "TrDelay/I");
}

void TofTriggerCalib::MakeHistograms()
{
  TH1F range("rangeTest","rangeTest",500,-2.5e5,2.5e5);

  int n_entries = dataTree->GetEntries();
  int i = 0;

  while(i<n_entries)
  {
    dataTree->GetEntry(i);
    i++;
    if(_t0!=0)
    {
      range.Fill(_t0);
      range.Fill(_t1);
      range.Fill(_t2);
      range.Fill(_t3);
    }
  }

  if(0)
  {
    range.Draw();
    c1_tr->Update();
    bool theEnd = false;
    char yn;
    while ( ! theEnd )
    {
      cout << "Do you want to continue ? [y/n] " ;
      cin >> yn;
      theEnd = ( yn == 'y' ) ? true : false;
    }
  }

  double mean = range.GetBinCenter( range.GetMaximumBin() );
  //cout<<"mean = "<<mean<<"   historange ("<<mean-2e4<<","<<mean+2e4<<")"<<endl;
  tr_t0.resize(nSlabs);
  t_LmR.resize(nSlabs);

  for(int j=0;j<nSlabs;j++)
  {
    tr_t0[j].resize(nSlabs);
    t_LmR[j].resize(nSlabs);

    for(int k=0;k<nSlabs;k++)
    {
      stringstream name0,name1;
      name0<<"LmR_Pl0_Sl"<<j<<"_refPl1_Slab"<<k;
      name1<<"LmR_Pl1_Sl"<<k<<"_refPl0_Slab"<<j;

      t_LmR[j][k].push_back( new TH1F(name0.str().c_str(), name0.str().c_str(),4000,-5e4,5e4) );
      t_LmR[j][k].push_back( new TH1F(name1.str().c_str(), name1.str().c_str(),4000,-5e4,5e4) );
      for(int i=0;i<4;i++)
      {
        stringstream name3;
        name3<<"t0_Sl"<<j<<"_Sl"<<k<<"_Pmt"<<i;
        tr_t0[j][k].push_back( new TH1F(name3.str().c_str(), name3.str().c_str(),4000,mean-5e4,mean+5e4) );
      }
    }
  }
}

void TofTriggerCalib::FillHistograms()
{
  int n_entries = dataTree->GetEntries();
  int i = 0;
  cout<<endl<<"Process Trigger Calibration"<<endl;
  while(i<n_entries)
  {
    dataTree->GetEntry(i);
    i++;
    if(_slX!=-99)
    {
      if(TWCalib)
      {
        int tw0 = TWCalib->TWCorrection(_adc0, Station, 0, _slX, 0);
        int tw1 = TWCalib->TWCorrection(_adc1, Station, 0, _slX, 1);
        int tw2 = TWCalib->TWCorrection(_adc2, Station, 1, _slY, 0);
        int tw3 = TWCalib->TWCorrection(_adc3, Station, 1, _slY, 1);

        //cout<<"TDC = "<<_t0<<" "<<_t1<<" "<<_t2<<" "<<_t3<<endl;
        if(tw0!=-99999 && tw1!=-99999 &&  tw2!=-99999 && tw3!=-99999  )
        {
          _t0 += tw0;
          _t1 += tw1;
          _t2 += tw2;
          _t3 += tw3;
        } else continue;
      }

      tr_t0[_slX][_slY][0]->Fill(_t0);
      tr_t0[_slX][_slY][1]->Fill(_t1);
      tr_t0[_slX][_slY][2]->Fill(_t2);
      tr_t0[_slX][_slY][3]->Fill(_t3);
      t_LmR[_slX][_slY][0]->Fill(_t0-_t1);
      t_LmR[_slX][_slY][1]->Fill(_t2-_t3);
    }
  }
}

TH1F* TofTriggerCalib::GetPMTHistogram(int plane, int slab, int refslab, int pmt)
{
  TH1F *h = NULL;
  if(plane==0)
    h = tr_t0[slab][refslab][pmt];
  if(plane==1)
    h = tr_t0[refslab][slab][pmt+2];

  if(plane == MarkedSlab.plane && slab == MarkedSlab.slab)
    h->SetLineColor(2);
  return h;
}

TH1F* TofTriggerCalib::GetDeltaTHistogram(int plane, int slab, int refslab)
{
  if(plane==0)
    return t_LmR[slab][refslab][0];
  if(plane==1)
    return t_LmR[refslab][slab][1];

  return NULL;
}

void TofTriggerCalib::SetHistoRange(TH1F *hist)
{
  int histoMaxBib = hist->GetMaximumBin();
  hist->GetXaxis()->SetRange(histoMaxBib-100,histoMaxBib+100);
  const char* title = hist->GetTitle();
  if(title[0] == 't')
    hist->GetXaxis()->SetTitle("t(pmt) - t(trigger)  [ps]");
  if(title[0] == 'L')
    hist->GetXaxis()->SetTitle("t(pmt0) - t(pmt1)  [ps]");
}

int TofTriggerCalib::GetTriggerT0(int slX, int slY)
{
  int _slabX, _slabY, _trDelay, t0=-99999;

  triggerT0Tree->GetBranch("slabX")->SetAddress(&_slabX);
  triggerT0Tree->GetBranch("slabY")->SetAddress(&_slabY);
  triggerT0Tree->GetBranch("TrDelay")->SetAddress(&_trDelay);

  int entries = triggerT0Tree->GetEntries();
  int i = 0;

  while(i<entries)
  {
    triggerT0Tree->GetEntry(i);
    i++;
    if(slX==_slabX && slY==_slabY)
      t0 = _trDelay;
  }
  //cout << "GetTriggerT0 : slX " << slX << "  slY " << slY << "   t0 "<< t0 << endl;
  return t0;
}

int TofTriggerCalib::GetTriggerT0(TTree *tree, int slX, int slY)
{
  int _slabX, _slabY, _trDelay, t0=-99999;

  tree->GetBranch("slabX")->SetAddress(&_slabX);
  tree->GetBranch("slabY")->SetAddress(&_slabY);
  tree->GetBranch("TrDelay")->SetAddress(&_trDelay);

  int entries = tree->GetEntries();
  int i = 0;

  while(i<entries)
  {
    tree->GetEntry(i);
    i++;
    if(slX==_slabX && slY==_slabY)
      t0 = _trDelay;
  }

  return t0;
}

int TofTriggerCalib::GetChannelT0(int station, int plane, int slab, int pmt, int &refslab)
{
  double* T = GetParameters(station, plane, slab, pmt);
  if(!T)    return 0;
  refslab = int(T[1]);
  return int( T[0] );
}

int TofTriggerCalib::T0Correction(int station, int plane, int slab, int pmt, int crossSlab)
{
  int refslab=-99,t0_pmt=0, t0_tr=0, t0_tr_ref=0;
  t0_pmt = GetChannelT0(station, plane, slab, pmt, refslab);
  if(!plane)
  {
    t0_tr = GetTriggerT0(slab,crossSlab);
    t0_tr_ref = GetTriggerT0(slab,refslab);
  }
  else
  {
    t0_tr = GetTriggerT0(crossSlab,slab);
    t0_tr_ref = GetTriggerT0(refslab,slab);
  }

  if(!t0_pmt || t0_tr==-99999)
    return -99999;

  return - t0_pmt - t0_tr_ref + t0_tr;
}

void TofTriggerCalib::Test()
{
  cout<<"Trigger Calibration Test"<<endl;

  TH1F h1("tof_trigger_resol","time resolution",400,-2e3,2e3);
  //TH1F h1("tof0resol","tof0 resolution",2000,-1e3,1e3);
  //TH1F h1("tof0resol","tof0 resolution",800,-8,8);
  //TH1F h2("tof0resol_raw","tof0 resolution",800,-8,8);
  TH2I h20("tof_trigger_profile","profile",10,-0.5,9.5,10,-0.5,9.5);

  vector< vector<TH1F*> > h_test;
  h_test.resize(nSlabs);
  for(int j=0;j<nSlabs;j++)
  {
    for(int k=0;k<nSlabs;k++)
    {
      stringstream name_test;
      name_test<<"pix_h"<<j<<"v"<<k<<"res";
      h_test[j].push_back( new TH1F(name_test.str().c_str(),name_test.str().c_str(),400,-2e3,2e3));
    }
  }

  int n_entries = dataTree->GetEntries();
  int i = 0;
  while(i<n_entries)
  {
    dataTree->GetEntry(i);
    i++;
    bool T0correction=false, TWcorrection=false;
    if(_slX!=-99)
    {
		//double dt_tof0_raw = double(_t3+_t2)/2. - double(_t1+_t0)/2.;

      int t0_0 = T0Correction(Station, 0, _slX, 0, _slY);
      int t0_1 = T0Correction(Station, 0, _slX, 1, _slY);
      int t0_2 = T0Correction(Station, 1, _slY, 0, _slX);
      int t0_3 = T0Correction(Station, 1, _slY, 1, _slX);

        //cout<<"TDC = "<<t0_0<<" "<<t0_1<<" "<<t0_2<<" "<<t0_3<<endl;
      if(t0_0!=-99999 && t0_1!=-99999 &&  t0_2!=-99999 && t0_3!=-99999)
      {
        //if(i<100) cout<<"TDC = "<<_t2<<" "<<t0_2<<" "<<tw2<<endl;
        _t0 += t0_0;
        _t1 += t0_1;
        _t2 += t0_2;
        _t3 += t0_3;
        T0correction = true;
      }

      if(TWCalib)
      {
        int tw0 = TWCalib->TWCorrection(_adc0, Station, 0, _slX, 0);
        int tw1 = TWCalib->TWCorrection(_adc1, Station, 0, _slX, 1);
        int tw2 = TWCalib->TWCorrection(_adc2, Station, 1, _slY, 0);
        int tw3 = TWCalib->TWCorrection(_adc3, Station, 1, _slY, 1);

        //cout<<"TDC = "<<t0_0<<" "<<t0_1<<" "<<t0_2<<" "<<t0_3<<endl;
        if(tw0!=-99999 && tw1!=-99999 &&  tw2!=-99999 && tw3!=-99999  )
        {
          _t0 += tw0;
          _t1 += tw1;
          _t2 += tw2;
          _t3 += tw3;
          TWcorrection = true;
        }
      }
      double dt_tof0 = double(_t3+_t2)/2. - double(_t1+_t0)/2.;
      //double tof = double(_t3+_t2)/2. + double(_t1+_t0)/2.;
      h20.Fill(_slY,_slX);
      if(T0correction && TWCalib)
        if(TWcorrection)
        {
          h1.Fill(dt_tof0);
          //h2.Fill(dt_tof0_raw/1e3);
          h_test[_slX][_slY]->Fill(dt_tof0);
        }
      if( T0correction && (!TWCalib) )
      {
        h1.Fill(dt_tof0);
        h_test[_slX][_slY]->Fill(dt_tof0);
      }
    }
  }

  TF1 fgaus("fgaus","gaus",-400,400);
  fgaus.SetLineWidth(1);
  fgaus.SetLineColor(4);
  cout<<endl<<"Station "<<Station<<" Total Resolution Fit"<<endl;
  h1.Fit("fgaus","R");
  h1.GetXaxis()->SetTitle("#Delta t [ps]");
  h1.GetYaxis()->SetTitle("events");
  gPad->Update();
  h20.GetXaxis()->SetTitle("slab");
  h20.GetYaxis()->SetTitle("slab");
  h20.SetStats(0);

  cout<<endl<<"Pixel Resolutions :"<<endl;
  for(int i=0;i<nSlabs;i++)
  {
    for(int j=0;j<nSlabs;j++)
    {
      if( h_test[i][j]->GetEntries() )
      {
        h_test[i][j]->Fit("fgaus","Q");
        cout<<"Pixel h"<<i<<"v"<<j;
        cout<<"   Mean = "<<fgaus.GetParameter("Mean")<<"    Sigma = "<<fgaus.GetParameter("Sigma");
        cout<<"    Entries = "<<h_test[i][j]->GetEntries()<<endl;
        gPad->Update();
        h_test[i][j]->Write();
      }
    }
  }
  h_test.erase( h_test.begin(), h_test.end() );
  h1.Write();
  //h2.Write();
  h20.Write();
}

void TofTriggerCalib::PrintTriggerDelays()
{
  int St, SlabX, SlabY, dt;

  //channelParTree->Scan();
  triggerT0Tree->GetBranch("station")->SetAddress(&St);
  triggerT0Tree->GetBranch("slabX")->SetAddress(&SlabX);
  triggerT0Tree->GetBranch("slabY")->SetAddress(&SlabY);
  triggerT0Tree->GetBranch("TrDelay")->SetAddress(&dt);

  int entries = triggerT0Tree->GetEntries();
  int i = 0;
  cout<<"=========  Trigger Delay  ==========="<<endl;

   while(i<entries)
  {
    triggerT0Tree->GetEntry(i);
    i++;

    cout<<"Station "<<St<<"  Pixel "<<SlabX<<SlabY<<"    Trigger Delay "<<dt<<endl;
  }
  cout<<"======================================"<<endl;
}

void TofTriggerCalib::PrintPMTParameters()
{
  cout<<"=========  PMTCable Delay  ==========="<<endl;
  MiceTofCalib::PrintPMTParameters(2);
}

void TofTriggerCalib::SaveTriggerDelay(int slX, int slY, int tr_t0)
{
  //(*Triggerfile) << setw(4) <<Station;
  //(*Triggerfile) << setw(4) <<slA;
  //(*Triggerfile) << setw(4) <<slB;
  stringstream ss;
  ss << "tof" << Station;
  MAUS::TOFPixelKey key(Station, slX, slY, ss.str());
  (*Triggerfile) << key;
  (*Triggerfile) << setw(5) <<tr_t0;
  (*Triggerfile) <<endl;
}

bool TofTriggerCalib::LoadTriggerDelays(string fileName)
{
  ifstream inf( fileName.c_str() );
  if( ! inf )
  {
    std::cerr << "Can't open calibration file " << fileName << std::endl;
    exit(1);
  }

  if( triggerT0Tree->GetEntries() != 0 )
  {
    delete triggerT0Tree;
    MakeTriggerTree();
  }

  while ( ! inf.eof() )
  {
    //inf >> _trT0.st;
    //inf >> _trT0.slabA;
    //inf >> _trT0.slabB;
    
    MAUS::TOFPixelKey key;
    inf >> key >> _trT0.TrDelay;
    _trT0.st = key.station();
    _trT0.slabX = key.slabX();
    _trT0.slabY = key.slabY();
//     if( _trT0.st || _trT0.slabX || _trT0.slabY || _trT0.TrDelay )
      triggerT0Tree->Fill();
  }
  return true;
}
