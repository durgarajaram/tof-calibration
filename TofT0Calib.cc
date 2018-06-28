#include "TofT0Calib.hh"
#include "TofTWCalib.hh"
#include "TofTriggerCalib.hh"

TofT0Calib::TofT0Calib(int station, TofTWCalib* tw, TofTriggerCalib* tr)
:TWCalib( tw ), TriggerCalib( tr )
{
  NHistoBins = 200;
  minEntries = 1000;
  nSlabs = (station==0 || station==2)? 10 : 7;
  Station = station;
  stringstream stname;
  stname<<"tof"<<Station;
  c1_t0 = new TCanvas( (stname.str()+"Canvas").c_str(),"T0Canvas",800,600);
}

TofT0Calib::TofT0Calib(int station, TofTWCalib* tw, TofTriggerCalib* tr, TTree *tree)
:TWCalib( tw ), TriggerCalib(tr)
{
  NHistoBins = 1000;
  minEntries = 1000;
  nSlabs = (station==0 || station==2)? 10 : 7;
  Station = station;
  stringstream stname;
  stname<<"tof"<<Station;
  c1_t0 = new TCanvas( (stname.str()+"Canvas").c_str(),"T0Canvas",800,600);
  SetDataTree( tree );
}

TofT0Calib::~TofT0Calib()
{
  hists.erase( hists.begin(), hists.end() );
}

void TofT0Calib::FullCalib(int refSlA, int refSlB, bool verboseMode)
{
  SetAddresses(Station);
  MakeHistograms();
  FillHistograms();

  PlaneCalib(0,refSlB,verboseMode);
  PlaneCalib(1,refSlA,verboseMode);
  //PrintParameters();
}

void TofT0Calib::PlaneCalib(int plane, int refSl, bool verboseMode)
{
  TSpectrum *s1 = new TSpectrum(10);

  double NoCalib[5] = {0., -9999, 0., 0., 0.};
  _pmtPar.st = Station;
  _pmtPar.plane = plane;

  for(int i=0;i<nSlabs;i++)
  {
    _pmtPar.slab = i;
    int Slab_A, Slab_B;
    if(!plane){ Slab_A = i; Slab_B =refSl;}
    else { Slab_A = refSl; Slab_B = i;}

    //cout<<Slab_A<<" "<<Slab_B<<" "<<hists[Slab_A][Slab_B][0]->GetEntries()<<endl;
    TF1 fgaus("fgaus","gaus");
    fgaus.SetLineWidth(1);
    fgaus.SetLineColor(4);
    int k;
    if(!plane) k = 0;
    else k = 2;
    int max = k+2;

    while(k<max)
    {
      if(!plane)
        _pmtPar.pmt = k;
      else _pmtPar.pmt = k-2;
      int histoMeanBin = hists[Slab_A][Slab_B][k]->GetXaxis()->FindBin( hists[Slab_A][Slab_B][k]->GetMean() );
      hists[Slab_A][Slab_B][k]->GetXaxis()->SetRange(histoMeanBin-120,histoMeanBin+120);
      if(hists[Slab_A][Slab_B][k]->GetEntries()>minEntries)
      {
        fgaus.SetParameters(0.,0.,0.);
        //int nfound = s1->Search( hists[Slab_A][Slab_B][k] );
        int nfound = s1->Search( hists[Slab_A][Slab_B][k],2,"goff");
        // search for peaks; if < 3; rebin & search again; if still < 3, skip
        //if(nfound != 3) { 
        if(!nfound) {
          hists[Slab_A][Slab_B][k]->Rebin();
          nfound = s1->Search( hists[Slab_A][Slab_B][k],2,"goff");
          int histoMeanBin = hists[Slab_A][Slab_B][k]->GetXaxis()->FindBin( hists[Slab_A][Slab_B][k]->GetMean() );
          hists[Slab_A][Slab_B][k]->GetXaxis()->SetRange(histoMeanBin-60,histoMeanBin+60);
        }
        //if(nfound!=3) {
        if(!nfound) {
          //cout << "not found: " << k << " " << Slab_A << " " << Slab_B << endl;
          hists[Slab_A][Slab_B][k]->Write();
          SaveToFile(_pmtPar.st, _pmtPar.plane, _pmtPar.slab, _pmtPar.pmt, NoCalib, 2);
          k++;
          continue;
        }

        float *posX = s1->GetPositionX();
        //float *posY = s1->GetPositionY();

        int first = 0, last = 0;
        for(int j=1;j<nfound;j++)
        {
          if(posX[j]<posX[first])
            first = j;
          if(posX[j]>posX[last])
            last = j;
        }
        double t0;
        if( Station > TriggerCalib->GetStation() )
        {
          if(verboseMode)
            cout<<"plane "<<plane<<"   slab "<<i<<"  pmt "<<k<<"   nfound "<<nfound<<" peaks"<<endl;

          hists[Slab_A][Slab_B][k]->Fit("fgaus","Q","",posX[first]-1e3,posX[first]+1e3);
        }else
        {
          if(verboseMode)
            cout<<"plane "<<plane<<"   slab "<<i<<"  pmt "<<k<<"   nfound "<<nfound<<" peaks"<<endl;

          hists[Slab_A][Slab_B][k]->Fit("fgaus","Q","",posX[last]-1e3,posX[last]+1e3);
        }
        gPad->Update();
        hists[Slab_A][Slab_B][k]->Write();
        t0 = fgaus.GetParameter("Mean");
        //cout <<"plane "<<plane<<"   slab "<<i<<"  pmt "<<k<<"  t0 = "<<t0<<" - "<<dT<<endl;
        _pmtPar.par[0] = t0 - dT;
        _pmtPar.par[1] = refSl;
        channelParTree->Fill();
        SaveToFile(_pmtPar.st, _pmtPar.plane, _pmtPar.slab, _pmtPar.pmt, _pmtPar.par, 2);
      }
      else {
        SaveToFile(_pmtPar.st, _pmtPar.plane, _pmtPar.slab, _pmtPar.pmt, NoCalib, 2);
        hists[Slab_A][Slab_B][k]->Write();
      }

      k++;
    }
  }
}

void TofT0Calib::MakeHistograms()
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
    c1_t0->Update();
    bool theEnd = false;
    char yn;
    while ( ! theEnd )
    {
      cout << "Do you want to continue ? [y/n] " ;
      cin >> yn;
      theEnd = ( yn == 'y' ) ? true : false;
    }
    delete c1_t0;
  }

  double mean = range.GetBinCenter( range.GetMaximumBin() );

  hists.resize(nSlabs);
  for(int j=0;j<nSlabs;j++)
  {
    hists[j].resize(nSlabs);
    for(int k=0;k<nSlabs;k++)
    {
      for(int i=0;i<4;i++)
      {
        stringstream name3;
        name3<<"St"<<Station<<"Pixel_h"<<j<<"v"<<k<<"_Pmt"<<i<<"_t0";
        //tr1[j][k].push_back( new TH1F(name3.str().c_str(), name3.str().c_str(),2400,-11e4,-5e4) );
        hists[j][k].push_back( new TH1F(name3.str().c_str(), name3.str().c_str(),NHistoBins,mean-5e4,mean+5e4) );
      }
    }
  }
}

void TofT0Calib::FillHistograms()
{
  cout<<endl<<"Process Station"<<Station<<"  calibration"<<endl;
  cout<<"Distance to the Trigger Station = "<<GetL()<<" cm    Time shift = "<<GetdT()/1e3<<" ns"<<endl;
  int trigger_slX, trigger_slY, triggerStation = TriggerCalib->GetStation();
  if(triggerStation == 0)
  {
    dataTree->GetBranch("slabA")->SetAddress(&trigger_slX);
    dataTree->GetBranch("slabB")->SetAddress(&trigger_slY);
  }
  if(triggerStation == 1)
  {
    dataTree->GetBranch("slabC")->SetAddress(&trigger_slX);
    dataTree->GetBranch("slabD")->SetAddress(&trigger_slY);
  }
  if(triggerStation == 2)
  {
    dataTree->GetBranch("slabE")->SetAddress(&trigger_slX);
    dataTree->GetBranch("slabF")->SetAddress(&trigger_slY);
  }

  int n_entries = dataTree->GetEntries();
  int i = 0;

  //dataTree->Print();
  while(i<n_entries)
  {
    dataTree->GetEntry(i);
    i++;

    if(_slX!=-99 && _slY != -99)
    {
      double tw0, tw1, tw2, tw3;
      if(TWCalib)
      {
        tw0 = TWCalib->TWCorrection(_adc0, Station, 0, _slX, 0);
        tw1 = TWCalib->TWCorrection(_adc1, Station, 0, _slX, 1);
        tw2 = TWCalib->TWCorrection(_adc2, Station, 1, _slY, 0);
        tw3 = TWCalib->TWCorrection(_adc3, Station, 1, _slY, 1);

        if(tw0!=-99999 && tw1!=-99999 &&  tw2!=-99999 && tw3!=-99999  )
        {
          _t0 += tw0;
          _t1 += tw1;
          _t2 += tw2;
          _t3 += tw3;
        }
      }

      double dt = 0.;
      //double dL = GetdL(trigger_slA, trigger_slB, _slA, _slB)*cm;
      //dt = dL/(c_light*picosecond);
      //cout<<"dL = "<< dL/cm <<" dt = "<<dt/picosecond<<"  "<<c_light*picosecond<<endl;
      double t0_tr = TriggerCalib->GetTriggerT0(trigger_slX,trigger_slY);
      if(t0_tr != -99999)
      {
        //if(i<50)cout<<"TDC = "<<tw0<<" "<<tw1<<" "<<tw2<<" "<<tw3<<" "<<dt0<<endl;
        hists[_slX][_slY][0]->Fill(_t0+dt+t0_tr);
        hists[_slX][_slY][1]->Fill(_t1+dt+t0_tr);
        hists[_slX][_slY][2]->Fill(_t2+dt+t0_tr);
        hists[_slX][_slY][3]->Fill(_t3+dt+t0_tr);
      }
    }
  }
}

int TofT0Calib::T0Correction(int station, int plane, int slab, int pmt)
{
  double* T = GetParameters(station, plane, slab, pmt);
  if(!T) return -99999;

  return - int( T[0] );
}

void TofT0Calib::Test()
{
  cout<<"Station"<<Station<<"  Calibration Test "<<endl;
  int trStation = TriggerCalib->GetStation();
  int _slA, _slB, _adc0, _adc1, _adc2, _adc3;
  double _t0, _t1, _t2, _t3, _t4, _t5, _t6, _t7;
  int _slC, _slD, _adc4, _adc5, _adc6, _adc7;
  if(Station==1 && trStation==0)
  {
    dataTree->GetBranch("slabA")->SetAddress(&_slA);
    dataTree->GetBranch("slabB")->SetAddress(&_slB);
    dataTree->GetBranch("adc0")->SetAddress(&_adc0);
    dataTree->GetBranch("adc1")->SetAddress(&_adc1);
    dataTree->GetBranch("adc2")->SetAddress(&_adc2);
    dataTree->GetBranch("adc3")->SetAddress(&_adc3);
    dataTree->GetBranch("t0")->SetAddress(&_t0);
    dataTree->GetBranch("t1")->SetAddress(&_t1);
    dataTree->GetBranch("t2")->SetAddress(&_t2);
    dataTree->GetBranch("t3")->SetAddress(&_t3);

    dataTree->GetBranch("slabC")->SetAddress(&_slC);
    dataTree->GetBranch("slabD")->SetAddress(&_slD);
    dataTree->GetBranch("adc4")->SetAddress(&_adc4);
    dataTree->GetBranch("adc5")->SetAddress(&_adc5);
    dataTree->GetBranch("adc6")->SetAddress(&_adc6);
    dataTree->GetBranch("adc7")->SetAddress(&_adc7);
    dataTree->GetBranch("t4")->SetAddress(&_t4);
    dataTree->GetBranch("t5")->SetAddress(&_t5);
    dataTree->GetBranch("t6")->SetAddress(&_t6);
    dataTree->GetBranch("t7")->SetAddress(&_t7);
  }
  if(Station==0 && trStation==1)
  {
    dataTree->GetBranch("slabA")->SetAddress(&_slC);
    dataTree->GetBranch("slabB")->SetAddress(&_slD);
    dataTree->GetBranch("adc0")->SetAddress(&_adc4);
    dataTree->GetBranch("adc1")->SetAddress(&_adc5);
    dataTree->GetBranch("adc2")->SetAddress(&_adc6);
    dataTree->GetBranch("adc3")->SetAddress(&_adc7);
    dataTree->GetBranch("t0")->SetAddress(&_t4);
    dataTree->GetBranch("t1")->SetAddress(&_t5);
    dataTree->GetBranch("t2")->SetAddress(&_t6);
    dataTree->GetBranch("t3")->SetAddress(&_t7);

    dataTree->GetBranch("slabC")->SetAddress(&_slA);
    dataTree->GetBranch("slabD")->SetAddress(&_slB);
    dataTree->GetBranch("adc4")->SetAddress(&_adc0);
    dataTree->GetBranch("adc5")->SetAddress(&_adc1);
    dataTree->GetBranch("adc6")->SetAddress(&_adc2);
    dataTree->GetBranch("adc7")->SetAddress(&_adc3);
    dataTree->GetBranch("t4")->SetAddress(&_t0);
    dataTree->GetBranch("t5")->SetAddress(&_t1);
    dataTree->GetBranch("t6")->SetAddress(&_t2);
    dataTree->GetBranch("t7")->SetAddress(&_t3);
  }

   if(Station==2 && trStation==0)
  {
    dataTree->GetBranch("slabA")->SetAddress(&_slA);
    dataTree->GetBranch("slabB")->SetAddress(&_slB);
    dataTree->GetBranch("adc0")->SetAddress(&_adc0);
    dataTree->GetBranch("adc1")->SetAddress(&_adc1);
    dataTree->GetBranch("adc2")->SetAddress(&_adc2);
    dataTree->GetBranch("adc3")->SetAddress(&_adc3);
    dataTree->GetBranch("t0")->SetAddress(&_t0);
    dataTree->GetBranch("t1")->SetAddress(&_t1);
    dataTree->GetBranch("t2")->SetAddress(&_t2);
    dataTree->GetBranch("t3")->SetAddress(&_t3);

    dataTree->GetBranch("slabE")->SetAddress(&_slC);
    dataTree->GetBranch("slabF")->SetAddress(&_slD);
    dataTree->GetBranch("adc8")->SetAddress(&_adc4);
    dataTree->GetBranch("adc9")->SetAddress(&_adc5);
    dataTree->GetBranch("adc10")->SetAddress(&_adc6);
    dataTree->GetBranch("adc11")->SetAddress(&_adc7);
    dataTree->GetBranch("t8")->SetAddress(&_t4);
    dataTree->GetBranch("t9")->SetAddress(&_t5);
    dataTree->GetBranch("t10")->SetAddress(&_t6);
    dataTree->GetBranch("t11")->SetAddress(&_t7);
  }
  if(Station==2 && trStation==1)
  {
    dataTree->GetBranch("slabE")->SetAddress(&_slC);
    dataTree->GetBranch("slabF")->SetAddress(&_slD);
    dataTree->GetBranch("adc8")->SetAddress(&_adc4);
    dataTree->GetBranch("adc9")->SetAddress(&_adc5);
    dataTree->GetBranch("adc10")->SetAddress(&_adc6);
    dataTree->GetBranch("adc11")->SetAddress(&_adc7);
    dataTree->GetBranch("t8")->SetAddress(&_t4);
    dataTree->GetBranch("t9")->SetAddress(&_t5);
    dataTree->GetBranch("t10")->SetAddress(&_t6);
    dataTree->GetBranch("t11")->SetAddress(&_t7);

    dataTree->GetBranch("slabC")->SetAddress(&_slA);
    dataTree->GetBranch("slabD")->SetAddress(&_slB);
    dataTree->GetBranch("adc4")->SetAddress(&_adc0);
    dataTree->GetBranch("adc5")->SetAddress(&_adc1);
    dataTree->GetBranch("adc6")->SetAddress(&_adc2);
    dataTree->GetBranch("adc7")->SetAddress(&_adc3);
    dataTree->GetBranch("t4")->SetAddress(&_t0);
    dataTree->GetBranch("t5")->SetAddress(&_t1);
    dataTree->GetBranch("t6")->SetAddress(&_t2);
    dataTree->GetBranch("t7")->SetAddress(&_t3);
  }
   if(Station==0 && trStation==2)
  {
    dataTree->GetBranch("slabE")->SetAddress(&_slA);
    dataTree->GetBranch("slabF")->SetAddress(&_slB);
    dataTree->GetBranch("adc8")->SetAddress(&_adc0);
    dataTree->GetBranch("adc9")->SetAddress(&_adc1);
    dataTree->GetBranch("adc10")->SetAddress(&_adc2);
    dataTree->GetBranch("adc11")->SetAddress(&_adc3);
    dataTree->GetBranch("t8")->SetAddress(&_t0);
    dataTree->GetBranch("t9")->SetAddress(&_t1);
    dataTree->GetBranch("t10")->SetAddress(&_t2);
    dataTree->GetBranch("t11")->SetAddress(&_t3);

    dataTree->GetBranch("slabA")->SetAddress(&_slC);
    dataTree->GetBranch("slabB")->SetAddress(&_slD);
    dataTree->GetBranch("adc0")->SetAddress(&_adc4);
    dataTree->GetBranch("adc1")->SetAddress(&_adc5);
    dataTree->GetBranch("adc2")->SetAddress(&_adc6);
    dataTree->GetBranch("adc3")->SetAddress(&_adc7);
    dataTree->GetBranch("t0")->SetAddress(&_t4);
    dataTree->GetBranch("t1")->SetAddress(&_t5);
    dataTree->GetBranch("t2")->SetAddress(&_t6);
    dataTree->GetBranch("t3")->SetAddress(&_t7);
  }
  if(Station==1 && trStation==2)
  {
    dataTree->GetBranch("slabC")->SetAddress(&_slC);
    dataTree->GetBranch("slabD")->SetAddress(&_slD);
    dataTree->GetBranch("adc4")->SetAddress(&_adc4);
    dataTree->GetBranch("adc5")->SetAddress(&_adc5);
    dataTree->GetBranch("adc6")->SetAddress(&_adc6);
    dataTree->GetBranch("adc7")->SetAddress(&_adc7);
    dataTree->GetBranch("t4")->SetAddress(&_t4);
    dataTree->GetBranch("t5")->SetAddress(&_t5);
    dataTree->GetBranch("t6")->SetAddress(&_t6);
    dataTree->GetBranch("t7")->SetAddress(&_t7);

    dataTree->GetBranch("slabE")->SetAddress(&_slA);
    dataTree->GetBranch("slabF")->SetAddress(&_slB);
    dataTree->GetBranch("adc8")->SetAddress(&_adc0);
    dataTree->GetBranch("adc9")->SetAddress(&_adc1);
    dataTree->GetBranch("adc10")->SetAddress(&_adc2);
    dataTree->GetBranch("adc11")->SetAddress(&_adc3);
    dataTree->GetBranch("t8")->SetAddress(&_t0);
    dataTree->GetBranch("t9")->SetAddress(&_t1);
    dataTree->GetBranch("t10")->SetAddress(&_t2);
    dataTree->GetBranch("t11")->SetAddress(&_t3);
  }

  TH1F h1;
  stringstream name;
  name<<"tof"<<trStation<<"_"<<Station;
  h1.SetNameTitle(name.str().c_str(),"time of flight");
  double e_time = GetdT()/1e3;
  if( e_time>0. ) h1.SetBins(300,e_time-2,e_time+13);
  else h1.SetBins(300,e_time-13,e_time+2);
  stringstream stname;
  stname<<"tof"<<Station;
  //TH1F h1("tof","time of flight",2000,-50.,50.);
  TH1F h2( (stname.str()+"resol").c_str(),"resolution",400,-2e3,2e3);
  //TH1F h2("tof1resol","tof1 resolution",20000,-1e4,1e4);
  //TH1F h2("tof1resol","tof1 resolution",2000,-1e3,1e3);
  TH2I h21( (stname.str()+"profile").c_str(),"tof1 profile",10,-0.5,9.5,10,-0.5,9.5);
  //TH2I h22("tof0profile","tof2 profile",10,-0.5,9.5,10,-0.5,9.5);

  vector< vector<TH1F*> > h_test;
  h_test.resize(nSlabs);
  for(int j=0;j<nSlabs;j++)
  {
    for(int k=0;k<nSlabs;k++)
    {
      stringstream name_test;
      name_test<<"st"<<Station<<"pixel_h"<<j<<"v"<<k<<"Resol";
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
    if(_slC!=-99)
    {
      int t0_0 = TriggerCalib->T0Correction(trStation, 0, _slA, 0, _slB);
      int t0_1 = TriggerCalib->T0Correction(trStation, 0, _slA, 1, _slB);
      int t0_2 = TriggerCalib->T0Correction(trStation, 1, _slB, 0, _slA);
      int t0_3 = TriggerCalib->T0Correction(trStation, 1, _slB, 1, _slA);
      int dt0 = TriggerCalib->GetTriggerT0(_slA,_slB);

      int t0_4 = this->T0Correction(Station, 0, _slC, 0);
      int t0_5 = this->T0Correction(Station, 0, _slC, 1);
      int t0_6 = this->T0Correction(Station, 1, _slD, 0);
      int t0_7 = this->T0Correction(Station, 1, _slD, 1);
      //cout<<"TDC = "<<t0_0<<" "<<t0_1<<" "<<t0_4<<" "<<dt0<<endl;
      if( t0_0!=-99999 && t0_1!=-99999 &&  t0_2!=-99999 && t0_3!=-99999 &&
          t0_4!=-99999 && t0_5!=-99999 &&  t0_6!=-99999 && t0_7!=-99999 && dt0!=-99999 )
      {
        //if(i<1000) cout<<"TDC = "<<_t2<<" "<<t0_2<<endl;
        _t0 += t0_0;
        _t1 += t0_1;
        _t2 += t0_2;
        _t3 += t0_3;

        _t4 += t0_4 + dt0;
        _t5 += t0_5 + dt0;
        _t6 += t0_6 + dt0;
        _t7 += t0_7 + dt0;
        T0correction = true;
      }

      if(TWCalib)
      {
        int tw0 = TWCalib->TWCorrection(_adc0, trStation, 0, _slA, 0);
        int tw1 = TWCalib->TWCorrection(_adc1, trStation, 0, _slA, 1);
        int tw2 = TWCalib->TWCorrection(_adc2, trStation, 1, _slB, 0);
        int tw3 = TWCalib->TWCorrection(_adc3, trStation, 1, _slB, 1);
        int tw4 = TWCalib->TWCorrection(_adc4, Station, 0, _slC, 0);
        int tw5 = TWCalib->TWCorrection(_adc5, Station, 0, _slC, 1);
        int tw6 = TWCalib->TWCorrection(_adc6, Station, 1, _slD, 0);
        int tw7 = TWCalib->TWCorrection(_adc7, Station, 1, _slD, 1);

        if(tw0!=-99999 && tw1!=-99999 &&  tw2!=-99999 && tw3!=-99999 &&
           tw4!=-99999 && tw5!=-99999 &&  tw6!=-99999 && tw7!=-99999  )
        {
          _t0 += tw0;
          _t1 += tw1;
          _t2 += tw2;
          _t3 += tw3;
          _t4 += tw4;
          _t5 += tw5;
          _t6 += tw6;
          _t7 += tw7;
          TWcorrection = true;
        }
      }
      double tof0 = double(_t3 + _t2 + _t1 + _t0)/4.;
      double tof1 = double(_t6 + _t7 + _t4 +_t5)/4.;
      //double dt_tof0 = double(_t3+_t2)/2. - double(_t1+_t0)/2.;
      double dt_tof1 = double(_t7+_t6)/2. - double(_t4+_t5)/2.;
      h21.Fill(_slC,_slD);

      if(T0correction && TWCalib)
        if(TWcorrection)
        {
          h1.Fill( double(tof1-tof0)/1e3);
          h2.Fill(dt_tof1);
          //h20.Fill(_slA,_slB);
          h_test[_slC][_slD]->Fill(dt_tof1);
        }
      if( T0correction && (!TWCalib) )
      {
        h1.Fill( double(tof1-tof0)/1e3);
        h2.Fill(dt_tof1);
        //h20.Fill(_slA,_slB);
        h_test[_slC][_slD]->Fill(dt_tof1);
      }
    }
  }

  //h1.Draw();
  h1.GetXaxis()->SetTitle("time of flight [ns]");
  gPad->Update();
  h1.Write();
  h21.GetXaxis()->SetTitle("slab");
  h21.GetYaxis()->SetTitle("slab");
  h21.SetStats(0);
  gPad->Update();
  h21.Write();
  TF1 fgaus("fgaus","gaus");
  fgaus.SetLineWidth(1);
  fgaus.SetLineColor(4);
  cout<<endl<<"Station "<<Station<<" Total Resolution Fit"<<endl;
  h2.Fit("fgaus");
  h2.GetXaxis()->SetTitle("#Delta t  [ps]");
  gPad->Update();
  h2.Write();
  cout<<endl<<"Pixel Resolutions :"<<endl;
  for(int i=0;i<nSlabs;i++)
  {
    for(int j=0;j<nSlabs;j++)
    {
      if( h_test[i][j]->GetEntries() )
      {
        h_test[i][j]->Fit("fgaus","Q");
        cout<<"   Mean = "<<fgaus.GetParameter("Mean")<<"    Sigma = "<<fgaus.GetParameter("Sigma");
        cout<<"    Entries = "<<h_test[i][j]->GetEntries()<<endl;
        h_test[i][j]->GetXaxis()->SetTitle("#Delta t  [ps]");
        gPad->Update();
        h_test[i][j]->Write();
      }
    }
  }
  h_test.erase( h_test.begin(), h_test.end() );
}

void TofT0Calib::DevTest()
{
  TH1F adc1, adc2, adc3, adc4, adc5;
  TH2F tw1, tw2, tw3;
  TH1F dt("dt","dt",100,-1e3,1e3);
  adc1.SetNameTitle("adc1","adc1");
  adc1.SetBins(120, 200, 6000);

  adc2.SetNameTitle("adc2","adc2");
  adc2.SetBins(120, 200, 6000);

  adc3.SetNameTitle("adc3","adc3");
  adc3.SetBins(120, 200, 6000);

  adc4.SetNameTitle("adc4","adc4");
  adc4.SetBins(120, 200, 6000);

  adc5.SetNameTitle("adc5","adc5");
  adc5.SetBins(120, 200, 6000);

  tw1.SetNameTitle("tw1","tw1");
  tw1.SetBins(120, 200, 6000, 4000,-5000,5000);

  tw2.SetNameTitle("tw2","tw2");
  tw2.SetBins(120, 200, 6000, 4000,-5000,5000);

  tw3.SetNameTitle("tw3","tw3");
  tw3.SetBins(120, 200, 6000, 4000,-5000,5000);

  int trigger_slX, trigger_slY, triggerStation = TriggerCalib->GetStation();
  if(triggerStation == 0)
  {
    dataTree->GetBranch("slabA")->SetAddress(&trigger_slX);
    dataTree->GetBranch("slabB")->SetAddress(&trigger_slY);
  }
  if(triggerStation == 1)
  {
    dataTree->GetBranch("slabC")->SetAddress(&trigger_slX);
    dataTree->GetBranch("slabD")->SetAddress(&trigger_slY);
  }
  if(triggerStation == 2)
  {
    dataTree->GetBranch("slabE")->SetAddress(&trigger_slX);
    dataTree->GetBranch("slabF")->SetAddress(&trigger_slY);
  }

  int n_entries = dataTree->GetEntries();
  int i = 0;

  //dataTree->Print();
  while(i<n_entries)
  {
    dataTree->GetEntry(i);
    i++;

    if(_slX==5)
    {
      double tw_0 = TWCalib->TWCorrection(_adc0, Station, 0, _slX, 0);
      double tw_1 = TWCalib->TWCorrection(_adc1, Station, 0, _slX, 1);
      double tw_2 = TWCalib->TWCorrection(_adc2, Station, 1, _slY, 0);
      double tw_3 = TWCalib->TWCorrection(_adc3, Station, 1, _slY, 1);

      int t0_0 = this->T0Correction(Station, 0, _slX, 0);
      int t0_1 = this->T0Correction(Station, 0, _slX, 1);
      int t0_2 = this->T0Correction(Station, 1, _slY, 0);
      int t0_3 = this->T0Correction(Station, 1, _slY, 1);

      int t0_tr = TriggerCalib->GetTriggerT0(trigger_slX,trigger_slY);

      double t0 = _t0 + t0_0 + tw_0 + t0_tr;
      double t1 = _t1 + t0_1 + tw_1 + t0_tr;
      double t2 = _t2 + t0_2 + tw_2 + t0_tr;
      double t3 = _t3 + t0_3 + tw_3 + t0_tr;
      double deltaT = double(t0+t1)/2. - double(t2+t3)/2.;
      if(_slY==1){
        adc1.Fill(_adc1);
      }
      if(_slY==3) {
        //cout<<deltaT<<endl;
        dt.Fill(deltaT);
        tw1.Fill(_adc1, (t2+t3)/2. - (t1-tw_1) + 1208.49/2.);
        adc2.Fill(_adc1);
      }
      if(_slY==4)
        adc3.Fill(_adc1);
      if(_slY==5){
        tw2.Fill(_adc1, (t2+t3)/2.-(t1-tw_1));
        adc4.Fill(_adc1);
      }
      if(_slY==7){
        tw3.Fill(_adc1, (t2+t3)/2.-(t1-tw_1) - 1153.1/2.);
        adc5.Fill(_adc1);
      }
	/*
      int tw0, tw1, tw2, tw3;
      if(TWCalib)
      {
        tw0 = TWCalib->TWCorrection(_adc0, Station, 0, _slA, 0);
        tw1 = TWCalib->TWCorrection(_adc1, Station, 0, _slA, 1);
        tw2 = TWCalib->TWCorrection(_adc2, Station, 1, _slB, 0);
        tw3 = TWCalib->TWCorrection(_adc3, Station, 1, _slB, 1);

        if(tw0!=-99999 && tw1!=-99999 &&  tw2!=-99999 && tw3!=-99999  )
        {
          _t0 += tw0;
          _t1 += tw1;
          _t2 += tw2;
          _t3 += tw3;
        }
      }

      double dt = 0.;
      //double dL = GetdL(trigger_slA, trigger_slB, _slA, _slB)*cm;
      //dt = dL/(c_light*picosecond);
      //cout<<"dL = "<< dL/cm <<" dt = "<<dt/picosecond<<"  "<<c_light*picosecond<<endl;
      int t0_tr = TriggerCalib->GetTriggerT0(trigger_slA,trigger_slB);
      //cout<<"TDC = "<<_t0<<" "<<_t1<<" "<<_t2<<" "<<_t3<<" "<<t0_tr<<endl;
      if(t0_tr != -99999)
      {

      }*/
    }
  }

  adc1.Write();
  adc2.Write();
  adc3.Write();
  adc4.Write();
  adc5.Write();
  dt.Write();
  tw1.Write();
  tw2.Write();
  tw3.Write();
}

void TofT0Calib::PrintPMTParameters()
{
  cout<<"=========  PMTCable Delay  ==========="<<endl;
  MiceTofCalib::PrintPMTParameters(2);
}
