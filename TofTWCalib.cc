//Y.Karadzhov 02.2009

#include "TofTWCalib.hh"
#include "TofTriggerCalib.hh"
#include "TofT0Calib.hh"

vector<double> coords;
vector<double> values;
vector<double> errors;

double fitf(double *x, double *par)
{
   double x0 = x[0] + par[0];
   double y = par[1] + par[2]/x0 + par[3]/pow(x0,2);
   //double y = par[1] + par[2]/x0 + par[3]/(x0*x0) + par[4]/(x0*x0*x0);
   return y;
}

void myFcn(Int_t &npar, Double_t *grad , Double_t &fval, Double_t *p, Int_t iflag  )
{
  int n = coords.size();
  double chi2 = 0;
  double tmp, x[1];
  for (int i = 0; i <n; ++i ) {
    x[0] = coords[i];
    if(errors[i]) tmp = ( values[i] - fitf(x,p))/errors[i];
    chi2 += tmp*tmp;
  }
  fval = chi2;
  //cout<<"chi2 = "<< chi2<<endl;
}

TofTWCalib::TofTWCalib(string n)
{
  _name = n;
  triggerCalib = NULL;
  t0Calib = NULL;
  nPar = 4;
  adcCut = 0.1;
  f1 = new TF1("f1",fitf,0.,6000.,nPar);
  f1->SetLineWidth(1);
  f1->SetLineColor(2);
  nPar = f1->GetNumberFreeParameters();
  string canvName =  _name + "TWCancas";
  c1_tw = new TCanvas(canvName.c_str(),"TWCancas",1000,900);
  InitializeMinuit();
  MakeHistograms();
}

TofTWCalib::TofTWCalib(TTree *tree, string n)
{
  SetDataTree( tree );
  //dataTree->Print();

  _name = n;
  triggerCalib = NULL;
  t0Calib = NULL;
  nPar = 4;
  adcCut = 0.1;
  string fName =  _name + "f1";
  f1 = new TF1(fName.c_str(),fitf,0.,6000.,nPar);
  f1->SetLineWidth(1);
  f1->SetLineColor(2);
  string canvName =  _name + "TWCancas";
  c1_tw = new TCanvas(canvName.c_str(),"TWCancas",1000,900);
  InitializeMinuit();
  MakeHistograms();
}

TofTWCalib::~TofTWCalib()
{
  ClearHistograms();
  delete f1;
  delete minuit;
}

void TofTWCalib::InitializeMinuit()
{
  minuit = new TMinuit(nPar);
  minuit->SetFCN(myFcn);
  ierflg = 0;

  arglist[0] = -1;
  minuit->mnexcm("SET PRINT",arglist,1,ierflg);
  minuit->mnexcm("SET NOW",arglist,1,ierflg);
}

void TofTWCalib::FullCalib(int st, int refSlA, int refSlB, bool verboseMode, bool Minuit)
{
  cout<<endl<<"Process Station"<<st<<" TW Calibration"<< " " << refSlA << " " << refSlB << endl;
  SetAddresses(st);
  if( this->ReferenceSlabCalib(st, refSlA, refSlB, verboseMode, Minuit) )
  {
    FillHistograms(st, refSlA, refSlB);
    int calibrated = 0;
    int nSlabs = (st==0 || st==2)? 10 : 7;

    for(int i=0;i<nSlabs;i++)
    {
      calibrated += this->SlabCalib(st, 0, i, refSlB, verboseMode);
      calibrated += this->SlabCalib(st, 1, i, refSlA, verboseMode);
    }
      stringstream cut;
      cut<<"station=="<<st;
      cout<<"The number of calibrated slabs is : "<< channelParTree->GetEntries( cut.str().c_str() )/2 <<endl;
  }
}

bool TofTWCalib::ReferenceSlabCalib(int station, int slabX, int slabY, bool verboseMode, bool Minuit)
{
  vector<TH1F*> h1d;
  for(int i=0;i<4;i++)
  {
    stringstream name;
    name<<"st"<<station<<"adcReference"<<i;
    h1d.push_back( new TH1F(name.str().c_str(), name.str().c_str(), 300, 500, 6000) );
  }

  vector<TH2F*> h2d;
  for(int i=0;i<4;i++)
  {
    stringstream name;
    name<<"st"<<station<<"h2dReference"<<i;
    h2d.push_back( new TH2F(name.str().c_str(), name.str().c_str(), 120, 500, 6000, 4000,-50000,50000) );
  }

  vector<TProfile*> hpf;
  for(int i=0;i<4;i++)
  {
    stringstream name;
    name<<"st"<<station<<"hpfReference"<<i;
    hpf.push_back( new TProfile() );
  }

  TH1F tof0("tof0","tof0",400,-50000,50000);

  int n_entries = dataTree->GetEntries();
  int i = 0;
  while(i<n_entries)
  {
    dataTree->GetEntry(i);
    i++;
    //if(station == _st)
    {
      if(slabX == _slX && slabY == _slY)
      {
        h1d[0]->Fill(_adc0);
        h1d[1]->Fill(_adc1);
        h1d[2]->Fill(_adc2);
        h1d[3]->Fill(_adc3);
        tof0.Fill( double(_t3+_t2)/2. - double(_t1+_t0)/2 );
      }
    }
  }
  tof0.Draw();
  
  double adc_mean[4], adc_rms[4];
  adc_mean[0] = h1d[0]->GetBinCenter( h1d[0]->GetMaximumBin() );
  adc_rms[0] = h1d[0]->GetRMS();
  //cout<<"ReferenceSlabCalib "<<h1d[0]->GetTitle()<<" "<<adc_mean[0]<<" +-"<<adc_rms[0]<<endl;

  adc_mean[1] = h1d[1]->GetBinCenter( h1d[1]->GetMaximumBin() );
  adc_rms[1] = h1d[1]->GetRMS();
  //cout<<"ReferenceSlabCalib "<<h1d[1]->GetTitle()<<" "<<adc_mean[1]<<" +-"<<adc_rms[1]<<endl;

  adc_mean[2] = h1d[2]->GetBinCenter( h1d[2]->GetMaximumBin() );
  adc_rms[2] = h1d[2]->GetRMS();
  //cout<<"ReferenceSlabCalib "<<h1d[2]->GetTitle()<<" "<<adc_mean[2]<<" +-"<<adc_rms[2]<<endl;

  adc_mean[3] = h1d[3]->GetBinCenter( h1d[3]->GetMaximumBin() );
  adc_rms[3] = h1d[3]->GetRMS();
  //cout<<"ReferenceSlabCalib "<<h1d[3]->GetTitle()<<" "<<adc_mean[3]<<" +-"<<adc_rms[3]<<endl;

  int RefPMT_A, RefPMT_B;

  if(adc_mean[3]>adc_mean[2]) RefPMT_B = 3;
  else RefPMT_B = 2;

  if(adc_mean[1]>adc_mean[0]) RefPMT_A = 1;
  else RefPMT_A = 0;


  double mean_tof = tof0.GetBinCenter( tof0.GetMaximumBin() );

  TH1F tof1("tof1","tof1", 200, mean_tof-2500, mean_tof+2500);
  TH1F tof2("tof2","tof2", 200, mean_tof-2500, mean_tof+2500);

  i = 0;
  int goodEntries=0;
  while(i<n_entries)
  {
    dataTree->GetEntry(i);
    i++;
    if(slabX == _slX && slabY == _slY)
    {
      if( RefPMT_B == 3 && fabs(_adc3-adc_mean[3])/adc_mean[3] < adcCut )
      {
        h2d[0]->Fill(_adc0, _t3 - _t0);
        h2d[1]->Fill(_adc1, _t3 - _t1);
      }
      if( RefPMT_B == 2 && fabs(_adc2-adc_mean[2])/adc_mean[2] < adcCut )
      {
        h2d[0]->Fill(_adc0, _t2 - _t0);
        h2d[1]->Fill(_adc1, _t2 - _t1);
      }
      if( RefPMT_A == 0 && fabs(_adc0-adc_mean[0])/adc_mean[0] < adcCut )
      {
        h2d[2]->Fill(_adc2, _t0 - _t2);
        h2d[3]->Fill(_adc3, _t0 - _t3);
      }
      if( RefPMT_A == 1 && fabs(_adc1-adc_mean[1])/adc_mean[1] < adcCut )
      {
        h2d[2]->Fill(_adc2, _t1 - _t2);
        h2d[3]->Fill(_adc3, _t1 - _t3);
      }
      goodEntries++;
    }
  }
  if(goodEntries<1000)
  {
    cout<<"Insufficient number of events in the referent pixel "<<goodEntries<<endl;
    return false;
  }

  for(unsigned int i=0;i<4;i++)
  {
    hpf[i] = MakeProfile(h2d[i]);
    TF1 *func;
    if(!Minuit)
      func = RootFit(hpf[i], verboseMode);
    else
      func = MinuitFit(hpf[i], adc_mean[i], adc_rms[i], verboseMode);

    double adc0 = func->GetParameter(1) - func->Eval(adc_mean[i]);
    func->GetParameters(ReferenceSlabPar[i]);
    ReferenceSlabPar[i][1] = adc0;

  }

  c1_tw->Clear();
  
  c1_tw->Divide(2,2);
  for(unsigned int i=0;i<4;i++)
  {
    c1_tw->cd(i+1);
    hpf[i]->SetLineColor(26);
    h2d[i]->Draw();
    hpf[i]->Draw("sames");
  }
  c1_tw->SetName("RefPixFit");
  c1_tw->SetTitle("RefPixFit");
  c1_tw->Update();
  c1_tw->Write();

  i = 0;
  while(i<n_entries)
  {
    dataTree->GetEntry(i);
    i++;
    //if(station == _st)
    {
      if(slabX == _slX && slabY == _slY)
        tof1.Fill( double(_t3+_t2)/2. - double(_t1+_t0)/2. );

      if(slabX == _slX)
      {
        _t0 += TWCorrection(_adc0, ReferenceSlabPar[0]);
        _t1 += TWCorrection(_adc1, ReferenceSlabPar[1]);
      }
      if(slabY == _slY)
      {
        _t2 += TWCorrection(_adc2, ReferenceSlabPar[2]);
        _t3 += TWCorrection(_adc3, ReferenceSlabPar[3]);
      }

      if(slabX == _slX && slabY == _slY)
        tof2.Fill( double(_t3+_t2)/2. - double(_t1+_t0)/2. );
    }
  }

  c1_tw->Clear();
  tof2.SetLineColor(4);
  tof2.Draw();
  tof2.Fit("gaus","Q");

  tof1.SetLineColor(2);
  tof1.Draw("same");
  c1_tw->Update();

  c1_tw->SetName("RefPixResol");
  c1_tw->SetTitle("RefPixResol");
  c1_tw->Write();

  for(unsigned int i=0;i<h1d.size();i++)
  {
    //h1d[i]->Write();
    delete h1d[i];
  }
  h1d.resize(0);

  for(unsigned int i=0;i<h2d.size();i++)
  {
    //h2d[i]->Write();
    delete h2d[i];
  }
  h2d.resize(0);

  for(unsigned int i=0;i<hpf.size();i++)
  {
    //hpf[i]->Write();
    delete hpf[i];
  }
  hpf.resize(0);

  return true;
}

bool TofTWCalib::SlabCalib(int Station, int Plane, int Slab, int refslab, bool verboseMode)
{
  cout<<"Plane "<<Plane<<"  Slab "<<Slab;
  TH2F *Ph2d_0 = GetPMTHistogram(Station, Plane, Slab, 0);
  TH2F *Ph2d_1 = GetPMTHistogram(Station, Plane, Slab, 1);
  TProfile *pf0, *pf1;
  double NoCalib[5] = {0.,0.,0.,0.,0.};
  int events = Ph2d_0->GetEntries();
  int events1 = Ph2d_1->GetEntries();
  if( events < minEntries )
  {
    SaveToFile(Station, Plane, Slab, 0, NoCalib, nPar);
    SaveToFile(Station, Plane, Slab, 1, NoCalib, nPar);
// write even if fail
    stringstream pixelname;
    if(Plane==0) pixelname<<"St"<<Station<<"_Pixel"<<Slab<<refslab<<"r";
    if(Plane==1) pixelname<<"St"<<Station<<"_Pixel"<<refslab<<"r"<<Slab;
    c1_tw->SetName( pixelname.str().c_str() );
    c1_tw->SetTitle( pixelname.str().c_str() );
    c1_tw->Write();
    Ph2d_0->Write();
    Ph2d_1->Write();
    cout<<"   Entries = "<< events << " " << events1 << " -- failed to calibrate" << endl;
    return false;
  }

  pf0 = MakeProfile(Station, Plane, Slab, 0);
  RootFit(Station, Plane, Slab, 0, verboseMode);
  MinuitFit(Station, Plane, Slab, 0, verboseMode);

  pf1 = MakeProfile(Station, Plane, Slab, 1);
  RootFit(Station, Plane, Slab, 1, verboseMode);
  MinuitFit(Station, Plane, Slab, 1, verboseMode);

  CompareFits(Station, Plane, Slab, refslab);
  cout<<"   Entries = "<< events <<" -- calibrated"<<endl;
  TH2F* h1, *h2;
  c1_tw->Clear();
  c1_tw->Divide(2,2);
  c1_tw->cd(1);
  h1 = GetPMTHistogram(Station, Plane, Slab, 0);
  h1->Draw();
  DrawMaxAdc( h1 );
  GetPMTProfile(Station, Plane, Slab, 0)->Draw("same");
  c1_tw->cd(2);
  h2 = GetPMTHistogram(Station, Plane, Slab, 1);
  h2->Draw();
  DrawMaxAdc( h2 );
  GetPMTProfile(Station, Plane, Slab, 1)->Draw("same");
  c1_tw->cd(3);
  GetResolHistogram(Station, Plane, Slab,0)->Draw();
  GetRawTHistogram(Station, Plane, Slab,1)->Draw("sames");
  c1_tw->cd(4);
  GetResolHistogram(Station, Plane, Slab,1)->Draw();
  GetRawTHistogram(Station, Plane, Slab,1)->Draw("sames");

  c1_tw->Update();
  stringstream pixelname;
  if(Plane==0) pixelname<<"St"<<Station<<"_Pixel"<<Slab<<refslab<<"r";
  if(Plane==1) pixelname<<"St"<<Station<<"_Pixel"<<refslab<<"r"<<Slab;
  c1_tw->SetName( pixelname.str().c_str() );
  c1_tw->SetTitle( pixelname.str().c_str() );
  c1_tw->Write();

  Ph2d_0->Write();
  Ph2d_1->Write();
  pf0->Write();
  pf1->Write();

  return true;
}

void TofTWCalib::DrawMaxAdc(TH2F* h2)
{
  double adc_mean = h2->GetBinCenter( h2->ProjectionX()->GetMaximumBin() );
  double max = h2->GetYaxis()->GetXmax();
  double min = h2->GetYaxis()->GetXmin();
  TLine *adc0_L = new TLine(adc_mean,min,adc_mean,max);
  adc0_L->SetLineColor(2);
  adc0_L->SetLineWidth(1);
  adc0_L->Draw("same");
}

void TofTWCalib::ClearHistograms()
{
  hists.erase( hists.begin(), hists.end() );
  profiles.erase( profiles.begin(), profiles.end() );
  raw_t.erase( raw_t.begin(), raw_t.end() );
  t1.erase( t1.begin(), t1.end() );
  t2.erase( t2.begin(), t2.end() );
  t2m.erase( t2m.begin(), t2m.end() );
}

void TofTWCalib::MakeHistograms()
{
  hists.resize(3);
  profiles.resize(3);
  raw_t.resize(3);
  t1.resize(3);
  t2.resize(3);
  t2m.resize(3);
  for(int st=0;st<3;st++)
  {
    hists[st].resize(2);
    profiles[st].resize(2);
    raw_t[st].resize(2);
    t1[st].resize(2);
    t2[st].resize(2);
    t2m[st].resize(2);
    for(int pl=0;pl<2;pl++)
    {
      hists[st][pl].resize(10);
      profiles[st][pl].resize(10);
      t2[st][pl].resize(10);
      t2m[st][pl].resize(10);
      for(int sl=0;sl<10;sl++)
      {
        profiles[st][pl][sl].resize(2);
        stringstream name;
        name<<_name<<"St"<<st<<"Pl"<<pl<<"Sl"<<sl;
        TH1F *tof0 = new TH1F((name.str()+"_rawT").c_str(),(name.str()+"_rawT").c_str(), 4000, -50000, 50000);
        TH1F *tof1 = new TH1F((name.str()+"_cslib1").c_str(),(name.str()+"_calib1").c_str(), 4000, -50000, 50000);

        raw_t[st][pl].push_back(tof0);
        t1[st][pl].push_back(tof1);

        stringstream name0;
        name0<<_name<<"h2d"<<name.str()<<"Pmt0";
        TH2F *Ph2d_0 = new TH2F(name0.str().c_str(), name0.str().c_str(), 120, 500, 6000, 10000,-50000,50000);
        hists[st][pl][sl].push_back( Ph2d_0 );

        stringstream name1;
        name1<<_name<<"h2d"<<name.str()<<"Pmt1";
        TH2F *Ph2d_1 = new TH2F(name1.str().c_str(), name1.str().c_str(), 120, 500, 6000, 10000,-50000,50000);
        hists[st][pl][sl].push_back( Ph2d_1 );
      }
    }
  }
}

void TofTWCalib::ResetHistograms()
{
  ClearHistograms();
  MakeHistograms();
}

void TofTWCalib::FillHistograms(int st, int refSlX, int refSlY)
{
  int n_entries = dataTree->GetEntries();
  int i = 0;
  while(i<n_entries)
  {
    dataTree->GetEntry(i);
    i++;
    if(refSlX==_slX && _slX != -99 && _slY != -99)
    {
        raw_t[st][1][_slY]->Fill( double(_t1+_t0)/2. - double(_t3+_t2)/2. );

        double _T0 = _t0 + TWCorrection(_adc0, ReferenceSlabPar[0]);
        double _T1 = _t1 + TWCorrection(_adc1, ReferenceSlabPar[1]);
        //if( _adc0>200 && _adc1>200 && _adc2>200 && _adc3>200 )
        {
          t1[st][1][_slY]->Fill( double(_T0+_T1)/2. - double(_t3+_t2)/2. );
          hists[st][1][_slY][0]->Fill(_adc2, double(_T0+_T1)/2. - double(_t2));
          hists[st][1][_slY][1]->Fill(_adc3, double(_T0+_T1)/2. - double(_t3));
        }
    }
    if(refSlY==_slY && _slY != -99 && _slX != -99)
    {
        raw_t[st][0][_slX]->Fill( double( _t3+_t2)/2. - double(_t1+_t0)/2. );

        double _T2 = _t2 + TWCorrection(_adc2, ReferenceSlabPar[2]);
        double _T3 = _t3 + TWCorrection(_adc3, ReferenceSlabPar[3]);
        //if( _adc0>200 && _adc1>200 && _adc2>200 && _adc3>200 )
        {
          t1[st][0][_slX]->Fill( double(_T2+_T3)/2. - double(_t1+_t0)/2. );
          hists[st][0][_slX][0]->Fill(_adc0, double(_T2+_T3)/2. - double(_t0));
          hists[st][0][_slX][1]->Fill(_adc1, double(_T2+_T3)/2. - double(_t1));
        }
    }
  }

  cout<<"TofTWCalib::FillHistograms - done"<<endl;
}

void TofTWCalib::MinuitFit(int Station, int Plane, int Slab, int Pmt, bool verboseMode)
{
  TProfile* pf = GetPMTProfile(Station, Plane, Slab, Pmt);
  TH2F *h2 = GetPMTHistogram(Station, Plane, Slab, Pmt);
  //double adc_mean = h2->GetBinCenter( h2->ProjectionX()->GetMaximumBin() );
  double adc_mean = h2->ProjectionX()->GetMean();
  double adc_rms = h2->GetRMS(1);
  TF1* func = MinuitFit(pf, adc_mean, adc_rms, verboseMode);
  double adc0m = func->GetParameter(1) - func->Eval(adc_mean);
  func->GetParameters(PmtPar[2+Pmt]);
  PmtPar[2+Pmt][1] = adc0m;
}

TF1* TofTWCalib::MinuitFit(TProfile* pfh, double adc_mean, double adc_rms, bool verboseMode)
{
  int nBinsX = pfh->GetNbinsX();
  coords.resize(0);
  values.resize(0);
  errors.resize(0);
  for(int bin=1;bin<nBinsX;bin++)
  {
    coords.push_back( pfh->GetBinCenter(bin) );
    values.push_back( pfh->GetBinContent(bin) );
    errors.push_back( pfh->GetBinError(bin) );
    double k = fabs(coords[bin-1] - adc_mean)/(2.*adc_rms);
    if(k<1) errors[bin-1] +=  5*cos(1.55*k)*errors[bin-1];
    //cout<<" adc "<<coords[bin-1]<<"   error "<<errors[bin-1]<<endl;
  }
    // Set starting values and step sizes for parameters
  static double vstart[4] = {0., 0. , 0. , 0.};
  static double step[4] = {0.1 , 0.1 , 0.01 , 0.001};
  minuit->mnparm(0, "a1", vstart[0], step[0], 0,0,ierflg);
  minuit->mnparm(1, "a2", vstart[1], step[1], 0,0,ierflg);
  minuit->mnparm(2, "a3", vstart[2], step[2], 0,0,ierflg);
  minuit->mnparm(3, "a4", vstart[3], step[3], 0,0,ierflg);
	    // minimize
  arglist[0] = 1000; // number of function calls
  arglist[1] = 0.01; // tolerance
  //minuit->mnseek();
  minuit->mnexcm("MIGRAD",arglist,2,ierflg);

  double params[5];
  double par_err[5];
  for (int p=0; p < nPar; ++p)
	 minuit->GetParameter(p, params[p], par_err[p]);
  double chi2, edm, errdef;
  int nvpar, nparx,stat;
  minuit->mnstat(chi2,edm,errdef,nvpar,nparx,stat);
  f1->SetParameters(params);
  f1->SetParErrors(par_err);
  f1->SetChisquare(chi2);
  int ndf = coords.size()-nvpar;
  f1->SetNDF(ndf);
  TF1 *func = new TF1(*f1);
  func->SetLineColor(4);
  pfh->GetListOfFunctions()->Add(func);
  //cout<<"Minuit Fit  chi2/ndf = "<<chi2/ndf<<endl;
  return func;
}
void TofTWCalib::RootFit(int Station, int Plane, int Slab, int Pmt, bool verboseMode)
{
  TProfile* pf = GetPMTProfile(Station, Plane, Slab, Pmt);
  TH2F *h2 = GetPMTHistogram(Station, Plane, Slab, Pmt);
  //double adc_mean = h2->GetBinCenter( h2->ProjectionX()->GetMaximumBin() );
  double adc_mean = h2->ProjectionX()->GetMean();
  TF1* func = RootFit(pf, verboseMode);
  double adc0 = func->GetParameter(1) - func->Eval(adc_mean);
  func->GetParameters(PmtPar[Pmt]);
  PmtPar[Pmt][1] = adc0;
  delete func;
}

TF1* TofTWCalib::RootFit(TProfile* pfh, bool verboseMode)
{
  f1->SetParameters(0.,0.,0.,0.);
  f1->SetLineColor(2);
  if(verboseMode) pfh->Fit("f1","R");
  else  pfh->Fit("f1","RQ");
  TF1 *func = new TF1(*f1);
  return func;
}

void TofTWCalib::CompareFits(int Station, int Plane, int Slab, int refslab)
{
  TH1F *rawT = GetRawTHistogram(Station,Plane,Slab,0);
  TH1F *t1 = GetRawTHistogram(Station,Plane,Slab,1);
  double mean_tof = rawT->GetMean();
  int meanbin = t1->GetXaxis()->FindBin(mean_tof);
  rawT->GetXaxis()->SetRange(meanbin-70, meanbin+70);
  t1->GetXaxis()->SetRange(meanbin-70, meanbin+70);
  rawT->GetXaxis()->SetTitle("#Delta t [ps]");
  t1->GetXaxis()->SetTitle("#Delta t [ps]");
  TF1 fgaus("fgaus","gaus");
  fgaus.SetLineWidth(1);
  t1->Fit("fgaus","Q");

  stringstream pixelname;
  if(Plane==0) pixelname<<"St"<<Station<<"_Pixel"<<Slab<<refslab<<"r";
  if(Plane==1) pixelname<<"St"<<Station<<"_Pixel"<<refslab<<"r"<<Slab;

  TH1F *TOF2 = new TH1F( ("resolNoM"+pixelname.str() ).c_str(),"Pixel Resol Normal Fit", 120, mean_tof-1500, mean_tof+1500);
  TH1F *TOF2m = new TH1F( ("resolM"+pixelname.str() ).c_str(),"Pixel Resol Minuit", 120, mean_tof-1500, mean_tof+1500);

  int n_entries = dataTree->GetEntries();
  int i = 0;
  while(i<n_entries)
  {
    dataTree->GetEntry(i);
    i++;
    if(Plane==0 && Slab==_slX && refslab==_slY)
    {
      _t2 += TWCorrection(_adc2, ReferenceSlabPar[2]);
      _t3 += TWCorrection(_adc3, ReferenceSlabPar[3]);

      double t0 = _t0 + TWCorrection(_adc0, PmtPar[0]);
      double t0m = _t0 + TWCorrection(_adc0, PmtPar[2]);

      double t1 = _t1 + TWCorrection(_adc1, PmtPar[1]);
      double t1m = _t1 + TWCorrection(_adc1, PmtPar[3]);
      if( _adc0>500 && _adc1>500 && _adc2>500 && _adc3>500 )
      {
        TOF2->Fill( double(_t3+_t2)/2. - double(t1+t0)/2. );
        TOF2m->Fill( double(_t3+_t2)/2. - double(t1m+t0m)/2. );
      }
    }

    if(Plane==1 && Slab==_slY && refslab==_slX)
    {
      _t0 += TWCorrection(_adc0, ReferenceSlabPar[0]);
      _t1 += TWCorrection(_adc1, ReferenceSlabPar[1]);

      double t2 = _t2 + TWCorrection(_adc2, PmtPar[0]);
      double t2m = _t2 + TWCorrection(_adc2, PmtPar[2]);

      double t3 = _t3 + TWCorrection(_adc3, PmtPar[1]);
      double t3m = _t3 + TWCorrection(_adc3, PmtPar[3]);
      if( _adc0>500 && _adc1>500 && _adc2>500 && _adc3>500 )
      {
        TOF2->Fill( double(_t1+_t0)/2. - double(t3+t2)/2. );
        TOF2m->Fill( double(_t1+_t0)/2. - double(t3m+t2m)/2. );
      }
    }
  }

  fgaus.SetLineColor(2);
  TOF2->SetLineColor(2);
  TOF2->GetXaxis()->SetTitle("#Delta t [ps]");
  TOF2->Fit("fgaus","Q");
  gPad->Update();
  TPaveStats *st_tof2 = (TPaveStats*)TOF2->FindObject("stats");
  st_tof2->SetTextColor(2);
  st_tof2->SetY1NDC(0.24); st_tof2->SetY2NDC(0.6);

  fgaus.SetLineColor(4);
  TOF2m->SetLineColor(4);
  TOF2m->GetXaxis()->SetTitle("#Delta t [ps]");
  TOF2m->Fit("fgaus","Q");
  gPad->Update();
  TPaveStats *st_tof2m = (TPaveStats*)TOF2m->FindObject("stats");
  st_tof2m->SetTextColor(4);
  st_tof2m->SetY1NDC(0.24); st_tof2m->SetY2NDC(0.6);

  t2[Station][Plane][Slab] = TOF2;
  t2m[Station][Plane][Slab] = TOF2m;

  double sigma, sigma_m;
  sigma = TOF2->GetFunction("fgaus")->GetParameter("Sigma");
  sigma_m = TOF2m->GetFunction("fgaus")->GetParameter("Sigma");

   _pmtPar.st = Station;
  _pmtPar.plane = Plane;
  _pmtPar.slab = Slab;

  cout << "s:m " << sigma << " " << sigma_m << endl;
  if(sigma<sigma_m)
  {
    _pmtPar.pmt = 0;
    for(int j=0;j<5;j++)
      _pmtPar.par[j] = PmtPar[0][j];
    channelParTree->Fill();
    SaveToFile(Station, Plane, Slab, 0, _pmtPar.par, nPar);

    _pmtPar.pmt = 1;
    for(int j=0;j<5;j++)
      _pmtPar.par[j] = PmtPar[1][j];
    channelParTree->Fill();
    SaveToFile(Station, Plane, Slab, 1, _pmtPar.par, nPar);
  } else
  {
    _pmtPar.pmt = 0;
    for(int j=0;j<5;j++)
      _pmtPar.par[j] = PmtPar[2][j];

    channelParTree->Fill();
    SaveToFile(Station, Plane, Slab, 0, _pmtPar.par, nPar);

    _pmtPar.pmt = 1;
    for(int j=0;j<5;j++)
      _pmtPar.par[j] = PmtPar[3][j];

    channelParTree->Fill();
    SaveToFile(Station, Plane, Slab, 1, _pmtPar.par, nPar);
  }
}

TProfile* TofTWCalib::MakeProfile(TH2F* h2)
{
  if(h2->GetEntries()<5e3)
  h2->RebinX();

  double dt_mean = h2->GetMean(2);
  int maxbin = h2->GetYaxis()->FindBin(dt_mean);
  string name( h2->GetTitle() );
  TProfile* hpf = h2->ProfileX( (name+"_pfh").c_str(),maxbin-350,maxbin+350,"");
  h2->GetYaxis()->SetRange(maxbin-350,maxbin+350);
  h2->GetXaxis()->SetTitle("adc");
  h2->GetYaxis()->SetTitle("#Delta t [ps]");
  h2->GetYaxis()->SetTitleOffset(1.3);

  hpf->GetXaxis()->SetTitle("adc");
  hpf->GetYaxis()->SetTitle("#Delta t [ps]");
  hpf->GetYaxis()->SetTitleOffset(1.3);
  hpf->SetLineColor(46);
  hpf->SetLineWidth(2);
  return hpf;
}

TProfile* TofTWCalib::MakeProfile(int Station, int Plane, int Slab, int Pmt)
{
  TH2F *h2 = GetPMTHistogram(Station, Plane, Slab, Pmt);
  TProfile* hpf = MakeProfile(h2);

  profiles[Station][Plane][Slab][Pmt] = hpf;
  return hpf;
}

TH2F* TofTWCalib::GetPMTHistogram(int station, int plane, int slab, int pmt)
{
  return hists[station][plane][slab][pmt];
}

TProfile* TofTWCalib::GetPMTProfile(int station, int plane, int slab, int pmt)
{
  return profiles[station][plane][slab][pmt];
}

TH1F* TofTWCalib::GetResolHistogram(int station, int plane, int slab, bool Minuit)
{
  if(Minuit)
    return t2m[station][plane][slab];

  else
    return t2[station][plane][slab];
}

TH1F* TofTWCalib::GetRawTHistogram(int station, int plane, int slab, bool RefCalib)
{
  if(RefCalib)
    return t1[station][plane][slab];

  else
    return raw_t[station][plane][slab];
}

double TofTWCalib::GetResolution(int st, int slabX, int slabY, bool verboseMode, int& Entris)
{
  if( !(this->GetParameters(st, 0, slabX, 0)) || !(this->GetParameters(st, 1, slabY, 0)) )
    return -99999;

  int Station, Plane, Slab, Pmt;
  double Par[5], SlabPar[4][5];
  for(int i=0;i<4;i++)
    for(int j=0;j<5;j++)
      SlabPar[i][j] = 0.;

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

    if(st==Station && Plane == 0 && slabX==Slab && Pmt==0)
    for(int i=0;i<5;i++)
      SlabPar[0][i] = Par[i];

    if(st==Station && Plane == 0 && slabX==Slab && Pmt==1)
    for(int i=0;i<5;i++)
      SlabPar[1][i] = Par[i];

    if(st==Station && Plane == 1 && slabY==Slab && Pmt==0)
    for(int i=0;i<5;i++)
      SlabPar[2][i] = Par[i];

    if(st==Station && Plane == 1 && slabY==Slab && Pmt==1)
    for(int i=0;i<5;i++)
      SlabPar[3][i] = Par[i];
  }

  if(verboseMode)
  for(int i=0;i<4;i++)
  {
    cout<<"Pmt"<<i<<" --> ";
    for(int j=0;j<4;j++)
      cout<<SlabPar[i][j]<<" ";
    cout<<endl;
  }

  TH1F tof0("tof0","tof0",400,-50000,50000);

  int n_entries = dataTree->GetEntries();
  i = 0;

  while(i<n_entries)
  {
    dataTree->GetEntry(i);
    i++;
    if(slabX==_slX && slabY==_slY)
      tof0.Fill(double(_t3+_t2)/2. - double(_t1+_t0)/2.);
  }

  double mean_tof = tof0.GetBinCenter( tof0.GetMaximumBin() );

  stringstream pixelname;
  pixelname<<"St"<<st<<"Sl"<<slabX<<"Sl"<<slabY;
  TH1F tof1(("resol_raw"+pixelname.str()).c_str(),"resol_raw", 200, mean_tof-2500, mean_tof+2500);
  TH1F tof2(("resol"+pixelname.str()).c_str(),"resol", 200, mean_tof-2500, mean_tof+2500);

  int n = 0;
  i = 0;
  while(i<n_entries)
  {
    dataTree->GetEntry(i);
    i++;
    if(slabX==_slX && slabY==_slY)
    {
      tof1.Fill(double(_t3+_t2)/2. - double(_t1+_t0)/2.);
      //if(i<50) cout<<_t0<<" "<<_t1<<" "<<_t2<<" "<<_t3<<" "<< (double(_t3+_t2)/2. - double(_t1+_t0)/2. )<<endl;
      _t0 += TWCorrection(_adc0, SlabPar[0]);
      _t1 += TWCorrection(_adc1, SlabPar[1]);
      _t2 += TWCorrection(_adc2, SlabPar[2]);
      _t3 += TWCorrection(_adc3, SlabPar[3]);
      n++;
      tof2.Fill(double(_t3+_t2)/2. - double(_t1+_t0)/2.);
    }
  }

  TF1 fgaus("fgaus","gaus");
  double sigma;
  if(verboseMode)
  {
    c1_tw->Clear();
    fgaus.SetLineWidth(2);
    fgaus.SetLineColor(4);
    tof2.Fit("fgaus");
    sigma = fgaus.GetParameter("Sigma");
    tof1.Fit("fgaus","","same");
    c1_tw->SetName( (pixelname.str()+"Resol").c_str() );
    c1_tw->SetTitle( (pixelname.str()+"Resol").c_str() );
    c1_tw->Update();
    c1_tw->Write();
  }
  else
  {
    tof2.Fit("fgaus","Q");
    sigma = fgaus.GetParameter("Sigma");
  }
  tof1.Write();
  tof2.Write();
  Entris = n;
  return sigma;
}

int TofTWCalib::TWCorrection(double adc, double *par)
{
  f1->SetParameters(par);
  int dt = (int)( f1->Eval(adc) );
  return dt;
}

int TofTWCalib::TWCorrection(double adc, int station, int plane, int slab, int pmt)
{
  double *par = GetParameters(station, plane, slab, pmt);
  if(!par) return -99999;

  f1->SetParameters(par);
  int dt = (int)( f1->Eval(adc) );
  return dt;
}

void TofTWCalib::PrintPMTParameters()
{
  cout<<"=========  TW Correction Parameters  ==========="<<endl;
  MiceTofCalib::PrintPMTParameters(4);
}
