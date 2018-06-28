{
  TChain chain("dataTree");
  chain.Add("tofcalibdata_3245.root");
  chain.Add("tofcalibdata_3246.root");
  chain.Add("tofcalibdata_3247.root");
  //chain.Add("tofcalibdata_3248.root");
  //chain.Draw("slabB");
  
  TFile file("tofcalibdata_merged.root", "recreate");
  TTree* mergeTree = chain.CloneTree();
  mergeTree->Write();
  file.Close();
}