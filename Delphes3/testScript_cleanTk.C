/*
root -l examples/myScript.C\(\"QCDoutput5.root\"\)

for clean tracks ie efficiency = 1, no smearing
*/
using namespace std;
// void myScript(const char *inputFile="")
void testScript_cleanTk()
{
  TH1::SetDefaultSumw2();

  gSystem->Load("libDelphes");

  bool doSignal = false;
  bool doMu = true;

  // Create chain of root trees
  TChain chain("Delphes");
  if (doSignal){
    // chain.Add("GG_H_aa.root");
    // chain.Add("Signal_cleanTk/signal_clean.root");
    chain.Add("Signal_1prong_cleanTk/signal_1prong_cleanTk.root");
    // chain.Add("Signal_3prong_cleanTk/signal_3prong_cleanTk.root");
    cout << "Doing signal" << endl;
  } else {
    if (doMu){
      cout << "Doing QCDb_mu" << endl;
      chain.Add("QCDb_mu_cleanTk/QCDb_mu_1.root");
      chain.Add("QCDb_mu_cleanTk/QCDb_mu_10.root");
      chain.Add("QCDb_mu_cleanTk/QCDb_mu_11.root");
      chain.Add("QCDb_mu_cleanTk/QCDb_mu_12.root");
      chain.Add("QCDb_mu_cleanTk/QCDb_mu_13.root");
      chain.Add("QCDb_mu_cleanTk/QCDb_mu_14.root");
      chain.Add("QCDb_mu_cleanTk/QCDb_mu_15.root");
      chain.Add("QCDb_mu_cleanTk/QCDb_mu_16.root");
      chain.Add("QCDb_mu_cleanTk/QCDb_mu_17.root");
      chain.Add("QCDb_mu_cleanTk/QCDb_mu_18.root");
      chain.Add("QCDb_mu_cleanTk/QCDb_mu_19.root");
      chain.Add("QCDb_mu_cleanTk/QCDb_mu_2.root");
      chain.Add("QCDb_mu_cleanTk/QCDb_mu_20.root");
      chain.Add("QCDb_mu_cleanTk/QCDb_mu_3.root");
      chain.Add("QCDb_mu_cleanTk/QCDb_mu_4.root");
      chain.Add("QCDb_mu_cleanTk/QCDb_mu_5.root");
      chain.Add("QCDb_mu_cleanTk/QCDb_mu_6.root");
    } else {
      cout << "Doing QCDb" << endl;
      chain.Add("QCDb_cleanTk/QCDb_10.root");
      chain.Add("QCDb_cleanTk/QCDb_2.root");
      chain.Add("QCDb_cleanTk/QCDb_3.root");
      chain.Add("QCDb_cleanTk/QCDb_4.root");
      chain.Add("QCDb_cleanTk/QCDb_5.root");
      chain.Add("QCDb_cleanTk/QCDb_6.root");
      chain.Add("QCDb_cleanTk/QCDb_7.root");
      chain.Add("QCDb_cleanTk/QCDb_8.root");
      chain.Add("QCDb_cleanTk/QCDb_9.root");
    }
  }

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  // Get pointers to branches used in this analysis
  // Use the data_flow.png and tcl file to figure out what branches are available, and what class they are
  // and use https://cp3.irmp.ucl.ac.be/projects/delphes/wiki/WorkBook/RootTreeDescription
  // TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchTracks = treeReader->UseBranch("Track");
  TClonesArray *branchGenMuons = treeReader->UseBranch("OnlyGenMuons");
  TClonesArray *branchStable = treeReader->UseBranch("StableParticle");

  // Book histograms
  TH1D *histNTracks1OS = new TH1D("hNTracks1OS" ,"Number of tracks about mu1, OS, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); N_{trk} about muon1 / N (muon1)", 10,0,1.0);
  TH1D *histNTracks1 = new TH1D("hNTracks1" ,"Number of tracks about mu1, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); N_{trk} about muon1 / N (muon1)", 10,0,1.0);
  TH1D *histNTracks2OS = new TH1D("hNTracks2OS" ,"Number of tracks about mu2, OS, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); N_{trk} about muon2 / N (muon2)", 10,0,1.0);
  TH1D *histNTracks2 = new TH1D("hNTracks2" ,"Number of tracks about mu2, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); N_{trk} about muon2 / N (muon2)", 10,0,1.0);
  TH1D *histNTracksCum1OS = new TH1D("hNTracksCum1OS" ,"Cumu Number of tracks about mu1, OS,p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); N_{trk} about muon1 / N (muon1)", 10,0,1.0);
  TH1D *histNTracksCum1 = new TH1D("hNTracksCum1" ,"Cumu Number of tracks about mu1, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{1}-track); N_{trk} about muon1 / N (muon1)", 10,0,1.0);
  TH1D *histNTracksCum2OS = new TH1D("hNTracksCum2OS" ,"Cumu Number of tracks about mu2, OS, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); N_{trk} about muon2 / N (muon2)", 10,0,1.0);
  TH1D *histNTracksCum2 = new TH1D("hNTracksCum2" ,"Cumu Number of tracks about mu2, p_{T}(trk)>2.5 GeV, muon selection;#Delta R (#mu_{2}-track); N_{trk} about muon2 / N (muon2)", 10,0,1.0);
  
  TH1D *histMu1Pt = new TH1D("hMu1Pt", "#mu_{1} p_{T}, no selection ;#mu_{1} p_{T}; N_{events}", 50,0,50.);
  TH1D *histMu2Pt = new TH1D("hMu2Pt", "#mu_{2} p_{T}, no selection;#mu_{1} p_{T}; N_{events}", 50,0,50.);
  
  TH1D *histMu1PtSel = new TH1D("hMu1PtSel", "#mu_{1} p_{T}, selection ;#mu_{1} p_{T}; N_{events}", 50,0,50.);
  TH1D *histMu2PtSel = new TH1D("hMu2PtSel", "#mu_{2} p_{T}, selection;#mu_{1} p_{T}; N_{events}", 50,0,50.);
  
  TH1D *histNMu = new TH1D("hNMu", "No. muons;N mu; N_{events}", 5,0,5);

  TH1D *histNTk25 = new TH1D("hNTk25", "No. tracks, p_{T} > 2.5 GeV;N_{tk}; N_{events}", 25,0,50);
  TH1D *histNTk1 = new TH1D("hNTk1", "No. tracks, p_{T} > 1 GeV;N_{tk}; N_{events}", 25,0,100);
  TH1D *histNTk = new TH1D("hNTk", "No. tracks, p_{T} > 0 GeV;N_{tk}; N_{events}", 25,50,250);
  
  TH1D *histDRMuMu = new TH1D("hDRMuMu", "#Delta R(#mu-#mu), muon selection;#Delta R(#mu_{1}-#mu_{2}); N_{events}", 20,0,TMath::Pi());
  
  int nMu(0);
    int n1(0), n2(0), nMuPass(0);

  // Loop over all events
  // for(Int_t entry = 0; entry < 500; ++entry){
  // 
  cout << "Nevts : " << numberOfEntries <<endl;

  for(Int_t entry = 0; entry < numberOfEntries; ++entry){
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
      
    // Do at gen particle level.
    // First, find the highest pt muons.
    // Note this method works fine for highest pt muon but is slow for highest two.
    // Really we want an ordered list
    histNMu->Fill(branchGenMuons->GetEntries());
    if (branchGenMuons->GetEntries() < 2) continue;
    
    GenParticle *cand(0),*mu1(0), *mu2(0);
    Track *candTk(0);

    // cout << "Event" << endl;
    double mu1PT(0.), mu2PT(0.);
    for (int i = 0; i < branchGenMuons->GetEntries(); i++){
      cand = (GenParticle*) branchGenMuons->At(i);
      if (cand->PT > mu1PT) {
        mu1 = cand;
        mu1PT = cand->PT;
      }
    }

    for(int j = 0; j < branchGenMuons->GetEntries(); j++){
      cand = (GenParticle*) branchGenMuons->At(j);
      if ((cand->PT > mu2PT) && (cand->PT != mu1->PT)) {
        mu2 = cand;
        mu2PT = cand->PT;
      }
    }

    histMu1Pt->Fill(mu1PT);
    histMu2Pt->Fill(mu2PT);

    TLorentzVector mu1Mom, mu2Mom;
    mu1Mom = mu1->P4();
    mu2Mom = mu2->P4();

    ////////////////////
    // Muon selection //
    ////////////////////
    if ( (mu1PT > 17.) 
      && (mu2PT > 10.) 
      && ((mu1->Charge) == (mu2->Charge))
      && (abs(mu1->Eta) < 2.1)
      && (abs(mu2->Eta) < 2.4)
      && ((mu1Mom.DeltaR(mu2Mom)) > 2.) 
      ){

      histDRMuMu->Fill(mu1Mom.DeltaR(mu2Mom));

      if (mu1Mom.DeltaR(mu2Mom)>1.){
        
        /////////////////////////////////
        // Look at tracks around muons //
        /////////////////////////////////
        
        int nAroundMu1 (0), nAroundMu2(0);

        int nTk1(0), nTk25(0);
              
        n1++;
        n2++;
        
        // cout << "Track mult: " << branchTracks->GetEntries() << endl;
        histNTk->Fill(branchTracks->GetEntries());
        
        for(int a = 0; a < branchTracks->GetEntries(); a++){
          candTk = (Track*) branchTracks->At(a);
          
          if ( (candTk->PT != mu1->PT) // Check it isn't the same object as the muons!
            && (candTk->PT != mu2->PT) 
            && (candTk->PT > 1.)
            && (abs(candTk->Z) < 1.) //dz < 1mm
            && ((pow(candTk->X,2)+pow(candTk->Y,2)) < 1.) //dxy < 1mm
            ){

            nTk1++;
            if (candTk->PT > 2.5){
            nTk25++;
              double dR1 = (candTk->P4()).DeltaR(mu1Mom);
              double dR2 = (candTk->P4()).DeltaR(mu2Mom);
              
              // if ((candTk->Charge) * (mu1->Charge) < 0){ // only need one if statement because SS muons
                hNTracks1->Fill(dR1);
                hNTracks2->Fill(dR2);
              // } 
             
             if ((candTk->Charge) * (mu1->Charge) < 0){ // only need one if statement because SS muons
                hNTracks1OS->Fill(dR1);
                hNTracks2OS->Fill(dR2);
              }

              // Count number of tracks with pT > 1 within a cone of 0.5 about each muon
              if (dR1 < 0.5)
                nAroundMu1++;
              if (dR2 < 0.5)
                nAroundMu2++;
            } //end of 2.5 cut
          } // End of track selection
        } // End of track loop
        histNTk1->Fill(nTk1);
        histNTk25->Fill(nTk25);
        histMu1PtSel->Fill(mu1PT);
        histMu2PtSel->Fill(mu2PT);

        if (nAroundMu1==1 && nAroundMu2==1){
          nMuPass++;
        }
      }// end of deltaR cut
    } // end of muon selection
  
  } // end of event loop

  histNTracks1->Scale(1./n1);
  histNTracks1OS->Scale(1./n1);
  histNTracks2->Scale(1./n2);
  histNTracks2OS->Scale(1./n2);

  histNTracks1Cum = (TH1F*)histNTracks1->Clone();
  histNTracks1CumOS = (TH1F*)histNTracks1OS->Clone();
  histNTracks2Cum = (TH1F*)histNTracks2->Clone();
  histNTracks2CumOS = (TH1F*)histNTracks2OS->Clone();

  for (int i = 1; i <= histNTracks1->GetNbinsX(); i++){
    histNTracks1Cum->SetBinContent(i,histNTracks1Cum->GetBinContent(i-1) + histNTracks1->GetBinContent(i));
    histNTracks1CumOS->SetBinContent(i,histNTracks1CumOS->GetBinContent(i-1) + histNTracks1OS->GetBinContent(i));
    histNTracks2Cum->SetBinContent(i,histNTracks2Cum->GetBinContent(i-1) + histNTracks2->GetBinContent(i));
    histNTracks2CumOS->SetBinContent(i,histNTracks2CumOS->GetBinContent(i-1) + histNTracks2OS->GetBinContent(i));
  }

  cout << "n1: " << n1 << endl;
  cout << "n2: " << n2 << endl;
  cout << "nMuPass: " << nMuPass << endl;

  TCanvas c;
  std::string name("");
  std::string app("");
  if (doSignal) {
    // name = "Signal_";
    name = "Signal_1prong_";
    // name = "Signal_3prong_";
    app = "_sig";
  } else {
    app = "_bg";
    if (doMu)
      name = "QCDb_mu_";
    else
      name = "QCDb_";
  }
  histNMu->Draw("HISTE");
  c.SaveAs((name+"cleanTk/NMu_clean"+app+".pdf").c_str());

  histMu1Pt->Draw("HISTE");
  c.SaveAs((name+"cleanTk/Mu1Pt_clean"+app+".pdf").c_str());
  histMu2Pt->Draw("HISTE");
  c.SaveAs((name+"cleanTk/Mu2Pt_clean"+app+".pdf").c_str());
  
  histMu1PtSel->Draw("HISTE");
  c.SaveAs((name+"cleanTk/Mu1PtSel_clean"+app+".pdf").c_str());
  histMu2PtSel->Draw("HISTE");
  c.SaveAs((name+"cleanTk/Mu2PtSel_clean"+app+".pdf").c_str());
  
  histNTracks1->Draw("HISTE");
  c.SaveAs((name+"cleanTk/NTracks1_NS_clean"+app+".pdf").c_str());
  histNTracks2->Draw("HISTE");
  c.SaveAs((name+"cleanTk/NTracks2_NS_clean"+app+".pdf").c_str());

  histNTracks1OS->Draw("HISTE");
  c.SaveAs((name+"cleanTk/NTracks1_OS_clean"+app+".pdf").c_str());
  histNTracks2OS->Draw("HISTE");
  c.SaveAs((name+"cleanTk/NTracks2_OS_clean"+app+".pdf").c_str());

  histNTracks1Cum->Draw("HISTE");
  c.SaveAs((name+"cleanTk/NTracks1Cum_NS_clean"+app+".pdf").c_str());
  histNTracks2Cum->Draw("HISTE");
  c.SaveAs((name+"cleanTk/NTracks2Cum_NS_clean"+app+".pdf").c_str());

  histNTracks1CumOS->Draw("HISTE");
  c.SaveAs((name+"cleanTk/NTracks1Cum_OS_clean"+app+".pdf").c_str());
  histNTracks2CumOS->Draw("HISTE");
  c.SaveAs((name+"cleanTk/NTracks2Cum_OS_clean"+app+".pdf").c_str());

  histDRMuMu->Draw("HISTE");
  c.SaveAs((name+"cleanTk/DRMuMu_clean"+app+".pdf").c_str());

  histNTk->Draw("HISTE");
  c.SaveAs((name+"cleanTk/NTk_clean"+app+".pdf").c_str());
  histNTk1->Draw("HISTE");
  c.SaveAs((name+"cleanTk/NTk1_clean"+app+".pdf").c_str());
  histNTk25->Draw("HISTE");
  c.SaveAs((name+"cleanTk/NTk25_clean"+app+".pdf").c_str());

  TFile* outFile = TFile::Open((name+"cleanTk/output"+app+".root").c_str(),"RECREATE");
  histNMu->Write("",TObject::kOverwrite);
  histMu1Pt->Write("",TObject::kOverwrite);
  histMu2Pt->Write("",TObject::kOverwrite);
  histMu1PtSel->Write("",TObject::kOverwrite);
  histMu2PtSel->Write("",TObject::kOverwrite);
  histNTracks1->Write("",TObject::kOverwrite);
  histNTracks2->Write("",TObject::kOverwrite);
  histNTracks1Cum->Write("",TObject::kOverwrite);
  histNTracks2Cum->Write("",TObject::kOverwrite);
  outFile->Close();

}
