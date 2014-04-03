// Basically takes two hists of same name from two files (fSig & fBg),
// plots on same canvas and saves to PDF
// 
// histName = name of histogram in ROOT file (assumed to be same for sig & bg)
// plotOpt  = plotting options for Draw() e.g. HISTE
// outputName  = name for output (needs .pdf, .png etc)
// 
void combineHists( TFile* fSig, TFile* fBg, std::string histName, std::string plotOpt, std::string outputName){

	TCanvas c1;
	TH1D* hSig = fSig->Get(histName.c_str());
	TH1D* hBg  = fBg->Get(histName.c_str());
	hSig->SetLineColor(kRed);
	THStack st("h","");
	st.Add(hSig);
	st.Add(hBg);
	// hSig->Draw(plotOpt.c_str());
	// hBg->Draw((plotOpt+"SAME").c_str());
	st.Draw((plotOpt+"NOSTACK").c_str());
	TLegend leg(0.75,0.75,0.88,0.88);
	leg.SetFillColor(kWhite);
	leg.SetLineColor(kWhite);
	leg.AddEntry(hSig,"Signal","l");
	leg.AddEntry(hBg,"Bg","l");
	leg.Draw();
	c1.SaveAs(outputName.c_str());

	if (!hSig) delete hSig;
	if (!hBg) delete hBg;
}

void combinePlots(){
	
	////////////////////////
	// Setup, open files //
	////////////////////////
	
	gStyle->SetOptStat("");

	TFile fSig("Signal_1prong_500K_bare/output_bare_sig.root","READ");
	TFile fBg("QCDb_mu_pthatmin20_Mu17_Mu8_bare/output_bare_bg.root","READ");
	
	TFile fSigRand("Signal_1prong_500K_bare/output_bare_sig_muRand.root","READ");
	TFile fBgRand("QCDb_mu_pthatmin20_Mu17_Mu8_bare/output_bare_bg_muRand.root","READ");

	//////////////////
	// Plot things //
	//////////////////

	// Cumulative track distr., for pT-ordered and random-ordered muons
	combineHists(&fSig, &fBg, "hNTracksCum1", "HISTE", "Combined/combined_NTrackCum1.pdf");
	combineHists(&fSig, &fBg, "hNTracksCum2", "HISTE", "Combined/combined_NTrackCum2.pdf");
	combineHists(&fSigRand, &fBgRand, "hNTracksCum1", "HISTE", "Combined/combined_NTrackCum1_muRand.pdf");
	combineHists(&fSigRand, &fBgRand, "hNTracksCum2", "HISTE", "Combined/combined_NTrackCum2_muRand.pdf");

	// Absolute track distr., for pT-ordered and random-ordered muons
	combineHists(&fSig, &fBg, "hNTracksAbs1", "HISTE", "Combined/combined_NTrackAbs1.pdf");
	combineHists(&fSig, &fBg, "hNTracksAbs2", "HISTE", "Combined/combined_NTrackAbs2.pdf");
	combineHists(&fSigRand, &fBgRand, "hNTracksAbs1", "HISTE", "Combined/combined_NTrackAbs1_muRand.pdf");
	combineHists(&fSigRand, &fBgRand, "hNTracksAbs2", "HISTE", "Combined/combined_NTrackAbs2_muRand.pdf");

	// Normalised track distr., for pT-ordered and random-ordered muons
	combineHists(&fSig, &fBg, "hNTracks1", "HISTE", "Combined/combined_NTrack1.pdf");
	combineHists(&fSig, &fBg, "hNTracks2", "HISTE", "Combined/combined_NTrack2.pdf");
	combineHists(&fSigRand, &fBgRand, "hNTracks1", "HISTE", "Combined/combined_NTrack1_muRand.pdf");
	combineHists(&fSigRand, &fBgRand, "hNTracks2", "HISTE", "Combined/combined_NTrack2_muRand.pdf");

	//////////////////
	// Close files //
	//////////////////
	
	fSig.Close();
	fBg.Close();
	fSigRand.Close();
	fBgRand.Close();
}