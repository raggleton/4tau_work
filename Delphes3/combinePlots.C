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
	st.GetXaxis()->SetTitle(hSig->GetXaxis()->GetTitle());
	st.GetYaxis()->SetTitle(hSig->GetYaxis()->GetTitle());
	st.SetTitle(hSig->GetTitle());
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

	TFile fSig("Signal_1prong_500K_bare/output_bare_sig_HLT.root","READ");
	TFile fBg("QCDb_mu_pthatmin20_Mu17_Mu8_bare/output_bare_bg_HLT.root","READ");
	
	TFile fSigRand("Signal_1prong_500K_bare/output_bare_sig_muRand_HLT.root","READ");
	TFile fBgRand("QCDb_mu_pthatmin20_Mu17_Mu8_bare/output_bare_bg_muRand_HLT.root","READ");

	//////////////////
	// Plot things //
	//////////////////

	// Cumulative track distr., for pT-ordered and random-ordered muons - no sign requirement
	combineHists(&fSig, &fBg, "hNTracksCum1", "HISTE", "Combined/combined_NTrackCum1.pdf");
	combineHists(&fSig, &fBg, "hNTracksCum2", "HISTE", "Combined/combined_NTrackCum2.pdf");
	combineHists(&fSigRand, &fBgRand, "hNTracksCum1", "HISTE", "Combined/combined_NTrackCum1_muRand.pdf");
	combineHists(&fSigRand, &fBgRand, "hNTracksCum2", "HISTE", "Combined/combined_NTrackCum2_muRand.pdf");

	// Cumulative track distr., for pT-ordered and random-ordered muons - opposite sign requirement
	combineHists(&fSig, &fBg, "hNTracksCum1OS", "HISTE", "Combined/combined_NTrackCum1OS.pdf");
	combineHists(&fSig, &fBg, "hNTracksCum2OS", "HISTE", "Combined/combined_NTrackCum2OS.pdf");
	combineHists(&fSigRand, &fBgRand, "hNTracksCum1OS", "HISTE", "Combined/combined_NTrackCum1OS_muRand.pdf");
	combineHists(&fSigRand, &fBgRand, "hNTracksCum2OS", "HISTE", "Combined/combined_NTrackCum2OS_muRand.pdf");

	// Absolute track distr., for pT-ordered and random-ordered muons - no sign requirement
	combineHists(&fSig, &fBg, "hNTracksAbs1", "HISTE", "Combined/combined_NTrackAbs1.pdf");
	combineHists(&fSig, &fBg, "hNTracksAbs2", "HISTE", "Combined/combined_NTrackAbs2.pdf");
	combineHists(&fSigRand, &fBgRand, "hNTracksAbs1", "HISTE", "Combined/combined_NTrackAbs1_muRand.pdf");
	combineHists(&fSigRand, &fBgRand, "hNTracksAbs2", "HISTE", "Combined/combined_NTrackAbs2_muRand.pdf");

	// Absolute track distr., for pT-ordered and random-ordered muons - opposite sign requirement
	combineHists(&fSig, &fBg, "hNTracksAbs1OS", "HISTE", "Combined/combined_NTrackAbs1OS.pdf");
	combineHists(&fSig, &fBg, "hNTracksAbs2OS", "HISTE", "Combined/combined_NTrackAbs2OS.pdf");
	combineHists(&fSigRand, &fBgRand, "hNTracksAbs1OS", "HISTE", "Combined/combined_NTrackAbs1OS_muRand.pdf");
	combineHists(&fSigRand, &fBgRand, "hNTracksAbs2OS", "HISTE", "Combined/combined_NTrackAbs2OS_muRand.pdf");

	// Normalised track distr., for pT-ordered and random-ordered muons - no sign requirement
	combineHists(&fSig, &fBg, "hNTracks1", "HISTE", "Combined/combined_NTrack1.pdf");
	combineHists(&fSig, &fBg, "hNTracks2", "HISTE", "Combined/combined_NTrack2.pdf");
	combineHists(&fSigRand, &fBgRand, "hNTracks1", "HISTE", "Combined/combined_NTrack1_muRand.pdf");
	combineHists(&fSigRand, &fBgRand, "hNTracks2", "HISTE", "Combined/combined_NTrack2_muRand.pdf");
	
	// Normalised track distr., for pT-ordered and random-ordered muons - opposite sign requirement
	combineHists(&fSig, &fBg, "hNTracks1OS", "HISTE", "Combined/combined_NTrackOS1.pdf");
	combineHists(&fSig, &fBg, "hNTracks2OS", "HISTE", "Combined/combined_NTrackOS2.pdf");
	combineHists(&fSigRand, &fBgRand, "hNTracks1OS", "HISTE", "Combined/combined_NTrackOS1_muRand.pdf");
	combineHists(&fSigRand, &fBgRand, "hNTracks2OS", "HISTE", "Combined/combined_NTrackOS2_muRand.pdf");

	//////////////////////
	// For soft tracks //
	//////////////////////
	// Cumulative track distr., for pT-ordered and random-ordered muons - no sign requirement
	combineHists(&fSigRand, &fBgRand, "hNSoftTracksCum1", "HISTE", "Combined/combined_NSoftTrackCum1_muRand.pdf");
	combineHists(&fSigRand, &fBgRand, "hNSoftTracksCum2", "HISTE", "Combined/combined_NSoftTrackCum2_muRand.pdf");

	// Cumulative track distr., for pT-ordered and random-ordered muons - opposite sign requirement
	combineHists(&fSigRand, &fBgRand, "hNSoftTracksCum1OS", "HISTE", "Combined/combined_NSoftTrackCum1OS_muRand.pdf");
	combineHists(&fSigRand, &fBgRand, "hNSoftTracksCum2OS", "HISTE", "Combined/combined_NSoftTrackCum2OS_muRand.pdf");

	// Absolute track distr., for pT-ordered and random-ordered muons - no sign requirement
	combineHists(&fSigRand, &fBgRand, "hNSoftTracksAbs1", "HISTE", "Combined/combined_NSoftTrackAbs1_muRand.pdf");
	combineHists(&fSigRand, &fBgRand, "hNSoftTracksAbs2", "HISTE", "Combined/combined_NSoftTrackAbs2_muRand.pdf");

	// Absolute track distr., for pT-ordered and random-ordered muons - opposite sign requirement
	combineHists(&fSigRand, &fBgRand, "hNSoftTracksAbs1OS", "HISTE", "Combined/combined_NSoftTrackAbs1OS_muRand.pdf");
	combineHists(&fSigRand, &fBgRand, "hNSoftTracksAbs2OS", "HISTE", "Combined/combined_NSoftTrackAbs2OS_muRand.pdf");

	// Normalised track distr., for pT-ordered and random-ordered muons - no sign requirement
	combineHists(&fSigRand, &fBgRand, "hNSoftTracks1", "HISTE", "Combined/combined_NSoftTrack1_muRand.pdf");
	combineHists(&fSigRand, &fBgRand, "hNSoftTracks2", "HISTE", "Combined/combined_NSoftTrack2_muRand.pdf");
	
	// Normalised track distr., for pT-ordered and random-ordered muons - opposite sign requirement
	combineHists(&fSigRand, &fBgRand, "hNSoftTracks1OS", "HISTE", "Combined/combined_NSoftTrackOS1_muRand.pdf");
	combineHists(&fSigRand, &fBgRand, "hNSoftTracks2OS", "HISTE", "Combined/combined_NSoftTrackOS2_muRand.pdf");

	//////////////////
	// Close files //
	//////////////////
	
	fSig.Close();
	fBg.Close();
	fSigRand.Close();
	fBgRand.Close();
}