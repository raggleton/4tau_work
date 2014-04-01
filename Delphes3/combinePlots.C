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
	hSig->Draw(plotOpt.c_str());
	hBg->SetLineColor(kRed);
	hBg->Draw((plotOpt+"SAME").c_str());
	TLegend leg(0.7,0.7,0.9,0.9);
	leg.SetFillColor(kWhite);
	leg.SetLineWidth(0);
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
	TFile fBg("QCDb_mu_pthatmin20_bare/output_bare_bg.root","READ");
	
	TFile fSigRand("Signal_1prong_500K_bare/output_bare_sig_muRand.root","READ");
	TFile fBgRand("QCDb_mu_pthatmin20_bare/output_bare_bg_muRand.root","READ");

	//////////////////
	// Plot things //
	//////////////////

	// Cumulative track distr., for pT-ordered and 	random-ordered muons
	combineHists(&fSig, &fBg, "hNTracksCum1", "HISTE", "combined_NTrackCum1.pdf");
	combineHists(&fSig, &fBg, "hNTracksCum2", "HISTE", "combined_NTrackCum2.pdf");
	combineHists(&fSigRand, &fBgRand, "hNTracksCum1", "HISTE", "combined_NTrackCum1_muRand.pdf");
	combineHists(&fSigRand, &fBgRand, "hNTracksCum2", "HISTE", "combined_NTrackCum2_muRand.pdf");

	//////////////////
	// Close files //
	//////////////////
	
	fSig.Close();
	fBg.Close();
	fSigRand.Close();
	fBgRand.Close();
}