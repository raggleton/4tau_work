void combinePlots(){

	gStyle->SetOptStat("");

	TCanvas c1;
	TFile fSig("Signal_1prong_500K_bare/output_bare_sig.root","READ");
	TFile fBg("QCDb_mu_pthatmin20_bare/output_bare_bg.root","READ");
	TH1D* hSig = fSig.Get("hNTracksCum1");
	TH1D* hBg = fBg.Get("hNTracksCum1");
	hSig->Draw("HISTE");
	hBg->SetLineColor(kRed);
	hBg->Draw("HISTESAME");
	TLegend leg(0.7,0.7,0.9,0.9);
	leg.SetFillColor(kWhite);
	leg.SetLineWidth(0);
	leg.AddEntry(hSig,"Signal","l");
	leg.AddEntry(hBg,"Bg","l");
	leg.Draw();
	c1.SaveAs("combined_NTrackCum1.pdf");

	TH1D* hSig2 = fSig.Get("hNTracksCum2");
	TH1D* hBg2 = fBg.Get("hNTracksCum2");
	hSig2->Draw("HISTE");
	hBg2->SetLineColor(kRed);
	hBg2->Draw("HISTESAME");
	leg.Draw();
	c1.SaveAs("combined_NTrackCum2.pdf");

	fSig.Close();
	fBg.Close();

	TFile fSigRand("Signal_1prong_500K_bare/output_bare_sig_muRand.root","READ");
	TFile fBgRand("QCDb_mu_pthatmin20_bare/output_bare_bg_muRand.root","READ");
	TH1D* hSigRand = fSigRand.Get("hNTracksCum1");
	TH1D* hBgRand = fBgRand.Get("hNTracksCum1");
	hSigRand->Draw("HISTE");
	hBgRand->SetLineColor(kRed);
	hBgRand->Draw("HISTESAME");
	leg.Draw();
	c1.SaveAs("combined_NTrackCum1_muRand.pdf");

	TH1D* hSigRand2 = fSigRand.Get("hNTracksCum2");
	TH1D* hBgRand2 = fBgRand.Get("hNTracksCum2");
	hSigRand2->Draw("HISTE");
	hBgRand2->SetLineColor(kRed);
	hBgRand2->Draw("HISTESAME");
	leg.Draw();
	c1.SaveAs("combined_NTrackCum2_muRand.pdf");

	fSigRand.Close();
	fBgRand.Close();
}