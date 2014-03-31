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

}