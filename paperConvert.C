// Convert plots to paper format - no title, bigger fonts etc
void combineHists( TFile* fSig, TFile* fBg, TFile* fBg2, std::string histName, std::string plotOpt, std::string outputName, int rebin=1){
    // for 3 hists - sig and 2 bg
    TCanvas c1;
    TH1D* hSig = fSig->Get(histName.c_str());
    TH1D* hBg  = fBg->Get(histName.c_str());
    TH1D* hBg2  = fBg2->Get(histName.c_str());
    hSig->SetLineColor(kRed);
    hSig->SetMarkerSize(0);
    hSig->Rebin(rebin);
    hBg->SetMarkerSize(0);
    hBg->Rebin(rebin);
    hBg2->SetMarkerSize(0);
    hBg2->SetLineColor(kGreen+2);
    hBg2->Rebin(rebin);
    THStack st("h","");
    st.Add(hSig);
    st.Add(hBg);
    st.Add(hBg2);
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
    leg.AddEntry(hBg2,"Bg2","l");
    leg.Draw();
    c1.SaveAs(outputName.c_str());

    if (!hSig) delete hSig;
    if (!hBg) delete hBg;
}

void combineHists( TFile* fSig, TFile* fBg, std::string histName, std::string plotOpt, std::string outputName, std::string label=""){

    TCanvas c1;
    TH1D* hSig = fSig->Get(histName.c_str());
    TH1D* hBg  = fBg->Get(histName.c_str());
    hSig->SetLineColor(kRed);
    hSig->SetMarkerSize(0);
    hBg->SetMarkerSize(0);
    THStack st("h","");
    st.Add(hSig);
    st.Add(hBg);
    // hSig->Draw(plotOpt.c_str());
    // hBg->Draw((plotOpt+"SAME").c_str());
    st.Draw((plotOpt+"NOSTACK").c_str());
    st.GetXaxis()->SetTitle(hSig->GetXaxis()->GetTitle());
    st.GetXaxis()->SetTitleSize(0.04);
    st.GetXaxis()->SetLabelSize(0.04);
    // st.GetXaxis()->SetTitleOffset();
    st.GetYaxis()->SetTitle(hSig->GetYaxis()->GetTitle());
    st.GetYaxis()->SetTitleSize(0.04);
    st.GetYaxis()->SetLabelSize(0.04);
    st.GetYaxis()->SetTitleOffset(1.2);
    st.SetTitle("");
    st.Draw((plotOpt+"NOSTACK").c_str());
    TLegend leg(0.65,0.6, 0.85,0.85);
    leg.SetFillColor(kWhite);
    leg.SetLineColor(kWhite);
    leg.AddEntry(hSig,"#splitline{Signal MC}{m_{#varphi} = 8 GeV}","l");
    leg.AddEntry(hBg,"QCD b#bar{b} MC","l");
    leg.Draw();

    TPaveText t(0.15, 0.7, 0.5, 0.8, "NDC");
    t.AddText(label.c_str());
    t.SetFillColor(kWhite);
    t.SetBorderSize(0);
    if (label != "") {
        t.Draw();
    }
    c1.SaveAs(outputName.c_str());

    if (!hSig) delete hSig;
    if (!hBg) delete hBg;
}

void paperConvert() {

    gStyle->SetOptStat("");
    TH1::SetDefaultSumw2();

    TFile *f_sig_main = TFile::Open("Signal_1prong_HLT_bare/output_bare_sig_muRand_HLT.root", "READ");
    TFile *f_sig_mass1 = TFile::Open("Signal_1prong_HLT_bare/output_bare_sig_muRand_HLT_dR1.root", "READ");
    TFile *f_sig_mass2 = TFile::Open("Signal_1prong_HLT_bare/output_bare_sig_muRand_HLT_dR2.root", "READ");
    TFile *f_bg_main = TFile::Open("QCDb_HLT_bare/output_bare_bg_muRand_HLT.root", "READ");
    TFile *f_bg_mass1 = TFile::Open("QCDb_HLT_bare/output_bare_bg_muRand_HLT_dR1.root", "READ");
    TFile *f_bg_mass2 = TFile::Open("QCDb_HLT_bare/output_bare_bg_muRand_HLT_dR2.root", "READ");


    ////////////////
    // Signal region plots
    ////////////////
    
    TH1D* hM1_bare_bg_muRand_HLT_dR2 = (TH1D*) f_bg_mass2->Get("hM1");
    TH1D* hM2_bare_bg_muRand_HLT_dR2 = (TH1D*) f_bg_mass2->Get("hM2");
    TH1D* hM_bare_bg_muRand_HLT_dR2 = (TH1D*)hM1_bare_bg_muRand_HLT_dR2->Clone();
    hM_bare_bg_muRand_HLT_dR2->Add(hM2_bare_bg_muRand_HLT_dR2);
    hM_bare_bg_muRand_HLT_dR2->SetTitle("");
    hM_bare_bg_muRand_HLT_dR2->SetXTitle("m(#mu-tk) [GeV]");
    hM_bare_bg_muRand_HLT_dR2->SetYTitle("A.U.");
    hM_bare_bg_muRand_HLT_dR2->SetTitleOffset(1.2, "Y");
    hM_bare_bg_muRand_HLT_dR2->SetTitleSize(0.04, "X");
    hM_bare_bg_muRand_HLT_dR2->SetTitleSize(0.04, "Y");
    hM_bare_bg_muRand_HLT_dR2->Draw("HISTE");
    hM_bare_bg_muRand_HLT_dR2->Scale(0.5);
    c1->SaveAs("M_10bins_bare_bg_muRand_HLT_dR2.pdf");

    TH1D* hM1_bare_sig_muRand_HLT_dR2 = (TH1D*) f_sig_mass2->Get("hM1");
    TH1D* hM2_bare_sig_muRand_HLT_dR2 = (TH1D*) f_sig_mass2->Get("hM2");
    TH1D* hM_bare_sig_muRand_HLT_dR2 = (TH1D*)hM1_bare_sig_muRand_HLT_dR2->Clone();
    hM_bare_sig_muRand_HLT_dR2->Add(hM2_bare_sig_muRand_HLT_dR2);
    hM_bare_sig_muRand_HLT_dR2->SetTitle("");
    hM_bare_sig_muRand_HLT_dR2->SetXTitle("m(#mu-tk) [GeV]");
    hM_bare_sig_muRand_HLT_dR2->SetYTitle("A.U.");
    hM_bare_sig_muRand_HLT_dR2->SetTitleOffset(1.2, "Y");
    hM_bare_sig_muRand_HLT_dR2->SetTitleSize(0.04, "X");
    hM_bare_sig_muRand_HLT_dR2->SetTitleSize(0.04, "Y");
    hM_bare_sig_muRand_HLT_dR2->SetLineColor(kRed);
    hM_bare_sig_muRand_HLT_dR2->SetMarkerColor(kRed);
    hM_bare_sig_muRand_HLT_dR2->Draw("HISTE");
    hM_bare_sig_muRand_HLT_dR2->Scale(0.5);
    c1->SaveAs("M_10bins_bare_sig_muRand_HLT_dR2.pdf");

    THStack st("st","");
    st.Add(hM_bare_bg_muRand_HLT_dR2);
    st.Add(hM_bare_sig_muRand_HLT_dR2);
    TLegend l(0.6,0.6, 0.85,0.85);
    l.AddEntry(hM_bare_sig_muRand_HLT_dR2, "#splitline{Signal MC}{m_{#varphi} = 8 GeV}", "lp");
    l.AddEntry(hM_bare_bg_muRand_HLT_dR2, "QCD b#bar{b} MC", "lp");
    l.SetFillColor(kWhite);
    l.SetLineColor(kWhite);
    st.Draw("NOSTACK HISTE");
    st.GetHistogram()->SetXTitle("m(#mu-tk) [GeV]");
    st.GetHistogram()->SetYTitle("A.U.");
    st.GetHistogram()->SetTitleOffset(1.2, "Y");
    st.GetHistogram()->SetTitleSize(0.04, "X");
    st.GetHistogram()->SetTitleSize(0.04, "Y");
    
    l.Draw();
    c1->SaveAs("M_10bins_bare_both_muRand_HLT_dR2.pdf");

    ////////////
    // sideband plots
    ////////////

    TH1D* hM1_side_1to2p5_bg = (TH1D*) f_bg_mass1->Get("hM1_side_1to2p5");
    hM1_side_1to2p5_bg->Rebin(5);
    TH1D* hM2_side_1to2p5_bg = (TH1D*) f_bg_mass1->Get("hM2_side_1to2p5");
    hM2_side_1to2p5_bg->Rebin(5);
    TH1D* hM_side_1to2p5_bg = (TH1D*)hM1_side_1to2p5_bg->Clone();
    hM_side_1to2p5_bg->Add(hM2_side_1to2p5_bg);
    hM_side_1to2p5_bg->SetTitle("");
    hM_side_1to2p5_bg->SetXTitle("m(#mu-tk) [GeV]");
    hM_side_1to2p5_bg->SetYTitle("A.U.");
    hM_side_1to2p5_bg->SetTitleOffset(1.2, "Y");
    hM_side_1to2p5_bg->SetTitleSize(0.04, "X");
    hM_side_1to2p5_bg->SetTitleSize(0.04, "Y");
    hM_side_1to2p5_bg->Draw("HISTE");
    hM_side_1to2p5_bg->Scale(0.5);
    c1->SaveAs("M_10bins_side_bg_muRand_HLT_dR1.pdf");

    TH1D* hM1_side_1to2p5_sig = (TH1D*) f_sig_mass1->Get("hM1_side_1to2p5");
    hM1_side_1to2p5_sig->Rebin(5);
    TH1D* hM2_side_1to2p5_sig = (TH1D*) f_sig_mass1->Get("hM2_side_1to2p5");
    hM2_side_1to2p5_sig->Rebin(5);
    TH1D* hM_side_1to2p5_sig = (TH1D*)hM1_side_1to2p5_sig->Clone();
    hM_side_1to2p5_sig->Add(hM2_side_1to2p5_sig);
    hM_side_1to2p5_sig->SetTitle("");
    hM_side_1to2p5_sig->SetXTitle("m(#mu-tk) [GeV]");
    hM_side_1to2p5_sig->SetYTitle("A.U.");
    hM_side_1to2p5_sig->SetTitleOffset(1.2, "Y");
    hM_side_1to2p5_sig->SetTitleSize(0.04, "X");
    hM_side_1to2p5_sig->SetTitleSize(0.04, "Y");
    hM_side_1to2p5_sig->SetLineColor(kRed);
    hM_side_1to2p5_sig->SetMarkerColor(kRed);
    hM_side_1to2p5_sig->Draw("HISTE");
    hM_side_1to2p5_sig->Scale(0.5);
    c1->SaveAs("M_10bins_side_sig_muRand_HLT_dR1.pdf");

    THStack st_side("st_side","");
    st_side.Add(hM_side_1to2p5_bg);
    st_side.Add(hM_side_1to2p5_sig);
    // TLegend l(0.5,0.72, 0.85,0.88);
    // l.AddEntry(hM_side_1to2p5_sig, "#splitline{Signal MC}{m_{#varphi} = 8 GeV}", "lp");
    // l.AddEntry(hM_side_1to2p5_bg, "QCD b#bar{b} MC", "lp");
    // l.SetFillColor(kWhite);
    // l.SetLineColor(kWhite);
    st_side.Draw("NOSTACK HISTE");
    st_side.GetHistogram()->SetXTitle("m(#mu-tk) [GeV]");
    st_side.GetHistogram()->SetYTitle("A.U.");
    st_side.GetHistogram()->SetTitleOffset(1.2, "Y");
    st_side.GetHistogram()->SetTitleSize(0.04, "X");
    st_side.GetHistogram()->SetTitleSize(0.04, "Y");
    
    l.Draw();
    c1->SaveAs("M_10bins_side_both_muRand_HLT_dR1.pdf");

    // Correlation plots
    TH1D* hCorr_bare_sig = (TH1D*) f_sig_mass2->Get("hCorr1D");
    TH1D* hCorr_side_sig = (TH1D*) f_sig_mass1->Get("hCorr1D_side_1to2p5");
    TH1D* hCorr_bare_bg = (TH1D*) f_bg_mass2->Get("hCorr1D");
    TH1D* hCorr_side_bg = (TH1D*) f_bg_mass1->Get("hCorr1D_side_1to2p5");
    
    hCorr_bare_sig->SetMaximum(2);
    hCorr_side_sig->SetMaximum(2);
    hCorr_bare_bg->SetMaximum(2);
    hCorr_side_bg->SetMaximum(2);

    hCorr_bare_sig->SetLabelSize(0.07, "X");
    hCorr_bare_sig->SetLabelSize(0.05, "Y");
    hCorr_bare_sig->SetTitleSize(0.05, "X");
    // hCorr_bare_sig->SetTitleOffset(0.05, "X");
    hCorr_bare_sig->SetTitleSize(0.07, "Y");
    hCorr_bare_sig->SetTitleSize(0.05, "Y");
    // hCorr_bare_sig->SetTitleOffset(0.05, "Y");
    hCorr_side_sig->SetLabelSize(0.07, "X");
    hCorr_side_sig->SetLabelSize(0.05, "Y");
    hCorr_side_sig->SetTitleSize(0.05, "X");
    // hCorr_side_sig->SetTitleOffset(0.05, "X");
    hCorr_side_sig->SetTitleSize(0.05, "Y");
    // hCorr_side_sig->SetTitleOffset(0.05, "Y");
    hCorr_bare_bg->SetLabelSize(0.07, "X");
    hCorr_bare_bg->SetLabelSize(0.05, "Y");
    hCorr_bare_bg->SetTitleSize(0.05, "X");
    // hCorr_bare_bg->SetTitleOffset(0.05, "X");
    hCorr_bare_bg->SetTitleSize(0.05, "Y");
    // hCorr_bare_bg->SetTitleOffset(0.05, "Y");
    hCorr_side_bg->SetLabelSize(0.07, "X");
    hCorr_side_bg->SetLabelSize(0.05, "Y");
    hCorr_side_bg->SetTitleSize(0.05, "X");
    // hCorr_side_bg->SetTitleOffset(0.05, "X");
    hCorr_side_bg->SetTitleSize(0.05, "Y");
    // hCorr_side_bg->SetTitleOffset(0.05, "Y");
    
    TLine line(0,1,10,1);
    line.SetLineStyle(2);
    line.SetLineColor(12);
    line.SetLineWidth(2);
    
    hCorr_bare_sig->SetLineColor(kRed);
    hCorr_bare_sig->SetMarkerColor(kRed);
    hCorr_bare_sig->Draw();
    line.Draw();
    c1->SaveAs("Corr_bare_sig.pdf");
    
    hCorr_bare_bg->Draw();
    line.Draw();
    c1->SaveAs("Corr_bare_bg.pdf");
    
    hCorr_bare_bg->Draw();
    hCorr_bare_sig->Draw("SAME");
    line.Draw();
    l.Draw();
    c1->SaveAs("Corr_bare.pdf");

    hCorr_side_sig->SetLineColor(kRed);
    hCorr_side_sig->SetMarkerColor(kRed);
    hCorr_side_sig->Draw();
    line.Draw();
    c1->SaveAs("Corr_side_sig.pdf");
    
    hCorr_side_bg->Draw();
    line.Draw();
    c1->SaveAs("Corr_side_bg.pdf");
    
    hCorr_side_sig->Draw();
    hCorr_side_bg->Draw("SAME");
    line.Draw();
    l.Draw();
    c1->SaveAs("Corr_side.pdf");

    // track distributions
    combineHists(f_sig_main, f_bg_main, "hNTracks1", "HISTE", "combined_NTrack1_muRand.pdf", "Tracks with p_{T} > 2.5 GeV");
    combineHists(f_sig_main, f_bg_main, "hNTracksAbs1", "HISTE", "combined_NTrackAbs1_muRand.pdf", "Tracks with p_{T} > 2.5 GeV");
    combineHists(f_sig_main, f_bg_main, "hNTracksAll1", "HISTE", "combined_NTrackAll1_muRand.pdf", "Tracks with p_{T} > 1 GeV");
    combineHists(f_sig_main, f_bg_main, "hNTracksAllAbs1", "HISTE", "combined_NTrackAllAbs1_muRand.pdf", "Tracks with p_{T} > 1 GeV");
    combineHists(f_sig_main, f_bg_main, "hNSoftTracks1", "HISTE", "combined_NSoftTrack1_muRand.pdf", "Tracks with 2.5 > p_{T} > 1 GeV");
    combineHists(f_sig_main, f_bg_main, "hNSoftTracksAbs1", "HISTE", "combined_NSoftTrackAbs1_muRand.pdf", "Tracks with 2.5 > p_{T} > 1 GeV");


    // cleanup
    f_sig_main->Close();
    f_sig_mass2->Close();
    f_sig_mass2->Close();
    f_bg_main->Close();
    f_bg_mass1->Close();
    f_bg_mass2->Close();
}