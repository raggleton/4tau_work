#include <vector>
// Convert plots to paper format - no title, bigger fonts etc
void combineHists( TFile* fSig, TFile* fBg, TFile* fBg2, std::string histName, std::string plotOpt, std::string outputName, int rebin=1){
    // for 3 hists - sig and 2 bg
    TCanvas c1;
    TH1D* hSig = fSig->Get(histName.c_str());
    TH1D* hBg  = fBg->Get(histName.c_str());
    TH1D* hBg2  = fBg2->Get(histName.c_str());
    hSig->SetLineColor(kRed);
    hSig->SetMarkerSize(0);
    doSignalHist(hSig);
    hSig->Rebin(rebin);
    hBg->SetMarkerSize(0);
    hBg->SetMarkerStyle(21);
    hBg->SetLineStyle(2);
    hBg->Rebin(rebin);
    hBg2->SetMarkerSize(0);
    hBg2->SetLineColor(kGreen+2);
    hBg2->SetMarkerStyle(22);
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
    hSig->SetMarkerSize(0);
    doSignalHist(hSig);
    hBg->SetMarkerSize(0);
    doAltBGHist(hBg);
    
    THStack st("h","");
    st.Add(hSig);
    st.Add(hBg);
    st.Draw((plotOpt+"NOSTACK").c_str());
    st.GetXaxis()->SetTitle(hSig->GetXaxis()->GetTitle());
    st.GetYaxis()->SetTitle(hSig->GetYaxis()->GetTitle());
    setAltOffsetSizes(&st.GetHistogram());

    st.SetTitle("");
    st.Draw((plotOpt+"NOSTACK").c_str());
    TLegend leg(0.65,0.6, 0.85,0.85);
    doStandardLegend(&leg);
    leg.AddEntry(hSig,"#splitline{Signal MC}{m_{#varphi} = 8 GeV}","l");
    leg.AddEntry(hBg,"QCD b#bar{b} MC","l");
    leg.Draw();

    TPaveText t(0.15, 0.7, 0.5, 0.8, "NDC");
    t.AddText(label.c_str());
    doStandardText(&t);
    if (label != "") {
        t.Draw();
    }
    c1.SaveAs(outputName.c_str());

    if (!hSig) delete hSig;
    if (!hBg) delete hBg;
}

void doStandardHist(TH1 *h) {
    h->SetLineWidth(2);
    h->SetTitle("");
}

void doSignalHist(TH1 *h) {
    doStandardHist(h);
    h->SetLineColor(kRed);
    h->SetMarkerColor(kRed);
}

void doBGHist(TH1 *h) {
    doStandardHist(h);
    h->SetMarkerStyle(21);
}

void doAltBGHist(TH1 *h) {
    doStandardHist(h);
    h->SetMarkerStyle(21);
    h->SetLineStyle(2);
}

void setOffsetSizes(TH1 *h) {
    h->SetLabelSize(0.07,"X");
    h->SetTitleSize(0.05,"X");
    h->SetLabelSize(0.05,"Y");
    h->SetTitleSize(0.05,"Y");
}

void setAltOffsetSizes(TH1 *h) {
    h->SetLabelSize(0.04,"X");
    h->SetTitleSize(0.04,"X");
    h->SetLabelSize(0.04,"Y");
    h->SetTitleSize(0.04,"Y");
    h->SetTitleOffset(1.2,"Y");
}

void doStandardLegend(TLegend *leg) {
    leg->SetFillStyle(0);
    leg->SetLineColor(kWhite);
}

void doStandardText(TPaveText *t) {
    t->SetFillColor(kWhite);
    t->SetBorderSize(0);
}

void setMassAUTitles(TH1 *h) {
    h->SetXTitle("m(#mu-tk) [GeV]");
    h->SetYTitle("A.U.");
}

TH1D* addScaleRebin(std::vector<TH1D*> plots, std::vector<double> scalingFactors) {
    TH1::SetDefaultSumw2();

    TH1D* h = (TH1D*)plots[0]->Clone(plots[0]->GetName());
    h->Scale(scalingFactors[0]);
    for (unsigned i = 1; i < plots.size(); i++) {
        TH1D* hTmp = (TH1D*)plots[i]->Rebin(nBinsX, plots[0]->GetTitle(), &massBins[0]);
        h->Add(hTmp, scalingFactors[i]);
    }
    return h;
}

TH1D* addScale(std::vector<TH1D*> plots, std::vector<double> scalingFactors) {
    TH1::SetDefaultSumw2();

    TH1D* h = (TH1D*)plots[0]->Clone(plots[0]->GetName());
    h->Scale(scalingFactors[0]);
    for (unsigned i = 1; i < plots.size(); i++) {
        // TH1D* hTmp = (TH1D*)plots[i]->Rebin(nBinsX, plots[0]->GetTitle(), &massBins[0]);
        h->Add(hTmp, scalingFactors[i]);
    }
    return h;
}

void paperConvert() {
    gStyle->SetHistLineWidth(2);
    gStyle->SetOptStat("");
    TH1::SetDefaultSumw2();

    TFile *f_sig_main = TFile::Open("Signal_1prong_HLT_bare/output_bare_sig_muRand_HLT.root", "READ");
    TFile *f_sig_mass1 = TFile::Open("Signal_1prong_HLT_bare/output_bare_sig_muRand_HLT_dR1.root", "READ");
    TFile *f_sig_mass2 = TFile::Open("Signal_1prong_HLT_bare/output_bare_sig_muRand_HLT_dR2.root", "READ");
    TFile *f_bg_main = TFile::Open("QCDb_HLT_bare/output_bare_bg_muRand_HLT.root", "READ");
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
    setMassAUTitles(hM_bare_bg_muRand_HLT_dR2);
    doBGHist(hM_bare_bg_muRand_HLT_dR2);
    setAltOffsetSizes(hM_bare_bg_muRand_HLT_dR2);
    hM_bare_bg_muRand_HLT_dR2->Draw("HISTE");
    hM_bare_bg_muRand_HLT_dR2->Scale(0.5); // as we added m1 + m2
    c1->SaveAs("Combined/M_10bins_bare_bg_muRand_HLT_dR2.pdf");

    TH1D* hM1_bare_sig_muRand_HLT_dR2 = (TH1D*) f_sig_mass2->Get("hM1");
    TH1D* hM2_bare_sig_muRand_HLT_dR2 = (TH1D*) f_sig_mass2->Get("hM2");
    TH1D* hM_bare_sig_muRand_HLT_dR2 = (TH1D*)hM1_bare_sig_muRand_HLT_dR2->Clone();
    hM_bare_sig_muRand_HLT_dR2->Add(hM2_bare_sig_muRand_HLT_dR2);
    setMassAUTitles(hM_bare_sig_muRand_HLT_dR2);
    doSignalHist(hM_bare_sig_muRand_HLT_dR2);
    setAltOffsetSizes(hM_bare_sig_muRand_HLT_dR2);
    hM_bare_sig_muRand_HLT_dR2->Draw("HISTE");
    hM_bare_sig_muRand_HLT_dR2->Scale(0.5);
    c1->SaveAs("Combined/M_10bins_bare_sig_muRand_HLT_dR2.pdf");

    THStack st("st","");
    st.Add(hM_bare_bg_muRand_HLT_dR2);
    st.Add(hM_bare_sig_muRand_HLT_dR2);
    TLegend l(0.67,0.67, 0.88,0.88);
    l.AddEntry(hM_bare_bg_muRand_HLT_dR2, "QCD b#bar{b} MC", "lp");
    l.AddEntry(hM_bare_sig_muRand_HLT_dR2, "#splitline{Signal MC}{m_{#varphi} = 8 GeV}", "lp");
    doStandardLegend(&l);
    st.Draw("NOSTACK HISTE");
    st.GetHistogram()->SetXTitle("m(#mu-tk) [GeV]");
    st.GetHistogram()->SetYTitle("A.U.");
    setAltOffsetSizes(st.GetHistogram());    
    l.Draw();
    c1->SaveAs("Combined/M_10bins_bare_both_muRand_HLT_dR2.pdf");

    ////////////
    // sideband plots
    ////////////
    
    TH1D* hM1_side_bg_muRand_HLT_dR2 = (TH1D*) f_bg_mass1->Get("hM1_side_1to2p5");
    hM1_side_bg_muRand_HLT_dR2->Rebin(5);
    TH1D* hM2_side_bg_muRand_HLT_dR2 = (TH1D*) f_bg_mass1->Get("hM2_side_1to2p5");
    hM2_side_bg_muRand_HLT_dR2->Rebin(5);
    TH1D* hM_side_bg_muRand_HLT_dR2 = (TH1D*)hM1_side_bg_muRand_HLT_dR2->Clone();
    hM_side_bg_muRand_HLT_dR2->Add(hM2_side_bg_muRand_HLT_dR2);
    setMassAUTitles(hM_side_bg_muRand_HLT_dR2);
    doBGHist(hM_side_bg_muRand_HLT_dR2);
    setAltOffsetSizes(hM_side_bg_muRand_HLT_dR2);
    hM_side_bg_muRand_HLT_dR2->Draw("HISTE");
    hM_side_bg_muRand_HLT_dR2->Scale(0.5); // as we added m1 + m2
    c1->SaveAs("Combined/M_10bins_side_bg_muRand_HLT_dR1.pdf");

    TH1D* hM1_side_sig_muRand_HLT_dR2 = (TH1D*) f_sig_mass1->Get("hM1_side_1to2p5");
    hM1_side_sig_muRand_HLT_dR2->Rebin(5);
    TH1D* hM2_side_sig_muRand_HLT_dR2 = (TH1D*) f_sig_mass1->Get("hM2_side_1to2p5");
    hM2_side_sig_muRand_HLT_dR2->Rebin(5);
    TH1D* hM_side_sig_muRand_HLT_dR2 = (TH1D*)hM1_side_sig_muRand_HLT_dR2->Clone();
    hM_side_sig_muRand_HLT_dR2->Add(hM2_side_sig_muRand_HLT_dR2);
    setMassAUTitles(hM_side_sig_muRand_HLT_dR2);
    doSignalHist(hM_side_sig_muRand_HLT_dR2);
    setAltOffsetSizes(hM_side_sig_muRand_HLT_dR2);
    hM_side_sig_muRand_HLT_dR2->Draw("HISTE");
    hM_side_sig_muRand_HLT_dR2->Scale(0.5);
    c1->SaveAs("Combined/M_10bins_side_sig_muRand_HLT_dR1.pdf");

    THStack st("st","");
    st.Add(hM_side_bg_muRand_HLT_dR2);
    st.Add(hM_side_sig_muRand_HLT_dR2);
    st.Draw("NOSTACK HISTE");
    st.GetHistogram()->SetXTitle("m(#mu-tk) [GeV]");
    st.GetHistogram()->SetYTitle("A.U.");
    setAltOffsetSizes(st.GetHistogram());    
    l.Draw();
    c1->SaveAs("Combined/M_10bins_side_both_muRand_HLT_dR1.pdf");

    ///////////////////////
    // Correlation plots //
    ///////////////////////
    TH1D* hCorr_bare_sig = (TH1D*) f_sig_mass2->Get("hCorr1D");
    TH1D* hCorr_side_sig = (TH1D*) f_sig_mass1->Get("hCorr1D_side_1to2p5");
    TH1D* hCorr_bare_bg = (TH1D*) f_bg_mass2->Get("hCorr1D");
    TH1D* hCorr_side_bg = (TH1D*) f_bg_mass1->Get("hCorr1D_side_1to2p5");
    
    TPaveText t(0.15, 0.7, 0.4, 0.8, "NDC");
    t.AddText("Signal region");
    t.SetFillColor(kWhite);
    t.SetBorderSize(0);

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
    hCorr_bare_sig->SetLineWidth(2);
    hCorr_bare_sig->SetMarkerColor(kRed);
    hCorr_bare_sig->Draw();
    line.Draw();
    t.Draw();
    c1->SetTicks(1,1);
    c1->SaveAs("Combined/Corr_bare_sig.pdf");
    
    hCorr_bare_bg->SetLineWidth(2);
    // hCorr_bare_bg->SetLineStyle(2);
    hCorr_bare_bg->SetMarkerStyle(21);
    hCorr_bare_bg->Draw();
    line.Draw();
    t.Draw();
    c1->SetTicks(1,1);
    c1->SaveAs("Combined/Corr_bare_bg.pdf");
    
    hCorr_bare_bg->Draw();
    hCorr_bare_sig->Draw("SAME");
    line.Draw();
    l.Draw();
    t.Draw();
    c1->SetTicks(1,1);
    c1->SaveAs("Combined/Corr_bare.pdf");

    hCorr_side_sig->SetLineColor(kRed);
    hCorr_side_sig->SetMarkerColor(kRed);
    hCorr_side_sig->Draw();
    line.Draw();
    c1->SaveAs("Combined/Corr_side_sig.pdf");
    
    hCorr_side_bg->Draw();
    line.Draw();
    c1->SaveAs("Combined/Corr_side_bg.pdf");
    
    hCorr_side_sig->Draw();
    hCorr_side_bg->Draw("SAME");
    line.Draw();
    l.Draw();
    c1->SaveAs("Combined/Corr_side.pdf");

    // track distributions
    combineHists(f_sig_main, f_bg_main, "hNTracks1", "HISTE", "Combined/combined_NTrack1_muRand.pdf", "Tracks with p_{T} > 2.5 GeV");
    combineHists(f_sig_main, f_bg_main, "hNTracksAbs1", "HISTE", "Combined/combined_NTrackAbs1_muRand.pdf", "Tracks with p_{T} > 2.5 GeV");
    combineHists(f_sig_main, f_bg_main, "hNTracksAll1", "HISTE", "Combined/combined_NTrackAll1_muRand.pdf", "Tracks with p_{T} > 1 GeV");
    combineHists(f_sig_main, f_bg_main, "hNTracksAllAbs1", "HISTE", "Combined/combined_NTrackAllAbs1_muRand.pdf", "Tracks with p_{T} > 1 GeV");
    combineHists(f_sig_main, f_bg_main, "hNSoftTracks1", "HISTE", "Combined/combined_NSoftTrack1_muRand.pdf", "Tracks with 2.5 > p_{T} > 1 GeV");
    combineHists(f_sig_main, f_bg_main, "hNSoftTracksAbs1", "HISTE", "Combined/combined_NSoftTrackAbs1_muRand.pdf", "Tracks with 2.5 > p_{T} > 1 GeV");


    // cleanup
    f_sig_main->Close();
    f_sig_mass2->Close();
    f_sig_mass2->Close();
    f_bg_main->Close();
    f_bg_mass1->Close();
    f_bg_mass2->Close();
}