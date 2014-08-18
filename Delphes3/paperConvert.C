#include <vector>



// Convert plots to paper format - no title, bigger fonts etc
void combineHists( TFile* fSig, TFile* fBg, TFile* fBg2, std::string histName, std::string plotOpt, std::string outputName, std::vector<double> scalingFactors, std::string label, std::string yTitle="DEFAULT"){
    // for 3 hists - sig and 2 bg
    TCanvas c1;
    TH1D* hSig = fSig->Get(histName.c_str());
    hSig->SetLineColor(kRed);
    hSig->SetMarkerSize(0);
    doSignalHist(hSig);
    
    // Make combined BG hist
    // Need to rescale carefully
    std::vector<TFile*> files;
    files.push_back(fBg);
    files.push_back(fBg2);
    // TH1D* hBg = combine(files, histName, scalingFactors);
    TH1D* hBgA = fBg->Get(histName.c_str());
    TH1D* hBgB =  fBg2->Get(histName.c_str());
    TH1D* hBg = (TH1D*) hBgA->Clone();
    double total = scalingFactors[0]+scalingFactors[1];
    hBg->Scale(scalingFactors[0]/total);
    hBg->Add(hBgB, scalingFactors[1]/total);
    hBg->SetMarkerSize(0);
    doAltBGHist(hBg);
    
    THStack st("h","");
    st.Add(hSig);
    st.Add(hBg);
    st.Draw((plotOpt+"NOSTACK").c_str());
    st.GetXaxis()->SetTitle(hSig->GetXaxis()->GetTitle());
    if (yTitle == "DEFAULT") {
        st.GetYaxis()->SetTitle(hSig->GetYaxis()->GetTitle());
    } else {
        st.GetYaxis()->SetTitle(yTitle.c_str());
    }
    setAltOffsetSizes(&st.GetHistogram());
    st.SetTitle("");
    st.Draw((plotOpt+"NOSTACK").c_str());

    TLegend* l_all = new TLegend(0.65,0.6,0.89,0.89);
    l_all->AddEntry(hBg,"Gen. level QCD MC","lp");
    l_all->AddEntry((TObject*)0,"(b#bar{b} + q-g scatter,",""); //null pointers for blank entries
    l_all->AddEntry((TObject*)0,"q = b, #bar{b}, c, #bar{c})","");
    l_all->AddEntry(hSig, "Signal MC", "lp");
    l_all->AddEntry((TObject*)0,"m_{#varphi} = 8 GeV", "");
    doStandardLegend(l_all);
    l_all->Draw();

    TPaveText t(0.15, 0.75, 0.5, 0.85, "NDC");
    t.AddText(label.c_str());
    doStandardText(&t);
    if (label != "") {
        t.Draw();
    }
    c1.SaveAs(outputName.c_str());

    if (!hSig) delete hSig;
    if (!hBg) delete hBg;
}

void combineHists( TFile* fSig, TFile* fBg, std::string histName, std::string plotOpt, std::string outputName, std::string label=""){

    TCanvas c1;
    TH1D* hSig = fSig->Get(histName.c_str());
    TH1D* hBg  = fBg->Get(histName.c_str());
    cout << "signal: " << hSig->Integral()  << " bg: " << hBg->Integral() << endl;
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
    h->SetMarkerColor(h->GetLineColor());
}

void doAltBGHist(TH1 *h) {
    doStandardHist(h);
    h->SetMarkerStyle(21);
    // h->SetMarkerColor(kBlue);
    h->SetLineStyle(3);
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
    leg->SetLineWidth(0);
}

void doStandardText(TPaveText *t) {
    t->SetFillColor(kWhite);
    t->SetBorderSize(0);
}

void setMassAUTitles(TH1 *h) {
    h->SetXTitle("m(#mu-tk) [GeV]");
    h->SetYTitle("A.U.");
    // h->SetMaximum(1);
    h->SetMinimum(0);
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

void normaliseHist(TH1* h) {
    TH1::SetDefaultSumw2();

    if (h->Integral() != 0) {
        h->Scale(1./h->Integral());
    }
}

TH1D* combineRebin4Bins(std::vector<TFile*> files, std::string histName, std::vector<double> scalingFactors) {
    TH1::SetDefaultSumw2();
    std::vector<double> massBins;
    massBins.push_back(0);
    massBins.push_back(1);
    massBins.push_back(2);
    massBins.push_back(3);
    massBins.push_back(10);
    int nBinsX = massBins.size()-1;
    // Get 1st hist in list
    TH1D* h = (TH1D*)files[0]->Get(histName.c_str())->Clone(files[0]->Get(histName.c_str())->GetName());
    h->Scale(scalingFactors[0]);
    
    for (unsigned i = 1; i < files.size(); i++) {
        TH1D* hTmp = (TH1D*) (files[i]->Get(histName.c_str()));
        h->Add(hTmp, scalingFactors[i]);
    }
    // /files[0]->Get(histName.c_str())->GetTitle()
    TH1D* hNew = h->Rebin(nBinsX, files[0]->Get(histName.c_str())->GetTitle(), &massBins[0]);
    return hNew;
}

TH1D* combineRebin10bins(std::vector<TFile*> files, std::string histName, std::vector<double> scalingFactors) {
    TH1::SetDefaultSumw2();
    std::vector<double> massBins;
    for (int a = 0; a <= 10; a++) {
        massBins.push_back(a);
    }
    int nBinsX = massBins.size()-1;
    // Get 1st hist in list
    TH1D* h = (TH1D*)files[0]->Get(histName.c_str())->Clone(files[0]->Get(histName.c_str())->GetName());
    h->Scale(scalingFactors[0]);
    
    for (unsigned i = 1; i < files.size(); i++) {
        TH1D* hTmp = (TH1D*) (files[i]->Get(histName.c_str()));
        h->Add(hTmp, scalingFactors[i]);
    }
    // /files[0]->Get(histName.c_str())->GetTitle()
    TH1D* hNew = h->Rebin(nBinsX, files[0]->Get(histName.c_str())->GetTitle(), &massBins[0]);
    return hNew;
}


TH1D* combine(std::vector<TFile*> files, std::string histName, std::vector<double> scalingFactors) {
    TH1::SetDefaultSumw2();
    // Get 1st hist in list
    TH1D* h = (TH1D*)files[0]->Get(histName.c_str())->Clone(files[0]->Get(histName.c_str())->GetName());
    h->Scale(scalingFactors[0]);
    for (unsigned i = 1; i < files.size(); i++) {
        TH1D* hTmp = (TH1D*) (files[i]->Get(histName.c_str()));
        h->Add(hTmp, scalingFactors[i]);
    }
    return h;
}

/////////////////
// MAIN SCRIPT //
/////////////////
void paperConvert() {
    gStyle->SetHistLineWidth(2);
    gStyle->SetOptStat("");
    TH1::SetDefaultSumw2();

    TFile *f_sig_main1 = TFile::Open("Signal_1prong_HLT_bare/output_main_bare_sig_muRand_HLT_dR1.root", "READ");
    TFile *f_sig_main2 = TFile::Open("Signal_1prong_HLT_bare/output_main_bare_sig_muRand_HLT_dR2.root", "READ");
    TFile *f_sig_mass1 = TFile::Open("Signal_1prong_HLT_bare/output_bare_sig_muRand_HLT_dR1.root", "READ");
    TFile *f_sig_mass2 = TFile::Open("Signal_1prong_HLT_bare/output_bare_sig_muRand_HLT_dR2.root", "READ");
    TFile *f_bg_main1 = TFile::Open("QCDb_HLT_bare/output_main_bare_bg_muRand_HLT_dR1.root", "READ");
    TFile *f_bg_main2 = TFile::Open("QCDb_HLT_bare/output_main_bare_bg_muRand_HLT_dR2.root", "READ");
    TFile *f_bg_mass1 = TFile::Open("QCDb_HLT_bare/output_bare_bg_muRand_HLT_dR1.root", "READ");
    TFile *f_bg_mass2 = TFile::Open("QCDb_HLT_bare/output_bare_bg_muRand_HLT_dR2.root", "READ");
    TFile *f_scatter_main1 = TFile::Open("QCDbcScatter_HLT_bare/output_main_bare_bg_muRand_HLT_dR1.root", "READ");
    TFile *f_scatter_main2 = TFile::Open("QCDbcScatter_HLT_bare/output_main_bare_bg_muRand_HLT_dR2.root", "READ");
    TFile *f_scatter_mass1 = TFile::Open("QCDbcScatter_HLT_bare/output_bare_bg_muRand_HLT_dR1.root", "READ");
    TFile *f_scatter_mass2 = TFile::Open("QCDbcScatter_HLT_bare/output_bare_bg_muRand_HLT_dR2.root", "READ");


    ////////////////
    // Signal region plots
    ////////////////
    std::vector<TFile*> filesdR2;
    filesdR2.push_back(f_bg_mass2);
    filesdR2.push_back(f_scatter_mass2);

    // bbbar only
    TH1D* hM1_bare_bg_muRand_HLT_dR2 = (TH1D*) f_bg_mass2->Get("hM1");
    TH1D* hM2_bare_bg_muRand_HLT_dR2 = (TH1D*) f_bg_mass2->Get("hM2");
    TH1D* hM_bare_bg_muRand_HLT_dR2 = (TH1D*)hM1_bare_bg_muRand_HLT_dR2->Clone();
    hM_bare_bg_muRand_HLT_dR2->Add(hM2_bare_bg_muRand_HLT_dR2);
    normaliseHist(hM_bare_bg_muRand_HLT_dR2); // as we added m1 + m2
    setMassAUTitles(hM_bare_bg_muRand_HLT_dR2);
    doBGHist(hM_bare_bg_muRand_HLT_dR2);
    setAltOffsetSizes(hM_bare_bg_muRand_HLT_dR2);
    hM_bare_bg_muRand_HLT_dR2->Draw("HISTE");
    c1->SaveAs("Combined/M_10bins_bare_bg_muRand_HLT_dR2.pdf");

    std::vector<double> scalingFactors;
    scalingFactors.push_back(2.9475);
    scalingFactors.push_back(2.6577);

    // // bbbar + scatter
    TH1D* hM1_bg_dR2 = combineRebin10bins(filesdR2, "hM1_unnormalised", scalingFactors);
    TH1D* hM2_bg_dR2 = combineRebin10bins(filesdR2, "hM2_unnormalised", scalingFactors);
    TH1D* hM_bg_dR2 = (TH1D*) hM1_bg_dR2->Clone() ;
    hM_bg_dR2->Add(hM2_bg_dR2);
    normaliseHist(hM_bg_dR2);
    setMassAUTitles(hM_bg_dR2);
    doBGHist(hM_bg_dR2);
    setAltOffsetSizes(hM_bg_dR2);
    hM_bg_dR2->Draw("HISTE");
    c1->SaveAs("Combined/M_10bins_bare_bg_both_muRand_HLT_dR2.pdf");

    // signal
    TH1D* hM1_bare_sig_muRand_HLT_dR2 = (TH1D*) f_sig_mass2->Get("hM1");
    TH1D* hM2_bare_sig_muRand_HLT_dR2 = (TH1D*) f_sig_mass2->Get("hM2");
    TH1D* hM_bare_sig_muRand_HLT_dR2 = (TH1D*)hM1_bare_sig_muRand_HLT_dR2->Clone();
    hM_bare_sig_muRand_HLT_dR2->Add(hM2_bare_sig_muRand_HLT_dR2);
    normaliseHist(hM_bare_sig_muRand_HLT_dR2);
    setMassAUTitles(hM_bare_sig_muRand_HLT_dR2);
    doSignalHist(hM_bare_sig_muRand_HLT_dR2);
    setAltOffsetSizes(hM_bare_sig_muRand_HLT_dR2);
    hM_bare_sig_muRand_HLT_dR2->Draw("HISTE");
    c1->SaveAs("Combined/M_10bins_bare_sig_muRand_HLT_dR2.pdf");


    // bbbar + signal
    THStack* st = new THStack("st","");
    st->Add(hM_bare_bg_muRand_HLT_dR2);
    st->Add(hM_bare_sig_muRand_HLT_dR2);
    TLegend* l = new TLegend(0.67, 0.67, 0.88, 0.88);
    l->AddEntry(hM_bare_bg_muRand_HLT_dR2, "QCD b#bar{b} MC", "lp");
    l->AddEntry(hM_bare_sig_muRand_HLT_dR2, "#splitline{Signal MC}{m_{#varphi} = 8 GeV}", "lp");
    doStandardLegend(l);
    st->Draw("NOSTACK HISTE");
    setMassAUTitles(st->GetHistogram());
    setAltOffsetSizes(st->GetHistogram());    
    l->Draw();
    c1->SaveAs("Combined/M_10bins_bare_both_muRand_HLT_dR2.pdf");

    // bbbar + scatter + signal
    THStack* st_all = new THStack("st_all","");
    st_all->Add(hM_bg_dR2);
    st_all->Add(hM_bare_sig_muRand_HLT_dR2);
    TLegend* l_all = new TLegend(0.55,0.56,0.89,0.89);
    l_all->AddEntry(hM_bg_dR2,"Gen. level QCD MC","lp");
    l_all->AddEntry((TObject*)0,"(b#bar{b} + q-g scatter,",""); //null pointers for blank entries
    l_all->AddEntry((TObject*)0,"q = b, #bar{b}, c, #bar{c})","");
    l_all->AddEntry(hM_bare_sig_muRand_HLT_dR2, "Signal MC", "lp");
    l_all->AddEntry((TObject*)0,"m_{#varphi} = 8 GeV", "");
    doStandardLegend(l_all);
    st_all->Draw("NOSTACKHISTE");
    setMassAUTitles(st_all->GetHistogram());
    setAltOffsetSizes(st_all->GetHistogram());    
    l_all->Draw();
    c1->SaveAs("Combined/M_10bins_bare_all_muRand_HLT_dR2.pdf");


    ////////////
    // sideband plots
    ////////////
    
    std::vector<TFile*> filesdR1;
    filesdR1.push_back(f_bg_mass1);
    filesdR1.push_back(f_scatter_mass1);

    // bbbar
    TH1D* hM1_side_bg_muRand_HLT_dR1 = (TH1D*) f_bg_mass1->Get("hM1_side_1to2p5");
    hM1_side_bg_muRand_HLT_dR1->Rebin(5);
    TH1D* hM2_side_bg_muRand_HLT_dR1 = (TH1D*) f_bg_mass1->Get("hM2_side_1to2p5");
    hM2_side_bg_muRand_HLT_dR1->Rebin(5);
    TH1D* hM_side_bg_muRand_HLT_dR1 = (TH1D*)hM1_side_bg_muRand_HLT_dR1->Clone();
    hM_side_bg_muRand_HLT_dR1->Add(hM2_side_bg_muRand_HLT_dR1);
    normaliseHist(hM_side_bg_muRand_HLT_dR1);
    setMassAUTitles(hM_side_bg_muRand_HLT_dR1);
    doBGHist(hM_side_bg_muRand_HLT_dR1);
    setAltOffsetSizes(hM_side_bg_muRand_HLT_dR1);
    hM_side_bg_muRand_HLT_dR1->Draw("HISTE");
    c1->SaveAs("Combined/M_10bins_side_bg_muRand_HLT_dR1.pdf");

    // // bbbar + scatter
    TH1D* hM1_side_bg_dR1 = combineRebin10bins(filesdR1, "hM1_side_1to2p5_unnormalised", scalingFactors);
    TH1D* hM2_side_bg_dR1 = combineRebin10bins(filesdR1, "hM1_side_1to2p5_unnormalised", scalingFactors);
    TH1D* hM_side_bg_dR1 = (TH1D*) hM1_side_bg_dR1->Clone();
    hM_side_bg_dR1->Add(hM2_side_bg_dR1);
    normaliseHist(hM_side_bg_dR1);
    setMassAUTitles(hM_side_bg_dR1);
    doBGHist(hM_side_bg_dR1);
    setAltOffsetSizes(hM_side_bg_dR1);
    hM_side_bg_dR1->Draw("HISTE");
    c1->SaveAs("Combined/M_10bins_side_bg_both_muRand_HLT_dR1.pdf");

    // signal
    TH1D* hM1_side_sig_muRand_HLT_dR1 = (TH1D*) f_sig_mass1->Get("hM1_side_1to2p5");
    hM1_side_sig_muRand_HLT_dR1->Rebin(5);
    TH1D* hM2_side_sig_muRand_HLT_dR1 = (TH1D*) f_sig_mass1->Get("hM2_side_1to2p5");
    hM2_side_sig_muRand_HLT_dR1->Rebin(5);
    TH1D* hM_side_sig_muRand_HLT_dR1 = (TH1D*)hM1_side_sig_muRand_HLT_dR1->Clone();
    hM_side_sig_muRand_HLT_dR1->Add(hM2_side_sig_muRand_HLT_dR1);
    normaliseHist(hM_side_sig_muRand_HLT_dR1);
    setMassAUTitles(hM_side_sig_muRand_HLT_dR1);
    doSignalHist(hM_side_sig_muRand_HLT_dR1);
    setAltOffsetSizes(hM_side_sig_muRand_HLT_dR1);
    hM_side_sig_muRand_HLT_dR1->Draw("HISTE");
    c1->SaveAs("Combined/M_10bins_side_sig_muRand_HLT_dR1.pdf");

    // bbar + signal
    THStack* st_side = new THStack("st_side","");
    st_side->Add(hM_side_bg_muRand_HLT_dR1);
    st_side->Add(hM_side_sig_muRand_HLT_dR1);
    st_side->Draw("NOSTACK HISTE");
    setMassAUTitles(st_side->GetHistogram());
    setAltOffsetSizes(st_side->GetHistogram());    
    l->Draw();
    c1->SaveAs("Combined/M_10bins_side_both_muRand_HLT_dR1.pdf");
    
    // bbar + scatter + signal
    THStack* st_side_all = new THStack("st_side_all","");
    st_side_all->Add(hM_side_bg_dR1);
    st_side_all->Add(hM_side_sig_muRand_HLT_dR1);
    st_side_all->Draw("NOSTACK HISTE");
    setMassAUTitles(st_side_all->GetHistogram());
    setAltOffsetSizes(st_side_all->GetHistogram());    
    l_all->Draw();
    c1->SaveAs("Combined/M_10bins_side_all_muRand_HLT_dR1.pdf");

    ///////////////////////
    // Correlation plots //
    ///////////////////////
    // TH1D* hCorr_bare_sig = (TH1D*) f_sig_mass2->Get("hCorr1D");
    // TH1D* hCorr_side_sig = (TH1D*) f_sig_mass1->Get("hCorr1D_side_1to2p5");
    // TH1D* hCorr_bare_bg = (TH1D*) f_bg_mass2->Get("hCorr1D");
    // TH1D* hCorr_side_bg = (TH1D*) f_bg_mass1->Get("hCorr1D_side_1to2p5");
    
    // TPaveText t(0.15, 0.7, 0.4, 0.8, "NDC");
    // t.AddText("Signal region");
    // t.SetFillColor(kWhite);
    // t.SetBorderSize(0);

    // hCorr_bare_sig->SetMaximum(2);
    // hCorr_side_sig->SetMaximum(2);
    // hCorr_bare_bg->SetMaximum(2);
    // hCorr_side_bg->SetMaximum(2);

    // setOffsetSizes(hCorr_bare_sig);
    // setOffsetSizes(hCorr_side_sig);
    // setOffsetSizes(hCorr_bare_bg);
    // setOffsetSizes(hCorr_side_bg);
    
    // TLine line(0,1,10,1);
    // line.SetLineStyle(2);
    // line.SetLineColor(12);
    // line.SetLineWidth(2);
    
    // doSignalHist(hCorr_bare_sig);
    // hCorr_bare_sig->Draw();
    // line.Draw();
    // t.Draw();
    // c1->SetTicks(1,1);
    // c1->SaveAs("Combined/Corr_bare_sig.pdf");
    
    // doBGHist(hCorr_bare_bg);
    // hCorr_bare_bg->Draw();
    // line.Draw();
    // t.Draw();
    // c1->SetTicks(1,1);
    // c1->SaveAs("Combined/Corr_bare_bg.pdf");
    
    // hCorr_bare_bg->Draw();
    // hCorr_bare_sig->Draw("SAME");
    // line.Draw();
    // l.Draw();
    // t.Draw();
    // c1->SetTicks(1,1);
    // c1->SaveAs("Combined/Corr_bare.pdf");

    // doSignalHist(hCorr_side_sig);
    // hCorr_side_sig->Draw();
    // line.Draw();
    // c1->SaveAs("Combined/Corr_side_sig.pdf");
    
    // doBGHist(hCorr_side_bg);
    // hCorr_side_bg->Draw();
    // line.Draw();
    // c1->SaveAs("Combined/Corr_side_bg.pdf");
    
    // hCorr_side_sig->Draw();
    // hCorr_side_bg->Draw("SAME");
    // line.Draw();
    // l.Draw();
    // c1->SaveAs("Combined/Corr_side.pdf");

    // track distributions
    combineHists(f_sig_main2, f_bg_main2, f_scatter_main2, "hNTracks1", "HISTE", "Combined/combined_NTrack1_muRand.pdf", scalingFactors, "Tracks with p_{T} > 2.5 GeV");
    combineHists(f_sig_main2, f_bg_main2, f_scatter_main2, "hNTracksAbs1", "HISTE", "Combined/combined_NTrackAbs1_muRand.pdf", scalingFactors, "Tracks with p_{T} > 2.5 GeV", "Average number of tracks per #mu_{1} / bin");
    combineHists(f_sig_main2, f_bg_main2, f_scatter_main2, "hNTracksAll1", "HISTE", "Combined/combined_NTrackAll1_muRand.pdf", scalingFactors, "Tracks with p_{T} > 1 GeV");
    combineHists(f_sig_main2, f_bg_main2, f_scatter_main2, "hNTracksAllAbs1", "HISTE", "Combined/combined_NTrackAllAbs1_muRand.pdf", scalingFactors, "Tracks with p_{T} > 1 GeV", "Average number of tracks per #mu_{1} / bin");
    combineHists(f_sig_main2, f_bg_main2, f_scatter_main2, "hNSoftTracks1", "HISTE", "Combined/combined_NSoftTrack1_muRand.pdf", scalingFactors, "Tracks with 2.5 > p_{T} > 1 GeV");
    combineHists(f_sig_main2, f_bg_main2, f_scatter_main2, "hNSoftTracksAbs1", "HISTE", "Combined/combined_NSoftTrackAbs1_muRand.pdf", scalingFactors, "Tracks with 2.5 > p_{T} > 1 GeV", "Average number of tracks per #mu_{1} / bin");


    // cleanup
    f_sig_main1->Close();
    f_sig_main2->Close();
    f_sig_mass2->Close();
    f_sig_mass2->Close();
    f_bg_main1->Close();
    f_bg_main2->Close();
    f_bg_mass1->Close();
    f_bg_mass2->Close();
    f_scatter_main1->Close();
    f_scatter_main2->Close();
    f_scatter_mass1->Close();
    f_scatter_mass2->Close();

}