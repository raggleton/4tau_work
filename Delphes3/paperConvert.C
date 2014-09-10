#include <vector>



// Convert plots to paper format - no title, bigger fonts etc
void combineHists( TFile* fSig, TFile* fBg, TFile* fBg2, std::string histName, std::string plotOpt, std::string outputName, std::vector<double> scalingFactors, std::string label, std::string yTitle="DEFAULT"){
    TH1::SetDefaultSumw2();

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
    setAltTitleLabelSizes(&st.GetHistogram());
    st.SetTitle("");
    st.Draw((plotOpt+"NOSTACK").c_str());

    TLegend* l_all = new TLegend(0.65,0.6,0.89,0.89);
    l_all->AddEntry(hBg,"Gen. level QCD MC","lp");
    l_all->AddEntry((TObject*)0,"(b#bar{b} + q-g scatter,",""); //null pointers for blank entries
    l_all->AddEntry((TObject*)0,"q = b, #bar{b}, c, #bar{c})","");
    l_all->AddEntry(hSig, "Signal MC", "lp");
    l_all->AddEntry((TObject*)0,"m_{#phi} = 8 GeV", "");
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
    TH1::SetDefaultSumw2();

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
    setAltTitleLabelSizes(&st.GetHistogram());

    st.SetTitle("");
    st.Draw((plotOpt+"NOSTACK").c_str());
    TLegend leg(0.65,0.6, 0.85,0.85);
    doStandardLegend(&leg);
    leg.AddEntry(hSig,"#splitline{Signal MC}{m_{#phi} = 8 GeV}","l");
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

std::string intToString(int n) {
    std::ostringstream convert;
    convert << n;
    return convert.str();
}

void doStandardHist(TH1 *h) {
    h->SetLineWidth(2);
    h->SetTitle("");
}

void doCustomHist(TH1* h, int color, int style = 1) {
    doStandardHist(h);
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(style);
}

void doSignalHist(TH1 *h) {
    doStandardHist(h);
    h->SetLineColor(kRed);
    h->SetMarkerColor(kRed);
}

void doDataHist(TH1* h) {
    h->SetLineColor(kGreen+2);
    h->SetLineWidth(2);
    h->SetMarkerColor(kGreen+2);
    h->SetMarkerStyle(22);
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

void setTitleLabelSizes(TH1 *h) {
    h->SetLabelSize(0.07,"X");
    h->SetTitleSize(0.05,"X");
    h->SetLabelSize(0.05,"Y");
    h->SetTitleSize(0.05,"Y");
    h->SetTitleOffset(1.2,"Y");
}

void setAltTitleLabelSizes(TH1 *h) {
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

std::vector<double> generate4Bins() {
    std::vector<double> massBins;
    massBins.push_back(0);
    massBins.push_back(1);
    massBins.push_back(2);
    massBins.push_back(3);
    massBins.push_back(10);
    return massBins;
}

std::vector<double> generate10Bins() {
    std::vector<double> massBins;
    for (int i = 0; i <= 10; i++) {
        massBins.push_back(i);
    }
    return massBins;
}

TH1D* combineRebin4Bins(std::vector<TFile*> files, std::string histName, std::vector<double> scalingFactors) {
    TH1::SetDefaultSumw2();
    std::vector<double> massBins = generate4Bins();
    // massBins.push_back(0);
    // massBins.push_back(1);
    // massBins.push_back(2);
    // massBins.push_back(3);
    // massBins.push_back(10);
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
    std::vector<double> massBins = generate10Bins();
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

/**
 * Combine various hists according to their weights
 * @param  hists          [description]
 * @param  scalingFactors [description]
 * @return                [description]
 */
TH1D* combineScale(std::vector<TH1D*> hists, std::vector<double> scalingFactors) {
    TH1::SetDefaultSumw2();
    // get total
    double total = 0;
    for (int i = 0; i < scalingFactors.size(); i++) {
        total += scalingFactors[i];
    }
    TH1D* h = (TH1D*) hists[0]->Clone(hists[0]->GetName());
    h->Scale(scalingFactors[0]/total);

    for (i = 1; i < scalingFactors.size(); i++) {
        h->Add(hists[i], scalingFactors[i]/total);
    }

    return h;
}

TH2D* combineScale(std::vector<TH2D*> hists, std::vector<double> scalingFactors) {
    TH1::SetDefaultSumw2();

    // get total
    double total = 0;
    for (int i = 0; i < scalingFactors.size(); i++) {
        total += scalingFactors[i];
    }
    TH2D* h = (TH2D*) hists[0]->Clone(hists[0]->GetName());
    h->Scale(scalingFactors[0]/total);

    for (i = 1; i < scalingFactors.size(); i++) {
        h->Add(hists[i], scalingFactors[i]/total);
    }

    return h;
}


///////////////////////////
// For correlation plots //
///////////////////////////

TH1D* dataCorr_side() {
    std::vector<double> massBins = generate4Bins();
    TH1D* histCorr1D_side_1to2p5_Data = new TH1D("hCorr_side_data","",10,0,10);
    double arrContents[] = {0.95, 1.02, 1.02, 1.10, 1.01, 0.96, 0.93, 1.07, 0.98, 0.91};
    double arrErrors[]   = {0.04, 0.03, 0.05, 0.07, 0.03, 0.05, 0.06, 0.12, 0.10, 0.19};
    for (unsigned i = 1; i <= histCorr1D_side_1to2p5_Data->GetNbinsX(); i++) {
        histCorr1D_side_1to2p5_Data->SetBinContent(i,arrContents[i-1]);
        histCorr1D_side_1to2p5_Data->SetBinError(i,arrErrors[i-1]);
    }
    return histCorr1D_side_1to2p5_Data;
}

TH2D* create2Dfrom1D(std::vector<double> bins, TH1D* h) {
    int nBins = bins.size()-1;
    TH2D* hist2D = new TH2D("", "", nBins,&bins[0],nBins,&bins[0]);
    for(int a = 1; a <= nBins; a++){
        for (int b = 1; b <=nBins; b++){
            hist2D->SetBinContent(a,b,h->GetBinContent(a)*h->GetBinContent(b));
            hist2D->SetBinError(a,b,sqrt(pow(h->GetBinContent(b)*h->GetBinError(a),2)
                                                            +pow(h->GetBinContent(a)*h->GetBinError(b),2)));
        }
    }
    return hist2D;
}

TH1D* unique1DBinsFrom2D(TH2D* h2D, int nUniqueBins) {
    TH1D* h1D = new TH1D("","",nUniqueBins,0,nUniqueBins);
    int counter = 1;
    int nBinsX = generate4Bins().size()-1;
    for (int i = 1; i <= nBinsX; i++) {
        for (int j = i; j <= nBinsX; j++) {
            std::string binLabel = "(" + intToString(i) + "," + intToString(j) + ")";
            h1D->SetBinContent(counter,h2D->GetBinContent(i,j));
            h1D->SetBinError(counter,h2D->GetBinError(i,j));
            h1D->GetXaxis()->SetBinLabel(counter,binLabel.c_str());
            counter++;
        }
    }
    return h1D;
}

void setCorrTitlesMaxMin(TH1D* h) {
    h->SetXTitle("Bin");
    h->SetYTitle("Correlation coefficiant C_{(i,j)}");
    h->SetMaximum(2.0);
    h->SetMinimum(0);
}

TH1D* createCorrelationPlot(TH2D* numerator2D, TH1D* denominator1D) {
    TH1::SetDefaultSumw2();
    normaliseHist(denominator1D);
    std::vector<double> massBins = generate4Bins();
    TH2D* denominator2D = create2Dfrom1D(massBins, denominator1D);
    normaliseHist(numerator2D);
    TH2D* corr2D = (TH2D*) numerator2D->Clone();
    corr2D->Divide(denominator2D);
    int nBinsX = massBins.size()-1;
    int nUniqueBins = (nBinsX+1)*nBinsX/2.;
    TH1D* corr1D = unique1DBinsFrom2D(corr2D, nUniqueBins);
    setCorrTitlesMaxMin(corr1D);
    return corr1D;
}


/////////////////
// MAIN SCRIPT //
/////////////////
void paperConvert() {
    gStyle->SetHistLineWidth(2);
    gStyle->SetOptStat("");
    TH1::SetDefaultSumw2();
    TCanvas *c1 = new TCanvas("c1", "c1",10,32,700,500);
    gStyle->SetOptStat(0);

    // For centralised hist with no title, legend inside plot
    c1->Range(-1.705566,0.8486743,10.78995,1.08077);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetTickx(1);
    c1->SetTicky(1);
    c1->SetLeftMargin(0.1364942);
    c1->SetRightMargin(0.06321839);
    c1->SetTopMargin(0.07415254);
    c1->SetBottomMargin(0.1271186);
    c1->SetFrameBorderMode(0);
    c1->SetFrameBorderMode(0);

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

    ///////////////////
    // 1D MASS PLOTS //
    ///////////////////

    ////////////////
    // Signal region plots
    ////////////////
    std::vector<TFile*> filesdR2;
    filesdR2.push_back(f_bg_mass2);
    filesdR2.push_back(f_scatter_mass2);

    TPaveText* t_signal = new TPaveText(0.6, 0.47, 0.85, 0.58, "NDC");
    t_signal->AddText("Signal region");
    t_signal->SetFillColor(kWhite);
    t_signal->SetBorderSize(0);
    
    TPaveText* t_side = new TPaveText(0.6, 0.47, 0.84, 0.57, "NDC");
    t_side->AddText("Control region B");
    t_side->SetFillColor(kWhite);
    t_side->SetBorderSize(0);

    // bbbar only
    TH1D* hM1_bare_bg_muRand_HLT_dR2 = (TH1D*) f_bg_mass2->Get("hM1");
    TH1D* hM2_bare_bg_muRand_HLT_dR2 = (TH1D*) f_bg_mass2->Get("hM2");
    TH1D* hM_bare_bg_muRand_HLT_dR2 = (TH1D*)hM1_bare_bg_muRand_HLT_dR2->Clone();
    hM_bare_bg_muRand_HLT_dR2->Add(hM2_bare_bg_muRand_HLT_dR2);
    normaliseHist(hM_bare_bg_muRand_HLT_dR2); // as we added m1 + m2
    setMassAUTitles(hM_bare_bg_muRand_HLT_dR2);
    doBGHist(hM_bare_bg_muRand_HLT_dR2);
    setAltTitleLabelSizes(hM_bare_bg_muRand_HLT_dR2);
    hM_bare_bg_muRand_HLT_dR2->Draw("HISTE");
    t_signal->Draw();
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
    setAltTitleLabelSizes(hM_bg_dR2);
    hM_bg_dR2->Draw("HISTE");
    t_signal->Draw();
    c1->SaveAs("Combined/M_10bins_bare_bg_both_muRand_HLT_dR2.pdf");

    // signal
    TH1D* hM1_bare_sig_muRand_HLT_dR2 = (TH1D*) f_sig_mass2->Get("hM1");
    TH1D* hM2_bare_sig_muRand_HLT_dR2 = (TH1D*) f_sig_mass2->Get("hM2");
    TH1D* hM_bare_sig_muRand_HLT_dR2 = (TH1D*)hM1_bare_sig_muRand_HLT_dR2->Clone();
    hM_bare_sig_muRand_HLT_dR2->Add(hM2_bare_sig_muRand_HLT_dR2);
    normaliseHist(hM_bare_sig_muRand_HLT_dR2);
    setMassAUTitles(hM_bare_sig_muRand_HLT_dR2);
    doSignalHist(hM_bare_sig_muRand_HLT_dR2);
    setAltTitleLabelSizes(hM_bare_sig_muRand_HLT_dR2);
    hM_bare_sig_muRand_HLT_dR2->Draw("HISTE");
    t_signal->Draw();
    c1->SaveAs("Combined/M_10bins_bare_sig_muRand_HLT_dR2.pdf");


    // bbbar + signal
    THStack* st = new THStack("st","");
    st->Add(hM_bare_bg_muRand_HLT_dR2);
    st->Add(hM_bare_sig_muRand_HLT_dR2);
    TLegend* l = new TLegend(0.67, 0.67, 0.88, 0.88);
    l->AddEntry(hM_bare_bg_muRand_HLT_dR2, "QCD b#bar{b} MC", "lp");
    l->AddEntry(hM_bare_sig_muRand_HLT_dR2, "#splitline{Signal MC}{m_{#phi} = 8 GeV}", "lp");
    doStandardLegend(l);
    st->Draw("NOSTACK HISTE");
    setMassAUTitles(st->GetHistogram());
    setAltTitleLabelSizes(st->GetHistogram());    
    l->Draw();
    t_signal->Draw();

    c1->SaveAs("Combined/M_10bins_bare_both_muRand_HLT_dR2.pdf");

    // bbbar + scatter + signal
    THStack* st_all = new THStack("st_all","");
    st_all->Add(hM_bg_dR2);
    st_all->Add(hM_bare_sig_muRand_HLT_dR2);
    TLegend* l_all = new TLegend(0.52,0.62,0.86,0.89);
    l_all->AddEntry(hM_bg_dR2,"Gen. level QCD MC","lp");
    l_all->AddEntry((TObject*)0,"(b#bar{b} + q-g scatter,",""); //null pointers for blank entries
    l_all->AddEntry((TObject*)0,"q = b, #bar{b}, c, #bar{c})","");
    l_all->AddEntry(hM_bare_sig_muRand_HLT_dR2, "Signal MC", "lp");
    l_all->AddEntry((TObject*)0,"m_{#phi} = 8 GeV", "");
    doStandardLegend(l_all);
    st_all->Draw("NOSTACKHISTE");
    setMassAUTitles(st_all->GetHistogram());
    setAltTitleLabelSizes(st_all->GetHistogram());    
    l_all->Draw();
    t_signal->Draw();
    c1->SaveAs("Combined/M_10bins_bare_allMC_muRand_HLT_dR2.pdf");

    
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
    setAltTitleLabelSizes(hM_side_bg_muRand_HLT_dR1);
    hM_side_bg_muRand_HLT_dR1->Draw("HISTE");
    t_side->Draw();
    c1->SaveAs("Combined/M_10bins_side_bg_muRand_HLT_dR1.pdf");

    // // bbbar + scatter
    TH1D* hM1_side_bg_dR1 = combineRebin10bins(filesdR1, "hM1_side_1to2p5_unnormalised", scalingFactors);
    TH1D* hM2_side_bg_dR1 = combineRebin10bins(filesdR1, "hM1_side_1to2p5_unnormalised", scalingFactors);
    TH1D* hM_side_bg_dR1 = (TH1D*) hM1_side_bg_dR1->Clone();
    hM_side_bg_dR1->Add(hM2_side_bg_dR1);
    normaliseHist(hM_side_bg_dR1);
    setMassAUTitles(hM_side_bg_dR1);
    doBGHist(hM_side_bg_dR1);
    setAltTitleLabelSizes(hM_side_bg_dR1);
    hM_side_bg_dR1->Draw("HISTE");
    t_side->Draw();
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
    setAltTitleLabelSizes(hM_side_sig_muRand_HLT_dR1);
    hM_side_sig_muRand_HLT_dR1->Draw("HISTE");
    t_side->Draw();
    c1->SaveAs("Combined/M_10bins_side_sig_muRand_HLT_dR1.pdf");

    // bbar + signal
    THStack* st_side = new THStack("st_side","");
    st_side->Add(hM_side_bg_muRand_HLT_dR1);
    st_side->Add(hM_side_sig_muRand_HLT_dR1);
    st_side->Draw("NOSTACK HISTE");
    setMassAUTitles(st_side->GetHistogram());
    setAltTitleLabelSizes(st_side->GetHistogram());    
    l->Draw();
    t_side->Draw();
    c1->SaveAs("Combined/M_10bins_side_both_muRand_HLT_dR1.pdf");
    
    // bbar + scatter + signal
    THStack* st_side_all = new THStack("st_side_all","");
    st_side_all->Add(hM_side_bg_dR1);
    st_side_all->Add(hM_side_sig_muRand_HLT_dR1);
    st_side_all->Draw("NOSTACK HISTE");
    setMassAUTitles(st_side_all->GetHistogram());
    setAltTitleLabelSizes(st_side_all->GetHistogram());    
    l_all->Draw();
    t_side->Draw();
    c1->SaveAs("Combined/M_10bins_side_allMC_muRand_HLT_dR1.pdf");


    ///////////////////////
    // CORRELATION PLOTS //
    ///////////////////////
    //signal mc
    TH1D* hCorr_bare_sig = (TH1D*) f_sig_mass2->Get("hCorr1D");
    TH1D* hCorr_side_sig = (TH1D*) f_sig_mass1->Get("hCorr1D_side_1to2p5");
    //bbbar
    TH1D* hCorr_bare_b = (TH1D*) f_bg_mass2->Get("hCorr1D");
    TH1D* hCorr_side_b = (TH1D*) f_bg_mass1->Get("hCorr1D_side_1to2p5");
    
    TH1D* hM1D_bare_b_unnormalised = (TH1D*) f_bg_mass2->Get("hM_unnormalised"); // signal region
    TH2D* hM2D_bare_b_unnormalised = (TH2D*) f_bg_mass2->Get("hM1vsM2_unnormalised");
    
    TH1D* hM1D_side_b_unnormalised = (TH1D*) f_bg_mass1->Get("hM_side_1to2p5_unnormalised"); //sideband
    TH2D* hM2D_side_b_unnormalised = (TH2D*) f_bg_mass1->Get("hM1vsM2_side_1to2p5_unnormalised");
    
    // scatter
    TH1D* hCorr_bare_scatter = (TH1D*) f_scatter_mass2->Get("hCorr1D");
    TH1D* hCorr_side_scatter = (TH1D*) f_scatter_mass1->Get("hCorr1D_side_1to2p5");
    
    TH1D* hM1D_bare_scatter_unnormalised = (TH1D*) f_scatter_mass2->Get("hM_unnormalised"); //signal region
    TH2D* hM2D_bare_scatter_unnormalised = (TH2D*) f_scatter_mass2->Get("hM1vsM2_unnormalised");
    
    TH1D* hM1D_side_scatter_unnormalised = (TH1D*) f_scatter_mass1->Get("hM_side_1to2p5_unnormalised"); //sideband
    TH2D* hM2D_side_scatter_unnormalised = (TH2D*) f_scatter_mass1->Get("hM1vsM2_side_1to2p5_unnormalised");
    
    t_signal = new TPaveText(0.15, 0.7, 0.4, 0.8, "NDC");
    t_signal->AddText("Signal region");
    t_signal->SetFillColor(kWhite);
    t_signal->SetBorderSize(0);

    t_side = new TPaveText(0.15, 0.7, 0.4, 0.8, "NDC");
    t_side->AddText("Control region B");
    t_side->SetFillColor(kWhite);
    t_side->SetBorderSize(0);

    // set title sizes, titles, max, min, etc
    setTitleLabelSizes(hCorr_bare_sig);
    doSignalHist(hCorr_bare_sig);
    setCorrTitlesMaxMin(hCorr_bare_sig);
    
    setTitleLabelSizes(hCorr_side_sig);
    doSignalHist(hCorr_side_sig);
    setCorrTitlesMaxMin(hCorr_side_sig);

    setTitleLabelSizes(hCorr_bare_b);
    doBGHist(hCorr_bare_b);
    setCorrTitlesMaxMin(hCorr_bare_b);

    setTitleLabelSizes(hCorr_side_b);
    doBGHist(hCorr_side_b);
    setCorrTitlesMaxMin(hCorr_side_b);

    setTitleLabelSizes(hCorr_bare_scatter);
    doAltBGHist(hCorr_bare_scatter);
    setCorrTitlesMaxMin(hCorr_bare_scatter);
    
    setTitleLabelSizes(hCorr_side_scatter);
    doAltBGHist(hCorr_side_scatter);
    setCorrTitlesMaxMin(hCorr_side_scatter);
    
    TLine line(0,1,10,1);
    line.SetLineStyle(2);
    line.SetLineColor(12);
    line.SetLineWidth(2);

    c1->SetTicks(1,1);
    
    ///////////////////
    // signal region //
    ///////////////////

    // signal MC by itself
    hCorr_bare_sig->Draw();
    line.Draw();
    t_signal->Draw();
    c1->SaveAs("Combined/Corr_bare_sig.pdf");
    
    // bbbar qcd by itself
    hCorr_bare_b->Draw();
    line.Draw();
    t_signal->Draw();
    c1->SaveAs("Combined/Corr_bare_b.pdf");
    
    // bbbar + signal
    hCorr_bare_sig->Draw("SAME");
    l->Draw();
    c1->SaveAs("Combined/Corr_bare_b_sig.pdf");


    // make combined corr plot for both qcd
    std::vector<TH1D*> bare_bg_hists1D;
    bare_bg_hists1D.push_back(hM1D_bare_b_unnormalised);
    bare_bg_hists1D.push_back(hM1D_bare_scatter_unnormalised);
    TH1D* hM1D_bare_bg = combineScale(bare_bg_hists1D, scalingFactors);
    
    std::vector<TH2D*> bare_bg_hists2D;
    bare_bg_hists2D.push_back(hM2D_bare_b_unnormalised);
    bare_bg_hists2D.push_back(hM2D_bare_scatter_unnormalised);
    TH2D* hM2D_bare_bg = combineScale(bare_bg_hists2D, scalingFactors);
    TH1D* hCorr_bare_bg = createCorrelationPlot(hM2D_bare_bg, hM1D_bare_bg);

    // both QCD summed together
    doBGHist(hCorr_bare_bg);
    setTitleLabelSizes(hCorr_bare_bg);
    setCorrTitlesMaxMin(hCorr_bare_bg);
    hCorr_bare_bg->Draw();
    line.Draw();
    t_signal->Draw();
    c1->SaveAs("Combined/Corr_bare_bg.pdf");
    
    // all qcd + signal mc
    hCorr_bare_sig->Draw("SAME");
    l_all->Draw();
    c1->SaveAs("Combined/Corr_bare_bg_sig.pdf");

    //////////////////////
    // control region B //
    //////////////////////

    // signal MC by itself
    hCorr_side_sig->Draw();
    line.Draw();
    c1->SaveAs("Combined/Corr_side_sig.pdf");
    
    // bbbar qcd by itself
    hCorr_side_b->Draw();
    line.Draw();
    c1->SaveAs("Combined/Corr_side_b.pdf");

    // signal + bbbar together
    hCorr_side_sig->Draw();
    hCorr_side_b->Draw("SAME");
    line.Draw();
    l->Draw();
    c1->SaveAs("Combined/Corr_side_b_sig.pdf");

    // make combined corr plot for both qcd
    std::vector<TH1D*> side_bg_hists1D;
    side_bg_hists1D.push_back(hM1D_side_b_unnormalised);
    side_bg_hists1D.push_back(hM1D_side_scatter_unnormalised);
    TH1D* hM1D_side_bg = combineScale(side_bg_hists1D, scalingFactors);
    
    std::vector<TH2D*> side_bg_hists2D;
    side_bg_hists2D.push_back(hM2D_side_b_unnormalised);
    side_bg_hists2D.push_back(hM2D_side_scatter_unnormalised);
    TH2D* hM2D_side_bg = combineScale(side_bg_hists2D, scalingFactors);
    TH1D* hCorr_side_bg = createCorrelationPlot(hM2D_side_bg, hM1D_side_bg);

    // both QCD summed together
    doBGHist(hCorr_side_bg);
    setTitleLabelSizes(hCorr_side_bg);
    setCorrTitlesMaxMin(hCorr_side_bg);
    hCorr_side_bg->Draw();
    line.Draw();
    t_side->Draw();
    c1->SaveAs("Combined/Corr_side_bg.pdf");
    
    // signal + both QCD
    hCorr_side_sig->Draw("SAME");
    l_all->Draw();
    c1->SaveAs("Combined/Corr_side_all_MC.pdf");

    // QCD + data
    hCorr_side_bg->Draw();
    TH1D* hCorr_side_data = dataCorr_side();
    doDataHist(hCorr_side_data);
    setCorrTitlesMaxMin(hCorr_side_data);
    hCorr_side_data->Draw("SAME");
    line.Draw();
    
    TLegend* l_all_data = new TLegend(0.52,0.65,0.86,0.89);
    l_all_data->AddEntry(hM_bg_dR2,"Gen. level QCD MC","lp");
    l_all_data->AddEntry((TObject*)0,"(b#bar{b} + q-g scatter,",""); //null pointers for blank entries
    l_all_data->AddEntry((TObject*)0,"q = b, #bar{b}, c, #bar{c})","");
    l_all_data->AddEntry(hCorr_side_data, "Data", "lp");
    doStandardLegend(l_all_data);
    l_all_data->Draw();
    
    t_side->Draw();
    c1->SaveAs("Combined/Corr_side_bg_data.pdf");

    // bbar + data
    hCorr_side_b->Draw();
    hCorr_side_data->Draw("SAME");
    TLegend* l_b_data = new TLegend(0.52,0.65,0.86,0.89);
    l_b_data->AddEntry(hM_bg_dR2,"QCD b#bar{b} MC","lp");
    l_b_data->AddEntry(hCorr_side_data, "Data", "lp");
    doStandardLegend(l_b_data);
    l_b_data->Draw();
    line.Draw();
    t_side->Draw();
    c1->SaveAs("Combined/Corr_side_b_data.pdf");


    /////////////////////////
    // TRACK DISTRIBUTIONS //
    /////////////////////////
    combineHists(f_sig_main2, f_bg_main2, f_scatter_main2, "hNTracks1", "HISTE", "Combined/combined_NTrack1_muRand.pdf", scalingFactors, "Tracks with p_{T} > 2.5 GeV");
    combineHists(f_sig_main2, f_bg_main2, f_scatter_main2, "hNTracksAbs1", "HISTE", "Combined/combined_NTrackAbs1_muRand.pdf", scalingFactors, "Tracks with p_{T} > 2.5 GeV", "Average number of tracks per #mu_{1} / bin");
    combineHists(f_sig_main2, f_bg_main2, f_scatter_main2, "hNTracksAll1", "HISTE", "Combined/combined_NTrackAll1_muRand.pdf", scalingFactors, "Tracks with p_{T} > 1 GeV");
    combineHists(f_sig_main2, f_bg_main2, f_scatter_main2, "hNTracksAllAbs1", "HISTE", "Combined/combined_NTrackAllAbs1_muRand.pdf", scalingFactors, "Tracks with p_{T} > 1 GeV", "Average number of tracks per #mu_{1} / bin");
    combineHists(f_sig_main2, f_bg_main2, f_scatter_main2, "hNSoftTracks1", "HISTE", "Combined/combined_NSoftTrack1_muRand.pdf", scalingFactors, "Tracks with 2.5 > p_{T} > 1 GeV");
    combineHists(f_sig_main2, f_bg_main2, f_scatter_main2, "hNSoftTracksAbs1", "HISTE", "Combined/combined_NSoftTrackAbs1_muRand.pdf", scalingFactors, "Tracks with 2.5 > p_{T} > 1 GeV", "Average number of tracks per #mu_{1} / bin");


    ////////////////////////////////////////////
    // mass shape as fn of # tracks about mu2 //
    ////////////////////////////////////////////
    THStack* st_Ntk2_234 = new THStack("","");
    TH1D* histM1_Ntk2_2 = (TH1D*) f_bg_mass2->Get("hM1_Ntk2_2");
    TH1D* histM1_Ntk2_3 = (TH1D*) f_bg_mass2->Get("hM1_Ntk2_3");
    TH1D* histM1_Ntk2_4 = (TH1D*) f_bg_mass2->Get("hM1_Ntk2_4");
    doStandardHist(hM1_bare_bg_muRand_HLT_dR2);
    hM1_bare_bg_muRand_HLT_dR2->SetMarkerStyle(1);
    doCustomHist(histM1_Ntk2_2, kRed);
    doCustomHist(histM1_Ntk2_3, kBlack);
    doCustomHist(histM1_Ntk2_4, kGreen+3);
    st_Ntk2_234->Add(hM1_bare_bg_muRand_HLT_dR2);
    st_Ntk2_234->Add(histM1_Ntk2_2);
    st_Ntk2_234->Add(histM1_Ntk2_3);
    st_Ntk2_234->Add(histM1_Ntk2_4);
    st_Ntk2_234->Draw("NOSTACK E");
    setMassAUTitles(st_Ntk2_234->GetHistogram());
    st_Ntk2_234->GetHistogram()->SetXTitle("m_{1}(#mu-tk) [GeV]");
    setAltTitleLabelSizes(st_Ntk2_234->GetHistogram());
    TLegend* l_Ntk2_234 = new TLegend(0.56, 0.6, 0.88, 0.88);
    l_Ntk2_234->AddEntry((TObject*)0, "Gen. level QCD b#bar{b} MC", "");
    l_Ntk2_234->AddEntry(hM1_bare_bg_muRand_HLT_dR2, "N_{tk,2} = 1", "l");
    l_Ntk2_234->AddEntry(histM1_Ntk2_2, "N_{tk,2} = 2", "l");
    l_Ntk2_234->AddEntry(histM1_Ntk2_3, "N_{tk,2} = 3", "l");
    l_Ntk2_234->AddEntry(histM1_Ntk2_4, "N_{tk,2} = 4", "l");
    doStandardLegend(l_Ntk2_234);
    l_Ntk2_234->Draw();
    c1->SaveAs("Combined/M1_Ntk2.pdf");

    THStack* st_Ntk2_2or3 = new THStack("","");
    TH1D* histM1_Ntk2_2or3 = (TH1D*) f_bg_mass2->Get("hM1_Ntk2_2or3");
    doCustomHist(histM1_Ntk2_2or3, kRed);
    st_Ntk2_2or3->Add(hM1_bare_bg_muRand_HLT_dR2);
    st_Ntk2_2or3->Add(histM1_Ntk2_2or3);
    st_Ntk2_2or3->Draw("NOSTACK E");
    setMassAUTitles(st_Ntk2_2or3->GetHistogram());
    st_Ntk2_2or3->GetHistogram()->SetXTitle("m_{1}(#mu-tk) [GeV]");
    setAltTitleLabelSizes(st_Ntk2_2or3->GetHistogram());
    TLegend* l_Ntk2_2or3 = new TLegend(0.56, 0.6, 0.88, 0.88);
    l_Ntk2_2or3->AddEntry((TObject*)0, "Gen. level QCD b#bar{b} MC", "");
    l_Ntk2_2or3->AddEntry(hM1_bare_bg_muRand_HLT_dR2, "N_{tk,2} = 1", "l");
    l_Ntk2_2or3->AddEntry(histM1_Ntk2_2or3, "N_{tk,2} = 2, 3", "l");
    doStandardLegend(l_Ntk2_2or3);
    l_Ntk2_2or3->Draw();
    c1->SaveAs("Combined/M1_Ntk2_2or3.pdf");


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