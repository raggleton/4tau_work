void printBinContents(TH1D* h) {
    TAxis *x_axis = h->GetXaxis();
    for(int i = 1; i <= h->GetNbinsX(); ++i) {
        std::cout << "Bin " << i << " [ " << x_axis->GetBinLabel(i) << " ] " << ": " << h->GetBinContent(i) << " \\pm " << h->GetBinError(i) << std::endl;
    }
}


/**
 * This script redoes the 1D correlation coefficents, making sure that the denominator is
 * calulated with M(i) x M(j), where M = M_1 + M_2,
 * and we ignore any errors on the denominator as correlated with numerator.
 */
int new_corr() {
    TH1::SetDefaultSumw2();
    gStyle->SetOptStat("");

    double massbins[5] = {0,1,2,3,10};
    int nBins = 4;
    // TFile * fqcdb = TFile::Open("QCDb/output_bare_bg_muRand_HLT_dR2.root", "READ");
    TFile * fqcdb = TFile::Open("QCDb/output_bare_bg_muRand_HLT_IPopt2_dR2.root", "READ");
    TH1D* hM1_orig = (TH1D*) fqcdb->Get("hM1");
    TH1D* hM1 = hM1_orig->Rebin(4,"hM1", massbins);
    TH1D* hM2_orig = (TH1D*) fqcdb->Get("hM2");
    TH1D* hM2 = hM2_orig->Rebin(4,"hM2", massbins);
    TH1D* hM = (TH1D*) fqcdb->Get("hM");

    TH2D* hM1vsM2_orig = (TH2D*) fqcdb->Get("hM1vsM2");
    TH2D* hM1vsM2 = hM1vsM2_orig->Clone("hM1vsM2");
    hM1vsM2->SetTitle("");

    TH2D* hMxM = new TH2D("hMxM", "", 4,massbins,4, massbins);
    for (int i = 1; i <= nBins; ++i) {
        for (int j = 1; j <= nBins; ++j) {
            hMxM->SetBinContent(i,j, hM->GetBinContent(i)* hM->GetBinContent(j));
            hMxM->SetBinError(i,j,0);
        }
    }
    hM1vsM2->Scale(1./hM1vsM2->Integral());
    hMxM->Scale(1./hMxM->Integral());
    hMxM->Draw("TEXTE");
    hM1vsM2->Draw("TEXTE");
    hM1vsM2->Divide(hMxM);

    TH1D* h1D = new TH1D("h1D", "", 10, 0, 10);
    int counter = 1;

    for (int i = 1; i <= nBins; i++) {
        for (int j = i; j <= nBins; j++) {
            h1D->SetBinContent(counter, hM1vsM2->GetBinContent(i,j));
            h1D->SetBinError(counter, hM1vsM2->GetBinError(i,j));
            TString label;
            label.Form("(%d,%d)", i, j);
            h1D->GetXaxis()->SetBinLabel(counter, label);
            counter++;
        }
    }
    TCanvas *c1 = new TCanvas();
    c1->SetTicks(1,1);

    TLine *line = new TLine(0,1,10,1);
    line->SetLineStyle(2);
    line->SetLineColor(12);
    line->SetLineWidth(2);

    h1D->Draw("E1");
    h1D->GetXaxis()->SetLabelSize(0.06);
    h1D->GetXaxis()->SetTitle("Bin");
    h1D->GetXaxis()->SetTitleSize(0.04);
    h1D->GetYaxis()->SetRangeUser(0,2.0);
    h1D->GetYaxis()->SetTitle("Correlation coefficient");
    h1D->GetYaxis()->SetTitleSize(0.04);
    h1D->SetLineWidth(2);
    h1D->Draw("E1");
    line->Draw();

    TPaveText* t_signal = new TPaveText(0.6, 0.72, 0.85, 0.83, "NDC");
    t_signal->AddText("Signal region");
    t_signal->SetFillColor(kWhite);
    t_signal->SetBorderSize(0);

    t_signal->Draw();

    c1->SaveAs("corr1D_corrected_IPopt2.pdf");

    printBinContents(h1D);
    std::cout << std::endl;

    TH1D* h1D_OLD = (TH1D*) fqcdb->Get("hCorr1D");
    printBinContents(h1D_OLD);

}
