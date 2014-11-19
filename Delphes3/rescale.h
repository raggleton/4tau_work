#ifndef RESCALE_H
#define RESCALE_H

#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <iostream>
#include <exception>

#include "commonFunctions.h"

using std::cout;
using std::endl;

/**
 * @brief Class to rescale track eta and phi
 * @details [long description]
 *
 */
class Rescaler {
public:
    Rescaler();
    void rescaleTrack(Track* tk, Track* mu);
private:
    bool doRescaling;
    // Add files for rescaling
    TFile * fileData;
    TFile * fileMC;
    TFile * weights;
    TH1 * histMC;
    TH1 * histData;
};

Rescaler::Rescaler(): doRescaling(true)
{
    TH1::SetDefaultSumw2();

    fileData = TFile::Open("rescaling_files/Data_kine.root", "READ");
    fileMC   = TFile::Open("rescaling_files/output_reweight_bare_bg_muRand_HLT_dR2.root", "READ");
    weights  = TFile::Open("rescaling_files/weights.root", "READ");
    if (!fileData || !fileMC || !weights) {
        doRescaling = false;
        cout << "Couldn't open rescaling files, no rescaling for you." << endl;
    } else {
        // get hists from MC and data
        // have to combine hard + soft
        TH1* histData_hard = (TH1*) fileData->Get("dRHardmuonTrk_InvQH");
        histData_hard->Sumw2();
        TH1* histData_soft = (TH1*) fileData->Get("dRSoftmuonTrk_InvQH");
        histData_soft->Sumw2();
        histData = (TH1*) histData_hard->Clone("dRmuonTrk_InvQH");
        histData->Add(histData_soft);
        TH1* histMC_hard = (TH1*) fileMC->Get("hDRmutk_hard");
        TH1* histMC_soft = (TH1*) fileMC->Get("hDRmutk_soft");
        histMC = (TH1*) histMC_hard->Clone("hDRmutk");
        histMC->Add(histMC_soft);

        if (!histMC || ! histData) {
            doRescaling = false;
            cout << "Couldn't get rescaling histograms, no rescaling for you." << endl;
        } else {
            normaliseHist(histMC);
            normaliseHist(histData);
        }
    }

}

void Rescaler::rescaleTrack(Track* tk, Track* mu) {
    if (!doRescaling) return;

    double deltaR = (tk->P4()).DeltaR(mu->P4());
    double deltaEta = tk->Eta - mu->Eta;
    double deltaPhi = (tk->P4()).DeltaPhi(mu->P4());

    // Rescale based on quantiles from hsits
    int bin = histMC->FindBin(deltaR);
    float contentBin = histMC->GetBinContent(bin);
    float binWidth = histMC->GetBinWidth(bin);
    float lowEdge = histMC->GetBinLowEdge(bin);
    float integral = contentBin*(deltaR-lowEdge)/binWidth;

    if (bin>1) {
        for (int i=1; i<bin; ++i) {
           integral += histMC->GetBinContent(i);
        }
    }

    double quantile[1];
    double prob[1];
    prob[0] = integral;

    histData->GetQuantiles(1,quantile,prob);

    // Get new deltaR(mu-tk) and deltaEta and deltaPhi to be applied to the Track
    float deltaR_new = float(quantile[0]);
    cout << "deltaR: " << deltaR << " deltaR_new: " << deltaR_new << endl;
    float deltaEta_new = deltaEta * deltaR_new / deltaR;
    float deltaPhi_new = deltaPhi * deltaR_new / deltaR;
    cout << mu->Eta << " : " << mu->Phi << endl;
    cout << tk->Eta << " : " << tk->Phi << endl;
    tk->Eta = mu->Eta + deltaEta_new;
    tk->Phi = mu->Phi + deltaPhi_new;
    cout << tk->Eta << " : " << tk->Phi << endl;
}

#endif
