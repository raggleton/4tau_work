#include "commonFunctions.h"
#include "classes/SortableObject.h"
#include "cuts.h"
#include "rescale.h"
// #include "tdrstyle.C"

using std::cout;
using std::endl;


/**
 * NB ALL DISTANCES IN MM, ALL ENERGIES IN GeV
 * SEE DelphesHepMCReader.cc
 */

// int getMassBin(double m){
// 	std::vector<double> massBins {0,1,2,3,10};
// 	for(unsigned i = 1; i < massBins.size(); i++){
// 		if (m < massBins[i])
// 			return i;
// 	}
// 	return 5;
// }

// /**
//  * These checks to see if muA and muB satisfy all muon condtions
//  * apart from pT.
//  * @param  muA Higher pT muon
//  * @param  muB Lesser pT muon
//  * @param  deltaR dR(mu-mu) cut
//  * @return    TRUE if muA and muB pass cuts, FALSE otherwise
//  */
// bool checkMuons(Track* muA, Track* muB, double deltaR){
// 	if ((muA->Charge == muB->Charge)
// 		&& (fabs(muA->Eta) < 2.1)
// 		&& (fabs(muB->Eta) < 2.1)
// 		&& (fabs(muA->Zd) < 1.) // dZ < 0.1cm
// 		&& (fabs(muB->Zd) < 1.) // dZ < 0.1cm
// 		&& (fabs(muA->Dxy) < 0.3 ) // d0 < 0.03cm
// 		&& (fabs(muB->Dxy) < 0.3 ) // d0 < 0.03cm
// 		&& ((muA->P4().DeltaR(muB->P4())) > deltaR)
// 		){
// 		return true;
// 	} else {
// 		return false;
// 	}
// }
// /**
//  * These checks to see if muA and muB satisfy all muon condtions
//  * apart from pT and impact params.
//  * (Old version kept incase you use RawGenMuons branch instead of GenMuons)
//  * @param  muA Higher pT muon
//  * @param  muB Lesser pT muon
//  * @param  deltaR dR(mu-mu) cut
//  * @return    TRUE if muA and muB pass cuts, FALSE otherwise
//  */
// bool checkMuons(GenParticle* muA, GenParticle* muB, double deltaR){
// 	if ((muA->Charge == muB->Charge)
// 		&& (fabs(muA->Eta) < 2.1)
// 		&& (fabs(muB->Eta) < 2.1)
// 		&& ((muA->P4().DeltaR(muB->P4())) > deltaR)
// 		){
// 		return true;
// 	} else {
// 		return false;
// 	}
// }

int countParticle(TClonesArray* branchAll, int PID) {
	int nParticle = 0;
	// cout << "**** EVENT *****" << endl;
	for (int a = 0; a < branchAll->GetEntries(); a++){
		GenParticle* candP = (GenParticle*) branchAll->At(a);
		// cout << candP->PID << endl;
		if (fabs(candP->PID) == fabs(PID)) {
			nParticle++;
		}
	}
	// cout << "**** ENDE *****" << endl;
	return nParticle;
}

/**
 * makes clone fo hist, appends "suffix" to hist name,
 * write to currently open file
 */
template <typename T>
void makeCopySave(T* h, std::string suffix="_unnormalised") {
	T* h_clone = (T*) h->Clone((h->GetName()+suffix).c_str());
	h_clone->Write("", TObject::kOverwrite);
}


/**
 * @brief For adding together mass plots and rebinning
 * @details [long description]
 *
 * @param histM1 [description]
 * @param histM2 [description]
 * @param bins [description]
 * @param name [description]
 * @return [description]
 */
template<typename T>
T* makeSumRebin(T * histM1, T * histM2, std::vector<double> bins, std::string name) {
	TObject* old=gDirectory->GetList()->FindObject("Temp");
	if (old)
		gDirectory->GetList()->Remove(old);
	T * histM = new T("Temp", "Inv. Mass of system, signal selection; m(#mu-tk) [GeV];A.U.",
		histM1->GetNbinsX(),
		histM1->GetXaxis()->GetBinLowEdge(1),
		histM1->GetXaxis()->GetBinLowEdge(histM1->GetNbinsX()+1));
	histM->Add(histM1);
	histM->Add(histM2);
	return (T*) histM->Rebin(bins.size()-1, name.c_str(), &bins[0]);
}


void massPlots(int argc, char* argv[])
{
	// setTDRStyle();

	TH1::SetDefaultSumw2();

	gSystem->Load("libDelphes");

	// Get program options from user and store
	ProgramOpts pOpts(argc, argv);

	MCsource source      = pOpts.getSource(); // get MC source (signal, qcdb, qcdc)
	bool doSignal        = pOpts.getSignal(); // for signal or not
	// bool doMu         = pOpts.getQCDMu(); // for QCDb - either inclusive decays or mu only decays
	bool swapMuRandomly  = pOpts.getMuOrdering(); // if true, fills plots for mu 1 and 2 randomly from highest & 2nd highest pt muons. Otherwise, does 1 = leading (highest pt), 2 = subleading (2nd highest pt)
	bool doHLT           = pOpts.getHLT(); // whether to use MC that has HLT cuts already applied or not.
	// bool DEBUG        = pOpts.getVerbose(); // output debug statments
	double deltaR        = pOpts.getdR(); // dR(mu-mu) value to use
	double rescaleFactor = pOpts.getRescale(); // factor to rescaleFactor track eta & phi to match data. 1 = no rescaling
	bool doRescale       = pOpts.doRescaling(); // factor to rescaleFactor track eta & phi to match data. 1 = no rescaling
	bool do1to1p5        = false; // for additional sideband studies. Slower?

	// for rescaling signal tracks
	Rescaler r;

	// Create chain of root trees
	TChain chain("Delphes");
	addInputFiles(&chain, &pOpts);
	pOpts.printProgramOptions();

	// Create object of class ExRootTreeReader
	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);

	// Get pointers to branches used in this analysis
	// Use the data_flow.png and tcl file to figure out what branches are available, and what class they are
	// and use https://cp3.irmp.ucl.ac.be/projects/delphes/wiki/WorkBook/RootTreeDescription
	// TClonesArray* branchMuon = treeReader->UseBranch("Muon");
	TClonesArray* branchTracks   = treeReader->UseBranch("Track");
	// TClonesArray* branchGenMuons = treeReader->UseBranch("OnlyGenMuons"); // GenParticle object
	TClonesArray* branchGenMuons = treeReader->UseBranch("GenMuon"); // Track object
	// TClonesArray* branchStable   = treeReader->UseBranch("StableParticle");
	TClonesArray* branchAll      = treeReader->UseBranch("AllParticle");

	//////////////////////
	// Book histograms  //
	//////////////////////

	// mass bins for correlation plots
	std::vector<double> massBins {0,1,2,3,10};
	int nBinsX = massBins.size()-1;

	// ------------------------
	// m1 & m2 1D distributions
	// ------------------------

	int nMassBins = 10;
	double massMin = 0.;
	double massMax = 10.;

	TH1D* histM1                                = new TH1D("hM1", "Inv. Mass of 1st system, full selection; m(#mu_{1}-tk) [GeV];A.U.", nMassBins, massMin, massMax);
	TH1D* histM2                                = new TH1D("hM2", "Inv. Mass of 2st system, full selection; m(#mu_{2}-tk) [GeV];A.U.", nMassBins, massMin, massMax);

	TH1D* histM1_ModMass                        = new TH1D("hM1_ModMass", "Inv. Mass of 1st system, full selection; m(#mu_{1}-tk) [GeV];A.U.", nMassBins, massMin, massMax);
	TH1D* histM2_ModMass                        = new TH1D("hM2_ModMass", "Inv. Mass of 2st system, full selection; m(#mu_{2}-tk) [GeV];A.U.", nMassBins, massMin, massMax);

	TH1D* histM1_regionA                        = new TH1D("hM1_regionA", "Inv. Mass of 1st system, control region A; m(#mu_{1}-tk) [GeV];A.U.", nMassBins, massMin, massMax);
	TH1D* histM2_regionA                        = new TH1D("hM2_regionA", "Inv. Mass of 2st system, control region A; m(#mu_{2}-tk) [GeV];A.U.", nMassBins, massMin, massMax);

	TH1D* histM1_regionA_ModMass                = new TH1D("hM1_regionA_ModMass", "Inv. Mass of 1st system, control region A; m(#mu_{1}-tk) [GeV];A.U.", nMassBins, massMin, massMax);
	TH1D* histM2_regionA_ModMass                = new TH1D("hM2_regionA_ModMass", "Inv. Mass of 2st system, control region A; m(#mu_{2}-tk) [GeV];A.U.", nMassBins, massMin, massMax);

	TH1D* histM1_regionB                        = new TH1D("hM1_regionB", "Inv. Mass of 1st system, control region B; m(#mu_{1}-tk) [GeV];A.U.", nMassBins, massMin, massMax);
	TH1D* histM2_regionB                        = new TH1D("hM2_regionB", "Inv. Mass of 2st system, control region B; m(#mu_{2}-tk) [GeV];A.U.", nMassBins, massMin, massMax);

	TH1D* histM1_regionB_ModMass                = new TH1D("hM1_regionB_ModMass", "Inv. Mass of 1st system, control region B; m(#mu_{1}-tk) [GeV];A.U.", nMassBins, massMin, massMax);
	TH1D* histM2_regionB_ModMass                = new TH1D("hM2_regionB_ModMass", "Inv. Mass of 2st system, control region B; m(#mu_{2}-tk) [GeV];A.U.", nMassBins, massMin, massMax);

	// for diff # tracks about 2nd muon
	TH1D* histM1_Ntk2_2                         = new TH1D("hM1_Ntk2_2", "Inv. Mass of 1st system, full selection with N_{tk,2} = 2; m(#mu_{1}-tk) [GeV];A.U.", nMassBins, massMin, massMax);
	TH1D* histM1_Ntk2_3                         = new TH1D("hM1_Ntk2_3", "Inv. Mass of 1st system, full selection with N_{tk,2} = 3; m(#mu_{1}-tk) [GeV];A.U.", nMassBins, massMin, massMax);
	TH1D* histM1_Ntk2_4                         = new TH1D("hM1_Ntk2_4", "Inv. Mass of 1st system, full selection with N_{tk,2} = 4; m(#mu_{1}-tk) [GeV];A.U.", nMassBins, massMin, massMax);
	TH1D* histM1_Ntk2_2or3                      = new TH1D("hM1_Ntk2_2or3", "Inv. Mass of 1st system, full selection with N_{tk,2} = 2 or 3; m(#mu_{1}-tk) [GeV];A.U.", nMassBins, massMin, massMax);

	// TH1D* histM1_Ntk2_2_loosetau             = new TH1D("hM1_Ntk2_2_loosietau", "Inv. Mass of 1st system, full selection with N_{tk,2} = 2; m(#mu_{1}-tk) [GeV];A.U.", nMassBins, massMin, massMax);
	// TH1D* histM1_Ntk2_3_loosetau             = new TH1D("hM1_Ntk2_3_loosietau", "Inv. Mass of 1st system, full selection with N_{tk,2} = 3; m(#mu_{1}-tk) [GeV];A.U.", nMassBins, massMin, massMax);
	// TH1D* histM1_Ntk2_4_loosetau             = new TH1D("hM1_Ntk2_4_loosietau", "Inv. Mass of 1st system, full selection with N_{tk,2} = 4; m(#mu_{1}-tk) [GeV];A.U.", nMassBins, massMin, massMax);
	// TH1D* histM1_Ntk2_2or3_loosetau          = new TH1D("hM1_Ntk2_2or3_loosietau", "Inv. Mass of 1st system, full selection with N_{tk,2} = 2 or 3; m(#mu_{1}-tk) [GeV];A.U.", nMassBins, massMin, massMax);

	// TH1D* histM1_fine                        = new TH1D("hM1Fine", "Inv. Mass of 1st system, full selection; m(#mu_{1}-tk) [GeV];A.U.",50,0,10);
	// TH1D* histM2_fine                        = new TH1D("hM2Fine", "Inv. Mass of 2st system, full selection; m(#mu_{2}-tk) [GeV];A.U.",50,0,10);

	// TH1D* histM1_fine_regionA                = new TH1D("hM1_fine_regionA", "Inv. Mass of 1st system, control region A; m(#mu_{1}-tk) [GeV];A.U.",50,0,10);
	// TH1D* histM2_fine_regionA                = new TH1D("hM_fine_regionA", "Inv. Mass of 2st system, control region A; m(#mu_{2}-tk) [GeV];A.U.",50,0,10);

	// TH1D* histM1_fine_regionB                = new TH1D("hM1_fine_regionB", "Inv. Mass of 1st system, control region B; m(#mu_{1}-tk) [GeV];A.U.",50,0,10);
	// TH1D* histM2_fine_regionB                = new TH1D("hM_fine_regionB", "Inv. Mass of 2st system, control region B; m(#mu_{2}-tk) [GeV];A.U.",50,0,10);

	// TH1D* histM1_side_1to2p5                 = new TH1D("hM1_side_1to2p5","m(#mu_{1}-tk) in sideband (soft tk p_{T} = 1 - 2.5 GeV);m(#mu_{1}-tk) [GeV];A.U.",50,0,10);
	// TH1D* histM2_side_1to2p5                 = new TH1D("hM2_side_1to2p5","m(#mu_{1}-tk) in sideband (soft tk p_{T} = 1 - 2.5 GeV);m(#mu_{1}-tk) [GeV];A.U.",50,0,10);

	// TH1D* histM1_side_1to2p5_loosetau        = new TH1D("hM1_side_1to2p5_loosetau","m(#mu_{1}-tk) in sideband (soft tk p_{T} = 1 - 2.5 GeV);m(#mu_{1}-tk) [GeV];A.U.",50,0,10);
	// TH1D* histM2_side_1to2p5_loosetau        = new TH1D("hM2_side_1to2p5_loosetau","m(#mu_{1}-tk) in sideband (soft tk p_{T} = 1 - 2.5 GeV);m(#mu_{1}-tk) [GeV];A.U.",50,0,10);

	// TH1D* histM1_side_1to1p5                 = new TH1D("hM1_side_1to1p5","m(#mu_{1}-tk) in sideband (soft tk p_{T} = 1 - 1.5 GeV);m(#mu_{1}-tk) [GeV];A.U.",50,0,10);
	// TH1D* histM2_side_1to1p5                 = new TH1D("hM2_side_1to1p5","m(#mu_{1}-tk) in sideband (soft tk p_{T} = 1 - 1.5 GeV);m(#mu_{1}-tk) [GeV];A.U.",50,0,10);



	//------------------
	// m1 in bins of m2
	//------------------
	/*
	// MC truth - use actual mu-tk pairs from tau
	TH1D* histM1_truth_0to1                     = new TH1D("hM1_truth_0to1","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 0-1 GeV; m(#mu_{1}-tk) [GeV]; A.U.",nBinsX,&massBins[0]);
	TH1D* histM1_truth_1to2                     = new TH1D("hM1_truth_1to2","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 1-2 GeV; m(#mu_{1}-tk) [GeV]; A.U.",nBinsX,&massBins[0]);
	TH1D* histM1_truth_2to3                     = new TH1D("hM1_truth_2to3","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 2-3 GeV; m(#mu_{1}-tk) [GeV]; A.U.",nBinsX,&massBins[0]);
	TH1D* histM1_truth_3toInf                   = new TH1D("hM1_truth_3toInf","m(#mu_{1}-tk) for m(#mu_{2}-tk) > 3 GeV; m(#mu_{1}-tk) [GeV]; A.U.",nBinsX,&massBins[0]);

	// actual dist using signal selection
	TH1D* histM1_0to1                           = new TH1D("hM1_0to1","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 0-1 GeV; m(#mu_{1}-tk) [GeV]; A.U.",nBinsX,&massBins[0]);
	TH1D* histM1_1to2                           = new TH1D("hM1_1to2","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 1-2 GeV; m(#mu_{1}-tk) [GeV]; A.U.",nBinsX,&massBins[0]);
	TH1D* histM1_2to3                           = new TH1D("hM1_2to3","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 2-3 GeV; m(#mu_{1}-tk) [GeV]; A.U.",nBinsX,&massBins[0]);
	TH1D* histM1_3toInf                         = new TH1D("hM1_3toInf","m(#mu_{1}-tk) for m(#mu_{2}-tk) > 3 GeV; m(#mu_{1}-tk) [GeV]; A.U.",nBinsX,&massBins[0]);

	// actual dist using sideband selection (soft track pT 1-2.5)
	TH1D* histM1_side_1to2p5_0to1               = new TH1D("hM1_side_1to2p5_0to1","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 0-1 GeV (soft tk p_{T} = 1-2.5 GeV); m(#mu_{1}-tk) [GeV]; A.U.",nBinsX,&massBins[0]);
	TH1D* histM1_side_1to2p5_1to2               = new TH1D("hM1_side_1to2p5_1to2","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 1-2 GeV (soft tk p_{T} = 1-2.5 GeV); m(#mu_{1}-tk) [GeV]; A.U.",nBinsX,&massBins[0]);
	TH1D* histM1_side_1to2p5_2to3               = new TH1D("hM1_side_1to2p5_2to3","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 2-3 GeV (soft tk p_{T} = 1-2.5 GeV); m(#mu_{1}-tk) [GeV]; A.U.",nBinsX,&massBins[0]);
	TH1D* histM1_side_1to2p5_3toInf             = new TH1D("hM1_side_1to2p5_3toInf","m(#mu_{1}-tk) for m(#mu_{2}-tk) > 3 GeV (soft tk p_{T} = 1-2.5 GeV); m(#mu_{1}-tk) [GeV]; A.U.",nBinsX,&massBins[0]);

	// actual dist using sideband selection (soft track pT 1-1.5)
	TH1D* histM1_side_1to1p5_0to1               = new TH1D("hM1_side_1to1p5_0to1","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 0-1 GeV (soft tk p_{T} = 1-1.5 GeV); m(#mu_{1}-tk) [GeV]; A.U.",nBinsX,&massBins[0]);
	TH1D* histM1_side_1to1p5_1to2               = new TH1D("hM1_side_1to1p5_1to2","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 1-2 GeV (soft tk p_{T} = 1-1.5 GeV); m(#mu_{1}-tk) [GeV]; A.U.",nBinsX,&massBins[0]);
	TH1D* histM1_side_1to1p5_2to3               = new TH1D("hM1_side_1to1p5_2to3","m(#mu_{1}-tk) for m(#mu_{2}-tk) = 2-3 GeV (soft tk p_{T} = 1-1.5 GeV); m(#mu_{1}-tk) [GeV]; A.U.",nBinsX,&massBins[0]);
	TH1D* histM1_side_1to1p5_3toInf             = new TH1D("hM1_side_1to1p5_3toInf","m(#mu_{1}-tk) for m(#mu_{2}-tk) > 3 GeV (soft tk p_{T} = 1-1.5 GeV); m(#mu_{1}-tk) [GeV]; A.U.",nBinsX,&massBins[0]);

	// --------------------------
	// mu1 pT plots in bins of M2
	// --------------------------


	// mu1 pT plots in bins of m2
	TH1D* histMu1Pt_0to1                        = new TH1D("hMu1Pt_0to1","#mu_{1} p_{T} for m(#mu_{2}-tk) = 0-1 GeV; #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);
	TH1D* histMu1Pt_1to2                        = new TH1D("hMu1Pt_1to2","#mu_{1} p_{T} for m(#mu_{2}-tk) = 1-2 GeV; #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);
	TH1D* histMu1Pt_2to3                        = new TH1D("hMu1Pt_2to3","#mu_{1} p_{T} for m(#mu_{2}-tk) = 2-3 GeV; #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);
	TH1D* histMu1Pt_3toInf                      = new TH1D("hMu1Pt_3toInf","#mu_{1} p_{T} for m(#mu_{2}-tk) = 3-Inf GeV; #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);

	// mu1 pT plots in bins of m2 - MC truth
	TH1D* histMu1Pt_truth_0to1                  = new TH1D("hMu1Pt_truth_0to1","#mu_{1} p_{T} for m(#mu_{2}-tk) = 0-1 GeV; #mu_{1} p_{T} [GeV];A.U.",10,0.,50.);
	TH1D* histMu1Pt_truth_1to2                  = new TH1D("hMu1Pt_truth_1to2","#mu_{1} p_{T} for m(#mu_{2}-tk) = 1-2 GeV; #mu_{1} p_{T} [GeV];A.U.",10,0.,50.);
	TH1D* histMu1Pt_truth_2to3                  = new TH1D("hMu1Pt_truth_2to3","#mu_{1} p_{T} for m(#mu_{2}-tk) = 2-3 GeV; #mu_{1} p_{T} [GeV];A.U.",10,0.,50.);
	TH1D* histMu1Pt_truth_3toInf                = new TH1D("hMu1Pt_truth_3toInf","#mu_{1} p_{T} for m(#mu_{2}-tk) = 3-Inf GeV; #mu_{1} p_{T} [GeV];A.U.",10,0.,50.);

	// mu1 pT plots in bins of m2 - sideband (soft track pT 1-2.5)
	TH1D* histMu1Pt_side_1to2p5_0to1            = new TH1D("hMu1Pt_side_1to2p5_0to1","#mu_{1} p_{T} for m(#mu_{2}-tk) = 0-1 GeV (soft tk p_{T} = 1-2.5 GeV); #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);
	TH1D* histMu1Pt_side_1to2p5_1to2            = new TH1D("hMu1Pt_side_1to2p5_1to2","#mu_{1} p_{T} for m(#mu_{2}-tk) = 1-2 GeV (soft tk p_{T} = 1-2.5 GeV); #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);
	TH1D* histMu1Pt_side_1to2p5_2to3            = new TH1D("hMu1Pt_side_1to2p5_2to3","#mu_{1} p_{T} for m(#mu_{2}-tk) = 2-3 GeV (soft tk p_{T} = 1-2.5 GeV); #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);
	TH1D* histMu1Pt_side_1to2p5_3toInf          = new TH1D("hMu1Pt_side_1to2p5_3toInf","#mu_{1} p_{T} for m(#mu_{2}-tk) = 3-Inf GeV (soft tk p_{T} = 1-2.5 GeV); #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);

	// mu1 pT plots in bins of m2 - sideband (soft track pT 1-1.5)
	TH1D* histMu1Pt_side_1to1p5_0to1            = new TH1D("hMu1Pt_side_1to1p5_0to1","#mu_{1} p_{T} for m(#mu_{2}-tk) = 0-1 GeV (soft tk p_{T} = 1-1.5 GeV); #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);
	TH1D* histMu1Pt_side_1to1p5_1to2            = new TH1D("hMu1Pt_side_1to1p5_1to2","#mu_{1} p_{T} for m(#mu_{2}-tk) = 1-2 GeV (soft tk p_{T} = 1-1.5 GeV); #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);
	TH1D* histMu1Pt_side_1to1p5_2to3            = new TH1D("hMu1Pt_side_1to1p5_2to3","#mu_{1} p_{T} for m(#mu_{2}-tk) = 2-3 GeV (soft tk p_{T} = 1-1.5 GeV); #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);
	TH1D* histMu1Pt_side_1to1p5_3toInf          = new TH1D("hMu1Pt_side_1to1p5_3toInf","#mu_{1} p_{T} for m(#mu_{2}-tk) = 3-Inf GeV (soft tk p_{T} = 1-1.5 GeV); #mu_{1} p_{T} [GeV];A.U.",10,10.,50.);
	*/

	// -------------------
	// 2D plots of m1 V m2
	// -------------------

	// 2D plots of m1 Vs m2 - signal region
	TH2D* histM1vsM2                            = new TH2D("hM1vsM2","m(#mu_{1}-tk) vs m(#mu_{2}-tk) (signal region);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",nBinsX,&massBins[0],nBinsX,&massBins[0]);
	TH2D* histM1timesM1                         = new TH2D("hM1timesM2","m(sideband) #times m(sideband) (signal region);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",nBinsX,&massBins[0],nBinsX,&massBins[0]);
	TH2D* histM1vsM2_ModMass                    = new TH2D("hM1vsM2_ModMass","m(#mu_{1}-tk) vs m(#mu_{2}-tk) (signal region);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",nBinsX,&massBins[0],nBinsX,&massBins[0]);
	TH2D* histM1timesM1_ModMass                 = new TH2D("hM1timesM2_ModMass","m(sideband) #times m(sideband) (signal region);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",nBinsX,&massBins[0],nBinsX,&massBins[0]);

	// 2D plots of m1 Vs m2 - control region B
	TH2D* histM1vsM2_regionA                    = new TH2D("hM1vsM2_regionA","m(#mu_{1}-tk) vs m(#mu_{2}-tk) (control region A);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",nBinsX,&massBins[0],nBinsX,&massBins[0]);
	TH2D* histM1timesM1_regionA                 = new TH2D("hM1timesM2_regionA","m(sideband) #times m(sideband) (control region A);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",nBinsX,&massBins[0],nBinsX,&massBins[0]);
	TH2D* histM1vsM2_regionA_ModMass            = new TH2D("hM1vsM2_regionA_ModMass","m(#mu_{1}-tk) vs m(#mu_{2}-tk) (control region A);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",nBinsX,&massBins[0],nBinsX,&massBins[0]);
	TH2D* histM1timesM1_regionA_ModMass         = new TH2D("hM1timesM2_regionA_ModMass","m(sideband) #times m(sideband) (control region A);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",nBinsX,&massBins[0],nBinsX,&massBins[0]);

	// 2D plots of m1 Vs m2 - control region B
	TH2D* histM1vsM2_regionB                    = new TH2D("hM1vsM2_regionB","m(#mu_{1}-tk) vs m(#mu_{2}-tk) (control region B);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",nBinsX,&massBins[0],nBinsX,&massBins[0]);
	TH2D* histM1timesM1_regionB                 = new TH2D("hM1timesM2_regionB","m(sideband) #times m(sideband) (control region B);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",nBinsX,&massBins[0],nBinsX,&massBins[0]);
	TH2D* histM1vsM2_regionB_ModMass            = new TH2D("hM1vsM2_regionB_ModMass","m(#mu_{1}-tk) vs m(#mu_{2}-tk) (control region B);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",nBinsX,&massBins[0],nBinsX,&massBins[0]);
	TH2D* histM1timesM1_regionB_ModMass         = new TH2D("hM1timesM2_regionB_ModMass","m(sideband) #times m(sideband) (control region B);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",nBinsX,&massBins[0],nBinsX,&massBins[0]);

	// 2D plots of m1 Vs m2 - signal region where tau cand only have to pass looser iso IP cuts
	// TH2D* histM1vsM2_loosetau                = new TH2D("hM1vsM2_loosetau","m(#mu_{1}-tk) vs m(#mu_{2}-tk) (signal region looser IP on tau cand);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",nBinsX,&massBins[0],nBinsX,&massBins[0]);
	// TH2D* histM1timesM1_loosetau             = new TH2D("hM1timesM2_loosetau","m(sideband) #times m(sideband) (signal region looser IP on tau cand);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",nBinsX,&massBins[0],nBinsX,&massBins[0]);

	// 2D plots of m1 Vs m2 - sideband (soft track pT 1-2.5)
	// TH2D* histM1vsM2_side_1to2p5             = new TH2D("hM1vsM2_side_1to2p5","m(#mu_{1}-tk) vs m(#mu_{2}-tk) (soft tk p_{T} = 1-2.5 GeV);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",nBinsX,&massBins[0],nBinsX,&massBins[0]);
	// TH2D* histM1timesM1_side_1to2p5          = new TH2D("hM1timesM2_side_1to2p5","m(sideband) #times m(sideband) (soft tk p_{T} = 1-2.5 GeV);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",nBinsX,&massBins[0],nBinsX,&massBins[0]);
	// TH2D* histM1vsM2_side_1to2p5_ModMass     = new TH2D("hM1vsM2_side_1to2p5_ModMass","m(#mu_{1}-tk) vs m(#mu_{2}-tk) (soft tk p_{T} = 1-2.5 GeV);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",nBinsX,&massBins[0],nBinsX,&massBins[0]);
	// TH2D* histM1timesM1_side_1to2p5_ModMass  = new TH2D("hM1timesM2_side_1to2p5_ModMass","m(sideband) #times m(sideband) (soft tk p_{T} = 1-2.5 GeV);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",nBinsX,&massBins[0],nBinsX,&massBins[0]);

	// 2D plots of m1 Vs m2 - sideband (soft track pT 1-2.5) with looser tau cand IP
	// TH2D* histM1vsM2_side_1to2p5_loosetau    = new TH2D("hM1vsM2_side_1to2p5_loosetau","m(#mu_{1}-tk) vs m(#mu_{2}-tk) (soft tk p_{T} = 1-2.5 GeV);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",nBinsX,&massBins[0],nBinsX,&massBins[0]);
	// TH2D* histM1timesM1_side_1to2p5_loosetau = new TH2D("hM1timesM2_side_1to2p5_loosetau","m(sideband) #times m(sideband) (soft tk p_{T} = 1-2.5 GeV);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",nBinsX,&massBins[0],nBinsX,&massBins[0]);

	// 2D plots of m1 Vs m2 - sideband (soft track pT 1-1.5)
	// TH2D* histM1vsM2_side_1to1p5             = new TH2D("hM1vsM2_side_1to1p5","m(#mu_{1}-tk) vs m(#mu_{2}-tk) (soft tk p_{T} = 1-1.5 GeV);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",nBinsX,&massBins[0],nBinsX,&massBins[0]);
	// TH2D* histM1timesM1_side_1to1p5          = new TH2D("hM1timesM2_side_1to1p5","m(sideband) #times m(sideband) (soft tk p_{T} = 1-1.5 GeV);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]",nBinsX,&massBins[0],nBinsX,&massBins[0]);

	// ------------------
	// Testing histograms
	// ------------------
	// TH1D* histN_JPsi                         = new TH1D("hN_JPsi","Fraction of J/#Psi in events signal region, 3 < m_{1} < 3.2 GeV;Number of J/#Psi; Fraction", 3, 0, 3);
	// TH1D* histN_JPsi_side_1to2p5             = new TH1D("hN_JPsi_side_1to2p5","Fraction of J/#Psi in events sideband (soft tk p_{T} = 1-2.5 GeV), 3 < m_{1} < 3.2 GeV;Number of J/#Psi; Fraction", 3, 0, 3);

	// TH1D* histN_2p5_looseIP                  = new TH1D("hN_2p5_looseIP","Number of tracks with p_{T} > 2.5 GeV with loose IP;Number of tracks; Fraction", 5, 0, 5);
	// TH1D* histN_2p5_OS_looseIP               = new TH1D("hN_2p5_OS_looseIP","Number of OS tracks with p_{T} > 2.5 GeV with loose IP;Number of tracks; Fraction", 5, 0, 5);

	// TH1D* histN_2p5                          = new TH1D("hN_2p5","Number of tracks with p_{T} > 2.5 GeV ;Number of tracks; Fraction", 5, 0, 5);
	// TH1D* histN_2p5_OS                       = new TH1D("hN_2p5_OS","Number of OS tracks with p_{T} > 2.5 GeV ;Number of tracks; Fraction", 5, 0, 5);

	// int nMu(0);
	// int n1(0), n2(0), nMuPass(0);

	// ----------------------------------------------------
	// For ARC study - 2nd muon has no track requirements
	// (only pt>10, eta < 2.1, deltaR(mu-mu) > 2)
	// ----------------------------------------------------
	// Mu/TK1 here has signal requirements for tracks
	//
	TH1D* histMu1PT_noTk2 = new TH1D("hMu1PT_noTk2", "", 100, 0., 100.);
	TH1D* histMu2PT_noTk2 = new TH1D("hMu2PT_noTk2", "", 100, 0., 100.);
	TH1D* histMuPT_noTk2 = new TH1D("hMuPT_noTk2", "", 100, 0., 100.);
	TH1D* histMu1Eta_noTk2 = new TH1D("hMu1Eta_noTk2", "", 50, -2.5, 2.5);
	TH1D* histMu2Eta_noTk2 = new TH1D("hMu2Eta_noTk2", "", 50, -2.5, 2.5);
	TH1D* histMuEta_noTk2 = new TH1D("hMuEta_noTk2", "", 50, -2.5, 2.5);
	TH1D* histTk1PT_noTk2 = new TH1D("hTk1PT_noTk2", "", 100, 0., 100.);
	TH1D* histTk2PT_noTk2 = new TH1D("hTk2PT_noTk2", "", 100, 0., 100.);
	TH1D* histTkPT_noTk2 = new TH1D("hTkPT_noTk2", "", 100, 0., 100.);
	TH1D* histTk1Eta_noTk2 = new TH1D("hTk1Eta_noTk2", "", 50, -2.5, 2.5);
	TH1D* histTk2Eta_noTk2 = new TH1D("hTk2Eta_noTk2", "", 50, -2.5, 2.5);
	TH1D* histTkEta_noTk2 = new TH1D("hTkEta_noTk2", "", 50, -2.5, 2.5);
	TH1D* histDRMuTk1_noTk2 = new TH1D("hDRMuTk1_noTk2", "", 50, 0., 0.5);
	TH1D* histDRMuTk2_noTk2 = new TH1D("hDRMuTk2_noTk2", "", 50, 0., 0.5);
	TH1D* histDRMuTk_noTk2 = new TH1D("hDRMuTk_noTk2", "", 50, 0., 0.5);
	TH1D* histMassMuTk1_noTk2 = new TH1D("hMassMuTk1_noTk2", "", 10, 0., 10.);
	TH1D* histMassMuTk2_noTk2 = new TH1D("hMassMuTk2_noTk2", "", 10, 0., 10.);
	TH1D* histMassMuTk_noTk2 = new TH1D("hMassMuTk_noTk2", "", 10, 0., 10.);
	TH1D* histMassMuTkModMass_noTk2 = new TH1D("hMassMuTkModMass_noTk2", "", 10, 0., 10.);

	// Store masses for modifiying mass later on
	const double pionMass = 0.13957018;
	const double muonMass = 0.1056583715;
	const double electronMass = 0.000510998928;

	///////////////////////
	// Loop over events  //
	///////////////////////
	Long64_t numberOfEntries = getNumberEvents(treeReader, &pOpts);
	cout << "Running over " << numberOfEntries << " events" << endl;

	bool stop = false; // used to stop the loop, for debugging/testing
	int last_pc = -1;
	for(Int_t entry = 0; entry < numberOfEntries && !stop; ++entry){

		// output something every 5%
		int pc = int(100 * entry/double(numberOfEntries));
		if(pc % 5 == 0 && pc != last_pc) {
			std::cout << pc << "% progress (event " << entry << ")" << endl;
			last_pc = pc;
		}

		// Load selected branches with data from specified event
		treeReader->ReadEntry(entry);

		// cout << "*** Event" << endl;

		if (branchGenMuons->GetEntries() < 2) continue; // skip if <2 muons! (alhtough pointless for HLT samples)

		////////////////////////////////////////////////////////////
		// Get the hard interaction particles for signal MC truth //
		// No selection cuts applied (only >=2 muons)             //
		////////////////////////////////////////////////////////////

		GenParticle *charged1a(nullptr);
		GenParticle *charged1b(nullptr);
		GenParticle *charged2a(nullptr);
		GenParticle *charged2b(nullptr);
/*
		if (doSignal) {
			GenParticle *a1(nullptr), *a2(nullptr);
			// Get a0s
			for(int j = 0; j < branchAll->GetEntries(); j++){
				GenParticle* cand = (GenParticle*) branchAll->At(j);
				// cout << j << " ID: " << cand->PID << " status: " << cand->Status << endl;

				if ((fabs(cand->PID)==36) && (fabs(cand->Status)==62)){
					if (a1==0){
						a1=cand;
						// cout << "found first a1 at " << j << endl;
					} else {
						// cout << "found second a1 at " << j << endl;
						a2=cand;
					}
				}
			}

			// Get the tau daughters from a1 and a2
			GenParticle *tau1a(nullptr), *tau1b(nullptr), *tau2a(nullptr), *tau2b(nullptr);
			tau1a = (GenParticle*) branchAll->At(a1->D1);
			tau1b = (GenParticle*) branchAll->At(a1->D2);
			tau2a = (GenParticle*) branchAll->At(a2->D1);
			tau2b = (GenParticle*) branchAll->At(a2->D2);

			TLorentzVector tau1aMom,tau1bMom, tau2aMom, tau2bMom;
			tau1aMom = tau1a->P4();
			tau1bMom = tau1b->P4();
			tau2aMom = tau2a->P4();
			tau2bMom = tau2b->P4();

			charged1a = getChargedObject(branchAll, tau1a);
			charged1b = getChargedObject(branchAll, tau1b);
			charged2a = getChargedObject(branchAll, tau2a);
			charged2b = getChargedObject(branchAll, tau2b);

			if (charged1a && charged1b && charged2a && charged2b){

				// To hold mu and tracks from each tau
				GenParticle* muTruth1(nullptr);
				GenParticle* trackTruth1(nullptr);
				GenParticle* muTruth2(nullptr);
				GenParticle* trackTruth2(nullptr);

				// Assign charged products to be mu or track
				bool truth1HasMu = assignMuonAndTrack(muTruth1, trackTruth1, *charged1a, *charged1b);
				bool truth2HasMu = assignMuonAndTrack(muTruth2, trackTruth2, *charged2a, *charged2b);

				// NOTE: muons are NOT pT ordered

				if (!truth1HasMu || !truth2HasMu) {
					// cout << "Problem, no truth mu for 1 and/or 2!" << endl;
				} else {

					// Assign system "1" to higher pT muon
					// Swap obj if necessary
					if (muTruth1->PT < muTruth2->PT) {
						GenParticle* tempMu = muTruth1;
						GenParticle* tempTk = trackTruth1;
						muTruth1 = muTruth2;
						trackTruth1 = trackTruth2;
						muTruth2 = tempMu;
						trackTruth2 = tempTk;
					}

					// Randomly swap trk-mu pairs 1<->2 if desired
					if(swapMuRandomly){
						double randNum = (double)rand() / RAND_MAX;
						if (randNum > 0.5){
							GenParticle* tempMu = muTruth1;
							GenParticle* tempTk = trackTruth1;
							muTruth1 = muTruth2;
							trackTruth1 = trackTruth2;
							muTruth2 = tempMu;
							trackTruth2 = tempTk;
						}
					}

					// Do m1 distribution in bins of m2 - for MC truth (is it actually correlated?)
					double m1 = (muTruth1->P4()+trackTruth1->P4()).M();
					double m2 = (muTruth2->P4()+trackTruth2->P4()).M();


					// plot mu-tk system properties (MC truth)
					// cout << m1 << "     " << m2 << endl;
					if(m2 < 1.){
						histM1_truth_0to1->Fill(m1);
						histMu1Pt_truth_0to1->Fill(muTruth1->PT);
					} else if (m2 < 2.){
						histM1_truth_1to2->Fill(m1);
						histMu1Pt_truth_1to2->Fill(muTruth1->PT);
					} else if (m2 < 3.){
						histM1_truth_2to3->Fill(m1);
						histMu1Pt_truth_2to3->Fill(muTruth1->PT);
					} else{
						histM1_truth_3toInf->Fill(m1);
						histMu1Pt_truth_3toInf->Fill(muTruth1->PT);
					}
				}

				if (!muTruth1) delete muTruth1;
				if (!trackTruth1) delete trackTruth1;
				if (!muTruth2) delete muTruth2;
				if (!trackTruth2) delete trackTruth2;
			} // end if(charged1a...)
		} // end if(doSignal)
*/

		/////////////////////////////////////////////////////////////////////////
		// Now, get the two highest pT muons in the event that pass selection, //
		// store pointers to the Track particles and 4-momenta                 //
		// (Use tracks for muons as store more info about position)
		/////////////////////////////////////////////////////////////////////////

		// Track *candTk(nullptr);

		// Fill vectors with muons, based on pT
		std::vector<Track*> muons10to17;
		std::vector<Track*> muons17toInf;
		for (int i = 0; i < branchGenMuons->GetEntries(); i++){
			Track* cand = (Track*) branchGenMuons->At(i);
			if (cand->PT > 17) {
				muons17toInf.push_back(cand);
			} else if (cand->PT > 10) {
				muons10to17.push_back(cand);
			}
		}

		// Check to see if we can skip the event if not enough muons
		if (!(muons17toInf.size() >= 1 && (muons17toInf.size() + muons10to17.size()) >= 2)) continue;

		// Sort both vectors by descending pT
		std::sort(muons17toInf.begin(), muons17toInf.end(), sortByPT<Track>);
		std::sort(muons10to17.begin(), muons10to17.end(), sortByPT<Track>);

		// // For muon pT and eta, we want fewer cuts on muons - just SS, IP, Eta
		// std::pair<Track*, Track*> p_loose = testMuons(muons17toInf, muons10to17, &checkMuonsIPSSEta);
		// if (p_loose.first && p_loose.second) {
		// 	histMu1Pt_fine_IPSS->Fill(p_loose.first->PT);
		// 	histMu1Eta_fine_IPSS->Fill(p_loose.first->Eta);
		// 	histMu2Pt_fine_IPSS->Fill(p_loose.second->PT);
		// 	histMu2Eta_fine_IPSS->Fill(p_loose.second->Eta);
		// }

		// // For muon pT and eta, we use standard selection, no tracks
		// std::pair<Track*, Track*> p_loose_DR = testMuons(muons17toInf, muons10to17, &checkMuonsIPSSDREta, deltaR);
		// if (p_loose_DR.first && p_loose_DR.second) {
		// 	histMu1Pt_fine_IPSSDR->Fill(p_loose_DR.first->PT);
		// 	histMu1Eta_fine_IPSSDR->Fill(p_loose_DR.first->Eta);
		// 	histMu2Pt_fine_IPSSDR->Fill(p_loose_DR.second->PT);
		// 	histMu2Eta_fine_IPSSDR->Fill(p_loose_DR.second->Eta);
		// 	HardMuonPtSoftMuonPt_DimuonsH->Fill(p_loose_DR.first->PT);
		// }

		// For rest of plots we want our muons to pass signal selection
		// Make pairs, see if they pass all cuts (SS, eta, deltaR, dZ, d0)
		// If they do, store in mu1 and mu2 (mu1 has higher pT)
		std::pair<Track*, Track*> p = testMuons(muons17toInf, muons10to17, &checkMuons, deltaR);
		Track* mu1 = p.first;
		Track* mu2 = p.second;

		if (!(p.first && p.second)) continue;


		// cout << " >>>>> Muon pair details: " << endl;
		// cout << " >>>>> Mu1 pT " << mu1->PT << " charge: " << mu1->Charge << " eta " << mu1->Eta << endl;
		// cout << " >>>>> Mu2 pT " << mu2->PT << " charge: " << mu2->Charge << " eta " << mu2->Eta << endl;
		// cout << " >>>>> deltaR " << mu1->P4().DeltaR(mu2->P4()) << endl;

		// Now randomly swap mu1 - mu2
		Track *origMu1(nullptr), *origMu2(nullptr);
		origMu1 = mu1;
		origMu2 = mu2;
		bool swapped = false;
		if (swapMuRandomly){
			// TLorentzVector mu1MomTmp = mu1->P4();
			// TLorentzVector mu2MomTmp = mu2->P4();
			double randNum = (double)rand() / RAND_MAX;
			if (randNum > 0.5){
				swapped = true;
				mu1 = origMu2;
				mu2 = origMu1;
			}
		}

		TLorentzVector mu1Mom = mu1->P4();
		TLorentzVector mu2Mom = mu2->P4();

		// actually setting the mass to it's proper mass
		TLorentzVector mu1MomMass = mu1->P4();
		mu1MomMass.SetPtEtaPhiM(mu1MomMass.Pt(), mu1MomMass.Eta(), mu1MomMass.Phi(), muonMass);
		TLorentzVector mu2MomMass = mu2->P4();
		mu2MomMass.SetPtEtaPhiM(mu2MomMass.Pt(), mu2MomMass.Eta(), mu2MomMass.Phi(), muonMass);

		/////////////////////////////////
		// Look at tracks around muons //
		/////////////////////////////////

		// Vectors of tracks with pT > 1, within dR < 0.5 of respective muons + other cuts
		// so tk1 is the track nearest to muon1, (may or may not be highest pT, depends if random swapping is on)
		std::vector<Track*> tk1_1;
		std::vector<Track*> tk2_1;

		// same but with pT >2.5
		std::vector<Track*> tk1_2p5;
		std::vector<Track*> tk2_2p5;

		// same but with pT > 2.5, OS to muon
		std::vector<Track*> tk1_2p5_OS;
		std::vector<Track*> tk2_2p5_OS;

		// same but with pT >2.5, looser (iso) IP cuts
		std::vector<Track*> tk1_2p5_looseip;
		std::vector<Track*> tk2_2p5_looseip;

		// same but with pT > 2.5, OS to muon, , looser (iso) IP cuts
		std::vector<Track*> tk1_2p5_OS_looseip;
		std::vector<Track*> tk2_2p5_OS_looseip;

		// same but with 1 < pT < 2.5 (for sideband) (loose/iso IP cuts)
		std::vector<Track*> tk1_1to2p5;
		std::vector<Track*> tk2_1to2p5;

		// same but with 1 < pT < 2.5 (for sideband) with tighter IP cuts (tau candidate IP)
		std::vector<Track*> tk1_1to2p5_taucand;
		std::vector<Track*> tk2_1to2p5_taucand;

		// same but with 1 < pT < 2.5 (for sideband) with looser IP cuts (tau candidate IP)
		std::vector<Track*> tk1_1to2p5_taucand_looseip;
		std::vector<Track*> tk2_1to2p5_taucand_looseip;

		// same but with 1 < pT < 1.5 (for sideband)
		std::vector<Track*> tk1_1to1p5;
		std::vector<Track*> tk2_1to1p5;

		// 2.5 < pT < 6 (for control region B) fail 1 prong tau criteria
		std::vector<Track*> tk1_2p5to6_fail1prong;
		std::vector<Track*> tk2_2p5to6_fail1prong;

		///////////////////////////////
		// FILL TRACK VECTORS HERE //
		///////////////////////////////
		Track *candTk(nullptr);
		for(int a = 0; a < branchTracks->GetEntries(); a++){
			candTk = (Track*) branchTracks->At(a);

			// Check it isn't the same object as the muons!
			if ((candTk->PT != mu1->PT) && (candTk->PT != mu2->PT)){

				if (checkTrackLoose(candTk)) {

					// Store track in suitable vector
					fillTrackVectors(candTk, mu1, mu2, &tk1_1, &tk2_1);


					// 1 prong candiates,
					if (checkTrackIPTight(candTk)) {

						// for N_trk = 2or3
						fillTrackVectors(candTk, mu1, mu2, &tk1_1to2p5_taucand, &tk2_1to2p5_taucand);

						if (checkTrackPTTight(candTk)) {
							// tau candidate must have pT > 2.5 GeV
							fillTrackVectors(candTk, mu1, mu2, &tk1_2p5, &tk2_2p5);

							if (checkTkMuOS(candTk, mu1)) {
								fillTrackVectors(candTk, mu1, mu2, &tk1_2p5_OS, &tk2_2p5_OS);
							}
						}
					} else {

						// if (checkTrackPTTight(candTk)) {
							// tau candidate must have pT > 2.5 GeV
							// fillTrackVectors(candTk, mu1, mu2, &tk1_2p5_looseip, & tk2_2p5_looseip);

							// if (checkTkMuOS(candTk, mu1)) {
								// fillTrackVectors(candTk, mu1, mu2, &tk1_2p5_OS_looseip, &tk2_2p5_OS_looseip);
							// }
						// } else {
							fillTrackVectors(candTk, mu1, mu2, &tk1_1to2p5, &tk2_1to2p5);
						// }

						if (do1to1p5 && candTk->PT < 1.5){
							fillTrackVectors(candTk, mu1, mu2, &tk1_1to1p5, &tk2_1to1p5);
						}
					}
					// for loose IP on signal tracks
					// if (checkTrackIPTightAlt(candTk)){
					// 	// for N_trk = 2or3
					// 	fillTrackVectors(candTk, mu1, mu2, &tk1_1to2p5_taucand_looseip, &tk2_1to2p5_taucand_looseip);

					// 	if (checkTrackPTTight(candTk)){ // alternate IP cut on 1-prong tracks
					// 		fillTrackVectors(candTk, mu1, mu2, &tk1_2p5_looseip, &tk2_2p5_looseip);
					// 		if (checkTkMuOS(candTk, mu1)) {
					// 			fillTrackVectors(candTk, mu1, mu2, &tk1_2p5_OS_looseip, &tk2_2p5_OS_looseip);
					// 		}
					// 	}
					// }

					// control region B tracks
					// must fail at least one of the tau 1-prong criteria
					if (checkTrackPTTight(candTk)
						&& candTk->PT < 6.0
						&& ((fabs(candTk->Zd) > 0.4) || ((fabs(calcDxy(candTk))) > 0.2) || !checkTkMuOS(candTk, mu1))) {
						fillTrackVectors(candTk, mu1, mu2, &tk1_2p5to6_fail1prong, &tk2_2p5to6_fail1prong);
					}
				}
			} // End of track selection criteria
		} // End of track loop

		// if ((tk1_2p5to6_fail1prong.size() == 1 || tk1_2p5to6_fail1prong.size() == 2)
		// 		&& (tk2_2p5to6_fail1prong.size() == 1 || tk2_2p5to6_fail1prong.size() == 2)) {
		// 	cout << "region B tracks: " << tk1_2p5to6_fail1prong.size() << " - " << tk2_2p5to6_fail1prong.size() << endl;
		// }

		////////////////////////////////////////
		// FILL HISTS FOR DIFFERENT REGIONS //
		////////////////////////////////////////

		// sort track vectors by descending pT
		std::sort(tk1_2p5_OS.begin(), tk1_2p5_OS.end(), sortByPT<Track>);
		std::sort(tk2_2p5_OS.begin(), tk2_2p5_OS.end(), sortByPT<Track>);

		// Do some preselection here that's common to all regions of consideration
		// exactly 1 OS track with pT > 2.5 within ∆R < 0.5 of muon
		// Don't put in tk1_2p5.size() == 1 as it screws with Region B
		if (tk1_2p5_OS.size() == 1 && tk2_2p5_OS.size() == 1) {


			/////////////////////
			// COMMON STUFF //
			/////////////////////

			// Get 4-momenta and invariant masses, since we never change the
			// objects used
			TLorentzVector track1Mom = tk1_2p5_OS[0]->P4();
			TLorentzVector track2Mom = tk2_2p5_OS[0]->P4();

			// random since mu1Mom andmu2Mom are randomly assigned (if selected at top)
			double m1 = (mu1Mom+tk1_2p5_OS[0]->P4()).M();
			double m2 = (mu2Mom+tk2_2p5_OS[0]->P4()).M();

			// recalculate 4 momenta, setting non-zero masses
			// for tracks, use muonMass if muon, electornMass if electron,
			// or pion mass otherwise
			TLorentzVector tk1MomMass = tk1_2p5_OS[0]->P4();
			double mass = pionMass;
			if (abs(tk1_2p5_OS[0]->PID) == 13) {
				mass = muonMass;
			} else if (abs(tk1_2p5_OS[0]->PID) == 11) {
				mass = electronMass;
			}
			tk1MomMass.SetPtEtaPhiM(tk1MomMass.Pt(), tk1MomMass.Eta(), tk1MomMass.Phi(), mass);
			TLorentzVector tk2MomMass = tk2_2p5_OS[0]->P4();
			mass = pionMass;
			if (abs(tk1_2p5_OS[0]->PID) == 13) {
				mass = muonMass;
			} else if (abs(tk1_2p5_OS[0]->PID) == 11) {
				mass = electronMass;
			}
			tk2MomMass.SetPtEtaPhiM(tk2MomMass.Pt(), tk2MomMass.Eta(), tk2MomMass.Phi(), muonMass);

			double m1Mod = (mu1MomMass+tk1MomMass).M();
			double m2Mod = (mu2MomMass+tk2MomMass).M();

			/////////////////////////
			// SIGNAL SELECTION    //
			/////////////////////////
			// Only 1 track within ∆R < 0.5 of muon has pT > 1,
			// and that track must have pT > 2.5, and be oppsite charge to muon
			if (tk1_1.size() == 1 && tk2_1.size() == 1
				&& tk1_2p5.size() == 1 && tk2_2p5.size() == 1 )
			{

				histM1->Fill(m1);
				histM2->Fill(m2);
				// histM1_fine->Fill(m1);
				// histM2_fine->Fill(m2);

				histM1_ModMass->Fill(m1Mod);
				histM2_ModMass->Fill(m2Mod);

				// Fill symmetrically to increase stats
				histM1vsM2->Fill(m1,m2);
				histM1vsM2->Fill(m2,m1);
				histM1vsM2_ModMass->Fill(m1Mod,m2Mod);
				histM1vsM2_ModMass->Fill(m2Mod,m1Mod);

				// if(m2 < 1.){
				// 	histM1_0to1->Fill(m1);
				// 	histMu1Pt_0to1->Fill(mu1->PT);
				// } else if (m2 < 2.){
				// 	histM1_1to2->Fill(m1);
				// 	histMu1Pt_1to2->Fill(mu1->PT);
				// } else if (m2 < 3.){
				// 	histM1_2to3->Fill(m1);
				// 	histMu1Pt_2to3->Fill(mu1->PT);
				// } else {
				// 	histM1_3toInf->Fill(m1);
				// 	histMu1Pt_3toInf->Fill(mu1->PT);
				// }

			}


			////////////////////////////////////////////////////////////////////////
			// CONTROL REGION A - where each muon must have 1 track satisfying //
			// the signal selection, and 1 ot 2 extra tracks nearby, with      //
			// pT 1 - 2.5, with loose IP cuts (ie not 1-prong criteria)        //
			////////////////////////////////////////////////////////////////////////
			if ((tk1_1.size() == 2 || tk1_1.size() == 3)
				&& (tk2_1.size() == 2 || tk2_1.size() == 3)
				&& (tk1_1to2p5.size() == 1 || tk1_1to2p5.size() == 2)
				&& (tk2_1to2p5.size() == 1 || tk2_1to2p5.size() == 2)
				&& tk1_2p5.size() == 1 && tk2_2p5.size() == 1 ){

				histM1_regionA->Fill(m1);
				histM2_regionA->Fill(m2);

				// histM1_fine_regionA->Fill(m1);
				// histM2_fine_regionA->Fill(m2);

				histM1_regionA_ModMass->Fill(m1Mod);
				histM2_regionA_ModMass->Fill(m2Mod);

				// Fill symmetrically to increase stats
				histM1vsM2_regionA->Fill(m1,m2);
				histM1vsM2_regionA->Fill(m2,m1);

				histM1vsM2_regionA_ModMass->Fill(m1Mod,m2Mod);
				histM1vsM2_regionA_ModMass->Fill(m2Mod,m1Mod);
			}


			///////////////////////////////////////////////////////////////////////
			// CONTROL REGION B - where each muon must have 1 or 2 exra tracks //
			// closeby, which fail the 1-prong criteria., with 2.5 < pT < 6.   //
			// Still require normal 1 OS tight track with pT > 2.5
			///////////////////////////////////////////////////////////////////////
			if ((tk1_1.size() == 2 || tk1_1.size() == 3)
				&& (tk2_1.size() == 2 || tk2_1.size() == 3)
				&& tk1_1to2p5.size() == 0 && tk2_1to2p5.size() == 0
				&& (tk1_2p5to6_fail1prong.size() == 1 || tk1_2p5to6_fail1prong.size() == 2)
				&& (tk2_2p5to6_fail1prong.size() == 1 || tk2_2p5to6_fail1prong.size() == 2)) {

				histM1_regionB->Fill(m1);
				histM2_regionB->Fill(m2);

				// histM1_fine_regionB->Fill(m1);
				// histM2_fine_regionB->Fill(m2);

				histM1_regionB_ModMass->Fill(m1Mod);
				histM2_regionB_ModMass->Fill(m2Mod);

				// Fill symmetrically to increase stats
				histM1vsM2_regionB->Fill(m1,m2);
				histM1vsM2_regionB->Fill(m2,m1);

				histM1vsM2_regionB_ModMass->Fill(m1Mod,m2Mod);
				histM1vsM2_regionB_ModMass->Fill(m2Mod,m1Mod);
			}


			//////////////////////////////////////////////////////////////////////
			// SIDEBAND REGION - where mu2 has 1, 2, 3 additional soft tracks with tight IP cuts//
			// and mu1 satisfies signal selection
			//////////////////////////////////////////////////////////////////////
	/*		if (   tk1_1.size() == 1
				&& tk1_2p5.size() == 1 && tk1_2p5_OS.size() == 1
				&& tk2_2p5.size() == 1 && tk2_2p5_OS.size() == 1
				){

				double m1(0);

				m1 = (mu1Mom+tk1_2p5_OS[0]->P4()).M();

				if ( //tk2_1.size() == 2
					tk2_1to2p5_taucand.size() == 1
					&& tk2_1to2p5.size() == 1
					) {

						histM1_Ntk2_2->Fill(m1);
						histM1_Ntk2_2or3->Fill(m1);

					}

				if (// tk2_1.size() == 3
					tk2_1to2p5_taucand.size() == 2
					&& tk2_1to2p5.size() == 2
					) {
						// m1 = (mu1Mom+tk1_2p5_OS[0]->P4()).M();

						histM1_Ntk2_3->Fill(m1);
						histM1_Ntk2_2or3->Fill(m1);
					}

				if ( //tk2_1.size() == 4
					tk2_1to2p5_taucand.size() == 3
					&& tk2_1to2p5.size() == 3
					) {
						// m1 = (mu1Mom+tk1_2p5_OS[0]->P4()).M();

						histM1_Ntk2_4->Fill(m1);
					}
			}

			if (   tk2_1.size() == 1
				&& tk2_2p5.size() == 1 && tk2_2p5_OS.size() == 1
				&& tk1_2p5.size() == 1 && tk1_2p5_OS.size() == 1
				){

				double m1(0);

				m1 = (mu2Mom+tk2_2p5_OS[0]->P4()).M();

				if ( //tk1_1.size() == 2
					tk1_1to2p5_taucand.size() == 1
					&& tk1_1to2p5.size() == 1
					) {

						histM1_Ntk2_2->Fill(m1);
						histM1_Ntk2_2or3->Fill(m1);

					}

				if (// tk1_1.size() == 3
					tk1_1to2p5_taucand.size() == 2
					&& tk1_1to2p5.size() == 2
					) {
						// m1 = (mu1Mom+tk2_2p5_OS[0]->P4()).M();

						histM1_Ntk2_3->Fill(m1);
						histM1_Ntk2_2or3->Fill(m1);
					}

				if ( //tk1_1.size() == 4
					tk1_1to2p5_taucand.size() == 3
					&& tk1_1to2p5.size() == 3
					) {
						// m1 = (mu1Mom+tk2_2p5_OS[0]->P4()).M();

						histM1_Ntk2_4->Fill(m1);
					}
			}*/


		} // end of track preselection

		//////////////////////////////////////////////////////
		// ARC STUDY - one muon has no track requirements, other has signal selection //
		//////////////////////////////////////////////////////
		if (tk1_1.size() == 1 && tk1_2p5.size() == 1 && tk1_2p5_OS.size() == 1) {
			histMu1PT_noTk2->Fill(mu1->PT);
			histMuPT_noTk2->Fill(mu1->PT);
			histMu1Eta_noTk2->Fill(mu1->Eta);
			histMuEta_noTk2->Fill(mu1->Eta);

			histTk1PT_noTk2->Fill(tk1_2p5_OS[0]->PT);
			histTk1Eta_noTk2->Fill(tk1_2p5_OS[0]->Eta);
			histDRMuTk1_noTk2->Fill(mu1->P4().DeltaR(tk1_2p5_OS[0]->P4()));

			histTkPT_noTk2->Fill(tk1_2p5_OS[0]->PT);
			histTkEta_noTk2->Fill(tk1_2p5_OS[0]->Eta);
			histDRMuTk_noTk2->Fill(mu1->P4().DeltaR(tk1_2p5_OS[0]->P4()));

			histMassMuTk1_noTk2->Fill((mu1->P4()+tk1_2p5_OS[0]->P4()).M());
			histMassMuTk_noTk2->Fill((mu1->P4()+tk1_2p5_OS[0]->P4()).M());

			// for modified mass of muon & track
			TLorentzVector mu1MomMass = mu1->P4();
			mu1MomMass.SetPtEtaPhiM(mu1MomMass.Pt(), mu1MomMass.Eta(), mu1MomMass.Phi(), muonMass);
			TLorentzVector tk1MomMass = tk1_2p5_OS[0]->P4();
			double mass = pionMass;
			if (tk1_2p5_OS[0]->PID == 11) {
				mass = electronMass;
			} else if (tk1_2p5_OS[0]->PID == 13) {
				mass = muonMass;
			}
			tk1MomMass.SetPtEtaPhiM(tk1MomMass.Pt(), tk1MomMass.Eta(), tk1MomMass.Phi(), mass);
			histMassMuTkModMass_noTk2->Fill((mu1MomMass+tk1MomMass).M());

			if (tk2_1.size() == 1 && tk2_2p5.size() == 1 && tk2_2p5_OS.size() == 1) {
				histMu2PT_noTk2->Fill(mu2->PT);
				histMuPT_noTk2->Fill(mu2->PT);
				histMu2Eta_noTk2->Fill(mu2->Eta);
				histMuEta_noTk2->Fill(mu2->Eta);

				histTk2PT_noTk2->Fill(tk2_2p5_OS[0]->PT);
				histTk2Eta_noTk2->Fill(tk2_2p5_OS[0]->Eta);
				histDRMuTk2_noTk2->Fill(mu2->P4().DeltaR(tk2_2p5_OS[0]->P4()));

				histTkPT_noTk2->Fill(tk2_2p5_OS[0]->PT);
				histTkEta_noTk2->Fill(tk2_2p5_OS[0]->Eta);
				histDRMuTk_noTk2->Fill(mu2->P4().DeltaR(tk2_2p5_OS[0]->P4()));

				histMassMuTk2_noTk2->Fill((mu2->P4()+tk2_2p5_OS[0]->P4()).M());
				histMassMuTk_noTk2->Fill((mu2->P4()+tk2_2p5_OS[0]->P4()).M());

				// for modified mass of muon & track
				TLorentzVector mu2MomMass = mu2->P4();
				mu2MomMass.SetPtEtaPhiM(mu2MomMass.Pt(), mu2MomMass.Eta(), mu2MomMass.Phi(), muonMass);
				TLorentzVector tk2MomMass = tk2_2p5_OS[0]->P4();
				double mass = pionMass;
				if (tk2_2p5_OS[0]->PID == 11) {
					mass = electronMass;
				} else if (tk2_2p5_OS[0]->PID == 13) {
					mass = muonMass;
				}
				tk2MomMass.SetPtEtaPhiM(tk2MomMass.Pt(), tk2MomMass.Eta(), tk2MomMass.Phi(), mass);
				histMassMuTkModMass_noTk2->Fill((mu2MomMass+tk2MomMass).M());
			}
		} else if (tk2_1.size() == 1 && tk2_2p5.size() == 1 && tk2_2p5_OS.size() == 1) {
			histMu1PT_noTk2->Fill(mu2->PT);
			histMuPT_noTk2->Fill(mu2->PT);
			histMu1Eta_noTk2->Fill(mu2->Eta);
			histMuEta_noTk2->Fill(mu2->Eta);

			histTk1PT_noTk2->Fill(tk2_2p5_OS[0]->PT);
			histTk1Eta_noTk2->Fill(tk2_2p5_OS[0]->Eta);
			histDRMuTk1_noTk2->Fill(mu2->P4().DeltaR(tk2_2p5_OS[0]->P4()));

			histTkPT_noTk2->Fill(tk2_2p5_OS[0]->PT);
			histTkEta_noTk2->Fill(tk2_2p5_OS[0]->Eta);
			histDRMuTk_noTk2->Fill(mu1->P4().DeltaR(tk2_2p5_OS[0]->P4()));

			histMassMuTk1_noTk2->Fill((mu2->P4()+tk2_2p5_OS[0]->P4()).M());
			histMassMuTk_noTk2->Fill((mu2->P4()+tk2_2p5_OS[0]->P4()).M());

			// for modified mass of muon & track
			TLorentzVector mu2MomMass = mu2->P4();
			mu2MomMass.SetPtEtaPhiM(mu2MomMass.Pt(), mu2MomMass.Eta(), mu2MomMass.Phi(), muonMass);
			TLorentzVector tk2MomMass = tk2_2p5_OS[0]->P4();
			double mass = pionMass;
			if (tk2_2p5_OS[0]->PID == 11) {
				mass = electronMass;
			} else if (tk2_2p5_OS[0]->PID == 13) {
				mass = muonMass;
			}
			tk2MomMass.SetPtEtaPhiM(tk2MomMass.Pt(), tk2MomMass.Eta(), tk2MomMass.Phi(), mass);
			histMassMuTkModMass_noTk2->Fill((mu2MomMass+tk2MomMass).M());

			if (tk1_1.size() == 1 && tk1_2p5.size() == 1 && tk1_2p5_OS.size() == 1) {
				histMu2PT_noTk2->Fill(mu1->PT);
				histMuPT_noTk2->Fill(mu1->PT);
				histMu2Eta_noTk2->Fill(mu1->Eta);
				histMuEta_noTk2->Fill(mu1->Eta);

				histTk2PT_noTk2->Fill(tk1_2p5_OS[0]->PT);
				histTk2Eta_noTk2->Fill(tk1_2p5_OS[0]->Eta);
				histDRMuTk2_noTk2->Fill(mu1->P4().DeltaR(tk1_2p5_OS[0]->P4()));

				histTkPT_noTk2->Fill(tk1_2p5_OS[0]->PT);
				histTkEta_noTk2->Fill(tk1_2p5_OS[0]->Eta);
				histDRMuTk_noTk2->Fill(mu1->P4().DeltaR(tk1_2p5_OS[0]->P4()));

				histMassMuTk_noTk2->Fill((mu1->P4()+tk1_2p5_OS[0]->P4()).M());
				histMassMuTk2_noTk2->Fill((mu1->P4()+tk1_2p5_OS[0]->P4()).M());

				// for modified mass of muon & track
				TLorentzVector mu1MomMass = mu1->P4();
				mu1MomMass.SetPtEtaPhiM(mu1MomMass.Pt(), mu1MomMass.Eta(), mu1MomMass.Phi(), muonMass);
				TLorentzVector tk1MomMass = tk1_2p5_OS[0]->P4();
				double mass = pionMass;
				if (tk1_2p5_OS[0]->PID == 11) {
					mass = electronMass;
				} else if (tk1_2p5_OS[0]->PID == 13) {
					mass = muonMass;
				}
				tk1MomMass.SetPtEtaPhiM(tk1MomMass.Pt(), tk1MomMass.Eta(), tk1MomMass.Phi(), mass);
				histMassMuTkModMass_noTk2->Fill((mu1MomMass+tk1MomMass).M());
			}
		} // end of arc study block

	} // end of event loop

	// Create sum of m1 and m2 1D hists and rebin
	TH1D * histM = makeSumRebin<TH1D>(histM1, histM2, massBins, "hM");
	histM->SetTitle("m(#mu-tk) in signal region; m(#mu-tk) [GeV];A.U.");

	TH1D * histM_ModMass = makeSumRebin<TH1D>(histM1_ModMass, histM2_ModMass, massBins, "hM_ModMass");
	histM_ModMass->SetTitle("m(#mu-tk) in signal region; m(#mu-tk) [GeV];A.U.");

	TH1D* histM_regionA = makeSumRebin<TH1D>(histM1_regionA, histM2_regionA, massBins, "hM_regionA");
	histM_regionA->SetTitle("m(#mu-tk) in region A;m(#mu-tk) [GeV];A.U.");

	TH1D* histM_regionA_ModMass = makeSumRebin<TH1D>(histM1_regionA_ModMass, histM2_regionA_ModMass, massBins, "hM_regionA_ModMass");
	histM_regionA_ModMass->SetTitle("m(#mu-tk) in region A;m(#mu-tk) [GeV];A.U.");

	TH1D* histM_regionB = makeSumRebin<TH1D>(histM1_regionB, histM2_regionB, massBins, "hM_regionB");
	histM_regionB->SetTitle("m(#mu-tk) in region B;m(#mu-tk) [GeV];A.U.");

	TH1D* histM_regionB_ModMass = makeSumRebin<TH1D>(histM1_regionB_ModMass, histM2_regionB_ModMass, massBins, "hM_regionB_ModMass");
	histM_regionB_ModMass->SetTitle("m(#mu-tk) in region B;m(#mu-tk) [GeV];A.U.");

/*
	TH1D* histM_side_1to2p5 = new TH1D("hM_side_1to2p5","m(#mu-tk) in sideband (soft tk p_{T} = 1 - 2.5 GeV);m(#mu-tk) [GeV];A.U.",
										histM1_side_1to2p5->GetNbinsX(),
										histM1_side_1to2p5->GetXaxis()->GetBinLowEdge(1),
										histM1_side_1to2p5->GetXaxis()->GetBinLowEdge(histM1_side_1to2p5->GetNbinsX()+1));
	histM_side_1to2p5->Add(histM1_side_1to2p5);
	histM_side_1to2p5->Add(histM2_side_1to2p5);
	histM_side_1to2p5 = (TH1D*) histM_side_1to2p5->Rebin(nBinsX, "hM_side_1to2p5", &massBins[0]);

	TH1D* histM_side_1to2p5_loosetau = new TH1D("hM_side_1to2p5_loosetau","m(#mu-tk) in sideband (soft tk p_{T} = 1 - 2.5 GeV);m(#mu-tk) [GeV];A.U.",
										histM1_side_1to2p5_loosetau->GetNbinsX(),
										histM1_side_1to2p5_loosetau->GetXaxis()->GetBinLowEdge(1),
										histM1_side_1to2p5_loosetau->GetXaxis()->GetBinLowEdge(histM1_side_1to2p5_loosetau->GetNbinsX()+1));
	histM_side_1to2p5_loosetau->Add(histM1_side_1to2p5_loosetau);
	histM_side_1to2p5_loosetau->Add(histM2_side_1to2p5_loosetau);
	histM_side_1to2p5_loosetau = (TH1D*) histM_side_1to2p5_loosetau->Rebin(nBinsX, "hM_side_1to2p5_loosetau", &massBins[0]);

	TH1D* histM_side_1to1p5;
	if(do1to1p5) {
		histM_side_1to1p5 = new TH1D("hM_side_1to1p5","m(#mu-tk) in sideband (soft tk p_{T} = 1 - 1.5 GeV);m(#mu-tk) [GeV];A.U.",
										histM1_side_1to1p5->GetNbinsX(),
										histM1_side_1to1p5->GetXaxis()->GetBinLowEdge(1),
										histM1_side_1to1p5->GetXaxis()->GetBinLowEdge(histM1_side_1to1p5->GetNbinsX()+1));

		histM_side_1to1p5->Add(histM1_side_1to1p5);
		histM_side_1to1p5->Add(histM2_side_1to1p5);
		histM_side_1to1p5 = (TH1D*) histM_side_1to1p5->Rebin(nBinsX, "hM_side_1to1p5", &massBins[0]);
	}

	TH1D* histM_loosetau = new TH1D("hM_loosetau", "Inv. Mass of system, full selection; m(#mu-tk) [GeV];A.U.",nBinsX,&massBins[0]);
	TH1D* histM1_rebin_loosetau = (TH1D*) histM1_loosetau->Rebin(nBinsX,"hM1_rebin_loosetau",&massBins[0]);
	TH1D* histM2_rebin_loosetau = (TH1D*) histM2_loosetau->Rebin(nBinsX,"hM2_rebin_loosetau",&massBins[0]);
	histM_loosetau->Add(histM1_rebin_loosetau);
	histM_loosetau->Add(histM2_rebin_loosetau);
*/

	// Print out integrals before we normalise
	printIntegral(histM1);
	printIntegral(histM2);
	printIntegral(histM);
	printIntegral(histM1vsM2);

	printIntegral(histM1_ModMass);
	printIntegral(histM2_ModMass);
	printIntegral(histM_ModMass);
	printIntegral(histM1vsM2_ModMass);

	printIntegral(histM1_regionA);
	printIntegral(histM2_regionA);
	printIntegral(histM_regionA);
	printIntegral(histM1vsM2_regionA);

	printIntegral(histM1_regionA_ModMass);
	printIntegral(histM2_regionA_ModMass);
	printIntegral(histM_regionA_ModMass);
	printIntegral(histM1vsM2_regionA_ModMass);

	printIntegral(histM1_regionB);
	printIntegral(histM2_regionB);
	printIntegral(histM_regionB);
	printIntegral(histM1vsM2_regionB);

	printIntegral(histM1_regionB_ModMass);
	printIntegral(histM2_regionB_ModMass);
	printIntegral(histM_regionB_ModMass);
	printIntegral(histM1vsM2_regionB_ModMass);

	printIntegral(histMuPT_noTk2);
/*
	printIntegral(histM1_side_1to2p5);
	printIntegral(histM2_side_1to2p5);
	printIntegral(histM_side_1to2p5);
	printIntegral(histM1vsM2_side_1to2p5);

	printIntegral(histM1_side_1to2p5_loosetau);
	printIntegral(histM2_side_1to2p5_loosetau);
	printIntegral(histM_side_1to2p5_loosetau);
	printIntegral(histM1vsM2_side_1to2p5_loosetau);

	if (do1to1p5) {
		printIntegral(histM1_side_1to1p5);
		printIntegral(histM2_side_1to1p5);
		printIntegral(histM_side_1to1p5);
		printIntegral(histM1vsM2_side_1to1p5);
	}
	printIntegral(histM1_loosetau);
	printIntegral(histM2_loosetau);
	printIntegral(histM_loosetau);
	printIntegral(histM1vsM2_loosetau);

	printIntegral(histM1_Ntk2_2);
	printIntegral(histM1_Ntk2_3);
	printIntegral(histM1_Ntk2_2or3);
	printIntegral(histM1_Ntk2_4);

	printIntegral(histM1_Ntk2_2_loosetau);
	printIntegral(histM1_Ntk2_3_loosetau);
	printIntegral(histM1_Ntk2_2or3_loosetau);
	printIntegral(histM1_Ntk2_4_loosetau);
*/


	/////////////////
	// PLOT THINGS //
	/////////////////
	std::string app(""); // text to append on end of plot filenames

	if (doSignal) {
		app = "sig";
	} else {
		app = "bg";
	}
	if (swapMuRandomly){
		app += "_muRand";
	}
	if (doHLT) {
		app += "_HLT";
	} else {
		app += "_NoHLT";
	}

	if (doRescale) {
		app += "_rescaleQuantile";
	}

	app += "_2500_arc_testing2";

	app += "_dR";
	app += boost::lexical_cast<std::string>(deltaR);

	if (source == test)
		app += "_TEST";

	// Get directory that input file was in - put plots in there
	std::string directory = getDirectory(chain.GetFile());

	// Get Delphes file config used - last part of directory name
	std::string delph = getDelph(directory);

	TFile* outFile = TFile::Open((directory+"/output_"+delph+"_"+app+".root").c_str(),"UPDATE");
	cout << "Writing to " << outFile->GetName() << endl;

	// Do some normalizing, make copies beforehand for combination plots

	// SIGNAL REGION
	makeCopySave(histM1);
	normaliseHist(histM1);

	makeCopySave(histM1_ModMass);
	normaliseHist(histM1_ModMass);

	makeCopySave(histM2);
	normaliseHist(histM2);

	makeCopySave(histM2_ModMass);
	normaliseHist(histM2_ModMass);

	makeCopySave(histM);
	normaliseHist(histM);

	makeCopySave(histM_ModMass);
	normaliseHist(histM_ModMass);

	makeCopySave(histM1vsM2);
	normaliseHist(histM1vsM2);

	makeCopySave(histM1vsM2_ModMass);
	normaliseHist(histM1vsM2_ModMass);


	// REGION A
	makeCopySave(histM1_regionA);
	normaliseHist(histM1_regionA);

	makeCopySave(histM1_regionA_ModMass);
	normaliseHist(histM1_regionA_ModMass);

	makeCopySave(histM2_regionA);
	normaliseHist(histM2_regionA);

	makeCopySave(histM2_regionA_ModMass);
	normaliseHist(histM2_regionA_ModMass);

	// makeCopySave(histM1_fine_regionA);
	// normaliseHist(histM1_fine_regionA);

	// makeCopySave(histM2_fine_regionA);
	// normaliseHist(histM2_fine_regionA);

	makeCopySave(histM_regionA);
	normaliseHist(histM_regionA);

	makeCopySave(histM_regionA_ModMass);
	normaliseHist(histM_regionA_ModMass);

	makeCopySave(histM1vsM2_regionA);
	normaliseHist(histM1vsM2_regionA);

	makeCopySave(histM1vsM2_regionA_ModMass);
	normaliseHist(histM1vsM2_regionA_ModMass);

	// REGION B
	makeCopySave(histM1_regionB);
	normaliseHist(histM1_regionB);

	makeCopySave(histM1_regionB_ModMass);
	normaliseHist(histM1_regionB_ModMass);

	makeCopySave(histM2_regionB);
	normaliseHist(histM2_regionB);

	makeCopySave(histM2_regionB_ModMass);
	normaliseHist(histM2_regionB_ModMass);

	// makeCopySave(histM1_fine_regionB);
	// normaliseHist(histM1_fine_regionB);

	// makeCopySave(histM2_fine_regionB);
	// normaliseHist(histM2_fine_regionB);

	makeCopySave(histM_regionB);
	normaliseHist(histM_regionB);

	makeCopySave(histM_regionB_ModMass);
	normaliseHist(histM_regionB_ModMass);

	makeCopySave(histM1vsM2_regionB);
	normaliseHist(histM1vsM2_regionB);

	makeCopySave(histM1vsM2_regionB_ModMass);
	normaliseHist(histM1vsM2_regionB_ModMass);


/*	makeCopySave(histM1_side_1to2p5);
	normaliseHist(histM1_side_1to2p5);

	makeCopySave(histM2_side_1to2p5);
	normaliseHist(histM2_side_1to2p5);

	makeCopySave(histM_side_1to2p5);
	normaliseHist(histM_side_1to2p5);

	makeCopySave(histM1vsM2_side_1to2p5);
	normaliseHist(histM1vsM2_side_1to2p5);

	makeCopySave(histM1_side_1to2p5_loosetau);
	normaliseHist(histM1_side_1to2p5_loosetau);

	makeCopySave(histM2_side_1to2p5_loosetau);
	normaliseHist(histM2_side_1to2p5_loosetau);

	makeCopySave(histM_side_1to2p5_loosetau);
	normaliseHist(histM_side_1to2p5_loosetau);

	makeCopySave(histM1vsM2_side_1to2p5_loosetau);
	normaliseHist(histM1vsM2_side_1to2p5_loosetau);

	if (do1to1p5) {
		normaliseHist(histM1_side_1to1p5);
		normaliseHist(histM2_side_1to1p5);
		normaliseHist(histM_side_1to1p5);
		normaliseHist(histM1vsM2_side_1to1p5);
	}

	makeCopySave(histM1_loosetau);
	normaliseHist(histM1_loosetau);

	makeCopySave(histM1_Ntk2_2);
	normaliseHist(histM1_Ntk2_2);

	makeCopySave(histM1_Ntk2_3);
	normaliseHist(histM1_Ntk2_3);

	makeCopySave(histM1_Ntk2_2or3);
	normaliseHist(histM1_Ntk2_2or3);

	makeCopySave(histM1_Ntk2_4);
	normaliseHist(histM1_Ntk2_4);

	makeCopySave(histM1_Ntk2_2_loosetau);
	normaliseHist(histM1_Ntk2_2_loosetau);

	makeCopySave(histM1_Ntk2_3_loosetau);
	normaliseHist(histM1_Ntk2_3_loosetau);

	makeCopySave(histM1_Ntk2_2or3_loosetau);
	normaliseHist(histM1_Ntk2_2or3_loosetau);

	makeCopySave(histM1_Ntk2_4_loosetau);
	normaliseHist(histM1_Ntk2_4_loosetau);

	makeCopySave(histM1_fine);
	normaliseHist(histM1_fine);

	makeCopySave(histM1_rebin);
	normaliseHist(histM1_rebin);

	makeCopySave(histM2_loosetau);
	normaliseHist(histM2_loosetau);

	makeCopySave(histM2_fine);
	normaliseHist(histM2_fine);

	makeCopySave(histM2_rebin);
	normaliseHist(histM2_rebin);

	makeCopySave(histM_loosetau);
	normaliseHist(histM_loosetau);

	makeCopySave(histM1vsM2_loosetau);
	normaliseHist(histM1vsM2_loosetau);

	// Testing ones
	normaliseHist(histN_JPsi);
	normaliseHist(histN_JPsi_side_1to2p5);
*/

	// Do correlations plot by making m1*m1 first, then dividing m1vsm2 by m1*m1
	// Don't need to normalise m1timesm1, as histM1_side already normalised
	// turn off errors as correlated with numerator
	for(int a = 1; a <= nBinsX; a++){
		for (int b = 1; b <= nBinsX; b++){
			// signal region
			histM1timesM1->SetBinContent(a,b,histM->GetBinContent(a)*histM->GetBinContent(b));
			histM1timesM1->SetBinError(a,b,0);

			histM1timesM1_ModMass->SetBinContent(a,b,histM_ModMass->GetBinContent(a)*histM_ModMass->GetBinContent(b));
			histM1timesM1_ModMass->SetBinError(a,b,0);

			// control region A
			histM1timesM1_regionA->SetBinContent(a,b,histM_regionA->GetBinContent(a)*histM_regionA->GetBinContent(b));
			histM1timesM1_regionA->SetBinError(a,b,0);

			histM1timesM1_regionA_ModMass->SetBinContent(a,b,histM_regionA_ModMass->GetBinContent(a)*histM_regionA_ModMass->GetBinContent(b));
			histM1timesM1_regionA_ModMass->SetBinError(a,b,0);

			// control region B
			histM1timesM1_regionB->SetBinContent(a,b,histM_regionB->GetBinContent(a)*histM_regionB->GetBinContent(b));
			histM1timesM1_regionB->SetBinError(a,b,0);

			histM1timesM1_regionB_ModMass->SetBinContent(a,b,histM_regionB_ModMass->GetBinContent(a)*histM_regionB_ModMass->GetBinContent(b));
			histM1timesM1_regionB_ModMass->SetBinError(a,b,0);

/*			histM1timesM1_side_1to2p5->SetBinContent(a,b,histM_side_1to2p5->GetBinContent(a)*histM_side_1to2p5->GetBinContent(b));
			histM1timesM1_side_1to2p5->SetBinError(a,b,0);
			// histM1timesM1_side_1to2p5->SetBinError(a,b,sqrt(pow(histM_side_1to2p5->GetBinContent(b)*histM_side_1to2p5->GetBinError(a),2)
															// +pow(histM_side_1to2p5->GetBinContent(a)*histM_side_1to2p5->GetBinError(b),2)));

			// loose tau
			histM1timesM1_side_1to2p5_loosetau->SetBinContent(a,b,histM_side_1to2p5_loosetau->GetBinContent(a)*histM_side_1to2p5_loosetau->GetBinContent(b));
			histM1timesM1_side_1to2p5_loosetau->SetBinError(a,b,0);
			// histM1timesM1_side_1to2p5_loosetau->SetBinError(a,b,sqrt(pow(histM_side_1to2p5_loosetau->GetBinContent(b)*histM_side_1to2p5_loosetau->GetBinError(a),2)
															// +pow(histM_side_1to2p5_loosetau->GetBinContent(a)*histM_side_1to2p5_loosetau->GetBinError(b),2)));

			if (do1to1p5) {
				histM1timesM1_side_1to1p5->SetBinContent(a,b,histM_side_1to1p5->GetBinContent(a)*histM_side_1to1p5->GetBinContent(b));
				histM1timesM1_side_1to1p5->SetBinError(a,b,0);
				// histM1timesM1_side_1to1p5->SetBinError(a,b,sqrt(pow(histM_side_1to1p5->GetBinContent(b)*histM_side_1to1p5->GetBinError(a),2)
																// +pow(histM_side_1to1p5->GetBinContent(a)*histM_side_1to1p5->GetBinError(b),2)));
			}

			// histM1timesM1->SetBinError(a,b,sqrt(pow(histM->GetBinContent(b)*histM->GetBinError(a),2)
												// +pow(histM->GetBinContent(a)*histM->GetBinError(b),2)));
			// loose tau cand
			histM1timesM1_loosetau->SetBinContent(a,b,histM_loosetau->GetBinContent(a)*histM_loosetau->GetBinContent(b));
			histM1timesM1_loosetau->SetBinError(a,b,0);
			// histM1timesM1_loosetau->SetBinError(a,b,sqrt(pow(histM_loosetau->GetBinContent(b)*histM_loosetau->GetBinError(a),2)
												// +pow(histM_loosetau->GetBinContent(a)*histM_loosetau->GetBinError(b),2)));*/
		}
	}

	TH2D* histM1vsM2_correlations = (TH2D*)histM1vsM2->Clone("hM1vsM2_correlations");
	histM1vsM2_correlations->SetTitle("m(#mu_{1}-tk) vs m(#mu_{2}-tk) / m_{1} #times m_{1};m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]");
	histM1vsM2_correlations->Divide(histM1timesM1);

	TH2D* histM1vsM2_correlations_ModMass = (TH2D*)histM1vsM2_ModMass->Clone("hM1vsM2_correlations_ModMass");
	histM1vsM2_correlations_ModMass->SetTitle("m(#mu_{1}-tk) vs m(#mu_{2}-tk) / m_{1} #times m_{1};m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]");
	histM1vsM2_correlations_ModMass->Divide(histM1timesM1_ModMass);

	TH2D* histM1vsM2_correlations_regionA = (TH2D*)histM1vsM2_regionA->Clone("hM1vsM2_correlations_regionA");
	histM1vsM2_correlations_regionA->SetTitle("m(#mu_{1}-tk) vs m(#mu_{2}-tk) / m_{1} #times m_{1};m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]");
	histM1vsM2_correlations_regionA->Divide(histM1timesM1_regionA);

	TH2D* histM1vsM2_correlations_regionA_ModMass = (TH2D*)histM1vsM2_regionA_ModMass->Clone("hM1vsM2_correlations_regionA_ModMass");
	histM1vsM2_correlations_regionA_ModMass->SetTitle("m(#mu_{1}-tk) vs m(#mu_{2}-tk) / m_{1} #times m_{1};m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]");
	histM1vsM2_correlations_regionA_ModMass->Divide(histM1timesM1_regionA_ModMass);

	TH2D* histM1vsM2_correlations_regionB = (TH2D*)histM1vsM2_regionB->Clone("hM1vsM2_correlations_regionB");
	histM1vsM2_correlations_regionB->SetTitle("m(#mu_{1}-tk) vs m(#mu_{2}-tk) / m_{1} #times m_{1};m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]");
	histM1vsM2_correlations_regionB->Divide(histM1timesM1_regionB);

	TH2D* histM1vsM2_correlations_regionB_ModMass = (TH2D*)histM1vsM2_regionB_ModMass->Clone("hM1vsM2_correlations_regionB_ModMass");
	histM1vsM2_correlations_regionB_ModMass->SetTitle("m(#mu_{1}-tk) vs m(#mu_{2}-tk) / m_{1} #times m_{1};m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]");
	histM1vsM2_correlations_regionB_ModMass->Divide(histM1timesM1_regionB_ModMass);

/*	TH2D* histM1vsM2_correlations_side_1to2p5 = (TH2D*)histM1vsM2_side_1to2p5->Clone("hM1vsM2_correlations_side_1to2p5");
	histM1vsM2_correlations_side_1to2p5->SetTitle("m(#mu_{1}-tk) vs m(#mu_{2}-tk) / m(sideband) #times m(sideband), (soft tk p_{T} = 1 - 2.5 GeV);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]");
	histM1vsM2_correlations_side_1to2p5->Divide(histM1timesM1_side_1to2p5);

	// loose tau
	TH2D* histM1vsM2_correlations_side_1to2p5_loosetau = (TH2D*)histM1vsM2_side_1to2p5_loosetau->Clone("hM1vsM2_correlations_side_1to2p5_loosetau");
	histM1vsM2_correlations_side_1to2p5_loosetau->SetTitle("m(#mu_{1}-tk) vs m(#mu_{2}-tk) / m(sideband) #times m(sideband), (soft tk p_{T} = 1 - 2.5 GeV);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]");
	histM1vsM2_correlations_side_1to2p5_loosetau->Divide(histM1timesM1_side_1to2p5_loosetau);

	TH2D* histM1vsM2_correlations_side_1to1p5;
	if (do1to1p5) {
		histM1vsM2_correlations_side_1to1p5 = (TH2D*)histM1vsM2_side_1to1p5->Clone("hM1vsM2_correlations_side_1to1p5");
		histM1vsM2_correlations_side_1to1p5->SetTitle("m(#mu_{1}-tk) vs m(#mu_{2}-tk) / m(sideband) #times m(sideband), (soft tk p_{T} = 1 - 1.5 GeV);m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]");
		histM1vsM2_correlations_side_1to1p5->Divide(histM1timesM1_side_1to1p5);
	}


	// loose tau
	TH2D* histM1vsM2_correlations_loosetau = (TH2D*)histM1vsM2_loosetau->Clone("hM1vsM2_correlations_loosetau");
	histM1vsM2_correlations_loosetau->SetTitle("m(#mu_{1}-tk) vs m(#mu_{2}-tk) / m_{1} #times m_{1};m(#mu_{1}-tk) [GeV];m(#mu_{2}-tk) [GeV]");
	histM1vsM2_correlations_loosetau->Divide(histM1timesM1_loosetau);*/

	// Make 1D plots of unique bins from 2D correlation plot, label x axis
	int nUniqueBins = (nBinsX+1)*nBinsX/2.;

	TH1D* histCorr1D                      = new TH1D("hCorr1D",";Bin;Correlation coefficient", nUniqueBins, 0, nUniqueBins);
	TH1D* histCorr1D_ModMass              = new TH1D("hCorr1D_ModMass",";Bin;Correlation coefficient", nUniqueBins, 0, nUniqueBins);
	TH1D* histCorr1D_regionA              = new TH1D("hCorr1D_regionA",";Bin;Correlation coefficient", nUniqueBins, 0, nUniqueBins);
	TH1D* histCorr1D_regionA_ModMass      = new TH1D("hCorr1D_regionA_ModMass",";Bin;Correlation coefficient", nUniqueBins, 0, nUniqueBins);
	TH1D* histCorr1D_regionB              = new TH1D("hCorr1D_regionB",";Bin;Correlation coefficient", nUniqueBins, 0, nUniqueBins);
	TH1D* histCorr1D_regionB_ModMass      = new TH1D("hCorr1D_regionB_ModMass",";Bin;Correlation coefficient", nUniqueBins, 0, nUniqueBins);
	/*
	TH1D* histCorr1D_side_1to2p5          = new TH1D("hCorr1D_side_1to2p5",";Bin;Correlation coefficient",nUniqueBins,0,nUniqueBins);
	TH1D* histCorr1D_side_1to2p5_loosetau = new TH1D("hCorr1D_side_1to2p5_loosetau",";Bin;Correlation coefficient",nUniqueBins,0,nUniqueBins);
	TH1D* histCorr1D_side_1to1p5          = new TH1D("hCorr1D_side_1to1p5",";Bin;Correlation coefficient",nUniqueBins,0,nUniqueBins);
	TH1D* histCorr1D_loosetau             = new TH1D("hCorr1D_loosetau",";Bin;Correlation coefficient",nUniqueBins,0,nUniqueBins);
	*/

	int counter = 1;
	for (int i = 1; i <= nBinsX; i++) {
		for (int j = i; j <= nBinsX; j++) {
			std::string binLabel = "(" + boost::lexical_cast<std::string>(i) + "," + boost::lexical_cast<std::string>(j) + ")";

			histCorr1D->SetBinContent(counter,histM1vsM2_correlations->GetBinContent(i,j));
			histCorr1D->SetBinError(counter,histM1vsM2_correlations->GetBinError(i,j));
			histCorr1D->GetXaxis()->SetBinLabel(counter,binLabel.c_str());

			histCorr1D_ModMass->SetBinContent(counter,histM1vsM2_correlations_ModMass->GetBinContent(i,j));
			histCorr1D_ModMass->SetBinError(counter,histM1vsM2_correlations_ModMass->GetBinError(i,j));
			histCorr1D_ModMass->GetXaxis()->SetBinLabel(counter,binLabel.c_str());

			histCorr1D_regionA->SetBinContent(counter,histM1vsM2_correlations_regionA->GetBinContent(i,j));
			histCorr1D_regionA->SetBinError(counter,histM1vsM2_correlations_regionA->GetBinError(i,j));
			histCorr1D_regionA->GetXaxis()->SetBinLabel(counter,binLabel.c_str());

			histCorr1D_regionA_ModMass->SetBinContent(counter,histM1vsM2_correlations_regionA_ModMass->GetBinContent(i,j));
			histCorr1D_regionA_ModMass->SetBinError(counter,histM1vsM2_correlations_regionA_ModMass->GetBinError(i,j));
			histCorr1D_regionA_ModMass->GetXaxis()->SetBinLabel(counter,binLabel.c_str());

			histCorr1D_regionB->SetBinContent(counter,histM1vsM2_correlations_regionB->GetBinContent(i,j));
			histCorr1D_regionB->SetBinError(counter,histM1vsM2_correlations_regionB->GetBinError(i,j));
			histCorr1D_regionB->GetXaxis()->SetBinLabel(counter,binLabel.c_str());

			histCorr1D_regionB_ModMass->SetBinContent(counter,histM1vsM2_correlations_regionB_ModMass->GetBinContent(i,j));
			histCorr1D_regionB_ModMass->SetBinError(counter,histM1vsM2_correlations_regionB_ModMass->GetBinError(i,j));
			histCorr1D_regionB_ModMass->GetXaxis()->SetBinLabel(counter,binLabel.c_str());

/*
			histCorr1D_side_1to2p5->SetBinContent(counter,histM1vsM2_correlations_side_1to2p5->GetBinContent(i,j));
			histCorr1D_side_1to2p5->SetBinError(counter,histM1vsM2_correlations_side_1to2p5->GetBinError(i,j));
			histCorr1D_side_1to2p5->GetXaxis()->SetBinLabel(counter,binLabel.c_str());

			// loose tau
			histCorr1D_side_1to2p5_loosetau->SetBinContent(counter,histM1vsM2_correlations_side_1to2p5_loosetau->GetBinContent(i,j));
			histCorr1D_side_1to2p5_loosetau->SetBinError(counter,histM1vsM2_correlations_side_1to2p5_loosetau->GetBinError(i,j));
			histCorr1D_side_1to2p5_loosetau->GetXaxis()->SetBinLabel(counter,binLabel.c_str());

			if (do1to1p5) {
				histCorr1D_side_1to1p5->SetBinContent(counter,histM1vsM2_correlations_side_1to1p5->GetBinContent(i,j));
				histCorr1D_side_1to1p5->SetBinError(counter,histM1vsM2_correlations_side_1to1p5->GetBinError(i,j));
				histCorr1D_side_1to1p5->GetXaxis()->SetBinLabel(counter,binLabel.c_str());
			}


			// loose tau
			histCorr1D_loosetau->SetBinContent(counter,histM1vsM2_correlations_loosetau->GetBinContent(i,j));
			histCorr1D_loosetau->SetBinError(counter,histM1vsM2_correlations_loosetau->GetBinError(i,j));
			histCorr1D_loosetau->GetXaxis()->SetBinLabel(counter,binLabel.c_str());
*/

			counter++;
		}
	}



	// m1, mu1 pT in bins of m2
	// -------------------------------
/*	drawMassPlot("m(#mu_{1}-tk) in bins of m(#mu_{2}-tk) - signal region;m(#mu_{1}-tk) [GeV]; A.U.",
		histM1_0to1, histM1_1to2, histM1_2to3, histM1_3toInf, "M1_M2", directory, app);
	drawMassPlot("m(#mu_{1}-tk) in bins of m(#mu_{2}-tk) - sideband region (soft tk with p_{T} = 1 - 2.5 GeV);m(#mu_{1}-tk) [GeV]; A.U.",
		histM1_side_1to2p5_0to1, histM1_side_1to2p5_1to2, histM1_side_1to2p5_2to3, histM1_side_1to2p5_3toInf, "M1_M2_side_1to2p5", directory, app);
	if (do1to1p5)
		drawMassPlot("m(#mu_{1}-tk) in bins of m(#mu_{2}-tk) - sideband region (soft tk with p_{T} = 1 - 1.5 GeV);m(#mu_{1}-tk) [GeV]; A.U.",
			histM1_side_1to1p5_0to1, histM1_side_1to1p5_1to2, histM1_side_1to1p5_2to3, histM1_side_1to1p5_3toInf, "M1_M2_side_1to1p5", directory, app);

	if(doSignal){
		drawMassPlot("m(#mu_{1}-tk) in bins of m(#mu_{2}-tk) - MC truth;m(#mu_{1}-tk) [GeV]; A.U.",
			histM1_truth_0to1, histM1_truth_1to2, histM1_truth_2to3, histM1_truth_3toInf, "M1_M2_truth", directory, app);
		drawMassPlot("#mu_{1} p_{T} in bins of m(#mu_{2}-tk) - MC truth;#mu_{1} p_{T} [GeV]; A.U.",
			histMu1Pt_truth_0to1, histMu1Pt_truth_1to2, histMu1Pt_truth_2to3, histMu1Pt_truth_3toInf, "Mu1Pt_M2_truth", directory, app);
	}

	drawMassPlot("#mu_{1} p_{T} in bins of m(#mu_{2}-tk) - signal region;#mu_{1} p_{T} [GeV]; A.U.",
		histMu1Pt_0to1, histMu1Pt_1to2, histMu1Pt_2to3, histMu1Pt_3toInf, "Mu1Pt_M2", directory, app);
	drawMassPlot("#mu_{1} p_{T} in bins of m(#mu_{2}-tk) - sideband (soft tk p_{T} = 1 - 2.5 GeV);#mu_{1} p_{T} [GeV]; A.U.",
		histMu1Pt_side_1to2p5_0to1, histMu1Pt_side_1to2p5_1to2, histMu1Pt_side_1to2p5_2to3, histMu1Pt_side_1to2p5_3toInf, "Mu1Pt_M2_side_1to2p5", directory, app);
	drawMassPlot("#mu_{1} p_{T} in bins of m(#mu_{2}-tk) - sideband (soft tk p_{T} = 1 - 1.5 GeV);#mu_{1} p_{T} [GeV]; A.U.",
		histMu1Pt_side_1to1p5_0to1, histMu1Pt_side_1to1p5_1to2, histMu1Pt_side_1to1p5_2to3, histMu1Pt_side_1to1p5_3toInf, "Mu1Pt_M2_side_1to1p5", directory, app);
*/
	// 1D and 2D plots of masses & correlation coefficents:
	// ----------------------------------------------------

	// SIGNAL region
	drawHistAndSave(histM1, "HISTE", "M1", directory, app);
	drawHistAndSave(histM2, "HISTE", "M2", directory, app);
	drawHistAndSave(histM, "HISTE", "M", directory, app);
	drawHistAndSave(histM1vsM2, "colzTEXTE","M1vsM2", directory, app);
	drawHistAndSave(histM1timesM1, "colzTEXTE","M1timesM1", directory, app);
	drawHistAndSave(histM1vsM2_correlations, "colzTEXTE","M1vsM2_correlations", directory, app);
	drawHistAndSave(histCorr1D, "e1", "Correlations1D", directory, app);

	drawHistAndSave(histM1_ModMass, "HISTE", "M1_ModMass", directory, app);
	drawHistAndSave(histM2_ModMass, "HISTE", "M2_ModMass", directory, app);
	drawHistAndSave(histM_ModMass, "HISTE", "M_ModMass", directory, app);
	drawHistAndSave(histM1vsM2_ModMass, "colzTEXTE","M1vsM2_ModMass", directory, app);
	drawHistAndSave(histM1timesM1_ModMass, "colzTEXTE","M1timesM1_ModMass", directory, app);
	drawHistAndSave(histM1vsM2_correlations_ModMass, "colzTEXTE","M1vsM2_correlations_ModMass", directory, app);
	drawHistAndSave(histCorr1D_ModMass, "e1", "Correlations1D_ModMass", directory, app);

	// control region A
	drawHistAndSave(histM1_regionA, "HISTE", "M1_regionA", directory, app);
	drawHistAndSave(histM2_regionA, "HISTE", "M2_regionA", directory, app);
	drawHistAndSave(histM_regionA, "HISTE", "M_regionA", directory, app);
	drawHistAndSave(histM1vsM2_regionA, "colzTEXTE","M1vsM2_regionA", directory, app);
	drawHistAndSave(histM1timesM1_regionA, "colzTEXTE","M1timesM1_regionA", directory, app);
	drawHistAndSave(histM1vsM2_correlations_regionA, "colzTEXTE","M1vsM2_correlations_regionA", directory, app);
	drawHistAndSave(histCorr1D_regionA, "e1", "Correlations1D_regionA", directory, app);

	drawHistAndSave(histM1_regionA_ModMass, "HISTE", "M1_regionA_ModMass", directory, app);
	drawHistAndSave(histM2_regionA_ModMass, "HISTE", "M2_regionA_ModMass", directory, app);
	drawHistAndSave(histM_regionA_ModMass, "HISTE", "M_regionA_ModMass", directory, app);
	drawHistAndSave(histM1vsM2_regionA_ModMass, "colzTEXTE","M1vsM2_regionA_ModMass", directory, app);
	drawHistAndSave(histM1timesM1_regionA_ModMass, "colzTEXTE","M1timesM1_regionA_ModMass", directory, app);
	drawHistAndSave(histM1vsM2_correlations_regionA_ModMass, "colzTEXTE","M1vsM2_correlations_regionA_ModMass", directory, app);
	drawHistAndSave(histCorr1D_regionA_ModMass, "e1", "Correlations1D_regionA_ModMass", directory, app);

	// control region B
	drawHistAndSave(histM1_regionB, "HISTE", "M1_regionB", directory, app);
	drawHistAndSave(histM2_regionB, "HISTE", "M2_regionB", directory, app);
	drawHistAndSave(histM_regionB, "HISTE", "M_regionB", directory, app);
	drawHistAndSave(histM1vsM2_regionB, "colzTEXTE","M1vsM2_regionB", directory, app);
	drawHistAndSave(histM1timesM1_regionB, "colzTEXTE","M1timesM1_regionB", directory, app);
	drawHistAndSave(histM1vsM2_correlations_regionB, "colzTEXTE","M1vsM2_correlations_regionB", directory, app);
	drawHistAndSave(histCorr1D_regionB, "e1", "Correlations1D_regionB", directory, app);

	drawHistAndSave(histM1_regionB_ModMass, "HISTE", "M1_regionB_ModMass", directory, app);
	drawHistAndSave(histM2_regionB_ModMass, "HISTE", "M2_regionB_ModMass", directory, app);
	drawHistAndSave(histM_regionB_ModMass, "HISTE", "M_regionB_ModMass", directory, app);
	drawHistAndSave(histM1vsM2_regionB_ModMass, "colzTEXTE","M1vsM2_regionB_ModMass", directory, app);
	drawHistAndSave(histM1timesM1_regionB_ModMass, "colzTEXTE","M1timesM1_regionB_ModMass", directory, app);
	drawHistAndSave(histM1vsM2_correlations_regionB_ModMass, "colzTEXTE","M1vsM2_correlations_regionB_ModMass", directory, app);
	drawHistAndSave(histCorr1D_regionB_ModMass, "e1", "Correlations1D_regionB_ModMass", directory, app);

/*
	// SIDEBAND - 1 < tk pT < 2.5
	drawHistAndSave(histM1_side_1to2p5, "HISTE", "M1_side_1to2p5", directory, app);
	drawHistAndSave(histM2_side_1to2p5, "HISTE", "M2_side_1to2p5", directory, app);
	drawHistAndSave(histM_side_1to2p5, "HISTE", "M_side_1to2p5", directory, app);
	drawHistAndSave(histM1vsM2_side_1to2p5, "colzTEXTE","M1vsM2_side_1to2p5", directory, app);
	drawHistAndSave(histM1timesM1_side_1to2p5, "colzTEXTE","M1timesM1_side_1to2p5", directory, app);
	drawHistAndSave(histM1vsM2_correlations_side_1to2p5, "colzTEXTE","M1vsM2_correlations_side_1to2p5", directory, app);
	drawHistAndSave(histCorr1D_side_1to2p5, "e1", "Correlations1D_side_1to2p5", directory, app);

	// SIDEBAND - 1 < tk pT < 2.5 - loose tau cand
	drawHistAndSave(histM1_side_1to2p5_loosetau, "HISTE", "M1_side_1to2p5_loosetau", directory, app);
	drawHistAndSave(histM2_side_1to2p5_loosetau, "HISTE", "M2_side_1to2p5_loosetau", directory, app);
	drawHistAndSave(histM_side_1to2p5_loosetau, "HISTE", "M_side_1to2p5_loosetau", directory, app);
	drawHistAndSave(histM1vsM2_side_1to2p5_loosetau, "colzTEXTE","M1vsM2_side_1to2p5_loosetau", directory, app);
	drawHistAndSave(histM1timesM1_side_1to2p5_loosetau, "colzTEXTE","M1timesM1_side_1to2p5_loosetau", directory, app);
	drawHistAndSave(histM1vsM2_correlations_side_1to2p5_loosetau, "colzTEXTE","M1vsM2_correlations_side_1to2p5_loosetau", directory, app);
	drawHistAndSave(histCorr1D_side_1to2p5_loosetau, "e1", "Correlations1D_side_1to2p5_loosetau", directory, app);

	// SIDEBAND - 1 < tk pT < 1.5
	if (do1to1p5) {
		drawHistAndSave(histM1_side_1to1p5, "HISTE", "M1_side_1to1p5", directory, app);
		drawHistAndSave(histM2_side_1to1p5, "HISTE", "M2_side_1to1p5", directory, app);
		drawHistAndSave(histM_side_1to1p5, "HISTE", "M_side_1to1p5", directory, app);
		drawHistAndSave(histM1vsM2_side_1to1p5, "colzTEXTE","M1vsM2_side_1to1p5", directory, app);
		drawHistAndSave(histM1timesM1_side_1to1p5, "colzTEXTE","M1timesM1_side_1to1p5", directory, app);
		drawHistAndSave(histM1vsM2_correlations_side_1to1p5, "colzTEXTE","M1vsM2_correlations_side_1to1p5", directory, app);
		drawHistAndSave(histCorr1D_side_1to1p5, "e1", "Correlations1D_side_1to1p5", directory, app);
	}


	// SIGNAL region - loose tau
	drawHistAndSave(histM1_loosetau, "HISTE", "M1_loosetau", directory, app);
	drawHistAndSave(histM2_loosetau, "HISTE", "M2_loosetau", directory, app);
	drawHistAndSave(histM_loosetau, "HISTE", "M_loosetau", directory, app);
	drawHistAndSave(histM1vsM2_loosetau, "colzTEXTE","M1vsM2_loosetau", directory, app);
	drawHistAndSave(histM1timesM1_loosetau, "colzTEXTE","M1timesM1_loosetau", directory, app);
	drawHistAndSave(histM1vsM2_correlations_loosetau, "colzTEXTE","M1vsM2_correlations_loosetau", directory, app);
	drawHistAndSave(histCorr1D_loosetau, "e1", "Correlations1D_loosetau", directory, app);


	// Extra tracks around mu2
	drawHistAndSave(histM1_Ntk2_2, "HISTE", "M1_Ntk2_2", directory, app);
	drawHistAndSave(histM1_Ntk2_3, "HISTE", "M1_Ntk2_3", directory, app);
	drawHistAndSave(histM1_Ntk2_2or3, "HISTE", "M1_Ntk2_2or3", directory, app);
	drawHistAndSave(histM1_Ntk2_4, "HISTE", "M1_Ntk2_4", directory, app);

	// Extra tracks around mu2 - looser IP on 1-prong track (IP option 2)
	drawHistAndSave(histM1_Ntk2_2_loosetau, "HISTE", "M1_Ntk2_2_loosetau", directory, app);
	drawHistAndSave(histM1_Ntk2_3_loosetau, "HISTE", "M1_Ntk2_3_loosetau", directory, app);
	drawHistAndSave(histM1_Ntk2_2or3_loosetau, "HISTE", "M1_Ntk2_2or3_loosetau", directory, app);
	drawHistAndSave(histM1_Ntk2_4_loosetau, "HISTE", "M1_Ntk2_4_loosetau", directory, app);

	// Testing hists
	// -------------

	drawHistAndSave(histN_JPsi, "HISTE", "N_JPsi", directory, app);
	drawHistAndSave(histN_JPsi_side_1to2p5, "HISTE", "N_JPsi_side_1to2p5", directory, app);

	drawHistAndSave(histN_2p5, "HISTE", "N_2p5", directory, app);
	drawHistAndSave(histN_2p5_OS, "HISTE", "N_2p5_OS", directory, app);

	drawHistAndSave(histN_2p5_looseIP, "HISTE", "N_2p5_looseIP", directory, app);
	drawHistAndSave(histN_2p5_OS_looseIP, "HISTE", "N_2p5_OS_looseIP", directory, app);

*/
	// ARC study region
	drawHistAndSave(histMu1PT_noTk2, "HISTE", "Mu1PT_noTk2", directory, app);
	drawHistAndSave(histMu2PT_noTk2, "HISTE", "Mu2PT_noTk2", directory, app);
	drawHistAndSave(histMuPT_noTk2, "HISTE", "MuPT_noTk2", directory, app);
	drawHistAndSave(histMu1Eta_noTk2, "HISTE", "Mu1Eta_noTk2", directory, app);
	drawHistAndSave(histMu2Eta_noTk2, "HISTE", "Mu2Eta_noTk2", directory, app);
	drawHistAndSave(histMuEta_noTk2, "HISTE", "MuEta_noTk2", directory, app);
	drawHistAndSave(histTk1PT_noTk2, "HISTE", "Tk1PT_noTk2", directory, app);
	drawHistAndSave(histTk2PT_noTk2, "HISTE", "Tk2PT_noTk2", directory, app);
	drawHistAndSave(histTkPT_noTk2, "HISTE", "TkPT_noTk2", directory, app);
	drawHistAndSave(histTk1Eta_noTk2, "HISTE", "Tk1Eta_noTk2", directory, app);
	drawHistAndSave(histTk2Eta_noTk2, "HISTE", "Tk2Eta_noTk2", directory, app);
	drawHistAndSave(histTkEta_noTk2, "HISTE", "TkEta_noTk2", directory, app);
	drawHistAndSave(histDRMuTk1_noTk2, "HISTE", "DRMuTk1_noTk2", directory, app);
	drawHistAndSave(histDRMuTk2_noTk2, "HISTE", "DRMuTk2_noTk2", directory, app);
	drawHistAndSave(histDRMuTk_noTk2, "HISTE", "DRMuTk_noTk2", directory, app);
	drawHistAndSave(histMassMuTk1_noTk2, "HISTE", "MassMuTk1_noTk2", directory, app);
	drawHistAndSave(histMassMuTk2_noTk2, "HISTE", "MassMuTk2_noTk2", directory, app);
	drawHistAndSave(histMassMuTk_noTk2, "HISTE", "MassMuTk_noTk2", directory, app);
	drawHistAndSave(histMassMuTkModMass_noTk2, "HISTE", "MassMuTk_noTk2", directory, app);


	//////////////////////////
	// Write hists to file //
	//////////////////////////

	// signal region
	histM1->Write("",TObject::kOverwrite);
	histM2->Write("",TObject::kOverwrite);
	histM->Write("",TObject::kOverwrite);
	histM1timesM1->Write("",TObject::kOverwrite);
	histM1vsM2->Write("",TObject::kOverwrite);
	histM1vsM2_correlations->Write("",TObject::kOverwrite);
	histCorr1D->Write("",TObject::kOverwrite);

	histM1_ModMass->Write("",TObject::kOverwrite);
	histM2_ModMass->Write("",TObject::kOverwrite);
	histM_ModMass->Write("",TObject::kOverwrite);
	histM1timesM1_ModMass->Write("",TObject::kOverwrite);
	histM1vsM2_ModMass->Write("",TObject::kOverwrite);
	histM1vsM2_correlations_ModMass->Write("",TObject::kOverwrite);
	histCorr1D_ModMass->Write("",TObject::kOverwrite);

	// control region A
	histM1_regionA->Write("", TObject::kOverwrite);
	histM2_regionA->Write("", TObject::kOverwrite);
	histM_regionA->Write("", TObject::kOverwrite);
	histM1vsM2_regionA->Write("", TObject::kOverwrite);
	histM1timesM1_regionA->Write("", TObject::kOverwrite);
	histM1vsM2_correlations_regionA->Write("", TObject::kOverwrite);
	histCorr1D_regionA->Write("", TObject::kOverwrite);

	histM1_regionA_ModMass->Write("", TObject::kOverwrite);
	histM2_regionA_ModMass->Write("", TObject::kOverwrite);
	histM_regionA_ModMass->Write("", TObject::kOverwrite);
	histM1vsM2_regionA_ModMass->Write("", TObject::kOverwrite);
	histM1timesM1_regionA_ModMass->Write("", TObject::kOverwrite);
	histM1vsM2_correlations_regionA_ModMass->Write("", TObject::kOverwrite);
	histCorr1D_regionA_ModMass->Write("", TObject::kOverwrite);

	// control region B
	histM1_regionB->Write("", TObject::kOverwrite);
	histM2_regionB->Write("", TObject::kOverwrite);
	histM_regionB->Write("", TObject::kOverwrite);
	histM1vsM2_regionB->Write("", TObject::kOverwrite);
	histM1timesM1_regionB->Write("", TObject::kOverwrite);
	histM1vsM2_correlations_regionB->Write("", TObject::kOverwrite);
	histCorr1D_regionB->Write("", TObject::kOverwrite);

	histM1_regionB_ModMass->Write("", TObject::kOverwrite);
	histM2_regionB_ModMass->Write("", TObject::kOverwrite);
	histM_regionB_ModMass->Write("", TObject::kOverwrite);
	histM1vsM2_regionB_ModMass->Write("", TObject::kOverwrite);
	histM1timesM1_regionB_ModMass->Write("", TObject::kOverwrite);
	histM1vsM2_correlations_regionB_ModMass->Write("", TObject::kOverwrite);
	histCorr1D_regionB_ModMass->Write("", TObject::kOverwrite);

	//_ModMass arc study region
	histMu1PT_noTk2->Write("", TObject::kOverwrite);
	histMu2PT_noTk2->Write("", TObject::kOverwrite);
	histMuPT_noTk2->Write("", TObject::kOverwrite);
	histMu1Eta_noTk2->Write("", TObject::kOverwrite);
	histMu2Eta_noTk2->Write("", TObject::kOverwrite);
	histMuEta_noTk2->Write("", TObject::kOverwrite);
	histTk1PT_noTk2->Write("", TObject::kOverwrite);
	histTk2PT_noTk2->Write("", TObject::kOverwrite);
	histTkPT_noTk2->Write("", TObject::kOverwrite);
	histTk1Eta_noTk2->Write("", TObject::kOverwrite);
	histTk2Eta_noTk2->Write("", TObject::kOverwrite);
	histTkEta_noTk2->Write("", TObject::kOverwrite);
	histDRMuTk1_noTk2->Write("", TObject::kOverwrite);
	histDRMuTk2_noTk2->Write("", TObject::kOverwrite);
	histDRMuTk_noTk2->Write("", TObject::kOverwrite);
	histMassMuTk1_noTk2->Write("", TObject::kOverwrite);
	histMassMuTk2_noTk2->Write("", TObject::kOverwrite);
	histMassMuTk_noTk2->Write("", TObject::kOverwrite);
	histMassMuTkModMass_noTk2->Write("", TObject::kOverwrite);

/*
	histM1_Ntk2_2->Write("",TObject::kOverwrite);
	histM1_Ntk2_3->Write("",TObject::kOverwrite);
	histM1_Ntk2_2or3->Write("",TObject::kOverwrite);
	histM1_Ntk2_4->Write("",TObject::kOverwrite);
	histM1_fine->Write("",TObject::kOverwrite);
	histM2_fine->Write("",TObject::kOverwrite);
	if (doSignal){
		histM1_truth_0to1->Write("",TObject::kOverwrite);
		histM1_truth_1to2->Write("",TObject::kOverwrite);
		histM1_truth_2to3->Write("",TObject::kOverwrite);
		histM1_truth_3toInf->Write("",TObject::kOverwrite);
	}
	histM1_0to1->Write("",TObject::kOverwrite);
	histM1_1to2->Write("",TObject::kOverwrite);
	histM1_2to3->Write("",TObject::kOverwrite);
	histM1_3toInf->Write("",TObject::kOverwrite);

	histM1_side_1to2p5_0to1->Write("",TObject::kOverwrite);
	histM1_side_1to2p5_1to2->Write("",TObject::kOverwrite);
	histM1_side_1to2p5_2to3->Write("",TObject::kOverwrite);
	histM1_side_1to2p5_3toInf->Write("",TObject::kOverwrite);

	if (do1to1p5) {
		histM1_side_1to1p5_0to1->Write("",TObject::kOverwrite);
		histM1_side_1to1p5_1to2->Write("",TObject::kOverwrite);
		histM1_side_1to1p5_2to3->Write("",TObject::kOverwrite);
		histM1_side_1to1p5_3toInf->Write("",TObject::kOverwrite);
	}

	histMu1Pt_0to1->Write("",TObject::kOverwrite);
	histMu1Pt_1to2->Write("",TObject::kOverwrite);
	histMu1Pt_2to3->Write("",TObject::kOverwrite);
	histMu1Pt_3toInf->Write("",TObject::kOverwrite);

	if (doSignal){
		histMu1Pt_truth_0to1->Write("",TObject::kOverwrite);
		histMu1Pt_truth_1to2->Write("",TObject::kOverwrite);
		histMu1Pt_truth_2to3->Write("",TObject::kOverwrite);
		histMu1Pt_truth_3toInf->Write("",TObject::kOverwrite);
	}

	histMu1Pt_side_1to2p5_0to1->Write("",TObject::kOverwrite);
	histMu1Pt_side_1to2p5_1to2->Write("",TObject::kOverwrite);
	histMu1Pt_side_1to2p5_2to3->Write("",TObject::kOverwrite);
	histMu1Pt_side_1to2p5_3toInf->Write("",TObject::kOverwrite);

	if (do1to1p5) {
		histMu1Pt_side_1to1p5_0to1->Write("",TObject::kOverwrite);
		histMu1Pt_side_1to1p5_1to2->Write("",TObject::kOverwrite);
		histMu1Pt_side_1to1p5_2to3->Write("",TObject::kOverwrite);
		histMu1Pt_side_1to1p5_3toInf->Write("",TObject::kOverwrite);
	}

	histM1_side_1to2p5->Write("",TObject::kOverwrite);
	histM2_side_1to2p5->Write("",TObject::kOverwrite);
	histM_side_1to2p5->Write("",TObject::kOverwrite);
	histM1timesM1_side_1to2p5->Write("",TObject::kOverwrite);
	histM1vsM2_side_1to2p5->Write("",TObject::kOverwrite);
	histM1vsM2_correlations_side_1to2p5->Write("",TObject::kOverwrite);
	histCorr1D_side_1to2p5->Write("",TObject::kOverwrite);

	// side - loose tau
	histM1_side_1to2p5_loosetau->Write("",TObject::kOverwrite);
	histM2_side_1to2p5_loosetau->Write("",TObject::kOverwrite);
	histM_side_1to2p5_loosetau->Write("",TObject::kOverwrite);
	histM1vsM2_side_1to2p5_loosetau->Write("",TObject::kOverwrite);
	histM1timesM1_side_1to2p5_loosetau->Write("",TObject::kOverwrite);
	histM1vsM2_correlations_side_1to2p5_loosetau->Write("",TObject::kOverwrite);
	histCorr1D_side_1to2p5_loosetau->Write("",TObject::kOverwrite);


	if (do1to1p5) {
		histM1_side_1to1p5->Write("",TObject::kOverwrite);
		histM2_side_1to1p5->Write("",TObject::kOverwrite);
		histM_side_1to1p5->Write("",TObject::kOverwrite);
		histM1vsM2_side_1to1p5->Write("",TObject::kOverwrite);
		histM1vsM2_correlations_side_1to1p5->Write("",TObject::kOverwrite);
		histCorr1D_side_1to1p5->Write("",TObject::kOverwrite);
	}

	// signal - loose tau
	histM1_loosetau->Write("",TObject::kOverwrite);
	histM2_loosetau->Write("",TObject::kOverwrite);
	histM_loosetau->Write("",TObject::kOverwrite);
	histM1vsM2_loosetau->Write("",TObject::kOverwrite);
	histM1timesM1_loosetau->Write("",TObject::kOverwrite);
	histM1vsM2_correlations_loosetau->Write("",TObject::kOverwrite);
	histCorr1D_loosetau->Write("",TObject::kOverwrite);
	histM1_Ntk2_2_loosetau->Write("",TObject::kOverwrite);
	histM1_Ntk2_3_loosetau->Write("",TObject::kOverwrite);
	histM1_Ntk2_2or3_loosetau->Write("",TObject::kOverwrite);
	histM1_Ntk2_4_loosetau->Write("",TObject::kOverwrite);

	histN_JPsi->Write("", TObject::kOverwrite);
	histN_JPsi_side_1to2p5->Write("", TObject::kOverwrite);

	histN_2p5->Write("", TObject::kOverwrite);
	histN_2p5_OS->Write("", TObject::kOverwrite);

	histN_2p5_looseIP->Write("", TObject::kOverwrite);
	histN_2p5_OS_looseIP->Write("", TObject::kOverwrite);
*/
	outFile->Close();

	delete treeReader;
}
