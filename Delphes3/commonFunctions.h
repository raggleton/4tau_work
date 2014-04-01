#include <vector>
#include <algorithm>
#include <string>
#include <sstream>

using std::cout;
using std::endl;

/**
 * This header contains common functiosn for all my Delphes analysis scripts
 *
 * Robin Aggleton 2014
 */


// For string splitting
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

// For sorting track vectors by pT
// Ideally we'd use the templated methods in classes/SortableObject.h ...
bool sortTracksByPT(Track* a, Track* b){ 
	return (a->PT) > (b->PT); 
}

// From the daughters of 2 taus, decide which is track and which is mu
bool assignMuonAndTrack(GenParticle* &mu, GenParticle* &tk, GenParticle &a, GenParticle &b){
	tk = 0;
	mu = 0;

	bool aIsMu = fabs(a.PID) == 13;
	bool bIsMu = fabs(b.PID) == 13;

	// Now look at both a and b
	if (aIsMu && !bIsMu) {
		mu = &a;
		tk = &b;
		return true;
	}

	if (!aIsMu && bIsMu) {
		mu = &b;
		tk = &a;
		return true;
	}

	// the muon is the mu with the higest pT, the other mu becomes a track
	if (aIsMu && bIsMu){
		if (a.PT > b.PT) {
			mu = &a;
			tk = &b;
			return true;
		} else {
			mu = &b;
			tk = &a;
			return true;
		}
	}
	// If it gets to here, then neither a nor b is a muon - trouble!
	return false;
}

// Get the 3 correct daughters of the tau


std::vector<GenParticle*> getTauDaughters(TClonesArray *branchAll, GenParticle *tau) { 

	// cout << "get tau daughters" << endl;

	std::vector<GenParticle*> tauDaughters;

	bool foundProduct = false;
	while(!foundProduct){
		if(tau->D1 == tau->D2) { // tau -> tau
			tau = (GenParticle*) branchAll->At(tau->D1);
		} else if (tau->D2-tau->D1 == 1) { //tau->tau+gamma
			if ( fabs(((GenParticle*) branchAll->At(tau->D1))->PID) == 22) {
				tau = (GenParticle*) branchAll->At(tau->D2);
			} else{
				tau = (GenParticle*) branchAll->At(tau->D1);
			}
		} else {
			for(int i = tau->D1; i <= tau->D2; i++){
				tauDaughters.push_back((GenParticle*) branchAll->At(i));
				// tauDaughters.push_back(i);
			}
			if (tauDaughters.size() ==3){
				// cout << "3 products" << endl;
				foundProduct = true;
				return tauDaughters;
			}
		}
	}
}

/**
 * From tau, get the final, stable, charged decay product
 * @param  branchAll [the branch of ALL gen particles, to loop through decay chain]
 * @param  tau       [tau to get decay product form]
 * @return           [returns charged, stable, decay product of param tau]
 */
GenParticle* getChargedObject(TClonesArray* branchAll, GenParticle* tau) { 

	std::vector<GenParticle*> history; // hold all unique particles in the decay chain (stores event posiiton number)
	std::vector<GenParticle*> current = getTauDaughters(branchAll, tau);
	std::vector<GenParticle*> next; // holds decay products, not nec. all unique
	GenParticle *prong(0); // holds final charged product
	int nCharged = 0; // count charged particles - shouldn't have more than 1!
	bool foundOne = false;

	while (current.size()>0){ // if current > 0 we haven't exhausted all the particles
		for (unsigned a = 0; a < current.size(); a++){
			// Check 1 - is this already in current?
			// Could probably do more efficiently using the unique function on std::vector
			bool alreadyDone = false;

			for (unsigned b = 0; b < a; b++){
				if ((current[a] == current[b]) && (a != 0)) {
					// cout << "--Found a duplicate" << endl;
					alreadyDone = true;
				}
			}

			// Check 2 - is this already in history?
			if (!alreadyDone){
				// cout << "-not in current" << endl;
				for (unsigned c = 0; c < history.size(); c++){
					if ((current[a] == history[c]) && (c!=0)){
						// cout << "--Found a dupicate in history" << endl;
						alreadyDone = true;
					}
				}
			}

			// Check 3 - is this final state?
			// cout << alreadyDone << endl;
			if (!alreadyDone){
				// cout << "-not in history" << endl;

				// Check if final state. Either status 1, or D1 == D2 == -1
				if (current[a]->Status == 1) {

					// Check if charged
					if (current[a]->Charge != 0) {
						// cout << "FINAL PRONG " << current[a]->PID << endl;
						prong = current[a];
						foundOne = true;
						nCharged++;
					}

					history.push_back(current[a]);
					// cout << "pushed into history" << endl;
				} else {
					// Load its daughters no. into next
					for (int d = current[a]->D1; d <= current[a]->D2; d++) {
						// cout << "pushing" << endl;
						next.push_back((GenParticle*) branchAll->At(d));
					}

					// Load it into history
					history.push_back(current[a]);
				}

			}// end of alreadyDone
		} // end of loop over current

		// Clear current - don't need to do as = assignment auto does this
		// current.clear();
		// Copy next into current
		current = next;
		// Empty next
		next.clear();
	} // end of while(!done)

	// Just in case
	history.clear();
	current.clear();
	next.clear();

	if (!foundOne) 
		cout << "No prongs!!!" << endl;

	if (nCharged == 1){
		// cout << "PRONG: " << prong->PID << endl;
		return prong;
	} else {
		// cout << "Damn, more than 1 prong!!!! Got " << nCharged << endl;
		return NULL;
	}

}

/**
 * Draw a histogram and save it to PDF, saves file as directory/filename_<delphes setup from directory>_app.pdf
 * @param h         [hist to draw (TObject* to handle THStacks, which don't inherit from TH1)]
 * @param drawOpt   [options for drawing hist]
 * @param filename  [main filename (e.g. Mu1Pt)]
 * @param directory [directory for output PDF]
 * @param app       [appendage eg muRand, sig]
 */
void drawHistAndSave(TObject* h, std::string drawOpt, std::string filename, std::string directory, std::string app){
	TCanvas c;
	h->Draw(drawOpt.c_str());

	std::vector<std::string> elems2;
	split(directory, '_', elems2);
	std::string delph = elems2[elems2.size()-1];
	
	c.SaveAs((directory+"/"+filename+"_"+delph+"_"+app+".pdf").c_str());
}

/**
 * Draws mass correlation plot, saves file as directory/filename_<delphes setup from directory>_app.pdf
 * @param title         [description]
 * @param histM1_0to1   [description]
 * @param histM1_1to2   [description]
 * @param histM1_2to3   [description]
 * @param histM1_3toInf [description]
 * @param filename      [description]
 * @param directory     [description]
 * @param app           [description]
 */
void drawMassPlot(std::string title, TH1* histM1_0to1, TH1* histM1_1to2, TH1* histM1_2to3, TH1* histM1_3toInf, std::string filename, std::string directory, std::string app){
	TCanvas c;
	THStack histM1_M2("hM1_M2",title.c_str());
	histM1_0to1->SetLineColor(kBlack);
	if (histM1_0to1->Integral() != 0)
		histM1_0to1->Scale(1./histM1_0to1->Integral());
	histM1_M2.Add(histM1_0to1);
	
	histM1_1to2->SetLineColor(kRed);
	if (histM1_1to2->Integral() != 0)
		histM1_1to2->Scale(1./histM1_1to2->Integral());
	histM1_M2.Add(histM1_1to2);

	histM1_2to3->SetLineColor(kGreen);
	if (histM1_2to3->Integral() != 0)
		histM1_2to3->Scale(1./histM1_2to3->Integral());
	histM1_M2.Add(histM1_2to3);

	histM1_3toInf->SetLineColor(kBlue);
	if (histM1_3toInf->Integral() != 0)
		histM1_3toInf->Scale(1./histM1_3toInf->Integral());
	histM1_M2.Add(histM1_3toInf);
	histM1_M2.Draw("nostack,HISTE");

	TLegend leg(0.7,0.7,0.9,0.9);
	leg.AddEntry(histM1_0to1,"m_{2} = 0-1 GeV","l");
	leg.AddEntry(histM1_1to2,"m_{2} = 1-2 GeV","l");
	leg.AddEntry(histM1_2to3,"m_{2} = 2-3 GeV","l");
	leg.AddEntry(histM1_3toInf,"m_{2} > 3 GeV","l");
	leg.Draw();
	
	std::vector<std::string> elems2;
	split(directory, '_', elems2);
	std::string delph = elems2[elems2.size()-1];
	
	c.SaveAs((directory+"/"+filename+"_"+delph+"_"+app+".pdf").c_str());
}
