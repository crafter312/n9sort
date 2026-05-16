/* New and improved version of the wood TTree output class (haha, just got the joke)
 * written by Henry Webb (h.s.webb@wustl.edu), originally written by Robert Charity
 * and then modified by Johnathan Phillips.
 * 
 * This new version is highly flexible and modularized. Instead of having a million
 * separate cases for the different numbers of particles and different reconstructions,
 * one simply initializes the `frags` list of pointers of type `GenericOut`, an output
 * class for which there is a child for each detector type supported by the sort code.
 * This initialization is done in correl2.cpp. One then retireves the relevant array
 * and array size of `solution` objects for the fragments and then passes it to this
 * class for transfer to the output objects. Direct output of these classes is via
 * ROOT dictionaries generated at compile time.
 * 
 * When adding a new detector type, one must update five things: add a new enum
 * option to `parType::detType`, add a new child class of `GenericOut` in wood.h,
 * add a new fragment vector and map in wood.h, add a new branch initialization
 * in the wood constructor, and update the switch statements in `wood::initFrag` and
 * `wood::loadFrag`.
 */

#include "wood.h"

#include <stdexcept>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// NOTE: correl2 object's masks must be set correctly before instantiating this class
// in order for initialization to proceed correctly for the desired reconstruction!
wood::wood(correl2& correl, string name, TDirectory* dir, bool gamma0) : gamma(gamma0) {

#ifdef ENABLE_DEBUG
	cout << "wood::wood 1" << endl;
#endif

	correl.initWood(this);
	
#ifdef ENABLE_DEBUG
	cout << "wood::wood 2" << endl;
#endif
	
	dir->cd();
	t = new TTree(name.c_str(), name.c_str());
	if (gobbiFrags.size() > 0) t->Branch("gobbiFrags", &gobbiFrags);
	if (texNeutFrags.size() > 0) t->Branch("texNeutFrags", &texNeutFrags);
	if (s800Frags.size() > 0) t->Branch("s800Frags", &s800Frags);

#ifdef ENABLE_DEBUG
	cout << "wood::wood 3" << endl;
#endif

	// Non-fragment branches
	t->Branch("Erel", &Erel);
	t->Branch("Ex", &Ex);
	t->Branch("Vcm", &Vcm);
	t->Branch("thetaCM", &thetaCM);
	t->Branch("cos_thetaH", &cos_thetaH);

	if (gamma) {
		t->Branch("Ngamma", &Ngamma);
		t->Branch("Ngamma_Select", &Ngamma_Select);
		t->Branch("Egamma", Egamma, "Egamma[15]/F");
		t->Branch("Egamma_Select", Egamma_Select, "Egamma_Select[15]/F");
		t->Branch("Tgamma", Tgamma, "Tgamma[15]/F");
		t->Branch("Tgamma_Select", Tgamma_Select, "Tgamma_Select[15]/F");
		t->Branch("Chgamma", Chgamma, "Chgamma[15]/I");
		t->Branch("Chgamma_Select", Chgamma_Select, "Chgamma_Select[15]/I");
	}

	t->Branch("runnum", &runnum);
	t->Branch("beamZ", &beamZ);
	
#ifdef ENABLE_DEBUG
	cout << "wood::wood 4" << endl;
#endif
	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

wood::~wood() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void wood::initFrag(size_t i, parType::detType det) {

#ifdef ENABLE_DEBUG
	cout << "wood::initFrag 1" << endl;
#endif

	switch (det) {
		case parType::detType::Gobbi:
			gobbiFrags.emplace_back();
			frags.push_back(&gobbiFrags.back());
			gobbiIndexMap[i] = gobbiFrags.size() - 1;
			break;
		case parType::detType::TexNeut:
			texNeutFrags.emplace_back();
			frags.push_back(&texNeutFrags.back());
			texNeutIndexMap[i] = texNeutFrags.size() - 1;
			break;
		case parType::detType::S800:
			s800Frags.emplace_back();
			frags.push_back(&s800Frags.back());
			s800IndexMap[i] = s800Frags.size() - 1;
			break;
		default:
			throw invalid_argument(string("wood initialization case for detector "));// + parType::ToString(det) + string(" does not exist"));
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void wood::loadFrag(size_t i, parType::detType det, solution* frag) {

	// First, set the common variables
	loadFragCommon(i, frag);

	// Then, set the detector specific variables
	GobbiOut* gobbiOut;
	TexNeutOut* texNeutOut;
	S800Out* s800Out;
	switch (det) {
		case parType::detType::Gobbi:
			if (gobbiIndexMap.find(i) == gobbiIndexMap.end()) break;
			gobbiOut = &gobbiFrags[gobbiIndexMap[i]];
			gobbiOut->itele = frag->itele;
			gobbiOut->id = frag->iCsI;
			gobbiOut->ifront = frag->ifront;
			gobbiOut->iback = frag->iback;
			gobbiOut->time = frag->CsITime;
			gobbiOut->denergy_R = frag->denergyR;
			gobbiOut->energy_p_R = frag->energyR;
			return;
		case parType::detType::TexNeut:
			if (texNeutIndexMap.find(i) == texNeutIndexMap.end()) break;
			texNeutOut = &texNeutFrags[texNeutIndexMap[i]];
			return;
		case parType::detType::S800:
			if (s800IndexMap.find(i) == s800IndexMap.end()) break;
			s800Out = &s800Frags[s800IndexMap[i]];
			s800Out->theta_s800 = frag->theta_s800;
			s800Out->phi_s800 = frag->phi_s800;
			return;
		default:
			throw invalid_argument(string("wood fragment output for detector "));// + parType::ToString(det) + string(" is not defined"));
	}

	throw invalid_argument("Something went wrong! Fragment output vector of requested detector type does not contain fragment of requested index.");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// NOTE: This version of the code doesn't use CAESAR or any other gamma detector.
// As such, I do not know what this function should look like, and it will have
// to be completed in the future case where gammas are required.
void wood::setGammas() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void wood::Fill(float _Erel, float _Ex, float _Vcm, float _thetaCM, float _cos_thetaH, int _runnum, int _beamZ) {
	Erel = _Erel;
	Ex = _Ex;
	Vcm = _Vcm;
	thetaCM = _thetaCM;
	cos_thetaH = _cos_thetaH;
	runnum = _runnum;
	beamZ = _beamZ;

	t->Fill();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void wood::loadFragCommon(size_t i, solution* frag) {
	if (i >= frags.size()) throw invalid_argument("Provided fragment index larger than expected number of fragments");

	GenericOut& fragOut = *frags[i];
	fragOut.M[0] = frag->Mvect[0];
	fragOut.M[1] = frag->Mvect[1];
	fragOut.M[2] = frag->Mvect[2];
	fragOut.et = frag->energyTot;
	fragOut.energy_p = frag->Ekin;
	fragOut.theta = frag->theta;
	fragOut.phi = frag->phi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



