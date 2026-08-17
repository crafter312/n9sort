#ifndef wood_H
#define wood_H

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

#include <cmath>
#include <deque>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <TDirectory.h>
#include <TTree.h>

#include "correl2.h"
#include "parType.h"
#include "solution.h"

class wood {

public:

	/******** OUTPUT CLASSES ********/
	// Notes:
	// The GenericOut constructor being protected prevents (I think)
	// this class from being instantiated on its own. Thus, one must
	// define a dedicated child class of GenericOut in order to use
	// it for a specific detector.

	class GenericOut {

	protected:
		GenericOut(parType::detType det) : detector(det) {}

	public:
		virtual ~GenericOut() = default;

		parType::detType GetDetector() const { return detector; }

		// Output variables every particle type should have from any detector
		float M[3];
		float et;
		float Ekin;
		float theta;
		float phi;
		float cos_thetaH;

	private:
		parType::detType detector;

	};

	class GobbiOut : public GenericOut {

	public:
		GobbiOut() : GenericOut(parType::detType::Gobbi) {}

		// Output variables unique to Gobbi
		int iCsI;
		int itele;
		int ifront;
		int iback;
		float CsITime;
		float time;
		float btime;
		float timediff;
		float energyR;
		float denergyR;
		float energy;
		float denergy;
		float qdc;

	};

	class TexNeutOut : public GenericOut {

	public:
		TexNeutOut() : GenericOut(parType::detType::TexNeut) {}

	};

	class S800Out : public GenericOut {

	public:
		S800Out() : GenericOut(parType::detType::S800) {}

		// Output variables unique to S800 (two sets of angles, one
		// from SFA and one from S800, I think)
		float theta_s800;
		float phi_s800;

	};

	/********************************/

	wood(correl2& correl, std::string name, TDirectory* dir, bool gamma0);
	~wood();

	void initFrag(size_t i, parType::detType det);
	void loadFrag(size_t i, parType::detType det, solution* frag);

	// NOTE: This version of the code doesn't use CAESAR or any other gamma detector.
	// As such, I do not know what this function should look like, and it will have
	// to be completed in the future case where gammas are required.
	void setGammas();

	void Fill(float _Erel, float _Ex, float _Vcm, float _thetaCM, int _runnum, int _beamZ);

private:

	void loadFragCommon(size_t i, solution* frag);

	TTree* t;

	std::vector<GenericOut*> frags;

	// Separate vector for each detector type output, these are
	// what is actually output to the TTree
	std::deque<GobbiOut> gobbiFrags;
	std::deque<TexNeutOut> texNeutFrags;
	std::deque<S800Out> s800Frags;

	// Maps the original solution index to the new fragment
	// index in the above vectors. Each of the above vectors
	// should have a corresponding map below
	std::map<size_t, size_t> gobbiIndexMap;
	std::map<size_t, size_t> texNeutIndexMap;
	std::map<size_t, size_t> s800IndexMap;

	// Other member variables
	bool gamma;

	double Erel{NAN};
	double Ex{NAN};
	double Vcm{NAN};
	double thetaCM{NAN};

	int Ngamma;
	int Ngamma_Select;
	float Egamma[15];
	float Egamma_Select[15];
	float Tgamma[15];
	float Tgamma_Select[15];
	int Chgamma[15];
	int Chgamma_Select[15];

	int runnum;
	int beamZ;

};

#endif
