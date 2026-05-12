/**
 * Created by Henry Webb (h.s.webb@wustl.edu) 6 March 2026.
 * This implementation file contains the SortConfig class, which loads and
 * stores all configuration information like filepaths and various settings.
 */

#include "SortConfig.h"

#include <exception>
#include <fstream>
#include <iostream>

#include "constants.h"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SortConfig::SortConfig(string configFilePath) {
	cout << "Reading sort code config file..." << endl;

	// Open config file, check that it exists	
	ifstream configfile(configFilePath);
	if (configfile.fail()) throw invalid_argument(string(BOLDRED) + string("Config file ") + configFilePath + string(" does not exist or failed to open") + string(RESET));

	// Read config file
	string line;
	while (getline(configfile, line)) {
		if (line.find("outputDir") != string::npos)
			outputDir = line.substr(line.find('=') + 2);
		else if (line.find("dataDir") != string::npos)
			dataDir = line.substr(line.find('=') + 2);
		else if (line.find("runNumbersFile") != string::npos)
			runNumbersFile = line.substr(line.find('=') + 2);
		else if (line.find("itreeName") != string::npos)
			itreeName = line.substr(line.find('=') + 2);
		else if (line.find("adcBranchName") != string::npos)
			adcBranchName = line.substr(line.find('=') + 2);
		else if (line.find("qdcBranchName") != string::npos)
			qdcBranchName = line.substr(line.find('=') + 2);
		else if (line.find("ofileName") != string::npos)
			ofileName = line.substr(line.find('=') + 2);
		else if (line.find("otreeName") != string::npos)
			otreeName = line.substr(line.find('=') + 2);
		else if (line.find("configDir") != string::npos)
			configDir = line.substr(line.find('=') + 2);
		else if (line.find("lossDir") != string::npos)
			lossDir = line.substr(line.find('=') + 2);
		else if (line.find("PIDDir") != string::npos)
			PIDDir = line.substr(line.find('=') + 2);
		else if (line.find("targetSuffix") != string::npos)
			targetSuffix = line.substr(line.find('=') + 2);
		else if (line.find("calDir") != string::npos)
			calDir = line.substr(line.find('=') + 2);
		else if (line.find("frontEcalFile") != string::npos)
			frontEcalFile = line.substr(line.find('=') + 2);
		else if (line.find("backEcalFile") != string::npos)
			backEcalFile = line.substr(line.find('=') + 2);
		else if (line.find("deltaEcalFile") != string::npos)
			deltaEcalFile = line.substr(line.find('=') + 2);
		else if (line.find("CsIEcalFile") != string::npos)
			CsIEcalFile = line.substr(line.find('=') + 2);
		else if (line.find("diamondEcalFile") != string::npos)
			diamondEcalFile = line.substr(line.find('=') + 2);
		else if (line.find("frontTimecalFile") != string::npos)
			frontTimecalFile = line.substr(line.find('=') + 2);
		else if (line.find("backTimecalFile") != string::npos)
			backTimecalFile = line.substr(line.find('=') + 2);
		else if (line.find("deltaTimecalFile") != string::npos)
			deltaTimecalFile = line.substr(line.find('=') + 2);
		else if (line.find("CsITimecalFile") != string::npos)
			CsITimecalFile = line.substr(line.find('=') + 2);
		else if (line.find("CsIStripExtentsFile") != string::npos)
			CsIStripExtentsFile = line.substr(line.find('=') + 2);
		else if (line.find("CsIChannelMapFile") != string::npos)
			CsIChannelMapFile = line.substr(line.find('=') + 2);
		else if (line.find("targdist") != string::npos) {
			string temps = line.substr(line.find('=') + 2);
			try {
				targdist = stod(temps);
			}
			catch (...) {
				throw invalid_argument("targdist in config file " + configFilePath + " is not a valid double");
			}
		}
		else if (line.find("targthick") != string::npos) {
			string temps = line.substr(line.find('=') + 2);
			try {
				targthick = std::stod(temps);
			}
			catch (...) {
				throw invalid_argument("targthick in config file " + configFilePath + " is not a valid double");
			}
		}
		else if (line.find("updateRate") != string::npos) {
			string temps = line.substr(line.find('=') + 2);
			try {
				sscanf(temps.c_str(), "%zu", &updateRate); // Note that the `z` specifier is Linux only, and will have to be changed for this to work on Windows
			}
			catch (...) {
				throw invalid_argument("updateRate in config file " + configFilePath + " is not a valid size_t (unsigned integer)");
			}
		}
		else if (line.find("hinpboards") != string::npos) {
			string temps = line.substr(line.find('=') + 2);
			try {
				sscanf(temps.c_str(), "%zu", &hinpboards); // Note that the `z` specifier is Linux only, and will have to be changed for this to work on Windows
			}
			catch (...) {
				throw invalid_argument("hinpboards in config file " + configFilePath + " is not a valid size_t (unsigned integer)");
			}
		}
		else if (line.find("hinpchans") != string::npos) {
			string temps = line.substr(line.find('=') + 2);
			try {
				sscanf(temps.c_str(), "%zu", &hinpchans); // Note that the `z` specifier is Linux only, and will have to be changed for this to work on Windows
			}
			catch (...) {
				throw invalid_argument("hinpchans in config file " + configFilePath + " is not a valid size_t (unsigned integer)");
			}
		}
		else if (line.find("gobbiHoleSize") != string::npos) {
			string temps = line.substr(line.find('=') + 2);
			try {
				gobbiHoleSize = std::stod(temps);
			}
			catch (...) {
				throw invalid_argument("gobbiHoleSize in config file " + configFilePath + " is not a valid double");
			}
		}
	}
	configfile.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



