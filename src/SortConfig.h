/**
 * Created by Henry Webb (h.s.webb@wustl.edu) 6 March 2026.
 * This header file contains the SortConfig class, which loads and stores
 * all configuration information like filepaths and various settings.
 */

#ifndef SortConfig_H
#define SortConfig_H

#include <string>

class SortConfig {

private:
	std::string outputDir;
	std::string dataDir;
	std::string runNumbersFile;
	std::string itreeName;
	std::string hinpBranchName;
	std::string adcBranchName;
	std::string qdcBranchName;
	std::string tdcBranchName;
	std::string ofileName;
	std::string otreeName;
	std::string configDir;
	std::string lossDir;
	std::string PIDDir;
	std::string targetSuffix;
	std::string calDir;
	std::string frontEcalFile;
	std::string backEcalFile;
	std::string deltaEcalFile;
	std::string CsIEcalFile;
	std::string diamondEcalFile;
	std::string frontTimecalFile;
	std::string backTimecalFile;
	std::string deltaTimecalFile;
	std::string CsITimecalFile;
	std::string CsIStripExtentsFile;
	std::string CsIChannelMapFile;
	double targdist;      // in cm
	double targthick;     // in mg/cm^2
	double alThick;       // thickness of aluminum absorbers in front of Gobbi (mg/cm^2)
	size_t updateRate;
	size_t hinpboards;    // total number of HINP boards in experimental setup
	size_t hinpchans;     // number of channels per HINP board (this should always be 32)
	double gobbiHoleSize; // size of inner hole in Gobbi 28 (or new Gobbi) in mm

public:
	SortConfig(std::string configFilePath);

	// Getters
	std::string GetOutputDir() const { return outputDir; }
	std::string GetDataDir() const { return dataDir; }
	std::string GetRunNumbersFile() const { return runNumbersFile; }
	std::string GetItreeName() const { return itreeName; }
	std::string GetHinpBranchName() const { return hinpBranchName; }
	std::string GetAdcBranchName() const { return adcBranchName; }
	std::string GetQdcBranchName() const { return qdcBranchName; }
	std::string GetTdcBranchName() const { return tdcBranchName; }
	std::string GetOfileName() const { return ofileName; }
	std::string GetOtreeName() const { return otreeName; }
	std::string GetConfigDir() const { return configDir; }
	std::string GetLossDir() const { return lossDir; }
	std::string GetPIDDir() const { return PIDDir; }
	std::string GetTargetSuffix() const { return targetSuffix; }
	std::string GetCalDir() const { return calDir; }
	std::string GetFrontEcalFile() const { return frontEcalFile; }
	std::string GetBackEcalFile() const { return backEcalFile; }
	std::string GetDeltaEcalFile() const { return deltaEcalFile; }
	std::string GetCsIEcalFile() const { return CsIEcalFile; }
	std::string GetDiamondEcalFile() const { return diamondEcalFile; }
	std::string GetFrontTimecalFile() const { return frontTimecalFile; }
	std::string GetBackTimecalFile() const { return backTimecalFile; }
	std::string GetDeltaTimecalFile() const { return deltaTimecalFile; }
	std::string GetCsITimecalFile() const { return CsITimecalFile; }
	std::string GetCsIStripExtentsFile() const { return CsIStripExtentsFile; }
	std::string GetCsIChannelMapFile() const { return CsIChannelMapFile; }
	double GetTargDist() const { return targdist; }
	double GetTargThick() const { return targthick; }
	double GetAlThick() const { return alThick; }
	size_t GetUpdateRate() const { return updateRate; }
	size_t GetHinpboards() const { return hinpboards; }
	size_t GetHinpchans() const { return hinpchans; }
	double GetGobbiHoleSize() const { return gobbiHoleSize; }
};

#endif
