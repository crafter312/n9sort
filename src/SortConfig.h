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
	std::string ofileName;
	std::string otreeName;
	std::string lossDir;
	std::string PIDDir;
	std::string targetSuffix;
	std::string calDir;
	std::string frontEcalFile;
	std::string backEcalFile;
	std::string deltaEcalFile;
	std::string diamondEcalFile;
	std::string frontTimecalFile;
	std::string backTimecalFile;
	std::string deltaTimecalFile;
	float targdist;
	float targthick;
	size_t updateRate;

public:
	SortConfig(std::string configFilePath);

	// Getters
	std::string GetOutputDir() const { return outputDir; }
	std::string GetDataDir() const { return dataDir; }
	std::string GetRunNumbersFile() const { return runNumbersFile; }
	std::string GetItreeName() const { return itreeName; }
	std::string GetOfileName() const { return ofileName; }
	std::string GetOtreeName() const { return otreeName; }
	std::string GetLossDir() const { return lossDir; }
	std::string GetPIDDir() const { return PIDDir; }
	std::string GetTargetSuffix() const { return targetSuffix; }
	std::string GetCalDir() const { return calDir; }
	std::string GetFrontEcalFile() const { return frontEcalFile; }
	std::string GetBackEcalFile() const { return backEcalFile; }
	std::string GetDeltaEcalFile() const { return deltaEcalFile; }
	std::string GetDiamondEcalFile() const { return diamondEcalFile; }
	std::string GetFrontTimecalFile() const { return frontTimecalFile; }
	std::string GetBackTimecalFile() const { return backTimecalFile; }
	std::string GetDeltaTimecalFile() const { return deltaTimecalFile; }
	float GetTargDist() const { return targdist; }
	float GetTargThick() const { return targthick; }
	size_t GetUpdateRate() const { return updateRate; }
};

#endif
