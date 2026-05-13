#ifndef Input_H
#define Input_H

/* Written by Henry Webb (h.s.webb@wustl.edu), August 2025
 * This class handles reading information from a SpecTcl-
 * generated ROOT file and reformatting it into a more
 * easily usable form. The resulting variables are stored
 * in this class on a per-event basis and can be accessed
 * via getter functions.
 * 
 * Modified by Henry Webb (h.s.webb@wustl.edu) and Johnathan
 * Phillips (j.s.phillips@wustl.edu) March 2026 for experiment
 * at TAMU Cyclotron Institute
 */

#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <string>
#include <optional>
#include <vector>

#include "SortConfig.h"

// Ensure that these match the configuration of your hardware system!
#define HINP_BOARD_COUNT 8
#define HINP_CHAN_COUNT 32
#define PSD_CHIP_COUNT 24
#define PSD_CHAN_COUNT 8
#define QDCV965_CHAN_COUNT 2 // CAEN v965, 16 channels with h and l each (32 total parameters)
#define ADC_CHAN_COUNT 32
#define QDC_CHAN_COUNT 32
#define TDC_CHAN_COUNT 48 // CAEN v1190a, 128 channels can be configured to accept up to 16 hits each
#define TDC_HIT_COUNT 3 // CAEN v1190a, 128 channels can be configured to accept up to 16 hits each

// Number of columns per parameter type in the input file
#define HINP_NCOLUMNS HINP_BOARD_COUNT * HINP_CHAN_COUNT
#define PSD_NCOLUMNS PSD_CHIP_COUNT * PSD_CHAN_COUNT
#define TDC_NCOLUMNS TDC_CHAN_COUNT * TDC_HIT_COUNT

// TODO: add private "checked access" helper function to all structs for boundary checking in getter functions

class Input {

public:
	Input(TTreeReader&, SortConfig&);
	~Input();

	void ReadAndRefactor();

	int badevt = 0;

	/******** GOBBI INPUT ********/

	struct GobbiInput {
		// Vectors for input branch readers
		std::vector<TTreeReaderValue<double>> eRVs;
		std::vector<TTreeReaderValue<double>> eLoRVs;
		std::vector<TTreeReaderValue<double>> tRVs;

		// Vectors for Gobbi (HINP) hit values
		size_t Nhits{0};
		std::vector<size_t> board;
		std::vector<size_t> chan;
		std::vector<size_t> e;
		std::vector<size_t> eLo;
		std::vector<size_t> t;

		void clear() {
			Nhits = 0;
			board.clear();
			chan.clear();
			e.clear();
			eLo.clear();
			t.clear();
		}

		// Hit getter functions
		size_t GetNhits() const         { return Nhits; }
		size_t GetBoard(size_t i) const { return board[i]; }
		size_t GetChan(size_t i) const  { return chan[i]; }
		size_t GetE(size_t i) const     { return e[i]; }
		size_t GetELo(size_t i) const   { return eLo[i]; }
		size_t GetT(size_t i) const     { return t[i]; }
	};

	struct PSDInput {
		// Vectors for input branch readers
		std::vector<TTreeReaderValue<double>> aRVs;
		std::vector<TTreeReaderValue<double>> bRVs;
		std::vector<TTreeReaderValue<double>> cRVs;
		std::vector<TTreeReaderValue<double>> tRVs;

		// Vectors for TexNeut (PSD) hit values
		size_t Nhits{0};
		std::vector<size_t> chip;
		std::vector<size_t> chan;
		std::vector<size_t> a;
		std::vector<size_t> b;
		std::vector<size_t> c;
		std::vector<size_t> t;

		void clear() {
			Nhits = 0;
			chip.clear();
			chan.clear();
			a.clear();
			b.clear();
			c.clear();
			t.clear();
		}

		// Hit getter functions
		size_t GetNhits() const        { return Nhits; }
		size_t GetChip(size_t i) const { return chip[i]; }
		size_t GetChan(size_t i) const { return chan[i]; }
		size_t GetA(size_t i) const    { return a[i]; }
		size_t GetB(size_t i) const    { return b[i]; }
		size_t GetC(size_t i) const    { return c[i]; }
		size_t GetT(size_t i) const    { return t[i]; }
	};
	
	// Input class supporting variety of CAEN ADC and QDCs which all share a common data format
	struct ADCQDCInput {
		// Vectors for input branch readers
		std::vector<TTreeReaderValue<double>> aqRVs;
		
		// Vectors for hit values
		size_t Nhits{0};
		std::vector<size_t> chan;
		std::vector<size_t> aq;
		
		void clear() {
			Nhits = 0;
			chan.clear();
			aq.clear();
		}
		
		// Hit getter functions
		size_t GetNhits() const { return Nhits; }
		size_t GetChan(size_t i) const { return chan[i]; }
		size_t GetAQ(size_t i) const { return aq[i]; }
	};
	
	// Input class specifically for CAEN v965 QDC
	struct QDCInputV965 {
		// Vectors for input branch readers
		std::vector<TTreeReaderValue<double>> qhRVs; // high range Q
		std::vector<TTreeReaderValue<double>> qlRVs; // low range Q
		
		// Vectors for QDC hit values
		size_t Nhits{0};
		std::vector<size_t> chan;
		std::vector<size_t> qh;
		std::vector<size_t> ql;
		
		void clear() {
			Nhits = 0;
			chan.clear();
			qh.clear();
			ql.clear();
		}
		
		// Hit getter functions
		size_t GetNhits() const { return Nhits; }
		size_t GetChan(size_t i) const { return chan[i]; }
		size_t GetQH(size_t i) const { return qh[i]; }
		size_t GetQL(size_t i) const { return ql[i]; }
	};
	
	struct TDCInput {
		// Vector for input branch readers
		std::vector<TTreeReaderValue<double>> tRVs; // flattened 2D array for channels and hits
		
		// Vector for TDC hit values
		size_t Nhits[TDC_CHAN_COUNT]; // one for each channel
		std::vector<double> t[TDC_CHAN_COUNT]; // outermost array stores channels, innermost vector stores hits
		
		void clear() {
			for (int i = 0;i < TDC_CHAN_COUNT; i++) {
				Nhits[i] = 0;
				t[i].clear();
			}
		}
		
		// Hit getter functions
		size_t GetNhits(size_t ch) const { return Nhits[ch]; }
		std::optional<double> GetT(size_t ch, size_t i) const { return (i >= t[ch].size()) ? std::nullopt : std::optional<double>(t[ch][i]); }
		
		// Group channel hits into single vector for interface with TexNeut
		// Make sure that this properly accesses the TDC channels dedicated to TexNeut (4-16 in my case)
		// For now, only consider the first hit in each TDC channel
		void FillTexNeutHitVectors(std::vector<size_t>& chan, std::vector<double>& tdct) const {
			for (int ch = 0; ch < TDC_CHAN_COUNT; ch++) {
				if (t[ch].size() == 0) continue;
				chan.push_back(ch);
				tdct.push_back(t[ch][0]);
			}
		}
	};

	// Getter functions (modify for expected input format)
	const GobbiInput& GetGobbi() const { return gobbi; }
	const PSDInput& GetTexNeut() const { return texneut; }
	const QDCInputV965& GetQDCV965() const { return qdcv965; }
	const ADCQDCInput& GetADC() const { return adc; }
	const ADCQDCInput& GetQDC() const { return qdc; }
	const TDCInput& GetTDC() const { return tdc; }

private:
	// TTreeReader reference for input
	TTreeReader& reader;

	// Input classes (modify for expected input)
	GobbiInput gobbi;
	PSDInput texneut;
	QDCInputV965 qdcv965;
	ADCQDCInput adc;
	ADCQDCInput qdc;
	TDCInput tdc;

	/******** PRIVATE STATIC HELPER FUNCTIONS ********/

	static std::vector<std::string> GenerateColumnNamesHINP(const std::string&, const std::string&);
	static std::vector<std::string> GenerateColumnNamesPSD(const std::string&);
	static std::vector<std::string> GenerateColumnNamesADCQDC(const std::string&, const size_t&);
	static std::vector<std::string> GenerateColumnNamesQDCV965(const std::string&);
	static std::vector<std::string> GenerateColumnNamesTDC(const std::string&);

};

#endif
