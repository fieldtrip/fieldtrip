/** Simple C++ class for writing GDF 2.20 files.
	WARNING: This will only work on little-endian machines!
	
	(C) 2010 S. Klanke
*/

#ifndef __GdfWriter_h
#define __GdfWriter_h

#include <stdio.h>
#include <stdint.h>
#include <string.h>

enum GDF_Type {
	GDF_INT8 = 1,
	GDF_UINT8 = 2,
	GDF_INT16 = 3,
	GDF_UINT16 = 4,
	GDF_INT32 = 5,
	GDF_UINT32 = 6,
	GDF_INT64 = 7,
	GDF_UINT64 = 8,
	GDF_FLOAT32 = 16,
	GDF_FLOAT64 = 17
	// GDF_FLOAT128 = 18,
	// GDF_INT24 = 279,
	// GDF_UINT24 = 535
};

enum GDF_Dimension {
	GDF_UNKNOWN = 0,
	GDF_DIMLESS = 512,
	GDF_PERCENT = 544,
	GDF_DEGREE = 736,
	GDF_RADIAN = 768,
	GDF_HZ = 2496,
	GDF_MMHG = 3872,
	GDF_VOLT = 4256,
	GDF_OHM = 4288,
	GDF_KELVIN = 4384,
	GDF_CELSIUS = 6048,
	GDF_L_MIN = 3072,
	GDF_L_MIN_M2 = 2848,
	GDF_DYNS_CM5 = 4128,
	GDF_DYNS_M2_CM5 = 6016
};

enum GDF_Scale {
	GDF_YOTTA = 10, // 1e24
	GDF_ZETTA = 9,  // 1e21
	GDF_EXA = 8,    // 1e18
	GDF_PETA = 7,   // 1e15
	GDF_TERA = 6,   // 1e12
	GDF_GIGA = 5,   // 1e9
	GDF_MEGA = 4,   // 1e6
	GDF_KILO = 3,   // 1000
	GDF_HECTO = 2,  // 100
	GDF_DECA = 1,   // 10
	GDF_DECI = 16,  // 0.1
	GDF_CENTI = 17, // 0.01
	GDF_MILLI = 18, // 0.001
	GDF_MICRO = 19, // 1e-6
	GDF_NANO = 20,  // 1e-9
	GDF_PICO = 21,  // 1e-12
	GDF_FEMTO = 22, // 1e-15
	GDF_ATTO = 23,  // 1e-18
	GDF_ZEPTO = 24, // 1e-21
	GDF_YOCTO = 25  // 1e-24
};

#pragma pack(push,1)

union GDF_SensorDescription {
	float electrodeImpedance;
	float probeFrequency;
	char reserved[20];
};

struct GDF_Header {
	char version[8];
	char patient[66];
	uint8_t reserved[10];
	uint8_t smokingAlcoholDrugMedication;
	uint8_t weight;
	uint8_t height;
	uint8_t genderHandedVisualHeart;
	char recordingId[64];
	uint32_t recordingLocation[4];
	uint32_t recordingTime[2];
	uint32_t birthday[2];
	uint16_t headerLengthInBlocks;
	uint8_t patientClass[6];
	uint64_t equipmentId;
	uint8_t reserved2[6];
	uint16_t headsize[3];
	float posRefElec[3];
	float posGroundElec[3];
	int64_t numDataRecords;
	uint32_t durDataRecord[2];
	uint16_t numChannels;
	uint16_t reserved3;

	// compiling this function will always yield an error, which is what we want!
	template<typename T> static GDF_Type getType(T dummy) {
		return (GDF_Type)0;
	}

	// "specialised versions" of the above template function, returning the
	// right GDF_Type for the given input argument
	static GDF_Type getType(int8_t dummy)   { return GDF_INT8; }
	static GDF_Type getType(uint8_t dummy)  { return GDF_UINT8; }
	static GDF_Type getType(int16_t dummy)  { return GDF_INT16; }
	static GDF_Type getType(uint16_t dummy) { return GDF_UINT16; }
	static GDF_Type getType(int32_t dummy)  { return GDF_INT32; }
	static GDF_Type getType(uint32_t dummy) { return GDF_UINT32; }
	static GDF_Type getType(int64_t dummy)  { return GDF_INT64; }
	static GDF_Type getType(uint64_t dummy) { return GDF_UINT64; }
	static GDF_Type getType(float dummy)    { return GDF_FLOAT32; }
	static GDF_Type getType(double dummy)   { return GDF_FLOAT64; }

};

#pragma pack(pop)

class GDF_Writer {
	public:

	/** Create a GDF_Write object by specifying number of channels, samplerate, and data type.
		Data will always be written in a continuous fashion, that is, record length = 1.
		Reasonable defaults will be set up for all the Phys/Dig-Min/Max values etc.
		You can change these if you like, but do so *before* calling createAndWriteHeader.
	*/
	GDF_Writer(int nChans, int sampleRate, GDF_Type gdfType);
	~GDF_Writer();

	/** Creates the specified file and writes the header */
	int createAndWriteHeader(const char *filename);

	/** Add samples (pointed to by 'data') to the file. Example: 4 channels have been specified,
		the data type is int32, and you want to write 3 samples, then data needs to point to an array of
		12 integers, in this channel/sample order: 1/1 2/1 3/1 4/1 1/2 2/2 3/2 4/2 1/3 2/3 3/3 4/3
		and you would use nSamples = 3.
		Note that this memory layout corresponds to a record length of 1 sample in the GDF way of thinking.
	*/
	int addSamples(int nSamples, const void *data);

	/** Update the sample counter (=number of records) and close the file */
	int close();

	void setPhysicalLimits(int channel, double minV, double maxV) {
		mPhysMin[channel] = minV;
		mPhysMax[channel] = maxV;
	}

	void setDigitalLimits(int channel, double minV, double maxV) {
		mDigMin[channel] = minV;
		mDigMax[channel] = maxV;
	}

	/** Call this to set filter information for the specific channel */
	void setFilters(int channel, float lowpass, float highpass, float notch) {
		mLowpass[channel]=lowpass;
		mHighpass[channel]=highpass;
		mNotch[channel]=notch;
	}

	/** Call this to set the dimension code of a particular channel (0..N-1).
		Use a value from GDF_Dimension and optionally and one from GDF_Scale,
		Example:  setPhysDimCode(0, GDF_MILLI + GDF_VOLT)
		to set the first channel to mV units.
	*/
	void setPhysDimCode(int channel, int code) {
		mPhysDimCode[channel] = code;
	}

	/** Call this to set the label of the specified channel (0..N-1).
		Only the first 16 characters will be considered.
	*/
	void setLabel(int channel, const char *label) {
		if (channel<0 || channel>=nChans) return;
		strncpy(mLabels + 16*channel, label, 16);
	}

	static int getSizeAndRangeByType(GDF_Type gdfType, double& minV, double& maxV);

	GDF_Header hdr;

	protected:

	union {
		float asFloat;
		int32_t asInt;
	} nanValue;
	int nChans, bytesPerSample;
	FILE *fp;

	char *mLabels; // 16 bytes each
	char *mTypes;  // 80 bytes each
	char *mPhysDim; // 6 bytes each
	uint16_t *mPhysDimCode;
	double *mPhysMin;
	double *mPhysMax;
	double *mDigMin;
	double *mDigMax;
	char *mPreFiltering; // 68 bytes each
	float *mLowpass;
	float *mHighpass;
	float *mNotch;
	uint32_t *mSamplesPerRecord;
	uint32_t *mGdfType;
	float *mSensorPosition; // 3 floats each
	GDF_SensorDescription *mSensorDescr;

	int64_t nSamplesWritten;
};

#endif
