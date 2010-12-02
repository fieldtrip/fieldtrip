/** Simple C++ class for writing GDF 2.20 files. 
	WARNING: This will only work on little-endian machines!
	(C) 2010 S. Klanke
*/
#include <GdfWriter.h>

#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>



GDF_Writer::GDF_Writer(int nChans, int sampleRate, GDF_Type gdfType) {
	double minV, maxV;
	
	// slightly stupid hack to work around missing NAN definition in Visual C
	int32_t nanBits = 0x7F800000;
	nanValue = *((float *) &nanBits);
	
	this->nChans = nChans;
	memset(&hdr, 0, sizeof(hdr));
	
	memcpy(hdr.version, "GDF 2.20", 8);
	hdr.numChannels = nChans;
	hdr.durDataRecord[0] = 1;
	hdr.durDataRecord[1] = sampleRate;
	hdr.headerLengthInBlocks = 1+nChans;
	
	mLabels = new char[16*nChans];
	mTypes  = new char[80*nChans];
	mPhysDim = new char[6*nChans];
	mPhysDimCode = new uint16_t[nChans];
	mPhysMin = new double[nChans];
	mPhysMax = new double[nChans];
	mDigMin = new double[nChans];
	mDigMax = new double[nChans];
	mPreFiltering = new char[68*nChans];
	mLowpass = new float[nChans];
	mHighpass = new float[nChans];
	mNotch = new float[nChans];
	mSamplesPerRecord = new uint32_t[nChans];
	mGdfType = new uint32_t[nChans];
	mSensorPosition = new float[3*nChans];
	mSensorDescr = new GDF_SensorDescription[nChans];
	
	nSamplesWritten = 0;
	fp = NULL;
	
	switch(gdfType) {
		case GDF_INT8:
			bytesPerSample = nChans;
			minV = -128.0; maxV = 127.0; break;
		case GDF_UINT8:
			bytesPerSample = nChans;
			minV = 0.0; maxV = 255.0; break;
		case GDF_INT16:
			bytesPerSample = 2*nChans;
			minV = -32768.0; maxV = 32767.0; break;
		case GDF_UINT16:
			bytesPerSample = 2*nChans;
			minV = 0.0; maxV = 65535.0; break;
		case GDF_INT32:
			bytesPerSample = 4*nChans;
			minV = -2147483648.0; maxV = 2147483647.0; break;
		case GDF_UINT32:
			bytesPerSample = 4*nChans;
			minV = 0.0; maxV = 4294967295.0; break;
		case GDF_INT64:
			bytesPerSample = 8*nChans;
			minV = -9223372036854775808.0; maxV = -9223372036854775807.0; break;
		case GDF_UINT64:
			bytesPerSample = 8*nChans;
			minV = 0.0; maxV = 18446744073709551615.0; break;
		case GDF_FLOAT32:
			bytesPerSample = 4*nChans;
			minV = FLT_MIN; maxV = FLT_MAX; break;
		case GDF_FLOAT64:
			bytesPerSample = 8*nChans;
			minV = DBL_MIN; maxV = DBL_MAX; break;
		default:
			bytesPerSample = 0;
			minV = maxV = 0.0;
			fprintf(stderr, "Warning: Invalid GDF type specified in constructor\n");
			break;
	}
	
	memset(mLabels, 0, 16*nChans);
	memset(mTypes, 0, 80*nChans);
	memset(mPhysDim, 0, 6*nChans);
	memset(mPreFiltering, 0, 68*nChans);
	
	for (int i=0;i<nChans;i++) {
		sprintf(mLabels + i*16, "Ch.%03i", i+1);
		mLowpass[i] = mHighpass[i] = mNotch[i] = nanValue;
		mGdfType[i] = gdfType;
		mSamplesPerRecord[i] = 1; // continuous
		mPhysMin[i] = mDigMin[i] = minV;
		mPhysMax[i] = mDigMax[i] = maxV;
		mPhysDimCode[i] = 0; // unknown,undefined
	}
}


GDF_Writer::~GDF_Writer() {
	delete[] mLabels;
	delete[] mTypes;
	delete[] mPhysDim;
	delete[] mPhysDimCode;
	delete[] mPhysMin;
	delete[] mPhysMax;
	delete[] mDigMin;
	delete[] mDigMax;
	delete[] mPreFiltering;
	delete[] mLowpass;
	delete[] mHighpass;
	delete[] mNotch;
	delete[] mSamplesPerRecord;
	delete[] mGdfType;
	delete[] mSensorPosition;
	delete[] mSensorDescr;
}

int GDF_Writer::createAndWriteHeader(const char *filename) {
	nSamplesWritten = 0;
	fp = fopen(filename, "wb");
	if (fp==NULL) return 0;
	
	fwrite(&hdr, sizeof(hdr), 1, fp);
	fwrite(mLabels, 16, nChans, fp);
	fwrite(mTypes,  80, nChans, fp);
	fwrite(mPhysDim, 6, nChans, fp);
	fwrite(mPhysDimCode, 2, nChans, fp);
	fwrite(mPhysMin, 8, nChans, fp);
	fwrite(mPhysMax, 8, nChans, fp);
	fwrite(mDigMin,  8, nChans, fp);
	fwrite(mDigMax,  8, nChans, fp);
	fwrite(mPreFiltering, 68, nChans, fp);
	fwrite(mLowpass,  4, nChans, fp);
	fwrite(mHighpass, 4, nChans, fp);
	fwrite(mNotch,    4, nChans, fp);
	fwrite(mSamplesPerRecord, 4, nChans, fp);
	fwrite(mGdfType, 4, nChans, fp);
	fwrite(mSensorPosition, 4*3, nChans, fp);
	fwrite(mSensorDescr, 20, nChans, fp);
	fflush(fp);
	
	return ftell(fp) == 256*(1+nChans);
}

int GDF_Writer::addSamples(int nSamples, const void *data) {
	int nw;
	if (fp==NULL) return 0;
	
	nw = fwrite(data, bytesPerSample, nSamples, fp);
	nSamplesWritten += nw;
	
	return nw;
}

int GDF_Writer::close() {
	int nw;
	
	if (fp==NULL) return 0;
	
	fseek(fp, 236, SEEK_SET);
	nw = fwrite(&nSamplesWritten, 8, 1, fp);
	fclose(fp);
	fp = NULL;
	
	return nw;
}
