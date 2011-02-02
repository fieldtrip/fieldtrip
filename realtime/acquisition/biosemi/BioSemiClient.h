#ifndef __BioSemiClient_h
#define __BioSemiClient_h

#ifndef WIN32
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#endif

#include <Labview_DLL.h>

#define BUFFER_SIZE   32*1024*1024  // 32 MByte as suggested by BioSemi
#define BUFFER_LEN    8*1024*1024   // in int32 units
#define SYNC_BV       (int32_t) 0xFFFFFF00

struct BioSemiBlock {
	int startIndex;
	int numSamples;
	int numInSync;
	int stride;	
	bool batteryLow;
};

class BioSemiClient {
	public:
	
	BioSemiClient();
	~BioSemiClient();
	
	bool openDevice();
	void closeDevice();
	
	int  getNumChannels()  const { return numChannels; };
	int  getNumChanAIB()   const { return numChanAIB; };
	int  getSpeedMode()    const { return speedMode; };
	int  getSamplingFreq() const { return sampleFreq; };
	int  getCurPointer()  ;
	
	bool checkNewBlock(BioSemiBlock &block);
	
	void msleep(int millis) const {
		#ifdef WIN32
		Sleep(millis);
		#else
		usleep(millis*1000);
		#endif
	}
	
	
	int getValue(int index) const {
	    return ringBuffer[index & (BUFFER_LEN-1)]; //  ? index - BUFFER_LEN : index];
	}	
	
	double getCurrentTime() {
		#ifdef WIN32
		return timeGetTime() * 0.001;
		#else
		struct timeval tv;
		gettimeofday(&tv, NULL);
		return tv.tv_sec + tv.tv_usec*1e-6;
		#endif
	}
	
	const int32_t *getRingbuffer() {
		return ringBuffer;
	}

	bool isDriverOk() const { return driverOk; }
	
	protected:
	

	
	bool driverOk, deviceOpen;
	double timeHandshake;
	int numChannels, numChanAIB;
	int speedMode;
	int sampleFreq;
	bool deviceIsMk2;
	bool cmsInRange;
	int stride;
	int lastPointerRead;
	int lastStartPtr;
	
	OPEN_DRIVER_ASYNC_T    lv_open_driver_async;
	USB_WRITE_T            lv_usb_write;
	READ_MULTIPLE_SWEEPS_T lv_read_multiple_sweeps;
	READ_POINTER_T         lv_read_pointer;
	CLOSE_DRIVER_ASYNC_T   lv_close_driver_async;
	HANDLE deviceHandle;
	int32_t *ringBuffer;

#ifdef WIN32
	HINSTANCE hLib; // handle to Labview_DLL.dll
#else
	void *hLib;
#endif
};

#endif
