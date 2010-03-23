#ifndef __PixelDataGrabber_h
#define __PixelDataGrabber_h

#include <vector>
#include <string>

#include <windows.h>

#include <siemensap.h>

class FolderWatcher;

class SimpleBuffer {
	public:
	
	SimpleBuffer() {
		m_data = NULL;
		m_capacity = m_size = 0;
	}
	
	~SimpleBuffer() {
		clear();
	}
	
	bool resize(unsigned int size) {
		if (size <= m_capacity) {
			m_size = size;
			return true;
		}
		if (m_data != NULL) {
			void *nd = realloc(m_data, size);
			if (nd == NULL) return false;
			m_data = nd;
		} else {
			m_data = malloc(size);
			if (m_data == NULL) return false;
		}
		m_size = m_capacity = size;
		return true;
	}
	
	void clear() {
		if (m_data!=NULL) free(m_data);
		m_data = NULL;
		m_size = m_capacity = 0;
	}	
	
	void *data() { return m_data; }
	unsigned int size() { return m_size; }
	
	protected:
	
	unsigned int m_size;
	unsigned int m_capacity;
	void *m_data;
};

class PixelDataGrabber {
	public: 
	
	enum Action { Nothing, BadPixelData, OutOfMemory, PixelsTransmitted, ProtocolRead };
	
	PixelDataGrabber();
	~PixelDataGrabber();
	
	bool monitorDirectory(const char *directory);
	bool connectToFieldTrip(const char *hostname, int port);
	
	/** Runs the PixelDataGrabber by waiting of 'ms' milliseconds for changes in the directory,
		and in case something relevant happened, carrying out the corresponding actions.
		@param ms 		Maximum time to wait for directory events
		@return  	-1 	in case of the PixelDataGrabber is not properly set up
					0	in case nothing happened during the wait
					1   in case an event was handled (but getLastAction() might still yield Nothing)
	*/
	int run(unsigned int ms);
	
	int getNumSlices() const { return numSlices; }
	int getReadoutResolution() const { return readResolution; }
	int getPhaseResolution() const { return phaseResolution; }
	double getReadoutFOV() const { return readoutFOV; }
	double getPhaseFOV() const { return phaseFOV; }
	const char *getLastFilename() { return fullName.c_str(); }
	Action getLastAction() const { return lastAction; }

	protected:	
	
	void writeHeader();
	bool writePixelData();
	void writeTimestampEvent(const struct timeval &tv);
	void sendFrameToBuffer(const struct timeval &tv);
	void handleProtocol(const char *info, unsigned int sizeInBytes);
	bool tryReadFile(const char *filename, SimpleBuffer &sBuf);
	void tryReadProtocol();
	void tryFolderToBuffer();
	void reshapeToSlices();

	std::string sourceDir;
	std::string fullName;
	FolderWatcher *FW;
	HANDLE fwEventHandle;
	int ftbSocket;
	
	unsigned int samplesWritten;
	sap_item_t *protInfo;
	unsigned int readResolution, phaseResolution, numSlices;
	double phaseFOV, readoutFOV;
	Action lastAction;
	
	SimpleBuffer pixBuffer, sliceBuffer, protBuffer;
};

#endif
