/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */
#ifndef __PixelDataGrabber_h
#define __PixelDataGrabber_h

#include <vector>
#include <string>

#include <windows.h>

#include <siemensap.h>

class FolderWatcher;

/** Simple class for managing a contiguous piece of memory with convenient resizing. 
	Pretty close to what std::vector<char> provides, but this implementation better
	fits a truly typeless void *buffer.
*/
class SimpleBuffer {
	public:
	
	/** Constructor. Does not allocate anything, sets capacity and size to 0 */
	SimpleBuffer() {
		m_data = NULL;
		m_capacity = m_size = 0;
	}
	
	/** Destructor. Just calls clear() */
	~SimpleBuffer() {
		clear();
	}
	
	/** Resizes the buffer. If no memory has been allocated so far, exactly 'size' bytes will
		be allocated. Otherwise, if 'size' is bigger than the current capacity, memory will
		be re-allocated. If 'size' is smaller than or equal to capacity, the memory will
		be kept as it is.
		@param size		Desired buffer size
		@return true	On success (ie., the buffer now holds at least 'size' bytes)
				false 	On memory allocation errors (capacity and size are unchanged)
	*/
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
	
	/** Free any internally allocated memory and set capacity and size to 0 */
	void clear() {
		if (m_data!=NULL) free(m_data);
		m_data = NULL;
		m_size = m_capacity = 0;
	}	
	
	/** Returns pointer to internally allocated memory */
	void *data() { return m_data; }
	
	/** Returns (minimum guaranteed) size of the buffer */
	unsigned int size() { return m_size; }
	
	protected:
	
	unsigned int m_size;
	unsigned int m_capacity;
	void *m_data;
};

/** The main class for monitoring a directory on the scanner host, and for transmitting
	protocol information, pixel data, and timestamps to a FieldTrip buffer.
*/
class PixelDataGrabber {
	public: 
	
	/** Enumeration of action and error codes as returned by getLastAction() after a call to run()
	*/
	enum Action { Nothing, BadPixelData, OutOfMemory, PixelsTransmitted, ProtocolRead, TransmissionError };
	
	/** Construct a PixelDataGrabber. Just sets internal state variables */
	PixelDataGrabber();
	/** Destructor. Will automatically close open connections and stop monitoring */
	~PixelDataGrabber();
	
	/** Monitor the given directory 
		@param 	directory 	Pathname as a 0-terminated string
		@return true		If a directory is now being monitored	
				false		Otherwise
	*/
	bool monitorDirectory(const char *directory);
	
	/** Connect to FieldTrip buffer, or disconnect if hostname == NULL.
		If a connection is already open, it will be closed before trying
		to open the new connection.
		
		@param hostname  A 0-terminated string encoding the hostname
		@param port		 The port number (default = 1972)
		@return true  	if a new connection was made
				false 	in case of errors, but also when called with hostname==NULL for closing.
	*/
	bool connectToFieldTrip(const char *hostname, int port=1972);
	
	/** Returns the status of the FieldTrip connection
		@return  true	if a FieldTrip buffer is connected
				 false  otherwise
	*/
	bool isConnected() const { return (ftbSocket>0); }
	
	/** Returns the status of the FolderWatcher mechanism
		@return  true   if a directory is monitored
				 false  otherwise
	*/
	bool isListening() const { return fwActive; }
	
	/** Determines how much status messages should be printed to stderr */
	void setVerbosity(int howMuch) { verbosity = howMuch; }
	
	/** Runs the PixelDataGrabber by waiting for 'ms' milliseconds for changes in the directory,
		and in case something relevant happened, carrying out the corresponding actions. When
		this function returns 1, you can get more detailed information using getLastAction().
		
		@param ms 		Maximum time to wait for directory events
		@return  	-1 	in case of the PixelDataGrabber is not properly set up
					0	in case nothing happened during the wait
					1   in case an event was handled (but getLastAction() might still yield Nothing)
	*/
	int run(unsigned int ms);
	
	/** Returns the number of slices as determined from the latest protocol file read */
	int getNumSlices() const { return numSlices; }
	/** Returns the readout resolution (number of pixels in X direction) */
	int getReadoutResolution() const { return readResolution; }
	/** Returns the phase resolution (number of pixels in Y direction) */
	int getPhaseResolution() const { return phaseResolution; }
	/** Returns the readout field of view (slice sice in X direction, mm) */
	double getReadoutFOV() const { return readoutFOV; }
	/** Returns the phase field of view (slice sice in Z direction, mm) */
	double getPhaseFOV() const { return phaseFOV; }
	/** Returns the full pathname of the last file read, can be either
		protocol data (...\mrprot.txt) or pixel data (...\xxx.PixelData).
	*/
	const char *getLastFilename() { return fullName.c_str(); }
	/** Returns the number of scans transmitted to the currently connected
		FieldTrip buffer.
	*/
	int getNumScansWritten() { return samplesWritten; }
	/** Returns a code for errors or the action carried out during the last run() */
	Action getLastAction() const { return lastAction; }

	protected:	
	
	/** This function writes the protocol information into the header of the buffer, 
		which also clears any previous data in there. Only call this function
		if there is indeed protocol information (protBuffer.size() != 0) and
		a FieldTrip buffer is connected (ftbSocket != -1). The function will set
		lastAction to TransmissionError in case of errors, and headerWritten to true
		in case of success.
		@return true 	on success
		@return false 	in case of errors
	*/
	bool writeHeader();
	
	/** This function writes the latest pixel data into the FieldTrip buffer by
		appending one sample. Only call this if there really is pixel data
		(sliceBuffer.size()!=0), and if a FieldTrip buffer is connected (ftbSocket != -1).
		The function will always set lastAction to either OutOfMemory, TransmissionError, or 
		PixelsTransmitted, depending on the success of the operations.
		@return true 	on success
		@return false 	in case of errors
	*/
	bool writePixelData();
	
	/** This function writes the latest timestamp as an event to the FieldTrip buffer.
		Only call this if a FieldTrip buffer is connected (ftbSocket != -1).
		The function will set lastAction to TransmissionError in case of errors.
		@param  tv      timeval struct, containing seconds + microseconds since 1970 (unix time)
		@return true 	on success
		@return false 	in case of errors
	*/
	bool writeTimestampEvent(const struct timeval &tv);
	
	/** This function is used to transmit a new scan to the FieldTrip buffer. Only call this
		if a buffer is connected (ftbSocket != -1). If no header has been written since we 
		connected to the FieldTrip buffer, writeHeader() is called.
		
		With both protocol and headers present, writePixelData() and writeTimestampEvent()
		are called. This function returns false as soon as one of the aforementioned subfunctions
		return false.
		@param  tv      timeval struct, containing seconds + microseconds since 1970 (unix time)
		@return true 	on success
		@return false 	in case of errors
	*/
	bool sendFrameToBuffer(const struct timeval &tv);
	
	/** Handle new protocol information as detected from monitoring the directory, and previously
		read into memory using tryReadFile(). This function will parse the ASCII representation
		into the linked list pointed to by protInfo, and try to retrieve values for
		readResolution, phaseResolution, numSlices, phaseFOV and readoutFOV;
		On errors, the xxxResolution and numSlices are set to 0.
		@param  info   		ASCII representation, not necessarily 0-terminated
		@param  sizeInBytes	length of info
	*/
	void handleProtocol(const char *info, unsigned int sizeInBytes);
	
	/** Try to read the complete contents of a file into a given SimpleBuffer.
		This is used internally for pixel data (->pixBuffer) and protocol information
		(->protBuffer). If the file cannot be opened, the sBuf is not modified,
		otherwise, sBuf is resized to the number of bytes that could actually be read.
		@param filename		Name of the file to be read
		@param sBuf			SimpleBuffer to receive the contents
		@return	true		On success (the file could be read completely)
				false		In case of errors
	*/
	bool tryReadFile(const char *filename, SimpleBuffer &sBuf);
	
	/** Try to read the protocol from the default location,which is <watch_directory>/mrprot.txt.
		@return true 	on success (irrespective of the protocol contents)
				false	in case the protocol file is not present
	*/
	bool tryReadProtocol();
	
	/** This function handles the results of the directory monitoring and will
		eventually call sendFrameToBuffer() in case new pixel data could be read, or 
		writeHeader() in case new protocol data has been detected. If pixel data
		comes in without having read protocol information, this function will call 
		tryReadProtocol() before reshapeToSlices().
	*/
	void tryFolderToBuffer();
	
	/** This function reshapes the contents of pixBuffer (where new pixel data is read into,
		and where slices are tiles of a 2D mosaic) to sliceBuffer (where slices are contiguous 
		in memory). If errors occur, the sliceBuffer will be empty.		
	*/
	void reshapeToSlices();

	std::string sourceDir;	/**< Contains the path of the directory that is being monitored */
	std::string fullName;	/**< Contains the full path of the latest read file */
	FolderWatcher *FW;		/**< Points to a FolderWatcher object */
	HANDLE fwEventHandle;	/**< WIN32 event handle used for the FolderWatcher */
	int ftbSocket;			/**< Socket identifying the remote FieldTrip buffer, or -1 if unconnected */
	bool fwActive;			/**< Flag that determines whether the FolderWatcher is indeed monitoring */
	int verbosity;			/**< Determines how much information should be printed to the console */ 
	bool headerWritten;		/**< Flag that determines whether header information has been written to the currently connected buffer */
	
	unsigned int samplesWritten;	/**< Number of samples written to the currently connected buffer */
	sap_item_t *protInfo;			/**< Linked list of protocol information */
	unsigned int readResolution;	/**< Number of pixels in readout direction (X) */
	unsigned int phaseResolution;	/**< Number of pixels in phase direction (Y) */
	unsigned int numSlices;			/**< Number of slices */
	double phaseFOV, readoutFOV;	/**< Size of the field of view in mm */
	Action lastAction;				/**< Contains last action or error that occured */
	
	SimpleBuffer pixBuffer;		/**< Simple buffer that contains pixel data as read from file (e.g., mosaic) */
	SimpleBuffer sliceBuffer;	/**< Simple buffer that contains slice-shaped pixel data */
	SimpleBuffer protBuffer;	/**< Simple buffer that contains ASCII protocol information */
};

#endif
