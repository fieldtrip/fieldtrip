/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */

#ifndef __Gdf_BackGroundWriter_h
#define __Gdf_BackGroundWriter_h

#include <GdfWriter.h>
#include <stdio.h>
#include <pthread.h>
#include <LocalPipe.h>
#include <string>

template <typename To>
class GDF_BackgroundWriter {
	public:

	GDF_BackgroundWriter(int nChans,  int sampleRate, GDF_Type gdfType, int seconds=5) : gdfWriter(nChans, sampleRate, gdfType) {
		running = threadStarted = false;
		this->nChans = nChans;
		rbSize = seconds*sampleRate; // 5 seconds of data for saving ring buffer
		rbData = new To[rbSize*nChans];
		rbWritePos = rbReadPos = 0;
		maxFileSize = 1024*1024*1024;
		fileCounter = 0;
	}

	GDF_Writer& gdf() {
		return gdfWriter;
	}

	void setMaxFileSize(int64_t maxSize) {
		maxFileSize = maxSize;
	}

	virtual ~GDF_BackgroundWriter() {
		stopSync();

		// clean up variables
		delete[] rbData;
	}

	bool start(const char *name) {
		lenBasename = strlen(name);

		// strip .gdf suffix, if any
		if (lenBasename > 4) {
			if ((toupper(name[lenBasename-1]) == 'F') &&
				(toupper(name[lenBasename-2]) == 'D') &&
				(toupper(name[lenBasename-3]) == 'G') &&
				(name[lenBasename-4] == '.')) {

				lenBasename-=4;
			}
		}

		filename.clear();
		filename.reserve(lenBasename + 8); // some extra space for numbers + .gdf
		filename.append(name, lenBasename);

		fileCounter = 0;

		if (pthread_create(&savingThread, NULL, staticSavingThreadFunction, this)) {
			fprintf(stderr, "Could not spawn GDF saving thread.\n");
			return false;
		}
		threadStarted = true;
		return true;
	}

	void stopSync() {
		if (running) {
			// stop saving thread
			int64_t minusOne = -1;
			locPipe.write(sizeof(int64_t), &minusOne);
		}
		pthread_join(savingThread, 0);
		threadStarted = false;
	}

	void stopAsync() {
		if (running) {
			// stop saving thread asynch
			pthread_detach(savingThread);
			int64_t minusTwo = -2;
			locPipe.write(sizeof(int64_t), &minusTwo);
		}
	}

	bool checkFreeBlock(int nSamples) {
		return (rbWritePos - rbReadPos <= rbSize - nSamples);
	}

	To *getSampleSlot() {
		int wrappedPos = (rbWritePos++) % rbSize;
		return rbData + wrappedPos * nChans;
	}

	void commitBlock() {
		int n = locPipe.write(sizeof(int64_t), &rbWritePos);
		if (n!=sizeof(int64_t)) {
			fprintf(stderr, "GDF_BackgroundWriter.commitBlock: Error when writing to pipe.\n");
		}
	}

	bool isRunning() const {
		return running;
	}

	bool threadWasStarted() const {
		return threadStarted;
	}

	const std::string& getFilename() const {
		return filename;
	}

	protected:

	bool createAndWriteHeader() {
		filename.resize(lenBasename);
		if (fileCounter == 0) {
			filename.append(".gdf");
		} else {
			char suffix[12];
			snprintf(suffix,12,"_%i.gdf", fileCounter);
			filename.append(suffix);
		}
		if (!gdfWriter.createAndWriteHeader(filename.c_str())) {
			fprintf(stderr, "Could not open GDF file %s for writing\n", filename.c_str());
			return false;
		}
		return true;
	}


	void savingThreadFunc() {
		int64_t fileSize = 256*(1+nChans);

		if (!createAndWriteHeader()) return;

		running = true;

		while (1) {
			int newSamplesA, newSamplesB;
			To *rbPtr;
			int64_t newSize, newWritePos;

			int n = locPipe.read(sizeof(int64_t), &newWritePos);

			if (n!=sizeof(int64_t)) {
				// this should never happen for blocking sockets/pipes
				fprintf(stderr, "Unexpected error in pipe communication\n");
				gdfWriter.close();
				break;
			}
			if (newWritePos < 0) {
				gdfWriter.close();
				if (newWritePos == -2) {
					printf("Stopping GDF writing and killing myself...\n");
					gdfWriter.close();
					delete this;
					return;
				} else {
					printf("\nSaving thread received %i - exiting...\n", (int) newWritePos);
					break;
				}
			}

			int writePtr = newWritePos % rbSize;
			int readPtr  = rbReadPos % rbSize;

			//printf("Saving %8lli [%i; %i(\n", newWritePos, readPtr, writePtr);

			rbPtr = rbData + readPtr*nChans;
			if (writePtr > readPtr) {
				newSamplesA = writePtr - readPtr;
				newSamplesB = 0;
			} else {
				newSamplesA = rbSize - readPtr;
				newSamplesB = writePtr;
			}

			int64_t addSize = (newSamplesA+newSamplesB) * nChans * sizeof(To);

			newSize = fileSize + addSize;
			if (newSize > maxFileSize) {
				gdfWriter.close();
				fileCounter++;
				if (!createAndWriteHeader()) break;
				newSize = 256*(1+nChans) + addSize;
			}

			gdfWriter.addSamples(newSamplesA, rbPtr);
			if (newSamplesB > 0) {
				gdfWriter.addSamples(newSamplesB, rbData);
			}
			rbReadPos = newWritePos;
			readPtr   = writePtr;
			fileSize  = newSize;
		}
		running = false;
	}

	static void *staticSavingThreadFunction(void *arg) {
		if (arg == 0) return NULL;

		GDF_BackgroundWriter<To> *GOW = (GDF_BackgroundWriter<To> *) arg;
		GOW->savingThreadFunc();
		return NULL;
	}

	GDF_Writer gdfWriter;
	int nChans;
	bool running, threadStarted;

	To *rbData;
	int rbSize;
	int64_t rbReadPos, rbWritePos;

	LocalPipe locPipe;
	pthread_t savingThread;

	int64_t maxFileSize;
	std::string filename;
	int lenBasename;
	int fileCounter;
};

#endif
