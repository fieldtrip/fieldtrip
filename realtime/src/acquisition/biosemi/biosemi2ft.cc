/** Biosemi acqusition tool to stream (downsampled) data to a FieldTrip buffer,
	and write full-rate data to one or multiple GDF files.
	(C) 2010 S. Klanke
*/
#include <stdio.h>
#include <BioSemiClient.h>
#include <OnlineDataManager.h>
#include <ConsoleInput.h>
#include <StringServer.h>

StringServer ctrlServ;
ConsoleInput conIn;
char hostname[256];
int port, ctrlPort;
char gdfname[1024];
char cfgname[1024];
bool haveGDF;

void acquisition(BioSemiClient& BS, int numHwChan, int fSample, int fSampleSaving) {
	// number of samples grabbed from USB in total
	int sampleCounter = 0;
	// value of the trigger channel (middle 16 bits of that)
	int triggerState = 0;
    int triggerSize = 0;     /**< size of triggerStates, i.e., number of stored trigger states for each status channel */
    int *triggerStates = NULL;  /**< remember all triggers received and that still need to be processed to ODM */
    int deci = (int)(fSample/fSampleSaving);
    int copied_trigger;
    
	// "CMS in range" flag as reported most recently
	int cmsInRange = 0;
	// number of samples until we're due a report again
	int nsBattery = 0, nsCMS = 0;
	// ODM grabs original 32-bit integer data, scales it to microVolts for streaming
	OnlineDataManager<int, float> ODM(1, numHwChan, (float) fSample, (float) fSampleSaving);
	// just for reporting: print current time relative to time when we started this
	double T0 = BS.getCurrentTime();
	
	if (ODM.configureFromFile(cfgname) != 0) {
		fprintf(stderr, "Configuration file is invalid\n");
		return;
	}
	if (!strcmp(hostname, "-")) {
		if (!ODM.useOwnServer(port)) {
			fprintf(stderr, "Could not spawn buffer server on port %d.\n",port);
			return;
		}
	} else {
		if (!ODM.connectToServer(hostname, port)) {
			fprintf(stderr, "Could not connect to buffer server at %s:%d.\n",hostname, port);
			return;
		}
	}
	ODM.setStatusLabel(0, "STATUS");
	// min:  (-2^31)/8192.0  (int32 -> microVolt)
	// max: (2^31-1)/8192.0  (int32 -> microVolt)
	ODM.setPhysicalLimits(-262144.0, 262143.99987792969);
	// scale streamed data to microvolts
	ODM.setPhysicalDimCode(GDF_MICRO + GDF_VOLT);	
	ODM.setSlopeAndOffset(1.0/8192.0, 0);
	
	if (haveGDF) ODM.setFilename(gdfname);
	
	if (!ODM.enableStreaming()) return;
	
	if (haveGDF) {
		printf("\nPress <Esc> to quit, <S> to enable saving, <D> to disable saving\n\n");
	} else {
		printf("\nPress <Esc> to quit\n\n");
	}
	
	while (1) {
		BioSemiBlock block;
		int *ptr;
		const SignalConfiguration& config = ODM.getSignalConfiguration();
		double batt = 0.0;
		if (conIn.checkKey()) {
			int c = conIn.getKey();
            if (c==27) {
                printf("\n");
                break; // quit
            }
			if (haveGDF) {
				if (c=='s' || c=='S') {
					if (!ODM.enableSaving()) {
						fprintf(stderr, "Cannot enable saving to GDF\n");
						break;
					}
				}
			
				if (c=='d' || c=='D') {
                    printf("\n");
					ODM.disableSaving();
				}
			}
		}
		ctrlServ.checkRequests(ODM);
		
		if (!BS.checkNewBlock(block)) {
			BS.msleep(1);
			continue;
		}
		
		if (numHwChan >= 280) {
			batt = BS.getValue(block.startIndex + 279);
			batt=(batt/2097152)-175.5;
			if (batt < 0) batt = 0;
			if (batt > 100) batt = 100;
		}
		
		printf("T=%8.3f Ptr=%8i, samples=%3i, batt=%5.1f%%\r", BS.getCurrentTime() - T0, block.startIndex, block.numSamples, batt);
		if (block.numSamples != block.numInSync) {
			fprintf(stderr, "USB device out of sync (%i / %i) -- exiting\n", block.numInSync, block.numSamples);
			break;
		}
		
		ptr = ODM.provideBlock(block.numSamples);
		if (ptr == NULL) {
			fprintf(stderr, "Out of memory\n");
			break;
		}
		
		if (nsBattery >= 0) {
			if (nsBattery<block.numSamples) {
				// we're due a report of the battery level
				ODM.getEventList().add(0, "BATTERY", (float) batt);
				//printf("\n-!- Battery level = %4.1f\n", batt);
				// next report in getBatteryRefresh() seconds
				nsBattery = config.getBatteryRefresh() * fSample;
				// if this is zero, user doesn't want events,
				// so prevent further output by setting it to -1
				if (nsBattery == 0) nsBattery = -1; 
			} else {
				nsBattery -= block.numSamples;
			}
		}

        bool printnewline = true;
		for (int j=0;j<block.numSamples;j++) {
			// dest_j points to j-th sample in block provided by ODM
            
			int *dest_j = ptr + j*(1+numHwChan);
			// source offset, j-th sample
			int soff   = block.startIndex + j*block.stride; 
			// status word (incl. triggers)
			int value  = BS.getValue(soff + 1);
            int orig_trigger = value;
			// "CMS within range" flag in Biosemi logic
			int cmsBit = (value & 0x10000000) ? 0:1;
						
			if (nsCMS == 0 || cmsBit != cmsInRange) {
				ODM.getEventList().add(j, "CMS_IN_RANGE", cmsBit);
				//printf("\n-!- CMS in range: %i\n", cmsBit);
				cmsInRange = cmsBit;
				nsCMS = config.getStatusRefresh() * fSample;
			} else {
				--nsCMS;
			}
			
			value = (value & 0x00FFFF00) >> 8;
			if (value && value!=triggerState) {
				int sNow = sampleCounter + j;
				
				if (!config.useSplittedTrigger())  {
					ODM.getEventList().add(j, "TRIGGER", value);
					printf("\n-!- Trigger at sample %i => %i\n", sNow, value);
				} else {
					// find out which of the bytes has changed, and send according trigger
					int hv = value >> 8;
					int lv = value & 0x00FF;
					int ht = triggerState >> 8;
					int lt = triggerState & 0x00FF;
                    
                    if ((hv != ht || lv != lt) && printnewline) { printf("\n"); printnewline=false; }
					if (hv != ht) {
						ODM.getEventList().add(j, config.getHighTriggerName(), hv);
						printf("-!- %s at sample %i => %i\n", config.getHighTriggerName(), sNow, hv);
					}
					if (lv != lt) {
						ODM.getEventList().add(j, config.getLowTriggerName(), lv);
						printf("-!- %s at sample %i => %i\n", config.getLowTriggerName(), sNow, lv);
					}
				}
                
                if (deci>1) {
                    // copy possible STATUS information as many times as decimation will skip STATUS samples
                    // if (sampleCounter modulus deci) is zero, change STATUS sample to next trigger value in the queue
                    // The queue stores triggers that were sent too quickly and would get lost by the decimation
                    triggerSize = triggerSize + 2; // also add intermediate zero trigger value
                    if (triggerSize > 500) {
                        fprintf(stderr,"Cannot process and delay such a large number (%d) of received triggers\n",(int)(triggerSize/2));
                        return;
                    }
                    int *more_numbers = (int*) realloc (triggerStates, triggerSize * sizeof(int));
                    if (more_numbers!=NULL) {
                        triggerStates=more_numbers;
                        triggerStates[triggerSize-2] = orig_trigger; // store original STATUS value
                        triggerStates[triggerSize-1] = 0; // also add inbetween zero STATUS value
                    }
                    else {
                        free (triggerStates);
                        fprintf(stderr,"Error (re)allocating memory");
                        return;
                    }
                }
			}
			triggerState = value;

            if (deci > 1) {
                if ((sampleCounter+j) % deci == 0) {
                    // change trigger for next #deci samples
                    if (triggerSize > 0) {
                        // take from queued triggers
                        copied_trigger = triggerStates[0]; // take oldest
                        // shift/remove items in array
                        for(int p=1; p < triggerSize; p++) triggerStates[p-1] = triggerStates[p];
                        triggerSize = triggerSize - 1;
                        triggerStates = (int*) realloc (triggerStates, triggerSize * sizeof(int));
                        // verify return value != null
                        // ....
                    } else {
                        copied_trigger = 0;
                    }
                }
                dest_j[0] = copied_trigger;
            }
            else {
                // copy status value
                dest_j[0] = orig_trigger;
            }
            
			// copy continuous channels into memory provided by ODM
			// selection, scaling, streaming + saving will be handled from there
			for (int i=0;i<numHwChan;i++) {
				// a memcpy would be a bit faster, but then we'd need to handle 
				// the wrapping-around of the BS ringbuffer here
				dest_j[1+i] = BS.getValue(soff + i + 2);
			}
		}
			
		if (!ODM.handleBlock()) {
			fprintf(stderr, "Error in handling this data block - stopping\n");
			break;
		}
		
		sampleCounter += block.numSamples;
	}
    
    if (triggerStates != NULL) free(triggerStates);
}


int main(int argc, char *argv[]) {
	BioSemiClient BS;
    int decimationSaving;
	
	if (argc<2) {
		printf("Usage: biosemi2ft <config-file> [gdf-file] [hostname=localhost [port=1972 [ctrlPort=8000 [decimation=1]]]\n");
		return 0;
	}
	
	strncpy(cfgname, argv[1], sizeof(cfgname));
#ifdef WIN32
	if (GetFileAttributes(cfgname)==INVALID_FILE_ATTRIBUTES) {
#else
	if (access(cfgname, R_OK)) {
#endif
		fprintf(stderr, "Cannot open configuration file %s\n", cfgname);
		return 1;
	}
	
	if (argc>2) {
		strncpy(gdfname, argv[2], sizeof(gdfname));
		haveGDF = strcmp(gdfname, "-") != 0;
	} else {
		strcpy(gdfname, "-");
		haveGDF = false;
	}
		
	if (argc>3) {
		strncpy(hostname, argv[3], sizeof(hostname));
	} else {
		strcpy(hostname, "localhost");
	}
	
	if (argc>4) {
		port = atoi(argv[4]);
	} else {
		port = 1972;
	}
	
	if (argc>5) {
		ctrlPort = atoi(argv[5]);
	} else {
		ctrlPort = 8000;
	}
	if (!ctrlServ.startListening(ctrlPort)) {
		fprintf(stderr, "Cannot listen on port %d for configuration commands\n", ctrlPort);
		return 1;
	}

    if (argc>6) {
        decimationSaving = atoi(argv[6]);
    } else {
        decimationSaving = 1; // save data to file with original amplifiers sample frequency
    }
    if (decimationSaving < 1 || decimationSaving > 8) {
        fprintf(stderr,"The downsample decimation factor (%d) for saving to file should be in the range [1-8]\n",decimationSaving);
        return 1;
    }
        
	if (!BS.openDevice()) return 1;
	
	printf("Speed mode.......................: %i\n", BS.getSpeedMode());
	printf("Number of channels...............: %i\n", BS.getNumChannels());
	printf("Number of AIB channels...........: %i\n", BS.getNumChanAIB());
	printf("Sampling frequency (amplifier)...: %i Hz\n", BS.getSamplingFreq());
	
	acquisition(BS, BS.getNumChannels() + BS.getNumChanAIB(), BS.getSamplingFreq(), BS.getSamplingFreq()/decimationSaving);
	
	BS.closeDevice();
	return 0;
}
