#include <FtBuffer.h>
#include <socketserver.h>
#include <MultiChannelFilter.h>
#include <StringServer.h>
#include <GDF_BackgroundWriter.h>
#include <SignalConfiguration.h>
#include <assert.h>

#ifndef __OnlineDataManager_h
#define __OnlineDataManager_h

/** To is type of original data, e.g. as coming out of the AD-converter
    Ts is the data type used for streaming, e.g. float

    The GDF_Type for writing to disk and the FieldTrip data type for streaming
    will be automatically deduced from these template parameters.
 */

template <typename To, typename Ts>
class OnlineDataManager : public StringRequestHandler {
    typedef char labelType[16];

public:

    OnlineDataManager(int nStatus, int nCont, float fSample) {
      OnlineDataManager(nStatus, nCont, fSample, fSample);
    }

    OnlineDataManager(int nStatus, int nCont, float fSample, float fSampleSaving) {
        this->nStatus = nStatus;
        this->nCont   = nCont;
        this->fSample = fSample;
        this->fSampleSaving = fSampleSaving;

        // automatically determine GDF_Type and Fieldtrip data type from template parameters
        gdfType = GDF_Header::getType((To) 0);
        ftType  = FtDataType::getType((Ts) 0);

        assert(ftType != DATATYPE_UNKNOWN);

        sampleBlock = new FtSampleBlock(ftType);

        pBlock = 0;
        allocSizeBlock = 0;
        auxVec  = new Ts[nCont];
        auxVec2 = new Ts[nCont];

        gdfPhysMin = new double[nStatus + nCont];
        gdfPhysMax = new double[nStatus + nCont];
        gdfDigMin  = new double[nStatus + nCont];
        gdfDigMax  = new double[nStatus + nCont];
        gdfPhysDimCode = new uint16_t[nStatus + nCont];

        if (nStatus > 0) {
            statusLabels = new labelType[nStatus];
        } else {
            statusLabels = 0;
        }

        int sizeTo;
        double minV,maxV;

        sizeTo = GDF_Writer::getSizeAndRangeByType(gdfType, minV, maxV);

        if (sizeTo != sizeof(To)) {
            fprintf(stderr, "Warning: GDF type specified in constructor does not match template parameter.\n");
        }

        for (int i=0;i<nStatus+nCont;i++) {
            gdfPhysMin[i] = minV;
            gdfPhysMax[i] = maxV;
            gdfDigMin[i]  = minV;
            gdfDigMax[i]  = maxV;
            gdfPhysDimCode[i] = 0;
        }

        for (int i=0;i<nStatus;i++) {
            sprintf(statusLabels[i], "Status_%i", i+1);
            gdfPhysDimCode[i] = GDF_DIMLESS;
        }

        offset = new Ts[nCont];
        slope  = new Ts[nCont];
        for (int i=0;i<nCont;i++) {
            offset[i] = 0.0;
            slope[i]  = 1.0;
        }

        ftSocket = -1;
        ftServer = 0;

        sampleCounter = 0;
        skipSamples = 0;
        skipSamples2 = 0;

        curWriter = 0;
        lpFilter = 0;
        lpFilter2 = 0;

        streamingEnabled = false;
        savingEnabled = false;
    }

    virtual ~OnlineDataManager() {
        // TODO: some more of this

        // stop buffer server, if spawned
        if (ftServer) ft_stop_buffer_server(ftServer);
        // FtConnection is cleaned up automatically, if needed

        // clean up GDF writer, if needed
        if (curWriter) {
            curWriter->stopSync();
            delete curWriter;
        }
        // clean up variables
        delete lpFilter;
        delete lpFilter2;
        delete[] auxVec;
        delete[] auxVec2;
        delete[] pBlock;
        delete[] gdfPhysMin;
        delete[] gdfPhysMax;
        delete[] gdfDigMin;
        delete[] gdfDigMax;
        delete[] gdfPhysDimCode;
        delete[] offset;
        delete[] slope;
        if (nStatus > 0) {
            delete[] statusLabels;
        }
        delete sampleBlock;
    }

    virtual std::string handleStringRequest(const std::string& request) {
        static const std::string ok("OK\n");
        static const std::string unknown("ERROR: UNKNOWN COMMAND\n");
        static const std::string malform("ERROR: MALFORMED COMMAND\n");
        static const std::string stopFirst("ERROR: STOP FIRST\n");
        static const std::string emptyFilename("ERROR: NO FILENAME SET\n");
        static const std::string chanOutOfRange("ERROR: CHANNEL INDEX OUTSIDE HARDWARE LIMITS\n");
        static const std::string cannotSave("ERROR: COULD NOT START SAVING\n");
        char response[200];
        static const char trueStr[] = "true";
        static const char falseStr[] = "false";
        static const char noFile[] = "";
        const char *saveTrueFalse;
        const char *filenameStr;

        int target = 0; // 1 = STREAM, 2= SAVE
        unsigned int pos = 0;

        std::string token1 = StringServer::getNextToken(request, pos);
        if (token1.empty()) return ok; // empty command -> nothing to do

        if (token1.compare("STREAM") == 0) {
            target = 1;
        } else if (token1.compare("SAVE") == 0) {
            target = 2;
        } else if (token1.compare("STATUS") == 0) {
            if (savingEnabled && curWriter!=NULL && curWriter->isRunning()) {
                saveTrueFalse = trueStr;
                filenameStr   = curWriter->getFilename().c_str();
            } else {
                saveTrueFalse = falseStr;
                filenameStr   = noFile;
            }

            snprintf(response, sizeof(response)-1,
                     "numacquired=%d numsaved=%d numstreamed=%d saving=%s savingto=\"%s\" downsample=%d bandwidth=%f bworder=%d\n",
                     nCont, // number of continuous channels from hardware
                     signalConf.getSavingSelection().getSize(),
                     signalConf.getStreamingSelection().getSize(),
                     saveTrueFalse,
                     filenameStr,
                     signalConf.getDownsampling(),
                     signalConf.getBandwidth(),
                     signalConf.getOrder());

            return static_cast<std::string>(response);
        } else {
            return unknown;
        }

        std::string token2 = StringServer::getNextToken(request, pos);

        if (token2.compare("START") == 0) {
            if (!StringServer::getNextToken(request, pos).empty()) return malform;
            if (target == 1) {
                if (streamingEnabled) return ok; // silently ignore that it's already running
                enableStreaming();
            } else {
                if (savingEnabled) return ok; // silently ignore that it's already running
                if (curFilename.empty()) return emptyFilename;
                if (!enableSaving()) return cannotSave;
            }
            return ok;
        } else if (token2.compare("STOP") == 0) {
            if (!StringServer::getNextToken(request, pos).empty()) return malform;
            if (target == 1) {
                if (!streamingEnabled) return ok; // silently ignore that it's already stopped
                disableStreaming();

            } else {
                if (!savingEnabled) return ok; // silently ignore that it's already stopped
                disableSaving();
            }
            return ok;
        } else if (token2.compare("SELECT") == 0) {
            if ((target==1 && streamingEnabled) || (target==2 && savingEnabled)) {
                return stopFirst;
            }
            ChannelSelection cs;
            if (!cs.parseString(request.size() - pos, request.data() + pos)) return malform;
            if (cs.getMaxIndex() >= nCont) return chanOutOfRange;

            if (target==1) {
                signalConf.setStreamingSelection(cs);
            } else {
                signalConf.setSavingSelection(cs);
            }
            return ok;
        } else if (target == 2 && token2.compare("FILE") == 0) {
            if (savingEnabled) return stopFirst;

            std::string filename = StringServer::getNextToken(request, pos);
            if (!StringServer::getNextToken(request, pos).empty()) return malform;

            setFilename(filename);
            return ok;
        } else if (target == 1 && token2.compare("FILTER") == 0) {
            if (streamingEnabled) return stopFirst;

            double bandwidth = 0.0;
            int order = 0;
            int factor = 0;

            if (convertToDouble(StringServer::getNextToken(request, pos), bandwidth)
                && (bandwidth >= 0)
                && convertToInt(StringServer::getNextToken(request, pos), order)
                && (order >= 0)
                && convertToInt(StringServer::getNextToken(request, pos), factor)
                && (factor >= 1)
                && StringServer::getNextToken(request, pos).empty()) {

                signalConf.setDownsampling(factor);
                if (bandwidth < 0.5*fSample) {
                    signalConf.setBandwidth(bandwidth);
                    signalConf.setOrder(order);
                } else {
                    signalConf.setOrder(0);
                }
                configureStreaming();

                return ok;
            } else {
                return malform;
            }
        } else if (token2.compare("STATUS") == 0) {
            if (target == 1) {
                // STREAM STATUS
                snprintf(response, sizeof(response)-1,
                         "numacquired=%d numstreamed=%d downsample=%d bandwidth=%f bworder=%d\n",
                         nCont, // number of continuous channels from hardware
                         signalConf.getStreamingSelection().getSize(),
                         signalConf.getDownsampling(),
                         signalConf.getBandwidth(),
                         signalConf.getOrder());
            } else {
                // SAVE STATUS
                if (savingEnabled && curWriter!=NULL && curWriter->isRunning()) {
                    saveTrueFalse = trueStr;
                    filenameStr   = curWriter->getFilename().c_str();
                } else {
                    saveTrueFalse = falseStr;
                    filenameStr   = noFile;
                }

                snprintf(response, sizeof(response)-1,
                         "numacquired=%d numsaved=%d saving=%s savingto=\"%s\"\n",
                         nCont, // number of continuous channels from hardware
                         signalConf.getSavingSelection().getSize(),
                         saveTrueFalse,
                         filenameStr);
            }
            return static_cast<std::string>(response);
        }

        return unknown;
    }

    /** Set label of i-th status channel (starting with 0).
     The default is Status_1, Status_2, ...
     A maxiumum of 16 characters is allowed (as in GDF).
     */
    void setStatusLabel(int i, const char *name) {
        if (i<0 || i>=nStatus) return;
        strncpy(statusLabels[i], name, sizeof(labelType));
    }

    /** Set physical dimension code for all continuous channels
     for writing to GDF. See GdfWriter.h for an enumeration
     of possible values.
     */
    void setPhysicalDimCode(uint16_t code) {
        for (int i=0;i<nCont;i++) {
            gdfPhysDimCode[nStatus + i] = code;
        }
    }

    /** Set physical dimension code for the i-th continuous channel,
     starting from 0, for writing to GDF. See GdfWriter.h for an
     enumeration	of possible values.
     */
    void setPhysicalDimCode(int i, uint16_t code) {
        if (i>=0 && i<nCont) {
            gdfPhysDimCode[nStatus + i] = code;
        }
    }

    /** Set the physical limits for all continuous channels for writing	to GDF.
     */
    void setPhysicalLimits(double minV, double maxV) {
        for (int i=0;i<nCont;i++) {
            gdfPhysMin[nStatus + i] = minV;
            gdfPhysMax[nStatus + i] = maxV;
        }
    }

    /** Set the physical limits for the i-th continuous channel (starting from 0)
     for writing	to GDF.
     */
    void setPhysicalLimit(int i, double minV, double maxV) {
        if (i<0 || i>=nCont) return;
        gdfPhysMin[nStatus + i] = minV;
        gdfPhysMax[nStatus + i] = maxV;
    }

    /** Set the digital limits for all continuous channels for writing to GDF.
     */
    void setDigitalLimits(double minV, double maxV) {
        for (int i=0;i<nCont;i++) {
            gdfDigMin[nStatus + i] = minV;
            gdfDigMax[nStatus + i] = maxV;
        }
    }

    /** Set the digital limits for the i-th continuous channel (starting from 0)
     for writing	to GDF.
     */
    void setDigitalLimit(int i, double minV, double maxV) {
        if (i<0 || i>=nCont) return;
        gdfDigMin[nStatus + i] = minV;
        gdfDigMax[nStatus + i] = maxV;
    }

    /** Set the affine transformation of the i-th continuous channel. The signal
     will first be subtracted 'offset', then multiplied by 'slope'. Note that
     this transformation does not affect the data that is saved to GDF files,
     but only the streamed channels.
     */
    void setSlopeAndOffset(int i, Ts slope, Ts offset) {
        if (i<0 || i>=nCont) return;
        this->slope[i] = slope;
        this->offset[i] = offset;
    }

    /** Set the affine transformation for all continuous channels. The signals
     will first be subtracted 'offset', then multiplied by 'slope'. Note that
     this transformation does not affect the data that is saved to GDF files,
     but only the streamed channels.
     */
    void setSlopeAndOffset(Ts slope, Ts offset) {
        for (int i=0;i<nCont;i++) {
            this->slope[i] = slope;
            this->offset[i] = offset;
        }
    }

    /** Connect to a remove FieldTrip buffer server on the given address. TCP-IP
     address must be given in the form hostname:port. Otherwise, the address
     will be interpreted as a UNIX domain socket name. Returns true on success,
     false if errors occured.
     */
    bool connectToServer(const char *address) {
        if (ftSocket != -1) return false;

        if (!ftConnection.connect(address)) return false;
        ftSocket = ftConnection.getSocket();
        return true;
    }

    /** Connect to a remove FieldTrip buffer server on the given address over
     TCP-IP. Returns true on success, false if errors occured.
     */
    bool connectToServer(const char *hostname, int port) {
        if (ftSocket != -1) return false;

        if (!ftConnection.connectTcp(hostname, port)) return false;
        ftSocket = ftConnection.getSocket();
        return true;
    }

    /** Spawn an internal buffer server on the given TCP-IP port, and subsequently
     use "dma" requests to the buffer. Returns true on success, false if errors
     occured (i.e. port not available).
     */
    bool useOwnServer(int port) {
        if (ftSocket != -1) return false;

        ftServer = ft_start_buffer_server(port, NULL, NULL, NULL);
        if (ftServer == 0) return false;

        ftSocket = 0; // => dma
        return true;
    }

    /** Read channel selections and other parameters from the given configuration file,
     and configure the streaming of data on success. Returns the number of parsing
     errors (0 on success).
     */
    int configureFromFile(const char *filename) {
        int numErr = signalConf.parseFile(filename);
        if (numErr == 0) {
            configureStreaming();
        }
        return numErr;
    }

    /** This function should be called by the acquisition driver to ask the OnlineDataManager
     to prepare a new block of N samples. It will try to allocate the necessary space for
     it (all HW channels), and return a pointer to that memory.
     */
    To *provideBlock(int N) {
        int needed = N * (nStatus + nCont);
        if (needed > allocSizeBlock) {
            delete[] pBlock;
            pBlock = new To[needed];
            allocSizeBlock = needed;
        }
        nThisBlock = N;

        eventList.clear();
        return pBlock;
    }

    /** This function should be called by the acquisition driver after data has been filled
     into the provided block, in order to stream out and save the selected channels.
     Returns true on success, false if errors occured.
     */
    bool handleBlock() {
        if (streamingEnabled) {
            if (!handleStreaming()) return false;
        }
        if (savingEnabled) {
            return handleSaving();
        }
        return true;
    }

    /** Call this to enable saving to a GDF file. You must have called "setFilename"
     before. Once a filename is set, multiple calls to enableSaving() will increase
     a file counter which is appended to the name of the target GDF file.
     Returns true on success, false on errors.
     */
    bool enableSaving() {
        if (curFilename.empty()) return false;
        if (curWriter == 0) {
            configureSaving();
            if (curWriter == 0) return false;
            printf("Reconfigured saving\n");
        }
        if (fileCounter == 0) {
            curWriter->start(curFilename.c_str());
            printf("Saving to: %s\n",curFilename.c_str());
        } else {
            char *compName = new char[curFilename.size() + 20];
            sprintf(compName, "%s_S%i", curFilename.c_str(), fileCounter);
            // printf("COMPOSED: %s\n", compName);
            curWriter->start(compName);
            printf("Saving to: %s\n",compName);
            delete[] compName;
        }
        fileCounter++;
        savingEnabled = true;
        return true;
    }

    /** Call this to disable saving to GDF. */
    void disableSaving() {
        savingEnabled = false;
        if (!curWriter) return;
        curWriter->stopAsync();
        curWriter = 0;
    }

    /** Call this to enable streaming to a FieldTrip buffer. */
    bool enableStreaming() {
        if (streamingEnabled) return true; // silently ignore
        if (!writeHeader()) return false;
        streamingEnabled = true;
        return true;
    }

    /** Call this to disable streaming */
    void disableStreaming() {
        // if (!streamingEnabled) return;
        streamingEnabled = false;
    }

    /** Returns a reference to the event list, which the acquisition driver
     should add events to. The event list will be flushed after every call
     to handleBlock()
     */
    FtEventList& getEventList() {
        return eventList;
    }

    /** Returns a pointer to the current GDF_Writer object, or NULL if there is
     no instance at the moment. This should probably not be called by users
     without a very good cause.
     */
    GDF_Writer *getCurrentGDF() {
        if (curWriter) return curWriter->gdf();
        return 0;
    }

    /** Set the filename for writing to GDF. Changes will note take effect
     before calling enableSaving()
     */
    void setFilename(const std::string& filename) {
        curFilename = filename;
        fileCounter = 0;
    }


    /** Retrieve reference to SignalConfiguration object */
    const SignalConfiguration& getSignalConfiguration() { return signalConf; }

    /** Configure signals for streaming and saving from C++ object */

    bool setSignalConfiguration(const SignalConfiguration& newCfg) {
        if (savingEnabled) return false;
        if (streamingEnabled) return false;

        if (newCfg.getMaxSavingChannel() >= nCont) return false;
        if (newCfg.getMaxStreamingChannel() >= nCont) return false;

        signalConf = newCfg;

        // check if selected bandwidth is beyond Nyquist
        // and if so, disable filtering silently
        float bw = signalConf.getBandwidth();
        if (bw < 0 || bw >= 0.5*fSample) {
            signalConf.setOrder(0);
        }
        configureStreaming();
        return true;
    }

protected:

    /** Helper function for writing a header to the FieldTrip buffer that corresponds
     to the current channel selection for streaming. Returns true on success,
     false on error.
     */
    bool writeHeader() {
        const ChannelSelection& streamSel = signalConf.getStreamingSelection();
        FtBufferRequest req;
        char *chunk_data;
        int N=0,P;

        for (int n=0;n<streamSel.getSize();n++) {
            int Ln = strlen(streamSel.getLabel(n))+1;
            N+=Ln;
        }
        chunk_data = new char[N];

        P=0;
        for (int n=0;n<streamSel.getSize();n++) {
            int Ln = strlen(streamSel.getLabel(n))+1;
            memcpy(chunk_data + P, streamSel.getLabel(n), Ln);
            P+=Ln;
        }

        req.prepPutHeader(streamSel.getSize(), ftType, fSample / signalConf.getDownsampling());
        req.prepPutHeaderAddChunk(FT_CHUNK_CHANNEL_NAMES, N, chunk_data);

        delete[] chunk_data;

        int err = clientrequest(ftSocket, req.out(), resp.in());
        if (err || !resp.checkPut()) {
            fprintf(stderr, "Could not write header to FieldTrip buffer\n");
            return false;
        }
        sampleCounter = 0;
        skipSamples = 0;
        skipSamples2 = 0;
        return true;
    }

    /** Called by handleBlock() to deal with streaming out samples and events.
     The raw data are first transformed by subtracting offsets and multiplying
     slope factors. If selected, the signal will then be filtered and optionally
     downsampled.
     */
    bool handleStreaming() {
        int stride = nStatus + nCont;
        int err;
        const ChannelSelection& streamSel = signalConf.getStreamingSelection();
        int nStream  = streamSel.getSize();

        // write events, if any
        if (eventList.count() > 0) {
            eventList.transform(sampleCounter, signalConf.getDownsampling());
            err = clientrequest(ftSocket, eventList.asRequest(), resp.in());
            if (err || !resp.checkPut()) {
                fprintf(stderr, "Could not write events to FieldTrip buffer.\n");
                return false;
            }
        }

        sampleCounter += nThisBlock; // sampleCounter ticks at original speed

        // write samples, if channels selected
        if (nStream == 0) return true;

        int deci        = signalConf.getDownsampling();
        int numThisTime = (nThisBlock - skipSamples + deci - 1)/deci;

        Ts *dest = (Ts *) sampleBlock->getMatrix(nStream, numThisTime);

        if (lpFilter) {
            for (int j=0;j<nThisBlock;j++) {
                To *src = pBlock + nStatus + j*stride;
                for (int i=0;i<nStream;i++) {
                    int idx = streamSel.getIndex(i);
                    auxVec[i] = slope[idx]*(src[idx] - offset[idx]);
                }

                if (skipSamples == 0) {
                    lpFilter->process(dest, auxVec);
                    dest += nStream;
                } else {
                    lpFilter->process(auxVec);
                }
                if (--skipSamples < 0) skipSamples = deci-1;
            }
        } else {
            for (int j=0;j<nThisBlock;j++) {
                To *src = pBlock + nStatus + j*stride;
                if (skipSamples == 0) {
                    for (int i=0;i<nStream;i++) {
                        int idx = streamSel.getIndex(i);
                        dest[i] = slope[idx]*(src[idx] - offset[idx]);
                    }
                    dest += nStream;
                }
                if (--skipSamples < 0) skipSamples = deci-1;
            }
        }

        if (numThisTime > 0) { // only send packet if there is actual data.
            err = clientrequest(ftSocket, sampleBlock->asRequest(), resp.in());
            if (err || !resp.checkPut()) {
                fprintf(stderr, "Could not write samples to FieldTrip buffer\n");
                return false;
            }
        }
        return true;
    }


    /** Called by handleBlock to deal with saving data to GDF. Actually this function
     doesn't save to disk itself, but copies the relevant channels to the
     internal ring buffer of the current GDF_BackgroundWriter instance. */
    bool handleSaving() {
        int stride = nStatus + nCont;
        const ChannelSelection& saveSel = signalConf.getSavingSelection();
        int nSave  = saveSel.getSize();

        if (curWriter == 0) return false;

        int deci        = (int)(fSample/fSampleSaving);
        int numThisTime = (nThisBlock - skipSamples2 + deci - 1)/deci;

        Ts *dest_filt = NULL; // to hold the filtered data

        if (numThisTime > 0) {
            if (!curWriter->checkFreeBlock(numThisTime)) {
                fprintf(stderr, "Error: saving data thread does not keep up with load\n");
                return false;
            }
        }
        if (lpFilter2) {
            // without nStatus: STATUS channel not filtered!
            dest_filt = (Ts *) sampleBlock->getMatrix(nSave, numThisTime);

            for (int j=0;j<nThisBlock;j++) {
                To *src = pBlock + nStatus + j*stride;
                for (int i=0;i<nSave;i++) {
                    int idx = saveSel.getIndex(i);
                    auxVec2[i] = slope[idx]*(src[idx] - offset[idx]);
                }
                // STATUS channel(s) need(s) to be duplicated from calling client in order to make sure
                // no trigger information is missing after decimation
                if (skipSamples2 == 0) {
                    lpFilter2->process(dest_filt, auxVec2);
                    // save the data ***
                    src = pBlock + j*stride; // src2: points to STATUS channel for this sample in the data
                    To *dest = curWriter->getSampleSlot(); // dest for saving this data sample (for all channels)
                    for (int i=0;i<nStatus;i++) {
                        *dest++ = *src++;
                    }
                    for (int i=0;i<nSave;i++) {
                        int idx = saveSel.getIndex(i);
                        // convert EEG sample back to original digital values
                        dest[i] = (To)((Ts)dest_filt[i]/slope[idx] + offset[idx]);
                    }
                    dest_filt += nSave;
                } else {
                    lpFilter2->process(auxVec2);
                }
                if (--skipSamples2 < 0) skipSamples2 = deci-1;
            }
        } else {
            // no filtering required
            for (int j=0;j<nThisBlock;j++) {
                To *dest = curWriter->getSampleSlot();
                To *src = pBlock + j*stride;
                if (skipSamples2 == 0) {
                    for (int i=0;i<nStatus;i++) {
                        *dest++ = *src++;
                    }
                    for (int i=0;i<nSave;i++) {
                        dest[i] = src[saveSel.getIndex(i)];
                    }
                }
                if (--skipSamples2 < 0) skipSamples2 = deci-1;
            }
        }
        if (numThisTime > 0) { // save if there is actual data.
            curWriter->commitBlock();
        }
        return true;
    }


    /** Create a new GDF_BackgroundWriter with the proper channels selected
     The old one (if any) will automatically terminate and delete itself
     */
    void configureSaving() {
        const ChannelSelection& saveSel = signalConf.getSavingSelection();
        int nSave = saveSel.getSize();
        if (nStatus + nSave == 0) return;

        // NOTE: filter settings fixed to downsample to sr=512Hz
        delete lpFilter2;
        // filter (sr=512Hz, decimation=4)

        if (fSample != fSampleSaving) {
            lpFilter2 = new MultiChannelFilter<Ts,Ts>(signalConf.getSavingSelection().getSize(), 4);
            lpFilter2->setButterLP(fSampleSaving/fSample);
        } else {
            lpFilter2 = 0;
        }

        curWriter = new GDF_BackgroundWriter<To>(nStatus + nSave, fSampleSaving, gdfType);

        for (int i=0;i<nStatus;i++) {
            curWriter->gdf().setLabel(i, statusLabels[i]);
            curWriter->gdf().setPhysicalLimits(i, gdfPhysMin[i], gdfPhysMax[i]);
            curWriter->gdf().setDigitalLimits(i, gdfDigMin[i], gdfDigMax[i]);
            curWriter->gdf().setPhysDimCode(i, gdfPhysDimCode[i]);
        }
        for (int i=0;i<nSave;i++) {
            int idx = nStatus + saveSel.getIndex(i);
            curWriter->gdf().setLabel(nStatus+i, saveSel.getLabel(i));
            curWriter->gdf().setPhysicalLimits(nStatus+i, gdfPhysMin[idx], gdfPhysMax[idx]);
            curWriter->gdf().setDigitalLimits(nStatus+i, gdfDigMin[idx], gdfDigMax[idx]);
            curWriter->gdf().setPhysDimCode(nStatus+i, gdfPhysDimCode[idx]);
        }
        printf("\nSampling frequency (saving)......: %.0f Hz\n", fSampleSaving);
        printf("Number of saved channels.........: %d\n",nSave+nStatus);
    }

    /** Set up lowpass filter with the correct characteristics and number of channels.
     Delete the current filter first, if any.
     */
    void configureStreaming() {
        const ChannelSelection& streamSel = signalConf.getStreamingSelection();
        int nStream = streamSel.getSize();
        delete lpFilter;

        if (signalConf.getOrder() > 0) {
            lpFilter = new MultiChannelFilter<Ts,Ts>(signalConf.getStreamingSelection().getSize(), signalConf.getOrder());
            lpFilter->setButterLP(signalConf.getBandwidth() / (0.5*fSample));
        } else {
            lpFilter = 0;
        }
        printf("Sampling frequency (streamed)....: %.0f Hz\n", fSample/signalConf.getDownsampling());
        printf("Number of streamed channels......: %d\n", nStream);
    }

    ////////////////////////////////////////////////////////////////
    // Class members
    ////////////////////////////////////////////////////////////////
    GDF_BackgroundWriter<To> *curWriter;	/**< currently active GDF writer (in background thread) */
    MultiChannelFilter<Ts,Ts> *lpFilter;	/**< currently active low-pass filter for streamed data */
    MultiChannelFilter<Ts,Ts> *lpFilter2;	/**< currently active low-pass filter for saved data */

    UINT32_T ftType;	/**< FieldTrip buffer data type */
    GDF_Type gdfType;	/**< GDF data type */
    int nStatus, nCont; /**< Number of status + continuously sampled channels */
    float fSample;		/**< Sampling rate */
    float fSampleSaving;/**< Sampling rate for saved data */

    int nThisBlock;		/**< Number of samples in this block (as requested by provideBlock) */
    int allocSizeBlock; /**< Size of buffer that is allocated for providing blocks */
    To *pBlock;			/**< Points to buffer that is allocated for providing blocks */
    Ts *auxVec;			/**< Big enough for keeping one sample of data (streaming format) */
    Ts *auxVec2;		/**< Big enough for keeping one sample of data (saving format) */
    Ts *offset;         /**< Offset subtracted from raw data before streaming */
    Ts *slope;          /**< Factor to multiply data with before streaming (after subtracting offset) */

    int ftSocket;		/**< The FT buffer socket identifier or 0 for dmarequests, -1 for none */
    int sampleCounter;	/**< Number of samples streamed out since last writeHeader */
    int skipSamples;	/**< Helper variable to keep track of downsampling operation for steamed data */
    int skipSamples2;	/**< Helper variable to keep track of downsampling operation for saved data */

    double *gdfPhysMin;		/**< Physical minimum for status + continuous channels, written to GDF */
    double *gdfPhysMax;		/**< Physical maximum for status + continuous channels, written to GDF */
    double *gdfDigMin;		/**< Digital minimum for status + continuous channels, written to GDF */
    double *gdfDigMax;      /**< Digital maximum for status + continuous channels, written to GDF */
    uint16_t *gdfPhysDimCode; 	/**< Dimension code for status + continuous channels, written to GDF */
    labelType *statusLabels;	/**< Labels of the status channels, up to 16 characters each, written to GDF */

    FtConnection ftConnection;	/**< Handles the connection to the FieldTrip buffer (either socket or dma) */
    FtBufferResponse resp;		/**< Receives responses from the buffer server */
    FtEventList eventList;		/**< Used for writing events to the buffer server, is flushed after each handleBlock() */
    FtSampleBlock *sampleBlock;	/**< Used for writing data to the buffer server */
    ft_buffer_server_t *ftServer;	/**< Handles the server sockets and background threads in case an own server is spawned */

    SignalConfiguration signalConf;	/**< Maintains the channel selection for streaming and saving, as well as a few other parameters */

    bool savingEnabled, streamingEnabled;	/**< Flags that determine the current mode of operation */

    std::string curFilename;	/**< Current filename for writing GDF to disk */
    int fileCounter;			/**< Current running ID for the file name, auto-incremented after every disableSaving() call */
};

#endif
