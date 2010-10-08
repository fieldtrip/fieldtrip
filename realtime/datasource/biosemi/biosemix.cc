#include <BioSemiClient.h>
#include <mex.h>

BioSemiClient *BS = NULL;
double tLastAccess;
int triggerState = 0;
int fSample, speedMode;
int numEEG, numAIB;

void exitFun() {
	if (BS) {
		mexPrintf("Shutting down BIOSEMI driver\n");
		delete BS;
	}
}

void initialize() {
	BS = new BioSemiClient();
	if (!BS->isDriverOk()) {
		delete BS;
		BS = NULL;
		mexErrMsgTxt("Could not initialize BIOSEMI driver\n");
	}
	if (!BS->openDevice()) {
		delete BS;
		BS = NULL;
		mexErrMsgTxt("Loaded driver, but cannot open BIOSEMI device\n");
	}
	mexAtExit(exitFun);
	tLastAccess = BS->getCurrentTime();
	numEEG      = BS->getNumChannels();
	numAIB      = BS->getNumChanAIB();
	speedMode   = BS->getSpeedMode();
	fSample     = BS->getSamplingFreq();
}

mxArray *createInfo() {
	const char *field_names[] = {"speedMode","numEEG","numAIB","fSample"};
	mxArray *A;
	
    A = mxCreateStructMatrix(1, 1, 4, field_names);
	                
	mxSetFieldByNumber(A, 0, 0, mxCreateDoubleScalar((double)speedMode));
	mxSetFieldByNumber(A, 0, 1, mxCreateDoubleScalar((double)numEEG));
	mxSetFieldByNumber(A, 0, 2, mxCreateDoubleScalar((double)numAIB));
	mxSetFieldByNumber(A, 0, 3, mxCreateDoubleScalar((double)fSample));
	return A;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	BioSemiBlock block;
	double tNow;
	
	if (BS==NULL) initialize();
	
	tNow = BS->getCurrentTime();
	
	if (nrhs==0) {
		plhs[0] = createInfo();
		BS->checkNewBlock(block); // this is for resetting the blocks
		tLastAccess = tNow;
		return;
	} else {
		bool newBlock;
		int dNumEEG, dNumAIB, dChan;
		int *destPtr;
		const double *pn;
		
		if (!mxIsDouble(prhs[0]) || mxGetM(prhs[0])*mxGetN(prhs[0])!=2) {
			mexErrMsgTxt("Argument, if given, must be 2-element vector [numEEG numAIB]");
		}
		pn = mxGetPr(prhs[0]);
		dNumEEG = (int) pn[0];
		dNumAIB = (int) pn[1];
		if (dNumEEG < 0 || dNumEEG > numEEG) mexErrMsgTxt("Desired number of EEG channels is out of range");
		if (dNumAIB < 0 || dNumAIB > numAIB) mexErrMsgTxt("Desired number of AIB channels is out of range");
		
		newBlock = BS->checkNewBlock(block);
		if (!newBlock) {
			plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
			return;
		}
		
		int predSamples = (int) ((tNow-tLastAccess)*fSample);
		tLastAccess = tNow;

		// simple heuristic to detect timeouts due to wrapped around ring buffer
		if (predSamples > 2*block.numSamples) {
			mexErrMsgTxt("BIOSEMI device: timeout detected - type 'clear biosemix' to force re-loading the driver.");
		}
		if (block.numSamples != block.numInSync) {
			mexErrMsgTxt("BIOSEMI device out of sync!\n");
		}
		plhs[0] = mxCreateNumericMatrix(1+dNumEEG+dNumAIB, block.numSamples, mxINT32_CLASS, mxREAL);
		destPtr = (int *) mxGetData(plhs[0]);
		dChan = 1+dNumEEG+dNumAIB;
		for (int j=0;j<block.numSamples;j++) {
			int idx = block.startIndex + 1 + j*block.stride;
			int status = BS->getValue(idx++);
			
			// only spit out the 16-bit trigger channel pattern
			*destPtr++ = (status & 0x00FFFF00) >> 8;
			
			for (int i=0;i<dNumEEG;i++) {
				*destPtr++ = BS->getValue(idx+i);
			}
			idx+=numEEG;
			for (int i=0;i<dNumAIB;i++) {
				*destPtr++ = BS->getValue(idx+i);
			}
		}
	}
}
