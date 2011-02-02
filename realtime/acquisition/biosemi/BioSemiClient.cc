#include <BioSemiClient.h>
#include <stdio.h>

#ifndef WIN32
#include <dlfcn.h>
#endif


#define DEVICE_TIMEOUT   1.0		// in seconds

#define MK2_BV        0x80000000	// bit 23/31
#define BATTERY_BV    0x40000000    // bit 22/30
#define SPEED_BIT3    0x20000000	// bit 21/29
#define CMS_RANGE_BV  0x10000000    // bit 20/28
#define SPEED_BIT2    0x08000000	// bit 19/27
#define SPEED_BIT1    0x04000000	// bit 18/26
#define SPEED_BIT0    0x02000000	// bit 17/25
#define EPOCH_BV	  0x01000000    // bit 16/24

// plus 32 AIB channels for speedmode 9
static const int sNumChannelsMk1[9] = { 256, 128,  64,  32, 256, 138, 64, 32, 256};
static const int sNumChannelsMk2[9] = { 608, 608, 608, 608, 280, 152, 88, 56, 280};
static const int sSamplingFreq[9]   = {2048, 4096, 8192, 16384, 2048, 4096, 8192, 16384, 2048};


BioSemiClient::BioSemiClient() {
	driverOk = false;
	deviceOpen = false;
	const char *dllName;
	
#ifdef WIN32
	timeBeginPeriod(1);
	if (sizeof(void *)==8) {
		dllName = "Labview_DLL64.dll";
	} else {
		dllName = "Labview_DLL.dll";
	}
	
	hLib = LoadLibrary(dllName);
	if (hLib == NULL) {
		fprintf(stderr, "Cannot load %s\n", dllName);
		return;
	}
		
	if (!(lv_open_driver_async    = (OPEN_DRIVER_ASYNC_T)    GetProcAddress(hLib,"OPEN_DRIVER_ASYNC"))) return;
	if (!(lv_usb_write            = (USB_WRITE_T)            GetProcAddress(hLib,"USB_WRITE"))) return;
	if (!(lv_read_multiple_sweeps = (READ_MULTIPLE_SWEEPS_T) GetProcAddress(hLib,"READ_MULTIPLE_SWEEPS"))) return;
	if (!(lv_read_pointer         = (READ_POINTER_T)         GetProcAddress(hLib,"READ_POINTER"))) return;
	if (!(lv_close_driver_async   = (CLOSE_DRIVER_ASYNC_T)   GetProcAddress(hLib,"CLOSE_DRIVER_ASYNC"))) return;
	
	// Use VirtualAlloc on Windows as done in Biosemi's synctest code
	ringBuffer = (int32_t *) VirtualAlloc(NULL, BUFFER_SIZE, MEM_COMMIT, PAGE_READWRITE);
	if (ringBuffer != NULL) ZeroMemory(ringBuffer, BUFFER_SIZE);
#else
   
    #ifdef __APPLE__
	if (sizeof(void *)==8) {
		dllName = "./liblabview_dll64.dylib";
	} else {
		dllName = "./liblabview_dll.dylib";
	}
    #else // should be Linux then
	if (sizeof(void *)==8) {
		dllName = "./liblabview_dll64.so";
	} else {
		dllName = "./liblabview_dll.so";
	}
    #endif
	
	hLib = dlopen(dllName, RTLD_LAZY);
	if (hLib == NULL) {
		fprintf(stderr, "Cannot load %s\n", dllName);
		return;
	}
		
	if (!(lv_open_driver_async    = (OPEN_DRIVER_ASYNC_T) 	 dlsym(hLib,"OPEN_DRIVER_ASYNC"))) return;
	if (!(lv_usb_write            = (USB_WRITE_T)            dlsym(hLib,"USB_WRITE"))) return;
	if (!(lv_read_multiple_sweeps = (READ_MULTIPLE_SWEEPS_T) dlsym(hLib,"READ_MULTIPLE_SWEEPS"))) return;
	if (!(lv_read_pointer         = (READ_POINTER_T)         dlsym(hLib,"READ_POINTER"))) return;
	if (!(lv_close_driver_async   = (CLOSE_DRIVER_ASYNC_T)   dlsym(hLib,"CLOSE_DRIVER_ASYNC"))) return;
	
	ringBuffer = (int32_t *) calloc(BUFFER_SIZE/sizeof(int32_t), sizeof(int32_t));
#endif
	if (ringBuffer == NULL) {
		fprintf(stderr, "Cannot allocate %i bytes for internal ringbuffer\n", BUFFER_SIZE);
		return;
	}
	
	driverOk = true;
}

BioSemiClient::~BioSemiClient() {
	if (deviceOpen) closeDevice();
	#ifdef WIN32
	if (hLib != NULL) FreeLibrary(hLib);
	timeEndPeriod(1);
	if (ringBuffer != NULL) VirtualFree(ringBuffer, 0, MEM_RELEASE);
	#else
	// TODO: dlclose
	if (hLib!=NULL) dlclose(hLib);
	if (ringBuffer != NULL) free(ringBuffer);
	#endif
}

bool BioSemiClient::openDevice() {
	char bytes[64];
	intptr_t pointer;
	
	if (!driverOk) return false;
	if (deviceOpen) closeDevice();
	
	deviceHandle = lv_open_driver_async();
	if (deviceHandle == NULL) {
		fprintf(stderr, "Cannot open the USB driver!\n");
		return false;
	}
	
	//printf("Ok so far\n");
	
	memset(bytes, 0, sizeof(bytes));
	
	//printf("Before lv_usb_write\n");
	
	if (!lv_usb_write(deviceHandle, bytes)) {
		fprintf(stderr, "Cannot write to/initialize USB device!\n");
		lv_close_driver_async(deviceHandle);
		return false;
	}
	
	//printf("Before lv_read_multiple_sweeps\n");
	
	lv_read_multiple_sweeps(deviceHandle, (char *) ringBuffer, BUFFER_SIZE);
	
	memset(bytes, 0, sizeof(bytes));
	bytes[0] = (char) 0xFF;
   
    //usleep(500);
	
	if (!lv_usb_write(deviceHandle, bytes)) {
		fprintf(stderr, "Cannot enable handshake!\n");
		lv_close_driver_async(deviceHandle);
		return false;
	}
	
	timeHandshake = getCurrentTime();
	while (1) {
		
		if (!lv_read_pointer(deviceHandle, &pointer)) {
			fprintf(stderr, "Cannot read ring buffer pointer!\n");
			closeDevice();
			return false;
		}
		
		if (pointer > 2) break;
		
		msleep(10);
		if (getCurrentTime() - timeHandshake > DEVICE_TIMEOUT) {
			fprintf(stderr, "Timeout for reading initial samples - power down?\n");
			closeDevice();
			return false;
		}
	}
	
	int32_t sync = ringBuffer[0];
	if (sync != SYNC_BV) {
		fprintf(stderr, "Device is not in sync\n");
		closeDevice();
		return false;
	}	
	
	int32_t status = ringBuffer[1];
	
	deviceIsMk2 = (status & MK2_BV) ? true : false;
	speedMode = 0;
	if (status & SPEED_BIT3) speedMode+=8;
	if (status & SPEED_BIT2) speedMode+=4;
	if (status & SPEED_BIT1) speedMode+=2;
	if (status & SPEED_BIT0) speedMode+=1;
	cmsInRange = (status & CMS_RANGE_BV) ? true : false;
	
	if (speedMode < 0 || speedMode > 8) {
		fprintf(stderr, "Unsupported speed mode discovered\n");
		closeDevice();
		return false;
	}
	if (deviceIsMk2) {
		numChannels = sNumChannelsMk2[speedMode];
	} else {
		numChannels = sNumChannelsMk1[speedMode];
	}
	numChanAIB = (speedMode == 8) ? 32 : 0;
	stride = numChannels + numChanAIB + 2;
	
	switch(speedMode) {
		case 0:
        case 4:
        case 8: 
            sampleFreq = 2048;
            break;
        case 1:
        case 5:
            sampleFreq = 4096;
            break;
        case 2:
        case 6:
            sampleFreq = 8192;
            break;
        case 3:
        case 7:
            sampleFreq = 16384;
            break;
	}

	lastPointerRead  = 0;
	lastStartPtr  = 0;
	deviceOpen = true;
	/*
	for (int i=0;i<pointer;i+=4) {
		int val = *((int *) (ringBuffer+i));
		if (val==SYNC_BV) printf("Sync at %i\n", i);
	}
	*/
	return true;
}


void BioSemiClient::closeDevice() {
	char bytes[64];
	
	if (!deviceOpen) return;
	
	memset(bytes, 0, sizeof(bytes));
	if (!lv_usb_write(deviceHandle, bytes)) {
		fprintf(stderr, "Cannot disable handshake!\n");
	}
	if (!lv_close_driver_async(deviceHandle)) {
		fprintf(stderr, "Error when closing device!\n");
		// well, what can we do?
	}
	deviceOpen = false;
}

int BioSemiClient::getCurPointer() {
	intptr_t pointer;
	
	if (!deviceOpen) return -1;

	if (!lv_read_pointer(deviceHandle, &pointer)) {
		return -1;
	}
	
	return pointer/sizeof(int32_t); 
}

bool BioSemiClient::checkNewBlock(BioSemiBlock &block) {
	int pointer, numSamples, nextStartPtr;
	
	pointer = getCurPointer();
	
	// If the pointer hasn't changed, there is no new data
	if (pointer==lastPointerRead) return false;	
	
	// Otherwise, check if the pointer has wrapped around,
	// and calculate how many new (complete) samples there are.
	// Apparently, the pointers we retrieve do not have to fall
	// on whole sample boundaries, so we maintain a lastStartPtr
	// as well.
	//printf("now/last: %6i / %6i -- ", pointer, lastPointerRead);
	if (pointer > lastStartPtr) {
		// no wrap around happened
		numSamples = (pointer - lastStartPtr) / stride;
	} else {
		// has wrapped at BUFFER_SIZE bytes
		numSamples = (pointer + BUFFER_LEN - lastStartPtr) / stride;
	}
	
	// new bytes, but not a complete sample - this is unlikely...
	if (numSamples == 0) return false;
	
	// Calculate the next start ptr and wrap around, if needed
	nextStartPtr = lastStartPtr + numSamples * stride;
	if (nextStartPtr >= BUFFER_LEN) nextStartPtr -= BUFFER_LEN;
	
	block.startIndex = lastStartPtr;
	block.numSamples = numSamples;
	block.numInSync  = 0;
	block.stride     = stride;
	block.batteryLow = false;
	
	int idx = lastStartPtr;
	for (int i=0;i<numSamples;i++) {
		int sync   = getValue(idx);
		int status = getValue(idx+1);
		
		if (sync==SYNC_BV)     block.numInSync++;
		if (status&BATTERY_BV) block.batteryLow=true;
		
		idx+=stride;
	}
	
	lastPointerRead = pointer;
	lastStartPtr = nextStartPtr;
	return true;
}

