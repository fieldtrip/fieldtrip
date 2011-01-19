/** TMSI2FT acquisition tool, based on OnlineDataManager and code from
    Bart Nienhuis' tmsidriver, as well as TMSi sample code.
	Hacked together by S. Klanke, 2010
*/

/*++

Parts of the code (TMSi examples)
 
Copyright (c) 2000 TMS International

Author:
    Rob Hemstede

Environment:
    User mode

Revision History: 
1 - Version 1,00:	first release
2 - Version 1,01:	Support for multiple devices 
3 - Version 1.02:	Support for Calibration and Impedance measurement 
4 - Version 1.03:	Support for device specific features 
5 - Version 1.04;   Using Device Instance Id

--*/

/*

  General:
		
		Whenever I refer to a 'sample' in this software I mean 
		a data structure representing a conversion result for all 
		channels/signal for one single sample moment.  
	
*/
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <Windows.h>
#include <stdio.h>
#include <wchar.h>
#include <conio.h>
#include <tchar.h> 
#include <iostream>
#include <fstream>
#include <string>


#include "RtDevice.h"
#include "Feature.h" 
#include <math.h>

//#include "buffer.h" 
//#include "pthread.h" 

#include <OnlineDataManager.h>
#include <StringServer.h>
#include <ConsoleInput.h>


using namespace std;
  
#define MAX_DEVICE			1	//Max number of devices supported by this demo
#define USE_MASTER_SLAVE	FALSE

#define MY_BUFFER_SIZE      10000


// The following function is from Bart Niehuises tmsidriver application,
// or from the TMSi example code.
ULONG UseMasterSlave( RTDeviceEx **Devices , ULONG Max )
{	ULONG NrOfDevices; 
	ULONG x;

	for(x = 0 ; x < Max ;x++ )
	{	Devices[x] = new RTDeviceEx(x);
		if( Devices[x] == NULL ) break; 
		
		if( !Devices[x]->InitOk) 
		{	delete Devices[x]; 
			Devices[x] = NULL; 
			break; 
		}		
	}

	NrOfDevices = x; 

	for( x = 0 ; x < NrOfDevices ; x++ ) 
	{	TCHAR DeviceName[40] = _T("Unknown Device"); 
		ULONG SerialNumber = 0; 
		ULONG Size;
	
		Size = sizeof( SerialNumber );
		RegQueryValueEx( Devices[x]->DeviceRegKey , _T("DeviceSerialNumber"), NULL , NULL , (PBYTE)&SerialNumber , &Size  );

		Size = sizeof( DeviceName );
		RegQueryValueEx( Devices[x]->DeviceRegKey , _T("DeviceDescription"), NULL , NULL , (PBYTE)&DeviceName[0] , &Size  ); 
		
		_tprintf( _T("%lud . %s %lud\n"), x + 1, DeviceName , SerialNumber	); 

		
		if( x!= 0 ) 
		{	HANDLE SlaveHandle = Devices[x]->GetSlaveHandle(); 
			if(SlaveHandle == 0 )
				_tprintf( _T("Unable to get a handle from device %lud\n"), x + 1 );
			break;

			Devices[0]->AddSlave(SlaveHandle);
		}
	}
	return NrOfDevices;
}

// The following function is from Bart Niehuises tmsidriver application,
// or from the TMSi example code.
RTDeviceEx *SelectDevice( IN BOOLEAN Present )
{	ULONG Count = 0; 
	ULONG Max = 0; 
	RTDeviceEx *Device; 

	Device = new RTDeviceEx; 

	if( Device == NULL )
		return NULL;

	if( !Device->InitOk  )
	{	delete Device;
		return NULL; 
	}
	
	PSP_DEVICE_PATH Id; 

	while(1)
	{	TCHAR DeviceName[40] = _T("Unknown Device"); 
		ULONG SerialNumber = 0; 

	
		HKEY hKey; 
		
		Id = Device->GetInstanceId( Count++ , Present , &Max ); 
		if( !Id ) break; 

		hKey = Device->OpenRegistryKey( Id ); 

		if( hKey != INVALID_HANDLE_VALUE ) 
		{	ULONG Size; 
	
			Size = sizeof( SerialNumber );
			RegQueryValueEx( hKey , _T("DeviceSerialNumber"), NULL , NULL , (PBYTE)&SerialNumber , &Size  );

			Size = sizeof( DeviceName );
			if( RegQueryValueEx( hKey , _T("DeviceDescription"), NULL , NULL , (PBYTE)&DeviceName[0] , &Size  ) 
				== ERROR_SUCCESS )
			{
				_tprintf( "%lud . %s %lud\n" , Count , DeviceName , SerialNumber	); 
			}
				
			RegCloseKey( hKey ); 
		}

		Device->Free( Id ); 		
	}

	if( Max == 0 ) 
	{	printf("There are no device connected to the PC\n");
		return NULL; 
	}

	if( Max == 1 )
	{	Id = Device->GetInstanceId( 0 , Present );
	
	}
	else
	{	printf("Please select device ...\n\n"); 
		while( _kbhit() ){}
		while( !_kbhit() ){}
		int key = _getch() - '0';
		Id = Device->GetInstanceId( key - 1 , Present );
	}

	if( !Device->Open( Id ) ) 
	{	Device->Free( Id );	
		delete Device; 
		return NULL; 
	}

	return Device; 
}

// The following function is from Bart Niehuises tmsidriver application,
// or from the TMSi example code.
// Get handle to my device and init the device 
RTDeviceEx *InitDevice( ULONG SampRate) {
	ULONG Index;	

	RTDeviceEx *Device[MAX_DEVICE];
	for(Index=0;Index < MAX_DEVICE;Index++)
		Device[Index] = NULL;
	
	RTDeviceEx *MasterL;
	
	if( USE_MASTER_SLAVE )
		UseMasterSlave( Device , MAX_DEVICE ); 
	else 
		Device[0] = SelectDevice( TRUE );  

	MasterL = Device[0]; 
		
	if( MasterL == NULL ) 
	{	_getch();
		return 0; 
	
	}

	MasterL->Reset();		

	return MasterL;
}

/** Retrieve number of channels and the index of the trigger channel.
    If there is more than one trigger channnel (digital input), only
	the last one is used. TODO: maybe change this in a later version.
*/
int getTotalNumberOfChannels(RTDeviceEx *Master, int& triggerChannel) {
	int numChan;
    PSIGNAL_FORMAT psf;

	psf = Master->GetSignalFormat(NULL);
	if (psf == NULL) return 0;
	
	// printf("%i x %i\n", (int) psf[0].Size, (int) psf[0].Elements);
	numChan = psf[0].Elements;
	triggerChannel = -1;
	
	for (int i=0;i<numChan;i++) {
		printf("Channel %i: %i %i ", i+1, (int) psf[i].Type, (int) psf[i].SubType);
		wprintf(psf[i].Name); // .Name field is of WCHAR type (unicode)
		printf("\n");
		
		// the documentation gives 0x13 for the type, but at least
		// for the Porti we need a "4".
		if (psf[i].Type == 4) triggerChannel = i;
	}
	// do we need to free this? or did we get static memory?
	// LocalFree(psf);   
	return numChan;
}

int main(int argc, char* argv[]) {	
	StringServer ctrlServ;
	ConsoleInput conIn;
	OnlineDataManager<int32_t, float> *ODM;
	SignalConfiguration sigConf;
	
	char hostname[256];
	int port;
	int ctrlPort;
	
	int triggerChannel = -1;
	int triggerStatus  = 0;
	
	RTDeviceEx *Master;
	ULONG SampleRate = MAX_SAMPLE_RATE;
	ULONG BufferSize = MAX_BUFFER_SIZE;
	
	ULONG PercentFull,Overflow;
	ULONG BytesPerSample=0;
	ULONG BytesReturned;
	ULONG numHwChans;

	// Buffer for storing the samples; 
	ULONG SignalBuffer[MY_BUFFER_SIZE];
	
	if (argc<2) {
		fprintf(stderr, "Usage: tmsi2ft configfile [hostname=localhost [port=1972 [ctrlPort=8000]]]\n");
		fprintf(stderr, "\nUse a minus (-) for the hostname parameter to spawn a local buffer server.\n");
		return 1;
	}
	
	if (argc>2) {
		strncpy(hostname, argv[2], sizeof(hostname));
	} else {
		strcpy(hostname, "localhost");
	}
	
	if (argc>3) {
		port = atoi(argv[3]);
	} else {
		port = 1972;
	}
	
	if (argc>4) {
		ctrlPort = atoi(argv[4]);
	} else {
		ctrlPort = 8000;
	}
	
	if (sigConf.parseFile(argv[1]) != 0) {
		fprintf(stderr, "Configuration file %s is invalid\n", argv[1]);
		return 1;
	}	
	
	if (!ctrlServ.startListening(ctrlPort)) {
		fprintf(stderr, "Cannot listen on port %d for configuration commands\n", ctrlPort);
		return 1;
	}

	Master=InitDevice(SampleRate);
	if (Master == 0) {
		fprintf(stderr, "Cannot initialise device\n");
		return 1;
	}
	Master->SetSignalBuffer(&SampleRate,&BufferSize);

	numHwChans = getTotalNumberOfChannels(Master, triggerChannel);
	BytesPerSample = 4*numHwChans;

	if( BytesPerSample == 0 ) {	
		fprintf(stderr, "Device returns no samples\n"); 
		return 1; 
	}

	printf("Maximum sample rate   = %.3f Hz\n",(float) SampleRate / 1000.0);
	printf("Maximum Buffer size   = %d Samples\n",(int) BufferSize);
	printf("Number of HW channels = %d\n", (int) numHwChans);        

	BufferSize = MY_BUFFER_SIZE; 
	if (sigConf.getSampleRate() != 0.0) {
		// override SampleRate if configured this way
		SampleRate = (int) (sigConf.getSampleRate() * 1000.0);
	}
	
	Master->SetSignalBuffer(&SampleRate ,&BufferSize);

	printf("Selected sample rate = %.3f Hz\n", (float) SampleRate / 1000.0);
	printf("Selected Buffer size = %d Samples\n", (int) BufferSize);

	/* these represent the acquisition system properties */
	float fSample      = SampleRate/1000.0;
	int nBufferSamp	   = 0;
	int nTotalSamp	   = 0;
	
	ODM = new OnlineDataManager<int32_t, float>(0, numHwChans, fSample);

	if (!Master->Start()) {
		fprintf(stderr, "Unable to start the Device\n");
		return 1;
	}
	if (!strcmp(hostname, "-")) {
		if (!ODM->useOwnServer(port)) {
			fprintf(stderr, "Could not spawn buffer server on port %d.\n",port);
			goto cleanup;
		}
	} else {
		if (!ODM->connectToServer(hostname, port)) {
			fprintf(stderr, "Could not connect to buffer server at %s:%d.\n",hostname, port);
			goto cleanup;
		}
	}
	
	if (!ODM->setSignalConfiguration(sigConf)) {
		fprintf(stderr, "Could not set OnlineDataManager configuration. Did you specify more channels than the HW provides?\n");
		goto cleanup;
	}
		
	ODM->enableStreaming();
	
	printf("\nPress [Escape] to quit...\n");
	
	while (1) {
		if (conIn.checkKey()) {	
			int c = conIn.getKey();
			if (c==27) break; // quit
		}
		// Process any incoming request on the control port
		ctrlServ.checkRequests(*ODM);
		
		//Get Signal buffer information
		Master->GetBufferInfo(&Overflow,&PercentFull);
			
		if (PercentFull > 0) {
			// If there is data available, get samples from the device
			// GetSamples returns the number of bytes written in the signal buffer
			// This will always be a multiple op BytesPerSample. 
			
			// Divide the result by BytesPerSamples to get the number of samples returned
			BytesReturned = Master->GetSamples((PULONG)SignalBuffer,sizeof(SignalBuffer));
				
			if (BytesReturned != 0) {
				// wprintf(L"\rSampleCounter = %8d, ,Sampwritten=%8d,Samp[%d]=%d, %d ,Buffer  = %d, Overflow = %d      " , Total /BytesPerSample ,sample, ShowChannel,SignalBuffer[ShowChannel], SignalBuffer[ShowChannel + 1],PercentFull,Overflow);
				nBufferSamp=BytesReturned/BytesPerSample;
				nTotalSamp+=nBufferSamp;
				
				int32_t *data = ODM->provideBlock(nBufferSamp);
				if (data==0) {
					fprintf(stderr, "Out of memory\n");
					break;
				}
				memcpy(data, SignalBuffer, nBufferSamp * BytesPerSample);
				
				// TODO: allow for multiple trigger channels
				if (triggerChannel >= 0) {
					for (int j=0;j<nBufferSamp;j++) {
						int trigVal = SignalBuffer[triggerChannel + j*numHwChans];
						
						if (trigVal != triggerStatus && trigVal != 0) {
							// TODO: maybe use the channel label as the event type instead of "Digi"
							ODM->getEventList().add(j, "Digi", trigVal);
						}
						triggerStatus = trigVal;
					}
				}
				
				ODM->handleBlock();
			}
		} else {
			Sleep(1);
		}
	}

cleanup:
	//Stop the device
	if (Master->Stop()) {	
		printf("Device stopped.\n");
	} else {	
		fprintf(stderr, "Unable to stop the device!\n");
	}
		
	//Unchain any slave devices 
	Master->GetSlaveHandle();
	delete ODM;

	return 0;
}
