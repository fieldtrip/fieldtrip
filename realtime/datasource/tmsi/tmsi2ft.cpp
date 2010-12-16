// tmsidriver.cpp : Defines the entry point for the console application.
//




/*++
 
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


using namespace std;
  
#define MAX_DEVICE			1	//Max number of devices supported by this demo
#define USE_MASTER_SLAVE	FALSE




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


ULONG getTotalNumChannels(RTDeviceEx *Master) {
	ULONG numChan;
	
	SIGNAL_FORMAT sf;
    PSIGNAL_FORMAT psf;
	
	memset(&sf, 0, sizeof(SIGNAL_FORMAT));
	sf.Elements = 1; 
	sf.Size = sizeof(SIGNAL_FORMAT); 

	psf = (SIGNAL_FORMAT *) LocalAlloc( LMEM_FIXED | LMEM_ZEROINIT , sf.Size * sf.Elements); 
	if( psf == NULL ) {
		fprintf(stderr, "Cannot allocate memory for SignalFormat structure\n");
		exit(1);
	}
	psf[0].Size = sf.Size; 
	psf[0].Elements = sf.Elements;
	Master->GetSignalFormat(psf);

	numChan = psf[0].Elements;
	LocalFree(psf); 
	
	return numChan;
}


int main(int argc, char* argv[]) {	
	StringServer ctrlServ;
	OnlineDataManager<int32_t, float> *ODM;
	
	RTDeviceEx *Master;
	ULONG SampleRate = MAX_SAMPLE_RATE;
	ULONG BufferSize = MAX_BUFFER_SIZE;
	
	ULONG PercentFull,Overflow;
	ULONG BytesPerSample=0;
	ULONG BytesReturned;
	ULONG numHwChans;

	// Buffer for storing the samples; 
	ULONG SignalBuffer[1000];

	Master=InitDevice(SampleRate);
	if (Master == 0) {
		fprintf(stderr, "Cannot initialise device\n");
		return 0;
	}
	Master->SetSignalBuffer(&SampleRate,&BufferSize);

	numHwChans = getTotalNumChannels(Master);
	BytesPerSample = 4*numHwChans;

	if( BytesPerSample == 0 ) 
	{	printf( "\nDevice returns no samples" ); 
		_getch();
		return 0; 
	}

	wprintf(L"\nMaximum sample rate = %d Hz",SampleRate  / 1000 );
	wprintf(L"\nMaximum Buffer size = %d Samples",BufferSize);
	        
	//Set sample rate 
	//SampleRate = 500000 ;     
	//Set buffer size;
	BufferSize = 1000; 

	Master->SetSignalBuffer(&SampleRate ,&BufferSize);

	wprintf(L"\nSelected sample rate = %d Hz",SampleRate  / 1000 );
	wprintf(L"\nSelected Buffer size = %d Samples",BufferSize);

	/* these represent the acquisition system properties */
	float fSample        = SampleRate/1000.0;
	int nBufferSamp	   = 0;
	int nTotalSamp	   = 0;
	
	ctrlServ.startListening(8000);
	ODM = new OnlineDataManager<int32_t, float>(0, numHwChans, fSample, GDF_INT32, DATATYPE_FLOAT32);

	if (!Master->Start()) {
		fprintf(stderr, "Unable to start the Device\n");
		return 1;
	}

	if (!ODM->useOwnServer(1972)) {
		fprintf(stderr, "Could not spawn buffer server.\n");
		return 0;
	}
	if (ODM->configureFromFile("config.txt") != 0) {
		fprintf(stderr, "Configuration file is invalid\n");
		return 0;
	}	
	ODM->enableStreaming();
	
	wprintf(L"\nPress any key to quit \n");

	// Stop the program if we hit any key
	while(!_kbhit()) {
		ctrlServ.checkRequests(*ODM);
		//Get Signal buffer information
		Master->GetBufferInfo(&Overflow,&PercentFull);
			
		if (PercentFull > 0) {
			// if there is data available get samples from the device
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
				ODM->handleBlock();
			} else {
				//allow other applications some extra process time 
				Sleep(10);
			}
		}//end PercentFull	
	}//endwhile

	//Stop the device

	if (Master->Stop()) {	
		wprintf(L"\nDevice stop\n");
	} else {	
		wprintf(L"\nUnable to stop the device\n");
	}
		
	//Unchain any slave devices 
	Master->GetSlaveHandle();
	delete ODM;

	return 0;
}
