/*++

Copyright (c) 2000 TMS International

Module Name:
    RTDevice.h 

Abstract:  
	Real time device, header file 
	
Version:  1.03 

Author:
    Rob Hemstede

Environment:
    User mode


--*/
#ifndef __RTDEVICE_H__
#define __RTDEVICE_H__

//Name of the DLL needed to 
#define RTLOADER "\\RTINST.Dll"

//Header file is used by appliaction and by the WDM Driver
#include "SADIO.H"


#define MAX_SAMPLE_RATE 0xFFFFFFFF
#define MAX_BUFFER_SIZE 0xFFFFFFFF

//classes 
class RTDevice{
		//Pointers to functions in the RTInst DLL
		POPEN fpOpen;
		PCLOSE fpClose; 
		PSTART fpStart;
		PSTOP fpStop;	
		PRESETDEVICE fpReset;
		PGETDEVICESTATE fpGetDeviceState;
		PSETSIGNALBUFFER fpSetSignalBuffer;
		PGETBUFFERINFO fpGetBufferInfo;
		PGETSAMPLES	fpGetSamples;
		PGETSLAVEHANDLE fpGetSlaveHandle;
		PADDSLAVE fpAddSlave; 
		PDEVICEFEATURE fpDeviceFeature; 
		PGETSIGNALFORMAT fpGetSignalFormat; 
		PGETINSTANCEID fpGetInstanceId; 
		POPENREGKEY fpOpenRegKey; 
		PFREE fpFree;
		
		HANDLE Handle;			//Device Handle
		HINSTANCE LibHandle;	//Liberary Handle
		BOOLEAN DidAttach; 
		BOOLEAN Init();	

	public:
	
		HKEY DeviceRegKey; 

		BOOLEAN OpenHandle(ULONG DeviceIndex);
		
		BOOLEAN Open( PSP_DEVICE_PATH InstanceId )
		{	
			if( (fpOpen == NULL) )
				return FALSE; 

			if( Handle != INVALID_HANDLE_VALUE ) 
			{	Close( Handle ); 
				Handle = INVALID_HANDLE_VALUE; 
			}			

			Handle = fpOpen( InstanceId  );
			
			if( Handle != INVALID_HANDLE_VALUE ) 
				return TRUE; 

			return FALSE; 
		}

		HKEY OpenRegistryKey( PSP_DEVICE_PATH InstanceId )
		{	if( fpOpenRegKey == NULL ) 
				return (HKEY)INVALID_HANDLE_VALUE; 
			
			return fpOpenRegKey( InstanceId );
		}


		PSP_DEVICE_PATH GetInstanceId( LONG DeviceIndex , BOOLEAN Present, OUT ULONG  *MaxDevices = NULL )
 		{	ULONG lMaxDevices = 0;
			PSP_DEVICE_PATH DeviceInstance = NULL; 
			
			if( fpGetInstanceId != NULL )
				DeviceInstance = fpGetInstanceId( DeviceIndex , Present , &lMaxDevices );
						
			if( MaxDevices ) 
				*MaxDevices = lMaxDevices; 

			return DeviceInstance; 
 		}

		PSP_DEVICE_PATH dp;		// Device path 

		BOOLEAN InitOk;			// True if device has initialized ok
		ULONG NrOfDevices;		// Number of devices on this PC

		RTDevice(); 
		RTDevice(ULONG Index);	
		RTDevice(HANDLE Attach); 
		~RTDevice();			//Destructor 

		//Functions for getting information about the device
		PWCHAR GetDescription(PWCHAR Destination,ULONG Size);
		PWCHAR GetManufacturer(PWCHAR Destination,ULONG Size);
		ULONG GetId();

		//Funtion for getting information about the signals of a device
		//SIGNAL_INFO defined in SADIO.H
		PSIGNAL_INFO GetSignalInfo(PSIGNAL_INFO pSignalInfo,PULONG NrOfSignals);
			
		BOOLEAN Reset();		//Reset device 
		BOOLEAN Start();		//Start sampling
		BOOLEAN Stop();			//Stop sampling
		//Get state ( Sampling or not sampling of a devive )
		//States defined in SADIO.H 
		ULONG GetDeviceState(); 
			
		//Initailize the signal buffer
		BOOLEAN SetSignalBuffer(IN OUT PULONG SampleRate,IN OUT PULONG BufferSize);
		//Get info about the status of the signal buffer
		BOOLEAN GetBufferInfo(OUT PULONG OverflowCounter,OUT PULONG PercentFull); 

		//Get Samples
		ULONG GetSamples(PULONG SampleBuffer,ULONG Size); 

		//Get slave handle 
		HANDLE GetSlaveHandle();
		
		//Add Slave 
		BOOLEAN AddSlave(HANDLE SlaveHandle); 

		
		BOOLEAN SetMode(ULONG Mode,ULONG Info); 
		BOOLEAN GetMode(ULONG *Mode,ULONG *Info); 

		//Device Specific feature 
		BOOLEAN DeviceFeature( PVOID InputBuffer , ULONG InLength , PVOID OutputBuffer , ULONG OutLength );
		
		//Used for getting device information from the registry 
		BOOLEAN GetDeviceKey( PWCHAR Name , PKEY_VALUE_INFORMATION Key );

		void Close(HANDLE DevHandle );  

		PSIGNAL_FORMAT GetSignalFormat( PSIGNAL_FORMAT InfoBuffer );

		BOOLEAN Free( void *Memory );
	
};

#endif  // __RTDEVICE_H__



