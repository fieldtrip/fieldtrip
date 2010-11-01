
/*++

Copyright (c) 2000 TMS International 

Module Name:
    RTDevice.cpp

Abstract:  
	RTDevice class implementation
	
Version:  1.3

Author:
    Rob Hemstede

Environment:
    User mode

Revision History:
1 - Version 1.00:	first release
2 - Version 1.01:	Support for multiple devices 
3 - Version 1.02:	Support for Calibration and Impedance measurement 
4 - Version 1.03:	Support for device specific features 

--*/

#include<Windows.h>
#include <wchar.h>

#include "RTDevice.h"


BOOLEAN RTDevice::Init()
{	TCHAR Path[	MAX_PATH ]; 

	//If a device is still in use close it
	if(Handle != INVALID_HANDLE_VALUE ) 
	{	//TO DO Maybe stop the device ? 
		fpClose(Handle);
		Handle = INVALID_HANDLE_VALUE;
	}
	
	//Close if there is an open lib 
	if(LibHandle != INVALID_HANDLE_VALUE)
	{	FreeLibrary((struct HINSTANCE__ *)LibHandle); 
		LibHandle = NULL;	
	}

	if( dp != NULL )
	{	Free(dp); 
		dp = NULL; 
	}

	//Open libratry
	GetSystemDirectory(Path, sizeof(Path) / sizeof(TCHAR) );
	lstrcat(Path,RTLOADER);
	LibHandle = LoadLibrary(Path); 
	
	//if the can not be opend return FALSE;
	if( LibHandle == NULL ) return FALSE; 
	
	//Get pointers to the funtions in the DLL 
	fpOpen				= (POPEN)			GetProcAddress(LibHandle,"Open");
	fpClose				= (PCLOSE)			GetProcAddress(LibHandle,"Close");
	fpStart				= (PSTART)			GetProcAddress(LibHandle,"Start");
	fpStop				= (PSTOP)			GetProcAddress(LibHandle,"Stop");
	fpReset				= (PRESETDEVICE)	GetProcAddress(LibHandle,"ResetDevice");
	fpGetDeviceState	= (PGETDEVICESTATE)	GetProcAddress(LibHandle,"GetDeviceState");
	fpSetSignalBuffer	= (PSETSIGNALBUFFER)GetProcAddress(LibHandle,"SetSignalBuffer");
	fpGetBufferInfo		= (PGETBUFFERINFO)	GetProcAddress(LibHandle,"GetBufferInfo");
	fpGetSamples		= (PGETSAMPLES)		GetProcAddress(LibHandle,"GetSamples");
	fpGetSlaveHandle	= (PGETSLAVEHANDLE) GetProcAddress(LibHandle,"GetSlaveHandle");
	fpAddSlave			= (PADDSLAVE)		GetProcAddress(LibHandle,"AddSlave");
	fpDeviceFeature		= (PDEVICEFEATURE) 	GetProcAddress(LibHandle,"DeviceFeature");
	fpGetSignalFormat	= (PGETSIGNALFORMAT)GetProcAddress(LibHandle,"GetSignalFormat"); 
	fpGetInstanceId		= (PGETINSTANCEID)	GetProcAddress(LibHandle, "GetInstanceId" ); 
	fpOpenRegKey		= (POPENREGKEY)		GetProcAddress(LibHandle, "OpenRegKey" );
	fpFree				= (PFREE)			GetProcAddress(LibHandle, "Free" ); 

	// All other DLL exported function are for obsolete but
	// available for compatibility reasons

	DeviceRegKey = NULL; 

	return TRUE;
}

BOOLEAN RTDevice::Free( void *Memory )
{	if( fpFree ) return fpFree( Memory );
	//If this is not supported we can do it manually
	if( LocalFree( Memory ) == NULL ) return TRUE; 
	return FALSE;
}


RTDevice::RTDevice()
{ 
	LibHandle = NULL ; 
	Handle=INVALID_HANDLE_VALUE; 
	NrOfDevices=0;
	DidAttach = FALSE; 
	InitOk = FALSE;
	dp = NULL;
	DeviceRegKey = NULL; 

	if( !Init() )
		InitOk = FALSE; 	
	else
		InitOk = TRUE; 
	return;
}


RTDevice::RTDevice(ULONG DeviceIndex)
{	LibHandle = NULL ; 
	Handle=INVALID_HANDLE_VALUE; 
	NrOfDevices=0;
	DidAttach = FALSE; 
	InitOk = FALSE;
	dp = NULL;
	DeviceRegKey = NULL; 

	// Initialize selected device 
	if(!Init()) return; 

	if( OpenHandle(DeviceIndex) )
	{	InitOk = TRUE; 
		return ;
	
	}
	return;
}

RTDevice::RTDevice(HANDLE NewHandle)
{	
	LibHandle = NULL; 
	Handle=INVALID_HANDLE_VALUE; 
	NrOfDevices=0;
	DidAttach = TRUE; 
	InitOk = FALSE; 
	dp = NULL;
	DeviceRegKey = NULL; 

	if(!Init() ) return; 
	
	Handle = NewHandle; 
	InitOk = TRUE; 
	return; 
} 


BOOLEAN RTDevice::OpenHandle(ULONG DeviceIndex)
{
	//If we have a handle then 
	if( Handle != INVALID_HANDLE_VALUE ) 
		return TRUE; 

	if( fpGetInstanceId != NULL ) 
	{	//This is the total number of devices connected to the system
		
		//the unique device name 
		dp =  fpGetInstanceId(DeviceIndex, TRUE , &NrOfDevices); 

		if(( dp != NULL ) && ( fpOpenRegKey != NULL )) 
			DeviceRegKey = fpOpenRegKey( dp ); 

		//Open a handle to the requested device 
		Handle = fpOpen(dp);
		//Check if the device is opened 
		if( Handle == INVALID_HANDLE_VALUE ) return FALSE; 
		return TRUE; 
	}
	
	return FALSE; 
}
	

//Destructor 
RTDevice::~RTDevice()
{
	if( dp != NULL ) Free(dp);

	// if we are holding a handle, close it.
	if( !DidAttach ) Close(Handle); 
	
	//if we are holding a Handle, close it.
	if( LibHandle != NULL ) FreeLibrary(LibHandle); 

	if( DeviceRegKey ) 
		RegCloseKey( DeviceRegKey ); 

	return;
}


void RTDevice::Close(HANDLE DevHandle )
{	if( DevHandle == INVALID_HANDLE_VALUE ) return; 
	if( fpClose == NULL ) return; 
	fpClose( DevHandle ); 
	return; 
}


/* Start
	Start sampling
	
	Output:
		if successful:	TRUE
		if unsuccessful:	FALSE

*/

BOOLEAN RTDevice::Start()
{	if(fpStart == NULL) return FALSE;
	return  fpStart(Handle); 
}

/* Stop
	Stop sampling
	
	Output:
		if successful:	TRUE
		if unsuccessful:	FALSE

*/
BOOLEAN RTDevice::Stop()
{	if(fpStop == NULL) return FALSE;
	return fpStop(Handle); 
}

/* Reset
	Reset device 
	
	Output:
		if successful:	TRUE
		if unsuccessful:	FALSE

*/
BOOLEAN RTDevice::Reset()
{	if(fpReset == NULL) return FALSE;
	return fpReset(Handle); 
}


/*
	GetDeviceState
		returns either:	
			DEVICE_STATE_STOP of device is not sampling 
			DEVICE_STATE_START of device is sampling 

			or DEVICE_STATE_UNKNOWN if unsuccessful. 


*/
ULONG RTDevice::GetDeviceState()
{	if(fpGetDeviceState == NULL) return DEVICE_STATE_UNKNOWN; 
	return fpGetDeviceState(Handle);
}



/*	
	SetSignalBuffer
		Initializes the signal buffer and sets the sample rate 
	Input
		SampleRate:	the requested sample rate in milli Hz (mHz) 
		BufferSize: the requested buffer size in samples 

	Output:
		SampleRate: if the device cannot satisfy the requested sample rate 
					it will return the closed math
		BufferSize:	if the device cannot satisfy the requested sample rate 
					it will return the closed math
	
	returns:
		TRUE if successfull
		FALSE if unsuccessfull


*/
BOOLEAN RTDevice::SetSignalBuffer(IN OUT PULONG SampleRate,IN OUT PULONG BufferSize)
{	if(fpSetSignalBuffer == NULL) return FALSE;
	return fpSetSignalBuffer(Handle,SampleRate,BufferSize);
}


/*
	GetBufferInfo
		Get information about the status of the signal buffer

	input:
	
	return:
		returns TRUE if successful / FALSE if unsuccessful 

*/

BOOLEAN RTDevice::GetBufferInfo(OUT PULONG Overflow,OUT PULONG PercentFull)
{	if(fpGetBufferInfo == NULL) return FALSE;
	return fpGetBufferInfo(Handle,Overflow,PercentFull);
}



/*
	GetSamples:
		Returns the samples in the SamplesBuffer if there are any available

	Input:
		OUT SampleBuffer. Pointer to the buffer were the samples will be stored
		IN Size. Size of the buffer in bytes. if the size permits the function will return 
		multiple samples 

	Output:
		Number of bytes returned in the buffer


*/

ULONG RTDevice::GetSamples(OUT PULONG SampleBuffer,IN ULONG Size)
{	if(fpGetSamples==NULL) return 0; 
	return fpGetSamples(Handle,SampleBuffer,Size);
}


/* 
	Get slave handle 

  	OUTPUT 
		Handle to the slave device, This handle is needed by the add slave function

*/


HANDLE RTDevice::GetSlaveHandle()
{	if(fpGetSlaveHandle==NULL) return NULL; 
	return fpGetSlaveHandle(Handle);
}
		


/*	
	Add Slave 

	INPUT
		Handle to the slave device obtained by the GetSlaveHandle function 

	OUTPUT
		returns TRUE if successful 
		
*/

BOOLEAN RTDevice::AddSlave(HANDLE SlaveHandle)
{	if(fpAddSlave==NULL) return FALSE; 
	return fpAddSlave(Handle,SlaveHandle);
}




/*
	Device Feature 
	  
	General function for setting device specific features 
	revert to your device specific software manual for 
	information abuot the parameter

*/


BOOLEAN RTDevice::DeviceFeature( PVOID InputBuffer , ULONG InLength , PVOID OutputBuffer , ULONG OutLength )
{	if( fpDeviceFeature == NULL ) return FALSE; 
	return fpDeviceFeature( Handle , InputBuffer , InLength , OutputBuffer , OutLength );
}



/*********************************************************************/

PSIGNAL_FORMAT RTDevice::GetSignalFormat( PSIGNAL_FORMAT InfoBuffer )
{	if( fpGetSignalFormat == NULL ) return NULL; 	
	return fpGetSignalFormat( Handle , InfoBuffer ); 
}


