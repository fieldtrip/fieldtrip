/*++

Copyright (c) 2001 TMS International

Module Name:
    Feature.h

Abstract:  
	Implementaion example for device specific features 
	
Version:  1.3

Author:
    Rob Hemstede

Environment:
    User mode

--*/

#ifndef __FEATURE_H__
#define __FEATURE_H__

#include"RTDevice.h" 

// added by Stefan Klanke for compilation using MinGW
#ifndef PSYSTEMTIME
#define PSYSTEMTIME LPSYSTEMTIME
#endif


class RTDeviceEx : public RTDevice
{	public:
		RTDeviceEx():RTDevice() {};	
		RTDeviceEx(ULONG Index):RTDevice(Index) {};	
		RTDeviceEx(HANDLE Attach):RTDevice(Attach){};	 
	

		//BOOLEAN FrontendControl(IN OUT PFRONTENDCONTROL fc );
		BOOLEAN Reference( BOOL RefOn );
	
		BOOLEAN GetRtcTime(IN PSYSTEMTIME Time );
		BOOLEAN GetRtcAlarmTime(IN PSYSTEMTIME Time );
		BOOLEAN SetRtcTime(IN PSYSTEMTIME Time ); 
		BOOLEAN SetRtcAlarmTime(IN PSYSTEMTIME Time );
		
		BOOLEAN MeasuringMode( ULONG Mode, ULONG Info );  
		BOOLEAN SetChanType( IN PULONG Type, IN ULONG NrOfChan ); 
		BOOLEAN SetHighPass( IN float *HighPass, IN ULONG NrOfChan ); 
		BOOLEAN SetLowPass(  IN float *LowPass , IN ULONG NrOfChan ); 
		BOOLEAN SetGain(     IN float *Gain    , IN ULONG NrOfChan ); 
		BOOLEAN SetOffset(   IN float *Offset ,  IN ULONG NrOfChan ); 
		
		BOOLEAN WriteId( IN ULONG StartAddress , IN ULONG Length , OUT VOID *Data );
		BOOLEAN ReadId( IN ULONG StartAddress , IN ULONG Length , OUT VOID *Data ); 
		
		//Additional features added for the Uro screener project 
		BOOLEAN SetIo( IN ULONG StartAddress , IN ULONG Length , OUT VOID *Data );
		BOOLEAN WriteMemory( IN ULONG StartAddress , IN ULONG Length , OUT VOID *Data );
		BOOLEAN ReadMemory( IN ULONG StartAddress , IN ULONG Length , OUT VOID *Data ); 
		BOOLEAN SetCorrection( IN ULONG Channel , IN float Gain , IN float Offset ); 
		BOOLEAN GetCorrection( IN ULONG Channel , OUT float *Gain , OUT float *Offset );
		BOOLEAN Storage( IN ULONG NewState , OUT ULONG *OldState);
		
};

void ShowError( const LPTSTR String ) ;

void MemCopy( LPVOID Dest , LPVOID Source , UINT uBytes ) ;
BOOLEAN FreeMem(LPCVOID Mem);
LPVOID AllocateMem( UINT uBytes ) ;

#endif // __FEATURE_H__
