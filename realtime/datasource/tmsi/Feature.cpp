/*++

Copyright (c) 2001 TMS International

Module Name:
    Feature.cpp

Abstract:  
	Implementaion example for device specific features 
	
Version:  1.3

Author:
    Rob Hemstede

Environment:
    User mode

--*/

#include <Windows.h>
#include <stdio.h> 

#include "Feature.h" 


/*******************************************************************************/
// Set Channel Type 
// Example for overwriting the channel type settings 
/*******************************************************************************/

BOOLEAN RTDeviceEx::SetChanType( IN PULONG Type , IN ULONG NrOfChan )
{	BOOLEAN Result;
	ULONG Size = sizeof(FEATURE_DATA) + sizeof( ULONG ) * NrOfChan;
	PFEATURE_TYPE ChanType = ( PFEATURE_TYPE )AllocateMem( Size );
	if( ChanType == NULL ) return FALSE; 

	ChanType->Feature.Id = DEVICE_FEATURE_TYPE; 
	ChanType->Feature.Info = 0;  //Member is not used  
	MemCopy( ChanType->Type , Type , sizeof( ULONG ) * NrOfChan ); 

	Result = DeviceFeature( ChanType , Size , NULL , 0  ) ;
	FreeMem( ChanType ); 
	return Result; 

}

/********************************************************************************/
// Impedance and calibration example 
/*******************************************************************************/
/* 

	Decription 
		Set the measuring mode of the device 
	
	INPUT 
		Mode :	
				MEASURE_MODE_NORMAL	     = for normal measument 
				MEASURE_MODE_IMPEDANCE   = for impedance measurement 
				MEASURE_MODE_CALIBRATION = for calibration measurment 

		Info :	depending on the selected operating mode this can be  
			for MEASURE_MODE_NORMAL   , 0
			for MEASURE_MODE_IMPEDANCE 
					IC_OHM_002	// 2K Impedance limit	
					IC_OHM_005	// 5K Impedance limit
					IC_OHM_010	// 10K Impedance limit
					IC_OHM_020	// 20K Impedance limit
					IC_OHM_050	// 50K Impedance limit
					IC_OHM_100	// 100K Impedance limit
			for MEASURE_MODE_CALIBRATION
					IC_VOLT_050 0	//50 uV t-t Calibration voltage
					IC_VOLT_100 1	//100 uV t-t
					IC_VOLT_200 2	//200 uV t-t
					IC_VOLT_500 3	//500 uV t-t
			
	OUTPUT
		returns TRUE if successful
				
*/





BOOLEAN RTDeviceEx::MeasuringMode(ULONG Mode,ULONG Info)
{	FEATURE_MODE FMode;
	FMode.Feature.Id = DEVICE_FEATURE_MODE; 
	FMode.Feature.Info = Info; 
	FMode.Mode = Mode; 
	
	return DeviceFeature( &FMode, sizeof(FEATURE_MODE) , NULL , 0  ) ;
}


/********************************************************************/
/* Demo implementaion for handling RTC 
/********************************************************************/



BOOLEAN RTDeviceEx::GetRtcTime(OUT PSYSTEMTIME Time )
{	FEATURE_RTC RtcData; 
	ULONG InSize; 
	//Check input parameters 
	if( Time == NULL ) return FALSE; 

	InSize = sizeof( FEATURE_DATA ); 

	RtcData.Feature.Info = 0; //Set info member to '0' for handling normal time
	RtcData.Feature.Id = DEVICE_FEATURE_RTC; 

	if(DeviceFeature( &RtcData, sizeof(FEATURE_DATA) , &RtcData , sizeof( FEATURE_RTC ) ) ) 
	{	MemCopy( Time , &RtcData.Time , sizeof( SYSTEMTIME )  ); 
		return TRUE; 
	}
	else
	{	return FALSE; 
	}
}

// 
BOOLEAN RTDeviceEx::GetRtcAlarmTime(OUT PSYSTEMTIME Time )
{	FEATURE_RTC RtcData; 
	//Check input parameters 
	if( Time == NULL ) return FALSE; 
	
	RtcData.Feature.Info = 1; //Set info member to '1' for handling alarm time
	RtcData.Feature.Id = DEVICE_FEATURE_RTC; 

	if(DeviceFeature( &RtcData, sizeof(FEATURE_DATA) , &RtcData , sizeof( FEATURE_RTC ) ) ) 
	{	MemCopy( Time , &RtcData.Time , sizeof( SYSTEMTIME )  ); 
		return TRUE; 
	}
	else
	{	return FALSE; 
	}
}

BOOLEAN RTDeviceEx::SetRtcTime(IN PSYSTEMTIME Time )
{	FEATURE_RTC RtcData; 
	//Check input parameters 
	if( Time == NULL ) return FALSE; 

	RtcData.Feature.Info = 0; //Set info member to '1' for handling normal time
	RtcData.Feature.Id = DEVICE_FEATURE_RTC; 
	MemCopy( &RtcData.Time ,Time ,  sizeof( SYSTEMTIME )  ); 
	return DeviceFeature( &RtcData, sizeof(FEATURE_RTC) , NULL , 0  ) ;
}

BOOLEAN RTDeviceEx::SetRtcAlarmTime(IN PSYSTEMTIME Time )
{	FEATURE_RTC RtcData; 
	//Check input parameters 
	if( Time == NULL ) return FALSE; 

	RtcData.Feature.Info = 1; //Set info member to '1' for handling normal time
	RtcData.Feature.Id = DEVICE_FEATURE_RTC; 
	MemCopy( &RtcData.Time ,Time ,  sizeof( SYSTEMTIME )  ); 
	return DeviceFeature( &RtcData, sizeof(FEATURE_RTC) , NULL , 0  ) ;
}

/**********************************************************************/


BOOLEAN RTDeviceEx::SetHighPass( IN float *HighPass, IN ULONG NrOfChan )
{	BOOLEAN Result; 
	ULONG size = sizeof( FEATURE_DATA ) + sizeof( float ) * NrOfChan; 
	PFEATURE_HIGHPASS Data = (PFEATURE_HIGHPASS) AllocateMem( size ); 
	if( Data == NULL ) return FALSE; 
	MemCopy( &Data->Highpass[0] , HighPass ,  sizeof( float ) * NrOfChan ); 
	Data->Feature.Id = DEVICE_FEATURE_HIGHPASS;
	Data->Feature.Info = 0 ; // Info is not used 
	Result = DeviceFeature( Data, size , NULL , 0  ) ;
	FreeMem( Data ); 
	return Result;
} 

/**********************************************************************/


BOOLEAN RTDeviceEx::SetLowPass(  IN float *LowPass , IN ULONG NrOfChan )
{	BOOLEAN Result; 
	ULONG size = sizeof( FEATURE_DATA ) + sizeof( float ) * NrOfChan; 
	PFEATURE_LOWPASS Data = (PFEATURE_LOWPASS) AllocateMem( size ); 
	if( Data == NULL ) return FALSE; 
	MemCopy( &Data->Lowpass[0] , LowPass ,  sizeof( float ) * NrOfChan ); 
	Data->Feature.Id = DEVICE_FEATURE_LOWPASS;
	Data->Feature.Info = 0 ; // Info is not used 
	Result = DeviceFeature( Data, size , NULL , 0  ) ;
	FreeMem( Data ); 
	return Result;
} 

/**********************************************************************/

 

BOOLEAN RTDeviceEx::SetGain( IN float *Gain , IN ULONG NrOfChan )
{	BOOLEAN Result; 
	ULONG size = sizeof( FEATURE_DATA ) + sizeof( float ) * NrOfChan; 
	PFEATURE_GAIN Data = (PFEATURE_GAIN) AllocateMem( size ); 
	if( Data == NULL ) return FALSE; 
	MemCopy( &Data->Gain[0] , Gain ,  sizeof( float ) * NrOfChan ); 
	Data->Feature.Id = DEVICE_FEATURE_GAIN;
	Data->Feature.Info = 0 ; // Info is not used 
	Result = DeviceFeature( Data, size , NULL , 0  ) ;
	FreeMem( Data ); 
	return Result;
}
 

/**********************************************************************/



BOOLEAN RTDeviceEx::SetOffset(   IN float *Offset ,  IN ULONG NrOfChan )
{	BOOLEAN Result; 
	ULONG size = sizeof( FEATURE_DATA ) + sizeof( float ) * NrOfChan; 
	PFEATURE_OFFSET Data = (PFEATURE_OFFSET) AllocateMem( size ); 
	if( Data == NULL ) return FALSE; 
	MemCopy( &Data->Offset[0] , Offset ,  sizeof( float ) * NrOfChan ); 
	Data->Feature.Id = DEVICE_FEATURE_OFFSET;
	Data->Feature.Info = 0 ; // Info is not used 
	Result = DeviceFeature( Data, size , NULL , 0  ) ;
	FreeMem( Data ); 
	return Result;
} 


/**********************************************************************/

/*	
	SetIo  

	description : This fucntion is used to set device specific io data, 
	these bits are used. What device features are controled by these data 
	dependes on the IoFeature parameter.  

	note: 
		Use IoFeature 0x0001 and bit 0 of the data to contole the Leakage led 
		on the UroScreener devices 

*/ 

BOOLEAN RTDeviceEx::SetIo( IN ULONG IoFeature , IN ULONG Length , OUT VOID *Data )
{	BOOLEAN ReturnCode; 
	
	PFEATURE_MEMORY pIo = (PFEATURE_MEMORY) AllocateMem( sizeof( FEATURE_DATA ) + Length ); 
	
	if( pIo == NULL ) return FALSE; 
	
	if( (Data != NULL) && (Length != 0) ) 
		MemCopy( &(pIo->Data) , Data , Length ); 

	pIo->Feature.Id = DEVICE_FEATURE_IO; 
	pIo->Feature.Info = IoFeature;

	ReturnCode = DeviceFeature( pIo , sizeof( FEATURE_DATA ) + Length , NULL , 0 ) ; 
	
	FreeMem( pIo ); 

	return ReturnCode; 
}
 

/**********************************************************************/



/* 
	ReadMemory 
	Description 
		Function to read memory from a device
	Parameters 
		StartAddress:	start address of the data transfer 
		Length : Length of the transfer in bytes 
		* Data : pointer to the destination buffer 


*/

BOOLEAN RTDeviceEx::ReadMemory( IN ULONG StartAddress , IN ULONG Length , OUT VOID *Data )
{	BOOLEAN ReturnCode; 

	PFEATURE_MEMORY pIo = (PFEATURE_MEMORY) AllocateMem( sizeof( FEATURE_DATA ) + Length ); 
	if( pIo == NULL ) return FALSE; 
	
	pIo->Feature.Id = DEVICE_FEATURE_MEMORY; 
	pIo->Feature.Info = StartAddress;

	ReturnCode = DeviceFeature( pIo , sizeof( FEATURE_DATA ) , pIo , ( sizeof( FEATURE_DATA ) + Length )  ) ; 
	if( ReturnCode ) 
	{	if(( Data != NULL ) && (Length!=0)) 
		{
			MemCopy( Data , &(pIo->Data) , Length ); 
		}
		FreeMem( pIo ); 
	}	return TRUE; 
	FreeMem( pIo ); 

	return FALSE; 
}

/*
	WriteMemory 
	Description 
		Function to write the memory of a device
	Parameters 
		StartAddress:	start address of the data transfer 
		Length : Length of the transfer in bytes 
		* Data : pointer to the source buffer 
*/

BOOLEAN RTDeviceEx::WriteMemory( IN ULONG StartAddress , IN ULONG Length , IN VOID *Data )
{	BOOLEAN ReturnCode; 
	ULONG Size = sizeof( FEATURE_DATA ) + Length; 
		
	PFEATURE_MEMORY pIo = (PFEATURE_MEMORY) AllocateMem( Size ); 
	if( pIo == NULL ) return FALSE; 
	
	if( (Length!=0) && (Data!=NULL) ) 
		MemCopy( &(pIo->Data) , Data , Length ); 

	pIo->Feature.Id = DEVICE_FEATURE_MEMORY; 
	pIo->Feature.Info = StartAddress;

	ReturnCode = DeviceFeature( pIo , Size , NULL , 0);
	
	FreeMem( pIo ); 

	return ReturnCode; 
}
 
/***********************************************************/
  


/*
	Storage 
	Description: 
		function to control the local storage of a device, 
	Parameters 
		NewState : this paramter is used to set the new storage state 
		OldState : returns the previous storage state 

*/ 


BOOLEAN RTDeviceEx::Storage( IN ULONG NewState ,  OUT ULONG *OldState )
{	FEATURE_DATA Feature;
	BOOLEAN ReturnCode; 

	Feature.Id = DEVICE_FEATURE_STORAGE; 
	Feature.Info = NewState;

	ReturnCode = DeviceFeature( &Feature , sizeof( FEATURE_DATA ) , &Feature ,  sizeof( FEATURE_DATA ) ) ; 
	
	if( ReturnCode ) if( OldState != NULL ) *OldState = Feature.Info; 
	
	return ReturnCode; 
}

/***********************************************************************/
// TODO test functions below when implemeted by embedded software

/*
	SetCorrection 
	Description 
		Some devices support signal correction for calibaration. This fucntion is 
		used to set the correction parameters 

	Parameters 
		Channel :	Channel index for which to set the new correction parameters 
		Gain	:	New gain to apply for this channel 
		Offset	:	new offset to apply for the channel 
		
		When correcting the signal first the gain and second the offset is applied 

*/ 


BOOLEAN RTDeviceEx::SetCorrection( IN ULONG Channel , IN float Gain , IN float Offset )
{	FEATURE_CORRECTION corr; 
	corr.Feature.Id = DEVICE_FEATURE_CORRECTION; 
	corr.Feature.Info = Channel; 
	corr.Gain = Gain;
	corr.Offset = Offset; 
	return DeviceFeature( &corr , sizeof( FEATURE_CORRECTION ) , NULL , 0 ); 
}

BOOLEAN RTDeviceEx::GetCorrection( IN ULONG Channel , OUT float *Gain , OUT float *Offset )
{	FEATURE_CORRECTION corr; 
	memset( &corr, 0, sizeof( FEATURE_CORRECTION )); 
	corr.Feature.Id = DEVICE_FEATURE_CORRECTION; 
	corr.Feature.Info = Channel; 
	if( DeviceFeature( &corr , sizeof( FEATURE_DATA ) , &corr , sizeof(FEATURE_CORRECTION) ) ) 
	{	if( Gain != NULL ) *Gain = corr.Gain; 
		if( Offset != NULL ) * Offset = corr.Offset; 
		return TRUE; 
	}
	return FALSE; 
}


/**********************************************************************/



/* 
	This function is for internal use by TMS only 

	ReadMemory 
	Description 
		Function to read id memory structure from a device
	Parameters 
		StartAddress:	start address of the data transfer 
		Length : Length of the transfer in bytes 
		* Data : pointer to the destination buffer 


*/

BOOLEAN RTDeviceEx::ReadId( IN ULONG StartAddress , IN ULONG Length , OUT VOID *Data )
{	BOOLEAN ReturnCode; 

	PFEATURE_MEMORY pIo = (PFEATURE_MEMORY) AllocateMem( sizeof( FEATURE_DATA ) + Length ); 
	if( pIo == NULL ) return FALSE; 
	
	pIo->Feature.Id = DEVICE_FEATURE_ID; 
	pIo->Feature.Info = StartAddress;

	ReturnCode = DeviceFeature( pIo , sizeof( FEATURE_DATA ) , pIo , ( sizeof( FEATURE_DATA ) + Length )  ) ; 
	if( ReturnCode ) 
	{	if( Data != NULL )  
		{	MemCopy( Data , &(pIo->Data) , Length ); 
		}
		FreeMem( pIo ); 
	}	return TRUE; 
	FreeMem( pIo ); 
	return FALSE; 
}

/*
	WriteId
	Description 
		Function ID structure of a device
	Parameters 
		StartAddress:	start address of the data transfer 
		Length : Length of the transfer in bytes 
		* Data : pointer to the source buffer 
*/

BOOLEAN RTDeviceEx::WriteId( IN ULONG StartAddress , IN ULONG Length , OUT VOID *Data )
{	BOOLEAN ReturnCode; 
	
	PFEATURE_MEMORY pIo = (PFEATURE_MEMORY) AllocateMem( sizeof( FEATURE_DATA ) + Length ); 
	if( pIo == NULL ) return FALSE; 
	
	if( Data != NULL ) MemCopy( &(pIo->Data) , Data , Length ); 

	pIo->Feature.Id = DEVICE_FEATURE_ID; 
	pIo->Feature.Info = StartAddress;

	ReturnCode = DeviceFeature( pIo , sizeof( FEATURE_DATA ) , pIo , ( sizeof( FEATURE_DATA ) + Length )  ) ; 
	
	FreeMem( pIo ); 

	return ReturnCode; 
}
 
/***********************************************************/



/********************************************************************/
// Example for debugging by TMS International only. 
// Het luikje van tms , alleen voor intern gebruik 
/********************************************************************/
/*
#define DEVICE_FEATURE_DSP_COMMAND  0x0100

//DSP Funtion codes 
#define  FRONTEND_CONTROL_ID 0x0003

typedef struct FrontendControl
{	WORD NuberOfChan;
	WORD ConversionType;
	WORD BlockSize; 
	WORD BufferSize;
	WORD SampleRateDiv;
	WORD Mode;
}FRONTENDCONTROL, *PFRONTENDCONTROL;

typedef struct _F_FrontendControl
{	FEATURE_DATA Feature;
	FRONTENDCONTROL Control; 
}F_FRONTENDCONTROL, *PF_FRONTENDCONTROL; 


BOOLEAN DeviceSpecificFeature::FrontendControl( IN OUT PFRONTENDCONTROL fc )
{	BOOLEAN Result; 
	F_FRONTENDCONTROL FControl; 
	
	if( fc == NULL ) return FALSE;
	MemCopy( &FControl.Control , fc , sizeof( FRONTENDCONTROL ) ); 

	//Initialize the feature data structure 
	FControl.Feature.Id = DEVICE_FEATURE_DSP_COMMAND;
	FControl.Feature.Info = FRONTEND_CONTROL_ID;
	
	Result = DeviceFeature( &FControl , sizeof(FControl) , &FControl , sizeof(FControl) ); 

	return Result;
}

*/

#define DEVICE_FEATURE_REF 0x0101

BOOLEAN RTDeviceEx::Reference( BOOL RefOn )
{	FEATURE_DATA fd; 
	
	fd.Id = DEVICE_FEATURE_REF; 
	fd.Info = RefOn; 
	
	return DeviceFeature( &fd, sizeof( FEATURE_DATA ) , NULL , 0); 

}




/*******************************************************************************/
// Simple error analyses example */ 
/*******************************************************************************/


void ShowError( const LPTSTR String ) 
{	LPTSTR ErrorMessage = NULL; 
	DWORD Error = 0;
	/* Get the error code from the system*/ 
	Error = GetLastError() ;
	
	/* Let the system convert this code to test*/ 
	FormatMessage( FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_ALLOCATE_BUFFER , 
				   NULL , Error , 0 , (char *)&ErrorMessage , 10 , NULL ); 
	
	if( String != NULL ) printf("\n%s : ", String );
	
	if( ErrorMessage != NULL )
	{	printf("0x%x %s", Error ,ErrorMessage ); 
		LocalFree( ErrorMessage );
	}
	else
	{	printf("Cause is unknown"); 
	}
	return;
}


/******************************************************************/
/* Some general memory util. functions 
/******************************************************************/ 

LPVOID AllocateMem( UINT uBytes ) 
{
	return (LPVOID) LocalAlloc( LMEM_FIXED | LMEM_ZEROINIT , uBytes ); 

}

BOOLEAN FreeMem(LPCVOID Mem)
{	HLOCAL Handle = NULL; 
	if( Mem == NULL) return FALSE; 
	Handle = LocalHandle( Mem ); 
	if( Handle == NULL ) return FALSE; 
	Handle = LocalFree( Handle ); 
	if( Handle != NULL ) return FALSE; 
	return TRUE;
}


void MemCopy( LPVOID Dest , LPVOID Source , UINT uBytes ) 
{
	UINT i; 
	PCHAR d = (PCHAR)Dest; 
	PCHAR s = (PCHAR)Source; 
	
	if( (d != NULL) && (s!=NULL) ) 
		for(i=0; i<uBytes ;i++ ) 
			d[i] = s[i];
	
	return; 
}

/*************************************************************/

