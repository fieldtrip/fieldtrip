/*++

   Copyright (c) 2000 TMS International

   Module Name:

    SADIO.H

   Abstract:
   Header for all Signal acquisition devices


   Version:  1.01

   Author:

    Rob Hemstede

   Environment:

    Kernel mode / User mode

   Revision History:

   1 - Version 1,00: first release
   Vesrion 1.01: Addional features


   --*/



#ifndef __SADIO_H__
#define __SADIO_H__




/*
 #include <initguid.h>

   DEFINE_GUID(GUID_CLASS_SAD,
   0x323e3f01, 0x45ce, 0x4753, 0xb0, 0xa4, 0xc9, 0x2f, 0x35, 0x15, 0x62, 0x30);

 */


/*
 #ifndef CTL_CODE
 #pragma message ("CTL_CODE undefined. Include winioctl.h or wdm.h")
 #endif
 */



//#define SADIO_DEVICE_INFO     CTL_CODE(FILE_DEVICE_UNKNOWN,0x0,	METHOD_BUFFERED,FILE_ANY_ACCESS)
//#define SADIO_ADD_SLAVE			CTL_CODE(FILE_DEVICE_UNKNOWN,0x1,	METHOD_BUFFERED,FILE_ANY_ACCESS)
//#define SADIO_CONFIG_INFO			CTL_CODE(FILE_DEVICE_UNKNOWN,0x2,	METHOD_BUFFERED,FILE_ANY_ACCESS)
//#define SADIO_SET_CONFIG      CTL_CODE(FILE_DEVICE_UNKNOWN,0x4,	METHOD_BUFFERED,FILE_ANY_ACCESS)
//#define SADIO_RESET_DEVICE		CTL_CODE(FILE_DEVICE_UNKNOWN,0x3,	METHOD_BUFFERED,FILE_ANY_ACCESS)
//#define SADIO_SIGNAL_INFO			CTL_CODE(FILE_DEVICE_UNKNOWN,0x5,   METHOD_BUFFERED,FILE_ANY_ACCESS)
//#define SADIO_SET_DEVICE_STATE	CTL_CODE(FILE_DEVICE_UNKNOWN,0x6,   METHOD_BUFFERED,FILE_ANY_ACCESS)
//#define SADIO_GET_SAMPLE      CTL_CODE(FILE_DEVICE_UNKNOWN,0x7,   METHOD_BUFFERED,FILE_ANY_ACCESS)
//#define SADIO_BUFFER_STATE		CTL_CODE(FILE_DEVICE_UNKNOWN,0x8,   METHOD_BUFFERED,FILE_ANY_ACCESS)




/*******************************************************************************/
#ifdef CTL_CODE
#define SADIO_DEVICE_INFO     CTL_CODE(FILE_DEVICE_UNKNOWN,0x0, METHOD_BUFFERED,FILE_ANY_ACCESS)
#endif

// Get device info
// Input (DEVICE_INFO_PROPERTY)
// OUTPUT (WCHAR)

typedef enum {
        DeviceDescription,
        HardwareID,
        CompatibleIDs,
        BootConfiguration,
        BootConfigurationTranslated,
        ClassName,
        ClassGuid,
        DriverKeyName,
        Manufacturer,
        FriendlyName,
        LocationInformation,
        PhysicalDeviceObjectName,
        BusTypeGuid,
        LegacyBusType,
        BusNumber,
        EnumeratorName,
        Address,
        UINumber
} DEVICE_INFO_PROPERTY;


/*******************************************************************************/
#ifdef CTL_CODE
#define SADIO_ADD_SLAVE     CTL_CODE(FILE_DEVICE_UNKNOWN,0x1, METHOD_BUFFERED,FILE_ANY_ACCESS)
#endif
// Get device info


/*******************************************************************************/


#ifdef CTL_CODE
#define SADIO_CONFIG_INFO   CTL_CODE(FILE_DEVICE_UNKNOWN,0x2, METHOD_BUFFERED,FILE_ANY_ACCESS)
#endif
// Get configuration info
// Input (ULONG) Index.
// Output  (CONFIG_INFO) Config

#ifdef CTL_CODE
#define SADIO_SET_CONFIG    CTL_CODE(FILE_DEVICE_UNKNOWN,0x4, METHOD_BUFFERED,FILE_ANY_ACCESS)
#endif
// Set configuration
// Input (ULONG) index

#define CONFIG_SET_MODE         ((ULONG)0x02420001)

typedef struct _MEASURE_MODE
{ ULONG ConfigId;       //Must be CONFIG_SET_
  ULONG Size;         //Must be sizeof(MEASURE_MODE)
  ULONG Mode;
  ULONG Info;}MEASURE_MODE,*PMEASURE_MODE;

#define MEASURE_MODE_NORMAL     ((ULONG)0x0)
#define MEASURE_MODE_IMPEDANCE    ((ULONG)0x1)
#define MEASURE_MODE_CALIBRATION  ((ULONG)0x2)
#define MEASURE_MODE_IMPEDANCE_EX ((ULONG)0x3)
#define MEASURE_MODE_CALIBRATION_EX ((ULONG)0x4)

#define MeasureModeDefault(m)    \
        m.ConfigId = CONFIG_SET_MODE;  \
        m.Size=sizeof(MEASURE_MODE);   \
        m.Mode=MEASURE_MODE_NORMAL;    \
        m.Info=0;


/*
   typedef struct _CONFIG_INFO
   {
   ULONG Size;
   ULONG Index;
   ULONG Type;
   WCHAR Name[CONFIGNAME];

   }CONFIG_INFO, *PCONFIG_INFO;
 */


/*******************************************************************************/

#ifdef CTL_CODE
#define SADIO_RESET_DEVICE    CTL_CODE(FILE_DEVICE_UNKNOWN,0x3, METHOD_BUFFERED,FILE_ANY_ACCESS)
#endif

// Get reset the hardware (for debugging only)
// Input (VOID)
// Output (VOID)

/*******************************************************************************/

#ifdef CTL_CODE
#define SADIO_SIGNAL_INFO   CTL_CODE(FILE_DEVICE_UNKNOWN,0x5,   METHOD_BUFFERED,FILE_ANY_ACCESS)
#endif

// Get Signal info

#define SIGNAL_NAME 40


//Signal types en sub types

#define ST_UNKNOWN  0x0

#define ST_NODE       0xF
#define SST_SINGLE      0x0
#define SST_MASTER      0x1
#define SST_SLAVE     0x2

#define ST_TIME   0x10        //Signal Type
#define SST_YEARS   0x00000006    //Signal sub type
#define SST_MONTH   0x00000005
#define SST_DAYS    0x00000004
#define SST_HOURS   0x00000002
#define SST_MINUTES   0x00000001
#define SST_SECONDS   0x00000000    // Default subtype
#define SST_FRAMES    0x80000000

#define ST_EXG    0x11
#define SST_EXG     0x00000000    // Default sub type
#define SST_EEG     0x00000001
#define SST_EMG     0x00000002
#define SST_ECG     0x00000003
#define SST_EGG     0x00000004

#define ST_AUX    0x12
#define SST_AUX     0x00000000
#define SST_PRESURE   0x00000001

#define ST_DIG    0x13
#define SST_DIG     0x00000000


//** Signal Formats

#define SF_UNKNOWN  0x0   // undefined
#define SF_INTEGER  0x1   // signed integer
#define SF_FLOAT    0x2   // IEEE floating point
#define SF_ASCII  0x3   // std 8-bit ascii sting
#define SF_UNICODE  0x4   // 16bit UNICODE string
#define SF_UNSIGNED 0x5   // Unsigned integer


typedef struct _SCALE
{ FLOAT UnitsHigh;
  FLOAT UnitsLow;
  FLOAT AdcHigh;
  FLOAT AdcLow;}SCALE,*PSCALE;


typedef struct _SIGNAL_INFO
{
        ULONG Size; // Size of this structure including sub-structures
        ULONG Type; // One of the signal types above
        ULONG SubType; // One of the signal sub-types above
        ULONG Format; // Float / Integer / Asci / Ect..
        ULONG Bytes; // Number of bytes per sample including subsignals

        SCALE Scale; // Unit to Bit conversion

        WCHAR SignalName[SIGNAL_NAME];
        WCHAR UnitName[SIGNAL_NAME];

        ULONG IndexNext; // Point to next signal in the list;
        ULONG IndexDown; // Point to 1st sub-signal

}SIGNAL_INFO,*PSIGNAL_INFO;


/*******************************************************************************/

#ifdef CTL_CODE
#define SADIO_DEVICE_STATE  CTL_CODE(FILE_DEVICE_UNKNOWN,0x6,   METHOD_BUFFERED,FILE_ANY_ACCESS)
#endif

//	Input (ULONG) State
//  Output (ULONG) OldState


#define DEVICE_MAX_KNOWN_STATES  2

#define DEVICE_STATE_UNKNOWN  0xFFFFFFFF
#define DEVICE_STATE_STOP   0x0
#define DEVICE_STATE_START    0x1


/*******************************************************************************/

#ifdef CTL_CODE
#define SADIO_GET_SAMPLE    CTL_CODE(FILE_DEVICE_UNKNOWN,0x7,   METHOD_BUFFERED,FILE_ANY_ACCESS)
#endif

//	Input (ULONG) SampleNumber



/*******************************************************************************/

#ifdef CTL_CODE
#define SADIO_BUFFER_STATE    CTL_CODE(FILE_DEVICE_UNKNOWN,0x8,   METHOD_BUFFERED,FILE_ANY_ACCESS)
#endif


typedef struct _BUF_STATE
{ ULONG Size;       //Out		Size of this structure
  ULONG SampleRate;   //in/out	milli Hz (mHz)
  ULONG BufferSize;   //in/out	number of samples
  ULONG Reserved;     //in/out	Must be 1
  ULONG Available;    //out		Number of available samples
  ULONG Overflow;     //out		incremented if an overflow
  ULONG BytesPerSample; //out
}BUF_STATE, *PBUF_STATE;

/********************************************************************************/

#ifdef CTL_CODE
#define SADIO_SET_SYNC    CTL_CODE(FILE_DEVICE_UNKNOWN,0x9,   METHOD_BUFFERED,FILE_ANY_ACCESS)
#endif


typedef VOID (*SADIO_PSYNC_SIGNALS)(IN PVOID Params);

typedef struct _SADIO_SYNC_SIGNAL_PARMAMS
{
        ULONG Size;
        PVOID Fdo;
        SADIO_PSYNC_SIGNALS Function;
        BOOLEAN Master;
        ULARGE_INTEGER SampleNumber;

}SADIO_SYNC_SIGNAL_PARAMS,*PSADIO_SYNC_SIGNAL_PARAMS;


//*******************************************************************************/

/********************************************************************************/

#ifdef CTL_CODE
#define SADIO_FEATURE     CTL_CODE(FILE_DEVICE_UNKNOWN,0x10,   METHOD_BUFFERED,FILE_ANY_ACCESS)
#endif



/*****************************************************************************/

#ifdef CTL_CODE
#define SADIO_QUERY_DEVICE_KEY  CTL_CODE( FILE_DEVICE_UNKNOWN,0x11,   METHOD_BUFFERED,FILE_ANY_ACCESS)
#endif

typedef struct _KEY_VALUE_INFORMATION {
        ULONG Length;
        ULONG Type;
        ULONG DataLength;
        UCHAR Data[1];        // Variable size
} KEY_VALUE_INFORMATION, *PKEY_VALUE_INFORMATION;

#ifndef REG_NONE // added by Stefan Klanke for compilation using MinGW
#define REG_NONE                    ( 0 )   // No value type
#define REG_SZ                      ( 1 )   // Unicode nul terminated string
#define REG_EXPAND_SZ               ( 2 )   // Unicode nul terminated string
                                            // (with environment variable references)
#define REG_BINARY                  ( 3 )   // Free form binary
#define REG_DWORD                   ( 4 )   // 32-bit number
#define REG_DWORD_LITTLE_ENDIAN     ( 4 )   // 32-bit number (same as REG_DWORD)
#define REG_DWORD_BIG_ENDIAN        ( 5 )   // 32-bit number
#define REG_LINK                    ( 6 )   // Symbolic Link (unicode)
#define REG_MULTI_SZ                ( 7 )   // Multiple Unicode strings
#define REG_RESOURCE_LIST           ( 8 )   // Resource list in the resource map
#define REG_FULL_RESOURCE_DESCRIPTOR ( 9 )  // Resource list in the hardware description
#define REG_RESOURCE_REQUIREMENTS_LIST ( 10 )
#define REG_QWORD                   ( 11 )  // 64-bit number
#define REG_QWORD_LITTLE_ENDIAN     ( 11 )  // 64-bit number (same as REG_QWORD)
#endif // added by Stefan Klanke for compilation using MinGW

/*************************************************************************/

#ifdef CTL_CODE
#define SADIO_SIGNAL_FORMAT CTL_CODE( FILE_DEVICE_UNKNOWN,0x12,   METHOD_BUFFERED,FILE_ANY_ACCESS)
#endif

typedef struct _SIGNAL_FORMAT
{
        ULONG Size; // Size of this structure
        ULONG Elements; // Number of elements in list

        ULONG Type; // One of the signal types above
        ULONG SubType; // One of the signal sub-types above
        ULONG Format; // Float / Integer / Asci / Ect..
        ULONG Bytes; // Number of bytes per sample including subsignals

        FLOAT UnitGain;
        FLOAT UnitOffSet;
        ULONG UnitId;
        LONG UnitExponent;

        WCHAR Name[SIGNAL_NAME];

        ULONG Port;
        WCHAR PortName[SIGNAL_NAME];
        ULONG SerialNumber;

}SIGNAL_FORMAT, *PSIGNAL_FORMAT;


//DLL ProtoTypes;


typedef struct _SP_DEVICE_PATH {
        DWORD cbSize;
        TCHAR DevicePath[1];
} SP_DEVICE_PATH, *PSP_DEVICE_PATH;


#ifndef NO_DLL_PROTO

HANDLE APIENTRY Open     (PSP_DEVICE_PATH DevicePath);
BOOL APIENTRY Close      (HANDLE hHandle);
PSP_DEVICE_PATH APIENTRY GetDevicePath  (IN ULONG InterfaceNumber,OUT ULONG  *MaxInterfaces);
ULONG APIENTRY GetNrDevices ();
HANDLE APIENTRY GetSADHandle (IN ULONG InterfaceNumber,OUT ULONG  *MaxInterfaces);
PWCHAR APIENTRY GetDescription (IN HANDLE Handle,OUT PWCHAR Destination,IN ULONG Size);
ULONG APIENTRY GetId      (IN HANDLE Handle);
PWCHAR APIENTRY GetManufacturer(IN HANDLE Handle,OUT PWCHAR Destination,IN ULONG Size);
BOOLEAN APIENTRY Start      (IN HANDLE Handle);
BOOLEAN APIENTRY Stop     (IN HANDLE Handle);
ULONG APIENTRY GetDeviceState (IN HANDLE Handle);
PSIGNAL_INFO APIENTRY GetSignalInfo  (IN HANDLE Handle,OUT PSIGNAL_INFO pSignalInfo,IN OUT PULONG NrOfSignals);
BOOLEAN APIENTRY SetSignalBuffer(IN HANDLE Handle,IN OUT PULONG SampleRate,IN OUT PULONG BufferSize);
BOOLEAN APIENTRY GetBufferInfo  (IN HANDLE Handle,OUT PULONG Overflow,OUT PULONG PercentFull);
ULONG APIENTRY GetSamples   (IN HANDLE Handle,OUT PULONG SampleBuffer,IN ULONG Size);
ULONG APIENTRY BytesPerSample (IN HANDLE Handle);
BOOLEAN APIENTRY ResetDevice  (IN HANDLE Handle);
HANDLE APIENTRY GetSlaveHandle (IN HANDLE Handle);
BOOLEAN APIENTRY AddSlave   (IN HANDLE Handle,IN HANDLE SlaveHandle);
BOOLEAN APIENTRY DeviceFeature(IN HANDLE Handle,IN LPVOID DataIn, IN DWORD InSize,IN LPVOID DataOut, IN DWORD OutSize );
BOOLEAN APIENTRY GetDeviceKey(IN HANDLE Handle, IN PWCHAR Name, IN OUT PKEY_VALUE_INFORMATION Info );
PSIGNAL_FORMAT APIENTRY GetSignalFormat (IN HANDLE Handle,IN OUT PSIGNAL_FORMAT pSignalFormat);

PSP_DEVICE_PATH APIENTRY GetInstanceId  ( IN LONG DeviceIndex, IN BOOLEAN Present, OUT ULONG  *MaxDevices );
HKEY APIENTRY OpenRegKey   ( IN PSP_DEVICE_PATH Path );
BOOL APIENTRY Free       ( IN VOID *Memory );



typedef HANDLE ( __stdcall * POPEN )(PSP_DEVICE_PATH DevicePath);
typedef BOOL ( __stdcall * PCLOSE )(HANDLE hHandle);
typedef PSP_DEVICE_PATH ( __stdcall * PGETDEVICEPATH)(IN ULONG InterfaceNumber,OUT PULONG MaxInterfaces);
typedef ULONG ( __stdcall * PGETNRDEVICES)();
typedef HANDLE ( __stdcall * PGETSADHANDLE)(IN ULONG InterfaceNumber,OUT PULONG MaxInterfaces);
typedef PWCHAR ( __stdcall * PGETDESCRIPTION)(IN HANDLE Handle,OUT PWCHAR Destination,IN ULONG Size);
typedef ULONG ( __stdcall * PGETID)(IN HANDLE Handle);
typedef PWCHAR ( __stdcall * PGETMANUFACTURER)(IN HANDLE Handle,OUT PWCHAR Destination,IN ULONG Size);
typedef BOOLEAN ( __stdcall * PSTART)(IN HANDLE Handle);
typedef BOOLEAN ( __stdcall * PSTOP)(IN HANDLE Handle);
typedef ULONG ( __stdcall * PGETDEVICESTATE)(IN HANDLE Handle);
typedef PSIGNAL_INFO ( __stdcall * PGETSIGNALINFO)(IN HANDLE Handle,OUT PSIGNAL_INFO pSignalInfo,IN OUT PULONG NrOfSignals);
typedef BOOLEAN ( __stdcall * PSETSIGNALBUFFER)(IN HANDLE Handle,IN OUT PULONG SampleRate,IN OUT PULONG BufferSize);
typedef BOOLEAN ( __stdcall * PGETBUFFERINFO)(IN HANDLE Handle,OUT PULONG Overflow,OUT PULONG PercentFull);
typedef ULONG ( __stdcall * PGETSAMPLES)(IN HANDLE Handle,OUT PULONG SampleBuffer,IN ULONG Size);
typedef ULONG ( __stdcall * PBYTESPERSAMPLE)(IN HANDLE Handle);
typedef BOOLEAN ( __stdcall * PRESETDEVICE)(IN HANDLE Handle);
typedef HANDLE ( __stdcall * PGETSLAVEHANDLE)(IN HANDLE Handle);
typedef BOOLEAN ( __stdcall * PADDSLAVE)(IN HANDLE Handle,IN HANDLE SlaveHandle);
typedef BOOLEAN ( __stdcall * PDEVICEFEATURE)(IN HANDLE Handle,IN LPVOID DataIn, IN DWORD InSize,OUT LPVOID DataOut, IN DWORD OutSize );
typedef BOOLEAN ( __stdcall * PGETDEVICEKEY)(IN HANDLE Handle,IN PWCHAR Name, IN OUT PKEY_VALUE_INFORMATION Info );
typedef PSIGNAL_FORMAT ( __stdcall * PGETSIGNALFORMAT)(IN HANDLE Handle,IN OUT PSIGNAL_FORMAT pSignalFormat);

typedef PSP_DEVICE_PATH ( __stdcall * PGETINSTANCEID)(IN LONG DeviceIndex, IN BOOLEAN Present, OUT ULONG  *MaxDevices );
typedef HKEY ( __stdcall * POPENREGKEY)(IN PSP_DEVICE_PATH Path );
typedef BOOL ( __stdcall * PFREE)(IN VOID *Memory);

#endif //NO_DLL_PROTO


/* Additional types and constants used by feature.cpp */

typedef struct _FeatureData
{ ULONG Id;   /* Feature ID */
  ULONG Info;   /* Feature dependend information */
}FEATURE_DATA, *PFEATURE_DATA;

//----------- TYPE ---------------------


//constants used by set chantype
#define EXG (ULONG) 0x0001
#define AUX (ULONG) 0x0002

#define DEVICE_FEATURE_TYPE     0x0303

typedef struct _FeatureChanType
{ FEATURE_DATA Feature;
  ULONG Type[1];}FEATURE_TYPE,*PFEATURE_TYPE;

//------------ MODE ---------------------

#define DEVICE_FEATURE_MODE       0x0302

typedef struct _FeatureMode
{ FEATURE_DATA Feature;
  ULONG Mode;}FEATURE_MODE,*PFEATURE_MODE;

//--------------- RTC ----------------------

#define DEVICE_FEATURE_RTC        0x0301

typedef struct _FeatureRtc
{ FEATURE_DATA Feature;
  SYSTEMTIME Time;}FEATURE_RTC,*PFEATURE_RTC;

//---------- HIGHPASS ------------------

#define DEVICE_FEATURE_HIGHPASS   0x0401

typedef struct _FeatureHighpass
{ FEATURE_DATA Feature;
  float Highpass[1];}FEATURE_HIGHPASS, *PFEATURE_HIGHPASS;

//------------- LOWPASS ----------------

#define DEVICE_FEATURE_LOWPASS  0x0402

typedef struct _FeatureLowpass
{ FEATURE_DATA Feature;
  float Lowpass[1];}FEATURE_LOWPASS, *PFEATURE_LOWPASS;

//--------------- GAIN ------------------

#define DEVICE_FEATURE_GAIN   0x0403

typedef struct _FeatureGain
{ FEATURE_DATA Feature;
  float Gain[1];}FEATURE_GAIN, *PFEATURE_GAIN;

//--------------- OFFSET -------------------

#define DEVICE_FEATURE_OFFSET 0x0404

typedef struct _F_Offset
{ FEATURE_DATA Feature;
  float Offset[1];}FEATURE_OFFSET, *PFEATURE_OFFSET;


//------------------ IO ----------------------

#define DEVICE_FEATURE_IO 0x0500

typedef struct _FeatureMemory
{ FEATURE_DATA Feature;
  ULONG Data[1];}FEATURE_MEMORY,  *PFEATURE_MEMORY;

//----------------- MEMORY ----------------------

#define DEVICE_FEATURE_MEMORY 0x0501

//------------------- STORAGE ---------------------

#define DEVICE_FEATURE_STORAGE 0x0502

#define DFM_STORAGE_NOCHANGE  0xFFFF
#define DFM_STORAGE_OFF     0x0000
#define DFM_STORAGE_ON      0x0001

//------------------ CORRECTION --------------------

#define DEVICE_FEATURE_CORRECTION 0x0503

typedef struct _Feature_Correction
{ FEATURE_DATA Feature;
  float Gain;
  float Offset;}FEATURE_CORRECTION, *PFEATURE_CORRECTION;

//--------------------- ID ---------------------------

#define DEVICE_FEATURE_ID 0x0504


#endif //__SADIO_H__
