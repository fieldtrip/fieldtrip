///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2003-2004 Neuroshare Project
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// A copy of the GNU Lesser General Public License can be obtained by writing to:
//  Free Software Foundation, Inc.,
//  59 Temple Place, Suite 330,
//  Boston, MA  02111-1307
//  USA
//
// Contact information:
//  Kirk Korver
//  CyberKinetics, Inc.,
//  391 G Chipeta Way
//  Salt Lake City,  UT  84108
//  USA
//  kkorver at cyberkineticsinc.com
//
// Website:
//  sourceforge.net/projects/neuroshare
//
// All other copyrights on this material are replaced by this license agreeement.
//
///////////////////////////////////////////////////////////////////////////////////////////////////
//
// $Author: putermutt $
// $Date: 2007/01/24 00:50:10 $
// $Revision: 1.2 $
// $Source: /cvsroot/neuroshare/NSLibraries/ns/ns.h,v $
//
///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Specification : API definitions for the Neuroshare file format standard.
//
// Description   : Neuroshare API Glue
// Use the following files to access Neuroshare files:
//   ns.h        --  API header file
//   ns.c        --  API interface module
//
// A Neuroshare library is required in conjuction with this module. 
// For example:
//   ns_NEV.so   --  Linux shared object
//      or
//   ns_NEV.dll  --  Windows dynamic link library
//
// Authors       : Gopal Santhanam
//
// $Date: 2007/01/24 00:50:10 $
//
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef NS_H_INCLUDED   // Include guards
#define NS_H_INCLUDED


/*=========================================================================
| VERSION
 ========================================================================*/
#define NS_HEADER_VERSION  0x0103  /* v1.03 */


/*=========================================================================
| C++ SUPPORT
 ========================================================================*/
#ifdef __cplusplus
extern "C" {
#endif


/*=========================================================================
| TYPES
 ========================================================================*/
#if defined(__GNUC__)
    typedef char           INT8;
    typedef unsigned char  UINT8;
    typedef short          INT16;
    typedef unsigned short UINT16;
    typedef int            INT32;
    typedef unsigned int   UINT32;
    typedef unsigned long long  UINT64;

    // Specify 4 byte structure alignment
    #define GCC_ALIGN_4     __attribute__((aligned(4)))
    
    
#elif defined(_MSC_VER)
    #include <windows.h>
    
    // Copied from the Feb2003 SDK, just in case it isn't installed
    typedef signed char         INT8, *PINT8;
    typedef signed short        INT16, *PINT16;
    typedef signed int          INT32, *PINT32;
    typedef signed __int64      INT64, *PINT64;
    typedef unsigned char       UINT8, *PUINT8;
    typedef unsigned short      UINT16, *PUINT16;
    typedef unsigned int        UINT32, *PUINT32;
    typedef unsigned __int64    UINT64, *PUINT64;


    // specify 4-byte structure alignment
    #pragma pack(push, 4)       
    #define GCC_ALIGN_4

#endif



///////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Library Return Code Definitions
//
///////////////////////////////////////////////////////////////////////////////////////////////////

typedef    INT32 ns_RESULT;

typedef enum {
    ns_OK        =  0,  //Function Successful
    ns_LIBERROR  = -1,  //Linked Library Error
    ns_TYPEERROR = -2,  //Library unable to open file type
    ns_FILEERROR = -3,  //File access or read Error
    ns_BADFILE   = -4,  //Invalid file handle passed to function
    ns_BADENTITY = -5,  //Invalid or inappropriate entity identifier specified
    ns_BADSOURCE = -6,  //Invalid source identifier specified
    ns_BADINDEX  = -7   //Invalid entity index or index range specified
} ns_STATUS;
    
    
///////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Definitions of constants and flags 
//
///////////////////////////////////////////////////////////////////////////////////////////////////

// Library description flags
typedef enum {
    ns_LIBRARY_DEBUG         = 0x01,  // includes debug info linkage
    ns_LIBRARY_MODIFIED      = 0x02,  // file was patched or modified
    ns_LIBRARY_PRERELEASE    = 0x04,  // pre-release or beta version
    ns_LIBRARY_SPECIALBUILD  = 0x08,  // different from release version
    ns_LIBRARY_MULTITHREADED = 0x10,  // library is multithread safe
} ns_LIBRARY_FLAGS;

//Definitions of Event Entity types
typedef enum {
    ns_EVENT_TEXT  = 0,  // null-terminated ascii text string
    ns_EVENT_CSV   = 1,  // comma separated ascii text values
    ns_EVENT_BYTE  = 2,  // 8-bit value
    ns_EVENT_WORD  = 3,  // 16-bit value
    ns_EVENT_DWORD = 4,  // 32-bit value
} ns_EVENT_TYPE;

//Definitions of entity types in the structure ns_ENTITYINFO
typedef enum {
    ns_ENTITY_UNKNOWN     = 0,     // unknown entity type
    ns_ENTITY_EVENT       = 1,     // Event entity
    ns_ENTITY_ANALOG      = 2,     // Analog entity
    ns_ENTITY_SEGMENT     = 3,     // Segment entity
    ns_ENTITY_NEURALEVENT = 4,     // Sorted Neural entity
} ns_ENTITY_TYPE;

//Flags used for locating data entries
typedef enum {
    ns_BEFORE  = -1,  // less than or equal to specified time    
    ns_CLOSEST =  0,  // closest time 
    ns_AFTER   = +1,  // greater than or equal to specified time
} ns_LOCATION_FLAGS;


///////////////////////////////////////////////////////////////////////////////////////////////////
// 
//             DLL library version information functions
//
///////////////////////////////////////////////////////////////////////////////////////////////////

//File descriptor structure
typedef struct
{
    char szDescription[32];    // Text description of the file type or file family
    char szExtension[8];       // Extension used on PC, Linux, and Unix Platforms
    char szMacCodes[8];        // Application and Type Codes used on Mac Platforms
    char szMagicCode[16];      // Null-terminated code used at the file beginning
} GCC_ALIGN_4  ns_FILEDESC;

// Library information structure
typedef struct 
{
    UINT32 dwLibVersionMaj;    // Major version number of library
    UINT32 dwLibVersionMin;    // Minor version number of library
    UINT32 dwAPIVersionMaj;    // Major version number of API
    UINT32 dwAPIVersionMin;    // Minor version number of API
    char   szDescription[64];  // Text description of the library
    char   szCreator[64];      // Name of library creator
    UINT32 dwTime_Year;        // Year of last modification date
    UINT32 dwTime_Month;       // Month (1-12; January = 1) of last modification date
    UINT32 dwTime_Day;         // Day of the month (1-31) of last modification date
    UINT32 dwFlags;            // Additional library flags
    UINT32 dwMaxFiles;         // Maximum number of files library can simultaneously open
    UINT32 dwFileDescCount;    // Number of valid description entries in the following array
    ns_FILEDESC FileDesc[16];  // Text descriptor of files that the DLL can interpret
} GCC_ALIGN_4  ns_LIBRARYINFO;

// File information structure (the time of file creation should be reported in GMT)
typedef struct
{
    char   szFileType[32];         // Manufacturer's file type descriptor
    UINT32 dwEntityCount;          // Number of entities in the data file.
    double dTimeStampResolution;   // Minimum timestamp resolution
    double dTimeSpan;              // Time span covered by the data file in seconds
    char   szAppName[64];          // Name of the application that created the file
    UINT32 dwTime_Year;            // Year the file was created
    UINT32 dwTime_Month;           // Month (1-12; January = 1)
    UINT32 dwTime_DayofWeek;       // Day of the week (0-6; Sunday = 0)
    UINT32 dwTime_Day;             // Day of the month (1-31)
    UINT32 dwTime_Hour;            // Hour since midnight (0-23)
    UINT32 dwTime_Min;             // Minute after the hour (0-59)
    UINT32 dwTime_Sec;             // Seconds after the minute (0-59)
    UINT32 dwTime_MilliSec;        // Milliseconds after the second (0-1000)
    char   szFileComment[256];     // Comments embedded in the source file
} GCC_ALIGN_4  ns_FILEINFO;

// General entity information structure
typedef struct
{
    char   szEntityLabel[32];  // Specifies the label or name of the entity
    UINT32 dwEntityType;       // One of the ns_ENTITY_* types defined above
    UINT32 dwItemCount;        // Number of data items for the specified entity in the file
} ns_ENTITYINFO;
    
// Event entity information structure
typedef struct
{
    UINT32 dwEventType;      // One of the ns_EVENT_* types defined above
    UINT32 dwMinDataLength;  // Minimum number of bytes that can be returned for an Event
    UINT32 dwMaxDataLength;  // Maximum number of bytes that can be returned for an Event
    char   szCSVDesc[128];   // Description of the data fields for CSV Event Entities
} GCC_ALIGN_4  ns_EVENTINFO;

// Analog information structure
typedef struct
{
    double    dSampleRate;           // The sampling rate in Hz used to digitize the analog values
    double    dMinVal;               // Minimum possible value of the input signal
    double    dMaxVal;               // Maximum possible value of the input signal
    char    szUnits[16];             // Specifies the recording units of measurement
    double    dResolution;           // Minimum resolvable step (.0000305 for a +/-1V 16-bit ADC)  
    double    dLocationX;            // X coordinate in meters
    double    dLocationY;            // Y coordinate in meters
    double    dLocationZ;            // Z coordinate in meters
    double    dLocationUser;         // Additional position information (e.g. tetrode number)
    double    dHighFreqCorner;       // High frequency cutoff in Hz of the source signal filtering
    UINT32    dwHighFreqOrder;       // Order of the filter used for high frequency cutoff
    char    szHighFilterType[16];    // Type of filter used for high frequency cutoff (text format)
    double    dLowFreqCorner;        // Low frequency cutoff in Hz of the source signal filtering
    UINT32    dwLowFreqOrder;        // Order of the filter used for low frequency cutoff
    char    szLowFilterType[16];     // Type of filter used for low frequency cutoff (text format)
    char    szProbeInfo[128];        // Additional text information about the signal source
} GCC_ALIGN_4  ns_ANALOGINFO;
    
//Segment Information structure
typedef struct
{
    UINT32 dwSourceCount;     // Number of sources in the Segment Entity, e.g. 4 for a tetrode
    UINT32 dwMinSampleCount;  // Minimum number of samples in each Segment data item
    UINT32 dwMaxSampleCount;  // Maximum number of samples in each Segment data item
    double dSampleRate;       // The sampling rate in Hz used to digitize source signals
    char   szUnits[32];       // Specifies the recording units of measurement
} GCC_ALIGN_4  ns_SEGMENTINFO;

// Segment source information structure
typedef struct
{
    double dMinVal;               // Minimum possible value of the input signal
    double dMaxVal;               // Maximum possible value of the input signal
    double dResolution;           // Minimum input step size that can be resolved
    double dSubSampleShift;       // Time diff btn timestamp and actual sampling time of source
    double dLocationX;            // X coordinate of source in meters
    double dLocationY;            // Y coordinate of source in meters
    double dLocationZ;            // Z coordinate of source in meters
    double dLocationUser;         // Additional position information (e.g tetrode number)
    double dHighFreqCorner;       // High frequency cutoff in Hz of the source signal filtering
    UINT32 dwHighFreqOrder;       // Order of the filter used for high frequency cutoff
    char   szHighFilterType[16];  // Type of filter used for high frequency cutoff (text format)
    double dLowFreqCorner;        // Low frequency cutoff in Hz of the source signal filtering
    UINT32 dwLowFreqOrder;        // Order of the filter used for low frequency cutoff
    char   szLowFilterType[16];   // Type of filter used for low frequency cutoff (text format)
    char   szProbeInfo[128];      // Additional text information about the signal source
} GCC_ALIGN_4  ns_SEGSOURCEINFO;

// Neural Information structure
typedef struct
{
    UINT32 dwSourceEntityID;  // Optional ID number of a source entity
    UINT32 dwSourceUnitID;    // Optional sorted unit ID number used in the source entity
    char   szProbeInfo[128];  // Additional probe text information or source entity label
} GCC_ALIGN_4  ns_NEURALINFO;

typedef INT32 ns_DLLHANDLE;

    
#if !defined(NS_COMPILING_LIB)
#  define NS_FIRST_ARG ns_DLLHANDLE nssDllHANDLE,
#  define ns_stdcall
#else
#  define NS_FIRST_ARG
#  if defined(WIN32) || defined(_WIN32)
#    define ns_stdcall __stdcall
#  else
#    define ns_stdcall
#  endif
#endif
    
    

/*=========================================================================
| PROTOTYPES
 ========================================================================*/
ns_DLLHANDLE ns_stdcall ns_LoadLibrary  (const char *libname);
ns_RESULT    ns_stdcall ns_CloseLibrary (ns_DLLHANDLE nsDllHandle);
    

///////////////////////////////////////////////////////////////////////////////////////////////////
// ns_GetLibraryInfo
//  
// Purpose:
//  Retrieves information about the loaded API library
//
// Parameters:
//  ns_LIBRARYINFO *pLibraryInfo       pointer to ns_LIBRARYINFO structure to receive information
//  UINT32 dwLibraryInfoSize           size in bytes of ns_LIBRARYINFO structure
//
// Return Values:
//  ns_OK                              ns_LIBIRARYINFO successfully retrieved
//  ns_LIBEERROR                       library error
//
///////////////////////////////////////////////////////////////////////////////////////////////////

ns_RESULT ns_stdcall ns_GetLibraryInfo
        (NS_FIRST_ARG ns_LIBRARYINFO *pLibraryInfo, UINT32 dwLibraryInfoSize);



///////////////////////////////////////////////////////////////////////////////////////////////////
// ns_OpenFile
//  
// Purpose:
//    Opens the data file and assigns a file handle for internal use by the library.
//
// Parameters:
//  char *pszFilename                  name of file to open
//  UINT32 *hFile                      pointer to a file handle
//
// Return Values:
//  ns_OK                              ns_LIBIRARYINFO successfully retrieved
//  ns_TYPEERROR                       library unable to open file type
//  ns_FILEERROR                       file access or read error 
//  ns_LIBEERROR                       library error
//
///////////////////////////////////////////////////////////////////////////////////////////////////

ns_RESULT ns_stdcall ns_OpenFile
    (NS_FIRST_ARG const char *pszFilename, UINT32 *hFile);



///////////////////////////////////////////////////////////////////////////////////////////////////
// ns_GetFileInfo
//  
// Purpose:
//  Retrieve general information about the data file
//
// Parameters:
//  UINT32 hFile                       handle to NS data file
//  ns_FILEINFO *pFileInfo             pointer to ns_FILEINFO structure that receives data
//  UINT32 dwFileInfoSize              number of bytes returned in ns_FILEINFO
//
// Return Values:
//  ns_OK                              function succeeded
//  ns_BADFILE                         invalid file handle
//  ns_FILEERROR                       file access or read error
//  ns_LIBERROR                        library error, null pointer
//
///////////////////////////////////////////////////////////////////////////////////////////////////

ns_RESULT ns_stdcall ns_GetFileInfo
    (NS_FIRST_ARG UINT32 hFile, ns_FILEINFO *pFileInfo, UINT32 dwFileInfoSize);



///////////////////////////////////////////////////////////////////////////////////////////////////
// ns_CloseFile
//  
// Purpose:
//  Close the open data file
//
// Parameters:
//  UINT32 hFile                       handle to NS data file
//
// Return Values:
//  ns_OK                              function succeeded
//  ns_BADFILE                         invalid file handle
//    
///////////////////////////////////////////////////////////////////////////////////////////////////

ns_RESULT ns_stdcall ns_CloseFile
    (NS_FIRST_ARG UINT32 hFile);



///////////////////////////////////////////////////////////////////////////////////////////////////
// ns_GetEntityInfo
//  
// Purpose:
//  Retrieve Entity information
//
// Parameters:
//  UINT32 hFile                       handle to NS data file
//  UINT32 dwEntityID                  entity ID
//  ns_ENTITYINFO *pEntityInfo         pointer to ns_ENTITYINFO structure that receives information 
//  UINT32 dwEntityInfoSize            number of bytes returned in ns_ENTITYINFO
//
// Return Values:
//  ns_OK                              function succeeded
//  ns_BADFILE                         invalid file handle
//  ns_LIBERROR                        library error, null pointer
//
///////////////////////////////////////////////////////////////////////////////////////////////////

ns_RESULT ns_stdcall ns_GetEntityInfo
    (NS_FIRST_ARG UINT32 hFile, UINT32 dwEntityID, ns_ENTITYINFO *pEntityInfo, UINT32 dwEntityInfoSize);



///////////////////////////////////////////////////////////////////////////////////////////////////
// ns_GetEventInfo
//  
// Purpose:
//  Retrieve information for Event entities.
//
// Parameters:
//  UINT32 hFile                       handle to NS data file
//  UINT32 dwEntityID                  Event entity ID
//  ns_EVENTINFO *pEventInfo           pointer to ns_EVENTINFO structure to receive information 
//  UINT32 dwEventInfoSize             number of bytes returned in ns_EVENTINFO
//
// Return Values:
//  ns_OK                function succeeded
//  ns_BADFILE            invalid file handle
//  ns_BADENTITY            inappropriate or invalid entity identifier
//  ns_LIBERROR            library error, null pointer
//
///////////////////////////////////////////////////////////////////////////////////////////////////

ns_RESULT ns_stdcall ns_GetEventInfo
    (NS_FIRST_ARG UINT32 hFile, UINT32 dwEntityID, ns_EVENTINFO *pEventInfo, UINT32 dwEventInfoSize);



///////////////////////////////////////////////////////////////////////////////////////////////////
// ns_GetEventData
//  
// Purpose:
//  Retrieve the timestamp and Event entity data items.
//
// Parameters:
//  UINT32 hFile                       handle to NS data file
//  UINT32 dwEntityID                  Event entity ID
//  UINT32 nIndex                      Event entity item number
//  double *pdTimeStamp                pointer to double timestamp (in seconds)
//  void   *pData                      pointer to data buffer to receive data
//  UINT32 *pdwDataSize                pointer to number of bytes of data retrieved into data buffer
//
// Return Values:
//  ns_OK                              function succeeded
//  ns_BADFILE                         invalid file handle
//  ns_BADENTITY                       inappropriate or invalie entity identifier
//  ns_LIBERROR                        library error, null pointer
//
///////////////////////////////////////////////////////////////////////////////////////////////////

ns_RESULT ns_stdcall ns_GetEventData
    (NS_FIRST_ARG UINT32 hFile, UINT32 dwEntityID, UINT32 nIndex, double *pdTimeStamp, void *pData,
     UINT32 dwDataSize, UINT32 *pdwDataRetSize);


///////////////////////////////////////////////////////////////////////////////////////////////////
// ns_GetAnalogInfo
//
// Purpose:
//  Retrieve information for Analog entities
//
// Parameters:
//
//  UINT32 hFile                       handle to NS data file
//  UINT32 dwEntityID                  Analog entity ID
//  ns_ANALOGINFO *pAnalogInfo         pointer to ns_ANALOGINFO structure to receive data 
//  UINT32 dwAnalogInfoSize            number of bytes returned in ns_ANALOGINFO
//
// Return Values:
//  ns_OK                              function succeeded
//  ns_BADFILE                         invalid file handle
//  ns_BADENTITY                       inappropriate or invalie entity identifier
//  ns_LIBERROR                        library error, null pointer
//
///////////////////////////////////////////////////////////////////////////////////////////////////

ns_RESULT ns_stdcall ns_GetAnalogInfo
    (NS_FIRST_ARG UINT32 hFile, UINT32 dwEntityID, ns_ANALOGINFO *pAnalogInfo, UINT32 dwAnalogInfoSize);



///////////////////////////////////////////////////////////////////////////////////////////////////
// ns_GetAnalogData
//  
// Purpose:
//  Retrieve analog data in the buffer at pData.  If possible, dwIndexCount, analog data values are
//  returned in the buffer.  As there may be time gaps in the sequential values, the number of
//  continuously sampled data items is returned in pdwContCount.
//
// Parameters:
//  UINT32 hFile                      handle to NS data file
//  UINT32 dwEntityID                 Analog entity ID
//  UINT32 dwStartIndex               starting index to search for timestamp
//  UINT32 dwIndexCount               number of timestamps to retrieve
//  UINT32 *pdwContCount              pointer to count of the first non-sequential analog item
//  double *pData                     pointer of data buffer to receive data values
//
// Return Values:
//  ns_OK                             function succeeded
//  ns_BADFILE                        invalid file handle
//  ns_BADENTITY                      inappropriate or invalid entity identifier
//  ns_LIBERROR                       library error, null pointer
//
///////////////////////////////////////////////////////////////////////////////////////////////////

ns_RESULT ns_stdcall ns_GetAnalogData
    (NS_FIRST_ARG UINT32 hFile, UINT32 dwEntityID, UINT32 dwStartIndex, UINT32 dwIndexCount, 
     UINT32 *pdwContCount, double *pData);



///////////////////////////////////////////////////////////////////////////////////////////////////
// ns_GetSegmentInfo
//  
// Purpose:
//  Retrieve information for Segment entities.
//
// Parameters:
//  UINT32 hFile                       handle to NS data file
//  UINT32 dwEntityID                  Segment entity ID
//  ns_SEGMENTINFO *pSegmentInfo       pointer to ns_SEGMENTINFO structure to receive information
//  UINT32 dwSegmentInfoSize           size in bytes retrieved in ns_SEGMENTINFO structure
//
// Return Values:
//  ns_OK                              function succeeded
//  ns_BADFILE                         invalid file handle
//  ns_BADENTITY                       invalid or inappropriate entity identifier specified
//  ns_FILEERROR                       file access or read error
//  ns_LIBERROR                        library error
//
///////////////////////////////////////////////////////////////////////////////////////////////////

ns_RESULT ns_stdcall ns_GetSegmentInfo
    (NS_FIRST_ARG UINT32 hFile, UINT32 dwEntityID, ns_SEGMENTINFO *pSegmentInfo, UINT32 dwSegmentInfoSize);


///////////////////////////////////////////////////////////////////////////////////////////////////
// ns_GetSegmentSourceInfo
//  
// Purpose:
//  Retrieve information on the source, dwSourceID, generating segment entity dwEntityID.
//
// Parameters:
//  UINT32 hFile                       handle to NS data file
//  UINT32 dwEntityID                  Segment entity ID
//  UINT32 dwSourceID                  entity ID of source
//  ns_SEGSOURCEINFO *pSourceInfo      pointer to ns_SEGSOURCEINFO structure to receive information
//  UINT32 dwSourceInfoSize            size in bytes retrieved in ns_SEGSOURCEINFO structure
//
// Return Values:
//  ns_OK                              function succeeded
//  ns_BADFILE                         invalid file handle
//  ns_BADENTITY                       invalid or inappropriate entity identifier specified
//  ns_FILEERROR                       file access or read error
//  ns_LIBERROR                        library error
//
///////////////////////////////////////////////////////////////////////////////////////////////////

ns_RESULT ns_stdcall ns_GetSegmentSourceInfo
    (NS_FIRST_ARG UINT32 hFile, UINT32 dwEntityID, UINT32 dwSourceID, ns_SEGSOURCEINFO *pSourceInfo,
     UINT32 dwSourceInfoSize);



///////////////////////////////////////////////////////////////////////////////////////////////////
// ns_GetSegmentData
//  
// Purpose:
//  Retrieve segment data waveform and its timestamp.
//  The number of data points read is returned at pdwSampleCount.
//
// Parameters:
//  UINT32 hFile                       handle to NS data file
//  UINT32 dwEntityID                  Segment entity ID
//  INT32 nIndex                       Segment item index to retrieve
//  double *pdTimeStamp                pointer to timestamp to retrieve
//  double *pData                      pointer to data buffer to receive data
//  UINT32 *pdwSampleCount             pointer to number of data items retrieved 
//  UINT32 *pdwUnitID                  pointer to unit ID of Segment data
//
// Return Values:
//  ns_OK                              function succeeded
//  ns_BADFILE                         invalid file handle
//  ns_BADENTITY                       invalid or inappropriate entity identifier specified
//  ns_FILEERROR                       file access or read error
//  ns_LIBERROR                        library error
//
///////////////////////////////////////////////////////////////////////////////////////////////////

ns_RESULT ns_stdcall ns_GetSegmentData
    (NS_FIRST_ARG UINT32 hFile, UINT32 dwEntityID, UINT32 nIndex, double *pdTimeStamp, double *pdData,
     UINT32 dwDataBufferSize, UINT32 *pdwSampleCount, UINT32 *pdwUnitID );



///////////////////////////////////////////////////////////////////////////////////////////////////
// ns_GetNeuralInfo
//  
// Purpose:
//  Retrieve information on Neural Events.
//
// Parameters:
//  UINT32 hFile                       handle to NS data file
//  UINT32 dwEntityID                  Neural entity ID
//  ns_NEURALINFO *pNeuralInfo         pointer to ns_NEURALINFO structure to receive information 
//  UINT32 dwNeuralInfoSize            number of bytes returned in ns_NEURALINFO
//
// Return Values:
//  ns_OK                              function succeeded
//  ns_BADFILE                         invalid file handle
//  ns_BADENTITY                       inappropriate or invalid entity identifier
//  ns_LIBERROR                        library error, null pointer
//
///////////////////////////////////////////////////////////////////////////////////////////////////


ns_RESULT ns_stdcall ns_GetNeuralInfo
    (NS_FIRST_ARG UINT32 hFile, UINT32 dwEntityID, ns_NEURALINFO *pNeuralInfo, UINT32 dwNeuralInfoSize);



///////////////////////////////////////////////////////////////////////////////////////////////////
// ns_GetNeuralData
//  
// Purpose:
//  Retrieve requested number of Neural event timestamps (in sec)
//
// Parameters:
//  UINT32 hFile                       handle to NS data file
//  UINT32 dwEntityID                  Neural event entity ID
//  UINT32 dwStartIndex                index of first Neural event item time to retrieve
//  UINT32 dwIndexCount                number of Neural event items to retrieve
//  double *pData                      pointer to buffer to receive times
//
// Return Values:
//  ns_OK                              function succeeded
//  ns_BADFILE                         invalid file handle
//  ns_BADENTITY                       inappropriate or invalie entity identifier
//  ns_LIBERROR                        library error, null pointer
//
///////////////////////////////////////////////////////////////////////////////////////////////////

ns_RESULT ns_stdcall ns_GetNeuralData
    (NS_FIRST_ARG UINT32 hFile, UINT32 dwEntityID, UINT32 dwStartIndex, UINT32 dwIndexCount, double *pdData);



///////////////////////////////////////////////////////////////////////////////////////////////////
// ns_GetIndexByTime
//  
// Purpose:
//  Given the time (sec), return the closest data item index, as specified by nFlag.
//  Finds the packet with the closest time to the requested time.  
//
// Parameters:
//  UINT32 hFile                       handle to NS data file
//  UINT32 dwEntityID                  entity ID to search for
//  UINT32 dwSearchTimeStamp           timestamp of item to search for
//  INT32 nFlag                        position of item relative to the requested timestamp
//  UINT32 *pdwIndex                   pointer to index of item to retrieve 
//
// Return Values:
//  ns_OK                              function succeeded
//  ns_BADFILE                         invalid file handle
//  ns_LIBERROR                        library error, null pointer
//
///////////////////////////////////////////////////////////////////////////////////////////////////

ns_RESULT ns_stdcall ns_GetIndexByTime
    (NS_FIRST_ARG UINT32 hFile, UINT32 dwEntityID, double dTime, INT32 nFlag, UINT32 *pdwIndex);



///////////////////////////////////////////////////////////////////////////////////////////////////
// ns_GetTimeByIndex
//  
// Purpose:
//  Given an index for an entity data item, return the time in seconds.
//
// Parameters:
//  UINT32 hFile                       handle to NS data file
//  UINT32 dwEntityID                  entity ID to search for
//  UINT32 dwIndex                     index of entity item to search for
//  double *pdTime                     time of entity to retrieve
//
// Return Values:
//  ns_OK                              function succeeded
//  ns_BADFILE                         invalid file handle
//  ns_LIBERROR                        library error, null pointer
//
///////////////////////////////////////////////////////////////////////////////////////////////////

ns_RESULT ns_stdcall ns_GetTimeByIndex
    (NS_FIRST_ARG UINT32 hFile, UINT32 dwEntityID, UINT32 dwIndex, double *pdTime);


///////////////////////////////////////////////////////////////////////////////////////////////////
// ns_GetLastErrorMsg
//  
// Purpose:  
//  Retrieve the most recent error text message
//
// Parameters:
//  char *pszMsgBuffer                 pointer to text buffer to receive error message  
//  UINT32 dwMsgBufferSize             size in bytes of text buffer
//
// Return Values:
//  ns_OK                              function succeeded
//  ns_LIBERROR                        library error
//
///////////////////////////////////////////////////////////////////////////////////////////////////

ns_RESULT ns_stdcall ns_GetLastErrorMsg
    (NS_FIRST_ARG char *pszMsgBuffer, UINT32 dwMsgBufferSize);



/*=========================================================================
| PACKING
 ========================================================================*/
#if defined(__GNUC__)
    #undef GCC_ALIGN_4
#elif defined(_MSC_VER)
    #pragma pack(pop)
#endif



/*=========================================================================
| C++ SUPPORT
 ========================================================================*/
#ifdef __cplusplus
}
#endif


#endif  // include guards
