///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2003-2005 Neuroshare Project
//
// This is free software; you can redistribute it and/or
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
// $Author: kirkkorver $
// $Date: 2005/01/24 05:10:47 $
// $Revision: 1.1 $
// $Source: /cvsroot/neuroshare/NSLibraries/ns/ns.c,v $
//
///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Specification : Glue to a neuroshare library
//
// Authors       : Gopal Santhanam
//
//////////////////////////////////////////////////////////////////////////////


/*=========================================================================
| Project : Neuroshare API Glue
| File    : ns.c
| Version : v1.01
|--------------------------------------------------------------------------
| Portions of this file are copyright (c) 2002-2003 Total Phase, Inc.
| All rights reserved.
| www.totalphase.com
|
| Redistribution and use in source and binary forms, with or without
| modification, are permitted provided that the following conditions
| are met:
|
| - Redistributions of source code must retain the above copyright
|   notice, this list of conditions and the following disclaimer.
|
| - Redistributions in binary form must reproduce the above copyright
|   notice, this list of conditions and the following disclaimer in the
|   documentation and/or other materials provided with the distribution.
|
| - Neither the name of Total Phase, Inc. nor the names of its
|   contributors may be used to endorse or promote products derived from
|   this software without specific prior written permission.
|
| THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
| "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
| LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
| FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE
| COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
| INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
| BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
| LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
| CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
| LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
| ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
| POSSIBILITY OF SUCH DAMAGE.
|--------------------------------------------------------------------------
|
| API definitions for the Neuroshare file format standard.
|
| Use the following files to access Neuroshare files:
|   ns.h        --  API header file
|   ns.c        --  API interface module
|
| A Neuroshare library is required in conjuction with this module. 
| For example:
|   ns_NEV.so   --  Linux shared object
|      or
|   ns_NEV.dll  --  Windows dynamic link library
|
========================================================================*/


/*=========================================================================
| INCLUDES
 ========================================================================*/
/* This #include can be customized to conform to the user's build paths. */
#include "ns.h"


/*=========================================================================
| VERSION CHECK
 ========================================================================*/
#define NS_CFILE_VERSION   0x0103  /* v1.03 */
#define NS_REQ_LIB_VERSION 0x0103  /* v1.03 */

/*
 * Make sure that the ns.h was included and that
 * the version numbers match.
 */
#ifndef NS_HEADER_VERSION
#  error Unable to include ns.h file. Please check include path.

#elif NS_HEADER_VERSION != NS_CFILE_VERSION
#  error Version mismatch between ns.c and ns.h files.

#endif


/*=========================================================================
| CONSTANTS
 ========================================================================*/
#define MAX_SO_PATH 256


/*=========================================================================
| LINUX SUPPORT
 ========================================================================*/
#ifdef linux

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <dlfcn.h>
#include <fcntl.h>

#define DLL_HANDLE  void *
#define SO_EXT      ".so"
#define stricmp              strcasecmp


/*
 * These functions allow the Linux behavior to emulate
 * the Windows behavior as specified below in the Windows
 * support section.
 * 
 * First, search for the shared object in the application
 * binary path, then in the current working directory.
 * 
 * Searching the application binary path requires /proc
 * filesystem support, which is standard in 2.4.x kernels.
 * 
 * If the /proc filesystem is not present, the shared object
 * will not be loaded from the execution path unless that path
 * is either the current working directory or explicitly
 * specified in LD_LIBRARY_PATH.
 */
static int _checkPath (const char *path, char *so_name) {
    char *filename = (char *)malloc(strlen(path) +1 + strlen(so_name) +1);
    int   fd;

    /* Check if the file is readable */
    sprintf(filename, "%s/%s", path, so_name);
    fd = open(filename, O_RDONLY);
    if (fd >= 0) {
        strncpy(so_name, filename, MAX_SO_PATH);
        close(fd);
    }

    /* Clean up and exit */
    free(filename);
    return (fd >= 0);
}

static void _setSearchPath (char *so_name) {
    char  path[MAX_SO_PATH+1];
    int   count;
    char *p;

    /* Make sure that SO_NAME is not an absolute path. */
    if (so_name[0] == '/')  return;

    /* Check the execution directory name. */
    count = readlink("/proc/self/exe", path, MAX_SO_PATH);
    if (count > 0) {
        path[count] = 0;
        if ((p = strrchr(path, '/')) != 0)  *p = '\0';
        if (path[0] == '\0')  strcpy(path, "/");

        /* If there is a match, return immediately. */
        if (_checkPath(path, so_name))  return;
    }

    /* Check the current working directory. */
    p = getcwd(path, MAX_SO_PATH);
    if (p != 0)  _checkPath(path, so_name);
}

#endif


/*=========================================================================
| WINDOWS SUPPORT
 ========================================================================*/
#if defined(WIN32) || defined(_WIN32)

#include <stdio.h>
#include <windows.h>
#include <string.h>

#define DLL_HANDLE           HINSTANCE
#define dlopen(name, flags)  LoadLibrary(name)
#define dlsym(handle, name)  GetProcAddress(handle, name)
#define dlclose(handle)      FreeLibrary(handle);
#define dlerror()            "Exiting program"
#define SO_EXT               ".dll"
#define snprintf             _snprintf
#define stricmp              _stricmp
#ifdef __GNUC__
    #undef stricmp
    #define stricmp          strcasecmp
#endif

/*
 * Use the default Windows DLL loading rules:
 *   1.  The directory from which the application binary was loaded.
 *   2.  The application's current directory.
 *   3a. [Windows NT/2000/XP only] 32-bit system directory
 *       (default: c:\winnt\System32)
 *   3b. 16-bit system directory
 *       (default: c:\winnt\System or c:\windows\system)
 *   4.  The windows directory
 *       (default: c:\winnt or c:\windows)
 *   5.  The directories listed in the PATH environment variable
 */
static void _setSearchPath (char *so_name) 
{
    /* Do nothing */
}

#endif


/*=========================================================================
| LIBRARY FUNCTION TABLE
 ========================================================================*/
#if defined(WIN32) || defined(_WIN32)
#  define ns_api_stdcall __stdcall
#else
#  define ns_api_stdcall
#endif

// Types for all of the function pointers
typedef ns_RESULT (ns_api_stdcall *ns_GetLibraryInfoType)
        (ns_LIBRARYINFO *pLibraryInfo, UINT32 dwLibraryInfoSize);

typedef ns_RESULT (ns_api_stdcall *ns_OpenFileType)
        (const char *pszFilename, UINT32 *hFile);
    
typedef ns_RESULT (ns_api_stdcall *ns_GetFileInfoType)
        (UINT32 hFile, ns_FILEINFO *pFileInfo, UINT32 dwFileInfoSize);
    
typedef ns_RESULT (ns_api_stdcall *ns_CloseFileType)
        (UINT32 hFile);
    
typedef ns_RESULT (ns_api_stdcall *ns_GetEntityInfoType)
        (UINT32 hFile, UINT32 dwEntityID, ns_ENTITYINFO *pEntityInfo, UINT32 dwEntityInfoSize);

    /* Event API */
typedef ns_RESULT (ns_api_stdcall *ns_GetEventInfoType)
        (UINT32 hFile, UINT32 dwEntityID, ns_EVENTINFO *pEventInfo, UINT32 dwEventInfoSize);
    
typedef ns_RESULT (ns_api_stdcall *ns_GetEventDataType)
        (UINT32 hFile, UINT32 dwEntityID, UINT32 nIndex, double *pdTimeStamp, void *pData,
         UINT32 dwDataSize, UINT32 *pdwDataRetSize);



    /* Analog API */
typedef ns_RESULT (ns_api_stdcall *ns_GetAnalogInfoType)
        (UINT32 hFile, UINT32 dwEntityID, ns_ANALOGINFO *pAnalogInfo,UINT32 dwAnalogInfoSize);
    
typedef ns_RESULT (ns_api_stdcall *ns_GetAnalogDataType)
        (UINT32 hFile, UINT32 dwEntityID, UINT32 dwStartIndex, UINT32 dwIndexCount, 
         UINT32 *pdwContCount, double *pData);


    /* Segment API */
typedef ns_RESULT (ns_api_stdcall *ns_GetSegmentInfoType)
        (UINT32 hFile, UINT32 dwEntityID, ns_SEGMENTINFO *pSegmentInfo,UINT32 dwSegmentInfoSize);

typedef ns_RESULT (ns_api_stdcall *ns_GetSegmentSourceInfoType)
        (UINT32 hFile, UINT32 dwEntityID, UINT32 dwSourceID, ns_SEGSOURCEINFO *pSourceInfo,
         UINT32 dwSourceInfoSize);

typedef ns_RESULT (ns_api_stdcall *ns_GetSegmentDataType)
        (UINT32 hFile, UINT32 dwEntityID, INT32 nIndex, double *pdTimeStamp, double *pdData,
         UINT32 dwDataBufferSize, UINT32 *pdwSampleCount, UINT32 *pdwUnitID );


    /* Neural API */
typedef ns_RESULT (ns_api_stdcall *ns_GetNeuralInfoType)
        (UINT32 hFile, UINT32 dwEntityID, ns_NEURALINFO *pNeuralInfo, UINT32 dwNeuralInfoSize);

typedef ns_RESULT (ns_api_stdcall *ns_GetNeuralDataType)
        (UINT32 hFile, UINT32 dwEntityID, UINT32 dwStartIndex, UINT32 dwIndexCount, double *pdData);


    /* More General API */
typedef ns_RESULT (ns_api_stdcall *ns_GetIndexByTimeType)
        (UINT32 hFile, UINT32 dwEntityID, double dTime, INT32 nFlag, UINT32 *pdwIndex);

typedef ns_RESULT (ns_api_stdcall *ns_GetTimeByIndexType)
        (UINT32 hFile, UINT32 dwEntityID, UINT32 dwIndex, double *pdTime);

typedef ns_RESULT (ns_api_stdcall *ns_GetLastErrorMsgType)
        (char *pszMsgBuffer, UINT32 dwMsgBufferSize);



struct NS_LIBRARY_FUNCS 
{
    char so_name[MAX_SO_PATH+1];
    int valid;
    int refCount;
    DLL_HANDLE dllHandle;
    
    /* General API */
    ns_GetLibraryInfoType ns_GetLibraryInfo;
    
    ns_OpenFileType ns_OpenFile;
    ns_GetFileInfoType ns_GetFileInfo;
    ns_CloseFileType ns_CloseFile;
    
    ns_GetEntityInfoType ns_GetEntityInfo;

    /* Event API */
    ns_GetEventInfoType ns_GetEventInfo;
    ns_GetEventDataType ns_GetEventData;

    /* Analog API */
    ns_GetAnalogInfoType ns_GetAnalogInfo;
    ns_GetAnalogDataType ns_GetAnalogData;

    /* Segment API */
    ns_GetSegmentInfoType ns_GetSegmentInfo;
    ns_GetSegmentSourceInfoType ns_GetSegmentSourceInfo;
    ns_GetSegmentDataType ns_GetSegmentData;

    /* Neural API */
    ns_GetNeuralInfoType ns_GetNeuralInfo;
    ns_GetNeuralDataType ns_GetNeuralData;
    
    /* More General API */
    ns_GetIndexByTimeType ns_GetIndexByTime;
    ns_GetTimeByIndexType ns_GetTimeByIndex;
    ns_GetLastErrorMsgType ns_GetLastErrorMsg;
};


#define MAX_LIBS 32

/* This is initialized to zero as per ANSI C specifications */
static struct NS_LIBRARY_FUNCS _libraries[MAX_LIBS];

#define CHECK_VALIDITY(nsDllHandle)                  \
    if (nsDllHandle < 1 || nsDllHandle > MAX_LIBS || \
        _libraries[nsDllHandle-1].valid == 0)        \
        return ns_LIBERROR;


/*=========================================================================
| SHARED LIBRARY LOADER
 ========================================================================*/
/* The error conditions can be customized depending on the application. */
/* a return value of 0 indicates an error */
ns_DLLHANDLE ns_stdcall ns_LoadLibrary (const char *libname)
{
    int i;
    DLL_HANDLE handle;

    ns_LIBRARYINFO lib_info;
    UINT32 lib_api_version;

    char so_name[MAX_SO_PATH+1];
    
    for (i=0; i < MAX_LIBS; ++i) 
    {
        if (_libraries[i].valid && !strcmp(libname, _libraries[i].so_name)) 
        {
            ++_libraries[i].refCount;
            return (i+1);
        }
    }

    /* Find an open slot */
    for (i=0; i < MAX_LIBS; ++i) 
    {
        if (!_libraries[i].valid)
            break;
    }

    /* Return an error if we don't have enough slots left */
    if (i == MAX_LIBS)
        return ns_LIBERROR;

    /* Add the extension if necessary */
    {
        const char * pExt;      /* point to the extension */
        pExt = libname + strlen(libname) - strlen(SO_EXT);
        if (pExt < libname) /* allow for really short library names) */
            pExt = libname;

        /* strcmp ignoring case */
        if (stricmp(pExt, SO_EXT) == 0)
            snprintf(so_name, MAX_SO_PATH, "%s", libname);              /* copy the name */
        else
            snprintf(so_name, MAX_SO_PATH, "%s%s", libname, SO_EXT);    /* copy the name, add extension */

        so_name[MAX_SO_PATH] = 0;   /* force termination */
    }
    
    _setSearchPath(so_name);
    handle = dlopen(so_name, RTLD_LAZY);
    if (handle == 0) 
    {
#if 1
        fprintf(stderr, "ns.c: Unable to load %s\n", so_name);
        fprintf(stderr, "ns.c: %s\n", dlerror());
#endif
        return 0;
    }
    
    _libraries[i].ns_GetLibraryInfo = (ns_GetLibraryInfoType)dlsym(handle,
                                                    "ns_GetLibraryInfo");
    
    if (_libraries[i].ns_GetLibraryInfo == 0) 
    {
#if 1
        fprintf(stderr, "ns.c: Unable to bind %s ns_GetLibraryInfo function\n",
                so_name);
        fprintf(stderr, "ns.c: %s\n", dlerror());
#endif
        dlclose(handle);
        return 0;
    }

    /* test for compatible versions */
    _libraries[i].ns_GetLibraryInfo(&lib_info, sizeof(lib_info));
    
    lib_api_version = ((lib_info.dwAPIVersionMaj << 8) & 0xff00) |
                      ((lib_info.dwAPIVersionMin << 0) & 0x00ff);
    
#if 0       /* ALL VERSIONS ARE COMPATIBLE */
    if (lib_api_version  < NS_REQ_LIB_VERSION ||
        NS_CFILE_VERSION < lib_api_version)
    {
#if 1
        fprintf(stderr, "\nns.c:Incompatible versions:\n");
        
        fprintf(stderr, "ns.c:  NS Glue version     = v%d.%02d  ",
                (NS_CFILE_VERSION >> 8) & 0xff, NS_CFILE_VERSION & 0xff);
        
        if (lib_api_version < NS_REQ_LIB_VERSION)
            fprintf(stderr, "(requires Library API version >= %d.%02d)\n",
                    (NS_REQ_LIB_VERSION >> 8) & 0xff,
                    (NS_REQ_LIB_VERSION >> 0) & 0xff);
        else
            fprintf(stderr, "(Library API version OK)\n");
        
        fprintf(stderr, "ns.c:  Library API version = v%d.%02d  ",
                (lib_api_version >> 8) & 0xff, lib_api_version  & 0xff);

        if (NS_CFILE_VERSION < NS_REQ_LIB_VERSION)
            fprintf(stderr, "(requires NS Glue version >= %d.%02d)\n",
                    (lib_api_version  >> 8) & 0xff,
                    (lib_api_version  >> 0) & 0xff);
        else
            fprintf(stderr, "(NS Glue version OK)\n");
        
#endif
        dlclose(handle);
        return 0;
    }
#endif


// This macro only works because several conventions are followed
// 1) A function named foo, has a type (typedef) called fooType
// 2) The name of the function in our structure is the same name as 3)
// 3) The exported DLL name is the same name as 2)
#define LOAD_LIBRARY(func)                                           \
    {                                                                \
        _libraries[i].func = (func##Type)dlsym(handle, #func);       \
        if (_libraries[i].func == 0) {                               \
            fprintf(stderr, "ns.c: Unable to bind %s %s function\n", \
                    so_name, #func);                                 \
            fprintf(stderr, "ns.c: %s\n", dlerror());                \
            goto error;                                              \
        }                                                            \
    }


    /* More General API */
    LOAD_LIBRARY(ns_OpenFile);
    LOAD_LIBRARY(ns_GetFileInfo);
    LOAD_LIBRARY(ns_CloseFile);
    LOAD_LIBRARY(ns_GetEntityInfo);

    /* Event API */
    LOAD_LIBRARY(ns_GetEventInfo);
    LOAD_LIBRARY(ns_GetEventData);

    /* Analog API */
    LOAD_LIBRARY(ns_GetAnalogInfo);
    LOAD_LIBRARY(ns_GetAnalogData);

    /* Segment API */
    LOAD_LIBRARY(ns_GetSegmentInfo);
    LOAD_LIBRARY(ns_GetSegmentSourceInfo);
    LOAD_LIBRARY(ns_GetSegmentData);

    /* Neural API */
    LOAD_LIBRARY(ns_GetNeuralInfo);
    LOAD_LIBRARY(ns_GetNeuralData);

    /* More General API */
    LOAD_LIBRARY(ns_GetIndexByTime);
    LOAD_LIBRARY(ns_GetTimeByIndex);
    LOAD_LIBRARY(ns_GetLastErrorMsg);

#undef LOAD_LIBRARY

    /* Setup the other variables in the library structure */
    _libraries[i].valid     = 1;
    _libraries[i].dllHandle = handle;
    _libraries[i].refCount  = 1;
    strncpy(_libraries[i].so_name, so_name, MAX_SO_PATH);
    _libraries[i].so_name[MAX_SO_PATH] = 0;
    
    return (i+1);

error:
    dlclose(handle);
    return 0;
}

ns_RESULT ns_stdcall ns_CloseLibrary (ns_DLLHANDLE nsDllHandle)
{
    CHECK_VALIDITY(nsDllHandle)
        
    --_libraries[nsDllHandle -1].refCount;
    if (_libraries[nsDllHandle -1].refCount == 0) {
        dlclose(_libraries[nsDllHandle -1].dllHandle);
        _libraries[nsDllHandle -1].valid = 0;
    }
    
    return ns_OK;
}

            
/*=========================================================================
| GENERAL API
 ========================================================================*/
ns_RESULT ns_stdcall ns_GetLibraryInfo (ns_DLLHANDLE nsDllHandle, ns_LIBRARYINFO *pLibraryInfo, UINT32 dwLibraryInfoSize)
{
    CHECK_VALIDITY(nsDllHandle)
    return _libraries[nsDllHandle-1].ns_GetLibraryInfo(pLibraryInfo, dwLibraryInfoSize);
}
    
ns_RESULT ns_stdcall ns_OpenFile (ns_DLLHANDLE nsDllHandle, const char *pszFilename, UINT32 *hFile)
{
    CHECK_VALIDITY(nsDllHandle)
    return _libraries[nsDllHandle-1].ns_OpenFile(pszFilename, hFile);
}
    
ns_RESULT ns_stdcall ns_GetFileInfo (ns_DLLHANDLE nsDllHandle, UINT32 hFile, ns_FILEINFO *pFileInfo, UINT32 dwFileInfoSize)
{
    CHECK_VALIDITY(nsDllHandle)
    return _libraries[nsDllHandle-1].ns_GetFileInfo(hFile, pFileInfo, dwFileInfoSize);
}
    
ns_RESULT ns_stdcall ns_CloseFile (ns_DLLHANDLE nsDllHandle, UINT32 hFile)
{
    CHECK_VALIDITY(nsDllHandle)
    return _libraries[nsDllHandle-1].ns_CloseFile(hFile);
}
    
ns_RESULT ns_stdcall ns_GetEntityInfo (ns_DLLHANDLE nsDllHandle, UINT32 hFile, UINT32 dwEntityID, ns_ENTITYINFO *pEntityInfo, UINT32 dwEntityInfoSize)
{
    CHECK_VALIDITY(nsDllHandle)
    return _libraries[nsDllHandle-1].ns_GetEntityInfo(hFile, dwEntityID, pEntityInfo, dwEntityInfoSize);
}


/*=========================================================================
| EVENT API
 ========================================================================*/
ns_RESULT ns_stdcall ns_GetEventInfo (ns_DLLHANDLE nsDllHandle, UINT32 hFile, UINT32 dwEntityID, ns_EVENTINFO *pEventInfo, UINT32 dwEventInfoSize)
{
    CHECK_VALIDITY(nsDllHandle)
    return _libraries[nsDllHandle-1].ns_GetEventInfo(hFile, dwEntityID, pEventInfo, dwEventInfoSize);
}
    
ns_RESULT ns_stdcall ns_GetEventData (ns_DLLHANDLE nsDllHandle, UINT32 hFile, UINT32 dwEntityID, UINT32 nIndex, double *pdTimeStamp, void *pData,
                           UINT32 dwDataSize, UINT32 *pdwDataRetSize)
{
    CHECK_VALIDITY(nsDllHandle)
    return _libraries[nsDllHandle-1].ns_GetEventData(hFile, dwEntityID, nIndex, pdTimeStamp, pData, dwDataSize, pdwDataRetSize);
}


/*=========================================================================
| ANALOG API
 ========================================================================*/
ns_RESULT ns_stdcall ns_GetAnalogInfo (ns_DLLHANDLE nsDllHandle, UINT32 hFile, UINT32 dwEntityID, ns_ANALOGINFO *pAnalogInfo, UINT32 dwAnalogInfoSize)
{
    CHECK_VALIDITY(nsDllHandle)
    return _libraries[nsDllHandle-1].ns_GetAnalogInfo(hFile, dwEntityID, pAnalogInfo, dwAnalogInfoSize);
}
    
ns_RESULT ns_stdcall ns_GetAnalogData (ns_DLLHANDLE nsDllHandle, UINT32 hFile, UINT32 dwEntityID, UINT32 dwStartIndex, UINT32 dwIndexCount, 
                            UINT32 *pdwContCount, double *pData)
{
    CHECK_VALIDITY(nsDllHandle)
    return _libraries[nsDllHandle-1].ns_GetAnalogData(hFile, dwEntityID, dwStartIndex, dwIndexCount, pdwContCount, pData);
}


/*=========================================================================
| SEGMENT API
 ========================================================================*/
ns_RESULT ns_stdcall ns_GetSegmentInfo (ns_DLLHANDLE nsDllHandle, UINT32 hFile, UINT32 dwEntityID, ns_SEGMENTINFO *pSegmentInfo, UINT32 dwSegmentInfoSize)
{
    CHECK_VALIDITY(nsDllHandle)
    return _libraries[nsDllHandle-1].ns_GetSegmentInfo(hFile, dwEntityID, pSegmentInfo, dwSegmentInfoSize);
}

ns_RESULT ns_stdcall ns_GetSegmentSourceInfo (ns_DLLHANDLE nsDllHandle, UINT32 hFile, UINT32 dwEntityID, UINT32 dwSourceID, ns_SEGSOURCEINFO *pSourceInfo,
                                   UINT32 dwSourceInfoSize)
{
    CHECK_VALIDITY(nsDllHandle)
    return _libraries[nsDllHandle-1].ns_GetSegmentSourceInfo(hFile, dwEntityID, dwSourceID, pSourceInfo, dwSourceInfoSize);
}

ns_RESULT ns_stdcall ns_GetSegmentData (ns_DLLHANDLE nsDllHandle, UINT32 hFile, UINT32 dwEntityID, UINT32 nIndex, double *pdTimeStamp, double *pdData,
                             UINT32 dwDataBufferSize, UINT32 *pdwSampleCount, UINT32 *pdwUnitID )
{
    CHECK_VALIDITY(nsDllHandle)
    return _libraries[nsDllHandle-1].ns_GetSegmentData(hFile, dwEntityID, nIndex, pdTimeStamp, pdData,
                                      dwDataBufferSize, pdwSampleCount, pdwUnitID );
}


/*=========================================================================
| NEURAL API
 ========================================================================*/
ns_RESULT ns_stdcall ns_GetNeuralInfo (ns_DLLHANDLE nsDllHandle, UINT32 hFile, UINT32 dwEntityID, ns_NEURALINFO *pNeuralInfo, UINT32 dwNeuralInfoSize)
{
    CHECK_VALIDITY(nsDllHandle)
    return _libraries[nsDllHandle-1].ns_GetNeuralInfo(hFile, dwEntityID, pNeuralInfo, dwNeuralInfoSize);
}

ns_RESULT ns_stdcall ns_GetNeuralData (ns_DLLHANDLE nsDllHandle, UINT32 hFile, UINT32 dwEntityID, UINT32 dwStartIndex, UINT32 dwIndexCount, double *pdData)
{
    CHECK_VALIDITY(nsDllHandle)
    return _libraries[nsDllHandle-1].ns_GetNeuralData(hFile, dwEntityID, dwStartIndex, dwIndexCount, pdData);
}

    
/*=========================================================================
| MORE GENERAL API
 ========================================================================*/
ns_RESULT ns_stdcall ns_GetIndexByTime (ns_DLLHANDLE nsDllHandle, UINT32 hFile, UINT32 dwEntityID, double dTime, INT32 nFlag, UINT32 *pdwIndex)
{
    CHECK_VALIDITY(nsDllHandle)
    return _libraries[nsDllHandle-1].ns_GetIndexByTime(hFile, dwEntityID, dTime, nFlag, pdwIndex);
}

ns_RESULT ns_stdcall ns_GetTimeByIndex (ns_DLLHANDLE nsDllHandle, UINT32 hFile, UINT32 dwEntityID, UINT32 dwIndex, double *pdTime)
{
    CHECK_VALIDITY(nsDllHandle)
    return _libraries[nsDllHandle-1].ns_GetTimeByIndex(hFile, dwEntityID, dwIndex, pdTime);
}

ns_RESULT ns_stdcall ns_GetLastErrorMsg (ns_DLLHANDLE nsDllHandle, char *pszMsgBuffer, UINT32 dwMsgBufferSize)
{
    CHECK_VALIDITY(nsDllHandle)
    return _libraries[nsDllHandle-1].ns_GetLastErrorMsg(pszMsgBuffer, dwMsgBufferSize);
}
