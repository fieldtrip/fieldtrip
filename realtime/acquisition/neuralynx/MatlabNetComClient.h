#ifndef __MatlabNetComClient_h
#define __MatlabNetComClient_h

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

//client functions
int ConnectToServer(char* serverName);
int DisconnectFromServer();
int OpenStream(char* cheetahObjectName);
int CloseStream(char* cheetahObjectName);
int SendCommand(char* command, char* &reply);

//setter functions
int SetApplicationName(char* myApplicationName);
//extern "C" __declspec(dllexport) bool SetLogFileName(string filename);

//getter functions
int GetCheetahObjectsAndTypes(char** cheetahObjects, char** cheetahTypes);

//getter status functions
char* GetServerPCName();
char* GetServerIPAddress();
char* GetServerApplicationName();
int AreWeConnected();

//setup functions for data retrieval
int GetRecordBufferSize(void);
bool SetRecordBufferSize(int numRecordsToBuffer);
int GetMaxCSCSamples(void);
int GetSpikeSampleWindowSize(void);
int GetMaxSpikeFeatures(void);
int GetMaxEventStringLength(void);


//data retrieval functions
int GetNewCSCData(char* acqEntName, uint64_t **timeStamps, int **channelNumbers, int **samplingFrequency, int **numValidSamples, short **samples, int *numRecordsReturned, int *numDroppedRecords);
int GetNewEventData(char* acqEntName, uint64_t **timeStamps, int **eventIDs, int **ttlValues, char** eventStrings, int *numRecordsReturned, int *numDroppedRecords);
/*
bool GetNewSEData(char* acqEntName, __int64* &timeStamps, int* &scNumbers, int* &cellNumbers, int* &featureValues, short* &samples, int &numRecordsReturned, int &numDroppedRecords); 
bool GetNewSTData(char* acqEntName, __int64* &timeStamps, int* &scNumbers, int* &cellNumbers, int* &featureValues, short* &samples, int &numRecordsReturned, int &numDroppedRecords); 
bool GetNewTTData(char* acqEntName, __int64* &timeStamps, int* &scNumbers, int* &cellNumbers, int* &featureValues, short* &samples, int &numRecordsReturned, int &numDroppedRecords); 
bool GetNewVTData(char* acqEntName, __int64* &timeStamps, int* &extractedLocations, int* &extractedAngles, int &numRecordsReturned, int &numDroppedRecords);
*/

#ifdef __cplusplus
}
#endif

#endif 
