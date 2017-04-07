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
#include "buffer.h"
#include "pthread.h"
#include <math.h>

using namespace std;

#define MAX_DEVICE			1	//Max number of devices supported by this demo
#define USE_MASTER_SLAVE	FALSE


/*
*	SignalInfo
*	function is used to get information about available signals of a device
*/


ULONG SignalInfo(RTDevice *Device)
{	BOOLEAN Allocate = FALSE;
	ULONG BytesPerSample = 0 ;
	ULONG i,NumberOfChannels = 0;


	PSIGNAL_FORMAT psf;

	if( Allocate )
	{
		SIGNAL_FORMAT sf;
		memset( &sf ,0, sizeof( SIGNAL_FORMAT ) );
		sf.Elements = 1;
		sf.Size = sizeof( SIGNAL_FORMAT );


		if( Device->GetSignalFormat( &sf ) != NULL )
		{	ULONG TotalSize = sf.Size * sf.Elements;

			psf = (SIGNAL_FORMAT *) LocalAlloc( LMEM_FIXED | LMEM_ZEROINIT , TotalSize );
			if( psf == NULL ) return 0;

			psf[0].Size = sf.Size;
			psf[0].Elements = sf.Elements;

			Device->GetSignalFormat( psf );
		}
	}
	else
		psf = Device->GetSignalFormat( NULL );

	if( psf != NULL )
	{
		UINT Size = LocalSize( psf );

		if( Size < sizeof( SIGNAL_FORMAT ) * psf->Elements )
			return 0;

		NumberOfChannels = psf->Elements;
		wprintf(L"\n");
		for( i = 0 ; i < NumberOfChannels ; i++ )
		{	BytesPerSample += psf[i].Bytes;
			wprintf(L"\nChannel %3d: %s on device %s:%d",i+1,psf[i].Name,psf[i].PortName,psf[i].SerialNumber );
		}
		//Remove the data from memory
		if( Allocate ) LocalFree( psf );
		else Device->Free( psf );

		wprintf(L"\n");
	}

	return BytesPerSample;
}




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

		_tprintf( _T("%d . %s %d\n"), x + 1, DeviceName , SerialNumber	);


		if( x!= 0 )
		{	HANDLE SlaveHandle = Devices[x]->GetSlaveHandle();
			if(SlaveHandle == 0 )
				_tprintf( _T("Unable to get a handle from device %d\n"), x + 1 );
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
				_tprintf( "%d . %s %d\n" , Count , DeviceName , SerialNumber	);
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



RTDeviceEx *InitDevice( ULONG SampRate)
// Get handle to my device and init the device
{	ULONG Index;
	ULONG NrOfSamples=0,Total=0;


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
//function to test if file exist
bool fexists(const char *filename)
{
  ifstream ifile(filename);
  return (!ifile.bad());
}


int main(int argc, char* argv[])
{
	ULONG SampleRate;
	int BufSizeFactor=5;
	int SaveData;
	ULONG BufferSize;


	ULONG NrOfSamples=0,Total=0;
	ULONG PercentFull,Overflow;
	ULONG BytesPerSample=72;
	ULONG BytesReturned;
	ULONG TotalNrChannelsInDevice;

// Buffer for storing the samples;
	ULONG SignalBuffer[1000];
	SampleRate = MAX_SAMPLE_RATE;
	BufferSize = MAX_BUFFER_SIZE;



	SIGNAL_FORMAT sf;
    PSIGNAL_FORMAT psf;
	memset( &sf ,0, sizeof( SIGNAL_FORMAT ) );
		sf.Elements = 1;
		sf.Size = sizeof( SIGNAL_FORMAT );


	RTDeviceEx *Master;

    Master=InitDevice(SampleRate );

	Master->SetSignalBuffer(&SampleRate,&BufferSize);

	ULONG TotalSize = sf.Size * sf.Elements;




	psf = (SIGNAL_FORMAT *) LocalAlloc( LMEM_FIXED | LMEM_ZEROINIT , TotalSize );
			if( psf == NULL ) return 0;

			psf[0].Size = sf.Size;
			psf[0].Elements = sf.Elements;

	Master->GetSignalFormat( psf );

	TotalNrChannelsInDevice=psf[0].Elements;
	//Remove the data from memory
	 LocalFree( psf );

// Get the size of each sample
//BytesPerSample = SignalInfo(Master);
//Leuk dit lijkt niet te kloppen volgens inf en documentatie geeft de refa inderdaad 220 bytes per samples af
// helaas wordt in de buffer vier bytes per channel * 74 =296 GEBRUIKT
// dus bij het lezen van de Buffer dit getal bebruiken
	BytesPerSample=4*(TotalNrChannelsInDevice);

	if( BytesPerSample == 0 )
	{	printf( "\nDevice returns no samples" );
		_getch();
		return 0;
	}

	wprintf(L"\nMaximum sample rate = %d Hz",SampleRate  / 1000 );
	wprintf(L"\nMaximum Buffer size = %d Samples",BufferSize);

	// Read sample freq and  blocksize factor from parameter file
	char dummy1[80],dummy2[80],dummy3[80];
	FILE *pFile;

	if( fexists("parameter.txt"))
	{


	pFile=fopen("parameter.txt","r");
	  fscanf (pFile, "%s %s %s", dummy1,dummy2,dummy3);
	  fscanf (pFile, "%u %u %u", &SampleRate,&BufSizeFactor,&SaveData);
      fclose (pFile);
	  printf("\n %s %s %s \n ", dummy1,dummy2,dummy3);
      printf("\n %u %u %u \n ",SampleRate,BufSizeFactor,SaveData);

	}
	//create datafile to store data from the aquisition device
	FILE *pOutputFile;
    char maand[2],dag[2],uur[4],minuut[2],seconde[4];
	char filename[22];
	time_t rawtime;
	struct tm * timeinfo;
    if (SaveData == 1)
	{
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
    itoa((timeinfo->tm_mon+1),maand,10);
    itoa((timeinfo->tm_mday),dag,10);
    itoa((timeinfo->tm_hour),uur,10);
    itoa((timeinfo->tm_min),minuut,10);
    itoa((timeinfo->tm_sec),seconde,10);

	strcpy(filename,"TmsiData");
    strcat(filename,maand);
	strcat(filename,dag);
    strcat(filename,uur);
	strcat(filename,minuut);
	strcat(filename,seconde);
	strcat(filename,".txt");
	pOutputFile=fopen(filename,"w");
	}
	//fclose (pOutputFile);

	//Set sample rate
	//SampleRate = 500000 ;
	//Set buffer size;
	BufferSize = 1000;

	Master->SetSignalBuffer(&SampleRate ,&BufferSize);

	wprintf(L"\nSelected sample rate = %d Hz",SampleRate  / 1000 );
	wprintf(L"\nSelected Buffer size = %d Samples",BufferSize);

	host_t host;
    int rc;
    pthread_t tid;
    int connection;

	check_datatypes();
	sprintf(host.name, DEFAULT_HOSTNAME);
	host.port = DEFAULT_PORT;

	/* start the buffer in a seperate thread */
	rc = pthread_create(&tid, NULL, tcpserver, (void *)(&host));
	if (rc) {
		fprintf(stderr, "main err1: return code from pthread_create() is %d\n", rc);
		exit(-1);
	}
    pthread_detach(tid);
	usleep(1000000);

    int i=0, j=0,k=0,sample = 0, status = 0, verbose = 0,datablockready=0;
	long tdif;


	/* these represent the acquisition system properties */




	int nchans         = TotalNrChannelsInDevice+1;//one channel extra for the sample nr
	int fsample        = SampleRate/1000;
	int buffersamp	   = 0;
	int Totalsamp	   = 0;
	int blocksize      = fsample/BufSizeFactor; // nr of samples that is written simultaniously to the buffer

	/* these are used in the communication and represent statefull information */
	int server             = -1;
	message_t    *request  = NULL;
	message_t    *response = NULL;
	header_t     *header   = NULL;
	data_t       *data     = NULL;
	event_t      *event    = NULL;

    if (verbose>0) fprintf(stderr, "Tmsidata host.name =  %s\n", host.name);
	if (verbose>0) fprintf(stderr, "Tmsidata host.port =  %d\n", host.port);

	/* allocate the elements that will be used in the communication */
	request      = static_cast<message_t*> (malloc(sizeof(message_t)));
	request->def = static_cast<messagedef_t*>(malloc(sizeof(messagedef_t)));
	request->buf = NULL;
	request->def->version = VERSION;
	request->def->bufsize = 0;

	header      = static_cast<header_t*>(malloc(sizeof(header_t)));
	header->def = static_cast<headerdef_t*>(malloc(sizeof(headerdef_t)));
	header->buf = NULL;
    data      = static_cast<data_t*>(malloc(sizeof(data_t)));
	data->def = static_cast<datadef_t*>(malloc(sizeof(datadef_t)));
	data->buf = NULL;

	event      = static_cast<event_t*>(malloc(sizeof(event_t)));
	event->def = static_cast<eventdef_t*>(malloc(sizeof(eventdef_t)));
	event->buf = NULL;

	/* define the header */
	header->def->nchans    = nchans;
	header->def->nsamples  = 0;
	header->def->nevents   = 0;
	header->def->fsample   = fsample;
	header->def->data_type = DATATYPE_INT32;
	header->def->bufsize   = 0;
	FREE(header->buf);

	/* define the constant part of the data and allocate space for the variable part */
	data->def->nchans    = nchans;
	data->def->nsamples  = blocksize;
	data->def->data_type = DATATYPE_INT32;
	data->def->bufsize   = WORDSIZE_INT32*nchans*blocksize;
	FREE(data->buf);
	data->buf            = malloc(WORDSIZE_INT32*nchans*blocksize);

	/* initialization phase, send the header */
	request->def->command = PUT_HDR;
	request->def->bufsize = append(&request->buf, request->def->bufsize, header->def, sizeof(headerdef_t));
	request->def->bufsize = append(&request->buf, request->def->bufsize, header->buf, header->def->bufsize);

    server = open_connection(host.name, host.port);
	status = clientrequest(server, request, &response);
    if (verbose>0) fprintf(stderr, "tmsidata: clientrequest returned %d\n", status);
    if (status) {
      fprintf(stderr, "tmsidata: err1\n");
      goto cleanup;
    }
	//status = close_connection(server);
    if (status) {
      fprintf(stderr, "tmsidata: err2\n");
      goto cleanup;
    }

	/* FIXME do someting with the response, i.e. check that it is OK */
	cleanup_message(reinterpret_cast<void**>(&request));
	cleanup_message(reinterpret_cast<void**>(&response));
    request = NULL;
    response = NULL;

	// Some example routines for testing device features

	/********** MEASURING MODE ****************/

	// Measuring mode can switch the device to
	// Calbartion mode
	// Impedance mode
	// and back to normal measurment

	/*
	#define MEASURE_MODE_NORMAL			((ULONG)0x0)
	#define MEASURE_MODE_IMPEDANCE		((ULONG)0x1)
	#define MEASURE_MODE_CALIBRATION	((ULONG)0x2)
	#define MEASURE_MODE_IMPEDANCE_EX	((ULONG)0x3)
	#define MEASURE_MODE_CALIBRATION_EX	((ULONG)0x4)


	if( !Master->MeasuringMode( MEASURE_MODE_IMPEDANCE_EX  ,  1   ) )
	{
		wprintf(L"\nCalibartion and Impedance measurement not supported by this device\n");
	}

	return 0;
	*/


	/************ RTC *****************/
	//Some portable devices has an internal clock

	/*
	SYSTEMTIME time ;

	GetLocalTime( &time );

	Master->SetRtcAlarmTime( &time );
	Master->SetRtcTime( &time );
	if( Master->GetRtcTime( &time ) )
	{	wprintf(L"\nDevice has a realtime clock ,");
		wprintf(L"Date = %d-%d-%d , Time = %d:%d:%d " ,
				time.wDay ,time.wMonth , time.wYear ,
				time.wHour , time.wMinute ,time.wSecond );
	}

	*/
	/********** OFFSET **********/
	/*
	float Offset[2] = { 1000, 1000 };
	if( !Master->SetOffset( Offset, sizeof(Offset) / sizeof(float)  ) )
	{	ShowError("Unable to set offset");
	}
	*/

	/********** GAIN **********/
	/*
	float Gain[2] = { 1 , 1 };
	if( !Master->SetGain( Gain , sizeof(Gain) / sizeof(float)  ) )
	{	ShowError("Unable to set gain");
	}
	*/
	/********** HIGHPASS **********/
	/*
	float HighPass[2];
	HighPass[0] = (float)10.0 * ( SampleRate / 1000 ) ; //10 second highpass filter for channel 0
	HighPass[1] = (float)0.1 * ( SampleRate / 1000 ) ;  //0.1 second highpass filter for channel 1
	if( !Master->SetHighPass( HighPass , sizeof(HighPass) / sizeof(float)  ) )
	{	ShowError("Unable to set highpass filter");
	}
	/********** LOWPASS **********/
	/*
	float LowPass[2];
	LowPass[0] = (float)10.0 * ( SampleRate / 1000 ) ;
	LowPass[1] = (float)0.1 * ( SampleRate / 1000 ) ;
	if( !Master->SetLowPass( LowPass , sizeof(LowPass) / sizeof(float)  ) )
	{	ShowError("Unable to set lowpass filter");
	}
	*/

//Set sample rate


	if(Master->Start() )
	{
		//ULONG ShowChannel = 1;
		wprintf(L"\nPress any key to quit \n");

		// Stop the program if we hit any key
		i=0; j=0; char intstr[5];
		while(!_kbhit())
		{	//beginwhile
			//Get Signal buffer information
			Master->GetBufferInfo(&Overflow,&PercentFull);

			if( PercentFull > 0)
			{	//beginPF
				// if there is data available get samples from the device
				// GetSamples returns the number of bytes written in the signal buffer
				// This will always be a multiple op BytesPerSample.

				// Divide the result by BytesPerSamples to get the number of samples returned


				BytesReturned = Master->GetSamples((PULONG)SignalBuffer,sizeof(SignalBuffer));

				if( BytesReturned != 0)
				{	Total += BytesReturned;
	//				wprintf(L"\rSampleCounter = %8d, ,Sampwritten=%8d,Samp[%d]=%d, %d ,Buffer  = %d, Overflow = %d      " , Total /BytesPerSample ,sample, ShowChannel,SignalBuffer[ShowChannel], SignalBuffer[ShowChannel + 1],PercentFull,Overflow);
				buffersamp=(BytesReturned/BytesPerSample);
				Totalsamp+=buffersamp;
				/*fill the data block */
				k=0;
				while (k<buffersamp)
				{
					while ((i<blocksize)&&(k<buffersamp))
						{
						for (j=0; j<nchans;j++)
						{
						if (j<1)
							{
							((INT32_T *)(data->buf))[i*nchans+0]=sample;
							if (SaveData == 1){
							fprintf(pOutputFile,"%12ld",sample);
							}//end if savedata
							}
						else
							{
							((INT32_T *)(data->buf))[(i*nchans)+j]=SignalBuffer[k*(nchans-1)+j-1];
							if (SaveData == 1){
							fprintf(pOutputFile,"%12ld",SignalBuffer[k*(nchans-1)+j-1]);
							}//end if savedata
							}
						}//end for
						if (SaveData == 1){
						fprintf(pOutputFile,"\n");
						}
						k++;i++;sample++;
						}//end while blocksize


				if (i==blocksize)
				{




		/* create the request */
		request      =static_cast<message_t*>(malloc(sizeof(message_t)));
		request->def =static_cast<messagedef_t*>(malloc(sizeof(messagedef_t)));
		request->buf = NULL;
		request->def->version = VERSION;
		request->def->bufsize = 0;
		request->def->command = PUT_DAT;
		request->def->bufsize = append(&request->buf, request->def->bufsize, data->def, sizeof(datadef_t));
		request->def->bufsize = append(&request->buf, request->def->bufsize, data->buf, data->def->bufsize);

	//	server = open_connection(host.name, host.port);
		status = clientrequest(server, request, &response);
        if (verbose>0) fprintf(stderr, "tmsidata: clientrequest returned %d\n", status);
    	if (status) {
			fprintf(stderr, "tmsidata: err3\n");
      		goto cleanup;
		}
	//	status = close_connection(server);
    	if (status) {
      		fprintf(stderr, "tmsidata: err4\n");
     		goto cleanup;
		}

		/* FIXME do someting with the response, i.e. check that it is OK */
		cleanup_message(reinterpret_cast<void**>(&request));
		cleanup_message(reinterpret_cast<void**>(&response));
    	request = NULL;
    	response = NULL;

		/* approximate delay in microseconds */
		//tdif = (long)(blocksize * 1000000 / fsample);
		//usleep(tdif);
		i=0;

		}//enddatablockready
		}//endK

		}//endBytesreturned
			else
			{	//allow other applications some extra process time
				Sleep(10);
			}
		}//end PercentFull


		}//endwhile

		//Stop the device

		status = close_connection(server);
    	if (status) {
      		fprintf(stderr, "tmsidata: err4\n");
     		goto cleanup;
		}
        cleanup:
		cleanup_event(reinterpret_cast<void**>(&event));
		cleanup_data(reinterpret_cast<void**>(&data));
		cleanup_header(reinterpret_cast<void**>(&header));
		cleanup_message(reinterpret_cast<void**>(&request));
		cleanup_message(reinterpret_cast<void**>(&response));

	//	pthread_exit(0);
//	return 0;
		if(Master->Stop())
		{	wprintf(L"\nDevice stop\n");
		}
		else
		{	wprintf(L"\nUnable to stop the device\n");

		}

	   //Unchain any slave devices
		Master->GetSlaveHandle();

	}
	else
	{	wprintf(L"\nUnable to start the Device\n");

	}

	//Clear memory before closing
//	for(Index=0;Index < MAX_DEVICE;Index++)
//	{	if( Device[Index] != NULL )
//		{	delete Device[Index];
//			Device[Index] = NULL;
//		}
//	}

	_getch();

    fclose (pOutputFile);
    pthread_exit(0);
	return 0;

}
