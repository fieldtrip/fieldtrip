/*
 * The low-level interface to the FieldTrip buffer is at the level of TCP sockets.
 * The medium-level interface to the FieldTrip buffer comprises request/response packets and the clientrequest function.
 * The high-level interface to the FieldTrip buffer is according to the following.
 *
 * #define DATATYPE_CHAR    
 * #define DATATYPE_UINT8   
 * #define DATATYPE_UINT16  
 * #define DATATYPE_UINT32  
 * #define DATATYPE_UINT64  
 * #define DATATYPE_INT8    
 * #define DATATYPE_INT16   
 * #define DATATYPE_INT32   
 * #define DATATYPE_INT64   
 * #define DATATYPE_FLOAT32 
 * #define DATATYPE_FLOAT64 
 * 
 * int start_server(int port);
 * int open_connection(const char *hostname, int port);
 * int close_connection(int s);
 * int read_header(int server, uint32_t *datatype, unsigned int *nchans, float *fsample, unsigned int *nsamples, unsigned int *nevents);
 * int read_data(int server, unsigned int begsample, unsigned int endsample, void *buffer);
 * int write_header(int server, uint32_t datatype, unsigned int nchans, float fsample);
 * int write_data(int server, uint32_t datatype, unsigned int nchans, unsigned int nsamples, void *buffer);
 * int wait_data(int server, unsigned int nsamples, unsigned int nevents, unsigned int milliseconds);
 */

#include <stdio.h>
#include <stdlib.h>

#include "platform_includes.h"
#include "interface.h"

#include "buffer.h"
#include "message.h"
#include "util.h"
#include "timestamp.h"

#define DIE(str) { fprintf(stderr, str " on line %d in file %s\n", __LINE__, __FILE__); exit(1); }

int main(int argc, char *argv[]) {
  int status, server = 0, nchans = 32, nsamples, datatype = DATATYPE_FLOAT32 ;
  float fsample = 250;
  void *buf;

  struct timespec before, after;

  if (argc<2) 
    DIE("incorrect number of arguments");

  if ((server = open_connection(argv[1], atoi(argv[2]))) < 1)
    DIE("cannot open connection");

  if ((status = write_header(server, datatype, nchans, fsample))!=0)
    DIE("cannot open connection");

  nsamples = fsample;
  if ((buf = malloc(nchans*nsamples*wordsize_from_type(datatype)))==0)
    DIE("cannot allocate memory");

  while (1) {
    get_monotonic_time(&before, TIMESTAMP_REF_BOOT);
    if ((status = write_data(server, datatype, nchans, nsamples, buf))!=0) {
      DIE("cannot write data");
    }
    else {
      fprintf(stderr, "wrote %d channels, %d samples, ", nchans, nsamples);
      get_monotonic_time(&after, TIMESTAMP_REF_BOOT);
      fprintf(stderr, "time = %f ms\n", 1000*get_elapsed_time(&before, &after));
    }
  }

  return 0;
}
