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

#define DIE(str) { fprintf(stderr, str " on line %d in file %s\n", __LINE__, __FILE__); exit(1); }

int main(int argc, char *argv[]) {
  int server;

  if (argc<2)
    DIE("incorrect number of arguments");

  /* the following starts the server in a pthread */
  if ((server = start_server(atoi(argv[1])))!=0) 
    DIE("cannot open server");
 
  while (1) {
    fprintf(stderr, "server running\n");
    usleep(1000000);
  }

  return 0;
}
