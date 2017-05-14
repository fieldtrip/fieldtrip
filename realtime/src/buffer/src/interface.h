/*
 * Copyright (C) 2017, Robert Oostenveld
 * Donders Institute for Brain, Cognition and Behaviour; Radboud University; NL
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef UINT32_T
	typedef uint32_t UINT32_T;
#endif

#define DATATYPE_CHAR    (UINT32_T)0
#define DATATYPE_UINT8   (UINT32_T)1
#define DATATYPE_UINT16  (UINT32_T)2
#define DATATYPE_UINT32  (UINT32_T)3
#define DATATYPE_UINT64  (UINT32_T)4
#define DATATYPE_INT8    (UINT32_T)5
#define DATATYPE_INT16   (UINT32_T)6
#define DATATYPE_INT32   (UINT32_T)7
#define DATATYPE_INT64   (UINT32_T)8
#define DATATYPE_FLOAT32 (UINT32_T)9
#define DATATYPE_FLOAT64 (UINT32_T)10

/* definition of simplified interface functions, see interface.c */
int start_server(int port);
int open_connection(const char *hostname, int port);
int close_connection(int s);
int read_header(int server, UINT32_T *datatype, unsigned int *nchans, float *fsample, unsigned int *nsamples, unsigned int *nevents);
int read_data(int server, unsigned int begsample, unsigned int endsample, void *buffer);
int write_header(int server, UINT32_T datatype, unsigned int nchans, float fsample);
int write_data(int server, UINT32_T datatype, unsigned int nchans, unsigned int nsamples, void *buffer);
int wait_data(int server, unsigned int nsamples, unsigned int nevents, unsigned int milliseconds);

#ifdef __cplusplus
}
#endif
