/*
 * Copyright (C) 2017, Robert Oostenveld
 * Donders Institute for Brain, Cognition and Behaviour; Radboud University; NL
 */

#ifdef __cplusplus
extern "C" {
#endif

#define DATATYPE_CHAR    (uint32_t)0
#define DATATYPE_UINT8   (uint32_t)1
#define DATATYPE_UINT16  (uint32_t)2
#define DATATYPE_UINT32  (uint32_t)3
#define DATATYPE_UINT64  (uint32_t)4
#define DATATYPE_INT8    (uint32_t)5
#define DATATYPE_INT16   (uint32_t)6
#define DATATYPE_INT32   (uint32_t)7
#define DATATYPE_INT64   (uint32_t)8
#define DATATYPE_FLOAT32 (uint32_t)9
#define DATATYPE_FLOAT64 (uint32_t)10

/* definition of simplified interface functions, see interface.c */
int start_server(int port);
int open_connection(const char *hostname, int port);
int close_connection(int s);
int read_header(int server, uint32_t *datatype, uint32_t *nchans, float *fsample, uint32_t *nsamples, uint32_t *nevents);
int read_data(int server, uint32_t begsample, uint32_t endsample, void *buffer);
int write_header(int server, uint32_t datatype, uint32_t nchans, float fsample);
int write_data(int server, uint32_t datatype, uint32_t nchans, uint32_t nsamples, void *buffer);
int wait_data(int server, uint32_t nsamples, uint32_t nevents, uint32_t milliseconds);

#ifdef __cplusplus
}
#endif
