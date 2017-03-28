/*
 * Copyright (C) 2017, Robert Oostenveld
 * Donders Institute for Brain, Cognition and Behaviour; Radboud University; NL
 */

#ifdef __cplusplus
extern "C" {
#endif

/* definition of simplified interface functions, see interface.c */
int close_connection(int s);
int open_connection(const char *hostname, int port);
int open_unix_connection(const char *name);
int write_header(int server, UINT32_T datatype, int nchans, int fsample);
int write_data(int server, UINT32_T datatype, int nchans, int nsamples, void *buffer);

#ifdef __cplusplus
}
#endif

