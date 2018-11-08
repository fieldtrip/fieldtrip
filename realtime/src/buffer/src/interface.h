/*
 * Copyright (C) 2017, Robert Oostenveld
 * Donders Institute for Brain, Cognition and Behaviour; Radboud University; NL
 */

#ifndef INTERFACE_H
#define INTERFACE_H

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#include "buffer.h"
#include "print.h"

#ifdef __cplusplus
extern "C" {
#endif

/* definition of simplified interface functions, see interface.c */
int start_server(int port);
int open_connection(const char *hostname, int port);
int close_connection(int s);
int read_header(int server, UINT32_T *datatype, UINT32_T *nchans, float *fsample, UINT32_T *nsamples, UINT32_T *nevents);
int read_data(int server, UINT32_T begsample, UINT32_T endsample, void *buffer);
int write_header(int server, UINT32_T datatype, UINT32_T nchans, float fsample);
int write_data(int server, UINT32_T datatype, UINT32_T nchans, UINT32_T nsamples, void *buffer);
int wait_data(int server, UINT32_T nsamples, UINT32_T nevents, UINT32_T milliseconds);

#ifdef __cplusplus
}
#endif

#endif
