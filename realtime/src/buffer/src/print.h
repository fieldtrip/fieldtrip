/*
 * Copyright (C) 2017, Robert Oostenveld
 * Donders Institute for Brain, Cognition and Behaviour; Radboud University; NL
 *
 */

#ifndef PRINT_H
#define PRINT_H

#include <stdio.h>
#include <stdlib.h>

#include "message.h"

#ifdef __cplusplus
extern "C" {
#endif

/* definition of functions for printing the content of various structures, see print.c */
void print_request(messagedef_t *);
void print_response(messagedef_t *);
void print_headerdef(headerdef_t *);
void print_datadef(datadef_t *);
void print_eventdef(eventdef_t *);
void print_datasel(datasel_t *);
void print_eventsel(eventsel_t *);
void print_buf(void *, int);

#ifdef __cplusplus
}
#endif

#endif
