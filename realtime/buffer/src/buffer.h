/*
 * Copyright (C) 2008, Robert Oostenveld & Christian Hesse
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * $Log: buffer.h,v $
 * Revision 1.31  2009/01/23 08:26:44  roboos
 * fixed a serious bug that caused a lot of memory to leak (in fact all packets that were sent over the socket would eventually leack away), both on the client and server side
 *
 * Revision 1.30  2008/11/14 15:08:38  roboos
 * renamed open_server into open_connection (idem for close)
 *
 * Revision 1.29  2008/10/29 20:42:46  roboos
 * added #define for backward compatibility support for open_remotehost
 *
 * Revision 1.28  2008/10/29 20:12:27  roboos
 * renamed open_remotehost into open_connection and added function close_connection
 * added some fprintf statements, no functional changes
 *
 * Revision 1.27  2008/07/09 10:37:53  roboos
 * moved some printing functions to seperate file
 *
 * Revision 1.26  2008/07/01 17:08:24  thohar
 * added cleanup_buffer function
 *
 * Revision 1.25  2008/06/19 21:07:04  roboos
 * switched back to more reasonable numbers for the data and event buffers
 *
 * Revision 1.24  2008/06/19 20:32:08  roboos
 * fixed bug in WRAP due to missing () and only when pointers were involved
 *
 * Revision 1.23  2008/06/19 19:21:52  roboos
 * ensure that wrapping also works for N+x where x>1
 *
 * Revision 1.22  2008/05/29 07:54:59  roboos
 * moved checks for the packing of structs and wordsize to seperate function
 *
 * Revision 1.21  2008/05/22 09:50:54  roboos
 * moved win32 specific into seperate header file
 *
 * Revision 1.20  2008/03/26 14:34:43  thohar
 * defines closesocket now as close on non-win32 platforms
 *
 * Revision 1.19  2008/03/19 09:21:00  thohar
 * added extern "C" statement for c++ builds
 *
 * Revision 1.18  2008/03/17 13:43:12  roboos
 * added client helper functions set/get_property
 *
 * Revision 1.17  2008/03/13 13:37:26  roboos
 * added DEFAULT_HOSTNAME, renamed PORT to DEFAULT_PORT
 * removed declaration of functions that are currently not used (sorry Christian)
 *
 * Revision 1.16  2008/03/13 12:31:11  thohar
 * Added defines for implementations of poll.h that do not have e.g. POLLRDNORM
 *
 * Revision 1.15  2008/03/10 10:02:42  roboos
 * renamed host.host into host.name
 *
 * Revision 1.14  2008/03/10 09:40:18  roboos
 * added property details, added host structure (for name and port)
 *
 * Revision 1.13  2008/03/08 10:37:40  roboos
 * some changes to reflect the renaming of the low-level functions and associated files
 *
 * Revision 1.12  2008/03/07 14:48:43  roboos
 * added declaration for write_request
 *
 * Revision 1.11  2008/03/02 13:24:09  roboos
 * changed SO_RCVBUF_SIZE and SO_SNDBUF_SIZE to 16k
 * added declaration of handle_request function
 *
 * Revision 1.10  2008/02/27 10:13:27  roboos
 * added print_buf declaration
 *
 * Revision 1.9  2008/02/26 21:43:25  roboos
 * renamed packet_t structure definition into messagedef_t, added message_t structure (contains def+buf)
 *
 * Revision 1.8  2008/02/20 13:35:54  roboos
 * changed comments to ansi style, needed for matlab
 * modified declaration of append function
 *
 * Revision 1.7  2008/02/20 07:10:25  roboos
 * tried out some low-level tcp specific details
 *
 * Revision 1.6  2008/02/19 10:24:26  roboos
 * added copyright statement
 *
 * Revision 1.5  2008/02/19 08:29:13  roboos
 * added print_datasel and print_eventsel declarations
 * added bufread and bufwrite declarations
 *
 * Revision 1.4  2008/02/18 17:04:08  roboos
 * lots of small changes, debugging for brainamp
 *
 * Revision 1.3  2008/02/18 12:13:46  roboos
 * moved executable from buffer to demo
 * fixed bugs in sinewave and socket for events
 * stripped down the eventdef_t fields
 * many small changes
 *
 * Revision 1.2  2008/02/18 10:20:48  roboos
 * merged bufferlib.h into buffer.h
 *
 * Revision 1.1  2008/02/18 10:05:25  roboos
 * restructured the directory layout, copied all code into src, added directory for external code
 *
 * Revision 1.4  2008/02/13 13:59:03  roboos
 * made a start with implementing GET_EVT
 *
 * Revision 1.3  2008/02/13 13:07:56  roboos
 * fixed numerous bugs
 * implemented append function, which after all does not seem to be neccessary (since the problem is in the TCP buffer size)
 *
 * Revision 1.2  2008/02/10 10:47:52  roboos
 * redefined max number of data and sampels, added define for wrapping
 *
 */

/* prevent double include */
#ifndef BUFFER_H
#define BUFFER_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "platform_includes.h"
#include "message.h"

#ifndef POLLRDNORM
#define POLLRDNORM POLLIN
#endif

#ifndef POLLRDBAND
#define POLLRDBAND POLLPRI
#endif

#ifndef POLLWRNORM
#define POLLWRNORM POLLOUT
#endif

#ifndef POLLWRBAND
#define POLLWRBAND POLLOUT
#endif

#define BACKLOG           16
#define DEFAULT_HOSTNAME "localhost"
#define DEFAULT_PORT      1972

#define SO_RCVBUF_SIZE 16384
#define SO_SNDBUF_SIZE 16384

/* this is because the function  has been renamed, but is perhaps already in use in other software */
#define open_remotehost open_connection

/* FIXME these should be variable */
#define MAXNUMSAMPLE    600000
#define MAXNUMEVENT     100
#define MAXNUMPROPERTY  100

#define WRAP(x,y) ((x) - ((int)((float)(x)/(y)))*(y))
#define FREE(x) {if (x) {free(x); x=NULL;}}

#ifdef __cplusplus
extern "C" {
	
#endif

/* declaration of "public" buffer API functions */
int read_header( const char *hostname, int port,           void **ppw);
int read_data(   const char *hostname, int port, int *pnw, void **ppw);
int read_event(  const char *hostname, int port, int *pnw, void **ppw);
int write_header(const char *hostname, int port,           void **ppw);
int write_data(  const char *hostname, int port, int *pnw, void **ppw);
int write_event( const char *hostname, int port, int *pnw, void **ppw);
int flush_header(const char *hostname, int port);
int flush_data(  const char *hostname, int port);
int flush_event( const char *hostname, int port);
void cleanup_buffer();

/* definition of the functions that implement the network transparent server */
void *tcpserver(void *);
void *tcpsocket(void *);

/* definition of test functions that emulate an acquisition system */
void *sinewave_thread(void *);
void *event_thread(void *);

/* definition of the functions used in thread cancelation, see cleanup.c */
void cleanup_message(void **arg);
void cleanup_header(void **arg);
void cleanup_data(void **arg);
void cleanup_event(void **arg);
void cleanup_property(void **arg);
void cleanup_buf(void **arg);
void cleanup_socket(int *);

/* definition of helper functions for debugging and printing the content of various structures */
void print_request(messagedef_t *);
void print_response(messagedef_t *);
void print_headerdef(headerdef_t *);
void print_datadef(datadef_t *);
void print_eventdef(eventdef_t *);
void print_propertydef(propertydef_t *);
void print_datasel(datasel_t *);
void print_eventsel(eventsel_t *);
void print_propertydef(propertydef_t *);
void print_buf(void *, int);

/* definition of even more helper functions, see util.c */
int open_connection(const char*, int);
int close_connection(int);
int append(void **, int, void *, int);
int bufread(int, void *, int);
int bufwrite(int, void *, int);
int clientrequest(int, message_t *, message_t**);
int dmarequest(message_t *, message_t**);
int tcprequest(int, message_t *, message_t**);
int find_property(property_t *);
int get_property(int, const char *, INT32_T *);
int set_property(int, const char *, INT32_T *);
void check_datatypes(void);

typedef struct {
	char name[256];
	int  port;
} host_t;

#ifdef __cplusplus
}
#endif

#endif

