/*
 * Copyright (C) 2017, Robert Oostenveld
 * Donders Institute for Brain, COgnition and Behaviour; Radboud University; NL
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "buffer.h"

/*******************************************************************************
* CLOSE CONNECTION
*******************************************************************************/
int close_connection(int s) {
  int status = 0, verbose = 0;
  if (verbose>0)
  fprintf(stderr, "close_connection: socket = %d\n", s);
  if (s>0)
  status = closesocket(s);	/* it is a TCP connection */
  if (status!=0)
  perror("close_connection");
  return status;
}

/*******************************************************************************
* OPEN CONNECTION
*******************************************************************************/
int open_connection(const char *hostname, int port) {
  int verbose = 0;
  int s, retry;
  struct sockaddr_in sa;
  struct hostent *host;
  #ifdef PLATFORM_WINDOWS
  static WSADATA wsa = {0,0}; /* check version fields to only initialise once */
  #endif

  if (port==0) {
    if (verbose>0)
    fprintf(stderr, "open_connection: using direct memory copy\n");
    return 0;
  }
  else {
    if (verbose>0)
    fprintf(stderr, "open_connection: server = %s, port = %d\n", hostname, port);
  }

  #ifdef PLATFORM_WINDOWS
  if (wsa.wVersion == 0) {
    /* We only need to do this once ... and actually have a corresponding WSACleanup call somewhere */
    if(WSAStartup(MAKEWORD(1, 1), &wsa)) {
      fprintf(stderr, "open_connection: cannot start sockets\n");
      /* FIXME should this exception be handled more explicitely?  */
    }
  }
  #endif

  if ((host = gethostbyname(hostname)) == NULL) {
    fprintf(stderr, "open_connection: nslookup1 failed on '%s'\n", hostname);
    return -1;
  }

  if (host->h_length == 0) {
    fprintf(stderr, "open_connection: nslookup2 failed on '%s'\n", hostname);
    return -1;
  }

  bzero(&sa, sizeof sa);
  sa.sin_family = AF_INET;
  sa.sin_port = htons(port);

  memcpy(&(sa.sin_addr.s_addr), host->h_addr_list[0], sizeof(sa.sin_addr.s_addr));

  if ((s = socket(PF_INET, SOCK_STREAM, 0)) < 0) {
    if (verbose>0)
    fprintf(stderr, "open_connection: socket = %d\n", s);
    perror("open_connection");
    return -1;
  }

  retry = 10;
  while (retry>0) {
    if (connect(s, (struct sockaddr *)&sa, sizeof sa)<0) {
      /* wait 5 miliseconds and try again */
      usleep(5000);
      retry--;
    }
    else {
      /* this signals that the connection has been made */
      retry = -1;
    }
  }
  if (retry==0) {
    /* close the socket */
    closesocket(s);
    /* it failed on mutliple attempts, give up */
    return -2;
  }

  /*
  while (connect(s, (struct sockaddr *)&sa, sizeof sa) < 0) {
  perror("open_connection connect");
  usleep(1000000);
}
*/

if (verbose>0)
fprintf(stderr, "open_connection: connected to %s:%d on socket %d\n", hostname, port, s);

#ifdef DISABLE_NAGLE
{
  int optval = 1;
  setsockopt(s, IPPROTO_TCP, TCP_NODELAY, &optval, sizeof(optval));
}
#endif

return s;
}

/*******************************************************************************
* OPEN CONNECTION
*******************************************************************************/
int open_unix_connection(const char *name) {
  #ifdef PLATFORM_WINDOWS
  return -1;
}
#else
int verbose = 0;
int s, retry;
struct sockaddr_un sa;

bzero(&sa, sizeof(sa));
sa.sun_family = AF_UNIX;
strncpy(sa.sun_path, name, sizeof(sa.sun_path));

s = socket(AF_UNIX, SOCK_STREAM, 0);
if (s < 0) {
  perror("open_unix_connection, socket");
  return -1;
}

retry = 10;
while (retry>0) {
  if (connect(s, (struct sockaddr *)&sa, sizeof(sa))<0) {
    /* wait 5 miliseconds and try again */
    perror("open_connection");
    usleep(5000);
    retry--;
  } else {
    /* this signals that the connection has been made */
    retry = -1;
  }
}
if (retry==0) {
  /* it failed on mutliple attempts, give up */
  return -2;
}

if (verbose>0)
fprintf(stderr, "open_unix_connection: connected to %s on socket %d\n", name,  s);

return s;
}
#endif

/*******************************************************************************
 * WRITE HEADER
 *******************************************************************************/
int write_header(int server, UINT32_T datatype, int nchans, int fsample) {
    int verbose = 2, status;

		/* these are used in the communication and represent statefull information */
		message_t    *request  = NULL;
		message_t    *response = NULL;
		header_t     *header   = NULL;

		/* define the header packet */
    header      = (header_t *)malloc(sizeof(header_t));
		header->def = (headerdef_t *)malloc(sizeof(headerdef_t));
		header->buf = NULL;
		header->def->nchans    = nchans;
		header->def->nsamples  = 0;
		header->def->nevents   = 0;
		header->def->fsample   = fsample;
		header->def->data_type = datatype;
		header->def->bufsize   = 0;

    /* create the request */
    request      = (message_t *)malloc(sizeof(message_t));
		request->def = (messagedef_t *)malloc(sizeof(messagedef_t));
		request->buf = NULL;
		request->def->version = VERSION;
		request->def->bufsize = 0;
		request->def->command = PUT_HDR;
		request->def->bufsize = append(&request->buf, request->def->bufsize, header->def, sizeof(headerdef_t));
		request->def->bufsize = append(&request->buf, request->def->bufsize, header->buf, header->def->bufsize);

    /* send the request */
		status = clientrequest(server, request, &response);
    cleanup_message((void **)&request);

		if (verbose>0) fprintf(stderr, "DEBUG: clientrequest returned %d\n", status);
		if (status) {
				fprintf(stderr, "DEBUG: err1\n");
				exit(1);
		}

    if (response == NULL || response->def == NULL || response->def->command!=PUT_OK) {
      fprintf(stderr, "Error when writing samples.\n");
    }
    cleanup_message((void **)&response);

    return 0;
}

/*******************************************************************************
 * WRITE DATA
 *******************************************************************************/
int write_data(int server, UINT32_T datatype, int nchans, int nsamples, void *buffer) {
  int verbose = 2, status;

  /* these are used in the communication and represent statefull information */
  message_t    *request  = NULL;
  message_t    *response = NULL;
  data_t       *data     = NULL;

  /* define the data packet */
  data      = (data_t *)malloc(sizeof(data_t));
  data->def = (datadef_t *)malloc(sizeof(datadef_t));
  data->def->nchans    = nchans;
  data->def->nsamples  = nsamples;
  data->def->data_type = nsamples;
  data->def->bufsize   = wordsize_from_type(datatype)*nchans*nsamples;
  data->buf            = buffer;

  /* create the request */
  request      = (message_t *)malloc(sizeof(message_t));
  request->def = (messagedef_t *)malloc(sizeof(messagedef_t));
  request->buf = NULL;
  request->def->version = VERSION;
  request->def->bufsize = 0;
  request->def->command = PUT_DAT;
  request->def->bufsize = append(&request->buf, request->def->bufsize, data->def, sizeof(datadef_t));
  request->def->bufsize = append(&request->buf, request->def->bufsize, data->buf, data->def->bufsize);

  /* send the request */
  status = clientrequest(server, request, &response);
  cleanup_message((void **)&request);

  if (verbose>0) fprintf(stderr, "DEBUG: clientrequest returned %d\n", status);
  if (status) {
      fprintf(stderr, "DEBUG: err3\n");
      exit(1);
  }

  if (response == NULL || response->def == NULL || response->def->command!=PUT_OK) {
    fprintf(stderr, "Error when writing samples.\n");
  }
  cleanup_message((void **)&response);

  return 0;
};
