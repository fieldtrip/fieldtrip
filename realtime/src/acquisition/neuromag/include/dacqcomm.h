/* Automatically generated file, chages will not survive!*/
/* Copyright (c) 2010 Elekta Neuromag Oy, all rights reserved.*/
/* Generated  Fri Jan 15 10:42:20 EET 2010  by  mjk @ mjk.neuromag.fi  */

#ifndef dacqcomm_h_included
#define dacqcomm_h_included
/*
 *
 * Copyright 1994-2003
 *
 * Lauri Parkkonen
 * Neuromag, Ltd.
 * Helsinki, Finland
 *
 * No part of this program may be photocopied, reproduced,
 * or translated to another program language without the
 * prior written consent of the author.
 *
 * $Id$
 */

#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>

typedef struct sockaddr_in dacq_inet_address;

#define DACQ_REPLY_PACKET 1
#define DACQ_REPLY_RFC    2
#define DACQ_REPLY_BINARY 4
#define DACQ_REPLY_ASCII  8

#define DACQ_DRAIN_INPUT  0
#define DACQ_KEEP_INPUT   MSG_PEEK

#define DACQ_REPLY_GOOD   1
#define DACQ_REPLY_BAD    0
#define DACQ_REPLY_ERROR -1

#define DACQ_CMD_PASSWORD "pass"
#define DACQ_CMD_NAME     "name"
#define DACQ_CMD_ABOUT    "abou"
#define DACQ_CMD_MONITOR  "moni"
#define DACQ_CMD_HELP     "help"
#define DACQ_CMD_QUIT     "quit"

/* Automagic section follows */
#if defined(__cplusplus)
extern "C" {
#endif
/* client.c */
extern int dacq_server_set_write_timeout(int s);
extern int dacq_server_set_read_timeout(int s);
extern int dacq_server_get_write_timeout(void);
extern int dacq_server_get_read_timeout(void);
extern int dacq_server_send(int *server, char *buf, int nbytes, int flags);
extern int dacq_server_recv(int *server, char *buf, int bufsiz, int flags);
extern int dacq_server_check_reply(const char *buf);
extern int dacq_server_command(int *server, char *command, ...);
extern char *dacq_server_query(int *server, char *command, ...);
extern int dacq_server_close(int *server, char *quit_message);
extern int dacq_server_connect_by_addresses(dacq_inet_address *hostaddr);
extern int dacq_server_connect_by_address(dacq_inet_address *hostaddr, const char *servicename);
extern int dacq_server_connect_by_name(const char *hostname, const char *servicename);
extern int dacq_server_connect_by_name_quiet(const char *hostname, const char *servicename);
extern int dacq_server_login(int *server, const char *password, const char *myname);

#if defined(__cplusplus)
}
#endif
#endif /* of file */

