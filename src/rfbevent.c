/*
 *
 * RFBEVENT sends a keyboard or mouse event to a VNC server
 *
 * Use as
 *   rfbevent(display, passwd, eventtype, eventvalue, ...)
 *
 * Some examples
 *  rfbevent('myserver:1', 'abc123', 'Text',   'xclock')          % type multiple characters
 *  rfbevent('myserver:1', 'abc123', 'Button', 'Return')          % single key event, press and release
 *  rfbevent('myserver:1', 'abc123', 'Button', 'Return',  0)      % single key event, press and release
 *  rfbevent('myserver:1', 'abc123', 'Button', 'Return',  1)      % single key event, press only
 *  rfbevent('myserver:1', 'abc123', 'Button', 'Return', -1)      % single key event, release only
 *  rfbevent('myserver:1', 'abc123', 'Pointer', [20 100])         % only mouse position
 *  rfbevent('myserver:1', 'abc123', 'Pointer', [20 100 1])       % mouse position and button 1, press and release
 *  rfbevent('myserver:1', 'abc123', 'Pointer', [20 100 1],  0)   % mouse position and button 1, press and release
 *  rfbevent('myserver:1', 'abc123', 'Pointer', [20 100 1],  1)   % mouse position and button 1, press only
 *  rfbevent('myserver:1', 'abc123', 'Pointer', [20 100 1], -1)   % mouse position and button 1, release only
 *
 * This code is based on rbfplaymacro version 0.2.2, see http://cyberelk.net/tim/software/rfbplaymacro
 * or http://people.redhat.com/twaugh/rfbplaymacro
 *
 * Copyright (C) 2000, 2001, 2002, 2003  Tim Waugh <twaugh@redhat.com>
 * Copyright (C) 2007 Robert Oostenveld <r.oostenveld@fcdonders.ru.nl>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * $Id$
 */

#include "platform.h"
#include "compiler.h"

#if defined(PLATFORM_WIN32) || defined(PLATFORM_WIN64)
#error The source code for this mex file has not yet been made compatible with Windows
/* the networking include files differ on windows and some of the functions are different */
/* a possibility to fix this would be to look into the fieldtrip buffer c-code, which is both unix and windows compatible */
#endif

#include <arpa/inet.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>
#include <getopt.h>
#include <netinet/in.h>
#include <unistd.h>
#include <wordexp.h>
#include <sys/time.h>
#include <sys/socket.h>
#include <sys/wait.h>
#include <netinet/in.h>
#include <netdb.h>

#include "d3des.h"
#include "mex.h"

#define VNC_BASE 5900
#define CHALLENGESIZE 16
#define CTS_BUF_SIZE 20    /* size of 'buf' used in connect_to_server, please make sure it's bigger than CHALLENGESIZE */
#define STRLEN 1024

static unsigned long int parse_keysym (const char *sym)
{
  size_t len = 0;
  char *end;
  unsigned long int key;
  
  /* Is it a single character?  That's easy. */
  if (*sym && !sym[1])
    return (unsigned long int)sym[0];
  
  /* Is it a number key code? */
  key = strtoul (sym, &end, 0);
  if (sym != end)
    return key;

       if (strcmp(sym, "Alt (left)"       )==0) { key = 0xffe9; }
  else if (strcmp(sym, "Alt (right)"      )==0) { key = 0xffea; }
  else if (strcmp(sym, "BackSpace"        )==0) { key = 0xff08; }
  else if (strcmp(sym, "Clear"            )==0) { key = 0xFF0B; }
  else if (strcmp(sym, "Control (left)"   )==0) { key = 0xffe3; }
  else if (strcmp(sym, "Control (right)"  )==0) { key = 0xffe4; }
  else if (strcmp(sym, "Delete"           )==0) { key = 0xFFFF; }
  else if (strcmp(sym, "Down"             )==0) { key = 0xff54; }
  else if (strcmp(sym, "End"              )==0) { key = 0xff57; }
  else if (strcmp(sym, "Enter"            )==0) { key = 0xff0d; }
  else if (strcmp(sym, "Esc"              )==0) { key = 0xFF1B; }
  else if (strcmp(sym, "Escape"           )==0) { key = 0xFF1B; }
  else if (strcmp(sym, "F1"               )==0) { key = 0xffbe; }
  else if (strcmp(sym, "F10"              )==0) { key = 0xffc7; }
  else if (strcmp(sym, "F11"              )==0) { key = 0xffc8; }
  else if (strcmp(sym, "F12"              )==0) { key = 0xffc9; }
  else if (strcmp(sym, "F2"               )==0) { key = 0xffbf; }
  else if (strcmp(sym, "F3"               )==0) { key = 0xffc0; }
  else if (strcmp(sym, "F4"               )==0) { key = 0xffc1; }
  else if (strcmp(sym, "F5"               )==0) { key = 0xffc2; }
  else if (strcmp(sym, "F6"               )==0) { key = 0xffc3; }
  else if (strcmp(sym, "F7"               )==0) { key = 0xffc4; }
  else if (strcmp(sym, "F8"               )==0) { key = 0xffc5; }
  else if (strcmp(sym, "F9"               )==0) { key = 0xffc6; }
  else if (strcmp(sym, "Home"             )==0) { key = 0xff50; }
  else if (strcmp(sym, "Insert"           )==0) { key = 0xff63; }
  else if (strcmp(sym, "Left"             )==0) { key = 0xff51; }
  else if (strcmp(sym, "Linefeed"         )==0) { key = 0xFF0A; }
  else if (strcmp(sym, "Meta (left)"      )==0) { key = 0xffe7; }
  else if (strcmp(sym, "Meta (right)"     )==0) { key = 0xffe8; }
  else if (strcmp(sym, "Page Down"        )==0) { key = 0xff56; }
  else if (strcmp(sym, "Page Up"          )==0) { key = 0xff55; }
  else if (strcmp(sym, "Pause"            )==0) { key = 0xFF13; }
  else if (strcmp(sym, "Return"           )==0) { key = 0xFF0D; }
  else if (strcmp(sym, "Right"            )==0) { key = 0xff53; }
  else if (strcmp(sym, "Scroll Lock"      )==0) { key = 0xFF14; }
  else if (strcmp(sym, "Shift (left)"     )==0) { key = 0xffe1; }
  else if (strcmp(sym, "Shift (right)"    )==0) { key = 0xffe2; }
  else if (strcmp(sym, "Sys Req"          )==0) { key = 0xFF15; }
  else if (strcmp(sym, "Tab"              )==0) { key = 0xff09; }
  else if (strcmp(sym, "Up"               )==0) { key = 0xff52; }
  else
  {
    char msg[STRLEN];
    snprintf(msg, (STRLEN-1), "Undefined keysym \"%s\"", sym);
    mexWarnMsgTxt(msg);
    return 0;
  }
  return key;
}

static ssize_t read_exact (int fd, void *buf, size_t len)
{
  ssize_t got;
  ssize_t need = len;
  do {
    got = read (fd, buf, need);
    if (got > 0) {
      buf += got;
      need -= got;
    }
    else return got;
  } while (need > 0);
  return len;
}

static ssize_t write_exact (int fd, const void *buf, size_t len)
{
  ssize_t wrote;
  ssize_t need = len;
  do {
    wrote = write (fd, buf, need);
    if (wrote > 0) {
      buf += wrote;
      need -= wrote;
    }
    else return wrote;
  } while (need > 0);
  return len;
}

/* This function ripped from vnc source as is (vncauth.c) */
void
vncEncryptBytes(unsigned char *bytes, char *passwd)
{
  unsigned char key[8];
  int i;
  
    /* key is simply password padded with nulls */
  
  for (i = 0; i < 8; i++) {
    if (i < strlen(passwd)) {
      key[i] = passwd[i];
    } else {
      key[i] = 0;
    }
  }
  
  deskey(key, EN0);
  
  for (i = 0; i < CHALLENGESIZE; i += 8) {
    des(bytes+i, bytes+i);
  }
}

static int connect_to_server (const char *display, int shared, const char *password)
{
  unsigned char buf[CTS_BUF_SIZE];
  uint32_t bit32;
  struct sockaddr_in sin;
  int s;
  char *server, *p, *end;
  unsigned long port;

  s = socket (PF_INET, SOCK_STREAM, IPPROTO_IP);
  if (s < 0)
    return s;

  /* Parse display as 'server:display' */
  server = strdup (display);
  p = strchr (server, ':');
  if (!p) {
    char msg[STRLEN];
    snprintf(msg, (STRLEN-1), "Invalid display: \"%s\"", display);
    mexWarnMsgTxt(msg);
    return -1;
  }
  *p++ = '\0';
  port = VNC_BASE + strtoul (p, &end, 10);
  if (p == end) {
    char msg[STRLEN];
    snprintf(msg, (STRLEN-1), "Invalid display number: %s", p);
    mexWarnMsgTxt(msg);
    return -1;
  }

  sin.sin_family = AF_INET;
  if (!server[0])
    sin.sin_addr.s_addr = htonl (INADDR_ANY);
  else {
    if ((sin.sin_addr.s_addr = inet_addr (server)) == -1) {
      struct hostent *hp;
      hp = gethostbyname (server);
      if (!hp) {
        char msg[STRLEN];
        snprintf(msg, (STRLEN-1), "Unknown host: %s", server);
        mexWarnMsgTxt(msg);
        exit (1);
      }
      memcpy (&sin.sin_addr.s_addr, hp->h_addr_list[0],
      hp->h_length);
    }
  }
  sin.sin_port = htons (port);

  if (connect (s, (struct sockaddr *) &sin, sizeof (sin))) {
    perror("connect");
    mexWarnMsgTxt ("Problem connecting to specified port");
    return -1;
  }

  /* ProtocolVersion from server */
  if (read_exact (s, buf, 12) < 12) {
    mexWarnMsgTxt("Couldn't read ProtocolVersion");
    close (s);
    return -1;
  }

  /* ProtocolVersion from client */
  write_exact (s, "RFB 003.003", 12);
  
  if (read_exact (s, &bit32, 4) < 4) {
    mexWarnMsgTxt("Couldn't read Authentication");
    close (s);
    return -1;
  }

  bit32 = ntohl (bit32);
  if (!bit32) {
    mexWarnMsgTxt("Connection failed");
    close (s);
    return -1;
  }

  if (bit32 == 2) {
    if (!password) {
      mexWarnMsgTxt("Need a password!");
      close (s);
      return -1;
    }

    /* Let's try to authenticate */
    if (read_exact (s, buf, CHALLENGESIZE) < CHALLENGESIZE) {
      mexWarnMsgTxt("Couldn't read DES Challenge");
      close(s);
      return -1;
    }

    /* buf now contains the server's 16 byte challenge */
    vncEncryptBytes(buf, (char *)password);

    write_exact(s, buf, CHALLENGESIZE);

    if (read_exact(s, &bit32, 4) < 4) {
      mexWarnMsgTxt( "DES Encrypted password sent, unable to read server response");
      close(s);
      return -1;
    }

    if (bit32 != 0) {
      mexWarnMsgTxt( "Authentication Failed");
      close(s);
      return -1;
    }

    /* Authentication successful!  wow!  let's set bit32 back to 1 so the next
       if statement doesn't bomb us out */
    bit32=1;
  }

  if (bit32 != 1) {
    char msg[STRLEN];
    snprintf(msg, (STRLEN-1), "Unknown authentication scheme: %d", bit32);
    mexWarnMsgTxt(msg);
    close (s);
    return -1;
  }

  /* ClientInitialisation from client */
  buf[0] = shared ? 1 : 0;
  write_exact (s, buf, 1);

  /* ServerInitialisation from server */
  if (read_exact (s, buf, 20) < 20) {
    mexWarnMsgTxt("Couldn't read ServerInitialisation");
    close (s);
    return -1;
  }

  if (read_exact (s, &bit32, 4) < 4) {
    mexWarnMsgTxt("Couldn't read ServerInitialisation");
    close (s);
    return -1;
  }

  bit32 = ntohl (bit32);
  while (bit32) {
    int size = 20;
    if (bit32 < size)
      size = bit32;
    read_exact (s, buf, size);
    bit32 -= size;
  }

  return s;
}

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  char display[STRLEN], password[STRLEN], type[STRLEN];
  int shared = 0;
  int socket, n, do_delay;
  struct timeval delay;

  delay.tv_sec  = 0;
  delay.tv_usec = 0;

  if (nrhs<3)
    mexErrMsgTxt ("Invalid number of input arguments");
  
  if (!mxIsChar(prhs[0]))
    mexErrMsgTxt ("Invalid input argument #1");
  if (!mxIsChar(prhs[1]))
    mexErrMsgTxt ("Invalid input argument #2");
  if (!mxIsChar(prhs[2]))
    mexErrMsgTxt ("Invalid input argument #3");

  n = mxGetNumberOfElements(prhs[0]); n = ( n>(STRLEN-1) ? (STRLEN-1) : n );
  mxGetString(prhs[0], display, n+1);

  n = mxGetNumberOfElements(prhs[1]); n = ( n>(STRLEN-1) ? (STRLEN-1) : n );
  mxGetString(prhs[1], password, n+1);

  n = mxGetNumberOfElements(prhs[2]); n = ( n>(STRLEN-1) ? (STRLEN-1) : n );
  mxGetString(prhs[2], type, n+1);

  /* 
   * mexPrintf("display  = %s\n", display);
   * mexPrintf("password = %s\n", password);
   * mexPrintf("type     = %s\n", type);
   */

  socket = connect_to_server (display, shared, password);
  if (socket < 0)
    mexErrMsgTxt ("Could not connect to remote display");

  if (strcmp(type, "Key")==0)
  {
    /*****************************************************************************
     *   rfbevent('myserver:1', 'abc123', 'Button', 'Return')          % single key event, press and release
     *   rfbevent('myserver:1', 'abc123', 'Button', 'Return',  0)      % single key event, press and release
     *   rfbevent('myserver:1', 'abc123', 'Button', 'Return',  1)      % single key event, press only
     *   rfbevent('myserver:1', 'abc123', 'Button', 'Return', -1)      % single key event, release only
     *****************************************************************************/
    int n, key, press = 0;
    char packet[20], str[STRLEN];
    uint32_t bit32;
    
    if (nrhs<4)
      mexErrMsgTxt ("Incorrect input arguments for Key event");
    if (!mxIsChar(prhs[3]))
      mexErrMsgTxt ("Incorrect input arguments for Key event");
    if (mxGetNumberOfElements(prhs[3])<1)
      mexErrMsgTxt ("Incorrect input arguments for Key event");
    
    n = mxGetNumberOfElements(prhs[3]); n = ( n>(STRLEN-1) ? (STRLEN-1) : n );
    mxGetString(prhs[3], str, n+1);
    key = parse_keysym (str);
    
    if (nrhs>4)
      press = mxGetScalar(prhs[4]);
    
    packet[0] = '\4';
    packet[2] = packet[3] = '\0';
    bit32 = htonl (key);
    memcpy (packet + 4, &bit32, 4);
    
    if (press>=0)
    {
      packet[1] = '\1'; /* key press */
      write_exact (socket, packet, 8);
    }
    if (press<=0)
    {
      packet[1] = '\0'; /* key release */
      write_exact (socket, packet, 8);
    }
  }

  else if (strcmp(type, "Pointer")==0)
  {
    /*****************************************************************************
     *   rfbevent('myserver:1', 'abc123', 'Pointer', [20 100])         % only mouse position
     *   rfbevent('myserver:1', 'abc123', 'Pointer', [20 100 1])       % mouse position and button 1, press and release
     *   rfbevent('myserver:1', 'abc123', 'Pointer', [20 100 1],  0)   % mouse position and button 1, press and release
     *   rfbevent('myserver:1', 'abc123', 'Pointer', [20 100 1],  1)   % mouse position and button 1, press only
     *   rfbevent('myserver:1', 'abc123', 'Pointer', [20 100 1], -1)   % mouse position and button 1, release only
     *****************************************************************************/
    char packet[20];
    double *p;
    int press = 0, number = 0, x, y;
    uint16_t xx, yy;
    uint8_t buttons = 0;

    if (nrhs<4)
      mexErrMsgTxt ("Incorrect input arguments for Pointer event");
    if (!mxIsDouble(prhs[3]))
      mexErrMsgTxt ("Incorrect input arguments for Pointer event");
    if (mxGetNumberOfElements(prhs[3])<2)
      mexErrMsgTxt ("Incorrect input arguments for Pointer event");

    p = mxGetPr(prhs[3]);
    x = (int)p[0];
    y = (int)p[1];
    if (x < 0) x = 0;
    if (y < 0) y = 0;

    /* get the optional button number */
    if (mxGetNumberOfElements(prhs[3])>2)
    {
      number = (int)p[2];
      if ((number<1) || (number>8))
        mexErrMsgTxt ("Incorrect button number for Pointer event (should be between 1 and 8)");
      buttons = (1<<(number-1));  /* set the appropriate bit for the button state */
    }

    /* get the press/release state of the button */
    if (nrhs>4)
      press = mxGetScalar(prhs[4]);

    /* construct the pointer event packet */
    packet[0] = '\5';
    xx = htons (x);
    yy = htons (y);
    memcpy (packet + 2, &xx, 2);
    memcpy (packet + 4, &yy, 2);

    if (press>=0)
    {
      packet[1] = buttons;      /* do the button press */
      write_exact (socket, packet, 6);
    }
    if (press<=0)
    {
      packet[1] = '\0';         /* do the button release, i.e. signal that no buttons are pressed */
      write_exact (socket, packet, 6);
    }
  }

  else if (strcmp(type, "Text")==0)
  {
    /*****************************************************************************
     * this is not a native RFB event type, but is included here for convenience
     * rfbevent('myserver:1', 'abc123', 'Text', 'xclock');
     *****************************************************************************/
    char packet[20], str[STRLEN];
    int i, n;   
    uint32_t bit32;

    if (nrhs<4)
      mexErrMsgTxt ("Incorrect input arguments for Text event");
    if (!mxIsChar(prhs[3]))
      mexErrMsgTxt ("Incorrect input arguments for Text event");
    if (mxGetNumberOfElements(prhs[3])<1)
      mexErrMsgTxt ("Incorrect input arguments for Text event");
    
    n = mxGetNumberOfElements(prhs[3]); n = ( n>(STRLEN-1) ? (STRLEN-1) : n );
    mxGetString(prhs[3], str, n+1);
    
    packet[0] = '\4';
    packet[2] = packet[3] = '\0';
    for (i=0; i<n; i++)
    {
      bit32 = htonl (str[i]);
      memcpy (packet + 4, &bit32, 4);
      packet[1] = '\1';  /* button press */
      write_exact (socket, packet, 8);
      packet[1] = '\0';  /* button release */
      write_exact (socket, packet, 8);
    }
  }
  
  if (do_delay)
    select (0, NULL, NULL, NULL, &delay);
  
  if (socket >= 0)
    close (socket);

  return;
}

