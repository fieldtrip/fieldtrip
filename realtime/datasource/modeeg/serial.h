/* Simple serial port library for Windows and Linux, written in plain C
 * (C) 2008 Stefan Klanke
 */

#ifndef __SERIAL_H
#define __SERIAL_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifdef WIN32

#include <windows.h>

typedef struct {
   HANDLE comPort;
   DCB oldDCB;
   COMMTIMEOUTS oldTimeOuts;
} SerialPort;

#else

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <termios.h>
#include <unistd.h>

typedef struct {
   int comPort;
   struct termios oldTermios;
} SerialPort;

#endif

int serialOpenByNumber(SerialPort *SP, int port);
int serialOpenByName(SerialPort *SP, const char *name);
int serialClose(SerialPort *SP);
int serialWrite(SerialPort *SP, int size, void *buffer);
int serialRead(SerialPort *SP, int size, void *buffer);
int serialSetParameters(SerialPort *SP, int baudrate, int bits, int parity, int stops, int timeout);
void serialFlushInput(SerialPort *SP);
void serialFlushOutput(SerialPort *SP);
int serialInputPending(SerialPort *SP);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
