/* Simple serial port library for Windows and Linux, written in plain C
 * (C) 2008 Stefan Klanke
 */

#include <stdio.h>
#include <string.h>

#include "serial.h"
#include "platform.h"
#include "compiler.h"
#include "platform_includes.h"

#define PRINTLINE { fprintf(stderr, "%s in %s at line %d\n", __FILE__, __func__, __LINE__); }

static const char serialErrOpen[]="Could not open the serial port.\n";
static const char serialErrGetP[]="Could not read serial port parameters.\n";
static const char serialErrGetT[]="Could not read serial port timeouts.\n";
static const char serialErrSetP[]="Could not set serial port parameters.\n";
static const char serialErrSetT[]="Could not set serial port timeout.\n";

#if defined(PLATFORM_WINDOWS)

#include <windows.h>

int serialOpenByNumber(SerialPort *SP, int port) {
  char device[32];

  sprintf(device,"\\\\.\\COM%d",port);

  return serialOpenByName(SP, device);
}

int serialOpenByName(SerialPort *SP, const char *device) {
  SP->comPort = CreateFile(device, GENERIC_READ|GENERIC_WRITE, 0, NULL, OPEN_EXISTING, 0, NULL);
  if (SP->comPort == INVALID_HANDLE_VALUE) {
    fputs(serialErrOpen,stderr);
    SP->comPort = NULL;
    return 0;
  }

  if (!GetCommState(SP->comPort, &(SP->oldDCB))) {
    fputs(serialErrGetP,stderr);
    CloseHandle(SP->comPort);
    SP->comPort=NULL;
    return 0;
  }
  if (!GetCommTimeouts(SP->comPort, &(SP->oldTimeOuts))) {
    fputs(serialErrGetP,stderr);
    CloseHandle(SP->comPort);
    SP->comPort=NULL;
    return 0;
  }
  PurgeComm(SP->comPort, PURGE_RXCLEAR | PURGE_TXCLEAR);
  return 1;
}


int serialSetParameters(SerialPort *SP, int baudrate, int bits, int parity, int stops, int timeout) {
  DCB newDCB;
  COMMTIMEOUTS newTO;

  memcpy(&newDCB, &(SP->oldDCB), sizeof(DCB));

  newDCB.BaudRate = baudrate;
  newDCB.ByteSize = bits;
  if (parity==1) {
    newDCB.Parity = EVENPARITY;
    newDCB.fParity = TRUE;
  } else {
    newDCB.Parity = NOPARITY;
    newDCB.fParity = FALSE;
  }
  /*
     newDCB.fOutX = FALSE;
     newDCB.fInX = FALSE;
     newDCB.fTXContinueOnXoff = TRUE;
     newDCB.fNull = FALSE;
   */
  switch(stops) {
    case 1:
      newDCB.StopBits = ONESTOPBIT;
      break;
    case 2:
      newDCB.StopBits = TWOSTOPBITS;
      break;
    default:
      newDCB.StopBits = ONESTOPBIT;
  }

  if (!SetCommState(SP->comPort, &newDCB)) {
    fprintf(stderr,"Couldn't change serial port settings, error = %li\n",GetLastError());
    return 0;
  }

  newTO.ReadIntervalTimeout = timeout*100; /* ms */
  newTO.ReadTotalTimeoutMultiplier = 0;
  newTO.ReadTotalTimeoutConstant = timeout*100; /* ms */
  newTO.WriteTotalTimeoutMultiplier = 0;
  newTO.WriteTotalTimeoutConstant = timeout*100; /* ms */

  if (!SetCommTimeouts(SP->comPort, &newTO)) {
    fputs("Couldn't set serial port timeouts\n",stderr);
    return 0;
  }
  return 1;
}

int serialClose(SerialPort *SP) {
  if (SP->comPort == NULL) return 0;

  PurgeComm(SP->comPort, PURGE_RXCLEAR | PURGE_TXCLEAR);

  SetCommState(SP->comPort, &(SP->oldDCB));
  SetCommTimeouts(SP->comPort, &(SP->oldTimeOuts));
  CloseHandle(SP->comPort);
  SP->comPort = NULL;
  /* TODO: think about error values */
  return 1;
}

int serialWrite(SerialPort *SP, int size, void *buffer) {
  DWORD numWritten;

  WriteFile(SP->comPort, buffer, size, &numWritten, NULL);
  return numWritten;
}

void serialFlushInput(SerialPort *SP) {
  PurgeComm(SP->comPort, PURGE_RXCLEAR);
}

void serialFlushOutput(SerialPort *SP) {
  PurgeComm(SP->comPort, PURGE_TXCLEAR);
}

int serialRead(SerialPort *SP, int size, void *buffer) {
  DWORD numRead;

  ReadFile(SP->comPort, buffer, size, &numRead, NULL);
  return numRead;
}

int serialInputPending(SerialPort *SP) {
  DWORD numErrors;
  COMSTAT stat;

  if (!ClearCommError(SP->comPort, &numErrors, &stat)) return -1;
  if (numErrors > 0) return -1;
  return stat.cbInQue;
}

#else /* Linux, POSIX */

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <termios.h>
#include <time.h>
#include <unistd.h>
#include <sys/ioctl.h>

int serialOpenByNumber(SerialPort *SP, int port) {
  char device[16];

  snprintf(device,16,"/dev/ttyS%d",port);
  return serialOpenByName(SP, device);
}


int serialOpenByName(SerialPort *SP, const char *device) {

  SP->comPort = open(device, O_RDWR | O_NOCTTY | O_NDELAY);

  /*
   * This used to be without the O_NDELAY for the modeeg implementation,
   * however that fails to work for the openbci implementation. It
   * is unclear whether the O_NDELAY is a problem for modeeg, but for
   * now lets assume that it can be added for both.
   *
   * SP->comPort = open(device, O_RDWR | O_NOCTTY);
   */

  if (!SP->comPort) {
    fputs(serialErrOpen, stderr);
    return 0;
  }

  if (tcgetattr(SP->comPort, &(SP->oldTermios))) {
    perror("tcgetattr");
    fputs(serialErrGetP, stderr);
    close(SP->comPort);
    SP->comPort=0;
    return 0;
  }
  tcflush(SP->comPort, TCIOFLUSH);
  return 1;
}


int serialSetParameters(SerialPort *SP, int baudrate, int bits, int parity, int stops, int timeout) {
  struct termios newtio;

  memset(&newtio, 0, sizeof(newtio)); /* clear struct for new port settings */

  newtio.c_cflag = CLOCAL | CREAD;

  if (stops>1) {
    newtio.c_cflag |= CSTOPB;
  }

  switch(baudrate) {
    case    300:
      cfsetispeed(&newtio, B300);
      cfsetospeed(&newtio, B300);
      break;
    case    600:
      cfsetispeed(&newtio, B600);
      cfsetospeed(&newtio, B600);
      break;
    case   1200:
      cfsetispeed(&newtio, B1200);
      cfsetospeed(&newtio, B1200);
      break;
    case   2400:
      cfsetispeed(&newtio, B2400);
      cfsetospeed(&newtio, B2400);
      break;
    case   4800:
      cfsetispeed(&newtio, B4800);
      cfsetospeed(&newtio, B4800);
      break;
    case   9600:
      cfsetispeed(&newtio, B9600);
      cfsetospeed(&newtio, B9600);
      break;
    case  19200:
      cfsetispeed(&newtio, B19200);
      cfsetospeed(&newtio, B19200);
      break;
    case  38400:
      cfsetispeed(&newtio, B38400);
      cfsetospeed(&newtio, B38400);
      break;
    case  57600:
      cfsetispeed(&newtio, B57600);
      cfsetospeed(&newtio, B57600);
      break;
    case 115200:
      cfsetispeed(&newtio, B115200);
      cfsetospeed(&newtio, B115200);
      break;
    case 230400:
      cfsetispeed(&newtio, B230400);
      cfsetospeed(&newtio, B230400);
      break;
    default: fputs("Unrecognized baudrate\n.",stderr); return 0;
  }


  switch(bits) {
    case 5: newtio.c_cflag |= CS5; break;
    case 6: newtio.c_cflag |= CS6; break;
    case 7: newtio.c_cflag |= CS7; break;
    case 8: newtio.c_cflag |= CS8; break;
    default: fputs("Unrecognized number of bits\n.",stderr); return 0;
  }

  if (parity) {
    newtio.c_cflag |= PARENB;
  }

  newtio.c_iflag = IGNBRK | IXOFF;
  newtio.c_oflag = 0;
  newtio.c_lflag = 0;

  newtio.c_cc[VTIME]    = timeout;     /* deciseconds */
  newtio.c_cc[VMIN ]    = 0;

  if (tcsetattr(SP->comPort, TCSANOW, &newtio)) {
    perror("tcsetattr");
    fputs("Couldn't change serial port settings\n",stderr);
    return 0;
  }
  return 1;
}

int serialClose(SerialPort *SP) {

  tcflush(SP->comPort, TCIOFLUSH);
  tcsetattr(SP->comPort,TCSANOW,&(SP->oldTermios));
  close(SP->comPort);
  return 1; /* TODO: think about error values */
}

int serialWrite(SerialPort *SP, int size, void *buffer) {
  int numWrite = write(SP->comPort, buffer, size);;

  return numWrite;
}

int serialRead(SerialPort *SP, int size, void *buffer) {
  return read(SP->comPort, buffer, size);
}

void serialFlushInput(SerialPort *SP) {
  tcflush(SP->comPort, TCIFLUSH);
}

void serialFlushOutput(SerialPort *SP) {
  tcflush(SP->comPort, TCOFLUSH);
}

int serialInputPending(SerialPort *SP) {
  int bytesWaiting, res;
  res = ioctl(SP->comPort, FIONREAD, &bytesWaiting);
  if (res < 0) return res;
  return bytesWaiting;
}

#endif
