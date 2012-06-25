#ifndef EXTERN_H
#define EXTERN_H

#include "buffer.h"

extern pthread_mutex_t mutexstatus;
extern int tcpserverStatus;

extern pthread_mutex_t mutexthreadcount;
extern int threadcount;

extern pthread_mutex_t mutexsocketcount;
extern int socketcount;

extern pthread_mutex_t mutexthreadcount;
extern int threadcount;

extern pthread_mutex_t mutexappendcount;
extern int appendcount;

#endif
