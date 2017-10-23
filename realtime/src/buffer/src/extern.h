/* 
 * Copyright (C) 2008, Robert Oostenveld & Christian Hesse
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 */

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
