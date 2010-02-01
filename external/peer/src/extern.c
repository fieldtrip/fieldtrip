#include "peer.h"

pthread_mutex_t mutexstatus = PTHREAD_MUTEX_INITIALIZER;
int tcpserverStatus = 0;
int announceStatus  = 0;
int discoverStatus  = 0;
int expireStatus    = 0;

pthread_mutex_t mutexappendcount = PTHREAD_MUTEX_INITIALIZER;
int appendcount = 0;

pthread_mutex_t mutexsocketcount = PTHREAD_MUTEX_INITIALIZER;
int socketcount = 0;

pthread_mutex_t mutexthreadcount = PTHREAD_MUTEX_INITIALIZER;
int threadcount = 0; 

pthread_mutex_t mutexhost = PTHREAD_MUTEX_INITIALIZER;
hostdef_t *host = NULL;

pthread_mutex_t mutexpeerlist = PTHREAD_MUTEX_INITIALIZER;
peerlist_t *peerlist = NULL;

pthread_mutex_t mutexjoblist = PTHREAD_MUTEX_INITIALIZER;
joblist_t *joblist = NULL;

pthread_mutex_t mutexuserlist = PTHREAD_MUTEX_INITIALIZER;
userlist_t *userlist = NULL;

pthread_mutex_t mutexgrouplist = PTHREAD_MUTEX_INITIALIZER;
grouplist_t *grouplist = NULL;

pthread_mutex_t mutexhostlist = PTHREAD_MUTEX_INITIALIZER;
hostlist_t *hostlist = NULL;

pthread_mutex_t mutexfairshare = PTHREAD_MUTEX_INITIALIZER;
struct {
  float a;
  float b;
  float c;
  float d;
  int   n;
  int prevhostcount;
  int prevhostid;
  time_t t0;
} param;

