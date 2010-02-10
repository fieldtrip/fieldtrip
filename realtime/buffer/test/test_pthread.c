#include <pthread.h>
#include <stdio.h>
#include <string.h>

#ifdef __BORLANDC__
#include <windows.h>
#include "gettimeofday.h"
#else
#include <unistd.h>
#include <strings.h>
#endif
#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#define bzero(b,len) (memset((b), '\0', (len)), (void) 0)
#define usleep(x)    (Sleep((x)/1000))
#endif
#if defined (__BORLANDC__) || defined (__WIN32__)
#define strcasecmp(a,b) (strcmpi(a,b))
#endif

pthread_mutex_t mutexthreadcount = PTHREAD_MUTEX_INITIALIZER;
int threadcount = 0;

void *hello_thread( void *arg ) {
		printf( "thread: threadcount = %d\n", threadcount );
		usleep(100000);
		pthread_mutex_lock(&mutexthreadcount);
		threadcount--;
		pthread_mutex_unlock(&mutexthreadcount);
		pthread_exit( 0 );
		return NULL;
}

int main( int argc, char *argv[] ) {
	int		n;
	pthread_t	tid;

	for (;;) {
		if ( n = pthread_create( &tid, NULL, hello_thread, NULL ) ) {
				fprintf( stderr, "pthread_create: %s\n", strerror( n ) );
				exit( 1 );
		}

		pthread_mutex_lock(&mutexthreadcount);
		threadcount++;
		pthread_mutex_unlock(&mutexthreadcount);
		printf("main: threadcount = %d\n", threadcount);

		if ( n = pthread_detach( tid ) ) {
			fprintf( stderr, "pthread_detach: %s\n", strerror( n ) );
			exit( 1 );
		}

		//if ( n = pthread_join( tid, NULL ) ) {
		//		fprintf( stderr, "pthread_join: %s\n", strerror( n ) );
		//		exit( 1 );
		//}

	}
	return( 0 );
}
