#include <stdio.h>
#include <pthread.h>
#include <sys/select.h>


pthread_mutex_t my_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t  my_cond  = PTHREAD_COND_INITIALIZER;

void *threadfunc(void *parm)
{
		struct timeval timeout;

		printf("Entered secondary thread\n");
		while (1) {

				printf("Secondary thread is looping\n");
				pthread_testcancel();

				pthread_mutex_lock(&my_mutex);
				pthread_cond_signal(&my_cond);
				pthread_mutex_unlock(&my_mutex);

				timeout.tv_sec = 1;
				timeout.tv_usec = 0;
				select( 0, NULL, NULL, NULL, & timeout );
		}
		return NULL;
}

int main(int argc, char **argv)
{
		pthread_t             thread;
		int                   rc=0;
		int i = 0;

		while (1) {
				rc = pthread_create(&thread, NULL, threadfunc, NULL);

				pthread_mutex_lock(&my_mutex);
				pthread_cond_wait(&my_cond, &my_mutex);
				pthread_mutex_unlock(&my_mutex);

				rc = pthread_cancel(thread);
				printf("%d\n", i++);
		}

		printf("Main completed\n");
		return 0;
}
