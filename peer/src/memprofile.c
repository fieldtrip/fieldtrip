/*
 * MEMPROFILE is a MATLAB mex file which can be used to sample the
 * memory useage while MATLAB is running arbitrary commands. A thread
 * is started in the background, which takes a sample of the memory
 * use every second.
 *
 * Copyright (C) 2010, Robert Oostenveld
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "mex.h"
#include "matrix.h"
#include "platform.h"
#include "platform_includes.h"

#define FREE(x) {if (x) {free(x); x=NULL;}}

/* how often to record a sample, start with 10x per second */
#define SAMPLINGDELAY 100000

/* the list will be thinned out with a factor 2x if it exceeds this length and the sampling delay will be doubled */
#define MAXLISTLENGTH 3600  

typedef struct memlist_s {
		uint64_t rss;  /* size of resident memory */
		uint64_t vs;   /* size of virtual memory */
		time_t   time; /* time at which the measurement was taken */
		struct memlist_s *next;
} memlist_t;

time_t reftime = 0;
memlist_t *memlist = NULL;
pthread_mutex_t mutexmemlist = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutexmemprofile  = PTHREAD_MUTEX_INITIALIZER;
pthread_t memprofileThread;
int memprofileStatus = 0;

void memprofile_cleanup(void *arg) {
		pthread_mutex_lock(&mutexmemprofile);	
		memprofileStatus = 0;
		pthread_mutex_unlock(&mutexmemprofile);	
}

/* this function takes a sample of the memory use and adds it to the list */
void memprofile_sample(void) {
		uint64_t rss, vs;
		memlist_t *memitem;

		if (getmem(&rss, &vs)!=0) {
				rss = 0;
				vs  = 0;
		}

		memitem = (memlist_t *)malloc(sizeof(memlist_t));
		memitem->rss  = rss;
		memitem->vs   = vs;
		memitem->time = time(NULL) - reftime;

		pthread_mutex_lock(&mutexmemlist);
		memitem->next = memlist;
		memlist       = memitem;
		pthread_mutex_unlock(&mutexmemlist);
}

/* this is the thread that records the memory usage in a linked list */
void *memprofile(void *arg) {
		int count;
		int pause = SAMPLINGDELAY;
		memlist_t *memitem, *next;

		pthread_mutex_lock(&mutexmemprofile);	
		if (memprofileStatus) {
				/* only a single instance should be running */
				pthread_mutex_unlock(&mutexmemprofile);	
				return NULL;
		}
		else {
				memprofileStatus = 1;
				pthread_mutex_unlock(&mutexmemprofile);	
		}

		pthread_cleanup_push(memprofile_cleanup, NULL);

		while (1) {
				/* count the number of items on the list */
				count = 0;
				pthread_mutex_lock(&mutexmemlist);
				memitem = memlist;
				while (memitem) {
						count++;
						memitem = memitem->next ;
				}
				pthread_mutex_unlock(&mutexmemlist);

				if (count>MAXLISTLENGTH) {
						/* remove every 2nd item from the list */
						pthread_mutex_lock(&mutexmemlist);
						memitem = memlist;
						next    = memitem->next;
						while (memitem && next) {
								memitem->next = next->next;
								FREE(next);
								memitem = memitem->next;
								next    = memitem->next;
						}
						pthread_mutex_unlock(&mutexmemlist);
						/* wait an extra delay before taking the next sample */
						usleep(pause);
						/* increment the subsequent sampling delay with a factor 2x */
						pause = pause*2;
				}

				memprofile_sample();
				pthread_testcancel();
				usleep(pause);
		}

		pthread_cleanup_pop(1);
		return NULL;
}

/* this function will be called upon unloading of the mex file */
void exitFun(void) {
		pthread_mutex_lock(&mutexmemprofile);	
		if (memprofileStatus) {
				memprofileStatus = 0;
				pthread_mutex_unlock(&mutexmemprofile);	
				pthread_cancel(memprofileThread);
				pthread_join(memprofileThread, NULL);
		}
		else {
				pthread_mutex_unlock(&mutexmemprofile);	
		}
}

void mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]) {
		char *command = NULL;
		int rc, i;
		memlist_t *memitem = NULL;
		const char *fieldnames[3] = {"time", "mem"};
		int numfields = 2;

		/* this function will be called upon unloading of the mex file */
		mexAtExit(exitFun);

		/* the first argument is the command string */
		if (nrhs<1)
				mexErrMsgTxt ("invalid number of input arguments");
		if (!mxIsChar(prhs[0]))
				mexErrMsgTxt ("invalid input argument #1");
		command = mxArrayToString(prhs[0]);

		/****************************************************************************/
		if      (strcasecmp(command, "on")==0) {

				/* clear the existing list */
				pthread_mutex_lock(&mutexmemlist);
				memitem = memlist;
				while (memitem) {
						memlist = memitem->next;
						FREE(memitem);
						memitem = memlist;
				}
				pthread_mutex_unlock(&mutexmemlist);

				/* set the reference time */
				reftime = time(NULL);

				/* include the current memory use in the list */
				memprofile_sample();

				pthread_mutex_lock(&mutexmemprofile);
				if (memprofileStatus) {
						pthread_mutex_unlock(&mutexmemprofile);
						return;
				}
				else {
						pthread_mutex_unlock(&mutexmemprofile);
				}

				/* start the thread */
				rc = pthread_create(&memprofileThread, NULL, memprofile, (void *)NULL);
				if (rc)
						mexErrMsgTxt("problem with return code from pthread_create()");
		}

		/****************************************************************************/
		else if (strcasecmp(command, "off")==0) {

				pthread_mutex_lock(&mutexmemprofile);
				if (!memprofileStatus) {
						pthread_mutex_unlock(&mutexmemprofile);
						return;
				}
				else {
						pthread_mutex_unlock(&mutexmemprofile);
				}

				/* stop the thread */
				rc = pthread_cancel(memprofileThread);
				if (rc)
						mexErrMsgTxt("problem with return code from pthread_cancel()");
				else
						memprofileStatus = 0;
		}

		/****************************************************************************/
		else if (strcasecmp(command, "resume")==0) {

				pthread_mutex_lock(&mutexmemprofile);
				if (memprofileStatus) {
						pthread_mutex_unlock(&mutexmemprofile);
						return;
				}
				else {
						pthread_mutex_unlock(&mutexmemprofile);
				}

				/* start the thread */
				rc = pthread_create(&memprofileThread, NULL, memprofile, (void *)NULL);
				if (rc)
						mexErrMsgTxt("problem with return code from pthread_create()");
		}

		/****************************************************************************/
		else if (strcasecmp(command, "clear")==0) {

				/* clear the existing list */
				pthread_mutex_lock(&mutexmemlist);
				memitem = memlist;
				while (memitem) {
						memlist = memitem->next;
						FREE(memitem);
						memitem = memlist;
				}
				pthread_mutex_unlock(&mutexmemlist);
		}

		/****************************************************************************/
		else if (strcasecmp(command, "info")==0) {

				/* include the current memory use in the list */
				memprofile_sample();

				pthread_mutex_lock(&mutexmemlist);
				/* count the number of items on the list */
				i = 0;
				memitem = memlist;
				while (memitem) {
						i++;
						memitem = memitem->next ;
				}

				plhs[0] = mxCreateStructMatrix(i, 1, numfields, fieldnames);

				/* return the items in a Matlab structure array */
				memitem = memlist;
				while (memitem) {
						i--;
						mxSetFieldByNumber(plhs[0], i, 0, mxCreateDoubleScalar((uint64_t)(memitem->time)));
						mxSetFieldByNumber(plhs[0], i, 1, mxCreateDoubleScalar((uint64_t)(memitem->rss)));
						memitem = memitem->next ;
				}
				pthread_mutex_unlock(&mutexmemlist);
		}

		/****************************************************************************/
		else if (strcasecmp(command, "report")==0) {
				uint64_t begmem, endmem, minmem, maxmem;
				time_t mintime, maxtime;
				int  num = 0;
				float summem = 0;
				memprofile_sample();

				/* include the current memory use in the list */
				memprofile_sample();

				pthread_mutex_lock(&mutexmemlist);
				if (memlist==NULL) {
						pthread_mutex_unlock(&mutexmemlist);
						mexErrMsgTxt ("memory profiling is disabled");
				}
				memitem = memlist;
				endmem = memitem->rss;
				minmem = memitem->rss;
				maxmem = memitem->rss;
				mintime = memitem->time;
				maxtime = memitem->time;
				while (memitem) {
						begmem  = memitem->rss; /* this will eventually contain the latest value on the list */
						minmem  = (memitem->rss < minmem ? memitem->rss : minmem);
						maxmem  = (memitem->rss > maxmem ? memitem->rss : maxmem);
						mintime = (memitem->time < mintime ? memitem->time : mintime);
						maxtime = (memitem->time > maxtime ? memitem->time : maxtime);
						summem += memitem->rss;
						num++;
						memitem = memitem->next ;
				}
				pthread_mutex_unlock(&mutexmemlist);
				mexPrintf("duration of the recording  = %6u s\n", (int)(maxtime-mintime+1));
				mexPrintf("memory in use at the begin = %6u MB\n", begmem/1048576);
				mexPrintf("memory in use at the end   = %6u MB\n", endmem/1048576);
				mexPrintf("minimum memory in use      = %6u MB\n", minmem/1048576);
				mexPrintf("maximum memory in use      = %6u MB\n", maxmem/1048576);
				mexPrintf("average memory in use      = %6u MB\n", (int)(summem/num)/1048576);
		}

		/****************************************************************************/
		else {
				mexErrMsgTxt ("unknown input option");
		}

		return;
}

