/*
 * Copyright (C) 2010, Robert Oostenveld
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/
 *
 */

#include "mex.h"
#include "matrix.h"
#include <pthread.h>

#include "peer.h"
#include "extern.h"
#include "platform_includes.h"

mxArray *mxSerialize(const mxArray*);
mxArray *mxDeserialize(const void*, size_t);

#define JOB_FIELDNUMBER 4
const char* job_fieldnames[JOB_FIELDNUMBER] = {"version", "jobid", "argsize", "optsize"};

#define JOBLIST_FIELDNUMBER 9
const char* joblist_fieldnames[JOBLIST_FIELDNUMBER] = {"version", "jobid", "argsize", "optsize", "hostid", "hostname", "user", "memreq", "timreq"};

#define PEERINFO_FIELDNUMBER 16
const char* peerinfo_fieldnames[PEERINFO_FIELDNUMBER] = {"hostid", "hostname", "user", "group", "socket", "port", "status", "timavail", "memavail", "cpuavail", "allowuser", "allowgroup", "allowhost", "refuseuser", "refusegroup", "refusehost"};

#define PEERLIST_FIELDNUMBER 11
const char* peerlist_fieldnames[PEERLIST_FIELDNUMBER] = {"hostid", "hostname", "user", "group", "socket", "port", "status", "timavail", "memavail", "cpuavail", "current"};

#define CURRENT_FIELDNUMBER 8
const char* current_fieldnames[CURRENT_FIELDNUMBER] = {"hostid", "jobid", "hostname", "user", "group", "timreq", "memreq", "cpureq"};

int peerInitialized = 0;

/* the thread IDs are needed for cancelation at cleanup */
pthread_t udsserverThread;
pthread_t tcpserverThread;
pthread_t announceThread;
pthread_t discoverThread;
pthread_t expireThread;

/* this is called the first time that the mex-file is loaded */
void initFun(void) {
		/* check whether the host has been initialized already */
		if (!peerInitialized) {

				/* open the logging facility and set the default level */
				#if SYSLOG ==1
				openlog("peer.mex", LOG_PID, LOG_USER);
				setlogmask(LOG_MASK(LOG_EMERG) | LOG_MASK(LOG_ALERT) | LOG_MASK(LOG_CRIT));
				#endif
				mexPrintf("peer: init\n");
				peerinit(NULL);
				peerInitialized = 1;
		}
		return;
}

/* this function is called upon unloading of the mex-file */
void exitFun(void) {
		mexPrintf("peer: exit\n");

		/* tell all threads to stop running */

		pthread_mutex_lock(&mutexstatus);
		if (udsserverStatus) {
				pthread_mutex_unlock(&mutexstatus);
				mexPrintf("peer: requesting cancelation of udsserver thread\n");
				pthread_cancel(udsserverThread);
				pthread_join(udsserverThread, NULL);
		}
		else {
				pthread_mutex_unlock(&mutexstatus);
		}

		pthread_mutex_lock(&mutexstatus);
		if (tcpserverStatus) {
				pthread_mutex_unlock(&mutexstatus);
				mexPrintf("peer: requesting cancelation of tcpserver thread\n");
				pthread_cancel(tcpserverThread);
				pthread_join(tcpserverThread, NULL);
		}
		else {
				pthread_mutex_unlock(&mutexstatus);
		}

		pthread_mutex_lock(&mutexstatus);
		if (announceStatus) {
				pthread_mutex_unlock(&mutexstatus);
				mexPrintf("peer: requesting cancelation of announce thread\n");
				pthread_cancel(announceThread);
				pthread_join(announceThread, NULL);
		}
		else {
				pthread_mutex_unlock(&mutexstatus);
		}

		pthread_mutex_lock(&mutexstatus);
		if (discoverStatus) {
				pthread_mutex_unlock(&mutexstatus);
				mexPrintf("peer: requesting cancelation of discover thread\n");
				pthread_cancel(discoverThread);
				pthread_join(discoverThread, NULL);
		}
		else {
				pthread_mutex_unlock(&mutexstatus);
		}

		pthread_mutex_lock(&mutexstatus);
		if (expireStatus) {
				pthread_mutex_unlock(&mutexstatus);
				mexPrintf("peer: requesting cancelation of expire thread\n");
				pthread_cancel(expireThread);
				pthread_join(expireThread, NULL);
		}
		else {
				pthread_mutex_unlock(&mutexstatus);
		}

		/* free the shared dynamical memory */
		peerexit(NULL);
		peerInitialized = 0;

		/* switch off the logging facility */
		#if SYSLOG == 1
		closelog();
		#endif

		return;
}

void mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]) {
		char command[STRLEN];
		char argument[STRLEN];
		char *ptr;
		int i, j, n, rc, t, found, handshake, success, count, server, status;
		UINT32_T peerid, jobid;
		UINT64_T memreq, cpureq, timreq;

		jobdef_t    *def;
		joblist_t   *job, *nextjob;
		peerlist_t  *peer;
		userlist_t  *allowuser, *refuseuser;
		grouplist_t *allowgroup, *refusegroup;
		hostlist_t  *allowhost, *refusehost;
		mxArray     *arg, *opt, *key, *val, *current;

		initFun();
		mexAtExit(exitFun);

		/* the first argument is always the command string */
		if (nrhs<1)
				mexErrMsgTxt ("invalid number of input arguments");

		if (!mxIsChar(prhs[0]))
				mexErrMsgTxt ("invalid input argument #1");
		if (mxGetString(prhs[0], command, STRLEN-1))
				mexErrMsgTxt ("invalid input argument #1");

		/* convert to lower case */
		ptr = command;
		while (*ptr) {
				*ptr = tolower(*ptr);
				ptr++;
		}

		ptr = argument;
		while (*ptr) {
				*ptr = tolower(*ptr);
				ptr++;
		}

		/****************************************************************************/
		if (strcmp(command, "announce")==0) {
				/* the input arguments should be "tcpserver <start|stop|status>" */
				if (nrhs<2)
						mexErrMsgTxt ("invalid number of input arguments");
				if (!mxIsChar(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				if (mxGetString(prhs[1], argument, STRLEN-1))
						mexErrMsgTxt ("invalid input argument #2");
				if (strcmp(argument, "start")==0) {
						if (announceStatus) {
								mexWarnMsgTxt("thread is already running");
								return;
						}
						mexPrintf("peer: spawning announce thread\n");
						rc = pthread_create(&announceThread, NULL, announce, (void *)NULL);
						if (rc)
								mexErrMsgTxt("problem with return code from pthread_create()");
						else {
								/* wait until the thread has properly started */
								pthread_mutex_lock(&mutexstatus);
								if (!announceStatus)
										pthread_cond_wait(&condstatus, &mutexstatus);
								pthread_mutex_unlock(&mutexstatus);
						}
				}
				else if (strcmp(argument, "stop")==0) {
						if (!announceStatus) {
								mexWarnMsgTxt("thread is not running");
								return;
						}
						mexPrintf("peer: requesting cancelation of announce thread\n");
						rc = pthread_cancel(announceThread);
						if (rc)
								mexErrMsgTxt("problem with return code from pthread_cancel()");
				}
				else if (strcmp(argument, "status")==0) {
						plhs[0] = mxCreateDoubleScalar(announceStatus);
				}
				else
						mexErrMsgTxt ("invalid input argument #2");
				return;
		}

		/****************************************************************************/
		else if (strcmp(command, "discover")==0) {
				/* the input arguments should be "discover <start|stop|status>" */
				if (nrhs<2)
						mexErrMsgTxt ("invalid number of input arguments");
				if (!mxIsChar(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				if (mxGetString(prhs[1], argument, STRLEN-1))
						mexErrMsgTxt ("invalid input argument #2");
				if (strcmp(argument, "start")==0) {
						if (discoverStatus) {
								mexWarnMsgTxt("thread is already running");
								return;
						}
						mexPrintf("peer: spawning discover thread\n");
						rc = pthread_create(&discoverThread, NULL, discover, (void *)NULL);
						if (rc)
								mexErrMsgTxt("problem with return code from pthread_create()");
						else {
								/* wait until the thread has properly started */
								pthread_mutex_lock(&mutexstatus);
								if (!discoverStatus)
										pthread_cond_wait(&condstatus, &mutexstatus);
								pthread_mutex_unlock(&mutexstatus);
						}
				}
				else if (strcmp(argument, "stop")==0) {
						if (!discoverStatus) {
								mexWarnMsgTxt("thread is not running");
								return;
						}
						mexPrintf("peer: requesting cancelation of discover thread\n");
						rc = pthread_cancel(discoverThread);
						if (rc)
								mexErrMsgTxt("problem with return code from pthread_cancel()");
				}
				else if (strcmp(argument, "status")==0) {
						plhs[0] = mxCreateDoubleScalar(discoverStatus);
				}
				else
						mexErrMsgTxt ("invalid input argument #2");
				return;
		}

		/****************************************************************************/
		else if (strcmp(command, "expire")==0) {
				/* the input arguments should be "expire <start|stop|status>" */
				if (nrhs<2)
						mexErrMsgTxt ("invalid number of input arguments");
				if (!mxIsChar(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				if (mxGetString(prhs[1], argument, STRLEN-1))
						mexErrMsgTxt ("invalid input argument #2");
				if (strcmp(argument, "start")==0) {
						if (expireStatus) {
								mexWarnMsgTxt("thread is already running");
								return;
						}
						mexPrintf("peer: spawning expire thread\n");
						rc = pthread_create(&expireThread, NULL, expire, (void *)NULL);
						if (rc)
								mexErrMsgTxt("problem with return code from pthread_create()");
						else {
								/* wait until the thread has properly started */
								pthread_mutex_lock(&mutexstatus);
								if (!expireStatus)
										pthread_cond_wait(&condstatus, &mutexstatus);
								pthread_mutex_unlock(&mutexstatus);
						}
				}
				else if (strcmp(argument, "stop")==0) {
						if (!expireStatus) {
								mexWarnMsgTxt("thread is not running");
								return;
						}
						mexPrintf("peer: requesting cancelation of expire thread\n");
						rc = pthread_cancel(expireThread);
						if (rc)
								mexErrMsgTxt("problem with return code from pthread_cancel()");
				}
				else if (strcmp(argument, "status")==0) {
						plhs[0] = mxCreateDoubleScalar(expireStatus);
				}
				else
						mexErrMsgTxt ("invalid input argument #2");
				return;
		}

		/****************************************************************************/
		else if (strcmp(command, "udsserver")==0) {
				/* the input arguments should be "udsserver <start|stop|status>" */
				if (nrhs<2)
						mexErrMsgTxt ("invalid number of input arguments");
				if (!mxIsChar(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				if (mxGetString(prhs[1], argument, STRLEN-1))
						mexErrMsgTxt ("invalid input argument #2");
				if (strcmp(argument, "start")==0) {
						if (udsserverStatus) {
								mexWarnMsgTxt("thread is already running");
								return;
						}
						mexPrintf("peer: spawning udsserver thread\n");
						rc = pthread_create(&udsserverThread, NULL, udsserver, (void *)NULL);
						if (rc)
								mexErrMsgTxt("problem with return code from pthread_create()");
						else {
								/* wait until the thread has properly started */
								pthread_mutex_lock(&mutexstatus);
								if (!udsserverStatus)
										pthread_cond_wait(&condstatus, &mutexstatus);
								pthread_mutex_unlock(&mutexstatus);
								/* inform the other peers of the updated unix domain socket */
								announce_once();
						}
				}
				else if (strcmp(argument, "stop")==0) {
						if (!udsserverStatus) {
								mexWarnMsgTxt("thread is not running");
								return;
						}
						mexPrintf("peer: requesting cancelation of udsserver thread\n");
						rc = pthread_cancel(udsserverThread);
						if (rc)
								mexErrMsgTxt("problem with return code from pthread_cancel()");
						else
								/* inform the other peers of the updated unix domain socket */
								announce_once();
				}
				else if (strcmp(argument, "status")==0) {
						plhs[0] = mxCreateDoubleScalar(udsserverStatus);
				}
				else
						mexErrMsgTxt ("invalid input argument #2");
				return;
		}

		/****************************************************************************/
		if (strcmp(command, "tcpserver")==0) {
				/* the input arguments should be "tcpserver <start|stop|status>" */
				if (nrhs<2)
						mexErrMsgTxt ("invalid number of input arguments");
				if (!mxIsChar(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				if (mxGetString(prhs[1], argument, STRLEN-1))
						mexErrMsgTxt ("invalid input argument #2");
				if (strcmp(argument, "start")==0) {
						if (tcpserverStatus) {
								mexWarnMsgTxt("thread is already running");
								return;
						}
						mexPrintf("peer: spawning tcpserver thread\n");
						rc = pthread_create(&tcpserverThread, NULL, tcpserver, (void *)NULL);
						if (rc)
								mexErrMsgTxt("problem with return code from pthread_create()");
						else {
								/* wait until the thread has properly started */
								pthread_mutex_lock(&mutexstatus);
								if (!tcpserverStatus)
										pthread_cond_wait(&condstatus, &mutexstatus);
								pthread_mutex_unlock(&mutexstatus);
								/* inform the other peers of the tcp server status */
								announce_once();
						}
				}
				else if (strcmp(argument, "stop")==0) {
						if (!tcpserverStatus) {
								mexWarnMsgTxt("thread is not running");
								return;
						}
						mexPrintf("peer: requesting cancelation of tcpserver thread\n");
						rc = pthread_cancel(tcpserverThread);
						if (rc)
								mexErrMsgTxt("problem with return code from pthread_cancel()");
						else
								/* inform the other peers of the updated tcp server status */
								announce_once();
				}
				else if (strcmp(argument, "status")==0) {
						plhs[0] = mxCreateDoubleScalar(tcpserverStatus);
				}
				else
						mexErrMsgTxt ("invalid input argument #2");
				return;
		}

		/****************************************************************************/
		else if (strcmp(command, "peerlist")==0) {
				/* it can be followed by an optional input argument "peerlist <status>" */
				if (nrhs>1) {
						/* only return the list for the desired status */
						if (!mxIsNumeric(prhs[1]))
								mexErrMsgTxt ("invalid input argument #2");
						else
								status = (UINT32_T)mxGetScalar(prhs[1]);
				}
				else {
						/* this is to be interpreted as "any status" */
						status = -1;
				}

				/* count the number of peers */
				i = 0;
				pthread_mutex_lock(&mutexpeerlist);
				peer = peerlist;
				while(peer) {
						i += (status==-1 ? 1 : (peer->host->status==status));
						peer = peer->next ;
				}
				plhs[0] = mxCreateStructMatrix(i, 1, PEERLIST_FIELDNUMBER, peerlist_fieldnames);
				i = 0;
				peer = peerlist;
				while(peer) {
						if (status!=-1 && peer->host->status!=status) {
								peer = peer->next;
								continue;
						}
						j = 0;
						mxSetFieldByNumber(plhs[0], i, j++, mxCreateDoubleScalar((UINT32_T)(peer->host->id)));
						mxSetFieldByNumber(plhs[0], i, j++, mxCreateString(peer->host->name));
						mxSetFieldByNumber(plhs[0], i, j++, mxCreateString(peer->host->user));
						mxSetFieldByNumber(plhs[0], i, j++, mxCreateString(peer->host->group));
						mxSetFieldByNumber(plhs[0], i, j++, mxCreateString(peer->host->socket));
						mxSetFieldByNumber(plhs[0], i, j++, mxCreateDoubleScalar((UINT32_T)(peer->host->port)));
						mxSetFieldByNumber(plhs[0], i, j++, mxCreateDoubleScalar((UINT32_T)(peer->host->status)));
						mxSetFieldByNumber(plhs[0], i, j++, mxCreateDoubleScalar((UINT64_T)(peer->host->timavail)));
						mxSetFieldByNumber(plhs[0], i, j++, mxCreateDoubleScalar((UINT64_T)(peer->host->memavail)));
						mxSetFieldByNumber(plhs[0], i, j++, mxCreateDoubleScalar((UINT64_T)(peer->host->cpuavail)));

						current = mxCreateStructMatrix(1, 1, CURRENT_FIELDNUMBER, current_fieldnames);
						mxSetFieldByNumber(current, 0, 0, mxCreateDoubleScalar(peer->host->current.hostid));
						mxSetFieldByNumber(current, 0, 1, mxCreateDoubleScalar(peer->host->current.jobid));
						mxSetFieldByNumber(current, 0, 2, mxCreateString(peer->host->current.name));
						mxSetFieldByNumber(current, 0, 3, mxCreateString(peer->host->current.user));
						mxSetFieldByNumber(current, 0, 4, mxCreateString(peer->host->current.group));
						mxSetFieldByNumber(current, 0, 5, mxCreateDoubleScalar(peer->host->current.timreq));
						mxSetFieldByNumber(current, 0, 6, mxCreateDoubleScalar(peer->host->current.memreq));
						mxSetFieldByNumber(current, 0, 7, mxCreateDoubleScalar(peer->host->current.cpureq));
						mxSetFieldByNumber(plhs[0], i, j++, current);

						i++;
						peer = peer->next ;
				}
				pthread_mutex_unlock(&mutexpeerlist);
				return;
		}

		/****************************************************************************/
		else if (strcmp(command, "joblist")==0) {
				/* count the number of jobs */
				i = 0;
				pthread_mutex_lock(&mutexjoblist);
				job = joblist;
				while(job) {
						i++;
						job = job->next ;
				}
				plhs[0] = mxCreateStructMatrix(i, 1, JOBLIST_FIELDNUMBER, joblist_fieldnames);
				i = 0;
				job = joblist;
				while(job) {
						mxSetFieldByNumber(plhs[0], i, 0, mxCreateDoubleScalar((UINT32_T)(job->job->version)));
						mxSetFieldByNumber(plhs[0], i, 1, mxCreateDoubleScalar((UINT32_T)(job->job->id)));
						mxSetFieldByNumber(plhs[0], i, 2, mxCreateDoubleScalar((UINT32_T)(job->job->argsize)));
						mxSetFieldByNumber(plhs[0], i, 3, mxCreateDoubleScalar((UINT32_T)(job->job->optsize)));
						mxSetFieldByNumber(plhs[0], i, 4, mxCreateDoubleScalar((UINT32_T)(job->host->id)));
						mxSetFieldByNumber(plhs[0], i, 5, mxCreateString(job->host->name));
						mxSetFieldByNumber(plhs[0], i, 6, mxCreateString(job->host->user));
						mxSetFieldByNumber(plhs[0], i, 7, mxCreateDoubleScalar((UINT64_T)(job->job->memreq)));
						mxSetFieldByNumber(plhs[0], i, 8, mxCreateDoubleScalar((UINT64_T)(job->job->timreq)));
						job = job->next;
						i++;
				}
				pthread_mutex_unlock(&mutexjoblist);
				return;
		}

		/****************************************************************************/
		else if (strcmp(command, "status")==0) {
				/* the input arguments should be "status <number>" */
				if (nrhs<2) {
						mexErrMsgTxt ("invalid number of input arguments");
				}
				else {
						if (!mxIsNumeric(prhs[1]))
								mexErrMsgTxt ("invalid input argument #2");
						if (!mxIsScalar(prhs[1]))
								mexErrMsgTxt ("invalid input argument #2");

						pthread_mutex_lock(&mutexhost);
						host->status = (UINT32_T)mxGetScalar(prhs[1]);
						pthread_mutex_unlock(&mutexhost);
						/* inform the other peers of the updated status */
						announce_once();
				}
		}

		/****************************************************************************/
		else if (strcmp(command, "info")==0) {
				/* the input arguments should be "info" */

				/* display  all possible information on screen */
				pthread_mutex_lock(&mutexstatus);
				mexPrintf("tcpserverStatus = %d\n", tcpserverStatus);
				mexPrintf("udsserverStatus = %d\n", udsserverStatus);
				mexPrintf("announceStatus  = %d\n", announceStatus);
				mexPrintf("discoverStatus  = %d\n", discoverStatus);
				mexPrintf("expireStatus    = %d\n", expireStatus);
				pthread_mutex_unlock(&mutexstatus);

				pthread_mutex_lock(&mutexhost);
				mexPrintf("host.id         = %u\n", host->id);
				mexPrintf("host.name       = %s\n", host->name);
				mexPrintf("host.port       = %u\n", host->port);
				mexPrintf("host.user       = %s\n", host->user);
				mexPrintf("host.group      = %s\n", host->group);
				mexPrintf("host.status     = %u\n", host->status);
				mexPrintf("host.timavail   = %llu\n", host->timavail);
				mexPrintf("host.memavail   = %llu\n", host->memavail);
				mexPrintf("host.cpuavail   = %llu\n", host->cpuavail);
				pthread_mutex_unlock(&mutexhost);

				pthread_mutex_lock(&mutexsmartmem);
				mexPrintf("smartmem.enabled = %d\n", smartmem.enabled);
				pthread_mutex_unlock(&mutexsmartmem);

				pthread_mutex_lock(&mutexsmartcpu);
				mexPrintf("smartcpu.enabled = %d\n", smartcpu.enabled);
				pthread_mutex_unlock(&mutexsmartcpu);

				pthread_mutex_lock(&mutexsmartshare);
				mexPrintf("smartshare.enabled = %d\n", smartshare.enabled);
				pthread_mutex_unlock(&mutexsmartshare);

				pthread_mutex_lock(&mutexallowuserlist);
				allowuser = allowuserlist;
				while (allowuser) {
						mexPrintf("allowuser = %s\n", allowuser->name);
						allowuser = allowuser->next;
				}
				pthread_mutex_unlock(&mutexallowuserlist);

				pthread_mutex_lock(&mutexallowgrouplist);
				allowgroup = allowgrouplist;
				while (allowgroup) {
						mexPrintf("allowgroup = %s\n", allowgroup->name);
						allowgroup = allowgroup->next;
				}
				pthread_mutex_unlock(&mutexallowgrouplist);

				pthread_mutex_lock(&mutexallowhostlist);
				allowhost = allowhostlist;
				while (allowhost) {
						mexPrintf("allowhost = %s\n", allowhost->name);
						allowhost = allowhost->next;
				}
				pthread_mutex_unlock(&mutexallowhostlist);

				pthread_mutex_lock(&mutexrefuseuserlist);
				refuseuser = refuseuserlist;
				while (refuseuser) {
						mexPrintf("refuseuser = %s\n", refuseuser->name);
						refuseuser = refuseuser->next;
				}
				pthread_mutex_unlock(&mutexrefuseuserlist);

				pthread_mutex_lock(&mutexrefusegrouplist);
				refusegroup = refusegrouplist;
				while (refusegroup) {
						mexPrintf("refusegroup = %s\n", refusegroup->name);
						refusegroup = refusegroup->next;
				}
				pthread_mutex_unlock(&mutexrefusegrouplist);

				pthread_mutex_lock(&mutexrefusehostlist);
				refusehost = refusehostlist;
				while (refusehost) {
						mexPrintf("refusehost = %s\n", refusehost->name);
						refusehost = refusehost->next;
				}
				pthread_mutex_unlock(&mutexrefusehostlist);

				i = 0;
				pthread_mutex_lock(&mutexpeerlist);
				peer = peerlist;
				while(peer) {
						mexPrintf("peerlist[%d] = \n", i);
						mexPrintf("  host.id       = %u\n", peer->host->id);
						mexPrintf("  host.name     = %s\n", peer->host->name);
						mexPrintf("  host.port     = %u\n", peer->host->port);
						mexPrintf("  host.user     = %s\n", peer->host->user);
						mexPrintf("  host.group    = %s\n", peer->host->group);
						mexPrintf("  host.status   = %u\n", peer->host->status);
						mexPrintf("  host.timavail = %llu\n", peer->host->timavail);
						mexPrintf("  host.memavail = %llu\n", peer->host->memavail);
						mexPrintf("  host.cpuavail = %llu\n", peer->host->cpuavail);
						mexPrintf("  ipaddr        = %s\n", peer->ipaddr);
						mexPrintf("  time          = %s",   ctime(&(peer->time)));
						peer = peer->next ;       
						i++;
				}
				pthread_mutex_unlock(&mutexpeerlist);

				i = 0;
				pthread_mutex_lock(&mutexjoblist);
				job = joblist;
				while(job) {
						mexPrintf("joblist[%d] = \n", i);
						mexPrintf("  job.version = %u\n", job->job->version);
						mexPrintf("  job.id      = %u\n", job->job->id);
						mexPrintf("  job.memreq  = %llu\n", job->job->memreq);
						mexPrintf("  job.cpureq  = %llu\n", job->job->cpureq);
						mexPrintf("  job.timreq  = %llu\n", job->job->timreq);
						mexPrintf("  job.argsize = %u\n", job->job->argsize);
						mexPrintf("  job.optsize = %u\n", job->job->optsize);
						mexPrintf("  job.host.id = %u\n", job->host->id);
						job = job->next ;
						i++;
				}
				pthread_mutex_unlock(&mutexjoblist);
				return;
		}

		/****************************************************************************/
		else if (strcmp(command, "smartmem")==0) {
				/* the input arguments should be "smartmem <0|1> */
				if (nrhs<2) {
						mexErrMsgTxt ("invalid number of input arguments");
				}
				else {

						if (!mxIsNumeric(prhs[1]))
								mexErrMsgTxt ("invalid input argument #2");
						if (!mxIsScalar(prhs[1]))
								mexErrMsgTxt ("invalid input argument #2");

						pthread_mutex_lock(&mutexsmartmem);
						smartmem.enabled = mxGetScalar(prhs[1]);
						pthread_mutex_unlock(&mutexsmartmem);
				}
		}

		/****************************************************************************/
		else if (strcmp(command, "smartcpu")==0) {
				/* the input arguments should be "smartcpu <0|1> */
				if (nrhs<2) {
						mexErrMsgTxt ("invalid number of input arguments");
				}
				else {

						if (!mxIsNumeric(prhs[1]))
								mexErrMsgTxt ("invalid input argument #2");
						if (!mxIsScalar(prhs[1]))
								mexErrMsgTxt ("invalid input argument #2");

						pthread_mutex_lock(&mutexsmartcpu);
						smartcpu.enabled = mxGetScalar(prhs[1]);
						pthread_mutex_unlock(&mutexsmartcpu);
				}
		}

		/****************************************************************************/
		else if (strcmp(command, "smartshare")==0) {
				/* the input arguments should be "smartshare <0|1> */
				if (nrhs<2) {
						mexErrMsgTxt ("invalid number of input arguments");
				}
				else {

						if (!mxIsNumeric(prhs[1]))
								mexErrMsgTxt ("invalid input argument #2");
						if (!mxIsScalar(prhs[1]))
								mexErrMsgTxt ("invalid input argument #2");

						pthread_mutex_lock(&mutexsmartshare);
						smartshare.enabled = mxGetScalar(prhs[1]);
						pthread_mutex_unlock(&mutexsmartshare);
				}
		}

		/****************************************************************************/
		else if (strcmp(command, "memavail")==0) {
				/* the input arguments should be "memavail <number>" */
				if (nrhs<2) {
						mexErrMsgTxt ("invalid number of input arguments");
				}
				else {
						if (!mxIsNumeric(prhs[1]))
								mexErrMsgTxt ("invalid input argument #2");
						if (!mxIsScalar(prhs[1]))
								mexErrMsgTxt ("invalid input argument #2");

						pthread_mutex_lock(&mutexhost);
						host->memavail = (UINT64_T)(mxGetScalar(prhs[1])+0.5);
						pthread_mutex_unlock(&mutexhost);
				}
		}

		/****************************************************************************/
		else if (strcmp(command, "cpuavail")==0) {
				/* the input arguments should be "cpuavail <number>" */
				if (nrhs<2) {
						mexErrMsgTxt ("invalid number of input arguments");
				}
				else {
						if (!mxIsNumeric(prhs[1]))
								mexErrMsgTxt ("invalid input argument #2");
						if (!mxIsScalar(prhs[1]))
								mexErrMsgTxt ("invalid input argument #2");

						pthread_mutex_lock(&mutexhost);
						host->cpuavail = (UINT64_T)(mxGetScalar(prhs[1])+0.5);
						pthread_mutex_unlock(&mutexhost);
				}
		}

		/****************************************************************************/
		else if (strcmp(command, "timavail")==0) {
				/* the input arguments should be "timavail <number>" */
				if (nrhs<2) {
						mexErrMsgTxt ("invalid number of input arguments");
				}
				else {
						if (!mxIsNumeric(prhs[1]))
								mexErrMsgTxt ("invalid input argument #2");
						if (!mxIsScalar(prhs[1]))
								mexErrMsgTxt ("invalid input argument #2");

						pthread_mutex_lock(&mutexhost);
						host->timavail = (UINT64_T)(mxGetScalar(prhs[1])+0.5);
						pthread_mutex_unlock(&mutexhost);
				}
		}

		/****************************************************************************/
		else if (strcmp(command, "tcpport")==0) {
				/* the input arguments should be "tcpport <number>" */
				if (nrhs<2)
						mexErrMsgTxt ("invalid number of input arguments");
				if (!mxIsNumeric(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				if (!mxIsScalar(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				if (tcpserverStatus) 
						mexErrMsgTxt ("cannot change the port while the tcpserver is running");
				pthread_mutex_lock(&mutexhost);
				host->port = (UINT32_T)mxGetScalar(prhs[1]);
				pthread_mutex_unlock(&mutexhost);
		}

		/****************************************************************************/
		else if (strcmp(command, "group")==0) {
				/* the input arguments should be "group <string>" */
				if (nrhs<2)
						mexErrMsgTxt ("invalid number of input arguments");

				pthread_mutex_lock(&mutexhost);
				if (mxIsEmpty(prhs[1])) {
						/* set the default group name */
						strncpy(host->group, DEFAULT_GROUP, STRLEN);
				}
				else {
						if (!mxIsChar(prhs[1]))
								/* FIXME this causes a deadlock */
								mexErrMsgTxt ("invalid input argument #2");
						if (mxGetString(prhs[1], host->group, STRLEN))
								/* FIXME this causes a deadlock */
								mexErrMsgTxt("FIXME: unexpected error");
				}
				pthread_mutex_unlock(&mutexhost);
		}

		/****************************************************************************/
		else if (strcmp(command, "hostname")==0) {
				/* the input arguments should be "group <string>" */
				if (nrhs<2)
						mexErrMsgTxt ("invalid number of input arguments");

				pthread_mutex_lock(&mutexhost);
				if (mxIsEmpty(prhs[1])) {
						/* set the default, i.e. the hostname of this computer */
						if (gethostname(host->name, STRLEN)) {
								pthread_mutex_unlock(&mutexhost);
								mexErrMsgTxt("could not get hostname");
						}
				}
				else {
						/* use the hostname specified by the user, e.g. localhost */
						if (!mxIsChar(prhs[1])) {
								pthread_mutex_unlock(&mutexhost);
								mexErrMsgTxt ("invalid input argument #2 for hostname");
						}
						if (mxGetString(prhs[1], host->name, STRLEN)) {
								pthread_mutex_unlock(&mutexhost);
								mexErrMsgTxt("could not copy string for hostname");
						}
				}
				/* the hostname was correctly updated */
				pthread_mutex_unlock(&mutexhost);
		}

		/****************************************************************************/
		else if (strcmp(command, "allowuser")==0) {
				/* the input arguments should be "allowuser {<string>, <string>, ...}" */
				if (nrhs<2)
						mexErrMsgTxt ("invalid number of input arguments");
				if (!mxIsCell(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				/* check that all elements of the cell array are strings */
				for (i=0; i<mxGetNumberOfElements(prhs[1]); i++) {
						arg = mxGetCell(prhs[1], i);
						if (!mxIsChar(arg))
								mexErrMsgTxt ("invalid input argument #2, the cell-array should contain strings");
				}

				/* erase the existing list */
				clear_allowuserlist();

				/* add all elements to the list */
				pthread_mutex_lock(&mutexallowuserlist);
				for (i=0; i<mxGetNumberOfElements(prhs[1]); i++) {
						arg = mxGetCell(prhs[1], i);
						allowuser = (userlist_t *)malloc(sizeof(userlist_t));
						allowuser->name = malloc(mxGetNumberOfElements(arg)+1);
						allowuser->next = allowuserlist;
						if(mxGetString(arg, allowuser->name, mxGetNumberOfElements(arg)+1)!=0)
								mexErrMsgTxt("FIXME: unexpected error, memory leak");
						allowuserlist = allowuser;
				}
				pthread_mutex_unlock(&mutexallowuserlist);

				/* erase the list of known peers, the updated filtering will be done upon discovery */
				clear_peerlist();
		}

		/****************************************************************************/
		else if (strcmp(command, "allowgroup")==0) {
				/* the input arguments should be "allowgroup {<string>, <string>, ...}" */
				if (nrhs<2)
						mexErrMsgTxt ("invalid number of input arguments");
				if (!mxIsCell(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				/* check that all elements of the cell array are strings */
				for (i=0; i<mxGetNumberOfElements(prhs[1]); i++) {
						arg = mxGetCell(prhs[1], i);
						if (!mxIsChar(arg))
								mexErrMsgTxt ("invalid input argument #2, the cell-array should contain strings");
				}

				/* erase the existing list */
				clear_allowgrouplist();

				/* add all elements to the list */
				pthread_mutex_lock(&mutexallowgrouplist);
				for (i=0; i<mxGetNumberOfElements(prhs[1]); i++) {
						arg = mxGetCell(prhs[1], i);
						allowgroup = (grouplist_t *)malloc(sizeof(grouplist_t));
						allowgroup->name = malloc(mxGetNumberOfElements(arg)+1);
						allowgroup->next = allowgrouplist;
						if(mxGetString(arg, allowgroup->name, mxGetNumberOfElements(arg)+1)!=0)
								mexErrMsgTxt("FIXME: unexpected error, memory leak");
						allowgrouplist = allowgroup;
				}
				pthread_mutex_unlock(&mutexallowgrouplist);

				/* flush the list of known peers, the updated filtering will be done upon discovery */
				clear_peerlist();
		}

		/****************************************************************************/
		else if (strcmp(command, "allowhost")==0) {
				/* the input arguments should be "allowhost {<string>, <string>, ...}" */
				if (nrhs<2)
						mexErrMsgTxt ("invalid number of input arguments");
				if (!mxIsCell(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				/* check that all elements of the cell array are strings */
				for (i=0; i<mxGetNumberOfElements(prhs[1]); i++) {
						arg = mxGetCell(prhs[1], i);
						if (!mxIsChar(arg))
								mexErrMsgTxt ("invalid input argument #2, the cell-array should contain strings");
				}

				/* erase the existing list */
				clear_allowhostlist();

				/* add all elements to the list */
				pthread_mutex_lock(&mutexallowhostlist);
				for (i=0; i<mxGetNumberOfElements(prhs[1]); i++) {
						arg = mxGetCell(prhs[1], i);
						allowhost = (hostlist_t *)malloc(sizeof(hostlist_t));
						allowhost->name = malloc(mxGetNumberOfElements(arg)+1);
						allowhost->next = allowhostlist;
						if(mxGetString(arg, allowhost->name, mxGetNumberOfElements(arg)+1)!=0)
								mexErrMsgTxt("FIXME: unexpected error, memory leak");
						allowhostlist = allowhost;
				}
				pthread_mutex_unlock(&mutexallowhostlist);

				/* flush the list of known peers, the updated filtering will be done upon discovery */
				clear_peerlist();
		}

		/****************************************************************************/
		else if (strcmp(command, "refuseuser")==0) {
				/* the input arguments should be "refuseuser {<string>, <string>, ...}" */
				if (nrhs<2)
						mexErrMsgTxt ("invalid number of input arguments");
				if (!mxIsCell(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				/* check that all elements of the cell array are strings */
				for (i=0; i<mxGetNumberOfElements(prhs[1]); i++) {
						arg = mxGetCell(prhs[1], i);
						if (!mxIsChar(arg))
								mexErrMsgTxt ("invalid input argument #2, the cell-array should contain strings");
				}

				/* erase the existing list */
				clear_refuseuserlist();

				/* add all elements to the list */
				pthread_mutex_lock(&mutexrefuseuserlist);
				for (i=0; i<mxGetNumberOfElements(prhs[1]); i++) {
						arg = mxGetCell(prhs[1], i);
						refuseuser = (userlist_t *)malloc(sizeof(userlist_t));
						refuseuser->name = malloc(mxGetNumberOfElements(arg)+1);
						refuseuser->next = refuseuserlist;
						if(mxGetString(arg, refuseuser->name, mxGetNumberOfElements(arg)+1)!=0)
								mexErrMsgTxt("FIXME: unexpected error, memory leak");
						refuseuserlist = refuseuser;
				}
				pthread_mutex_unlock(&mutexrefuseuserlist);

				/* erase the list of known peers, the updated filtering will be done upon discovery */
				clear_peerlist();
		}

		/****************************************************************************/
		else if (strcmp(command, "refusegroup")==0) {
				/* the input arguments should be "refusegroup {<string>, <string>, ...}" */
				if (nrhs<2)
						mexErrMsgTxt ("invalid number of input arguments");
				if (!mxIsCell(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				/* check that all elements of the cell array are strings */
				for (i=0; i<mxGetNumberOfElements(prhs[1]); i++) {
						arg = mxGetCell(prhs[1], i);
						if (!mxIsChar(arg))
								mexErrMsgTxt ("invalid input argument #2, the cell-array should contain strings");
				}

				/* erase the existing list */
				clear_refusegrouplist();

				/* add all elements to the list */
				pthread_mutex_lock(&mutexrefusegrouplist);
				for (i=0; i<mxGetNumberOfElements(prhs[1]); i++) {
						arg = mxGetCell(prhs[1], i);
						refusegroup = (grouplist_t *)malloc(sizeof(grouplist_t));
						refusegroup->name = malloc(mxGetNumberOfElements(arg)+1);
						refusegroup->next = refusegrouplist;
						if(mxGetString(arg, refusegroup->name, mxGetNumberOfElements(arg)+1)!=0)
								mexErrMsgTxt("FIXME: unexpected error, memory leak");
						refusegrouplist = refusegroup;
				}
				pthread_mutex_unlock(&mutexrefusegrouplist);

				/* flush the list of known peers, the updated filtering will be done upon discovery */
				clear_peerlist();
		}

		/****************************************************************************/
		else if (strcmp(command, "refusehost")==0) {
				/* the input arguments should be "refusehost {<string>, <string>, ...}" */
				if (nrhs<2)
						mexErrMsgTxt ("invalid number of input arguments");
				if (!mxIsCell(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				/* check that all elements of the cell array are strings */
				for (i=0; i<mxGetNumberOfElements(prhs[1]); i++) {
						arg = mxGetCell(prhs[1], i);
						if (!mxIsChar(arg))
								mexErrMsgTxt ("invalid input argument #2, the cell-array should contain strings");
				}

				/* erase the existing list */
				clear_refusehostlist();

				/* add all elements to the list */
				pthread_mutex_lock(&mutexrefusehostlist);
				for (i=0; i<mxGetNumberOfElements(prhs[1]); i++) {
						arg = mxGetCell(prhs[1], i);
						refusehost = (hostlist_t *)malloc(sizeof(hostlist_t));
						refusehost->name = malloc(mxGetNumberOfElements(arg)+1);
						refusehost->next = refusehostlist;
						if(mxGetString(arg, refusehost->name, mxGetNumberOfElements(arg)+1)!=0)
								mexErrMsgTxt("FIXME: unexpected error, memory leak");
						refusehostlist = refusehost;
				}
				pthread_mutex_unlock(&mutexrefusehostlist);

				/* flush the list of known peers, the updated filtering will be done upon discovery */
				clear_peerlist();
		}

		/****************************************************************************/
		else if (strcmp(command, "put")==0) {
				int hasuds, hastcp;
				/* the input arguments should be "put <peerid> <arg> <opt> ... "   */
				/* where additional options should be specified as key-value pairs */

				if (nrhs<2)
						mexErrMsgTxt("invalid argument #2");
				if (!mxIsNumeric(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				peerid = (UINT32_T)mxGetScalar(prhs[1]);

				if (nrhs<3)
						mexErrMsgTxt("invalid argument #3");

				if (nrhs<4)
						mexErrMsgTxt("invalid argument #4");

				jobid = rand();		/* assign a random jobid by default */
				memreq = 0; 		/* default assumption */
				cpureq = 0; 		/* default assumption */
				timreq = 0; 		/* default assumption */

				i = 4;
				while ((i+1)<nrhs) {
						key = (mxArray *)prhs[i++];
						val = (mxArray *)prhs[i++];

						if (!mxIsChar(key))
								mexErrMsgTxt ("optional arguments should come in key-value pairs");
						if (mxGetString(key, argument, STRLEN-1))
								mexErrMsgTxt ("optional arguments should come in key-value pairs");

						/* use the optional parameter value, beware of the typecasting to 32 and 64 bit integers */
						if      (strcmp(argument, "jobid")==0)
								jobid = (UINT32_T)mxGetScalar(val);
						else if (strcmp(argument, "memreq")==0)
								memreq = (UINT64_T)(mxGetScalar(val)+0.5);
						else if (strcmp(argument, "cpureq")==0)
								cpureq = (UINT64_T)(mxGetScalar(val)+0.5);
						else if (strcmp(argument, "timreq")==0)
								timreq = (UINT64_T)(mxGetScalar(val)+0.5);
				}

				found = 0;
				pthread_mutex_lock(&mutexpeerlist);
				peer = peerlist;
				while(peer) {
						if (peer->host->id==peerid) {
								found = 1;
								break;
						}
						peer = peer->next ;
				}

				if (!found) {
						pthread_mutex_unlock(&mutexpeerlist);
						mexErrMsgTxt("failed to locate specified peer\n");
				}

				arg = (mxArray *) mxSerialize(prhs[2]);
				if (!arg) {
						mexErrMsgTxt("could not serialize job arguments");
				}

				opt = (mxArray *) mxSerialize(prhs[3]);
				if (!opt) {
						mxDestroyArray(arg);
						arg = NULL;
						mexErrMsgTxt("could not serialize job options");
				}

				def = (jobdef_t *)malloc(sizeof(jobdef_t));
				if (!def) {
						mxDestroyArray(arg);
						arg = NULL;
						mxDestroyArray(opt);
						opt = NULL;
						mexErrMsgTxt("could not allocate memory");
				}

				pthread_mutex_lock(&mutexhost);
				hasuds = (strlen(peer->host->socket)>0 && strcmp(peer->host->name, host->name)==0);
				hastcp = (peer->host->port>0);
				pthread_mutex_unlock(&mutexhost);

				if (hasuds) {
						/* open the UDS socket */
						if ((server = open_uds_connection(peer->host->socket)) < 0) {
								pthread_mutex_unlock(&mutexpeerlist);
								mexErrMsgTxt("failed to create socket\n");
						}
				}
				else if (hastcp) {
						/* open the TCP socket */
						if ((server = open_tcp_connection(peer->ipaddr, peer->host->port)) < 0) {
								pthread_mutex_unlock(&mutexpeerlist);
								mexErrMsgTxt("failed to create socket\n");
						}
				}
				else {
						pthread_mutex_unlock(&mutexpeerlist);
						mexErrMsgTxt("failed to create socket\n");
				}

				/* the connection was opened without error */
				pthread_mutex_unlock(&mutexpeerlist);

				if ((n = bufread(server, &handshake, sizeof(int))) != sizeof(int)) {
						close_connection(server);
						mexErrMsgTxt("tcpsocket: could not write handshake");
				}

				if (!handshake) {
						close_connection(server);
						mexErrMsgTxt("failed to negociate connection");
				}

				/* the message that will be written consists of
				   message->host
				   message->job
				   message->arg
				   message->opt
				 */

				def->version  = VERSION;
				def->id       = jobid;
				def->memreq   = memreq;
				def->cpureq   = cpureq;
				def->timreq   = timreq;
				def->argsize  = mxGetNumberOfElements(arg);
				def->optsize  = mxGetNumberOfElements(opt);

				/* write the message  (hostdef, jobdef, arg, opt) with handshakes in between */
				/* the slave may close the connection between the message segments in case the job is refused */
				success = 1;

				pthread_mutex_lock(&mutexhost);
				if (success) 
						success = (bufwrite(server, host, sizeof(hostdef_t)) == sizeof(hostdef_t));
				pthread_mutex_unlock(&mutexhost);

				if ((n = bufread(server, &handshake, sizeof(int))) != sizeof(int)) {
						close_connection(server);
						mexErrMsgTxt("could not write handshake");
				}

				if (!handshake) {
						close_connection(server);
						mexErrMsgTxt("failed to write hostdef");
				}

				if (success)
						success = (bufwrite(server, def, sizeof(jobdef_t)) == sizeof(jobdef_t));

				if ((n = bufread(server, &handshake, sizeof(int))) != sizeof(int)) {
						close_connection(server);
						mexErrMsgTxt("could not write handshake");
				}

				if (!handshake) {
						close_connection(server);
						mexErrMsgTxt("failed to write jobdef");
				}

				if (success) 
						success = (bufwrite(server, (void *)mxGetData(arg), def->argsize) == def->argsize);

				if ((n = bufread(server, &handshake, sizeof(int))) != sizeof(int)) {
						close_connection(server);
						mexErrMsgTxt("could not write handshake");
				}

				if (!handshake) {
						close_connection(server);
						mexErrMsgTxt("failed to write arg");
				}

				if (success) 
						success = (bufwrite(server, (void *)mxGetData(opt), def->optsize) == def->optsize);

				if ((n = bufread(server, &handshake, sizeof(int))) != sizeof(int)) {
						close_connection(server);
						mexErrMsgTxt("could not write handshake");
				}

				if (!handshake) {
						close_connection(server);
						mexErrMsgTxt("failed to write opt");
				}

				close_connection(server);

				mxDestroyArray(arg);
				arg = NULL;
				mxDestroyArray(opt);
				opt = NULL;

				if (success) {
						/* return the job details */
						plhs[0] = mxCreateStructMatrix(1, 1, JOB_FIELDNUMBER, job_fieldnames);
						mxSetFieldByNumber(plhs[0], 0, 0, mxCreateDoubleScalar((UINT32_T)(def->version)));
						mxSetFieldByNumber(plhs[0], 0, 1, mxCreateDoubleScalar((UINT32_T)(def->id)));
						mxSetFieldByNumber(plhs[0], 0, 2, mxCreateDoubleScalar((UINT32_T)(def->argsize)));
						mxSetFieldByNumber(plhs[0], 0, 3, mxCreateDoubleScalar((UINT32_T)(def->optsize)));
						FREE(def);
				}
				else {
						FREE(def);
						mexErrMsgTxt ("failed to put job arguments");
				}

				return;
		}

		/****************************************************************************/
		else if (strcmp(command, "get")==0) {
				/* the input arguments should be "get <jobid>" */
				if (nrhs<2)
						mexErrMsgTxt ("invalid number of input arguments");
				if (!mxIsNumeric(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				if (!mxIsScalar(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				jobid = (UINT32_T)mxGetScalar(prhs[1]);

				found = 0;
				pthread_mutex_lock(&mutexjoblist);
				job = joblist;
				while(job) {
						found = (job->job->id==jobid);
						if (found) {
								plhs[0] = (mxArray *)mxDeserialize(job->arg, job->job->argsize);
								plhs[1] = (mxArray *)mxDeserialize(job->opt, job->job->optsize);
								break;
						}
						job = job->next ;
				}

				if (!found) {
						pthread_mutex_unlock(&mutexjoblist);
						mexErrMsgTxt("failed to locate specified job\n");
				}

				pthread_mutex_unlock(&mutexjoblist);
				return;
		}

		/****************************************************************************/
		else if (strcmp(command, "clear")==0) {
				/* the input arguments should be "clear <jobid>" */
				if (nrhs<2)
						mexErrMsgTxt ("invalid number of input arguments");
				if (!mxIsNumeric(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				if (!mxIsScalar(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				jobid = (UINT32_T)mxGetScalar(prhs[1]);

				found = 0;
				pthread_mutex_lock(&mutexjoblist);

				/* test the first item on the list */
				if (joblist)
				{
						if (joblist->job->id==jobid) {
								found = 1;
								/* delete the first item in the list */
								nextjob = joblist->next;
								FREE(joblist->job);
								FREE(joblist->host);
								FREE(joblist->arg);
								FREE(joblist->opt);
								FREE(joblist);
								joblist = nextjob;
						}
				}

				/* traverse the list */
				job = joblist;
				while(job) {
						/* test the next item on the list, remember the current one */
						nextjob = job->next;
						if (nextjob) {
								if (nextjob->job->id==jobid) {
										found = 1;
										/* skip the next item in the list */
										job->next = nextjob->next;
										/* delete the next item in the list */
										FREE(nextjob->job);
										FREE(nextjob->host);
										FREE(nextjob->arg);
										FREE(nextjob->opt);
										FREE(nextjob);
										break;
								}
						}
						job = job->next;
				}

				pthread_mutex_unlock(&mutexjoblist);

				if (!found) {
						mexWarnMsgTxt("failed to locate specified job\n");
				}

				return;
		}

		/****************************************************************************/
		else if (strcmp(command, "peerinfo")==0) {
				pthread_mutex_lock(&mutexhost);
				plhs[0] = mxCreateStructMatrix(1, 1, PEERINFO_FIELDNUMBER, peerinfo_fieldnames);
				j = 0;
				mxSetFieldByNumber(plhs[0], 0, j++, mxCreateDoubleScalar((UINT32_T)(host->id)));
				mxSetFieldByNumber(plhs[0], 0, j++, mxCreateString(host->name));
				mxSetFieldByNumber(plhs[0], 0, j++, mxCreateString(host->user));
				mxSetFieldByNumber(plhs[0], 0, j++, mxCreateString(host->group));
				mxSetFieldByNumber(plhs[0], 0, j++, mxCreateString(host->socket));
				mxSetFieldByNumber(plhs[0], 0, j++, mxCreateDoubleScalar((UINT32_T)(host->port)));
				mxSetFieldByNumber(plhs[0], 0, j++, mxCreateDoubleScalar((UINT32_T)(host->status)));
				mxSetFieldByNumber(plhs[0], 0, j++, mxCreateDoubleScalar((UINT64_T)(host->timavail)));
				mxSetFieldByNumber(plhs[0], 0, j++, mxCreateDoubleScalar((UINT64_T)(host->memavail)));
				mxSetFieldByNumber(plhs[0], 0, j++, mxCreateDoubleScalar((UINT64_T)(host->cpuavail)));
				pthread_mutex_unlock(&mutexhost);

				/* create a cell-array for allowuser */
				i = 0;
				pthread_mutex_lock(&mutexallowuserlist);
				allowuser = allowuserlist;
				while (allowuser) {
						/* count the number of items */
						i++;
						allowuser = allowuser->next;
				}
				if (i==0)
						val= mxCreateCellMatrix(0, 0);
				else 
						val= mxCreateCellMatrix(1, i);
				/* loop over the list to assign the items to a cell-array */
				allowuser = allowuserlist;
				while (allowuser) {
						/* start assigning at the back, to keep them in the original order */
						mxSetCell(val, --i, mxCreateString(allowuser->name));
						allowuser = allowuser->next;
				}
				mxSetFieldByNumber(plhs[0], 0, j++, val);
				pthread_mutex_unlock(&mutexallowuserlist);

				/* create a cell-array for allowgroup */
				i = 0;
				pthread_mutex_lock(&mutexallowgrouplist);
				allowgroup = allowgrouplist;
				while (allowgroup) {
						/* count the number of items */
						i++;
						allowgroup = allowgroup->next;
				}
				if (i==0)
						val= mxCreateCellMatrix(0, 0);
				else 
						val= mxCreateCellMatrix(1, i);
				/* loop over the list to assign the items to a cell-array */
				allowgroup = allowgrouplist;
				while (allowgroup) {
						/* start assigning at the back, to keep them in the original order */
						mxSetCell(val, --i, mxCreateString(allowgroup->name));
						allowgroup = allowgroup->next;
				}
				mxSetFieldByNumber(plhs[0], 0, j++, val);
				pthread_mutex_unlock(&mutexallowgrouplist);

				/* create a cell-array for allowhost */
				i = 0;
				pthread_mutex_lock(&mutexallowhostlist);
				allowhost = allowhostlist;
				while (allowhost) {
						/* count the number of items */
						i++;
						allowhost = allowhost->next;
				}
				if (i==0)
						val= mxCreateCellMatrix(0, 0);
				else 
						val= mxCreateCellMatrix(1, i);
				/* loop over the list to assign the items to a cell-array */
				allowhost = allowhostlist;
				while (allowhost) {
						/* start assigning at the back, to keep them in the original order */
						mxSetCell(val, --i, mxCreateString(allowhost->name));
						allowhost = allowhost->next;
				}
				mxSetFieldByNumber(plhs[0], 0, j++, val);
				pthread_mutex_unlock(&mutexallowhostlist);

				/* create a cell-array for refuseuser */
				i = 0;
				pthread_mutex_lock(&mutexrefuseuserlist);
				refuseuser = refuseuserlist;
				while (refuseuser) {
						/* count the number of items */
						i++;
						refuseuser = refuseuser->next;
				}
				if (i==0)
						val= mxCreateCellMatrix(0, 0);
				else 
						val= mxCreateCellMatrix(1, i);
				/* loop over the list to assign the items to a cell-array */
				refuseuser = refuseuserlist;
				while (refuseuser) {
						/* start assigning at the back, to keep them in the original order */
						mxSetCell(val, --i, mxCreateString(refuseuser->name));
						refuseuser = refuseuser->next;
				}
				mxSetFieldByNumber(plhs[0], 0, j++, val);
				pthread_mutex_unlock(&mutexrefuseuserlist);

				/* create a cell-array for refusegroup */
				i = 0;
				pthread_mutex_lock(&mutexrefusegrouplist);
				refusegroup = refusegrouplist;
				while (refusegroup) {
						/* count the number of items */
						i++;
						refusegroup = refusegroup->next;
				}
				if (i==0)
						val= mxCreateCellMatrix(0, 0);
				else 
						val= mxCreateCellMatrix(1, i);
				/* loop over the list to assign the items to a cell-array */
				refusegroup = refusegrouplist;
				while (refusegroup) {
						/* start assigning at the back, to keep them in the original order */
						mxSetCell(val, --i, mxCreateString(refusegroup->name));
						refusegroup = refusegroup->next;
				}
				mxSetFieldByNumber(plhs[0], 0, j++, val);
				pthread_mutex_unlock(&mutexrefusegrouplist);

				/* create a cell-array for refusehost */
				i = 0;
				pthread_mutex_lock(&mutexrefusehostlist);
				refusehost = refusehostlist;
				while (refusehost) {
						/* count the number of items */
						i++;
						refusehost = refusehost->next;
				}
				if (i==0)
						val= mxCreateCellMatrix(0, 0);
				else 
						val= mxCreateCellMatrix(1, i);
				/* loop over the list to assign the items to a cell-array */
				refusehost = refusehostlist;
				while (refusehost) {
						/* start assigning at the back, to keep them in the original order */
						mxSetCell(val, --i, mxCreateString(refusehost->name));
						refusehost = refusehost->next;
				}
				mxSetFieldByNumber(plhs[0], 0, j++, val);
				pthread_mutex_unlock(&mutexrefusehostlist);

				return;
		}

		/****************************************************************************/
		else {
				mexErrMsgTxt ("unknown command for peer");
				return;
		}

		return;
}

