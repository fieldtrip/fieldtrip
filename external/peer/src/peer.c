#include "mex.h"
#include "matrix.h"

#include "peer.h"
#include "extern.h"
#include "unix_includes.h"

#define NUMPEERSTRUCTFIELDS 5
const char* peerstructfieldnames[NUMPEERSTRUCTFIELDS] = {"hostid", "hostname", "hostaddr", "hostport", "hoststatus"};

#define NUMJOBSTRUCTFIELDS 5
const char* jobstructfieldnames[NUMJOBSTRUCTFIELDS] = {"version", "jobid", "hostsize", "argsize", "optsize"};

#define NUMJOBPEERSTRUCTFIELDS 10
const char* jobpeerstructfieldnames[NUMJOBPEERSTRUCTFIELDS] = {"version", "jobid", "hostsize", "argsize", "optsize", "hostid", "hostname", "hostaddr", "hostport", "hoststatus"}; 

int peerInitialized = 0;

/* the thread IDs are needed for cancelation at cleanup */
pthread_t tcpserverThread = 0;
pthread_t announceThread = 0;
pthread_t discoverThread = 0;
pthread_t expireThread = 0;

/* this is called the first time that the mex-file is loaded */
void initFun(void) {
		/* check whether the host has been initialized */
		if (!peerInitialized) {
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
		return;
}

void mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]) {
		char *command = NULL, *argument = NULL;
		int i, rc, t, found, success, jobid, peerid, count, server;
		jobdef_t    *def;
		joblist_t   *job, *nextjob;
		peerlist_t  *peer;
		userlist_t  *allowuser;
		grouplist_t *allowgroup;
		hostlist_t  *allowhost;
		mxArray     *arg, *opt;

		initFun();
		mexAtExit(exitFun);

		/* the first argument is always the command string */
		if (nrhs<1)
				mexErrMsgTxt ("invalid number of input arguments");
		if (!mxIsChar(prhs[0]))
				mexErrMsgTxt ("invalid input argument #1");
		/* FIXME should command be mxFree'd ? */
		command = mxArrayToString(prhs[0]);

		/****************************************************************************/
		if (strcasecmp(command, "tcpserver")==0) {
				/* the input arguments should be "tcpserver <start|stop>" */
				if (nrhs<2)
						mexErrMsgTxt ("invalid number of input arguments");
				if (!mxIsChar(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				argument = mxArrayToString(prhs[1]);
				if (strcasecmp(argument, "start")==0) {
						if (tcpserverStatus) {
								mexWarnMsgTxt("thread is already running");
								return;
						}
						mexPrintf("peer: spawning tcpserver thread\n");
						rc = pthread_create(&tcpserverThread, NULL, tcpserver, (void *)NULL);
						if (rc)
								mexErrMsgTxt("problem with return code from pthread_create()");
				}
				else if (strcasecmp(argument, "stop")==0) {
						if (!tcpserverStatus) {
								mexWarnMsgTxt("thread is not running");
								return;
						}
						mexPrintf("peer: requesting cancelation of tcpserver thread\n");
						rc = pthread_cancel(tcpserverThread);
						if (rc)
								mexErrMsgTxt("problem with return code from pthread_cancel()");
						else
								tcpserverThread = 0;
				}
				else
						mexErrMsgTxt ("invalid input argument #2");
				return;
		}

		/****************************************************************************/
		else if (strcasecmp(command, "announce")==0) {
				/* the input arguments should be "tcpserver <start|stop>" */
				if (nrhs<2)
						mexErrMsgTxt ("invalid number of input arguments");
				if (!mxIsChar(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				argument = mxArrayToString(prhs[1]);
				if (strcasecmp(argument, "start")==0) {
						if (announceStatus) {
								mexWarnMsgTxt("thread is already running");
								return;
						}
						mexPrintf("peer: spawning announce thread\n");
						rc = pthread_create(&announceThread, NULL, announce, (void *)NULL);
						if (rc)
								mexErrMsgTxt("problem with return code from pthread_create()");
				}
				else if (strcasecmp(argument, "stop")==0) {
						if (!announceStatus) {
								mexWarnMsgTxt("thread is not running");
								return;
						}
						mexPrintf("peer: requesting cancelation of announce thread\n");
						rc = pthread_cancel(announceThread);
						if (rc)
								mexErrMsgTxt("problem with return code from pthread_cancel()");
						else
								announceThread = 0;
				}
				else
						mexErrMsgTxt ("invalid input argument #2");
				return;
		}

		/****************************************************************************/
		else if (strcasecmp(command, "discover")==0) {
				/* the input arguments should be "discover <start|stop>" */
				if (nrhs<2)
						mexErrMsgTxt ("invalid number of input arguments");
				if (!mxIsChar(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				argument = mxArrayToString(prhs[1]);
				if (strcasecmp(argument, "start")==0) {
						if (discoverStatus) {
								mexWarnMsgTxt("thread is already running");
								return;
						}
						mexPrintf("peer: spawning discover thread\n");
						rc = pthread_create(&discoverThread, NULL, discover, (void *)NULL);
						if (rc)
								mexErrMsgTxt("problem with return code from pthread_create()");
				}
				else if (strcasecmp(argument, "stop")==0) {
						if (!discoverStatus) {
								mexWarnMsgTxt("thread is not running");
								return;
						}
						mexPrintf("peer: requesting cancelation of discover thread\n");
						rc = pthread_cancel(discoverThread);
						if (rc)
								mexErrMsgTxt("problem with return code from pthread_cancel()");
						else
								discoverThread = 0;
				}
				else
						mexErrMsgTxt ("invalid input argument #2");
				return;
		}

		/****************************************************************************/
		else if (strcasecmp(command, "expire")==0) {
				/* the input arguments should be "expire <start|stop>" */
				if (nrhs<2)
						mexErrMsgTxt ("invalid number of input arguments");
				if (!mxIsChar(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				argument = mxArrayToString(prhs[1]);
				if (strcasecmp(argument, "start")==0) {
						if (expireStatus) {
								mexWarnMsgTxt("thread is already running");
								return;
						}
						mexPrintf("peer: spawning expire thread\n");
						rc = pthread_create(&expireThread, NULL, expire, (void *)NULL);
						if (rc)
								mexErrMsgTxt("problem with return code from pthread_create()");
				}
				else if (strcasecmp(argument, "stop")==0) {
						if (!expireStatus) {
								mexWarnMsgTxt("thread is not running");
								return;
						}
						mexPrintf("peer: requesting cancelation of expire thread\n");
						rc = pthread_cancel(expireThread);
						if (rc)
								mexErrMsgTxt("problem with return code from pthread_cancel()");
						else
								expireThread = 0;
				}
				else
						mexErrMsgTxt ("invalid input argument #2");
				return;
		}

		/****************************************************************************/
		else if (strcasecmp(command, "info")==0) {
				/* the input arguments should be "info" */

				/* display  all possible information on screen */
				pthread_mutex_lock(&mutexstatus);
				mexPrintf("tcpserverStatus = %d\n", tcpserverStatus);
				mexPrintf("announceStatus  = %d\n", announceStatus);
				mexPrintf("discoverStatus  = %d\n", discoverStatus);
				mexPrintf("expireStatus    = %d\n", expireStatus);
				pthread_mutex_unlock(&mutexstatus);

				pthread_mutex_lock(&mutexhost);
				mexPrintf("host.name   = %s\n", host->name);
				mexPrintf("host.addr   = %s\n", host->addr);
				mexPrintf("host.port   = %d\n", host->port);
				mexPrintf("host.user   = %s\n", host->user);
				mexPrintf("host.group  = %s\n", host->group);
				mexPrintf("host.status = %d\n", host->status);
				mexPrintf("host.id     = %d\n", host->id);
				pthread_mutex_unlock(&mutexhost);

				pthread_mutex_lock(&mutexuserlist);
				allowuser = userlist;
				while (allowuser) {
						mexPrintf("allowuser = %s\n", allowuser->name);
						allowuser = allowuser->next;
				}
				pthread_mutex_unlock(&mutexuserlist);

				pthread_mutex_lock(&mutexgrouplist);
				allowgroup = grouplist;
				while (allowgroup) {
						mexPrintf("allowgroup = %s\n", allowgroup->name);
						allowgroup = allowgroup->next;
				}
				pthread_mutex_unlock(&mutexgrouplist);

				pthread_mutex_lock(&mutexpeerlist);
				i = 0;
				peer = peerlist;
				while(peer) {
						mexPrintf("peerlist[%d] = \n", i);
						mexPrintf("  host.name   = %s\n", peer->host->name);
						mexPrintf("  host.addr   = %s\n", peer->host->addr);
						mexPrintf("  host.port   = %d\n", peer->host->port);
						mexPrintf("  host.user   = %s\n", peer->host->user);
						mexPrintf("  host.group  = %s\n", peer->host->group);
						mexPrintf("  host.status = %d\n", peer->host->status);
						mexPrintf("  host.id     = %d\n", peer->host->id);
						mexPrintf("  time        = %s",   ctime(&(peer->time)));
						peer = peer->next ;       
						i++;
				}
				pthread_mutex_unlock(&mutexpeerlist);

				pthread_mutex_lock(&mutexjoblist);
				i = 0;
				job = joblist;
				while(job) {
						mexPrintf("joblist[%d] = \n", i);
						mexPrintf("  job.version = %d\n", job->job->version);
						mexPrintf("  job.id      = %d\n", job->job->id);
						mexPrintf("  job.argsize = %d\n", job->job->argsize);
						mexPrintf("  job.optsize = %d\n", job->job->optsize);
						mexPrintf("  host.name   = %s\n", job->host->name);
						mexPrintf("  host.id     = %d\n", job->host->id);
						job = job->next ;
						i++;
				}
				pthread_mutex_unlock(&mutexjoblist);
				return;
		}

		/****************************************************************************/
		else if (strcasecmp(command, "status")==0) {
				/* the input arguments should be "status <number>" */
				if (nrhs<2)
						mexErrMsgTxt ("invalid number of input arguments");
				if (!mxIsNumeric(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				if (!mxIsScalar(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				pthread_mutex_lock(&mutexhost);
				host->status = (UINT32_T)mxGetScalar(prhs[1]);
				pthread_mutex_unlock(&mutexhost);
		}

		/****************************************************************************/
		else if (strcasecmp(command, "tcpport")==0) {
				/* the input arguments should be "tcpport <number>" */
				if (tcpserverThread)
						mexErrMsgTxt ("cannot change the port while the tcpserver is running");
				if (nrhs<2)
						mexErrMsgTxt ("invalid number of input arguments");
				if (!mxIsNumeric(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				if (!mxIsScalar(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				pthread_mutex_lock(&mutexhost);
				host->port = (UINT32_T)mxGetScalar(prhs[1]);
				pthread_mutex_unlock(&mutexhost);
		}

		/****************************************************************************/
		else if (strcasecmp(command, "group")==0) {
				/* the input arguments should be "group <string>" */
				if (nrhs<2)
						mexErrMsgTxt ("invalid number of input arguments");
				if (!mxIsChar(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				pthread_mutex_lock(&mutexhost);
				if (mxGetString(prhs[1], host->group, STRLEN))
						mexErrMsgTxt("FIXME: unexpected error");
				pthread_mutex_unlock(&mutexhost);
		}

		/****************************************************************************/
		else if (strcasecmp(command, "allowuser")==0) {
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
				pthread_mutex_lock(&mutexuserlist);
				/* erase the existing list */
				allowuser = userlist;
				while (allowuser) {
						userlist = allowuser->next;
						FREE(allowuser->name);
						FREE(allowuser);
						allowuser = userlist;
				}
				/* add all elements to the list */
				for (i=0; i<mxGetNumberOfElements(prhs[1]); i++) {
						arg = mxGetCell(prhs[1], i);
						allowuser = (userlist_t *)malloc(sizeof(userlist_t));
						allowuser->name = malloc(mxGetNumberOfElements(arg)+1);
						allowuser->next = userlist;
						if(mxGetString(arg, allowuser->name, mxGetNumberOfElements(arg)+1)!=0)
								mexErrMsgTxt("FIXME: unexpected error, memory leak");
						userlist = allowuser;
				}
				pthread_mutex_unlock(&mutexuserlist);
		}

		/****************************************************************************/
		else if (strcasecmp(command, "allowgroup")==0) {
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
				pthread_mutex_lock(&mutexgrouplist);
				/* erase the existing list */
				allowgroup = grouplist;
				while (allowgroup) {
						grouplist = allowgroup->next;
						FREE(allowgroup->name);
						FREE(allowgroup);
						allowgroup = grouplist;
				}
				/* add all elements to the list */
				for (i=0; i<mxGetNumberOfElements(prhs[1]); i++) {
						arg = mxGetCell(prhs[1], i);
						allowgroup = (grouplist_t *)malloc(sizeof(grouplist_t));
						allowgroup->name = malloc(mxGetNumberOfElements(arg)+1);
						allowgroup->next = grouplist;
						if(mxGetString(arg, allowgroup->name, mxGetNumberOfElements(arg)+1)!=0)
								mexErrMsgTxt("FIXME: unexpected error, memory leak");
						grouplist = allowgroup;
				}
				pthread_mutex_unlock(&mutexgrouplist);
		}

		/****************************************************************************/
		else if (strcasecmp(command, "allowhost")==0) {
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
				pthread_mutex_lock(&mutexhostlist);
				/* erase the existing list */
				allowhost = hostlist;
				while (allowhost) {
						hostlist = allowhost->next;
						FREE(allowhost->name);
						FREE(allowhost);
						allowhost = hostlist;
				}
				/* add all elements to the list */
				for (i=0; i<mxGetNumberOfElements(prhs[1]); i++) {
						arg = mxGetCell(prhs[1], i);
						allowhost = (hostlist_t *)malloc(sizeof(hostlist_t));
						allowhost->name = malloc(mxGetNumberOfElements(arg)+1);
						allowhost->next = hostlist;
						if(mxGetString(arg, allowhost->name, mxGetNumberOfElements(arg)+1)!=0)
								mexErrMsgTxt("FIXME: unexpected error, memory leak");
						hostlist = allowhost;
				}
				pthread_mutex_unlock(&mutexhostlist);
		}

		/****************************************************************************/
		else if (strcasecmp(command, "put")==0) {
				/* the input arguments should be "put <peerid> <arg> <opt>"         */
				/* the input arguments should be "put <peerid> <arg> <opt> <jobid>" */

				if (nrhs<2)
						mexErrMsgTxt("invalid argument #2");
				if (!mxIsNumeric(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				peerid = (UINT32_T)mxGetScalar(prhs[1]);

				if (nrhs>4) {
						/* use the specified jobid */
						if (!mxIsNumeric(prhs[4]))
								mexErrMsgTxt ("invalid input argument #5");
						jobid = (UINT32_T)mxGetScalar(prhs[4]);
				}
				else {
						/* assign a random jobid */
						jobid = random();
				}

				pthread_mutex_lock(&mutexpeerlist);
				found = 0;

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

				/* open the TCP socket */
				if ((server = open_connection(peer->host->name, peer->host->port)) < 0) {
						pthread_mutex_unlock(&mutexpeerlist);
						mexErrMsgTxt("failed to create socket\n");
				}

				pthread_mutex_unlock(&mutexpeerlist);

				if (nrhs<3)
						mexErrMsgTxt("invalid argument #3");

				if (nrhs<4)
						mexErrMsgTxt("invalid argument #4");

				def = (jobdef_t *)malloc(sizeof(jobdef_t));
				def->version  = VERSION;
				def->id       = jobid;
				def->hostsize = sizeof(hostdef_t);

				arg = (mxArray *) mxSerialize(prhs[2]);
				opt = (mxArray *) mxSerialize(prhs[3]);

				def->argsize  = mxGetNumberOfElements(arg);
				def->optsize  = mxGetNumberOfElements(opt);

				/* write the message, the slave may close the connection if he does not accept the job */
				success = 1;
				if (success)
						success = (bufwrite(server, def, sizeof(jobdef_t)) == sizeof(jobdef_t));
				if (success) 
						success = (bufwrite(server, host, sizeof(hostdef_t)) == sizeof(hostdef_t));
				if (success) 
						success = (bufwrite(server, (void *)mxGetData(arg), def->argsize) == def->argsize);
				if (success) 
						success = (bufwrite(server, (void *)mxGetData(opt), def->optsize) == def->optsize);

				close_connection(server);

				mxDestroyArray(arg);
				mxDestroyArray(opt);

				if (success) {
						/* return the job details */
						plhs[0] = mxCreateStructMatrix(1, 1, NUMJOBSTRUCTFIELDS, jobstructfieldnames);
						mxSetFieldByNumber(plhs[0], 0, 0, mxCreateDoubleScalar((UINT32_T)(def->version)));
						mxSetFieldByNumber(plhs[0], 0, 1, mxCreateDoubleScalar((UINT32_T)(def->id)));
						mxSetFieldByNumber(plhs[0], 0, 2, mxCreateDoubleScalar((UINT32_T)(def->hostsize)));
						mxSetFieldByNumber(plhs[0], 0, 3, mxCreateDoubleScalar((UINT32_T)(def->argsize)));
						mxSetFieldByNumber(plhs[0], 0, 4, mxCreateDoubleScalar((UINT32_T)(def->optsize)));
						FREE(def);
				}
				else {
						FREE(def);
						mexErrMsgTxt ("failed to put job arguments");
				}

				return;
		}

		/****************************************************************************/
		else if (strcasecmp(command, "get")==0) {
				/* the input arguments should be "get <jobid>" */
				if (nrhs<2)
						mexErrMsgTxt ("invalid number of input arguments");
				if (!mxIsNumeric(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				if (!mxIsScalar(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				jobid = (UINT32_T)mxGetScalar(prhs[1]);

				pthread_mutex_lock(&mutexjoblist);
				found = 0;
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
		else if (strcasecmp(command, "clear")==0) {
				/* the input arguments should be "clear <jobid>" */
				if (nrhs<2)
						mexErrMsgTxt ("invalid number of input arguments");
				if (!mxIsNumeric(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				if (!mxIsScalar(prhs[1]))
						mexErrMsgTxt ("invalid input argument #2");
				jobid = (UINT32_T)mxGetScalar(prhs[1]);

				pthread_mutex_lock(&mutexjoblist);
				found = 0;

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

				if (!found) {
						mexWarnMsgTxt("failed to locate specified job\n");
				}

				pthread_mutex_unlock(&mutexjoblist);
				return;
		}

		/****************************************************************************/
		else if (strcasecmp(command, "peerlist")==0) {
				pthread_mutex_lock(&mutexpeerlist);
				/* count the number of peers */
				i = 0;
				peer = peerlist;
				while(peer) {
						i++;
						peer = peer->next ;
				}
				plhs[0] = mxCreateStructMatrix(i, 1, NUMPEERSTRUCTFIELDS, peerstructfieldnames);
				i = 0;
				peer = peerlist;
				while(peer) {
						mxSetFieldByNumber(plhs[0], i, 0, mxCreateDoubleScalar((UINT32_T)(peer->host->id)));
						mxSetFieldByNumber(plhs[0], i, 1, mxCreateString(peer->host->name));
						mxSetFieldByNumber(plhs[0], i, 2, mxCreateString(peer->host->addr));
						mxSetFieldByNumber(plhs[0], i, 3, mxCreateDoubleScalar((UINT32_T)(peer->host->port)));
						mxSetFieldByNumber(plhs[0], i, 4, mxCreateDoubleScalar((UINT32_T)(peer->host->status)));
						i++;
						peer = peer->next ;
				}
				pthread_mutex_unlock(&mutexpeerlist);
				return;
		}

		/****************************************************************************/
		else if (strcasecmp(command, "joblist")==0) {
				pthread_mutex_lock(&mutexjoblist);
				/* count the number of jobs */
				i = 0;
				job = joblist;
				while(job) {
						i++;
						job = job->next ;
				}
				plhs[0] = mxCreateStructMatrix(i, 1, NUMJOBPEERSTRUCTFIELDS, jobpeerstructfieldnames);
				i = 0;
				job = joblist;
				while(job) {
						mxSetFieldByNumber(plhs[0], i, 0, mxCreateDoubleScalar((UINT32_T)(job->job->version)));
						mxSetFieldByNumber(plhs[0], i, 1, mxCreateDoubleScalar((UINT32_T)(job->job->id)));
						mxSetFieldByNumber(plhs[0], i, 2, mxCreateDoubleScalar((UINT32_T)(job->job->hostsize)));
						mxSetFieldByNumber(plhs[0], i, 3, mxCreateDoubleScalar((UINT32_T)(job->job->argsize)));
						mxSetFieldByNumber(plhs[0], i, 4, mxCreateDoubleScalar((UINT32_T)(job->job->optsize)));
						mxSetFieldByNumber(plhs[0], i, 5, mxCreateDoubleScalar((UINT32_T)(job->host->id)));
						mxSetFieldByNumber(plhs[0], i, 6, mxCreateString(job->host->name));
						mxSetFieldByNumber(plhs[0], i, 7, mxCreateString(job->host->addr));
						mxSetFieldByNumber(plhs[0], i, 8, mxCreateDoubleScalar((UINT32_T)(job->host->port)));
						mxSetFieldByNumber(plhs[0], i, 9, mxCreateDoubleScalar((UINT32_T)(job->host->status)));
						job = job->next ;
						i++;
				}
				pthread_mutex_unlock(&mutexjoblist);
				return;
		}

		/****************************************************************************/
		else {
				mexErrMsgTxt ("unknown command for peer");
				return;
		}

		return;
}

