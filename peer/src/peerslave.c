#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "mex.h"
#include "engine.h"

#include "peer.h"
#include "extern.h"
#include "platform_includes.h"

mxArray *mxSerialize(const mxArray*);
mxArray *mxDeserialize(const void*, size_t);

#define panic(X) {fprintf(stderr, X); exit(1);}
#define MAXIDLETIME 30     /* in seconds */
#define SLEEPTIME   10000  /* in microseconds */

#ifndef STARTCMD
#define STARTCMD "/Applications/MATLAB_R2009a.app/bin/matlab"
#endif

int main(int argc, char *argv[]) {
		Engine *en;
		mxArray *argin, *argout, *options, *arg, *opt;
		joblist_t  *job  = NULL;
		peerlist_t *peer = NULL;
		jobdef_t   *def  = NULL;

		int matlabRunning = 0, matlabStart, matlabFinished;
		int n, rc, found, handshake, success, server, verbose = 0, jobnum = 0;
		unsigned int peerid, jobid;

		/* the thread IDs are needed for cancelation at cleanup */
		pthread_t tcpserverThread;
		pthread_t announceThread;
		pthread_t discoverThread;
		pthread_t expireThread;

		peerinit(NULL);

		if ((rc = pthread_create(&tcpserverThread, NULL, tcpserver, (void *)NULL))>0) {
				fprintf(stderr, "failed to start tcpserver thread\n");
				exit(1);
		}
		else {
				fprintf(stderr, "started tcpserver thread\n");
		}

		if ((rc = pthread_create(&announceThread, NULL, announce, (void *)NULL))>0) {
				fprintf(stderr, "failed to start announce thread\n");
				exit(1);
		}
		else {
				fprintf(stderr, "started announce thread\n");
		}

		if ((rc = pthread_create(&discoverThread, NULL, discover, (void *)NULL))>0) {
				fprintf(stderr, "failed to start discover thread\n");
				exit(1);
		}
		else {
				fprintf(stderr, "started discover thread\n");
		}

		if ((rc = pthread_create(&expireThread, NULL, expire, (void *)NULL))>0) {
				fprintf(stderr, "failed to start expire thread\n");
				exit(1);
		}
		else {
				fprintf(stderr, "started expire thread\n");
		}

		/* status = 0 means zombie mode, don't accept anything   */
		/* status = 1 means master mode, accept everything       */
		/* status = 2 means idle slave, accept only a single job */
		/* status = 3 means busy slave, don't accept a new job   */
		/* any other status is interpreted as zombie mode        */
		host->status = 2;

		while (1) {

				if (jobcount()>0) {

						/* there is a job to be executed */
						if (matlabRunning==0) {
								if ((en = engOpen(STARTCMD)) == NULL) {
										panic("failed to start MATLAB engine\n");
								}
								else {
										fprintf(stderr, "started MATLAB engine\n");
										matlabRunning = 1;
								}
						}

						/* switch the mode to busy slave */
						host->status = 3;
						matlabStart = time(NULL);

						/* get the first job input arguments and options */
						pthread_mutex_lock(&mutexjoblist);
						job = joblist;
						argin   = (mxArray *)mxDeserialize(job->arg, job->job->argsize);
						options = (mxArray *)mxDeserialize(job->opt, job->job->optsize);
						jobid   = job->job->id;
						peerid  = job->host->id;
						fprintf(stderr, "executing job %d from %s@%s (jobid=%d)\n", ++jobnum, job->host->user, job->host->name, jobid);
						pthread_mutex_unlock(&mutexjoblist);

						/* copy them over to the engine */
						engPutVariable(en, "argin", argin);
						engPutVariable(en, "options", options);
						mxDestroyArray(argin);
						mxDestroyArray(options);

						/* execute the job */
						engEvalString(en, "[argout, options] = peerexec(argin, options);");

						/* get the job output arguments and options */
						if ((argout = engGetVariable(en, "argout")) == NULL) {
								fprintf(stderr, "error getting argout");
								exit(1);
						}

						if ((options = engGetVariable(en, "options")) == NULL) {
								fprintf(stderr, "error getting options");
								exit(1);
						}

						/* send them back to the master */

						/*****************************************************************************
						 * the following code is largely shared with the get-option in the peer mex file
						 *****************************************************************************/
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
								panic("failed to locate specified peer\n");
						}

						/* open the TCP socket */
						if ((server = open_connection(peer->ipaddr, peer->host->port)) < 0) {
								pthread_mutex_unlock(&mutexpeerlist);
								panic("failed to create socket\n");
						}

						pthread_mutex_unlock(&mutexpeerlist);

						if ((n = bufread(server, &handshake, sizeof(int))) != sizeof(int)) {
								panic("tcpsocket: could not write handshake");
						}
						else if (!handshake) {
								close_connection(server);
								panic("failed to negociate connection");
						}

						/* the message that will be written consists of
						   message->host
						   message->job
						   message->arg
						   message->opt
						 */

						arg = (mxArray *) mxSerialize(argout);
						opt = (mxArray *) mxSerialize(options);

						if (!arg)
								panic("could not serialize job arguments");

						if (!opt)
								panic("could not serialize job options");

						def = (jobdef_t *)malloc(sizeof(jobdef_t));

						if (!def)
								panic("could not allocate memory");

						def->version  = VERSION;
						def->id       = jobid;
						def->memreq   = 0;
						def->cpureq   = 0;
						def->timreq   = 0;
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
								panic("could not write handshake");
						}
						else if (!handshake) {
								close_connection(server);
								panic("failed to write hostdef");
						}

						if (success)
								success = (bufwrite(server, def, sizeof(jobdef_t)) == sizeof(jobdef_t));

						if ((n = bufread(server, &handshake, sizeof(int))) != sizeof(int)) {
								panic("could not write handshake");
						}
						else if (!handshake) {
								close_connection(server);
								panic("failed to write jobdef");
						}

						if (success) 
								success = (bufwrite(server, (void *)mxGetData(arg), def->argsize) == def->argsize);

						if ((n = bufread(server, &handshake, sizeof(int))) != sizeof(int)) {
								panic("could not write handshake");
						}
						else if (!handshake) {
								close_connection(server);
								panic("failed to write arg");
						}

						if (success) 
								success = (bufwrite(server, (void *)mxGetData(opt), def->optsize) == def->optsize);

						if ((n = bufread(server, &handshake, sizeof(int))) != sizeof(int)) {
								panic("could not write handshake");
						}
						else if (!handshake) {
								close_connection(server);
								panic("failed to write opt");
						}

						close_connection(server);

						mxDestroyArray(arg);
						mxDestroyArray(opt);

						/*****************************************************************************/

						/* remove the first job from the joblist */
						pthread_mutex_lock(&mutexjoblist);
						joblist = job->next;
						FREE(job->job);
						FREE(job->host);
						FREE(job->arg);
						FREE(job->opt);
						FREE(job);
						job = joblist;
						pthread_mutex_unlock(&mutexjoblist);

						host->status = 2; /* idle slave */
						matlabFinished = time(NULL);
						fprintf(stderr, "executing job %d took %d seconds\n", jobnum, matlabFinished - matlabStart);
				}
				else {
						usleep(SLEEPTIME);
				} /* if jobcount */

				/* test that the matlab engine is not idle for too long */
				if ((matlabRunning==1) && (time(NULL)-matlabFinished)>MAXIDLETIME) {
						if (engClose(en)!=0) {
								panic("could not stop the MATLAB engine\n");
						}
						else {
								fprintf(stderr, "stopped idle MATLAB engine\n");
								matlabRunning = 0;
						}
				}

		} /* while */

		return 0;
} /* main */

