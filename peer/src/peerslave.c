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

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>

#include "engine.h"
#include "matrix.h"

#include "peer.h"
#include "extern.h"
#include "platform_includes.h"

mxArray *mxSerialize(const mxArray*);
mxArray *mxDeserialize(const void*, size_t);

#define ENGINETIMEOUT     30     /* in seconds */
#define ZOMBIETIMEOUT     300    /* in seconds */
#define SLEEPTIME         10000  /* in microseconds */

#define STARTCMD "matlab -nosplash"

void print_help(char *argv[]) {
		printf("\n");
		printf("This starts a FieldTrip peer-to-peer distributed computing peer, which\n");
		printf("will wait for an incoming job and subsequently start the MATLAB engine and\n");
		printf("evaluate the job. Use as\n");
		printf("  %s [options]\n", argv[0]);
		printf("where the options can include\n");
		printf("  --memavail    = number, amount of memory available       (default = inf)\n");
		printf("  --cpuavail    = number, speed of the CPU                 (default = inf)\n");
		printf("  --timavail    = number, maximum duration of a single job (default = inf)\n");
		printf("  --allowhost   = {...}\n");
		printf("  --allowuser   = {...}\n");
		printf("  --allowgroup  = {...}\n");
		printf("  --group       = string\n");
		printf("  --hostname    = string\n");
		printf("  --matlab      = string\n");
		printf("  --timeout     = number, time to keep the engine running after the job finished\n");
		printf("  --smartshare  = 0|1\n");
		printf("  --smartmem    = 0|1\n");
		printf("  --smartcpu    = 0|1\n");
		printf("  --daemon\n");
		printf("  --udsserver\n");
		printf("  --verbose\n");
		printf("  --help\n");
		printf("\n");
}

int help_flag;
int daemon_flag;
int tcpserver_flag = 1;
int udpserver_flag = 0;
int udsserver_flag = 0;

int main(int argc, char *argv[]) {
		Engine *en;
		mxArray *argin, *argout, *options, *arg, *opt;
		joblist_t  *job  = NULL;
		peerlist_t *peer = NULL;
		jobdef_t   *def  = NULL;
		pid_t childpid;

		int matlabRunning = 0, matlabStart, matlabFinished, engineFailed = 0;
		int c, rc, found, handshake, success, server, jobnum = 0, jobcount, jobfailed, timallow;
		unsigned int enginetimeout = ENGINETIMEOUT;
		unsigned int zombietimeout = ZOMBIETIMEOUT;
		unsigned int peerid, jobid;
		char *str = NULL, *startcmd = NULL;

		userlist_t  *allowuser;
		grouplist_t *allowgroup;
		hostlist_t  *allowhost;

		/* the thread IDs are needed for cancelation at cleanup */
		pthread_t udsserverThread;
		pthread_t tcpserverThread;
		pthread_t announceThread;
		pthread_t discoverThread;
		pthread_t expireThread;

		openlog("peerslave", LOG_PID | LOG_PERROR, LOG_USER);
		setlogmask(LOG_MASK(LOG_EMERG) | LOG_MASK(LOG_ALERT) | LOG_MASK(LOG_CRIT));
		peerinit(NULL);

		/* use GNU getopt_long for the command-line options */
		while (1)
		{
				static struct option long_options[] =
				{
						{"help",       no_argument, &help_flag, 1},
						{"daemon",     no_argument, &daemon_flag, 1},
						{"udsserver",  no_argument, &udsserver_flag, 1},
						{"memavail",   required_argument, 0, 'a'}, /* numeric argument */
						{"cpuavail",   required_argument, 0, 'b'}, /* numeric argument */
						{"timavail",   required_argument, 0, 'c'}, /* numeric argument */
						{"hostname",   required_argument, 0, 'd'}, /* single string argument */
						{"group",      required_argument, 0, 'e'}, /* single string argument */
						{"allowuser",  required_argument, 0, 'f'}, /* single or multiple string argument */
						{"allowhost",  required_argument, 0, 'g'}, /* single or multiple string argument */
						{"allowgroup", required_argument, 0, 'h'}, /* single or multiple string argument */
						{"matlab",     required_argument, 0, 'i'}, /* single string argument */
						{"smartmem",   required_argument, 0, 'j'}, /* numeric, 0 or 1 */
						{"smartcpu",   required_argument, 0, 'k'}, /* numeric, 0 or 1 */
						{"smartshare", required_argument, 0, 'l'}, /* numeric, 0 or 1 */
						{"timeout",    required_argument, 0, 'm'}, /* numeric argument */
						{"verbose",    required_argument, 0, 'n'}, /* numeric argument */
						{0, 0, 0, 0}
				};

				/* getopt_long stores the option index here. */
				int option_index = 0;

				c = getopt_long (argc, argv, "", long_options, &option_index); 

				/* detect the end of the options. */
				if (c == -1)
						break;

				switch (c)
				{
						case 0:
								/* if this option set a flag, do nothing else now. */
								if (long_options[option_index].flag != 0)
										break;
								break;

						case 'a':
								DEBUG(LOG_NOTICE, "option --memavail with value `%s'", optarg);
								pthread_mutex_lock(&mutexhost);
								host->memavail = atol(optarg);
								pthread_mutex_unlock(&mutexhost);
								pthread_mutex_lock(&mutexsmartmem);
								smartmem.enabled = 0;
								pthread_mutex_unlock(&mutexsmartmem);
								break;

						case 'b':
								DEBUG(LOG_NOTICE, "option --cpuavail with value `%s'", optarg);
								pthread_mutex_lock(&mutexhost);
								host->cpuavail = atol(optarg);
								pthread_mutex_unlock(&mutexhost);
								break;

						case 'c':
								DEBUG(LOG_NOTICE, "option --timavail with value `%s'", optarg);
								pthread_mutex_lock(&mutexhost);
								host->timavail = atol(optarg);
								pthread_mutex_unlock(&mutexhost);
								break;

						case 'd':
								DEBUG(LOG_NOTICE, "option --hostname with value `%s'", optarg);
								pthread_mutex_lock(&mutexhost);
								strncpy(host->name, optarg, STRLEN);
								pthread_mutex_unlock(&mutexhost);
								break;

						case 'e':
								DEBUG(LOG_NOTICE, "option --group with value `%s'", optarg);
								pthread_mutex_lock(&mutexhost);
								strncpy(host->group, optarg, STRLEN);
								pthread_mutex_unlock(&mutexhost);
								break;

						case 'f':
								DEBUG(LOG_NOTICE, "option --allowuser with value `%s'", optarg);
								str = strtok(optarg, ",");
								while (str) {
										allowuser = (userlist_t *)malloc(sizeof(userlist_t));
										allowuser->name = malloc(STRLEN);
										strncpy(allowuser->name, str, STRLEN);
										allowuser->next = userlist;
										userlist = allowuser;
										str = strtok(NULL, ",");
								}
								break;

						case 'g':
								DEBUG(LOG_NOTICE, "option --allowhost with value `%s'", optarg);
								str = strtok(optarg, ",");
								while (str) {
										allowhost = (hostlist_t *)malloc(sizeof(hostlist_t));
										allowhost->name = malloc(STRLEN);
										strncpy(allowhost->name, str, STRLEN);
										allowhost->next = hostlist;
										hostlist = allowhost;
										str = strtok(NULL, ",");
								}
								break;

						case 'h':
								DEBUG(LOG_NOTICE, "option --allowgroup with value `%s'", optarg);
								str = strtok(optarg, ",");
								while (str) {
										allowgroup = (grouplist_t *)malloc(sizeof(grouplist_t));
										allowgroup->name = malloc(STRLEN);
										strncpy(allowgroup->name, str, STRLEN);
										allowgroup->next = grouplist;
										grouplist = allowgroup;
										str = strtok(NULL, ",");
								}
								break;

						case 'i':
								DEBUG(LOG_NOTICE, "option --matlab with value `%s'", optarg);
								startcmd = malloc(STRLEN);
								strncpy(startcmd, optarg, STRLEN);
								break;

						case 'j':
								DEBUG(LOG_NOTICE, "option --smartmem with value `%s'", optarg);
								pthread_mutex_lock(&mutexsmartmem);
								smartmem.enabled = atol(optarg);
								pthread_mutex_unlock(&mutexsmartmem);
								break;

						case 'k':
								DEBUG(LOG_NOTICE, "option --smartcpu with value `%s'", optarg);
								pthread_mutex_lock(&mutexsmartcpu);
								smartcpu.enabled = atol(optarg);
								pthread_mutex_unlock(&mutexsmartcpu);
								break;

						case 'l':
								DEBUG(LOG_NOTICE, "option --smartshare with value `%s'", optarg);
								pthread_mutex_lock(&mutexsmartshare);
								smartshare.enabled = atol(optarg);
								pthread_mutex_unlock(&mutexsmartshare);
								break;

						case 'm':
								DEBUG(LOG_NOTICE, "option --timeout with value `%s'", optarg);
								enginetimeout = atol(optarg);
								break;

						case 'n':
								DEBUG(LOG_NOTICE, "option --verbose with value `%s'", optarg);
								syslog_level = atol(optarg);
								break;

						case '?':
								/* getopt_long already printed an error message. */
								break;

						default:
								PANIC("invalid command line options\n");
								break;
				}
		}

		switch (syslog_level) {
				case 0:
						setlogmask(LOG_MASK(LOG_EMERG) | LOG_MASK(LOG_ALERT) | LOG_MASK(LOG_CRIT) | LOG_MASK(LOG_ERR) | LOG_MASK(LOG_WARNING) | LOG_MASK(LOG_NOTICE) | LOG_MASK(LOG_INFO) | LOG_MASK(LOG_DEBUG));
						break;
				case 1:
						setlogmask(LOG_MASK(LOG_EMERG) | LOG_MASK(LOG_ALERT) | LOG_MASK(LOG_CRIT) | LOG_MASK(LOG_ERR) | LOG_MASK(LOG_WARNING) | LOG_MASK(LOG_NOTICE) | LOG_MASK(LOG_INFO));
						break;
				case 2:
						setlogmask(LOG_MASK(LOG_EMERG) | LOG_MASK(LOG_ALERT) | LOG_MASK(LOG_CRIT) | LOG_MASK(LOG_ERR) | LOG_MASK(LOG_WARNING) | LOG_MASK(LOG_NOTICE));
						break;
				case 3:
						setlogmask(LOG_MASK(LOG_EMERG) | LOG_MASK(LOG_ALERT) | LOG_MASK(LOG_CRIT) | LOG_MASK(LOG_ERR) | LOG_MASK(LOG_WARNING));
						break;
				case 4:
						setlogmask(LOG_MASK(LOG_EMERG) | LOG_MASK(LOG_ALERT) | LOG_MASK(LOG_CRIT) | LOG_MASK(LOG_ERR));
						break;
				case 5:
						setlogmask(LOG_MASK(LOG_EMERG) | LOG_MASK(LOG_ALERT) | LOG_MASK(LOG_CRIT));
						break;
				case 6:
						setlogmask(LOG_MASK(LOG_EMERG) | LOG_MASK(LOG_ALERT));
						break;
				case 7:
						setlogmask(LOG_MASK(LOG_EMERG));
						break;
		}

		if (help_flag) {
				/* display the help message and return to the command line */
				print_help(argv);
				exit(0);
		}

		if (daemon_flag) {
				/* now create new process */
				childpid = fork();

				if (childpid >= 0) /* fork succeeded */
				{
						if (childpid == 0) /* fork() returns 0 to the child process */
						{
								printf("CHILD: I am the child process!\n");
								printf("CHILD: Here's my PID: %d\n", getpid());
								printf("CHILD: My parent's PID is: %d\n", getppid());
								printf("CHILD: The value of my copy of childpid is: %d\n", childpid);
								/* the child continues as the actual executable */
						}
						else /* fork() returns new pid to the parent process */
						{
								printf("PARENT: I am the parent process!\n");
								printf("PARENT: Here's my PID: %d\n", getpid());
								printf("PARENT: The value of my copy of childpid is %d\n", childpid);
								/* parent exits */ 
								return 0;
						}
				}
				else /* fork returns -1 on failure */
				{
						perror("fork"); /* display error message */
						exit(0); 
				}
		}

		/* set the default command to start matlab */
		if (!startcmd) {
				startcmd = malloc(STRLEN);
				strncpy(startcmd, STARTCMD, STRLEN);
		}

		if (udsserver_flag) {
				if ((rc = pthread_create(&udsserverThread, NULL, udsserver, (void *)NULL))>0) {
						PANIC("failed to start udsserver thread\n");
				}
				else {
						DEBUG(LOG_NOTICE, "started udsserver thread");
				}
		}

		if (tcpserver_flag) {
				if ((rc = pthread_create(&tcpserverThread, NULL, tcpserver, (void *)NULL))>0) {
						PANIC("failed to start tcpserver thread\n");
				}
				else {
						DEBUG(LOG_NOTICE, "started tcpserver thread");
				}
		}

		if ((rc = pthread_create(&announceThread, NULL, announce, (void *)NULL))>0) {
				PANIC("failed to start announce thread\n");
		}
		else {
				DEBUG(LOG_NOTICE, "started announce thread");
		}

		if ((rc = pthread_create(&discoverThread, NULL, discover, (void *)NULL))>0) {
				PANIC("failed to start discover thread\n");
		}
		else {
				DEBUG(LOG_NOTICE, "started discover thread");
		}

		if ((rc = pthread_create(&expireThread, NULL, expire, (void *)NULL))>0) {
				DEBUG(LOG_NOTICE, "failed to start expire thread");
				exit(1);
		}
		else {
				DEBUG(LOG_NOTICE, "started expire thread");
		}

		/* start as idle slave */
		pthread_mutex_lock(&mutexhost);
		host->status = STATUS_IDLE;
		/* update the current job description */
		bzero(&(host->current), sizeof(current_t));
		pthread_mutex_unlock(&mutexhost);

		while (1) {

				pthread_mutex_lock(&mutexjoblist);
				jobcount = 0;
				job = joblist;
				while (job) {
						jobcount++;
						job = job->next;
				}
				pthread_mutex_unlock(&mutexjoblist);

				if (jobcount>0 && !engineFailed) {

						/* there is a job to be executed */
						if (matlabRunning==0) {
								DEBUG(LOG_CRIT, "starting MATLAB engine");
								if ((en = engOpen(startcmd)) == NULL) {
										/* this may be due to a licensing problem */
										/* do not attempt to start again during the timeout period */
										DEBUG(LOG_ERR, "failed to start MATLAB engine, deleting job and switching to zombie");
										engineFailed = time(NULL);
										pthread_mutex_lock(&mutexhost);
										host->status = STATUS_ZOMBIE;
										pthread_mutex_unlock(&mutexhost);

										/* remove the failed job from the joblist */
										pthread_mutex_lock(&mutexjoblist);
										job = joblist;
										joblist = job->next; 
										FREE(job->job);
										FREE(job->host);
										FREE(job->arg);
										FREE(job->opt);
										FREE(job);
										pthread_mutex_unlock(&mutexjoblist);
										continue;
								}
								else {
										matlabRunning = 1;
								}
						}

						/* get the first job input arguments and options */
						pthread_mutex_lock(&mutexjoblist);
						job = joblist;

						pthread_mutex_lock(&mutexhost);
						/* switch the mode to busy slave */
						host->status = STATUS_BUSY;
						/* determine the maximum allowed job duration */
						timallow = 2*(host->timavail+1);
						/* update the current job description */
						bzero(&(host->current), sizeof(current_t));
						// FIXME
						host->current.pid = getpid();
						host->current.id  = job->host->id;
						strncpy(host->current.name, job->host->name, STRLEN);
						strncpy(host->current.user, job->host->user, STRLEN);
						strncpy(host->current.group, job->host->group, STRLEN);
						host->current.timreq  = job->job->timreq;
						host->current.memreq  = job->job->memreq;
						host->current.cpureq  = job->job->cpureq;

						pthread_mutex_unlock(&mutexhost);

						matlabStart = time(NULL);

						/* inform the other peers of the updated status */
						announce_once();

						argin   = (mxArray *)mxDeserialize(job->arg, job->job->argsize);
						options = (mxArray *)mxDeserialize(job->opt, job->job->optsize);
						jobid   = job->job->id;
						peerid  = job->host->id;
						DEBUG(LOG_CRIT, "executing job %d from %s@%s (jobid=%u, memreq=%lu, timreq=%lu)", ++jobnum, job->host->user, job->host->name, job->job->id, job->job->memreq, job->job->timreq);
						pthread_mutex_unlock(&mutexjoblist);

						/* create a copy of the optin cell-array */
						int n, i;
						n = mxGetM(options) * mxGetN(options);
						mxArray *previous = options;
						options = mxCreateCellMatrix(1, n+4);
						for (i=0; i<n; i++)
								mxSetCell(options, i, mxGetCell(previous, i));
						/* add the masterid and timallow options, these are used by peerexec for the killswitch */
						mxSetCell(options, n+0, mxCreateString("masterid\0"));
						mxSetCell(options, n+1, mxCreateDoubleScalar(peerid));
						mxSetCell(options, n+2, mxCreateString("timallow\0"));
						mxSetCell(options, n+3, mxCreateDoubleScalar(timallow));

						/* copy them over to the engine */
						engPutVariable(en, "argin", argin);
						engPutVariable(en, "options", options);
						mxDestroyArray(argin);
						mxDestroyArray(options);

						/* execute the job */
						engEvalString(en, "[argout, options] = peerexec(argin, options);");
						jobfailed = 0;

						/* get the job output arguments and options */
						if ((argout = engGetVariable(en, "argout")) == NULL) {
								DEBUG(LOG_ERR, "error getting argout");
								jobfailed = 1;
						}

						if ((options = engGetVariable(en, "options")) == NULL) {
								DEBUG(LOG_ERR, "error getting options");
								jobfailed = 1;
						}

						if (jobfailed) {
								/* create an empty ouptut argument */
								argout = mxCreateCellMatrix(1,1);
								/* specify the error in the options */
								options = mxCreateCellMatrix(1,2);
								mxSetCell(options, 0, mxCreateString("lasterr\0"));
								mxSetCell(options, 1, mxCreateString("the matlab engine aborted\0"));
								/* try to close the engine, probably this will fail */
								engClose(en);
								matlabRunning = 0;
						};

						/* send the results back to the master */

						/*****************************************************************************
						 * the following code is largely shared with the put-option in the peer mex file
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
								DEBUG(LOG_ERR, "failed to locate specified peer");
								goto cleanup;
						}

						pthread_mutex_lock(&mutexhost);
						int hasuds = (strlen(peer->host->socket)>0 && strcmp(peer->host->name, host->name)==0);
						int hastcp = (peer->host->port>0);
						pthread_mutex_unlock(&mutexhost);

						/* these have to be initialized to NULL to ensure that the cleanup works */
						server = 0;
						def = NULL;
						arg = NULL;
						opt = NULL;

						if (hasuds) {
								/* open the UDS socket */
								if ((server = open_uds_connection(peer->host->socket)) < 0) {
										pthread_mutex_unlock(&mutexpeerlist);
										DEBUG(LOG_ERR, "failed to create socket");
										goto cleanup;
								}
						}
						else if (hastcp) {
								/* open the TCP socket */
								if ((server = open_tcp_connection(peer->ipaddr, peer->host->port)) < 0) {
										pthread_mutex_unlock(&mutexpeerlist);
										DEBUG(LOG_ERR, "failed to create socket");
										goto cleanup;
								}
						}
						else {
								pthread_mutex_unlock(&mutexpeerlist);
								DEBUG(LOG_ERR, "failed to create socket");
								goto cleanup;
						}

						/* the connection was opened without error */
						pthread_mutex_unlock(&mutexpeerlist);

						if ((n = bufread(server, &handshake, sizeof(int))) != sizeof(int)) {
								DEBUG(LOG_ERR, "could not write handshake");
								goto cleanup;
						}
						else if (!handshake) {
								DEBUG(LOG_ERR, "failed to negociate connection");
								goto cleanup;
						}

						/* the message that will be written consists of
						   message->host
						   message->job
						   message->arg
						   message->opt
						 */

						if ((arg = (mxArray *) mxSerialize(argout))==NULL) {
								DEBUG(LOG_ERR, "could not serialize job arguments");
								goto cleanup;
						}

						if ((opt = (mxArray *) mxSerialize(options))==NULL) {
								DEBUG(LOG_ERR, "could not serialize job options");
								goto cleanup;
						}

						if ((def = (jobdef_t *)malloc(sizeof(jobdef_t)))==NULL) {
								DEBUG(LOG_ERR, "could not allocate memory");
								goto cleanup;
						}

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
								DEBUG(LOG_ERR, "could not write handshake");
								goto cleanup;
						}
						else if (!handshake) {
								DEBUG(LOG_ERR, "failed to write hostdef");
								goto cleanup;
						}

						if (success)
								success = (bufwrite(server, def, sizeof(jobdef_t)) == sizeof(jobdef_t));

						if ((n = bufread(server, &handshake, sizeof(int))) != sizeof(int)) {
								DEBUG(LOG_ERR, "could not write handshake");
								goto cleanup;
						}
						else if (!handshake) {
								DEBUG(LOG_ERR, "failed to write jobdef");
								goto cleanup;
						}

						if (success) 
								success = (bufwrite(server, (void *)mxGetData(arg), def->argsize) == def->argsize);

						if ((n = bufread(server, &handshake, sizeof(int))) != sizeof(int)) {
								DEBUG(LOG_ERR, "could not write handshake");
								goto cleanup;
						}
						else if (!handshake) {
								DEBUG(LOG_ERR, "failed to write arg");
								goto cleanup;
						}

						if (success) 
								success = (bufwrite(server, (void *)mxGetData(opt), def->optsize) == def->optsize);

						if ((n = bufread(server, &handshake, sizeof(int))) != sizeof(int)) {
								DEBUG(LOG_ERR, "could not write handshake");
								goto cleanup;
						}
						else if (!handshake) {
								DEBUG(LOG_ERR, "failed to write opt");
								goto cleanup;
						}

cleanup:
						if (server>0) {
								close_connection(server);
								server = 0;
						}

						if (arg) {
								mxDestroyArray(arg);
								arg = NULL;
						}

						if (opt) {
								mxDestroyArray(opt);
								opt = NULL;
						}

						FREE(def);

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

						/* make the slave available again */
						pthread_mutex_lock(&mutexhost);
						host->status = STATUS_IDLE;
						/* update the current job description */
						bzero(&(host->current), sizeof(current_t));
						pthread_mutex_unlock(&mutexhost);

						/* inform the other peers of the updated status */
						announce_once();

						matlabFinished = time(NULL);
						DEBUG(LOG_CRIT, "executing job %d took %d seconds", jobnum, matlabFinished - matlabStart);
				}
				else {
						usleep(SLEEPTIME);
				} /* if jobcount */

				/* switch the engine off if it is idle for too long */
				if ((matlabRunning!=0) && ((time(NULL)-matlabFinished)>enginetimeout)) {
						if (engClose(en)!=0) {
								DEBUG(LOG_CRIT, "could not stop the MATLAB engine");
								matlabRunning = 0;
						}
						else {
								DEBUG(LOG_CRIT, "stopped idle MATLAB engine");
								matlabRunning = 0;
						}
				}

				/* switch back to the default state after having waiting some time */
				if ((engineFailed!=0) && ((time(NULL)-engineFailed)>zombietimeout)) {
						DEBUG(LOG_NOTICE, "switching back to idle mode");
						pthread_mutex_lock(&mutexhost);
						host->status = STATUS_IDLE;
						/* update the current job description */
						bzero(&(host->current), sizeof(current_t));
						pthread_mutex_unlock(&mutexhost);
						engineFailed = 0;
						continue;
				}

		} /* while */

		return 0;
} /* main */

