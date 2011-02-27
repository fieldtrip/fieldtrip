#include <stdio.h>
#include <string.h>

#include "peer.h"
#include "parser.h"

void backspace(char *s, int pos, int num) {
		int len_s, p;
		len_s = strlen(s);
		if (num > 0) {
				for (p = pos; p+num < len_s; p++) {
						s[p]=s[p+num];
				}
				s[len_s-num]='\0';
		}
}

void trimline(char *s) {
		/* trim spaces and tabs from beginning */
		int i = 0, len_s;
		len_s = strlen(s);
		while((s[i]==' ') || (s[i]=='\t')) {
				i++;
		}
		backspace(s,0,i);

		/* trim spaces and tabs from end */
		i = len_s - 1;
		while((s[i]==' ') || (s[i]=='\t') || (s[i]=='\r')  || (s[i]=='\n')) {
				i--;
		}
		if(i < (len_s - 1)) {
				s[i+1] = '\0';
		}
}

void initconfig(config_t *cconf) {
		if (cconf) {
				cconf->pid         = 0;
				cconf->next        = NULL;
				/* these should either be NULL or point to a string */
				cconf->memavail    = NULL;
				cconf->cpuavail    = NULL;
				cconf->timavail    = NULL;
				cconf->allowuser   = NULL;
				cconf->allowgroup  = NULL;
				cconf->allowhost   = NULL;
				cconf->refuseuser   = NULL;
				cconf->refusegroup  = NULL;
				cconf->refusehost   = NULL;
				cconf->group       = NULL;
				cconf->hostname    = NULL;
				cconf->matlab      = NULL;
				cconf->timeout     = NULL;
				cconf->smartshare  = NULL;
				cconf->smartmem    = NULL;
				cconf->smartcpu    = NULL;
				cconf->udsserver   = NULL;
				cconf->verbose     = NULL;
		}
}

/* parseline gets a textline and if it finds a key=value pair, 
 * the value will be put in newly allocated memory and 
 * the pointer to it returned, if not found, return NULL */
char* parseline(char* line, char* key) {
		char* keyp = NULL;
		char* keyis;
		char* valp;
		char* valcpy=NULL; /* will be return value */

		/* make a copy of the key combined with '=' */
		keyis = malloc((strlen(key) + 2) * sizeof(char));
		strcpy(keyis, key);
		strcat(keyis, "=");

		/* find the occurrance of the key in the line */
		/* copy the value after the '=' to a newly allocated string */
		if ( (keyp = strstr(line, keyis)) ) {
				valp = keyp + strlen(keyis);
				valcpy = malloc(sizeof(char)*strlen(valp));
				strcpy(valcpy, valp);
		}
		free(keyis);
		return valcpy;
}

/* parsefile takes a filename as input and returns a linked list
 * with the configuration for the peer slaves
 * config is a pointer to the start of the configuration list */
int parsefile(char* fname, config_t** config) {
		char line[LINELENGTH];
		int numread=0;
		FILE* cf;
		config_t *cconf=NULL, *pconf=NULL;

		if ( (cf = fopen(fname, "r")) ) {
				while ( fgets(line, LINELENGTH, cf) ) {

						trimline(line); /* remove leading/trailing spaces/tabs */
						if (*line == '#') continue; /* skip lines that start with a comment */
						if (*line == '%') continue; /* skip lines that start with a comment */

						/* make a new list item when a [peer] line is found */
						if ( strstr(line, "[peer]") ) {
								if (pconf == NULL) {
										cconf = (config_t *)malloc(sizeof(config_t));
										pconf = cconf;
								} else {
										cconf->next = (config_t *)malloc(sizeof(config_t));
										cconf = cconf->next;
								}
								initconfig(cconf);
						}
						else if (cconf) {
								/* fill the current configuration structure with specific key=value pairs */
								if (!cconf->memavail) 
										cconf->memavail    = parseline(line, "memavail");
								if (!cconf->cpuavail) 
										cconf->cpuavail    = parseline(line, "cpuavail");
								if (!cconf->timavail) 
										cconf->timavail    = parseline(line, "timavail");
								if (!cconf->allowhost) 
										cconf->allowhost   = parseline(line, "allowhost");
								if (!cconf->allowuser) 
										cconf->allowuser   = parseline(line, "allowuser");
								if (!cconf->allowgroup) 
										cconf->allowgroup  = parseline(line, "allowgroup");
								if (!cconf->group) 
										cconf->group       = parseline(line, "group");
								if (!cconf->hostname) 
										cconf->hostname    = parseline(line, "hostname");
								if (!cconf->matlab) 
										cconf->matlab      = parseline(line, "matlab");
								if (!cconf->timeout) 
										cconf->timeout     = parseline(line, "timeout");
								if (!cconf->smartshare) 
										cconf->smartshare  = parseline(line, "smartshare");
								if (!cconf->smartmem) 
										cconf->smartmem    = parseline(line, "smartmem");
								if (!cconf->smartcpu) 
										cconf->smartcpu    = parseline(line, "smartcpu");
								if (!cconf->udsserver) 
										cconf->udsserver   = parseline(line, "udsserver");
								if (!cconf->verbose) 
										cconf->verbose     = parseline(line, "verbose");
						}

				} /* while */
		} else {
				fprintf(stderr, "parser: cannot open file %s\n", fname);
				return 0;
		} /* if fopen */

		/* this will be returned to the calling function */
		*config = pconf;

		/* count the number of nodes */
		while (pconf) {
				numread++;
				pconf= pconf->next;
		}
		return numread;				
}

