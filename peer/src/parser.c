#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "peer.h"

/* 
 *  number      = number, number of slaves to start        (default = 1)
 *  memavail    = number, amount of memory available       (default = inf)
 *  cpuavail    = number, speed of the CPU                 (default = inf)
 *  timavail    = number, maximum duration of a single job (default = inf)
 *  allowhost   = {...}
 *  allowuser   = {...}
 *  allowgroup  = {...}
 *  group       = string
 *  hostname    = string
 *  matlab      = string
 *  timeout     = number, time to keep the engine running after the job finishes
 *  smartshare  = 0|1
 *  smartmem    = 0|1
 *  smartcpu    = 0|1
 *  udsserver   = 0|1
 *  verbose     = number
 */

typedef struct {
		int option;
		int memavail;
		int cpuavail;
		int timavail;
		int allowhost;
		int allowuser;
		int allowgroup;
		int group;
		int hostname;
		int matlab;
		int timeout;
		int smartshare;
		int smartmem;
		int smartcpu;
		int udsserver;
} setup_t;

void assign_default(setup_t *setup) {
		setup->option     = 0; 
		setup->memavail   = 0; 
		setup->cpuavail   = 0; 
		setup->timavail   = 0; 
		setup->allowhost  = 0; 
		setup->allowuser  = 0; 
		setup->allowgroup = 0;
		setup->group      = 0; 
		setup->hostname   = 0; 
		setup->matlab     = 0; 
		setup->timeout    = 0; 
		setup->smartshare = 0;
		setup->smartmem   = 0; 
		setup->smartcpu   = 0; 
		setup->udsserver  = 0; 
}

void cleanup_line(char *line) {
		char *ptr1 = NULL, *ptr2 = NULL;

		/* remove the newline and any comments*/
		for (ptr1=line; *ptr1; ptr1++) 
				if (*ptr1=='\r' || *ptr1=='\n' || *ptr1=='#' || *ptr1=='%')
						*ptr1 = 0;

		/* convert all characters up to the '=' to lower case */
		ptr1=line;
		while (*ptr1) {
				*ptr1 = tolower(*ptr1);
				ptr1++;
		}

		/* remove all spaces and tabs */
		ptr1 = line;
		ptr2 = line;
		for (ptr1=line; *ptr1; ptr1++) {
				while (*ptr1==' ')
						ptr1++;
				*ptr2 = *ptr1;
				ptr2++;
		}
		*ptr2=0;
}

int parser(const char *filename, setup_t **setup, int maxcount) {
		char line[256];
		int count = 0;
		FILE *fp;

		if ((fp = fopen(filename, "r"))==0) {
				perror(filename);
				return(-1);
		}

		while (!feof(fp) && count<maxcount) {
				if (fgets(line, 256, fp)!=NULL) {

						/* remove spaces etcetera */
						cleanup_line(line);

						if (strlen(line)==0)
								/* empty lines are not interesting */
								continue;

						if (strcmp(line, "[peer]")==0) {
								setup[count] = malloc(sizeof(setup_t));
								assign_default(setup[count]);
								count++;
						}

						if (count>0) {
								sscanf(line, "option=%d"     , &(setup[count-1]->option    ));
								sscanf(line, "memavail=%d"   , &(setup[count-1]->memavail  ));
								sscanf(line, "cpuavail=%d"   , &(setup[count-1]->cpuavail  ));
								sscanf(line, "timavail=%d"   , &(setup[count-1]->timavail  ));
								sscanf(line, "allowhost=%d"  , &(setup[count-1]->allowhost ));
								sscanf(line, "allowuser=%d"  , &(setup[count-1]->allowuser ));
								sscanf(line, "allowgroup=%d" , &(setup[count-1]->allowgroup));
								sscanf(line, "group=%d"      , &(setup[count-1]->group     ));
								sscanf(line, "hostname=%d"   , &(setup[count-1]->hostname  ));
								sscanf(line, "matlab=%d"     , &(setup[count-1]->matlab    ));
								sscanf(line, "timeout=%d"    , &(setup[count-1]->timeout   ));
								sscanf(line, "smartshare=%d" , &(setup[count-1]->smartshare));
								sscanf(line, "smartmem=%d"   , &(setup[count-1]->smartmem  ));
								sscanf(line, "smartcpu=%d"   , &(setup[count-1]->smartcpu  ));
								sscanf(line, "udsserver=%d"  , &(setup[count-1]->udsserver ));
						} /* if count */
				} /* if fgets */
		} /* while not feof */

		fclose(fp);
		return(count);
}

void main(void) {
		int count, i;
		setup_t *setup[128];

		count = parser("parser.conf", setup, 128);

		for (i=0; i<count; i++) {
				fprintf(stderr, "setup[%d].option  = %d\n", i, setup[i]->option);
				fprintf(stderr, "setup[%d].timeout = %d\n", i, setup[i]->timeout);
		}

		return;
}

