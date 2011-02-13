#include <stdio.h>
#include <stdlib.h>

#include "peer.h"
#include "parser.h"

int main () {
        config_t *config=NULL;
        int i=0;

        if (parsefile("test_parser.conf", &config) ) {
                while (config) {
                        i++;
                        if (config->memavail) printf("%d: memavail %s\n", i, config->memavail);
                        if (config->timavail) printf("%d: timavail %s\n", i, config->timavail);
                        if (config->cpuavail) printf("%d: cpuavail %s\n", i, config->cpuavail);
                        config=config->next;
                }
        } else {
                fprintf(stderr, "Error reading test_parser.conf\n");
                exit(1);
        }
        exit(0);
}

