#include "free.h"

void free_header() {
	int verbose = 0;
	if (verbose>0) fprintf(stderr, "free_header: freeing header buffer\n");
	if (header) {
		FREE(header->def);
		FREE(header->buf);
		FREE(header);
	}
}


void free_data() {
	int verbose = 0;
	if (verbose>0) fprintf(stderr, "free_data: freeing data buffer\n");
	if (data) {
		FREE(data->def);
		FREE(data->buf);
		FREE(data);
	}
	thissample = 0;
	if (header) header->def->nsamples = 0;
}


void free_event() {
	int verbose = 0;
	int i;
	if (verbose>0) fprintf(stderr, "free_event: freeing event buffer\n");
	if (event) {
		for (i=0; i<MAXNUMEVENT; i++) {
			FREE(event[i].def);
			FREE(event[i].buf);
		}
		FREE(event);
	}
	thisevent = 0;
	if (header) header->def->nevents = 0;
}
