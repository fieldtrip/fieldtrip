#include "extern.h"


void init_data(void) {
	int verbose = 0;
	if (verbose>0) fprintf(stderr, "init_data: creating data buffer\n");
	if (header) {
		unsigned int wordsize = wordsize_from_type(header->def->data_type);

		if (wordsize==0) {
			fprintf(stderr, "init_data: unsupported data type (%u)\n", header->def->data_type);
			return;
		}
		/* heuristic of choosing size of buffer:
			 set current_max_num_sample to MAXNUMSAMPLE if nchans <= 256
			 otherwise, allocate about MAXNUMBYTE and calculate current_max_num_sample from nchans + wordsize
		 */
		if (header->def->nchans <= 256) {
			current_max_num_sample = MAXNUMSAMPLE;
		} else {
			current_max_num_sample = MAXNUMBYTE / (wordsize * header->def->nchans);
		}
		data = (data_t*)malloc(sizeof(data_t));

		DIE_BAD_MALLOC(data);

		data->def = (datadef_t*)malloc(sizeof(datadef_t));

		DIE_BAD_MALLOC(data->def);

		data->def->nchans    = header->def->nchans;
		data->def->nsamples  = current_max_num_sample;
		data->def->data_type = header->def->data_type;
		data->buf = malloc(header->def->nchans*current_max_num_sample*wordsize);

		DIE_BAD_MALLOC(data->buf);
	}
}


void init_event(void) {
	int verbose = 0;
	int i;
	if (verbose>0) fprintf(stderr, "init_event: creating event buffer\n");
	if (header) {
		event = (event_t*)malloc(MAXNUMEVENT*sizeof(event_t));
		DIE_BAD_MALLOC(event);
		for (i=0; i<MAXNUMEVENT; i++) {
			event[i].def = NULL;
			event[i].buf = NULL;
		}
	}
}
