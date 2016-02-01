/*
 * Real-time proxy between the Elekta Neuromag MEG system and FieldTrip real-time buffer
 *
 * (C)opyright Gustavo Sudre and Lauri Parkkonen, 2010
 *
 * This software comes without warranty of any kind, and it may not be fit for any particular
 * purpose. Use at your own risk. This software is NOT by Elekta Oy, and Elekta Oy does not endorse
 * its use.
 */

#include <fiff.h>
#include <buffer.h>
#include "neuromag2ft.h"

/* --------------------------------------------------------------------------
 * Process one data buffer. 
 * 
 * Do your processing here. Note that if it takes too long, a buffer
 * overflow may occur.
 */

void process_data(float **scannerData, int nchan, int nsamp, fiffChInfo *ch_info)
{

  int ch = 0, samp;

  //printf("Data buffer: %d channels, %d samples\n", nchan, nsamp);

  // writing data to fiedltrip buffer
  data_t       *data     = NULL;
  message_t    *request  = NULL;
  message_t    *response = NULL;
  int		status = 0;

  // allocate the elements that will be used in the communication
  data      = malloc(sizeof(data_t));
  data->def = malloc(sizeof(datadef_t));
  data->buf = NULL;

  // define the constant part of the data and allocate space for the variable part 
  data->def->nchans    = nchan;
  data->def->nsamples  = nsamp;
  data->def->data_type = DATATYPE_FLOAT32;
  data->def->bufsize   = WORDSIZE_FLOAT32*nchan*nsamp;
  FREE(data->buf);
  data->buf            = malloc(WORDSIZE_FLOAT32*nchan*nsamp);

  // the buffer is organized with the first sample for all channels, then the second
  // sample for all channels, etc
  for (samp = 0; samp < nsamp; samp++) {
    for (ch = 0; ch < nchan; ch++)
      ((FLOAT32_T *)(data->buf))[samp*nchan+ch] = scannerData[ch][samp];
  }

  // create the request
  request      = malloc(sizeof(message_t));
  request->def = malloc(sizeof(messagedef_t));
  request->buf = NULL;
  request->def->version = VERSION;
  request->def->bufsize = 0;
  request->def->command = PUT_DAT;
  request->def->bufsize = append(&request->buf, request->def->bufsize, data->def, sizeof(datadef_t));
  request->def->bufsize = append(&request->buf, request->def->bufsize, data->buf, data->def->bufsize);

  status = clientrequest(fieldtrip_sock, request, &response);

  if (status) {
    fprintf(stderr, "Something wrong with FieldTrip buffer while writing.\n");
  }

  // FIXME do someting with the response, i.e. check that it is OK 
  free(request->def);
  request->def = NULL;
  free(request->buf);
  request->buf = NULL;
  free(request);
  request = NULL;
  free(data->def);
  data->def = NULL;
  free(data->buf);
  data->buf = NULL;
  free(data);	
  data = NULL;
  free(response);
  response = NULL;
}
