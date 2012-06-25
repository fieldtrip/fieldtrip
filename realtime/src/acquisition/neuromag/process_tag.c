/*
 * Real-time proxy between the Elekta Neuromag MEG system and FieldTrip real-time buffer
 *
 * (C)opyright Gustavo Sudre and Lauri Parkkonen, 2010
 *
 * This software comes without warranty of any kind, and it may not be fit for any particular
 * purpose. Use at your own risk. This software is NOT by Elekta Oy, and Elekta Oy does not endorse
 * its use.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* Elekta libraries */
#include <err.h>
#include <fiff.h>
#include <dacq.h>

/* FieldTrip real-time buffer */
#include <buffer.h>

/* Local */
#include "neuromag2ft.h"

static int nsamp             = 0;    /* Number of samples in a data buffer */
static int nchan             = 0;    /* Total number of channels */
static int ch_count          = 0;    /* General channel counter */
static fiffChInfo *ch_info   = NULL; /* Channel information */
static float sfreq           = -1.0; /* Sampling frequency (Hz) */
static int bufcnt            = 0;    /* How many buffers processed? */
static float **databuf       = NULL; /* Internal storage for one databuffer */
static int collect_data      = 1;    /* Really use the data? */

#ifdef FOO
static void set_data_filter(int state)
{
  static int filt_kind[] = {FIFF_DATA_BUFFER};
  if (state)
    dacq_set_data_filter (filt_kind,1);
  else 
    dacq_set_data_filter (NULL,0);
  return;
}
#endif

/* -------------------------------------------------------------------------
 * Get one data buffer, apply calibration coefficients and store the data
 */

static void catch_data(fiffTag tag,     /* Our data tag */ 
		       float **ch_data, /* (nchan x nsamp) data matrix */  
		       int nchan,       /* Total number of channels */
		       int nsamp)       /* Number of samples in data buffer */
{
  int ch, ns;
  float a;

  if (tag->type == FIFFT_DAU_PACK16) {
    fiff_dau_pack16_t *data16 = (fiff_dau_pack16_t *)tag->data;
    for (ch = 0; ch < nchan; ch++) {
      switch(ch_info[ch]->kind) {
      case FIFFV_MAGN_CH:
	if (ch_info[ch]->unit == FIFF_UNIT_T_M)
	  a = meg_grad_multiplier;
	else
	  a = meg_mag_multiplier;
	break;
      case FIFFV_EL_CH:
	a = eeg_multiplier;
	break;
      default:
	a = 1.0;
      }
      for (ns = 0; ns < nsamp; ns++)
	ch_data[ch][ns] = a * ch_info[ch]->cal * ch_info[ch]->range * data16[nchan*ns+ch];
    }
  }
  else {
    fiff_int_t *data32 = (fiff_int_t *)tag->data;
    for (ch = 0; ch < nchan; ch++) {
      switch(ch_info[ch]->kind) {
      case FIFFV_MAGN_CH:
	if (ch_info[ch]->unit == FIFF_UNIT_T_M)
	  a = meg_grad_multiplier;
	else
	  a = meg_mag_multiplier;
	break;
      case FIFFV_EL_CH:
	a = eeg_multiplier;
	break;
      default:
	a = 1.0;
      }
      for (ns = 0; ns < nsamp; ns++) 
	ch_data[ch][ns] = a * ch_info[ch]->cal * ch_info[ch]->range * data32[nchan*ns+ch];
    }
  }
}


/* -------------------------------------
 * Send the file header to the FT buffer
 */

int send_header_to_FT()
{
  int c, namelen = 0, status = 0;
  char *namevec = NULL, *name;
  message_t     *request  = NULL;
  message_t     *response = NULL;
  header_t      *header   = NULL;
  ft_chunk_t    *chunk    = NULL;

  dacq_log("Creating a header: %d channels, sampling rate %g Hz\n", nchan, sfreq);
  
  // Construct the channel name chunk
  for (c = 0; c < nchan; c++) {
    name = ch_info[c]->ch_name;
    namevec = realloc(namevec, namelen + strlen(name) + 1);
    strcpy(namevec + namelen, name);
    namelen = namelen + strlen(name) + 1;
  }
  chunk = malloc(sizeof(ft_chunkdef_t) + namelen * sizeof(char));
  chunk->def.type = FT_CHUNK_CHANNEL_NAMES;
  chunk->def.size = namelen;
  memcpy(chunk->data, namevec, namelen);  // I don't like this kind of assumptions on how structures are laid out in memory
  FREE(namevec);

  // Construct the header with the channel name chunk
  header      = malloc(sizeof(header_t));
  header->def = malloc(sizeof(headerdef_t));
  header->def->nchans    = nchan;
  header->def->nsamples  = 0;
  header->def->nevents   = 0;
  header->def->fsample   = sfreq;
  header->def->data_type = DATATYPE_FLOAT32;
  header->def->bufsize   = sizeof(ft_chunkdef_t) + namelen * sizeof(char);
  header->buf = chunk;

  // Construct the message with the header
  request      = malloc(sizeof(message_t));
  request->def = malloc(sizeof(messagedef_t));
  request->buf = NULL;
  request->def->version = VERSION;
  request->def->bufsize = 0;
  request->def->command = PUT_HDR;
  request->def->bufsize = append(&request->buf, request->def->bufsize, header->def, sizeof(headerdef_t));
  request->def->bufsize = append(&request->buf, request->def->bufsize, header->buf, header->def->bufsize);

  // Send the message to the buffer
  status = clientrequest(fieldtrip_sock, request, &response);
  if (status) {
    dacq_log("Something wrong with FieldTrip buffer during initialization, status %d. Exiting.\n", status);
    clean_up();
    exit(1);
  }
  dacq_log("Header sent to the FieldTrip buffer\n");
  
  // FIXME: Do someting with the response, i.e. check that it is OK 

  FREE(request->def);
  FREE(request);
  FREE(header->def);
  FREE(header);
  FREE(chunk);
  FREE(response->def);
  FREE(response);
  return(0);
}

/* -------------------------------------------------------------------------
 * Process one tag of the stream from the data server
 */

int process_tag (fiffTag tag)
{
  int block_kind;
  int c;
  fiffChInfo ch = NULL;
	
  switch (tag->kind) {
			
  case FIFF_NCHAN:               /* Number of channels */
    nchan = *(int *)tag->data;
    ch_info = calloc(nchan, sizeof(fiffChInfo));
    break;
			
  case FIFF_CH_INFO:             /* Information about one channel */
    ch = (fiffChInfo)(tag->data);
    ch_info[ch_count] = malloc(sizeof(fiffChInfoRec));
    memcpy(ch_info[ch_count], ch, sizeof(fiffChInfoRec));
    ch_count++;
    break;
			
  case FIFF_DATA_BUFFER:         /* One buffer of data */
    if (tag->type == FIFFT_DAU_PACK16)
      nsamp = nchan > 0 ? tag->size / (nchan * sizeof(fiff_dau_pack16_t)) : 0;
    else  /* FIFFT_INT */
      nsamp = nchan > 0 ? tag->size / (nchan * sizeof(fiff_int_t)) : 0;
    bufcnt++;
			
    if (collect_data) {
				
      if (databuf == NULL) {
	databuf = calloc(nchan, sizeof(float *));
	for (c = 0; c < nchan; c++)
	  databuf[c] = calloc(nsamp, sizeof(float));
      }
				
      catch_data(tag,databuf,nchan,nsamp);
      process_data(databuf,nchan,nsamp,ch_info);
				
    } else
      printf("Data buffer not sent\n");
			
    break;
			
  case FIFF_ERROR_MESSAGE : /* Message from the front-end */
    dacq_log("Error message from the data acquisition system: %s\n",tag->data);
    break;
			
  case FIFF_SFREQ :         /* Sampling frequency */
    if (tag->data)
      sfreq = *(float *)tag->data;
    break;
			
  case FIFF_BLOCK_START :
    block_kind = *(int *)(tag->data);
    switch (block_kind) {
					
    case FIFFB_MEAS :       /* Every file starts with this */
      dacq_log("New measurement is starting...\n");
      nchan = 0;
      sfreq = -1.0;
      /* set_data_filter (TRUE); */
      break;
					
    case FIFFB_RAW_DATA :   /* The data buffers start */
      dacq_log("Data buffers coming soon...\n");
      bufcnt = 0;
      send_header_to_FT();
      break;
    }
    break;
			
  case FIFF_BLOCK_END : 
    block_kind = *(int *)(tag->data);
    switch (block_kind) {
					
    case FIFFB_MEAS :
      dacq_log("Measurement ended (%d buffers).\n", bufcnt);
					
      // data was acquired, so we need to clear the data and the channel info
      if (databuf != NULL) {
	for (c = 0; c < nchan; c++) {
	  free(databuf[c]);
	  free(ch_info[c]);
	}
	free(databuf); databuf = NULL;
	free(ch_info); ch_info = NULL;
      }
      // channel info was created, but data was not acquired
      if (ch_info != NULL) {
	for (c = 0; c < nchan; c++)
	  free(ch_info[c]);
	free(ch_info); ch_info = NULL;
      }
					
      nchan = 0;
      nsamp = 0;
      bufcnt = 0;
      ch_count = 0;
      sfreq = -1.0;
      break;
    }
    break;
			
  case FIFF_CLOSE_FILE :
    dacq_log("File closed.\n");
    break;
			
  default:                  /* An unknown tag but it doesn't harm */
    break;
  }
	
  fflush(stdout);
  return(0);
}
