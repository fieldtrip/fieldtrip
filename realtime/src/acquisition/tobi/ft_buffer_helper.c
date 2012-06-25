#include "stdio.h"
#include "assert.h"
#include "buffer.h"


/* Starts FieldTrip buffer server on local host. This does *not* spawn an extra
 * thread.*/
void ft_buffer_serve(int port)
{
  host_t host;
  host.port = port;

  check_datatypes();  /* sanity check for sizes of different data types. */
  tcpserver((void *)(&host));  /* start server */
}


/* Construct a ft-buffer chunk containing channel names. Do not forget to free
 * the chunk. This frees the content of the chunk as well. 
 * FIXME: This depends on byte alignment of structs. Use proper serialization.
 */
ft_chunk_t *ft_create_chanlab_chunk(int n, const char **labels)
{
  ft_chunk_t *chunk = NULL;
  size_t payload_size = 0, chunk_size = 0;
  int i, l;
  char *p;

  /* Find out how big our payload and chunk should be. */
  for(i=0; i < n; ++i) 
    payload_size += strlen(labels[i]) + 1;
  chunk_size = payload_size + 8; /* sizeof(chunk) == 9 due to char .data[1]! */

  /* Allocate packet. */
  chunk = (ft_chunk_t *) malloc(chunk_size);
  if (!chunk) 
    return NULL;

  /* Fill final packet. */
  memset(chunk, 0, chunk_size);
  chunk->def.type = FT_CHUNK_CHANNEL_NAMES;
  chunk->def.size = payload_size;
  
  /* Copy label strings in buffer. */
  p = (char *) &(chunk->data);
  for (i = 0; i < n; ++i) {
    l = strlen(labels[i]) + 1;
    memcpy(p, labels[i], l);
    p += l;
  }

  return chunk;
}

/* Convenience function to call put header (PUT_DAT) in a FT buffer
 * FIXME: This depends on byte alignment of structs. Use proper serialization.
 */
int ft_put_hdr(int ft_buffer_handle, int nchann, int fsample, 
  ft_chunk_t *chunk)
{
  int status = 0;
  message_t req = {0}, *response = NULL;
  messagedef_t *req_hdr;
  headerdef_t *hdr;
  char *chunks, *packet;
  size_t packet_size = 0;

  /* Sanity check for sizes of different data types. */
  check_datatypes();  

  /* Find packet size */
  assert(sizeof(messagedef_t) == 8);
  assert(sizeof(headerdef_t) == 24);
  packet_size = sizeof(messagedef_t) + sizeof(headerdef_t);  /* fixed part */
  packet_size += 8 + chunk->def.size;  /* Note: sizeof(chunk) == 9! */

  /* Allocate packet */
  packet = (char *) malloc(packet_size);
  if (!packet) 
    return -1;
  memset(packet, 0, packet_size);

  /* Setup pointers */
  req_hdr = (messagedef_t *) packet;
  hdr = (headerdef_t *) ((char *) req_hdr + sizeof(messagedef_t));
  chunks = (char *) hdr + sizeof(headerdef_t);

  /* Fill content */
  req_hdr->version = VERSION;
  req_hdr->command = PUT_HDR;
  req_hdr->bufsize = packet_size - sizeof(messagedef_t);

  hdr->nchans = nchann;
  hdr->fsample = fsample;
  hdr->data_type = DATATYPE_FLOAT32;
  hdr->bufsize = req_hdr->bufsize - sizeof(headerdef_t);

  memcpy(chunks, chunk, hdr->bufsize);

  /* Request. */
  req.def = req_hdr;  /* bit of hack to get pointers into the same packet :/ */
  req.buf = hdr; 
  if (!clientrequest(ft_buffer_handle, &req, &response))
    status = -1;
  if (!response->def->command == PUT_OK)
    status = -2;

  /* Cleanup */
  free(packet);
  free(response->buf);
  free(response->def);
  free(response);

  return status;
}


/* Convenience wrapper to perform a PUT_DAT on ft_buffer. 
 * TODO: Add support for different datatypes?
 * FIXME: This depends on byte alignment of structs. Use proper serialization.
 */
int ft_put_data(int ft_buffer_handle, int nchannels, int nsamples, const float
  *chan_samp)
{
  int status = 0;
  message_t req = {0}, *response = NULL;
  messagedef_t req_hdr = {0};
  datadef_t data_hdr = {0};

  /* Create descriptor of raw data for FT-buffer: */
  data_hdr.nchans = nchannels;
  data_hdr.nsamples = nsamples;
  data_hdr.data_type = DATATYPE_FLOAT32;
  data_hdr.bufsize = nchannels * nsamples * sizeof(float);

  /* Create packet header for FT-buffer request */
  req_hdr.version = VERSION;
  req_hdr.command = PUT_DAT;
  req_hdr.bufsize = sizeof(data_hdr) + data_hdr.bufsize;

  /* Compose final FT-buffer request */
  req.def = &req_hdr;
  req.buf = malloc(req_hdr.bufsize);
  memset(req.buf, 0, req_hdr.bufsize);
  memcpy(req.buf, &data_hdr, sizeof(data_hdr));
  memcpy((char *) req.buf + sizeof(data_hdr), chan_samp, data_hdr.bufsize);

  if (!clientrequest(ft_buffer_handle, &req, &response))
    status = -1;
  if (response->def->command != PUT_OK)
    status = -2;

  free(req.buf);
  free(response->buf);
  free(response->def);
  free(response);

  return status;
}

