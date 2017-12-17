/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "buffer.h"

#define MERGE_THRESHOLD 4096 /* TODO: optimize this value? Maybe look at MTU size */

/*******************************************************************************
 * communicate with the buffer through TCP
 * returns 0 on success, -1 on error
 *******************************************************************************/
int tcprequest(int server, const message_t *request, message_t **response_ptr) {
  unsigned int n, total;

  /* this will hold the response */
  message_t *response;
  response      = (message_t*)malloc(sizeof(message_t));
  response->def = (messagedef_t*)malloc(sizeof(messagedef_t));
  DIE_BAD_MALLOC(response->def);
  response->buf = NULL;
  /* the response should be passed to the calling function, where it should be freed */
  *response_ptr = response;

  total = sizeof(messagedef_t) + request->def->bufsize;

  /* Check whether request->def and request->buf are already contiguous in memory,
     or whether request->buf is empty. If that's the case, we can write the request in one go.
   */
  if (request->def->bufsize == 0 || (request->def+1) == (messagedef_t *) request->buf) {
    if ((n = bufwrite(server, request->def, total)) != total) {
      fprintf(stderr, "write size = %u, should be %u\n", n, total);
      goto cleanup;
    }
  }
  /* Now check whether the total size is below the merge threshold, in which case
     we'll copy it to contiguous memory and again send it in one go
   */
  else if (total <= MERGE_THRESHOLD) {
    char merged[MERGE_THRESHOLD];

    memcpy(merged, request->def, sizeof(messagedef_t));
    memcpy(merged + sizeof(messagedef_t), request->buf, request->def->bufsize);

    if ((n = bufwrite(server, merged, total)) != total) {
      fprintf(stderr, "write size = %u, should be %u\n", n, total);
      goto cleanup;
    }
  }
  /* Otherwise, send "def" and "buf" in separate pieces. This might introduce latencies
     if the other end runs Windows :-(
   */
     else {
       /* send the request to the server, first the message definition */
       /* FIXME: bufwrite expects unsigned int, gets size_t. Similar for return. */
       if ((n = bufwrite(server, request->def, sizeof(messagedef_t)))!=sizeof(messagedef_t)) {
         fprintf(stderr, "write size = %u, should be %lu\n", n, sizeof(messagedef_t));
         goto cleanup;
       }

       /* send the request to the server, then the message payload */
       if ((n = bufwrite(server, request->buf, request->def->bufsize))!=request->def->bufsize) {
         fprintf(stderr, "write size = %d, should be %u\n", n, request->def->bufsize);
         goto cleanup;
       }
     }

     /* read the response from the server, first the message definition */
     if ((n = bufread(server, response->def, sizeof(messagedef_t))) != sizeof(messagedef_t)) {
       fprintf(stderr, "packet size = %d, should be %lu\n", n, sizeof(messagedef_t));
       goto cleanup;
     }

     if (response->def->version!=VERSION) {
       fprintf(stderr, "incorrect version\n");
       goto cleanup;
     }

     /* read the response from the server, then the message payload */
     if (response->def->bufsize>0) {
       response->buf = malloc(response->def->bufsize);
       if ((n = bufread(server, response->buf, response->def->bufsize)) != response->def->bufsize) {
         fprintf(stderr, "read size = %d, should be %d\n", n, response->def->bufsize);
         goto cleanup;
       }
     }

     /* everything went fine, return with the response */
     /* print_response(response->def); */
     return 0;

cleanup:
     /* there was a problem, clear the response and return */
     FREE(response->def);
     FREE(response->buf);
     FREE(response);
     *response_ptr = NULL;
     return -1;
}
