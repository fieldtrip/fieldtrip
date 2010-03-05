/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "buffer.h"

int set_property(int server, const char *name, INT32_T *value) {
	int status = 0, verbose = 0;

	message_t *request   = NULL;
	message_t *response  = NULL;
	property_t *property = NULL;

	property                   = (property_t*)malloc(sizeof(property_t));
	property->def              = (propertydef_t*)malloc(sizeof(propertydef_t));
	property->buf              = NULL;
	property->def->type_type   = DATATYPE_CHAR;
	property->def->type_numel  = (UINT32_T) strlen(name)+1;
	property->def->value_type  = DATATYPE_INT32;
	property->def->value_numel = 1;
	property->def->bufsize     = 0;
	property->def->bufsize     = append(&property->buf, property->def->bufsize, (void*)name, property->def->type_numel*WORDSIZE_CHAR);
	property->def->bufsize     = append(&property->buf, property->def->bufsize, value, property->def->value_numel*WORDSIZE_INT32);

	if (verbose>0) print_propertydef(property->def);

	request               = (message_t*)malloc(sizeof(message_t));
	request->def          = (messagedef_t*)malloc(sizeof(messagedef_t));
	request->buf          = NULL;
	request->def->version = VERSION;
	request->def->command = PUT_PRP;
	request->def->bufsize = 0;
	request->def->bufsize = append(&request->buf, request->def->bufsize, property->def, sizeof(propertydef_t));
	request->def->bufsize = append(&request->buf, request->def->bufsize, property->buf, property->def->bufsize);

	if (verbose>0) print_request(request->def);
	if (verbose>0) print_buf(request->buf, request->def->bufsize);

	if ((status = clientrequest(server, request, &response)))
      exit(1);
    
	if (verbose>0) print_response(response->def);
	if (verbose>0) print_buf(response->buf, response->def->bufsize);

	FREE(property->def);
	FREE(property->buf);
	FREE(property);

	FREE(request->def);
	FREE(request->buf);
	FREE(request);

	if (response) {
		if (response->def->command==PUT_OK)
			status = 0;
		else {
			fprintf(stderr, "set_property: unexpected PUT_ERR\n");
		    status = -1;
        }
		FREE(response->def);
		FREE(response->buf);
		FREE(response);
	}
	else {
		fprintf(stderr, "set_property: unexpected response==NULL\n");
		status = -1;
	}

	return status;
}

int get_property(int server, const char *name, INT32_T *value) {
	int status = 0, verbose = 0;
	unsigned int offset;
	void *ptr_type, *ptr_value;

	message_t *request  = NULL;
	message_t *response = NULL;
	propertydef_t *propertydef;

	/* allocate the elements that will be used in the communication */
	request      = (message_t*)malloc(sizeof(message_t));
	request->def = (messagedef_t*)malloc(sizeof(messagedef_t));
	request->buf = NULL;
	request->def->version = VERSION;
	request->def->command = GET_PRP;
	request->def->bufsize = 0;

	if (verbose>0) print_request(request->def);
	if (verbose>0) print_buf(request->buf, request->def->bufsize);

	if ((status = clientrequest(server, request, &response)))
      exit(1);

	if (verbose>0) print_response(response->def);
	if (verbose>0) print_buf(response->buf, response->def->bufsize);

	if (response) {
		if (response->def->command == GET_OK) {

			status = -1; /* requested property has not (yet) been found */
			offset = 0;  /* this represents the offset of the property in the buffer */

			while (status && (offset<response->def->bufsize)) {

				propertydef = (propertydef_t*)((char*)response->buf+offset);
				if (verbose>1) print_propertydef(propertydef);

				if ((propertydef->type_type   == DATATYPE_CHAR)  &&
					(propertydef->type_numel  == strlen(name)+1)   && 
					(propertydef->value_type  == DATATYPE_INT32) &&
					(propertydef->value_numel == 1)) {

					ptr_type  = (char*)response->buf+offset+sizeof(propertydef_t);
					ptr_value = (char*)response->buf+offset+sizeof(propertydef_t) + propertydef->type_numel;
					if (strcmp((char *)ptr_type, name)==0) {
						/* property has been found */
						(*value) = *((INT32_T *)ptr_value);
						status = 0;
					}
				}
				if (status==-1)
					offset += sizeof(propertydef_t) + propertydef->bufsize;
			}

            if (verbose>0) {
				if (status==0)
					fprintf(stderr, "property \"%s\"has been found at offset %d, its value is %d\n", name, offset, *value);
				else
					fprintf(stderr, "property \"%s\"has not been found\n", name);
            }

			FREE(response->def);
			FREE(response->buf);
			FREE(response);
		}
		else {
			fprintf(stderr, "get_property: unexpected PUT_ERR\n");
			status = -2;
		}
	}
	else {
		fprintf(stderr, "get_property: unexpected empty response\n");
		status = -3;
	}

	return status;
}

