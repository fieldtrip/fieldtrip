/*
 * Copyright (C) 2008, Robert Oostenveld
 * F.C. Donders Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * $Log: buffer_gethdr.c,v $
 * Revision 1.9  2009/06/30 11:35:48  roboos
 * pass more explicit error from buffer onto matlab
 *
 * Revision 1.8  2008/10/29 20:46:12  roboos
 * consistent use of open_connection and close_connection
 * there were some incorrect uses of close() that on windows did not actually close, resulting in the buffer running out of sockets/threads after prolonged use
 *
 * Revision 1.7  2008/07/09 13:34:21  roboos
 * small change in verbose output, using verbose=0|1
 *
 * Revision 1.6  2008/05/22 09:27:22  roboos
 * fixed some issues with Borland compiler, correct pointer arithmetic, declarations at beginning
 *
 * Revision 1.5  2008/03/09 22:36:38  roboos
 * implemented new buffer mex-function, which acts as dispatcher for all functionality
 * changed buffer_xxx code from mex-function into c-function to be included in dispatcher
 *
 * Revision 1.4  2008/02/26 21:43:25  roboos
 * renamed packet_t structure definition into messagedef_t, added message_t structure (contains def+buf)
 *
 * Revision 1.3  2008/02/20 13:49:17  roboos
 * changed somments into ansi style, needed for matlab
 *
 * Revision 1.2  2008/02/19 10:28:33  roboos
 * added copyright statement and log message
 *
 *
 */

#include "mex.h"
#include "matrix.h"
#include "buffer.h"

#define NUMBER_OF_FIELDS 6
#define MAX_NUM_BLOBS  32
#define SIZE_NIFTI_1   348

/* return field number on success, or -1 if already existing or invalid */
int addIfNew(mxArray *S, const char *name) {
	int field;
	if (mxGetFieldNumber(S, name) >= 0) {
		printf("Chunk '%s' already defined. Skipping.\n", name);
		return -1;
	}
	field = mxAddField(S, name);
	if (field<0) {
		printf("Can not add '%s' to output struct. Skipping.\n", name);
	}
	return field;
}

mxArray *channelNames2Cell(const char *str, int len, int numChannels) {
	int i,pe,ps = 0;
	mxArray *A = mxCreateCellMatrix(numChannels, 1);
	for (i=0;i<numChannels;i++) {
		mxArray *name;
		
		for (pe=ps; pe<len; pe++) {
			if (str[pe]==0) break;
		}
		if (pe>=len) {
			printf("Invalid name for channel %i. Skipping the rest.\n", i+1);
			break;
		}
		
		name = mxCreateString(str + ps);
		mxSetCell(A, i, name);
		/* next channel name begins after the 0 of the previous one */
		ps = pe+1;
	}
	return A;
}

void addChunksToMatrix(mxArray *S, const char *buf, int bufsize, int numChannels) {
	int bufpos = 0;
	int numBlobs = 0;
	mxArray *blobs[MAX_NUM_BLOBS];
	mxArray *keyval = NULL;
	mxArray *A;
	int field;

	while (bufpos + sizeof(ft_chunkdef_t) <= bufsize) {
		ft_chunk_t *chunk = (ft_chunk_t *) (buf + bufpos);
				
		/* "chunk" now points to the right location, make sure it has a valid size definition */
		if (bufpos + sizeof(ft_chunkdef_t) + chunk->def.size > bufsize) {
			printf("Invalid chunk size (%i) in Fieldtrip header detected. Stopping to parse.\n", chunk->def.size);
			break;
		}
				
		switch (chunk->def.type) {
			case FT_CHUNK_CHANNEL_NAMES:
				field = addIfNew(S, "channel_names");
				if (field < 0) break;
				A  = channelNames2Cell(chunk->data, chunk->def.size, numChannels);
				mxSetFieldByNumber(S, 0, field, A);
				break;
			case FT_CHUNK_NIFTI1:
				if (chunk->def.size != SIZE_NIFTI_1) {
					mexWarnMsgTxt("Invalid NIFTI-1 chunk detected. Skipping.");
					break;
				}
				field = addIfNew(S, "nifti_1");
				if (field < 0) break;
				/* pass on as 348 bytes (uint8), should be decoded on Matlab level (?) */
				A  = mxCreateNumericMatrix(1, SIZE_NIFTI_1, mxUINT8_CLASS, mxREAL);
				memcpy(mxGetData(A), chunk->data, SIZE_NIFTI_1);
				mxSetFieldByNumber(S, 0, field, A);
				break;
			case FT_CHUNK_SIEMENS_AP:
				field = addIfNew(S, "siemensap");
				if (field < 0) break;
				/* pass on as uint8, should be decoded on Matlab level (?) */
				A  = mxCreateNumericMatrix(1, chunk->def.size, mxUINT8_CLASS, mxREAL);
				memcpy(mxGetData(A), chunk->data, chunk->def.size);
				mxSetFieldByNumber(S, 0, field, A);
				break;
			case FT_CHUNK_CTF_RES4:
				field = addIfNew(S, "ctf_res4");
				if (field < 0) break;
				/* pass on as uint8, should be decoded on Matlab level (?) */
				A  = mxCreateNumericMatrix(1, chunk->def.size, mxUINT8_CLASS, mxREAL);
				memcpy(mxGetData(A), chunk->data, chunk->def.size);
				mxSetFieldByNumber(S, 0, field, A);
				break;
			case FT_CHUNK_RESOLUTIONS:
				field = addIfNew(S, "resolutions");
				if (field >=0) {
					int nc = chunk->def.size / sizeof(double);
					/*  If the chunk is buggy and there are less resolution values present,
						we only fill in those we have. If there are more, we only fill in numChannels.
						So the returned 'resolutions' field will always match the number of channels in the buffer.
					*/
					if (nc>numChannels) nc = numChannels;
					A = mxCreateDoubleMatrix(numChannels, 1, mxREAL);
					memcpy(mxGetPr(A), chunk->data, nc*sizeof(double));
				}
				break;
			case FT_CHUNK_UNSPECIFIED:
			default:
				if (numBlobs < MAX_NUM_BLOBS) {
					/* pass on the binary(?) blob as an uint8 matrix */
					A = mxCreateNumericMatrix(chunk->def.size, (chunk->def.size>0)?1:0, mxUINT8_CLASS, mxREAL);
					memcpy(mxGetData(A), chunk->data, chunk->def.size);
					blobs[numBlobs++] = A;
				} else {
					mexWarnMsgTxt("Encountered too many unspecified chunks in header. Skipping this one.");
				}
		}
		/* jump to next chunk */
		bufpos += chunk->def.size + sizeof(ft_chunkdef_t);
	}
	
	if (numBlobs > 0) {
		int i;
		
		field = addIfNew(S, "unspecified_blob");
		if (field < 0) return;
		
		A = mxCreateCellMatrix(numBlobs,1);
		for (i=0;i<numBlobs;i++) {
			mxSetCell(A, i, blobs[i]);
		}
		mxSetFieldByNumber(S, 0, field, A);
	}
}


int buffer_gethdr(int server, mxArray *plhs[], const mxArray *prhs[])
{
	int verbose = 0;
	int result  = 0;

	message_t *request  = NULL;
	message_t *response = NULL;

	/* this is for the Matlab specific output */
	const char *field_names[NUMBER_OF_FIELDS] = {"nchans", "nsamples", "nevents", "fsample", "data_type", "bufsize"};

	/* allocate the elements that will be used in the communication */
	request      = malloc(sizeof(message_t));
	request->def = malloc(sizeof(messagedef_t));
	request->buf = NULL;
	request->def->version = VERSION;
	request->def->command = GET_HDR;
	request->def->bufsize = 0;
	
	if (verbose) print_request(request->def);
	result = clientrequest(server, request, &response);
	
	if (result == 0) {
		if (verbose) print_response(response->def);

		if (response->def->command==GET_OK) {
			headerdef_t *headerdef = (headerdef_t *) response->buf;
			
			if (verbose) print_headerdef(headerdef);

			plhs[0] = mxCreateStructMatrix(1, 1, NUMBER_OF_FIELDS, field_names);
			mxSetFieldByNumber(plhs[0], 0, 0, mxCreateDoubleScalar((double)(headerdef->nchans)));
			mxSetFieldByNumber(plhs[0], 0, 1, mxCreateDoubleScalar((double)(headerdef->nsamples)));
			mxSetFieldByNumber(plhs[0], 0, 2, mxCreateDoubleScalar((double)(headerdef->nevents)));
			mxSetFieldByNumber(plhs[0], 0, 3, mxCreateDoubleScalar((double)(headerdef->fsample)));
			mxSetFieldByNumber(plhs[0], 0, 4, mxCreateDoubleScalar((double)(headerdef->data_type)));
			mxSetFieldByNumber(plhs[0], 0, 5, mxCreateDoubleScalar((double)(headerdef->bufsize)));
			
			addChunksToMatrix(plhs[0], (const char *) response->buf + sizeof(headerdef_t), headerdef->bufsize, headerdef->nchans);
		}
		else {
			result = response->def->command;
		}
	}

	if (response) {
		FREE(response->def);
		FREE(response->buf);
		FREE(response);
	}

	if (request) {
		FREE(request->def);
		FREE(request->buf);
		FREE(request);
	}

	return result;
}

