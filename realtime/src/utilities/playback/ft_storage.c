/*
 * Collection of routines for saving FieldTrip buffer data to disk.
 *
 * (C) 2010 S. Klanke
 */

#include <ft_storage.h>

#ifdef WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/time.h>
#endif

static char *datatype_names[]={"char","uint8","uint16","uint32","uint64","int8","int16","int32","int64","float32","float64"};

ft_storage_t *ft_storage_create(const char *directory, const headerdef_t *hdef, const void *chunks, int *errCode) {
	ft_storage_t *S;
	int r;
	char endianness[10];
	union {
		short word;
		char bytes[2];
	} endianTest;
	
	endianTest.bytes[0] = 1;
	endianTest.bytes[1] = 0;

	if (endianTest.word == 1) {
		strcpy(endianness, "little");
	} else {
		strcpy(endianness, "big");
	}

	S = (ft_storage_t *) calloc(1, sizeof(ft_storage_t));
	if (S==NULL) {
		if (errCode) *errCode=FT_OUT_OF_MEMORY;
		return NULL;
	}
	S->dirLen = strlen(directory);
	
	/* overallocate to have space for /header.txt and other names */
	S->dirName = (char *) malloc(S->dirLen + 16);
	if (S->dirName==NULL) {
		if (errCode) *errCode=FT_OUT_OF_MEMORY;
		goto cleanup;
	}
	strcpy(S->dirName, directory);

	if (errCode) *errCode = FT_FILE_ERROR;
	#ifdef WIN32
	r=mkdir(directory);
	#else
	r=mkdir(directory, 0700);
	#endif
	if (r==-1) {
		fprintf(stderr, "ERROR: cannot create directory %s\n", directory);
		goto cleanup;
	}

	strcpy(S->dirName + S->dirLen, "/header");
	S->fHeader = fopen(S->dirName, "wb");
	if (S->fHeader==NULL) goto cannotWrite;
	
	strcpy(S->dirName + S->dirLen, "/header.txt");
	S->fHeaderTxt = fopen(S->dirName, "w");
	if (S->fHeaderTxt==NULL) goto cannotWrite;
		
	strcpy(S->dirName + S->dirLen, "/samples");
	S->fSamples = fopen(S->dirName, "wb");
	if (S->fSamples == NULL) goto cannotWrite;
	
	strcpy(S->dirName + S->dirLen, "/events");
	S->fEvents = fopen(S->dirName, "wb");
	if (S->fEvents == NULL) goto cannotWrite;

	strcpy(S->dirName + S->dirLen, "/timing");
	S->fTime = fopen(S->dirName, "wb");
	if (S->fTime == NULL) goto cannotWrite;
	
	/* write binary header definition */
	fwrite(hdef, sizeof(headerdef_t), 1, S->fHeader);
	
	/* write ASCII header definition */
	fprintf(S->fHeaderTxt, "version=%i\n", VERSION);
	fprintf(S->fHeaderTxt, "endian=%s\n", endianness);
	fprintf(S->fHeaderTxt, "dataType=%s\n", datatype_names[hdef->data_type]);
	fprintf(S->fHeaderTxt, "fSample=%f\n", hdef->fsample);		
	fprintf(S->fHeaderTxt, "nChans=%i\n", hdef->nchans);

	if (hdef->bufsize>0 && chunks) {
		const ft_chunk_t *cnc;
		/* write all chunks to binary header */
		fwrite(chunks, 1, hdef->bufsize, S->fHeader);
		/* write channel name chunk to ASCII header */
		
		cnc = find_chunk(chunks, 0, hdef->bufsize, FT_CHUNK_CHANNEL_NAMES);
		if (cnc) {
			int i;
			const char *ni = (const char *) cnc->data;
			for (i=0;i<hdef->nchans;i++) {
				int n = strlen(ni);
				fprintf(S->fHeaderTxt, "%i:%s\n", i+1, ni);
				ni+=n+1;
			}
		}
	} 
	
	fflush(S->fHeaderTxt);
	
	S->sampleSize  = wordsize_from_type(hdef->data_type) * hdef->nchans;
	S->numChannels = hdef->nchans;
	S->numSamples  = 0;
	S->numEvents   = 0;
	S->created = 1;
	S->curSampleFile = 0;
	
	
	return S;
	
cannotWrite:
	fprintf(stderr, "ERROR: cannot create file %s\n", S->dirName);	
cleanup:
	if (S->fHeader)    fclose(S->fHeader);
	if (S->fHeaderTxt) fclose(S->fHeaderTxt);
	if (S->fSamples)   fclose(S->fSamples);
	if (S->fEvents)    fclose(S->fEvents);
	if (S->fTime)      fclose(S->fTime);
	if (S->dirName)    free(S->dirName);
	free(S);
	return NULL;
}

void ft_storage_close(ft_storage_t *S) {
	fclose(S->fSamples);
	fclose(S->fEvents);
	fclose(S->fTime);
	
	if (S->created) {
		/* update sample and event numbers */
		fseek(S->fHeader, 4, SEEK_SET);
		fwrite(&S->numSamples, 4, 1, S->fHeader);
		fwrite(&S->numEvents,  4, 1, S->fHeader);
		
		fprintf(S->fHeaderTxt, "nSamples=%i\n", S->numSamples);
		fprintf(S->fHeaderTxt, "nEvents=%i\n",  S->numEvents);
	}

	fclose(S->fHeader);
	fclose(S->fHeaderTxt);
	
	free(S->dirName);
	free(S);
}

int ft_storage_add_samples(ft_storage_t *S, int numSamples, const void *data) {
	long fPos = ftell(S->fSamples);
	long addSize = S->sampleSize * numSamples;
	
	if ((unsigned long) fPos + (unsigned long) addSize > (unsigned long) 2*1024*1024*1024) {
		/* would violate 2 GB boundary - create new file instead */
		fclose(S->fSamples);
		
		S->curSampleFile++;
		
		sprintf(S->dirName + S->dirLen, "/samples%i", S->curSampleFile);
		S->fSamples = fopen(S->dirName, "wb");
		if (S->fSamples == NULL) {
			fprintf(stderr, "ERROR: cannot create file %s\n", S->dirName);
			return FT_FILE_ERROR;
		}
	}
	
	if (fwrite(data, S->sampleSize, numSamples, S->fSamples) != numSamples) return FT_FILE_ERROR;
	S->numSamples += numSamples;
	fflush(S->fSamples);
	
	return 0;
}


int ft_storage_add_events(ft_storage_t *S, int size, const void *events) {
	if (fwrite(events, 1, size, S->fEvents) != size) return FT_FILE_ERROR;
	fflush(S->fEvents);
	return 0;
}

int ft_storage_add_event(ft_storage_t *S, const eventdef_t *event, const void *type, const void *value) {
	UINT32_T siz = sizeof(eventdef_t);
	if (fwrite(event, 1, siz, S->fEvents) != siz) return FT_FILE_ERROR;
	
	siz = wordsize_from_type(event->type_type) * event->type_numel;
	if (fwrite(type, 1, siz, S->fEvents) != siz) return FT_FILE_ERROR;
	
	siz = wordsize_from_type(event->value_type) * event->value_numel;
	if (fwrite(value, 1, siz, S->fEvents) != siz) return FT_FILE_ERROR;
	fflush(S->fEvents);
	return 0;
}

int ft_storage_add_timing(ft_storage_t *S, const ft_timing_element_t *te) {
	if (te->numEvents > 0) {
		if (fprintf(S->fTime, "E %i %f\n", te->numEvents, te->time) < 0) return FT_FILE_ERROR;
	}
	if (te->numSamples > 0) {
		if (fprintf(S->fTime, "S %i %f\n", te->numSamples, te->time) < 0) return FT_FILE_ERROR;
	}
	return 0;
}

	
