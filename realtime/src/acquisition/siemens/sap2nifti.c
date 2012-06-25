/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */

#include <siemensap.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <nifti1.h>
#include <math.h>

/* TODO: not everything get's translated yet. Look at spatial transformation and slice time */
int sap_to_nifti1(const sap_item_t *list, nifti_1_header *NH) {
	const sap_item_t *item;
	int numFound = 0;
	double phaseFOV = 0.0, readoutFOV = 0.0;

	if (list == NULL || NH == NULL) return -1;
	
	memset(NH, 0, sizeof(nifti_1_header));
	
	NH->sizeof_hdr = sizeof(nifti_1_header);
	/* number of dimensions = 3 (should be 4 for time series) */
	NH->dim[0] = 3;
	
	item = sap_search_deep(list, "sKSpace.lBaseResolution");
	if (item!=NULL && item->type == SAP_LONG) {
		long res = *((long *) item->value);
		if (res>0) {
			NH->dim[1] = (short) res;
			numFound++;
		}
	}
	
	item = sap_search_deep(list, "sSliceArray.lSize");
	if (item!=NULL && item->type == SAP_LONG) {
		long slices = *((long *) item->value);
		if (slices > 0) {
			NH->dim[3] = (short) slices;
			numFound++;
		}
	}
	
	item = sap_search_deep(list, "sSliceArray.asSlice[0].dPhaseFOV");
	if (item!=NULL && item->type == SAP_DOUBLE) {
		phaseFOV = *((double *) item->value);
	}
	
	item = sap_search_deep(list, "sSliceArray.asSlice[0].dReadoutFOV");
	if (item!=NULL && item->type == SAP_DOUBLE) {
		readoutFOV = *((double *) item->value);
	}
	
	if (phaseFOV > 0.0 && readoutFOV > 0.0) {
		float pixelSize = (float) (readoutFOV / NH->dim[1]);		
		
		/* number of pixels in phase direction */
		NH->dim[2] = (short) round(phaseFOV / pixelSize);
		NH->pixdim[1] = pixelSize;
		NH->pixdim[2] = pixelSize;
		numFound += 3;
	}
	
	item = sap_search_deep(list, "sSliceArray.asSlice[0].dThickness");
	if (item!=NULL && item->type == SAP_DOUBLE) {
		NH->pixdim[3] = (float) *((double *) item->value);
		numFound++;
	}
	
	NH->dim_info = FPS_INTO_DIM_INFO(1,2,3);
	NH->magic[0] = 'n';
	NH->magic[1] = 'i';
	NH->magic[2] = '1';
	NH->magic[3] = 0;	
	
	return numFound;
}
