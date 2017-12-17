/*
 * Copyright (C) 2010, Stefan Klanke
 * 	Modified by Tim van Mourik 2015
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */

#include <mex.h>
#include <siemensap.h>
#include <string.h>

mxArray *createStructFromSAP(sap_item_t *item) {
	mxArray *A;
	
	if (item==NULL) {
		return mxCreateStructMatrix(0,0,0,NULL);
	}
	
	A = mxCreateStructMatrix(1,1,0,NULL);
	
	while (item!=NULL) {
		mxArray *F = NULL;
		int nr;
		
		nr = mxAddField(A, item->fieldname);
		if (nr==-1) {
			printf("Invalid field name '%s' or out of memory -- skipping item...\n", item->fieldname);
			item = item->next;
			continue;
		}
		
		switch(item->type) {
			case SAP_DOUBLE:
				F = mxCreateDoubleMatrix(item->num_elements, 1, mxREAL);
				memcpy(mxGetPr(F), item->value, item->num_elements*sizeof(double));
				break;
			case SAP_LONG:
				if (sizeof(long) == 4) {
					/* 32 bit machine */
					F = mxCreateNumericMatrix(item->num_elements, 1, mxINT32_CLASS, mxREAL);
				} else {
					F = mxCreateNumericMatrix(item->num_elements, 1, mxINT64_CLASS, mxREAL);
				}
				memcpy(mxGetPr(F), item->value, item->num_elements*sizeof(long));
				break;
			case SAP_TEXT:
				F = mxCreateString((char *) item->value);
				break;
			case SAP_STRUCT:
				if (item->is_array) {
					sap_item_t **children = (sap_item_t **) item->value;
					int i;
					/*  We need to use a cell array here, since we're not guaranteed to have 
						the same fields in each element
					*/
					F = mxCreateCellMatrix(item->num_elements,1);
					for (i=0;i<item->num_elements;i++) {
						mxArray *Fi = createStructFromSAP(children[i]);
						mxSetCell(F,i,Fi);
					}
				} else {
					sap_item_t **children = (sap_item_t **) item->value;
					F = createStructFromSAP(children[0]);
				}
				break;
		}
		if (F!=NULL) mxSetFieldByNumber(A,0,nr,F);
		
		item = item->next;
	}
	
	return A;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	char *buffer;
	int size;
	sap_item_t *L = NULL;
	
	if (nrhs!=1) mexErrMsgTxt("This function needs exactly one (string or uint8) argument.");
	
	size = (int) mxGetNumberOfElements(prhs[0]);
	
	if ((sizeof(char) == sizeof(mxChar) && mxIsChar(prhs[0])) || mxIsUint8(prhs[0]))  {
		buffer = (char *) mxGetData(prhs[0]);
		L = sap_parse(buffer, size);
	} else if (mxIsChar(prhs[0])) {
		buffer = mxArrayToString(prhs[0]);
		L = sap_parse(buffer, size);
		mxFree(buffer);
	} else {
		mexErrMsgTxt("Argument must be either of type 'char' or 'uint8'.");
	}
	
	sap_reverse_in_place(&L);
	
	plhs[0] = createStructFromSAP(L);
	
	sap_destroy(L);
}
