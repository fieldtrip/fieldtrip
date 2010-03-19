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


sap_item_t *sap_alloc_field(int len_name, const char *name, int size_value) {
	sap_item_t *F = malloc(sizeof(sap_item_t));
	F->fieldname = malloc(len_name+1);
	memcpy(F->fieldname, name, len_name);
	F->num_elements = 0;
	F->fieldname[len_name]=0;
	if (size_value == 0) {
		F->value = NULL;
	} else {
		F->value = calloc(size_value,1);
	}
	return F;
}


sap_item_t *sap_search_field(sap_item_t *first, int len_name, const char *name) {
	while (first!=NULL) {
		if (0==strncmp(first->fieldname, name, len_name)) return first;
		first = first->next;
	}
	return NULL;
}

/* parse one line and add to the list pointed to by "*first" */
int sap_handle_line(sap_item_t **first, const char *line, const char *line_end) {
	const char *name, *end_name;
	sap_field_type_t typ = SAP_LONG;
	sap_item_t *item;
	int len_name, size_value;
	long arrayIndex = 0;
	int isArray=0;
	
	long longVal;
	double doubleVal;
	const char *val_start,*val_end,*aux;
	
	name = line;
	while (*name == ' ') name++;
	/* check for comments and empty lines */
	if (name[0]=='#' || name[0]=='\r' || name[0]=='\n') return 0;
	
	end_name = name;
	/* detect field name */
	while (1) {	
		if (*end_name == '[') {
			isArray = 1;
			break;
		}
		if (*end_name == '.') {
			typ = SAP_STRUCT;
			break;
		}
		/* otherwise, name is finished at the first space or = sign */
		if (*end_name == ' ') break;
		if (*end_name == '=') break;
		if (*end_name == 0) {
			fprintf(stderr,"Something went wrong in sap_parse -- this should not happen.\n");
			return -1; 
		}
		end_name++;
	}
	
	len_name = end_name - name;
	
	if (isArray) {
		char *end_index;
		arrayIndex = strtol(end_name+1, (char **) &end_index, 10);
		if (arrayIndex < 0 || end_index[0]!=']' ) { /* this LONG_MIN in case of parse error */
			fprintf(stderr, "Invalid array index description.\n");
			return 0;
		}
		if (end_index[1]=='.') {
			typ = SAP_STRUCT;
			end_name = end_index + 1;
		}
	} else {
		arrayIndex = 0;
	}
	
	val_start = strchr(end_name, '=');
	if (val_start == NULL) {
		fprintf(stderr, "Invalid definition: No = sign\n");
		return 0;
	}
	val_start++; /* now points just behind the = sign */

	/* If it's not a struct, try to parse the value */
	if (typ != SAP_STRUCT) {
		/* Investigate right hand side to see whether it's a string, 
			a hexadecimal number (->long), or contains a . (-> double) */
		for (aux = val_start; aux!=line_end; aux++) {
			if (*aux == '"' || *aux == 'x' || *aux == '.') break;
		}
		
		if (aux==line_end) {
			/* If . is not found, we'll still go for a double if indicated by the name */
			if (name[0]=='d' || (name[0]=='a' && name[1]=='d')) {
				doubleVal = strtod(val_start, (char **) &val_end);
				if (val_start == val_end) {
					fprintf(stderr, "Field name indicated double precision value, but could not be parsed.\n");
					return 0;
				}
				size_value = sizeof(double);
				typ = SAP_DOUBLE;
				/* printf("%f  -> ", doubleVal); */
			} else {
				/* everything else is treated as a decimal long */
				longVal = strtol(val_start, (char **) &val_end, 10);
				if (longVal == LONG_MIN) {
					fprintf(stderr, "Could not parse long integer value as guessed by default.\n");
					return 0;
				}
				typ = SAP_LONG;
				size_value = sizeof(long);
				/* printf("%li  -> ", longVal);		*/
			}
		} else if (*aux=='"') {
			/* must be a string */
			val_start = aux+1;
			val_end = strchr(val_start, '"');
			if (val_end == NULL) {
				fprintf(stderr, "RHS started to look like a string, but could not be parsed.\n");
				return 0;
			}
			typ = SAP_TEXT;
			size_value = val_end - val_start + 1; /* +1 for trailing 0 */
		} else if (*aux=='x') {
			/* must be hexadecimal */
			val_start = aux+1;
			longVal = strtol(val_start, (char **) &val_end, 16);
			if (longVal == LONG_MIN) {
				fprintf(stderr, "RHS looked like a hexadecimal number, but could not be parsed.\n");
				return 0;
			}
			size_value = sizeof(long);
			typ = SAP_LONG;
			/* printf("%lX  -> ", longVal); */
		} else if (*aux=='.') {
			doubleVal = strtod(val_start, (char **) &val_end);
			if (val_start == val_end) {
				fprintf(stderr, "RHS looked like double precision value, but could not be parsed.\n");
				return 0;
			}
			size_value = sizeof(double);
			typ = SAP_DOUBLE;
			/* printf("%f  -> ", doubleVal); */
		} else {
			fprintf(stderr, "Parsing error\n");
			return 0;
		}
	} else {
		size_value = sizeof(void *); /* values are pointers */
	}
					
	item = sap_search_field(*first, len_name, name);
	if (item == NULL) {
		/* create new field */
		item = sap_alloc_field(len_name, name, size_value);
		item->next = *first;
		item->num_elements = 1;
		item->is_array = isArray;
		item->type = typ;
		*first = item;
		/* printf("New field: %.*s  %i, %i\n",len_name,name, isArray, typ); */
	} else {
		if (isArray) {
			if (arrayIndex >= item->num_elements) { /* need to reallocate, but we don't use realloc since we need to initialize to zero in case we jump over an index */
				void *newBuf = calloc(arrayIndex+1, size_value);
				memcpy(newBuf, item->value, item->num_elements*size_value);
				free(item->value);
				item->value = newBuf;
				item->num_elements = arrayIndex+1;
			}
		}
	}
	
	
	switch(typ) {
		case SAP_STRUCT:
			sap_handle_line(((sap_item_t **) item->value)+arrayIndex, end_name +1, line_end);
			break;
		case SAP_DOUBLE:
			((double *) item->value)[arrayIndex]=doubleVal;
			break;
		case SAP_LONG:
			((long *) item->value)[arrayIndex]=longVal;
			break;
		case SAP_TEXT:
			memcpy(item->value, val_start, size_value-1); 	/* only size_value-1 in string, but... */
			((char *)item->value)[size_value-1]=0;			/* we allocated one more for the trailing 0 */
			break;
	} 
	return 0;
}
	
	
sap_item_t *sap_parse(const char *buffer, int length) {
	const char *curLine = buffer;
	const char *nextLine;
	sap_item_t *first = NULL;
	
	while (length>0) {
		for (nextLine = curLine; *nextLine!='\n'; nextLine++) {
			/* We've reached the end before finding a new line break -- return what we have */
			if (--length == 0) return first;
		}
		
		sap_handle_line(&first, curLine, nextLine);
		/* skip the \n for next time */
		curLine=nextLine+1;
		length--;
	} 
	
	return first;
}


void sap_destroy(sap_item_t *item) {
	while (item!=NULL) {
		sap_item_t *aux;
	
		if (item->fieldname!=NULL) free(item->fieldname);
		if (item->value!=NULL) {
			if (item->type == SAP_STRUCT) {
				int i;
				sap_item_t **children = (sap_item_t **) item->value;
				for (i=0;i<item->num_elements;i++) {
					sap_destroy(children[i]);
				}
			} 
			free(item->value);
		}
		aux = item;
		item = item->next;
		
		free(aux);
	}
}

void sap_print_items(int indent, sap_item_t *item) {
	while (item!=NULL) {
		int i;
		for (i=0;i<item->num_elements;i++) {
			printf("%*s%s",indent, "", item->fieldname);
			if (item->is_array) {
				printf("[%i]",i);
			}
			switch(item->type) {
				case SAP_STRUCT:
					printf("\n");
					sap_print_items(indent + 4, ((sap_item_t **) item->value)[i]);
					break;
				case SAP_DOUBLE:
					printf(" = %f\n", ((double *)item->value)[i]);
					break;
				case SAP_LONG:
					printf(" = %li\n", ((long *)item->value)[i]);
					break;
				case SAP_TEXT:
					printf(" = \"%s\"\n", (char *)item->value);
					break;
			}
		}
		item = item->next;
	}

}

void sap_print(sap_item_t *list) {
	printf("----------------------------------------------\r\n");
	sap_print_items(0, list);
	printf("----------------------------------------------\r\n");
}

