/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */


/* small library for parsing Siemens Protocol information as produced by their pMrProt->fwrite command */
 
#ifndef _SIEMENSAP_H
#define _SIEMENSAP_H

#ifdef __cplusplus
#extern "C" {
#endif

typedef enum { SAP_LONG, SAP_DOUBLE, SAP_STRUCT, SAP_TEXT } sap_field_type_t;

typedef struct sap_item {
	char *fieldname;			/* Name of the field, 0-terminated */
	void *value;				/* Pointer to value(s), or pointer to pointer(s) in case of SAP_STRUCT */
	sap_field_type_t type;		/* see above */
	int is_array;               /* flag */
	int num_elements;			/* must be 1 for scalar, >=0 for arrays */
	
	struct sap_item *next;	/* next element in list, NULL if last */
} sap_item_t;

sap_item_t *sap_alloc_field(int len_name, const char *name, int size_value);
sap_item_t *sap_search_field(sap_item_t *first, int len_name, const char *name);
int sap_handle_line(sap_item_t **first, const char *line, const char *line_end);
void sap_print_items(int indent, sap_item_t *item);

sap_item_t *sap_parse(const char *buffer, int length);
void sap_destroy(sap_item_t *item);
void sap_print(sap_item_t *list);

#ifdef __cplusplus
}
#endif

#endif
