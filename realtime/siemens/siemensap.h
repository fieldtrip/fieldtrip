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
extern "C" {
#endif

/** Possible types we can deal with. 
*/
typedef enum { SAP_LONG, SAP_DOUBLE, SAP_STRUCT, SAP_TEXT } sap_field_type_t;

/** Elements of the linked list of name/value pairs. 
	Values can also be arrays, sub-structs, or arrays of sub-structs 
*/
typedef struct sap_item {
	char *fieldname;			/* Name of the field, 0-terminated */
	void *value;				/* Pointer to value(s), or pointer to pointer(s) in case of SAP_STRUCT */
	sap_field_type_t type;		/* see above */
	int is_array;               /* flag */
	int num_elements;			/* must be 1 for scalar, >=0 for arrays */
	
	struct sap_item *next;	/* next element in list, NULL if last */
} sap_item_t;

/** This function allocates a new list item including space for the fieldname and the value.
	@param len_name     Length of the field name, must not be negative
	@param name  		Fieldname, must not be NULL, does not need to be 0-terminated
	@param size_value 	Size to be allocated for the value field, or 0 if no allocation should be done.
	@return Pointer to the new list item, or NULL if out of memory.
*/
sap_item_t *sap_alloc_field(int len_name, const char *name, int size_value);

/** This function searches for a given fieldname at the level of items pointed to by "first".
	Note that fieldnames are case-sensitive, and that you cannot use this directly 
	to retrieve nested fields like 'structA.fieldB'
	@param first	Pointer to first element of linked list, may be NULL
	@param len_name Length of fieldname, must not be negative
	@param name		Fieldname, must not be NULL, does not need to be 0-terminated
	@return 	Pointer to corresponding list item, or NULL in case no matching item has been found.
*/
sap_item_t *sap_search_field(sap_item_t *first, int len_name, const char *name);

/** This function is used internally for parsing a single line and adding the results to the
	linked list pointed to by 'first'. It will call itself recursively to deal with substructures.
	@param first	First element of linked list
	@param line		Start of line to be parsed
	@param line_end	End of line to be parsed
	@return		1 	in case the line could be succesfully parsed
				0	in case of errors (but also for empty lines and comments)
	The return value and its semantics might change in future implementations.
*/
int sap_handle_line(sap_item_t **first, const char *line, const char *line_end);

/** This functions is used internally for printing a substructure at a given indentation level.
	@param indent	Indentation level (in number of spaces), must not be negative
	@param item		Item to be printed
	If the 'item' contains a substructure, this will be printed recursively at
	an increased indentation.
*/
void sap_print_items(int indent, sap_item_t *item);

/** This function parses the string representation pointed to by 'buffer' and creates
	a nested linked list of (fieldname, value) pairs.
	@param buffer	String representation to be parsed
	@param length	Length of 'buffer'
	@return 	Pointer to first element of linked list, or NULL if nothing could be parsed
*/
sap_item_t *sap_parse(const char *buffer, int length);

/** Recursively reverses the list in place, and changes 'first' pointer accordingly.
	This function does not re-allocate any memory but only changes the pointers. Pointers
	to items retrieved previously (e.g., using sap_search_field) will still be valid,
	but the 'next' element will be different.
	@param first	Pointer to pointer to first list item
*/
void sap_reverse_in_place(sap_item_t **first);

/** Recursively destroys the linked list pointed to by 'item' and frees memory.
	@param item 	First element of list to be destroyed (can also be NULL)
	You cannot use the list or any of its elements after this call.
*/
void sap_destroy(sap_item_t *item);

/** Recursively print the linked list including indentation for substructures. Mainly for
	inspection and debugging.
	@param	list	First element of list (can also be NULL).
*/
void sap_print(sap_item_t *list);

#ifdef __cplusplus
}
#endif

#endif
