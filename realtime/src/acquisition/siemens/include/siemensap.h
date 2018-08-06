/*
 * This is a small library for parsing Siemens Protocol information as produced by their pMrProt->fwrite command 
 * 
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
*/
 
#ifndef _SIEMENSAP_H
#define _SIEMENSAP_H

#ifdef __cplusplus
extern "C" {
#endif

/** Possible types we can deal with. 
*/
typedef enum { 
	SAP_LONG, 			/**< Integer values */
	SAP_DOUBLE, 		/**< Real numbers 	*/
	SAP_STRUCT, 		/**< Substructures  */
	SAP_TEXT, 			/**< Strings		*/
   SAP_ARRAY_SIZE = 256   /**< Only temporarily use (Skyra) */
} sap_field_type_t;

/** Elements of the linked list of name/value pairs. 
	Values can also be arrays, sub-structs, or arrays of sub-structs 
*/
typedef struct sap_item {
	char *fieldname;			/**< Name of the field, 0-terminated */
	void *value;				/**< Pointer to value(s), or pointer to pointer(s) in case of SAP_STRUCT */
	sap_field_type_t type;		/**< Type of this field */
	int is_array;               /**< Flag that determines whether this field is an array */
	int num_elements;			/**< This must be 1 for scalar, >=0 for arrays */
	
	struct sap_item *next;		/**< next element in list, NULL if last */
} sap_item_t;

/** Essential information about the MR protocol provided for convenience.
*/
typedef struct sap_essentials {
	long TR;							/**< TR[0] */
	unsigned int readoutPixels;			/**< Image size (pixels) in readout (X) direction */
	unsigned int phasePixels;			/**< Image size (pixels) in phase (Y) direction */
	unsigned int numberOfSlices;		/**< Number of slices */
	unsigned int numberOfContrasts;		/**< Number of contrasts */
	double readoutFOV;					/**< Image size (mm) in readout (X) direction */
	double phaseFOV;					/**< Image size (mm) in phase (Y) direction */
	double sliceThickness;				/**< Slice thickness (mm) */
} sap_essentials_t;

/** This constant must be equal to the number of fields in sap_essentials_t */
#define SAP_NUM_ESSENTIALS        8

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
const sap_item_t *sap_search_field(const sap_item_t *first, int len_name, const char *name);

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
void sap_print_items(int indent, const sap_item_t *item);

/** This function parses the string representation pointed to by 'buffer' and creates
	a nested linked list of (fieldname, value) pairs.
	@param buffer	String representation to be parsed
	@param length	Length of 'buffer'
	@return 	Pointer to first element of linked list, or NULL if nothing could be parsed
*/
sap_item_t *sap_parse(const char *buffer, unsigned int length);

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
void sap_print(const sap_item_t *list);

/** Recursively search for a field such as "foo[13].bar[3].element". Array indices
    are only meaningful here if another subfield follows, that is, searching for
	"foo[12]", "foo[5]", and "foo" will return the same item.
	@param  list		Linked list of fieldname/value items
	@param  fieldname  	0-terminated string describing the desired field (as found in Siemens .pro files)
	@return 	The corresponding item within 'list', or NULL if not found
*/
const sap_item_t *sap_search_deep(const sap_item_t *list, const char *fieldname);

/** Try to retrieve "essential" protocol information from a given linked list.
	@param list		Linked list of protocol definition key/value pairs
	@param E		Pointer to essential information data structure
	@return 	The number of detected parameters, or -1 if one of the parameters is NULL
*/
int sap_get_essentials(const sap_item_t *list, sap_essentials_t *E);

void sap_remove_empty_arrays(sap_item_t **list);

#ifdef __cplusplus
}
#endif

#endif
