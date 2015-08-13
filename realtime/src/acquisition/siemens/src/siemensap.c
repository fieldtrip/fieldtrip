/*
 * Copyright (C) 2010, Stefan Klanke
 * 	Modified by Tim van Mourik 2015
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 *
 * For API documentation see siemensap.h
 */

#include <siemensap.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

sap_item_t *sap_alloc_field(int len_name, const char *name, int size_value) {
    sap_item_t *F = (sap_item_t *) calloc(1, sizeof(sap_item_t));
    /*sap_item_t *F = (sap_item_t *) malloc(sizeof(sap_item_t));*/
    if (F==NULL) return NULL;
    
#ifdef _DEBUG_
    printf("-allocating- %.*s (%i)\n", len_name, name, size_value);
#endif
    
    if ((F->fieldname = (char *) calloc(len_name+1, 1)) == NULL) {
        /*if ((F->fieldname = (char *) malloc(len_name+1)) == NULL) {*/
        free(F);
        return NULL;
    }
    memcpy(F->fieldname, name, len_name);
    F->fieldname[len_name]=0;
    if (size_value == 0) {
        F->num_elements = 0;
        F->value = NULL;
    } else {
        F->num_elements = 1;
        F->value = calloc(size_value,1);
        /*F->value = malloc(size_value);*/
        if (F->value == NULL) {
            free(F->fieldname);
            free(F);
            return NULL;
        }
    }
    return F;
}


const sap_item_t *sap_search_field(const sap_item_t *first, int len_name, const char *name) {
    while (first!=NULL) {
        if (0==strncmp(first->fieldname, name, len_name)) {
            /* the first len_name characters are the same
             * -> this means that fieldname is at least len_name characters long
             * -> return item if fieldname is exactly len_name characters long
             */
            if (first->fieldname[len_name] == 0) return first;
        }
        first = first->next;
    }
    return NULL;
}

/* parse one line and add to the list pointed to by "*first"
 * returns 0 on errors, empty lines, or comments
 * returns 1 if line was parsed and added to the list
 */
int sap_handle_line(sap_item_t **first, const char *line, const char *line_end) {
    /* looked-up or newly created item */
    sap_item_t *item;
    /* used during parsing for storing the determined type */
    sap_field_type_t typ = SAP_LONG;
    /* for detecting arrays and parsing their indicies */
    long arrayIndex = 0;
    int isArray = 0;
    int isArraySize = 0;
    /* these are for storing numeric values during parsing, before a new item is created */
    long longVal;
    double doubleVal;
    /* these are for storing noteworthy positions inside the line */
    const char *val_start, *val_end, *aux;
    const char *name, *name_end;
    const char *nn;
    /* length of field name and size of current value type (LONG / DOUBLE / VOID * / TEXT) */
    int len_name, size_value;
    
#ifdef _DEBUG_
    printf("**********\n");
    printf("%.*s\n", line_end-line, line);
    if (*first != NULL)
    {
        printf("[%s] Type %i  IsArray: %i  NumElements: %i\n",
                (*first)->fieldname, (*first)->type, (*first)->is_array, (*first)->num_elements);
    }
    else { printf("first = NULL"); }
    printf("**********]\n");
#endif
    
    name = line;
    /* skip whitespace at the beginning */
    while (*name == ' ' || *name == '\t') name++;
    /* check for comments and empty lines */
    if (name[0]=='#' || name[0]=='\r' || name[0]=='\n') return 0;
    
    /* detect end of field name */
    for (name_end = name; name_end != line_end; name_end++) {
        if (*name_end == '[') {
            isArray = 1;
            break;
        }
        if (*name_end == '.') {
            if (!strncmp(name_end+1, "__attribute__.size", 18))
            {
                len_name = (int) (name_end - name);
                nn = name_end + 19;
#ifdef _DEBUG_
                printf("ARRAY SIZE: (%.*s) %.*s\n", name_end - name, name, line_end - nn, nn);
#endif
                isArraySize = 1;
                name_end += 19;
            }
            else
            {
                typ = SAP_STRUCT;
                break;
            }
        }
        /* otherwise, name is finished at the first space, tab, or = sign */
        if (*name_end == ' ') break;
        if (*name_end == '\t') break;
        if (*name_end == '=') break;
        if (*name_end == 0) {
            fprintf(stderr,"Something went wrong in sap_parse -- this should not happen.\n");
            return -1;
        }
    }
    
    /* length of the fieldname is given by difference between pointers */
    if (!isArraySize)
    {
        len_name = (int) (name_end - name);
    }
    
#ifdef _DEBUG_
    printf("*NAME* = '%.*s'\n", name_end - name, name);
#endif
    
    if (isArray) {
        /* try to parse the index */
        char *end_index;
        arrayIndex = strtol(name_end+1, (char **) &end_index, 10);
        if (arrayIndex < 0 || end_index[0]!=']' ) { /* this includes LONG_MIN in case of parse error */
            fprintf(stderr, "Invalid array index description.\n");
            return 0;
        }
        /* we can also have an array of structs */
        if (end_index[1]=='.') {
            typ = SAP_STRUCT;
            name_end = end_index + 1;
        }
    } else {
        arrayIndex = 0;
    }
    
    
    for (aux = name_end; aux != line_end; aux++) {
        if (*aux == '=') break;
    }
    if (aux == line_end) {
        fprintf(stderr, "Invalid definition: No = sign\n");
        return 0;
    }
    val_start = aux+1; /* now points just behind the = sign */
    
    /* If it's not a struct, try to parse the value */
    if (isArraySize)
    {
        longVal = strtol(val_start, (char **) &val_end, 10);
        if (longVal == LONG_MIN) {
            fprintf(stderr, "Could not parse array size.\n");
            return 0;
        }
        
        /* New feature in Siemens Skyra protocols: array length is given beforehand
         We use this information to create an otherwise empty sap_item and leave for now */
        item = sap_alloc_field(len_name, name, 0);
#ifdef _DEBUG_
        printf("ARRAY SIZE: %s   [%i]\n", item->fieldname, arrayIndex);
#endif
        item->next = *first;
        item->num_elements = longVal;
        item->is_array = 1;
        item->type = SAP_ARRAY_SIZE;
        *first = item;
        return 1;
    }
    else if (typ != SAP_STRUCT)
    {
        /* Investigate right hand side to see whether it's a string,
         * a hexadecimal number (->long), or contains a . (-> double) */
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
            size_value = (int) (val_end - val_start + 1); /* +1 for trailing 0 */
        } else if (aux[0]=='x' && aux[-1]=='0') {
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
        size_value = sizeof(void *); /* values are pointers for SAP_STRUCT */
    }
    
#ifdef _DEBUG_
    printf("Still here. Typ=%i, size=%i longVal=%i\n", typ, size_value, longVal);
#endif
    
    /* the following cast is for removing the constness warning */
    item = (sap_item_t *) sap_search_field(*first, len_name, name);
#ifdef _DEBUG_
    printf("Still here. Typ=%i, size=%i longVal=%i\n", typ, size_value, longVal);
#endif
    if (item == NULL) {
        /* create new field */
        item = sap_alloc_field(len_name, name, size_value);
        if (item == NULL) {
            fprintf(stderr, "Could not allocate sap_item: Out of memory\n");
            return 0;
        }
        item->next = *first;
        item->num_elements = 1;
        item->is_array = isArray;
        item->type = typ;
        *first = item;
#ifdef _DEBUG_
        printf("New field: %.*s  %i, %i\n",len_name,name, isArray, typ);
#endif
    }
    if (isArray)
    {
        if (item->type == SAP_ARRAY_SIZE)
        {
            item->value = calloc(item->num_elements, size_value);
            if (item->value == NULL) {
                /* We can't allocate the buffer */
                fprintf(stderr, "Could not allocate memory for value array: Out of memory\n");
                return 0;
            }
            item->type = typ;
        }
        else if (arrayIndex >= item->num_elements)
        {
            int newSize = (arrayIndex+1)*size_value;
            int oldSize = item->num_elements*size_value;
            /* 	need to reallocate */
            char *newBuf = (char *) realloc(item->value, newSize);
            
            if (newBuf == NULL) {
                /* No harm done so far, we just can't add this value to the buffer */
                fprintf(stderr, "Could not allocate memory for value array: Out of memory\n");
                return 0;
            }
            /* zero-out the new memory */
            memset(newBuf + oldSize, 0, newSize - oldSize);
            item->value = (void *) newBuf;
            item->num_elements = arrayIndex+1;
        }
    }
    
    
    switch(typ) {
        case SAP_STRUCT:
            /* If there is an error further back in this line, this call will return 0, so we do as well
             * Regarding clean-up, we could think of removing this item, e.g., if it was newly created
             * But the default value will be filled with zeros, corresponding to an empty structure,
             * so this won't cause problems.
             */
            return sap_handle_line(((sap_item_t **) item->value)+arrayIndex, name_end +1, line_end);
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
    return 1;
}


void sap_reverse_in_place(sap_item_t **first) {
    sap_item_t *prev = NULL, *next;
    sap_item_t *item;
    
    if (first == NULL || *first == NULL) return;
    
    item = *first;
    while (1) {
        if (item->type == SAP_STRUCT) {
            sap_item_t **children = (sap_item_t **) item->value;
            int i;
            /* the array itself will of course not be reversed */
            for (i=0;i<item->num_elements;i++) {
                sap_reverse_in_place(children + i);
            }
        }
        
        next = item->next;
        item->next = prev;
        prev = item;
        if (next == NULL) break;
        item = next;
    }
    *first = item;
}

sap_item_t *sap_parse(const char *buffer, unsigned int length) {
    const char *curLine = buffer;
    const char *nextLine;
    sap_item_t *first = NULL;
    
    if (buffer==NULL || length==0) return NULL;
    
    while (length>0) {
        for (nextLine = curLine; *nextLine!='\n'; nextLine++) {
            /* We've reached the end before finding a new line break -- return what we have */
            if (--length == 0) return first;
        }
        
        sap_handle_line(&first, curLine, nextLine);
        /*
#ifdef _DEBUG_
      if (first != NULL)
      {
         printf("[%s] Type %i  IsArray: %i  NumElements: %i\n",
               first->fieldname, first->type, first->is_array, first->num_elements);
      }
#endif
         */
        /* skip the \n for next time */
        curLine=nextLine+1;
        length--;
    }
    
    sap_remove_empty_arrays(&first);
    
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

void sap_print_items(int indent, const sap_item_t *item) {
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

void sap_print(const sap_item_t *list) {
    sap_print_items(0, list);
}


const sap_item_t *sap_search_deep(const sap_item_t *list, const char *fieldname) {
    const char *dotPos, *brkPos;
    int len, dotLen, brkLen, slen, index;
    const sap_item_t *item;
    
    if (fieldname == NULL) return NULL;
    
    len = (int) strlen(fieldname);
    
    dotPos = strchr(fieldname, '.');
    dotLen = (dotPos==NULL) ? len : (int) (dotPos - fieldname);
    
    brkPos = strchr(fieldname, '[');
    brkLen = (brkPos==NULL) ? len : (int) (brkPos - fieldname);
    
    slen = (dotLen<brkLen) ? dotLen : brkLen;
    item = sap_search_field(list, slen, fieldname);
    
    if (item == NULL) return NULL; 	/* not found */
    if (dotLen == len) return item; /* user did not search for a sub-struct */
    
    if (brkLen < len) {
        /* user is looking for sub-struct of an array */
        long val = strtol(brkPos+1, NULL, 10);
        if (val < 0) return NULL;
        index = (int) val;
    } else {
        index = 0;
    }
    
    item = ((sap_item_t **) item->value)[index];
    
    return sap_search_deep(item, dotPos+1);
}


int sap_get_essentials(const sap_item_t *list, sap_essentials_t *E) {
    const sap_item_t *item;
    int numFound = 0;
    
    if (list == NULL || E == NULL) return -1;
    
    memset(E, 0, sizeof(sap_essentials_t));
    
    item = sap_search_deep(list, "alTR");
    if (item!=NULL && item->type == SAP_LONG) {
        E->TR = ((long *) item->value)[0];
        numFound++;
    }
    
    item = sap_search_deep(list, "lContrasts");
    if (item!=NULL && item->type == SAP_LONG) {
        long cs = *((long *) item->value);
        if (cs>=0) {
            E->numberOfContrasts = (unsigned int) cs;
            numFound++;
        }
    }
    
    item = sap_search_deep(list, "sKSpace.lBaseResolution");
    if (item!=NULL && item->type == SAP_LONG) {
        long res = *((long *) item->value);
        if (res>0) {
            E->readoutPixels = res;
            numFound++;
        }
    }
    
    item = sap_search_deep(list, "sSliceArray.lSize");
    if (item!=NULL && item->type == SAP_LONG) {
        long slices = *((long *) item->value);
        if (slices > 0) {
            E->numberOfSlices = slices;
            numFound++;
        }
    }
    
    item = sap_search_deep(list, "sSliceArray.asSlice[0].dPhaseFOV");
    if (item!=NULL && item->type == SAP_DOUBLE) {
        E->phaseFOV = *((double *) item->value);
        numFound++;
    }
    
    item = sap_search_deep(list, "sSliceArray.asSlice[0].dReadoutFOV");
    if (item!=NULL && item->type == SAP_DOUBLE) {
        E->readoutFOV = *((double *) item->value);
        numFound++;
    }
    
    if (E->phaseFOV > 0.0 && E->readoutFOV > 0.0) {
        E->phasePixels = (unsigned int) (E->readoutPixels * E->phaseFOV / E->readoutFOV + 0.5);
        numFound++;
    }
    
    item = sap_search_deep(list, "sSliceArray.asSlice[0].dThickness");
    if (item!=NULL && item->type == SAP_DOUBLE) {
        E->sliceThickness = *((double *) item->value);
        numFound++;
    }
    
    return numFound;
}

void sap_remove_empty_arrays(sap_item_t **list)
{
    sap_item_t *item = *list;
    
    if (item == NULL) return;
    
    if (item->type == SAP_STRUCT)
    {
        int i;
        sap_item_t **children = (sap_item_t **) item->value;
        for (i=0;i<item->num_elements;i++)
        {
            sap_remove_empty_arrays(&children[i]);
        }
    }
    if (item->type == SAP_ARRAY_SIZE)
    {
        *list = item->next;
        if (item->value != NULL) free(item->value);
        if (item->fieldname!=NULL) free(item->fieldname);
        free(item);
        
        sap_remove_empty_arrays(list);
    }
    else
    {
        sap_remove_empty_arrays(&item->next);
    }
}
