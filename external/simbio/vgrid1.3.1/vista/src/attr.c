/*
 *  Copyright 1993, 1994 University of British Columbia
 *
 *  Permission to use, copy, modify, distribute, and sell this software and its
 *  documentation for any purpose is hereby granted without fee, provided that
 *  the above copyright notice appears in all copies and that both that
 *  copyright notice and this permission notice appear in supporting
 *  documentation. UBC makes no representations about the suitability of this
 *  software for any purpose. It is provided "as is" without express or
 *  implied warranty.
 *
 *  Author: Arthur Pope, UBC Laboratory for Computational Intelligence
 */

/*
 *  Remodeled for SimBio by F. Kruggel (kruggel@cns.mpg.de) - 15/07/00
 */

#include <vista.h>

static VStringConst Encode (VDictEntry *dict, VRepnKind repn, va_list *args)
{
    VLong i_value = 0;
    VDouble f_value = 0.0;
    VString s_value = NULL;
    static char buf[40];

    /* Fetch the attribute value: */
    switch (repn) {
    case VBitRepn: i_value = va_arg (*args, VBitPromoted); break;
    case VUByteRepn: i_value = va_arg (*args, VUBytePromoted); break;
    case VSByteRepn: i_value = va_arg (*args, VSBytePromoted); break;
    case VShortRepn: i_value = va_arg (*args, VShortPromoted); break;
    case VLongRepn: i_value = va_arg (*args, VLongPromoted); break;
    case VFloatRepn: f_value = va_arg (*args, VFloatPromoted); break;
    case VDoubleRepn: f_value = va_arg (*args, VDoublePromoted); break;
    case VBooleanRepn: i_value = va_arg (*args, VBooleanPromoted); break;
    case VStringRepn: s_value = va_arg (*args, VString); break;

    default:
	VError ("VEncodeAttrValue: Can't encode from %s", VRepnName (repn));
    }

    /* If its numeric, convert it to a string: */
    switch (repn) {

    case VBitRepn:
    case VUByteRepn:
    case VSByteRepn:
    case VShortRepn:
    case VLongRepn:
    case VBooleanRepn:
	sprintf (s_value = buf, "%ld", (long) i_value);
	break;

    case VFloatRepn:
    case VDoubleRepn:
	sprintf (s_value = buf, "%.20g", (double) f_value);
	break;

    default:
	break;
    }

    /* If a dictionary was supplied, try to map the encoded value to
       a keyword: */
    if (dict)
	switch (repn) {

	case VBitRepn:
	case VUByteRepn:
	case VSByteRepn:
	case VShortRepn:
	case VLongRepn:
	case VBooleanRepn:
	    dict = VLookupDictValue (dict, VLongRepn, i_value);
	    break;

	case VFloatRepn:
	case VDoubleRepn:
	    dict = VLookupDictValue (dict, VDoubleRepn, f_value);
	    break;

	case VStringRepn:
	    dict = VLookupDictValue (dict, VStringRepn, s_value);
	    break;

	default:
	    break;
	}
    return dict ? dict->keyword : s_value;
}

static VAttrRec *NewAttr (VStringConst name, VDictEntry *dict,
			  VRepnKind repn, va_list *args)
{
    size_t new_value_size, name_size;
    VPointer value;
    VAttrRec *a;

    name_size = strlen (name);
    switch (repn) {

    case VBitRepn:
    case VUByteRepn:
    case VSByteRepn:
    case VShortRepn:
    case VLongRepn:
    case VFloatRepn:
    case VDoubleRepn:
    case VBooleanRepn:
    case VStringRepn:

	/* Determine the amount of storage needed to record the new value.
	   In some cases, this requires first encoding the new value as a
	   string. */
	if (repn == VStringRepn && ! dict)
	    value = (VPointer) va_arg (*args, VStringConst);
	else value = (VPointer) Encode (dict, repn, args);
	new_value_size = strlen (value) + 1;

	/* Allocate storage for the new attribute and copy in its value: */
	a = VMalloc (sizeof (VAttrRec) + name_size + new_value_size);
	a->repn = VStringRepn;
	a->value = (a->name + name_size + 1);
	strcpy (a->value, value);
	break;

    default:
	a = VMalloc (sizeof (VAttrRec) + name_size);
	a->repn = repn;
	a->value = va_arg (*args, VPointer);
    }
    strcpy (a->name, name);
    return a;
}

static void SetAttr (VAttrListPosn *posn, VDictEntry *dict,
		     VRepnKind repn, va_list *args)
{
    size_t old_value_size, new_value_size, name_size;
    VPointer value;
    VAttrRec *a = posn->ptr;

    /* Determine the amount of storage needed to record the new value. In some
       cases, this requires first encoding the new value as a string. */
    name_size = strlen (a->name);
    switch (repn) {

    case VBitRepn:
    case VUByteRepn:
    case VSByteRepn:
    case VShortRepn:
    case VLongRepn:
    case VFloatRepn:
    case VDoubleRepn:
    case VBooleanRepn:
    case VStringRepn:
	if (repn == VStringRepn && ! dict)
	    value = (VPointer) va_arg (*args, VStringConst);
	else value = (VPointer) Encode (dict, repn, args);
	new_value_size = strlen (value) + 1;
	break;

    default:
	value = va_arg (*args, VPointer);
	new_value_size = 0;
    }

    /* Is enough storage allocated for it in the existing attribute node? */
    switch (a->repn) {

    case VStringRepn:
	old_value_size = strlen (a->value) + 1;
	break;

    default:
	old_value_size = 0;
    }
    if (old_value_size < new_value_size) {

	/* It exists, but it's too small: */
	a = VMalloc (sizeof (VAttrRec) + name_size + new_value_size);
	a->next = posn->ptr->next;
	a->prev = posn->ptr->prev;
	if (a->next)
	    a->next->prev = a;
	else posn->list->prev = a;
	if (a->prev)
	    a->prev->next = a;
	else posn->list->next = a;
	strcpy (a->name, posn->ptr->name);
	VFree (posn->ptr);
	posn->ptr = a;
    }

    /* Copy in the attribute's new value: */
    switch (repn) {

    case VBitRepn:
    case VUByteRepn:
    case VSByteRepn:
    case VShortRepn:
    case VLongRepn:
    case VFloatRepn:
    case VDoubleRepn:
    case VBooleanRepn:
    case VStringRepn:
	a->repn = VStringRepn;
	a->value = a->name + name_size + 1;
	strcpy (a->value, value);
	break;

    default:
	a->repn  = repn;
	a->value = value;
    }
}

static void FreeAttrValue (VStringConst routine, VAttrRec *a)
{
    VTypeMethods *methods;

    switch (a->repn) {

    case VAttrListRepn:
  	VDestroyAttrList (a->value);
	break;

    case VBundleRepn:
	VDestroyBundle (a->value);
	break;

    case VPointerRepn:
    case VStringRepn:
	break;

    default:
	methods = VRepnMethods(a->repn);
	if (! methods) {
          VError ("%s: %s attribute has invalid repn %d",
		          routine, a->name, a->repn);
        }
        else {
	  (methods->destroy) (a->value);
        }
    }
}

void VAppendAttr (VAttrList list, VStringConst name,
		  VDictEntry *dict, VRepnKind repn, ...)
{
    va_list args;
    VAttrRec *a;

    /* Create the new attribute node: */
    va_start (args, repn);
    a = NewAttr (name, dict, repn, & args);
    va_end (args);

    /* Append it: */
    a->next = NULL; a->prev = list->prev;
    if (a->prev) a->prev->next = a;
    else list->next = a;
    list->prev = a;
}

VAttrList VCopyAttrList (VAttrList list)
{
    VAttrList new_list = VCreateAttrList ();
    size_t name_size, value_size;
    VAttrRec *old_a, *new_a;
    VBundle old_b, new_b;
    VTypeMethods *methods;

    /* For each node of the old list: */
    for (old_a = list->next; old_a; old_a = old_a->next) {

	/* Compute the amount of storage needed for a copy of the node: */
	name_size = strlen (old_a->name);
	value_size = (old_a->repn == VStringRepn) ?
	    strlen ((VStringConst) old_a->value) + 1 : 0;

	/* Allocate that size and fill in the node's value: */
	new_a = VMalloc (sizeof (VAttrRec) + name_size + value_size);
	strcpy (new_a->name, old_a->name);
	switch (new_a->repn = old_a->repn) {

	case VAttrListRepn:
	    new_a->value = VCopyAttrList (old_a->value);
	    break;

	case VBundleRepn:
	    old_b = old_a->value;
	    new_b = VCreateBundle (old_b->type_name,
				   VCopyAttrList (old_b->list),
				   old_b->length, NULL);
	    if (old_b->length > 0) {
		new_b->data = VMalloc (old_b->length);
		memcpy (new_b->data, old_b->data, old_b->length);
	    }
	    new_a->value = new_b;
	    break;

	case VPointerRepn:
	    new_a->value = old_a->value;
	    break;

	case VStringRepn:
	    new_a->value = (VPointer) (new_a->name + name_size + 1);
	    strcpy (new_a->value, old_a->value);
	    break;

	default:
	    methods = VRepnMethods (new_a->repn);
	    if (methods) new_a->value = (methods->copy) (old_a->value);
	    else VError ("VCopyAttrList: %s attribute has invalid repn %d",
			 old_a->name, old_a->repn);
	}

	/* Append it to the new list: */
	new_a->next = NULL;
	new_a->prev = new_list->prev;
	if (new_a->prev) new_a->prev->next = new_a;
	if (! new_list->next) new_list->next = new_a;
	new_list->prev = new_a;
    }
    return new_list;
}

VAttrList VCreateAttrList (void)
{
    VAttrList list;

    list = VNew (VAttrRec);
    list->next = list->prev = list->value = NULL;
    list->repn = VUnknownRepn;	/* not mistakable for an attribute */
    list->name[0] = 0;
    return list;
}

VBundle VCreateBundle (VStringConst type_name, VAttrList list,
		       size_t length, VPointer data)
{
    VBundle b;

    b = VMalloc (sizeof (VBundleRec) + strlen (type_name));
    strcpy (b->type_name, type_name);
    b->list = list;
    b->length = length;
    b->data = data;
    return b;
}

VBoolean VDecodeAttrValue (VStringConst str, VDictEntry *dict,
			   VRepnKind repn, VPointer value)
{
    VLong i_value = 0;
    VDouble f_value = 0.0;
    char *cp = NULL, buf[20];

    /* If a dict is provided, see if str maps to any dict entry keyword,
       substituting the associated value if found: */
    if (dict) {
	dict = VLookupDictKeyword (dict, str);

	/* If there's a dictionary entry, complete it: */
	if (dict && ! dict->svalue) {
	    str = NULL;
	    dict->icached = dict->fcached = TRUE;
	    sprintf (buf, "%ld", (long) dict->ivalue);
	    dict->svalue = VNewString (buf);
	    dict->fvalue = dict->ivalue;
	}
    }

    /* Convert to the internal representation: */
    switch (repn) {

    case VBitRepn:
    case VUByteRepn:
    case VSByteRepn:
    case VShortRepn:
    case VLongRepn:
    case VBooleanRepn:
	if (dict) {
	    if (dict->icached)
		i_value = dict->ivalue;
	    else {
		dict->ivalue = i_value = strtol (dict->svalue, & cp, 0);
		dict->icached = TRUE;
	    }
	} else i_value = strtol (str, & cp, 0);
	break;

    case VFloatRepn:
    case VDoubleRepn:
	if (dict) {
	    if (dict->fcached)
		f_value = dict->fvalue;
	    else {
		dict->fvalue = f_value = strtod (dict->svalue, & cp);
		dict->fcached = TRUE;
	    }
	} else f_value = strtod (str, & cp);
	break;

    case VStringRepn:
	if (dict)
	    str = dict->svalue;
	break;

    default:
	VError ("VDecodeAttrValue: Can't decode to %s", VRepnName (repn));
    }
    if (cp && *cp)
	return FALSE;

    /* Store at *value: */
    switch (repn) {
    case VBitRepn: * (VBit *) value = i_value; break;
    case VUByteRepn: * (VUByte *) value = i_value; break;
    case VSByteRepn: * (VSByte *) value = i_value; break;
    case VShortRepn: * (VShort *) value = i_value; break;
    case VLongRepn: * (VLong *) value = i_value; break;
    case VFloatRepn: * (VFloat *) value = f_value; break;
    case VDoubleRepn: * (VDouble *) value = f_value; break;
    case VBooleanRepn: * (VBoolean *) value = i_value; break;
    case VStringRepn: * (VStringConst *) value = str; break;

    default:
	break;
    }

    return TRUE;
}

void VDeleteAttr (VAttrListPosn *posn)
{
    VAttrRec *a = posn->ptr;

    /* Remove it from the list: */
    if (a->next)
	a->next->prev = a->prev;
    if (a->prev)
	a->prev->next = a->next;
    if (posn->list->next == a)
	posn->list->next = a->next;
    if (posn->list->prev == a)
	posn->list->prev = a->prev;

    /* Make posn point to the next attribute, or nothing: */
    posn->ptr = a->next;

    VFree (a);
}

void VDestroyAttrList (VAttrList list)
{
    VAttrRec *a, *a_next;

    if (! list) {
	VWarning ("VDestroyAttrList: called with NULL list");
	return;
    }

    /* For each attribute in the list: */
    for (a = list->next; a; a = a_next) {
	a_next = a->next;

	/* Free any storage used for the attribute's value: */
	FreeAttrValue ("VDestroyAttrList", a);

	/* Free the attribute record itself: */
  	VFree (a);
    }
    VFree (list);
}

void VDestroyBundle (VBundle b)
{
    VDestroyAttrList (b->list);
    if (b->length > 0)
	VFree (b->data);
    VFree (b);
}

VStringConst VEncodeAttrValue (VDictEntry *dict, VRepnKind repn, ...)
{
    va_list args;
    VStringConst str;

    va_start (args, repn);
    str = Encode (dict, repn, & args);
    va_end (args);
    return str;
}

VBoolean VExtractAttr (VAttrList list, VStringConst name,
		       VDictEntry *dict, VRepnKind repn, VPointer value,
		       VBooleanPromoted required)
{
    VAttrListPosn posn;

    /* If the attribute is in the list... */
    if (VLookupAttr (list, name, & posn)) {

	if (value) {

	    /* Get its value: */
	    if (! VGetAttrValue (& posn, dict, repn, value)) {
		VWarning ("VExtractAttr: %s attribute has bad value", name);
		return FALSE;
	    }

	    /* Clone or hide the value if we're about to delete it: */
	    if (repn == VStringRepn)
		* (VString *) value = VNewString (* (VString *) value);
	}

	/* Remove it from the list: */
	VDeleteAttr (& posn);
	return TRUE;
    }

    /* Otherwise complain if the attribute was a required one: */
    if (required)
	VWarning ("VExtractAttr: %s attribute missing", name);
    return ! required;
}

VGetAttrResult VGetAttr (VAttrList list, VStringConst name,
			 VDictEntry *dict, VRepnKind repn, VPointer value)
{
    VAttrListPosn posn;

    /* Look up the attribute name in the list: */
    if (! VLookupAttr (list, name, & posn))
	return VAttrMissing;

    /* Get its value in the specified representation: */
    return VGetAttrValue (& posn, dict, repn, value) ?
	VAttrFound : VAttrBadValue;
}

VBoolean VGetAttrValue (VAttrListPosn *posn, VDictEntry *dict,
			VRepnKind repn, VPointer value)
{
    /* Convert it to the requested representation: */
    switch (repn) {

    case VBitRepn:
    case VUByteRepn:
    case VSByteRepn:
    case VShortRepn:
    case VLongRepn:
    case VFloatRepn:
    case VDoubleRepn:
    case VBooleanRepn:
    case VStringRepn:
	return (VGetAttrRepn (posn) == VStringRepn &&
		VDecodeAttrValue (posn->ptr->value, dict, repn, value));

    default:
	if (VGetAttrRepn (posn) != repn)
	    return FALSE;
	* (VPointer *) value = posn->ptr->value;
	return TRUE;
    }
}

void VInsertAttr (VAttrListPosn *posn, VBooleanPromoted after,
		  VStringConst name, VDictEntry *dict, VRepnKind repn, ...)
{
    va_list args;
    VAttrRec *a;

    /* Create the new attribute node: */
    va_start (args, repn);
    a = NewAttr (name, dict, repn, & args);
    va_end (args);

    /* Insert it at the specified position: */
    if (! posn->ptr) {			/* the pointer points nowhere */
	a->next = posn->list->next;
	if (a->next) a->next->prev = a;
	a->prev = 0;
	posn->list->next = a;
	if (! posn->list->prev)
	    posn->list->prev = a;
    } else if (after) {
	a->next = posn->ptr->next;
	if (a->next) a->next->prev = a;
	else posn->list->prev = a;
	a->prev = posn->ptr;
	a->prev->next = a;
	if (posn->list->prev == a->prev)
	    posn->list->prev = a;
    } else {
	a->next = posn->ptr;
	a->prev = posn->ptr->prev;
	if (a->prev) a->prev->next = a;
	else posn->list->next = a;
	a->next->prev = a;
	if (posn->list->next == a->next)
	    posn->list->next = a;
    }
}

VBoolean VLookupAttr (VAttrList list, VStringConst name, VAttrListPosn *posn)
{
    for (VFirstAttr (list, posn); VAttrExists (posn); VNextAttr (posn))
	if (strcmp (VGetAttrName (posn), name) == 0)
	    return TRUE;
    return FALSE;
}

void VPrependAttr (VAttrList list, VStringConst name,
		   VDictEntry *dict, VRepnKind repn, ...)
{
    va_list args;
    VAttrRec *a;

    /* Create the new attribute node: */
    va_start (args, repn);
    a = NewAttr (name, dict, repn, & args);
    va_end (args);

    /* Prepend it: */
    a->next = list->next;
    if (a->next) a->next->prev = a;
    else list->prev = a;
    a->prev = NULL;
    list->next = a;
}

void VSetAttr (VAttrList list, VStringConst name,
	       VDictEntry *dict, VRepnKind repn, ...)
{
    va_list args;
    VAttrListPosn posn;
    VAttrRec *a;

    /* Locate any existing attribute of the specified name: */
    va_start (args, repn);
    if (VLookupAttr (list, name, & posn))
	SetAttr (& posn, dict, repn, & args);
    else {

	/* None exists -- append a new attribute of that name: */
	a = NewAttr (name, dict, repn, & args);
	a->next = NULL;
	a->prev = list->prev;
	if (a->prev) a->prev->next = a;
	else list->next = a;
	list->prev = a;
    }
    va_end (args);
}

void VSetAttrValue (VAttrListPosn *posn, VDictEntry *dict, VRepnKind repn, ...)
{
    va_list args;

    /* Locate any existing attribute of the specified name: */
    va_start (args, repn);
    SetAttr (posn, dict, repn, & args);
    va_end (args);
}

static VNodePtrType MakeNode (VPointer item, VNodePtrType prev,
			     VNodePtrType next)
{
    VNodePtrType result = VMalloc (sizeof (struct V_Node));

    result->item = item;
    result->prev = prev;
    result->next = next;

    return result;
}

VList VListCreate (void)
{
    VList vlist = VMalloc (sizeof (struct V_List));
    VNodePtrType dummy_head, dummy_tail;

    dummy_head = VMalloc (sizeof (struct V_Node));

    dummy_tail = VMalloc (sizeof (struct V_Node));

    dummy_head->item = NULL;
    dummy_head->prev = NULL;
    dummy_head->next = dummy_tail;

    dummy_tail->item = NULL;
    dummy_tail->prev = dummy_head;
    dummy_tail->next = NULL;

    vlist->head	   = dummy_head;
    vlist->tail	   = dummy_tail;
    vlist->current = dummy_head;
    vlist->count   = 0;

    return vlist;
}

VPointer VListFirst (VList vlist)
{
    if ( vlist->count == 0 )   /* empty vist, move beyond beginning */
	vlist->current = vlist->head;
    else		       /* vlist not empty, move to beginning */
	vlist->current = vlist->head->next;

    return vlist->current->item;
}

VPointer VListLast (VList vlist)
{
    if ( vlist->count == 0 )   /* empty vlist, move beyond end */
	vlist->current = vlist->tail;
    else		       /* vlist not empty, move to end */
	vlist->current = vlist->tail->prev;

    return vlist->current->item;
}

VPointer VListNext (VList vlist)
{
    if ( vlist->current == vlist->tail )
	/* already beyond end, no action */
	;
    else   /* move to next node */
	vlist->current = vlist->current->next;

    return vlist->current->item;
}

VPointer VListPrev (VList vlist)
{
    if ( vlist->current == vlist->head )
	/* already before beginning, no action */
	;
    else   /* move to previous node */
	vlist->current = vlist->current->prev;

    return vlist->current->item;
}

void VListAdd (VList vlist, VPointer item)
{
    VNodePtrType add_me;

    if ( vlist->current == vlist->tail )
	/* current pointer beyond end, add to end */
	vlist->current = vlist->tail->prev;

    add_me = MakeNode (item, vlist->current, vlist->current->next);

    add_me->prev->next = add_me;
    add_me->next->prev = add_me;

    vlist->current = add_me;
    vlist->count++;
}

void VListInsert (VList vlist, VPointer item)
{
    VNodePtrType add_me;

    if ( vlist->current == vlist->head )
	/* current pointer before beginning, add to beginning */
	vlist->current = vlist->head->next;

    add_me = MakeNode (item, vlist->current->prev, vlist->current);

    add_me->prev->next = add_me;
    add_me->next->prev = add_me;

    vlist->current = add_me;
    vlist->count++;
}

void VListAppend (VList vlist, VPointer item)
{
    vlist->current = vlist->tail;   /* move beyond end */
    VListAdd (vlist, item);
}

void VListPrepend (VList vlist, VPointer item)
{
    vlist->current = vlist->head;   /* move before beginning */
    VListAdd (vlist, item);
}

VPointer VListRemove (VList vlist)
{
    VPointer return_me;
    VNodePtrType free_me;

    return_me = vlist->current->item;

    if ((vlist->current == vlist->tail)
	|| (vlist->current == vlist->head) )
	/* current pointer before beginning or beyond end, no action */
	;
    else {  /* free current node */

	vlist->current->prev->next = vlist->current->next;
	vlist->current->next->prev = vlist->current->prev;
	free_me = vlist->current;
	vlist->current = vlist->current->next;

	VFree (free_me);
	vlist->count--;
    }

    return return_me;
}

void VListConcat (VList vlist1, VList vlist2)
{
    VNodePtrType free_me, free_me_too;

    free_me = vlist1->tail;
    free_me_too = vlist2->head;

    vlist1->tail->prev->next = vlist2->head->next;
    vlist2->head->next->prev = vlist1->tail->prev;

    if ( vlist1->current == vlist1->tail )
	/* current pointer of vlist1 points beyond end,
	   set it to first node of vlist2 */
	vlist1->current = vlist2->head->next;

    vlist1->tail = vlist2->tail;
    vlist1->count += vlist2->count;

    VFree (free_me);
    VFree (free_me_too);
    VFree (vlist2);
}

void VListDestroy (VList vlist, void (*item_free) ())
{
   VPointer free_me;

   vlist->current = vlist->head->next;
   while ( vlist->current != vlist->tail )
   {
       free_me = VListRemove (vlist);
       (*item_free)(free_me);
   }

   VFree (vlist->head);
   VFree (vlist->tail);
   VFree (vlist);
}

VPointer VListTrim (VList vlist)
{
    VPointer return_me;

    return_me = VListLast(vlist);
    VListRemove (vlist);
    VListLast (vlist);

    return return_me;
}

VPointer VListSearch (VList vlist, int (*comp) (), VPointer comp_arg)
{
    if ( vlist->current == vlist->head )
	/* before beginning, go to next node */
	vlist->current = vlist->current->next;

    while ( vlist->current != vlist->tail ) {
	if ( (*comp)(vlist->current->item, comp_arg) )
	    /* a match is found */
	    return (vlist->current->item);
	else
	    vlist->current = vlist->current->next;
    }

    /* no match */
    return vlist->current->item;
}
