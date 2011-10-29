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

/*
 *  Information about built-in types.
 */
extern VTypeMethods VImageMethods;	/* in ImageType.c */
extern VTypeMethods VGraphMethods;      /* in GraphType.c */

static VRepnInfoRec builtin_repn_info[] = {
    { "unknown" },

    /* Integer and floating-point numbers: */
    { "bit",	   sizeof (VBit),	1, 0.0, 1.0 },
    { "ubyte",	   sizeof (VUByte),	8, 0.0, 255.0 },
    { "sbyte",	   sizeof (VSByte),	8, -128.0, 127.0 },
    { "short",	   sizeof (VShort),    16, -32768.0, 32767.0 },
    { "long",	   sizeof (VLong),     32, -2147483648.0, 2147483647.0 },
    { "float",	   sizeof (VFloat),    32, -3.40282346638528860e+38, 3.40282346638528860e+38 },
    { "double",	   sizeof (VDouble),   64, -1.797693134862315708e+308, 1.797693134862315708e+308 },

    /* Miscellaneous representations: */
    { "attr-list", sizeof (VAttrList),  0, 0.0, 0.0 },
    { "boolean",   sizeof (VBoolean),	1, 0.0, 0.0 },
    { "bundle",	   sizeof (VPointer),	0, 0.0, 0.0 },
    { "list",	   sizeof (VList),	0, 0.0, 0.0 },
    { "pointer",   sizeof (VPointer),	0, 0.0, 0.0 },
    { "string",	   sizeof (VString),	0, 0.0, 0.0 },

    /* Standard object types: */
    { "image",	   sizeof (VPointer),	0, 0.0, 0.0, &VImageMethods },
    { "graph",     sizeof (VPointer),   0, 0.0, 0.0, &VGraphMethods },
};

/* Keywords for representing TRUE or FALSE: */
VDictEntry VBooleanDict[] = {
    { "false",		FALSE },
    { "true",		TRUE },
    { "no",		FALSE },
    { "yes",		TRUE },
    { "off",		FALSE },
    { "on",		TRUE },
    { NULL }
};

/* Keywords for representing kinds of numeric representation: */
VDictEntry VNumericRepnDict[] = {
    { "bit",		VBitRepn },
    { "double",		VDoubleRepn },
    { "float",		VFloatRepn },
    { "long",		VLongRepn },
    { "sbyte",		VSByteRepn },
    { "short",		VShortRepn },
    { "ubyte",		VUByteRepn },
    { NULL }
};

VBoolean V_RequiredOpt, V_OptionalOpt;
VRepnInfoRec *VRepnInfo = builtin_repn_info;
static VRepnKind nRepnKinds = VNRepnKinds;

/* Maximum length of an error message: */
#define maxErrorMessageLength	500

VPointer VMalloc (size_t size)
{
    VPointer p;

    if (size == 0)
	return NULL;
    if (! (p = (VPointer) malloc (size)))
	VSystemError ("VMalloc: Memory allocation failure");
    return p;
}

VPointer VRealloc (VPointer p, size_t size)
{
    if (size == 0) {
	VFree (p);
	return NULL;
    }
    if (! p)
	return VMalloc (size);
    if (! (p = (VPointer) realloc (p, size)))
	VSystemError ("VRealloc: Memory allocation failure");
    return p;
}

VPointer VCalloc (size_t n, size_t size)
{
    VPointer p;

    if (n == 0 || size == 0)
	return NULL;
    if (! (p = (VPointer) calloc (n, size)))
	VSystemError ("VCalloc: Memory allocation failure");
    return p;
}

void VFree (VPointer p)
{
    if (p != NULL) {
      free ((char *) p);
      p = NULL;
    }
}

static void FormatMsg (char *buf, VStringConst severity, VStringConst format,
		       va_list *args, VStringConst extra)
{
    sprintf (buf, "%s: ", severity);
    vsprintf (buf + strlen (buf), format, *args);
    if (extra)
	sprintf (buf + strlen (buf), ": %s", extra);
    strcat (buf, ".\n");
}

void VError (VStringConst format, ...)
{
    va_list args;
    char buf[maxErrorMessageLength];

    va_start (args, format);
    FormatMsg (buf, "Fatal", format, & args, 0);
    va_end (args);
    fprintf(stderr, "%s\n", buf);
    exit (EXIT_FAILURE);
}

void VWarning (VStringConst format, ...)
{
    va_list args;
    char buf[maxErrorMessageLength];

    va_start (args, format);
    FormatMsg (buf, "Warning", format, & args, 0);
    va_end (args);
    fprintf(stderr, "%s\n", buf);
}

void VSystemError (VStringConst format, ...)
{
    va_list args;
    char buf[maxErrorMessageLength];

    va_start (args, format);
    FormatMsg (buf, "Fatal", format, & args, strerror (errno));
    va_end (args);
    fprintf(stderr, "%s\n", buf);
    exit (EXIT_FAILURE);
}

void VSystemWarning (VStringConst format, ...)
{
    va_list args;
    char buf[maxErrorMessageLength];

    va_start (args, format);
    FormatMsg (buf, "Warning", format, & args, strerror (errno));
    va_end (args);
    fprintf(stderr, "%s\n", buf);
}

VRepnKind VRegisterType (VStringConst name, VTypeMethods *methods)
{
    VRepnInfoRec *p;

    /* Move the existing type information into a bigger table: */
    if (VRepnInfo == builtin_repn_info) {
	VRepnInfo = VMalloc ((VNRepnKinds + 1) * sizeof (VRepnInfoRec));
	memcpy(VRepnInfo, builtin_repn_info, VNRepnKinds * sizeof (VRepnInfoRec));
    } else
	VRepnInfo =
	    VRealloc (VRepnInfo, (nRepnKinds + 1) * sizeof (VRepnInfoRec));

    /* Write the new type's info into the last table entry: */
    p = VRepnInfo + nRepnKinds;
    p->name = VNewString (name);
    p->size = p->precision = p->min_value = p->max_value = 0.0;
    p->methods = methods;

    return nRepnKinds++;
}

VRepnKind VLookupType (VStringConst name)
{
    VRepnKind repn;

    for (repn = VUnknownRepn; repn < nRepnKinds; repn++)
	if (strcmp (VRepnInfo[repn].name, name) == 0)
	    return repn;
    return VUnknownRepn;
}

VDictEntry *VLookupDictKeyword (VDictEntry *dict, VStringConst keyword)
{
    if (dict)
	for ( ; dict->keyword; dict++)
	    if (strcmp (keyword, dict->keyword) == 0)
		return dict;
    return NULL;
}

VDictEntry *VLookupDictValue (VDictEntry *dict, VRepnKind repn, ...)
{
    va_list args;
    VLong i_value = 0;
    VDouble f_value = 0.0;
    VString s_value = NULL;
    VBoolean i_valid;

    /* Unravel the arguments passed: */
    if (! dict)
	return NULL;
    va_start (args, repn);
    switch (repn) {
    case VBitRepn: i_value = va_arg (args, VBitPromoted); break;
    case VUByteRepn: i_value = va_arg (args, VUBytePromoted); break;
    case VSByteRepn: i_value = va_arg (args, VSBytePromoted); break;
    case VShortRepn: i_value = va_arg (args, VShortPromoted); break;
    case VLongRepn: i_value = va_arg (args, VLongPromoted); break;
    case VFloatRepn: f_value = va_arg (args, VFloatPromoted); break;
    case VDoubleRepn: f_value = va_arg (args, VDoublePromoted); break;
    case VBooleanRepn: i_value = va_arg (args, VBooleanPromoted); break;
    case VStringRepn: s_value = va_arg (args, VString); break;
    default:
	VError ("VLookupDictValue: Can't lookup %s value", VRepnName (repn));
    }
    va_end (args);

    /* Search the dictionary by value: */
    switch (repn) {

    case VBitRepn:
    case VUByteRepn:
    case VSByteRepn:
    case VShortRepn:
    case VLongRepn:
    case VBooleanRepn:
	for ( ; dict->keyword; dict++) {

	    /* Is the entry's value only stored as a string? */
	    if (dict->svalue && ! dict->icached) {

		/* Yes -- try to convert the string to an integer, and
		   cache that value: */
		if (! VDecodeAttrValue (dict->svalue, NULL,
					VLongRepn, & dict->ivalue))
		    break;
		dict->icached = TRUE;
	    }

	    /* Test against the integer value stored in the entry: */
	    if (i_value == dict->ivalue)
		return dict;
	}
	break;

    case VFloatRepn:
    case VDoubleRepn:
	for ( ; dict->keyword; dict++) {

	    /* Does the entry include a cached floating point value? */
	    if (! dict->fcached) {

		/* No -- obtain it from an integer or string value: */
		if (dict->svalue) {
		    if (! VDecodeAttrValue (dict->svalue, NULL,
					    VDoubleRepn, & dict->fvalue))
			break;
		} else dict->fvalue = dict->ivalue;
		dict->fcached = TRUE;
	    }

	    /* Test against the cached float value now stored in the entry: */
	    if (f_value == dict->fvalue)
		return dict;
	}
	break;

    case VStringRepn:

	/* In case we're searching a dictionary with only integer values
	   stored, try to convert the supplied string value to an integer: */
	i_valid = VDecodeAttrValue (s_value, NULL, VLongRepn, & i_value);

	for ( ; dict->keyword; dict++) {

	    /* If the entry includes a string value, compare with it: */
	    if (dict->svalue) {
		if (strcmp (s_value, dict->svalue) == 0)
		    return dict;
	    }

	    /* Otherwise, compare with its integer value: */
	    else if (i_valid && i_value == dict->ivalue)
		return dict;
	}
	break;

    default:
	break;
    }
    return NULL;
}

static VPackOrder MachineByteOrder  ()
{
    union {
	short s;
	char c[sizeof (short)];
    } u;

    u.s = 1;
    if (u.c[0] == 1)
	return VLsbFirst;
    if (u.c[sizeof (short) - 1] != 1)
	VError ("VPackImage or VUnpackImage: Byte order not recognized");
    return VMsbFirst;
}

static void SwapBytes (size_t nels, size_t elsize, char *data)
{
    int i;
    char *pl, *pu, byte;

    for (i = 0; i < nels; i++, data += elsize)
	for (pl = data, pu = data + elsize - 1; pl < pu; pl++, pu--) {
	    byte = *pl;
	    *pl = *pu;
	    *pu = byte;
	}
}

VBoolean VPackData (VRepnKind repn,
		    size_t nels, VPointer unpacked, VPackOrder packed_order,
		    size_t *length, VPointer *packed, VBoolean *alloced)
{
    VPackOrder unpacked_order;
    size_t unpacked_elsize = VRepnSize (repn) * CHAR_BIT;
    size_t packed_elsize = VRepnPrecision (repn);
    size_t packed_length = (nels * packed_elsize + 7) / 8;

    /* If space for the packed data was supplied, ensure there's
       enough of it: */
    if (! alloced && packed_length > *length) {
	VWarning ("VPackData: Insufficient space for packed data");
	return FALSE;
    }
    *length = packed_length;

    /* Determine the present machine's internal byte order: */
    unpacked_order = MachineByteOrder ();

    /* If the desired byte order matches that of the present machine's, and
       the unpacked and packed data element sizes are identical,
       just return the unpacked data: */
    if (unpacked_order == packed_order && unpacked_elsize == packed_elsize) {
	if (alloced) {
	    *packed = unpacked;
	    *alloced = FALSE;
	} else if (unpacked != packed)
	    memcpy (*packed, unpacked, packed_length);
	return TRUE;
    }

    /* Allocate a buffer for the packed data if none was provided: */
    if (alloced) {
	*packed = VMalloc (packed_length);
	*alloced = TRUE;
    }

    /* Pack data elements into the buffer: */
    if (unpacked_elsize == packed_elsize) {
    
	/* If the packed and unpacked are the same size, do a straight copy: */
	if (unpacked != *packed)
	    memcpy (*packed, unpacked, packed_length);
	
	/* Swap bytes if necessary: */
	if (packed_order != unpacked_order && packed_elsize > 8)
	    SwapBytes (nels, packed_elsize / 8, (char *) *packed);

    } else if (packed_elsize == 1) {

	/* If the elements are VBits, this packs them: */
	VPackBits (nels, packed_order, (VBit *) unpacked, (char *) *packed);

    } else

	/* Packing multi-byte integers or floats is currently not supported: */
	VError ("VPackData: Packing %s from %d to %d bits is not supported",
		VRepnName (repn), unpacked_elsize, packed_elsize);

    return TRUE;
}

VBoolean VUnpackData (VRepnKind repn,
		      size_t nels, VPointer packed, VPackOrder packed_order,
		      size_t *length, VPointer *unpacked, VBoolean *alloced)
{
    VPackOrder unpacked_order;
    size_t unpacked_elsize = VRepnSize (repn) * CHAR_BIT;
    size_t packed_elsize = VRepnPrecision (repn);
    size_t unpacked_length = nels * VRepnSize (repn);

    /* If a space for the unpacked data was supplied, ensure there's
       enough of it: */
    if (! alloced && unpacked_length > *length) {
	VWarning ("VUnpackData: Insufficient space for unpacked data");
	return FALSE;
    }
    *length = unpacked_length;

    /* Determine the present machine's internal byte order: */
    unpacked_order = MachineByteOrder ();

    /* If the desired byte order matches that of the present machine's, and
       the unpacked and packed data element sizes are identical,
       just return the packed data: */
    if (unpacked_order == packed_order && unpacked_elsize == packed_elsize) {
	if (alloced) {
	    *unpacked = packed;
	    *alloced = FALSE;
	} else if (packed != *unpacked)
	    memcpy (*unpacked, packed, unpacked_length);
	return TRUE;
    }

    /* Unpack data elements into the buffer: */
    if (packed_elsize == unpacked_elsize) {

	/* If the packed and unpacked are the same size, do a straight copy: */
	if (packed != *unpacked)
	    memcpy (*unpacked, packed, unpacked_length);

	/* Swap bytes if necessary: */
	if (packed_order != unpacked_order && packed_elsize > 8)
	    SwapBytes (nels, packed_elsize / 8, (char *) *unpacked);

    } else if (packed_elsize == 1) {

	/* If the elements are VBits, this unpacks them: */
	VUnpackBits (nels, packed_order, (char *) packed, (char *) *unpacked);

    } else 
    
	/* Unpacking multi-byte integers or floats is currently not
	   supported: */
	VError ("VUnpackData: "
		"Unpacking %s from %d to %d bits is not supported",
		VRepnName (repn), packed_elsize, unpacked_elsize);

    return TRUE;
}

void VPackBits (size_t nels, VPackOrder packed_order,
		VBit *unpacked, char *packed)
{
    int bit;
    char byte;

    if (packed_order == VLsbFirst)
	while (nels > 0) {
	    byte = 0;
	    for (bit = 0; bit < 8 && nels > 0; nels--, bit++)
		if (*unpacked++)
		    byte |= (1 << bit);
	    *packed++ = byte;
	}
    else
	while (nels > 0) {
	    byte = 0;
	    for (bit = 7; bit >= 0 && nels > 0; nels--, bit--)
		if (*unpacked++)
		    byte |= (1 << bit);
	    *packed++ = byte;
	}
}

void VUnpackBits (size_t nels, VPackOrder packed_order,
		  char *packed, VBit *unpacked)
{
    int bit;
    char byte;

    /* Compute the position of the first bit to be unpacked, which is the
       last bit of the vector: */
    bit = (nels + 7) % 8;
    if (packed_order == VMsbFirst)
	bit = 7 - bit;

    /* Unpack bits from last to first. For each byte to be unpacked: */
    packed += (nels + 7) / 8;
    unpacked += nels;
    if (packed_order == VLsbFirst)
	while (nels > 0) {
	    byte = *--packed;

	    /* For each bit to be unpacked within that byte: */
	    for ( ; bit >= 0 && nels > 0; nels--, bit--) {
		
		/* Unpack a bit: */
		*--unpacked = (byte >> bit) & 1;
	    }
	    bit = 7;
	}
    else
	while (nels > 0) {
	    byte = *--packed;

	    /* For each bit to be unpacked within that byte: */
	    for ( ; bit < 8 && nels > 0; nels--, bit++) {
		
		/* Unpack a bit: */
		*--unpacked = (byte >> bit) & 1;
	    }
	    bit = 0;
	}
}

