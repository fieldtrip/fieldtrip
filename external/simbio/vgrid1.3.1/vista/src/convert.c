/*
 *
 * convert.c: converting an image from one pixel representation to another by 
 *            simply copying (not mapping) pixel values
 * BRIAN Software Package Version 1.14
 *
 * version 1.1 (12/12/00)
 *
 * Copyright © 1996 Max-Planck-Institute of Cognitive Neuroscience. 
 * All rights reserved.
 *
 * Vista library is copyrighted of University of British Columbia.
 * Copyright © 1993, 1994 University of British Columbia.
 *
 */

#include "vista.h"

/*
 *  Some macros for converting one type to another.
 */

#define FillTable(table, op)	for (i = 0; i < VNumber (table); i++) table[i] = op;

#define Convert(src_type, result_type, op)				\
    {									\
	src_type *src_pp = (src_type *) src_first;			\
	result_type *res_pp = VImageData (result);			\
									\
	for (i = 0; i < npixels; i++) {					\
	    op;								\
	}								\
    }

#define Clip(test, limit)	if (pixel test limit) { pixel = limit; clipped = TRUE; }
#define ClipNeg			if (t < 0) { t = 0; clipped = TRUE; }

#define Bracket(t, lower, upper)					\
    if (t < lower) { t = lower; clipped = TRUE; } 			\
    else if (t > upper) { t = upper; clipped = TRUE; }

#define Cast(src_type, result_type, op)					\
    {									\
	src_type pixel, *src_pp = (src_type *) src_first;		\
	result_type *res_pp = VImageData (result);			\
									\
	for (i = 0; i < npixels; i++) {					\
	    pixel = *src_pp++;						\
	    op;								\
	    *res_pp++ = pixel;						\
	}								\
    }
#define Nothing

/*
 *  VConvertImageCopy
 *
 *  Converts an image from one pixel representation to another by copying
 *  (and perhaps rounding) pixel values, not mapping them from one range
 *  to another.
 *
 *  Note: This uses rint() to round a floating point value to the nearest
 *	  integer (under the default rounding mode). It isn't part of
 *	  ANSI C, but it's widely available.
 */

VImage VConvertImageCopy (VImage src, VImage dest, VBand band,
			  VRepnKind pixel_repn)
{
	VImage result;
	int npixels, i;
	VPointer src_first;
	VBoolean clipped = FALSE;

	/* If src already has the requested representation, simply copy: */
	if (pixel_repn == VPixelRepn (src))
		return (src == dest) ? src : VCopyImage (src, dest, band);

	/* Prepare to iterate over all source pixels: */
	if (!VSelectBand
	    ("VConvertImageCopy", src, band, &npixels, &src_first))
		return NULL;

	/* Check dest if it exists; create it if it doesn't: */
	result = VSelectDestImage ("VConvertImageCopy", dest,
				   band == VAllBands ? VImageNBands (src) : 1,
				   VImageNRows (src), VImageNColumns (src),
				   pixel_repn);
	if (!result)
		return NULL;

	/* Copy pixels, casting or rounding them: */
	switch (VPixelRepn (src)) {

	case VBitRepn:
		switch (VPixelRepn (result)) {
		case VUByteRepn:
			Cast (VBit, VUByte, Nothing);
			break;
		case VSByteRepn:
			Cast (VBit, VSByte, Nothing);
			break;
		case VShortRepn:
			Cast (VBit, VShort, Nothing);
			break;
		case VLongRepn:
			Cast (VBit, VLong, Nothing);
			break;
		case VFloatRepn:
			Cast (VBit, VFloat, Nothing);
			break;
		case VDoubleRepn:
			Cast (VBit, VDouble, Nothing);
			break;
		default:
			break;
		}
		break;

	case VUByteRepn:
		switch (VPixelRepn (result)) {
		case VBitRepn:
			Cast (VUByte, VBit, Clip (>, 1));
			break;
		case VSByteRepn:
			Cast (VUByte, VSByte, Clip (>, 127));
			break;
		case VShortRepn:
			Cast (VUByte, VShort, Nothing);
			break;
		case VLongRepn:
			Cast (VUByte, VLong, Nothing);
			break;
		case VFloatRepn:
			Cast (VUByte, VFloat, Nothing);
			break;
		case VDoubleRepn:
			Cast (VUByte, VDouble, Nothing);
			break;
		default:
			break;
		}
		break;

	case VSByteRepn:
		switch (VPixelRepn (result)) {
		case VBitRepn:
			Cast (VSByte, VBit, Bracket (pixel, 0, 1));
			break;
		case VUByteRepn:
			Cast (VSByte, VUByte, Clip (<, 0));
			break;
		case VShortRepn:
			Cast (VSByte, VShort, Nothing);
			break;
		case VLongRepn:
			Cast (VSByte, VLong, Nothing);
			break;
		case VFloatRepn:
			Cast (VSByte, VFloat, Nothing);
			break;
		case VDoubleRepn:
			Cast (VSByte, VDouble, Nothing);
			break;
		default:
			break;
		}
		break;

	case VShortRepn:
		switch (VPixelRepn (result)) {
		case VBitRepn:
			Cast (VShort, VBit, Bracket (pixel, 0, 1));
			break;
		case VUByteRepn:
			Cast (VShort, VUByte, Bracket (pixel, 0, 255));
			break;
		case VSByteRepn:
			Cast (VShort, VSByte, Bracket (pixel, -128, 127));
			break;
		case VLongRepn:
			Cast (VShort, VLong, Nothing);
			break;
		case VFloatRepn:
			Cast (VShort, VFloat, Nothing);
			break;
		case VDoubleRepn:
			Cast (VShort, VDouble, Nothing);
			break;
		default:
			break;
		}
		break;

	case VLongRepn:
		switch (VPixelRepn (result)) {
		case VBitRepn:
			Cast (VLong, VBit, Bracket (pixel, 0, 1));
			break;
		case VUByteRepn:
			Cast (VLong, VUByte, Bracket (pixel, 0, 255));
			break;
		case VSByteRepn:
			Cast (VLong, VSByte, Bracket (pixel, -128, 127));
			break;
		case VShortRepn:
			Cast (VLong, VShort, Bracket (pixel, -32768, 32767));
			break;
		case VFloatRepn:
			Cast (VLong, VFloat, Nothing);
			break;
		case VDoubleRepn:
			Cast (VLong, VDouble, Nothing);
			break;
		default:
			break;
		}
		break;

	case VFloatRepn:
		switch (VPixelRepn (result)) {
		case VBitRepn:
			Cast (VFloat, VBit, Bracket (pixel, 0, 1);
			      pixel = rint (pixel));
			break;
		case VUByteRepn:
			Cast (VFloat, VUByte, Bracket (pixel, 0, 255);
			      pixel = rint (pixel));
			break;
		case VSByteRepn:
			Cast (VFloat, VSByte, Bracket (pixel, -128, 127);
			      pixel = rint (pixel));
			break;
		case VShortRepn:
			Cast (VFloat, VShort, Bracket (pixel, -32768, 32767);
			      pixel = rint (pixel));
			break;
		case VLongRepn:
			Cast (VFloat, VLong,
			      Bracket (pixel, VRepnMinValue (VLongRepn),
				       VRepnMaxValue (VLongRepn));
				);
			break;
		case VDoubleRepn:
			Cast (VFloat, VDouble, Nothing);
			break;
		default:
			break;
		}
		break;

	case VDoubleRepn:
		switch (VPixelRepn (result)) {
		case VBitRepn:
			Cast (VDouble, VBit, Bracket (pixel, 0, 1);
			      pixel = rint (pixel));
			break;
		case VUByteRepn:
			Cast (VDouble, VUByte, Bracket (pixel, 0, 255);
			      pixel = rint (pixel));
			break;
		case VSByteRepn:
			Cast (VDouble, VSByte, Bracket (pixel, -128, 127);
			      pixel = rint (pixel));
			break;
		case VShortRepn:
			Cast (VDouble, VShort, Bracket (pixel, -32768, 32767);
			      pixel = rint (pixel));
			break;
		case VLongRepn:
			Cast (VDouble, VLong,
			      Bracket (pixel, VRepnMinValue (VLongRepn),
				       VRepnMaxValue (VLongRepn)));
			break;
		case VFloatRepn:
			Cast (VDouble, VFloat, Nothing);
			break;
		default:
			break;
		}

	default:
		break;
	}
	if (clipped)
		VWarning ("VConvertImageCopy: "
			  "Source pixel value exceeds destination range");

	VCopyImageAttrs (src, result);

	return result;
}


/*
 *  VConvertImageLinear
 *
 *  Converts an image from one pixel representation to another using a
 *  linear mapping between source and destination pixel values:
 *
 *	dest_pixel = src_pixel * a + b
 *
 *  for a pair a, b supplied as parameters.
 */

VImage VConvertImageLinear (VImage src, VImage dest, VBand band,
			    VRepnKind pixel_repn, double a, double b)
{
	VImage result;
	int npixels, i;
	VPointer src_first;
	VDouble a_min, a_max, d0, d1, t;
	VDouble d_min =
		VIsFloatPtRepn (pixel_repn) ? -1.0 :
		VRepnMinValue (pixel_repn);
	VDouble d_max =
		VIsFloatPtRepn (pixel_repn) ? 1.0 :
		VRepnMaxValue (pixel_repn);
	VBoolean clipped = FALSE;

	/* Prepare to iterate over all source pixels: */
	if (!VSelectBand
	    ("VConvertImageLinear", src, band, &npixels, &src_first))
		return NULL;

	/* Check dest if it exists; create it if it doesn't: */
	result = VSelectDestImage ("VConvertImageLinear", dest,
				   band == VAllBands ? VImageNBands (src) : 1,
				   VImageNRows (src), VImageNColumns (src),
				   pixel_repn);
	if (!result)
		return NULL;

	/* Map pixels from one representation to the other: */
	switch (VPixelRepn (src)) {

	case VBitRepn:
		d0 = b;
		d1 = a + b;
		if (pixel_repn != VFloatRepn && pixel_repn != VDoubleRepn) {
			Bracket (d0, d_min, d_max);
			Bracket (d1, d_min, d_max);
		}
		switch (VPixelRepn (result)) {
		case VBitRepn:
			{
				VBit dd0 = rint (d0), dd1 = rint (d1);
				Convert (VBit, VBit, *res_pp++ =
					 *src_pp++ ? dd1 : dd0);
			}
			break;
		case VUByteRepn:
			{
				VUByte dd0 = rint (d0), dd1 = rint (d1);
				Convert (VBit, VUByte, *res_pp++ =
					 *src_pp++ ? dd1 : dd0);
			}
			break;
		case VSByteRepn:
			{
				VSByte dd0 = rint (d0), dd1 = rint (d1);
				Convert (VBit, VUByte, *res_pp++ =
					 *src_pp++ ? dd1 : dd0);
			}
			break;
		case VShortRepn:
			{
				VShort dd0 = rint (d0), dd1 = rint (d1);
				Convert (VBit, VUByte, *res_pp++ =
					 *src_pp++ ? dd1 : dd0);
			}
			break;
		case VLongRepn:
			{
				VLong dd0 = rint (d0), dd1 = rint (d1);
				Convert (VBit, VUByte, *res_pp++ =
					 *src_pp++ ? dd1 : dd0);
			}
			break;
		case VFloatRepn:
			{
				Convert (VBit, VFloat, *res_pp++ =
					 *src_pp++ ? d1 : d0);
			}
			break;
		case VDoubleRepn:
			{
				Convert (VBit, VDouble, *res_pp++ =
					 *src_pp++ ? d1 : d0);
			}
			break;
		default:
			break;
		}
		break;

	case VUByteRepn:
		switch (VPixelRepn (result)) {
		case VBitRepn:
			{
				VBit table[256];
				FillTable (table, t = i * a + b;
					   Bracket (t, d_min, d_max);
					   table[i] = rint (t));
				Convert (VUByte, VBit, *res_pp++ =
					 table[*src_pp++]);
			}
			break;
		case VUByteRepn:
			{
				VUByte table[256];
				FillTable (table, t = i * a + b;
					   Bracket (t, d_min, d_max);
					   table[i] = rint (t));
				Convert (VUByte, VUByte, *res_pp++ =
					 table[*src_pp++]);
			}
			break;
		case VSByteRepn:
			{
				VSByte table[256];
				FillTable (table, t = i * a + b;
					   Bracket (t, d_min, d_max);
					   table[i] = rint (t));
				Convert (VUByte, VSByte, *res_pp++ =
					 table[*src_pp++]);
			}
			break;
		case VShortRepn:
			{
				VShort table[256];
				FillTable (table, t = i * a + b;
					   Bracket (t, d_min, d_max);
					   table[i] = rint (t));
				Convert (VUByte, VShort, *res_pp++ =
					 table[*src_pp++]);
			}
			break;
		case VLongRepn:
			{
				VLong table[256];
				FillTable (table, t = i * a + b;
					   Bracket (t, d_min, d_max);
					   table[i] = rint (t));
				Convert (VUByte, VLong, *res_pp++ =
					 table[*src_pp++]);
			}
			break;
		case VFloatRepn:
			{
				VFloat table[256];
				FillTable (table, table[i] = i * a + b);
				Convert (VUByte, VFloat, *res_pp++ =
					 table[*src_pp++]);
			}
			break;
		case VDoubleRepn:
			{
				VDouble table[256];
				FillTable (table, table[i] = i * a + b);
				Convert (VUByte, VFloat, *res_pp++ =
					 table[*src_pp++]);
			}
			break;
		default:
			break;
		}
		break;

	case VSByteRepn:
		switch (VPixelRepn (result)) {
		case VBitRepn:
			{
				VBit table[256];
				FillTable (table, t = (i - 128) * a + b;
					   Bracket (t, d_min, d_max);
					   table[i] = rint (t));
				Convert (VSByte, VBit, *res_pp++ =
					 table[*src_pp++ + 128]);
			}
			break;
		case VUByteRepn:
			{
				VUByte table[256];
				FillTable (table, t = (i - 128) * a + b;
					   Bracket (t, d_min, d_max);
					   table[i] = rint (t));
				Convert (VSByte, VUByte, *res_pp++ =
					 table[*src_pp++ + 128]);
			}
			break;
		case VSByteRepn:
			{
				VSByte table[256];
				FillTable (table, t = (i - 128) * a + b;
					   Bracket (t, d_min, d_max);
					   table[i] = rint (t));
				Convert (VSByte, VSByte, *res_pp++ =
					 table[*src_pp++ + 128]);
			}
			break;
		case VShortRepn:
			{
				VShort table[256];
				FillTable (table, t = (i - 128) * a + b;
					   Bracket (t, d_min, d_max);
					   table[i] = rint (t));
				Convert (VSByte, VShort, *res_pp++ =
					 table[*src_pp++ + 128]);
			}
			break;
		case VLongRepn:
			{
				VLong table[256];
				FillTable (table, t = (i - 128) * a + b;
					   Bracket (t, d_min, d_max);
					   table[i] = rint (t));
				Convert (VSByte, VLong, *res_pp++ =
					 table[*src_pp++ + 128]);
			}
			break;
		case VFloatRepn:
			{
				VFloat table[256];
				FillTable (table, table[i] =
					   (i - 128) * a + b);
				Convert (VSByte, VFloat, *res_pp++ =
					 table[*src_pp++ + 128]);
			}
			break;
		case VDoubleRepn:
			{
				VDouble table[256];
				FillTable (table, table[i] =
					   (i - 128) * a + b);
				Convert (VSByte, VFloat, *res_pp++ =
					 table[*src_pp++ + 128]);
			}
			break;
		default:
			break;
		}
		break;

	case VShortRepn:
		switch (VPixelRepn (result)) {
		case VBitRepn:
			Convert (VShort, VBit, t = *src_pp++ * a + b;
				 Bracket (t, d_min, d_max);
				 *res_pp++ = rint (t));
			break;
		case VUByteRepn:
			Convert (VShort, VUByte, t = *src_pp++ * a + b;
				 Bracket (t, d_min, d_max);
				 *res_pp++ = rint (t));
			break;
		case VSByteRepn:
			Convert (VShort, VSByte, t = *src_pp++ * a + b;
				 Bracket (t, d_min, d_max);
				 *res_pp++ = rint (t));
			break;
		case VShortRepn:
			Convert (VShort, VShort, t = *src_pp++ * a + b;
				 Bracket (t, d_min, d_max);
				 *res_pp++ = rint (t));
			break;
		case VLongRepn:
			Convert (VShort, VLong, t = *src_pp++ * a + b;
				 Bracket (t, d_min, d_max);
				 *res_pp++ = rint (t));
			break;
		case VFloatRepn:
			Convert (VShort, VFloat, *res_pp++ =
				 *src_pp++ * a + b);
			break;
		case VDoubleRepn:
			Convert (VShort, VDouble, *res_pp++ =
				 *src_pp++ * a + b);
			break;
		default:
			break;
		}
		break;

	case VLongRepn:
		switch (VPixelRepn (result)) {
		case VBitRepn:
			Convert (VLong, VBit, t = *src_pp++ * a + b;
				 Bracket (t, d_min, d_max);
				 *res_pp++ = rint (t));
			break;
		case VUByteRepn:
			Convert (VLong, VUByte, t = *src_pp++ * a + b;
				 Bracket (t, d_min, d_max);
				 *res_pp++ = rint (t));
			break;
		case VSByteRepn:
			Convert (VLong, VSByte, t = *src_pp++ * a + b;
				 Bracket (t, d_min, d_max);
				 *res_pp++ = rint (t));
			break;
		case VShortRepn:
			Convert (VLong, VShort, t = *src_pp++ * a + b;
				 Bracket (t, d_min, d_max);
				 *res_pp++ = rint (t));
			break;
		case VLongRepn:
			Convert (VLong, VLong, t = *src_pp++ * a + b;
				 Bracket (t, d_min, d_max);
				 *res_pp++ = rint (t));
			break;
		case VFloatRepn:
			Convert (VLong, VFloat, *res_pp++ =
				 *src_pp++ * a + b);
			break;
		case VDoubleRepn:
			Convert (VLong, VDouble, *res_pp++ =
				 *src_pp++ * a + b);
			break;
		default:
			break;
		}
		break;

	case VFloatRepn:
		switch (VPixelRepn (result)) {
		case VBitRepn:
			Convert (VFloat, VBit, t = *src_pp++ * a + b;
				 Bracket (t, d_min, d_max);
				 *res_pp++ = rint (t));
			break;
		case VUByteRepn:
			Convert (VFloat, VUByte, t = *src_pp++ * a + b;
				 Bracket (t, d_min, d_max);
				 *res_pp++ = rint (t));
			break;
		case VSByteRepn:
			Convert (VFloat, VSByte, t = *src_pp++ * a + b;
				 Bracket (t, d_min, d_max);
				 *res_pp++ = rint (t));
			break;
		case VShortRepn:
			Convert (VFloat, VShort, t = *src_pp++ * a + b;
				 Bracket (t, d_min, d_max);
				 *res_pp++ = rint (t));
			break;
		case VLongRepn:
			Convert (VFloat, VLong, t = *src_pp++ * a + b;
				 Bracket (t, d_min, d_max);
				 *res_pp++ = rint (t));
			break;
		case VFloatRepn:
			Convert (VFloat, VFloat, *res_pp++ =
				 *src_pp++ * a + b);
			break;
		case VDoubleRepn:
			Convert (VFloat, VDouble, *res_pp++ =
				 *src_pp++ * a + b);
			break;
		default:
			break;
		}
		break;

	case VDoubleRepn:
		switch (VPixelRepn (result)) {
		case VBitRepn:
			Convert (VDouble, VBit, t = *src_pp++ * a + b;
				 Bracket (t, d_min, d_max);
				 *res_pp++ = rint (t));
			break;
		case VUByteRepn:
			Convert (VDouble, VUByte, t = *src_pp++ * a + b;
				 Bracket (t, d_min, d_max);
				 *res_pp++ = rint (t));
			break;
		case VSByteRepn:
			Convert (VDouble, VSByte, t = *src_pp++ * a + b;
				 Bracket (t, d_min, d_max);
				 *res_pp++ = rint (t));
			break;
		case VShortRepn:
			Convert (VDouble, VShort, t = *src_pp++ * a + b;
				 Bracket (t, d_min, d_max);
				 *res_pp++ = rint (t));
			break;
		case VLongRepn:
			Convert (VDouble, VLong, t = *src_pp++ * a + b;
				 Bracket (t, d_min, d_max);
				 *res_pp++ = rint (t));
			break;
		case VFloatRepn:
			Convert (VDouble, VFloat, *res_pp++ =
				 *src_pp++ * a + b);
			break;
		case VDoubleRepn:
			Convert (VDouble, VDouble, *res_pp++ =
				 *src_pp++ * a + b);
			break;
		default:
			break;
		}

	default:
		break;
	}

	/* In the cases were a table was used for the conversion, we don't know
	   for sure whether any pixel values were clipped. This determines
	   with certainty: */
	if (clipped && (pixel_repn == VBitRepn ||
			pixel_repn == VUByteRepn
			|| pixel_repn == VSByteRepn)) {
		clipped = FALSE;
		VImageStats (src, band, &a_min, &a_max, NULL, NULL);
		t = a_min * a + b;
		if (t < d_min || t > d_max)
			clipped = TRUE;
		t = a_max * a + b;
		if (t < d_min || t > d_max)
			clipped = TRUE;
	}

	if (clipped)
		VWarning ("VConvertImageLinear: "
			  "Pixel value exceeds destination range");

	VCopyImageAttrs (src, result);

	return result;
}

static VImage ConstantImage (VImage src, VImage dest, VBand band,
			     VRepnKind pixel_repn, VDouble c)
{
	dest = VSelectDestImage ("VConvertImageOpt", dest,
				 band == VAllBands ? VImageNBands (src) : 1,
				 VImageNRows (src), VImageNColumns (src),
				 pixel_repn);
	if (dest)
		VFillImage (dest, VAllBands, c);
	return dest;
}

/*
 *  VConvertImageOpt
 *
 *  Converts an image from one pixel representation to another using a
 *  mapping from the actual range of source pixel values to full range of
 *  destination pixel values, thus preserving as much as possible of the
 *  information in the source image. If preserve_zero is TRUE, the mapping
 *  will be such that 0 maps to 0.
 */

VImage VConvertImageOpt (VImage src, VImage dest, VBand band,
			 VRepnKind pixel_repn, int method)
{
	VDouble a_min, a_max, d_min, d_max, v, a, b, e;

	/* Determine range of pixel values in the source and destination images: */
	VImageStats (src, band, &a_min, &a_max, NULL, NULL);
	d_min = VIsFloatPtRepn (pixel_repn) ? -1.0 :
		VRepnMinValue (pixel_repn);
	d_max = VIsFloatPtRepn (pixel_repn) ? 1.0 :
		VRepnMaxValue (pixel_repn);

	/* Set v = max ( | a_min |, | a_max | ): */
	v = a_max > 0.0 ? a_max : -a_max;
	if (a_min > v)
		v = a_min;
	else if (-a_min > v)
		v = -a_min;

	/* Choose an appropriate linear mapping from [a_min, a_max] to
	   [d_min, d_max]: */
	switch (method) {

	case 1:

		/* A mapping that preserves sign and maps 0 -> 0: */
		if (a_min == 0.0 && a_max == 0.0)
			return ConstantImage (src, dest, band, pixel_repn,
					      0.0);
		if (a_min < 0 && d_min >= 0.0)
			VWarning ("VConvertImageOpt: Negative pixel value converted " "to unsigned destination");
		return VConvertImageLinear (src, dest, band, pixel_repn,
					    d_max / v, 0.0);

	case 2:

		/* A mapping that preserves the sign of the actual pixel values mapped
		   but doesn't necessarily map 0 -> 0: */
		if (a_min == a_max)
			return ConstantImage (src, dest, band, pixel_repn,
					      a_max < 0.0 ? d_min : a_min >
					      0.0 ? d_max : 0.0);
		if (a_min < 0 && d_min >= 0.0)
			VWarning ("VConvertImageOpt: Negative pixel value converted " "to unsigned destination");
		if (a_max < 0.0) {
			e = 1.0 / (a_max - a_min + 1.0);
			a = (-e - d_min) / (a_max - a_min);
			b = d_min - a_min * a;
		} else if (a_max == 0.0) {
			a = d_min / a_min;
			b = 0.0;
		} else if (a_min == 0.0) {
			a = d_max / a_max;
			b = 0.0;
		} else if (a_min > 0.0) {
			e = 1.0 / (a_max - a_min + 1.0);
			a = (d_max - e) / (a_max - a_min);
			b = d_max - a_max * a;
		} else {
			a = d_max / v;
			b = 0.0;
		}
		return VConvertImageLinear (src, dest, band, pixel_repn,
					    (double)a, (double)b);

	case 3:

		/* A mapping from the range of actual source pixel values to the
		   range of possible destination pixel values, not necessarily
		   preserving either sign or zero: */
		if (a_min == a_max)
			return ConstantImage (src, dest, band, pixel_repn,
					      a_max < 0.0 ? d_min : a_min >
					      0.0 ? d_max : 0.0);
		a = (d_max - d_min) / (a_max - a_min);
		b = d_max - a_max * a;
		return VConvertImageLinear (src, dest, band, pixel_repn,
					    (double)a, (double)b);

	default:
		VError ("VConvertImageopt: Unknown conversion method %d",
			method);
	}
	return NULL;
}


/*
 *  VConvertImageRange
 *
 *  Converts an image from one pixel representation to another using a
 *  mapping from all possible source pixel values to all destination pixel
 *  values.
 *
 *  Note: In ANSI C, left shifting a signed, negative value may or may
 *	  produce a negative result. Where we now the value's positive,
 *	  we use >>; where it could be negative, we use /.
 */

VImage VConvertImageRange (VImage src, VImage dest, VBand band,
			   VRepnKind pixel_repn)
{
	VImage result;
	int npixels, i;
	double v;
	VPointer src_first;
	VBoolean clipped = FALSE;

	/* If src already has the requested representation, simply copy: */
	if (pixel_repn == VPixelRepn (src))
		return (src == dest) ? src : VCopyImage (src, dest, band);

	/* Prepare to iterate over all source pixels: */
	if (!VSelectBand
	    ("VConvertImageRange", src, band, &npixels, &src_first))
		return NULL;

	/* Check dest if it exists; create it if it doesn't: */
	result = VSelectDestImage ("VConvertImageRange", dest,
				   band == VAllBands ? VImageNBands (src) : 1,
				   VImageNRows (src), VImageNColumns (src),
				   pixel_repn);
	if (!result)
		return NULL;

	/* Map pixels from one representation to the other: */
	switch (VPixelRepn (src)) {

	case VBitRepn:
		switch (VPixelRepn (result)) {
		case VUByteRepn:
			Convert (VBit, VUByte, *res_pp++ =
				 *src_pp++ ? 255 : 0);
			break;
		case VSByteRepn:
			Convert (VBit, VSByte, *res_pp++ =
				 *src_pp++ ? 127 : 0);
			break;
		case VShortRepn:
			Convert (VBit, VShort, *res_pp++ =
				 *src_pp++ ? 32767 : 0);
			break;
		case VLongRepn:
			Convert (VBit, VLong, *res_pp++ =
				 *src_pp++ ? 0x7FFFFFFFl : 0);
			break;
		case VFloatRepn:
			Convert (VBit, VFloat,
				 *res_pp++ =
				 (*src_pp++ ? VFloatConst (1.0) :
				  VFloatConst (0.0)));
			break;
		case VDoubleRepn:
			Convert (VBit, VDouble,
				 *res_pp++ =
				 (*src_pp++ ? VDoubleConst (1.0) :
				  VDoubleConst (0.0)));
			break;
		default:
			break;
		}
		break;

	case VUByteRepn:
		switch (VPixelRepn (result)) {
		case VBitRepn:
			Convert (VUByte, VBit, *res_pp++ = *src_pp++ >> 7);
			break;
		case VSByteRepn:
			Convert (VUByte, VSByte, *res_pp++ = *src_pp++ >> 1);
			break;
		case VShortRepn:
			Convert (VUByte, VShort, *res_pp++ = *src_pp++ << 7);
			break;
		case VLongRepn:
			Convert (VUByte, VLong, *res_pp++ = *src_pp++ << 23);
			break;
		case VFloatRepn:
			{
				VFloat table[256];
				FillTable (table, (i / VFloatConst (255.0)));
				Convert (VUByte, VFloat, *res_pp++ =
					 table[*src_pp++]);
			}
			break;
		case VDoubleRepn:
			{
				VDouble table[256];
				FillTable (table, (i / VDoubleConst (255.0)));
				Convert (VUByte, VDouble, *res_pp++ =
					 table[*src_pp++]);
			}
			break;
		default:
			break;
		}
		break;

	case VSByteRepn:
		switch (VPixelRepn (result)) {

		case VBitRepn:
			{
				VSByte t;
				Convert (VSByte, VBit, t = *src_pp++;
					 ClipNeg; *res_pp++ = t >> 6);
			}
			break;
		case VUByteRepn:
			{
				VSByte t;
				Convert (VSByte, VUByte, t = *src_pp++;
					 ClipNeg; *res_pp++ = t << 1);
			}
			break;
		case VShortRepn:
			Convert (VSByte, VShort, *res_pp++ =
				 ((VShort) * src_pp++ << 8));
			break;
		case VLongRepn:
			Convert (VSByte, VLong, *res_pp++ =
				 ((VLong) * src_pp++ << 24));
			break;
		case VFloatRepn:
			{
				VFloat table[256];
				FillTable (table,
					   (i -
					    VFloatConst (128.0)) /
					   VFloatConst (128.0));
				Convert (VSByte, VFloat, *res_pp++ =
					 table[*src_pp++ + 128]);
			}
			break;
		case VDoubleRepn:
			{
				VDouble table[256];
				FillTable (table,
					   (i -
					    VDoubleConst (128)) /
					   VDoubleConst (128.0));
				Convert (VSByte, VDouble, *res_pp++ =
					 table[*src_pp++ + 128]);
			}
			break;
		default:
			break;
		}
		break;

	case VShortRepn:
		v = -VPixelMinValue (src);
		switch (VPixelRepn (result)) {
		case VBitRepn:
			{
				VShort t;
				Convert (VShort, VBit, t = *src_pp++;
					 ClipNeg; *res_pp++ = t >> 14);
			}
			break;
		case VUByteRepn:
			{
				VShort t;
				Convert (VShort, VUByte, t = *src_pp++;
					 ClipNeg; *res_pp++ = t >> 7);
			}
			break;
		case VSByteRepn:
			Convert (VShort, VSByte, *res_pp++ =
				 (*src_pp++ / 256));
			break;
		case VLongRepn:
			Convert (VShort, VLong, *res_pp++ =
				 ((VLong) * src_pp++ << 8));
			break;
		case VFloatRepn:
			Convert (VShort, VFloat, *res_pp++ = *src_pp++ / v);
			break;
		case VDoubleRepn:
			Convert (VShort, VDouble, *res_pp++ = *src_pp++ / v);
			break;
		default:
			break;
		}
		break;

	case VLongRepn:
		v = -VPixelMinValue (src);
		switch (VPixelRepn (result)) {

		case VBitRepn:
			{
				VLong t;
				Convert (VLong, VBit, t = *src_pp++;
					 ClipNeg; *res_pp++ = t >> 30);
			}
			break;
		case VUByteRepn:
			{
				VLong t;
				Convert (VLong, VUByte, t = *src_pp++;
					 ClipNeg; *res_pp++ = t >> 23);
			}
			break;
		case VSByteRepn:
			Convert (VLong, VSByte, *res_pp++ =
				 (*src_pp++ / (1 << 23)));
			break;
		case VShortRepn:
			Convert (VLong, VShort, *res_pp++ =
				 ((VLong) * src_pp++ / (1 << 16)));
			break;
		case VFloatRepn:
			Convert (VLong, VFloat, *res_pp++ = *src_pp++ / v);
			break;
		case VDoubleRepn:
			Convert (VLong, VDouble, *res_pp++ = *src_pp++ / v);
			break;
		default:
			break;
		}
		break;

	case VFloatRepn:
		{
			VFloat t, m = VPixelMaxValue (result);
			switch (VPixelRepn (result)) {

			case VBitRepn:
				Convert (VFloat, VBit, t = *src_pp++;
					 Bracket (t, 0, 1);
					 *res_pp++ = t >= 0.5);
				break;
			case VUByteRepn:
				Convert (VFloat, VUByte, t = *src_pp++;
					 Bracket (t, 0, 1);
					 *res_pp++ = rint (t * m));
				break;
			case VSByteRepn:
				Convert (VFloat, VSByte, t = *src_pp++;
					 Bracket (t, -1, 1);
					 *res_pp++ = rint (t * m));
				break;
			case VShortRepn:
				Convert (VFloat, VShort, t = *src_pp++;
					 Bracket (t, -1, 1);
					 *res_pp++ = rint (t * m));
				break;
			case VLongRepn:
				Convert (VFloat, VLong, t = *src_pp++;
					 Bracket (t, -1, 1);
					 *res_pp++ = rint (t * m));
				break;
			case VDoubleRepn:
				Convert (VFloat, VDouble, *res_pp++ =
					 *src_pp++);
				break;

			default:
				break;
			}
		}
		break;

	case VDoubleRepn:
		{
			VDouble t, m = VPixelMaxValue (result);
			switch (VPixelRepn (result)) {
			case VBitRepn:
				Convert (VDouble, VBit, t = *src_pp++;
					 Bracket (t, 0, 1);
					 *res_pp++ = t >= 0.5);
				break;
			case VUByteRepn:
				Convert (VDouble, VUByte, t = *src_pp++;
					 Bracket (t, 0, 1);
					 *res_pp++ = rint (t * m));
				break;
			case VSByteRepn:
				Convert (VDouble, VSByte, t = *src_pp++;
					 Bracket (t, -1, 1);
					 *res_pp++ = rint (t * m));
				break;
			case VShortRepn:
				Convert (VDouble, VShort, t = *src_pp++;
					 Bracket (t, -1, 1);
					 *res_pp++ = rint (t * m));
				break;
			case VLongRepn:
				Convert (VDouble, VLong, t = *src_pp++;
					 Bracket (t, -1, 1);
					 *res_pp++ = rint (t * m));
				break;
			case VFloatRepn:
				Convert (VDouble, VFloat, *res_pp++ =
					 *src_pp++);
				break;
			default:
				break;
			}

		}
		break;
	default:
		break;
	}

	if (clipped)
		VWarning ("VConvertImageRange: "
			  "Pixel value exceeds destination range");

	VCopyImageAttrs (src, result);

	return result;
}
