/*
 *
 * binarize.c: binarize image
 * BRIAN Software Package Version 1.14
 *
 * version 1.0. (01/04/01)
 *
 * Copyright © 1996 Max-Planck-Institute of Cognitive Neuroscience. 
 * All rights reserved.
 *
 * Vista library is copyrighted of University of British Columbia.
 * Copyright © 1993, 1994 University of British Columbia.
 *
 */

#include "vista.h"

VImage VBinarizeImage (VImage src, VImage dest, VBand band, double xmin,
		       double xmax)
{
	int nbands, npixels, i;
	VPointer src_pixels;
	VBit *dest_pixels;

	/* Locate the block of source image pixels: */
	if (!VSelectBand ("BinarizeImage", src, band, &npixels, &src_pixels))
		return NULL;

	/* Grumble about a threshold outside the source pixel range: */
	if (xmin < VPixelMinValue (src) || xmax > VPixelMaxValue (src))
		VWarning ("BinarizeImage: Thresholds (%g %g) are outside source pixel range [%g,%g]", xmin, xmax, VPixelMaxValue (src), VPixelMaxValue (src));

	/* Check the destination image, which should have VBit pixels.
	   If it is NULL, then create one of the appropriate size and type. */
	nbands = (band == VAllBands) ? VImageNBands (src) : 1;
	dest = VSelectDestImage ("BinarizeImage", dest, nbands,
				 VImageNRows (src), VImageNColumns (src),
				 VBitRepn);
	if (!dest)
		return NULL;

	/* Locate the block of destination image pixels: */
	dest_pixels = &VPixel (dest, 0, 0, 0, VBit);

	/* Define a macro, parameterized by source pixel type,
	   that can be instantiated for each possible type: */
#define Binarize(type)						\
    {								\
      type *src_pp = (type *) src_pixels;			\
      type tmin = (type) xmin,tmax = (type) xmax;		\
								\
      for (i = 0; i < npixels; i++) {				\
	 *dest_pixels++ = ((*src_pp >= tmin) && (*src_pp <= tmax)); \
         src_pp++;                                            \
      }                                                       \
    }

	/* Instantiate the macro once for each source pixel type: */
	switch (VPixelRepn (src)) {
	case VBitRepn:
		Binarize (VBit);
		break;
	case VUByteRepn:
		Binarize (VUByte);
		break;
	case VSByteRepn:
		Binarize (VSByte);
		break;
	case VShortRepn:
		Binarize (VShort);
		break;
	case VLongRepn:
		Binarize (VLong);
		break;
	case VFloatRepn:
		Binarize (VFloat);
		break;
	case VDoubleRepn:
		Binarize (VDouble);
		break;
	default:;
	}

#undef Binarize

	/* Let the destination inherit any attributes of the source image: */
	VCopyImageAttrs (src, dest);
	return dest;
}
