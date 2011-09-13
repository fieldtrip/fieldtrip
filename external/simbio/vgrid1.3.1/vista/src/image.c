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

/* Keywords for representing band interpretation values: */
VDictEntry VBandInterpDict[] = {
  { "complex",		VBandInterpComplex },
  { "gradient",		VBandInterpGradient },
  { "intensity",	VBandInterpIntensity },
  { "orientation",	VBandInterpOrientation },
  { "rgb",		VBandInterpRGB },
  { "stereo_pair",	VBandInterpStereoPair },
  { NULL }
};

static VPointer VImageCopyMethod (VPointer value)
{
  return VCopyImage ((VImage) value, NULL, VAllBands);
}

static VPointer VImageDecodeMethod (VStringConst name, VBundle b)
{
  VImage image;
  VLong nbands, nrows, ncolumns, pixel_repn;
  VLong nframes, nviewpoints, ncolors, ncomponents;
  VAttrList list;
  size_t length;

#define Extract(name, dict, locn, required)	\
	VExtractAttr (b->list, name, dict, VLongRepn, & locn, required)

  /* Extract the number of bands, rows, columns, pixel repn, etc.: */
  nbands = nframes = nviewpoints = ncolors = ncomponents = 1;	/* defaults */
  if (! Extract (VNBandsAttr, NULL, nbands, FALSE) ||
      ! Extract (VNRowsAttr, NULL, nrows, TRUE) ||
      ! Extract (VNColumnsAttr, NULL, ncolumns, TRUE) ||
      ! Extract (VRepnAttr, VNumericRepnDict, pixel_repn, TRUE) ||
      ! Extract (VNFramesAttr, NULL, nframes, FALSE) ||
      ! Extract (VNViewpointsAttr, NULL, nviewpoints, FALSE) ||
      ! Extract (VNColorsAttr, NULL, ncolors, FALSE) ||
      ! Extract (VNComponentsAttr, NULL, ncomponents, FALSE))
    return NULL;

  /* Ensure that nbands == nframes * nviewpoints * ncolors * ncomponents.
     For backwards compatibility, set ncomponents to nbands if nbands != 1
     but nframes == nviewpoints == ncolors == ncomponents == 1. */
  if (nbands != nframes * nviewpoints * ncolors * ncomponents) {
    if (nbands != 1 && nframes == 1 && nviewpoints == 1 &&
	ncolors == 1 && ncomponents == 1)
      ncomponents = nbands;
    else {
      VWarning ("VImageDecodeMethod: %s image has inconsistent nbands",
		name);
      return NULL;
    }
  }

  /* Create an image with the specified properties: */
  if (! (image = VCreateImage ((int) nbands, (int) nrows, (int) ncolumns,
			       (VRepnKind) pixel_repn)))
    return NULL;
  VImageNFrames (image) = nframes;
  VImageNViewpoints (image) = nviewpoints;
  VImageNColors (image) = ncolors;
  VImageNComponents (image) = ncomponents;

  /* Give it whatever attributes remain: */
  list = VImageAttrList (image);
  VImageAttrList (image) = b->list;
  b->list = list;

  /* Check that the expected amount of binary data was read: */
  length = VImageNPixels (image);
  if (VPixelRepn (image) == VBitRepn)
    length = (length + 7) / 8;
  else length *= VPixelPrecision (image) / 8;
  if (length != b->length) {
    VWarning ("VImageDecodeMethod: %s image has wrong data length", name);
  Fail:   VDestroyImage (image);
    return NULL;
  }

  /* Unpack the binary pixel data: */
  length = VImageSize (image);
  if (! VUnpackData (VPixelRepn (image), VImageNPixels (image),
		     b->data, VMsbFirst, & length, & VImageData (image),
		     NULL))
    goto Fail;
  return image;

#undef Extract
}

static VAttrList VImageEncodeAttrMethod (VPointer value, size_t *lengthp)
{
  VImage image = value;
  VAttrList list;
  size_t length;

#define OptionallyPrepend(value, name)				\
	if (value != 1)							\
	    VPrependAttr (list, name, NULL, VLongRepn, (VLong) value)

  /* Temporarily prepend several attributes to the image's attribute list: */
  if ((list = VImageAttrList (image)) == NULL)
    list = VImageAttrList (image) = VCreateAttrList ();
  VPrependAttr (list, VRepnAttr, VNumericRepnDict,
		VLongRepn, (VLong) VPixelRepn (image));
  VPrependAttr (list, VNColumnsAttr, NULL,
		VLongRepn, (VLong) VImageNColumns (image));
  VPrependAttr (list, VNRowsAttr, NULL,
		VLongRepn, (VLong) VImageNRows (image));
  OptionallyPrepend (VImageNComponents (image), VNComponentsAttr);
  OptionallyPrepend (VImageNColors (image), VNColorsAttr);
  OptionallyPrepend (VImageNViewpoints (image), VNViewpointsAttr);
  OptionallyPrepend (VImageNFrames (image), VNFramesAttr);
  OptionallyPrepend (VImageNBands (image), VNBandsAttr);

  /* Compute the file space needed for the image's binary data: */
  length = VImageNPixels (image);
  if (VPixelRepn (image) == VBitRepn)
    length = (length + 7) / 8;
  else length *= VPixelPrecision (image) / 8;
  *lengthp = length;

  return list;

#undef OptionallyPrepend
}

static VPointer VImageEncodeDataMethod (VPointer value, VAttrList list,
					size_t length, VBoolean *free_itp)
{
  VImage image = value;
  VAttrListPosn posn;
  size_t len;
  VPointer ptr;

  /* Remove the attributes prepended by the VImageEncodeAttrsMethod: */
  for (VFirstAttr (list, & posn);
       strcmp (VGetAttrName (& posn), VRepnAttr) != 0;
       VDeleteAttr (& posn)) ;
  VDeleteAttr (& posn);

  /* Pack and return pixel data: */
  if (! VPackData (VPixelRepn (image), VImageNPixels (image),
		   VImageData (image), VMsbFirst, & len, & ptr, free_itp))
    return NULL;
  if (len != length)
    VError ("VImageEncodeDataMethod: Encoded data has unexpected length");
  return ptr;
}

/* Used in Type.c to register this type: */
VTypeMethods VImageMethods = {
  VImageCopyMethod,			/* copy a VImage */
  (VDestroyMethod *) VDestroyImage,	/* destroy a VImage */
  VImageDecodeMethod,			/* decode a VImage's value */
  VImageEncodeAttrMethod,		/* encode a VImage's attr list */
  VImageEncodeDataMethod		/* encode a VImage's binary data */
};

VImage VCreateImage (int nbands, int nrows, int ncolumns, VRepnKind pixel_repn)
{
  size_t row_size = ncolumns * VRepnSize (pixel_repn);
  size_t data_size = nbands * nrows * row_size;
  size_t row_index_size = nbands * nrows * sizeof (char *);
  size_t band_index_size = nbands * sizeof (char **);
  size_t pixel_size;
  char *p;
  VImage image;
  int band, row;

#define AlignUp(v, b) ((((v) + (b) - 1) / (b)) * (b))

  /* Check parameters: */
  if (nbands < 1) {
    VWarning ("VCreateImage: Invalid number of bands: %d", (int) nbands);
    return NULL;
  }
  if (nrows < 1) {
    VWarning ("VCreateImage: Invalid number of rows: %d", (int) nrows);
    return NULL;
  }
  if (ncolumns < 1) {
    VWarning ("VCreateImage: Invalid number of columns: %d",
	      (int) ncolumns);
    return NULL;
  }
  if (pixel_repn != VBitRepn && pixel_repn != VUByteRepn &&
      pixel_repn != VSByteRepn && pixel_repn != VShortRepn &&
      pixel_repn != VLongRepn && pixel_repn != VFloatRepn &&
      pixel_repn != VDoubleRepn) {
    VWarning ("VCreateImage: Invalid pixel representation: %d",
	      (int) pixel_repn);
    return NULL;
  }

  /* Allocate memory for the VImage, its indices, and pixel values, while
     padding enough to ensure pixel values are appropriately aligned: */
  pixel_size = VRepnSize (pixel_repn);
  p = VMalloc (AlignUp (sizeof (VImageRec) + row_index_size +
			band_index_size, pixel_size) + data_size);

  /* Initialize the VImage: */
  image = (VImage) p;
  image->nbands = nbands;
  image->nrows = nrows;
  image->ncolumns = ncolumns;
  image->flags = VImageSingleAlloc;
  image->pixel_repn = pixel_repn;
  image->attributes = VCreateAttrList ();
  image->band_index = (VPointer **) (p += sizeof (VImageRec));
  image->row_index = (VPointer *) (p += band_index_size);
  image->data = (VPointer) AlignUp ((long) p + row_index_size, pixel_size);
  image->nframes = nbands;
  image->nviewpoints = image->ncolors = image->ncomponents = 1;

  /* Initialize the indices: */
  for (band = 0; band < nbands; band++)
    image->band_index[band] = image->row_index + band * nrows;
  for (row = 0, p = image->data; row < nbands * nrows; row++, p += row_size)
    image->row_index[row] = p;

  return image;

#undef AlignUp
}

VImage VCreateImageLike (VImage src)
{
  return VCopyImageAttrs (src, NULL);
}	

void VDestroyImage (VImage image)
{
  if (image == NULL)
    return;

  if (! (image->flags & VImageSingleAlloc)) {
    if (image->data != NULL) {
      VFree (image->data);
    }
    if ((VPointer) image->row_index != NULL) {
      VFree ((VPointer) image->row_index);
    }
    if ((VPointer) image->band_index != NULL) {
      VFree ((VPointer) image->band_index);
    }
  }

  if (VImageAttrList(image) != NULL) {
    VDestroyAttrList (VImageAttrList (image));
  }
  if ((VPointer) image != NULL) {
    VFree ((VPointer) image);
  }
}

VDouble VGetPixel (VImage image, int band, int row, int column)
{
  VPointer p = VPixelPtr (image, band, row, column);
  switch (VPixelRepn (image)) {
  case VBitRepn:	return (VDouble) *(VBit *) p;
  case VUByteRepn:	return (VDouble) *(VUByte *) p;
  case VSByteRepn:	return (VDouble) *(VSByte *) p;
  case VShortRepn:	return (VDouble) *(VShort *) p;
  case VLongRepn:	return (VDouble) *(VLong *) p;
  case VFloatRepn:	return (VDouble) *(VFloat *) p;
  case VDoubleRepn:	return (VDouble) *(VDouble *) p;
  default:		VError("VGetPixel: %s images not supported", VPixelRepnName (image));
  }
  return 0.0;		/* to make lint happy */
}

void VSetPixel(VImage image, int band, int row, int column, VDoublePromoted value)
{
  VPointer p = VPixelPtr (image, band, row, column);
  switch (VPixelRepn (image)) {
  case VBitRepn:	*(VBit *) p = value; break;
  case VUByteRepn:	*(VUByte *) p = value; break;
  case VSByteRepn:	*(VSByte *) p = value; break;
  case VShortRepn:	*(VShort *) p = value; break;
  case VLongRepn:	*(VLong *) p = value; break;
  case VFloatRepn:	*(VFloat *) p = value; break;
  case VDoubleRepn:	*(VDouble *) p = value; break;
  default:		VError("VSetPixel: %s images not supported", VPixelRepnName (image));
  }
}

VImage VCopyImage (VImage src, VImage dest, VBand band)
{
  VImage result;

  if (src == dest &&
      (band == VAllBands || (band == 0 && VImageNBands (src) == 1)))
    return src;

  if ((result = VCopyImagePixels (src, dest, band)) != 0)
    VCopyImageAttrs (src, result);
  return result;
}

VImage VCopyImageAttrs (VImage src, VImage dest)
{
  VAttrList list;

  if (src == dest)
    return dest;
  if (! dest) {
    dest = VCreateImage (VImageNBands (src), VImageNRows (src),
			 VImageNColumns (src), VPixelRepn (src));
    if (! dest)
      return NULL;
  }

  /* Clone the source image's attribute list if it isn't empty: */
  if (! VAttrListEmpty (VImageAttrList (src))) {
    list = VImageAttrList (dest);
    VImageAttrList (dest) = VCopyAttrList (VImageAttrList (src));
  } else if (! VAttrListEmpty (VImageAttrList (dest))) {
    list = VImageAttrList (dest);
    VImageAttrList (dest) = VCreateAttrList ();
  } else list = NULL;
  if (list)
    VDestroyAttrList (list);

  /* Preserve band interpretation attributes only if the source and
     destination images have the same number of bands: */
  if (VImageNBands (src) > 1 && VImageNBands (dest) == VImageNBands (src)) {
    VImageNFrames (dest) = VImageNFrames (src);
    VImageNViewpoints (dest) = VImageNViewpoints (src);
    VImageNColors (dest) = VImageNColors (src);
    VImageNComponents (dest) = VImageNComponents (src);
  } else {
    VExtractAttr (VImageAttrList (dest), VFrameInterpAttr, NULL,
		  VBitRepn, NULL, FALSE);
    VExtractAttr (VImageAttrList (dest), VViewpointInterpAttr, NULL,
		  VBitRepn, NULL, FALSE);
    VExtractAttr (VImageAttrList (dest), VColorInterpAttr, NULL,
		  VBitRepn, NULL, FALSE);
    VExtractAttr (VImageAttrList (dest), VComponentInterpAttr, NULL,
		  VBitRepn, NULL, FALSE);
    VImageNComponents (dest) = VImageNColors (dest) = 
      VImageNViewpoints (dest) = 1;
    VImageNFrames (dest) = VImageNBands (dest);
  }
  return dest;
}

VImage VCopyImagePixels (VImage src, VImage dest, VBand band)
{
  int npixels;
  VPointer src_pixels;
  VImage result;

  /* Locate the source and destination of the copy: */
  if (! VSelectBand ("VCopyImagePixels", src, band, & npixels, & src_pixels))
    return NULL;
  result = VSelectDestImage ("VCopyImagePixels", dest,
			     band == VAllBands ? VImageNBands (src) : 1,
			     VImageNRows (src), VImageNColumns (src),
			     VPixelRepn (src));
  if (! result)
    return NULL;

  /* Copy pixel values from src to dest: */
  memcpy (VImageData (result), src_pixels, npixels * VPixelSize (src));

  return result;
}

VBoolean VCopyBand (VImage src, VBand src_band, VImage dest, VBand dest_band)
{
  int nbands, src_npixels, dest_npixels;
  VPointer src_pixels, dest_pixels;

  /* The destination image must exist: */
  if (! dest) {
    VWarning ("VCopyBand: No destination specified");
    return FALSE;
  }

  /* VAllBands not accepted for destination band: */
  if (dest_band < 0 || dest_band >= VImageNBands (dest)) {
    VWarning ("VCopyBand: Band %d referenced in image of %d bands",
	      (int) dest_band, (int) VImageNBands (dest));
    return FALSE;
  }

  /* Ensure that the destination image has the appropriate dimensions
     and pixel representation: */
  nbands = dest_band;
  if (src_band == VAllBands)
    nbands += VImageNBands (src) - 1;
  if (nbands < VImageNBands (dest))
    nbands = VImageNBands (dest);
  if (! VSelectDestImage ("VCopyBand", dest, nbands, VImageNRows (src),
			  VImageNColumns (src), VPixelRepn (src)))
    return FALSE;

  /* Locate the specified source and destination bands: */
  if (! VSelectBand ("VCopyBand", src, src_band,
		     & src_npixels, & src_pixels))
    return FALSE;
  if (! VSelectBand ("VCopyBand", dest, dest_band,
		     & dest_npixels, & dest_pixels))
    return FALSE;

  /* Copy from the source band to the destination band: */
  memcpy (dest_pixels, src_pixels, src_npixels * VPixelSize (src));

  return TRUE;
}

VImage VCombineBands (int nels, VImage src_images[], VBand src_bands[],
		      VImage dest)
{
  int n, i;
  VImage result, src = src_images[0];

  /* Count the number of bands needed in the destination image: */
  for (i = n = 0; i < nels; i++)
    n += (src_bands[i] == VAllBands) ? VImageNBands (src_images[i]) : 1;

  /* Check or allocate the destination image: */
  result = VSelectDestImage ("VCombineBands", dest, n,
			     VImageNRows (src), VImageNColumns (src), 
			     VPixelRepn (src));
  if (! result)
    return NULL;

  /* Copy each source band into the destination image: */
  for (i = n = 0; i < nels; i++) {
    if (! VCopyBand (src_images[i], src_bands[i], result, n)) {
      if (result != dest)
	VDestroyImage (result);
      return NULL;
    }
    n += (src_bands[i] == VAllBands) ? VImageNBands (src_images[i]) : 1;
  }
  return result;
}

VImage VCombineBandsVa (VImage dest, ...)
{
  va_list args;
  VImage src, result;
  int nbands;
  VBand src_band, dest_band;

  /* Count the number of bands to be combined: */
  va_start (args, dest);
  for (nbands = 0; (src = va_arg (args, VImage)) != 0; nbands +=
	 (va_arg (args, VBand) == VAllBands) ? VImageNBands (src) : 1) ;
  va_end (args);

  /* Check or allocate the destination image: */
  va_start (args, dest);
  src = va_arg (args, VImage);
  va_end (args);
  result = VSelectDestImage ("VCombineBandsVa", dest, nbands,
			     VImageNRows (src), VImageNColumns (src),
			     VPixelRepn (src));
  if (! result)
    return NULL;

  /* Copy each source band into the destination image: */
  va_start (args, dest);
  for (dest_band = 0; (src = va_arg (args, VImage)) != 0; ) {
    src_band = va_arg (args, VBand);
    if (! VCopyBand (src, src_band, result, dest_band)) {
      if (result != dest)
	VDestroyImage (result);
      return NULL;
    }
    dest_band += (src_band == VAllBands) ? VImageNBands (src) : 1;
  }
  va_end (args);
  return result;
}

VImage VSelectDestImage (VStringConst routine, VImage dest,
			 int nbands, int nrows, int ncolumns,
			 VRepnKind pixel_repn)
{
  /* If no destination image was specified, allocate one: */
  if (! dest)
    return VCreateImage (nbands, nrows, ncolumns, pixel_repn);

    /* Otherwise check that the destination provided has the appropriate
       characteristics: */
  if (VImageNBands (dest) != nbands) {
    VWarning ("%s: Destination image has %d bands; %d expected",
	      routine, VImageNBands (dest), nbands);
    return NULL;
  }
  if (VImageNRows (dest) != nrows) {
    VWarning ("%s: Destination image has %d rows; %d expected",
	      routine, VImageNRows (dest), nrows);
    return NULL;
  }
  if (VImageNColumns (dest) != ncolumns) {
    VWarning ("%s: Destination image has %d columns; %d expected",
	      routine, VImageNColumns (dest), ncolumns);
    return NULL;
  }
  if (VPixelRepn (dest) != pixel_repn) {
    VWarning ("%s: Destination image has %s pixels; %s expected", routine, 
	      VPixelRepnName (dest), VRepnName (pixel_repn));
    return NULL;
  }
  return dest;
}

VBoolean VSelectBand (VStringConst routine, VImage image, VBand band,
		      int *npixels, VPointer *first_pixel)
{
  if (band == VAllBands) {
    if (npixels)
      *npixels = VImageNPixels (image);
    if (first_pixel)
      *first_pixel = VImageData (image);
  } else if (band >= 0 && band < VImageNBands (image)) {
    if (npixels)
      *npixels = VImageNRows (image) * VImageNColumns (image);
    if (first_pixel)
      *first_pixel = image->band_index[band][0];
  } else {
    VWarning ("%s: Band %d referenced in image of %d band(s)",
	      routine, band, VImageNBands (image));
    return FALSE;
  }
  return TRUE;
}

VBandInterp VImageFrameInterp (VImage image)
{
  VLong interp;
  VGetAttrResult result;

  if (VImageNBands (image) !=
      (VImageNFrames (image) * VImageNViewpoints (image) *
       VImageNColors (image) * VImageNComponents (image)))
    VWarning ("VImageFrameInterp: No. bands (%d) conflicts with no. "
	      "of frames, etc. (%d %d %d %d)",
	      VImageNBands (image), VImageNFrames (image),
	      VImageNViewpoints (image), VImageNColors (image),
	      VImageNComponents (image));

  if (! VImageAttrList (image) ||
      (result =
       VGetAttr (VImageAttrList (image), VFrameInterpAttr,
		 VBandInterpDict, VLongRepn, & interp)) == VAttrMissing)
    return VImageNFrames (image) > 1 ? VBandInterpOther : VBandInterpNone;

  if (result == VAttrBadValue)
    return VBandInterpOther;

  switch (interp) {

  }
  return VBandInterpOther;
}

VBandInterp VImageViewpointInterp (VImage image)
{
  VLong interp;
  VGetAttrResult result;

  if (VImageNBands (image) !=
      (VImageNFrames (image) * VImageNViewpoints (image) *
       VImageNColors (image) * VImageNComponents (image)))
    VWarning ("VImageViewpointInterp: No. bands (%d) conflicts with no. "
	      "of frames, etc. (%d %d %d %d)",
	      VImageNBands (image), VImageNFrames (image),
	      VImageNViewpoints (image), VImageNColors (image),
	      VImageNComponents (image));
    
  if (! VImageAttrList (image) ||
      (result =
       VGetAttr (VImageAttrList (image), VViewpointInterpAttr,
		 VBandInterpDict, VLongRepn, & interp)) == VAttrMissing)
    return VImageNViewpoints (image) > 1 ?
      VBandInterpOther : VBandInterpNone;

  if (result == VAttrBadValue)
    return VBandInterpOther;

  switch (interp) {

  case VBandInterpStereoPair:
    if (VImageNViewpoints (image) != 2) {
      VWarning ("VBandViewpointInterp: "
		"Stereo-pair image has %d viewpoint dimension(s)",
		VImageNViewpoints (image));
      return VBandInterpOther;
    }
    return VBandInterpStereoPair;
  }
  return VBandInterpOther;
}

VBandInterp VImageColorInterp (VImage image)
{
  VLong interp;
  VGetAttrResult result;

  if (VImageNBands (image) !=
      (VImageNFrames (image) * VImageNViewpoints (image) *
       VImageNColors (image) * VImageNComponents (image)))
    VWarning ("VImageColorInterp: No. bands (%d) conflicts with no. "
	      "of frames, etc. (%d %d %d %d)",
	      VImageNBands (image), VImageNFrames (image),
	      VImageNViewpoints (image), VImageNColors (image),
	      VImageNComponents (image));
    
  if (! VImageAttrList (image) ||
      (result =
       VGetAttr (VImageAttrList (image), VColorInterpAttr,
		 VBandInterpDict, VLongRepn, & interp)) == VAttrMissing)
    return VImageNColors (image) > 1 ? VBandInterpOther : VBandInterpNone;

  if (result == VAttrBadValue)
    return VBandInterpOther;

  switch (interp) {

  case VBandInterpRGB:
    if (VImageNColors (image) != 3) {
      VWarning ("VBandColorInterp: RGB image has %d color dimension(s)",
		VImageNColors (image));
      return VBandInterpOther;
    }
    return VBandInterpRGB;
  }
  return VBandInterpOther;
}

VBandInterp VImageComponentInterp (VImage image)
{
  VLong interp;
  VGetAttrResult result;

  if (VImageNBands (image) !=
      (VImageNFrames (image) * VImageNViewpoints (image) *
       VImageNColors (image) * VImageNComponents (image)))
    VWarning ("VImageComponentInterp: No. bands (%d) conflicts with no. "
	      "of frames, etc. (%d %d %d %d)",
	      VImageNBands (image), VImageNFrames (image),
	      VImageNViewpoints (image), VImageNColors (image),
	      VImageNComponents (image));
    
  if (! VImageAttrList (image) ||
      (result =
       VGetAttr (VImageAttrList (image), VComponentInterpAttr,
		 VBandInterpDict, VLongRepn, & interp)) == VAttrMissing)
    return VImageNComponents (image) > 1 ?
      VBandInterpOther : VBandInterpNone;

  if (result == VAttrBadValue)
    return VBandInterpOther;

  switch (interp) {

  case VBandInterpComplex:
    if (VImageNComponents (image) != 2) {
      VWarning ("VBandColorInterp: Complex image has %d component(s)",
		VImageNComponents (image));
      return VBandInterpOther;
    }
    return VBandInterpComplex;

  case VBandInterpGradient:
    if (VImageNComponents (image) > 3) {
      VWarning ("VBandColorInterp: Gradient image has %d component(s)",
		VImageNComponents (image));
      return VBandInterpOther;
    }
    return VBandInterpGradient;

  case VBandInterpIntensity:
    if (VImageNComponents (image) > 1) {
      VWarning ("VBandColorInterp: Intensity image has %d component(s)",
		VImageNComponents (image));
      return VBandInterpOther;
    }
    return VBandInterpIntensity;

  case VBandInterpOrientation:
    if (VImageNComponents (image) > 1) {
      VWarning ("VBandColorInterp: "
		"Orientation image has %d component(s)",
		VImageNComponents (image));
      return VBandInterpOther;
    }
    return VBandInterpOrientation;
  }
  return VBandInterpOther;
}

VBoolean VSetBandInterp (VImage image,
			 VBandInterp frame_interp, int nframes,
			 VBandInterp viewpoint_interp, int nviewpoints,
			 VBandInterp color_interp, int ncolors,
			 VBandInterp component_interp, int ncomponents)
{
  VBoolean result = TRUE;
  VString str;

  if (VImageNBands (image) !=
      nframes * nviewpoints * ncolors * ncomponents) {
    VWarning ("VSetBandInterp: No. bands (%d) conflicts with no. "
	      "of frames, etc. (%d %d %d %d)", VImageNBands (image),
	      nframes, nviewpoints, ncolors, ncomponents);
    result = FALSE;
  }

  if (frame_interp == VBandInterpNone)
    result &= VExtractAttr (VImageAttrList (image), VFrameInterpAttr,
			    NULL, VStringRepn, & str, FALSE);
  else VSetAttr (VImageAttrList (image), VFrameInterpAttr,
		 VBandInterpDict, VLongRepn, (VLong) frame_interp);
  VImageNFrames (image) = nframes;

  if (viewpoint_interp == VBandInterpNone)
    result &= VExtractAttr (VImageAttrList (image), VViewpointInterpAttr,
			    NULL, VStringRepn, & str, FALSE);
  else VSetAttr (VImageAttrList (image), VViewpointInterpAttr,
		 VBandInterpDict, VLongRepn, (VLong) viewpoint_interp);
  VImageNViewpoints (image) = nviewpoints;

  if (color_interp == VBandInterpNone)
    result &= VExtractAttr (VImageAttrList (image), VColorInterpAttr,
			    NULL, VStringRepn, & str, FALSE);
  else VSetAttr (VImageAttrList (image), VColorInterpAttr,
		 VBandInterpDict, VLongRepn, (VLong) color_interp);
  VImageNColors (image) = ncolors;

  if (component_interp == VBandInterpNone)
    result &= VExtractAttr (VImageAttrList (image), VComponentInterpAttr,
			    NULL, VStringRepn, & str, FALSE);
  else VSetAttr (VImageAttrList (image), VComponentInterpAttr,
		 VBandInterpDict, VLongRepn, (VLong) component_interp);
  VImageNComponents (image) = ncomponents;

  return result;
}

int VReadImages (FILE *file, VAttrList *attributes, VImage **images)
{
  return VReadObjects (file, VImageRepn, attributes, (VPointer **) images);
}

VBoolean VWriteImages (FILE *file, VAttrList attributes,
		       int nimages, VImage images[])
{
  return VWriteObjects (file, VImageRepn, attributes, nimages,
			(VPointer *) images);
}

#define TallyStats(type)						\
    {									\
	type *pp = first_pixel;						\
									\
	tmin = tmax = *pp;						\
	for (i = 0; i < npixels; i++) {					\
	    tmean += pixel = *pp++;					\
	    tvar  += (pixel * pixel);					\
	    if (pixel < tmin)						\
		tmin = pixel;						\
	    else if (pixel > tmax)					\
		tmax = pixel;						\
	}								\
}

VBoolean VImageStats(VImage src, VBand band, VDouble *pmin, VDouble *pmax,
		      VDouble *pmean, VDouble *pvar)
{
    int npixels, i;
    VPointer first_pixel;
    VDouble pixel, tmin, tmax, tmean, tvar;

    /* Prepare to iterate over the specified band(s) of source pixels: */
    if (! VSelectBand ("VImageStats", src, band, & npixels, & first_pixel))
	return FALSE;

    /* Initialize accumulators: */
    tmin = tmax = tmean = tvar = 0.0;

    /* Tally pixels: */
    switch (VPixelRepn (src)) {
    case VBitRepn:	TallyStats (VBit); break;
    case VUByteRepn:	TallyStats (VUByte); break;
    case VSByteRepn:	TallyStats (VSByte); break;
    case VShortRepn:	TallyStats (VShort); break;
    case VLongRepn:	TallyStats (VLong); break;
    case VFloatRepn:	TallyStats (VFloat); break;
    case VDoubleRepn:	TallyStats (VDouble); break;
    default:		break;
    }

    /* Compute and return the final results: */
    tmean /= npixels;
    tvar = tvar / npixels - (tmean * tmean);
    if (pmin) *pmin = tmin;
    if (pmax) *pmax = tmax;
    if (pmean) *pmean = tmean;
    if (pvar) *pvar = tvar;
    return TRUE;
}

VImage VScaleIntensity(VImage src, double white, double black)
/* scale any input image to VUByte,
   mapping white (black) percent of the voxel to white (black) */
{
	int x, y, z, nx, ny, nz, i, range, maxshort;
	unsigned int lb, ub, limit, sum, *histo;
	double m, b, max, min, mean, var, v;
	VImage dst;
	
	maxshort = (int)(VRepnMaxValue(VShortRepn));
	histo = (unsigned int *)VCalloc(maxshort, sizeof(unsigned int));
	nx = VImageNColumns(src);
	ny = VImageNRows(src);	
	nz = VImageNBands(src);

	if (white < 0 || white > 100 || black < 0 || black > 100 || white+black >= 100)  {
		fprintf(stderr, "VScaleIntensity: illegal percentage given.\n");
		return 0;
	};
	
	/* first pass: find maximum and minimum values */
	VImageStats(src, VAllBands, &min, &max, &mean, &var);
	if (max == min) {
		fprintf(stderr, "VScaleIntensity: illegal data in image.\n");
		return 0;
	};
	b = min;
	m = (max-min) / (double)maxshort;

	/* second pass: build a histogram*/
	for (z = 0; z < nz; z++)  {
		for (y = 0; y < ny; y++)  {
			for (x = 0; x < nx; x++)  {
				v = VGetPixel(src, z, y, x);
				i = (int)((v-b)/m+0.5);
				histo[i]++;
			};
		};
	};

	/* throw away pc percent of the voxel below lb and pc percent above ub */
	limit = (black * nx * ny * nz) / 100;
        lb = 0; sum = 0;
        for (i = 0; i < maxshort; i++)  {
        	sum += histo[i];
        	if (sum >= limit) { lb = i; break; };
        };
	limit = (white * nx * ny * nz) / 100;
        ub = maxshort-1; sum = 0;
        for (i = maxshort-1; i >= 0; i--)  {
        	sum += histo[i];
        	if (sum >= limit) { ub = i; break; };
        };
	min = lb*m+b;
	max = ub*m+b;

	/* third pass: create and convert image */
	dst = VCreateImage(nz, ny, nx, VUByteRepn);
	if (dst == 0) return 0;
	
	range = 256;
        m = range / (max - min);
        b = range - (m * max);
	for (z = 0; z < nz; z++)  {
		for (y = 0; y < ny; y++)  {
			for (x = 0; x < nx; x++)  {
				v = VGetPixel(src, z, y, x);
				i = (int)(v * m + b + 0.5);
                        	if (i < 0) i = 0;
                        	else if (i >= range) i = range-1;
                        	VPixel(dst, z, y, x, VUByte) = i;
                	};
                };
        };
        VFree(histo);
	VCopyImageAttrs(src, dst);
        return dst;
}

VBoolean VFillImage (VImage image, VBand band, VDoublePromoted value)
{
  int i, npixels;
  VPointer first_pixel;

#define Fill(type)					\
	{						\
	    type d = value, *pp = first_pixel;		\
	    for (i = npixels; i > 0; i--)		\
		*pp++ = d;				\
	}

  /* Locate the specified band(s) in the image: */
  if (! VSelectBand ("VFillImage", image, band, & npixels, & first_pixel))
    return FALSE;

  /* To avoid surprises when filling integer pixels, round the fill
     value to the nearest integer: */
  if (VIsIntegerRepn (VPixelRepn (image)))
    value += (value > 0.0) ? 0.5 : -0.5;

  /* The memset() routine is probably the fastest way of filling memory, but
     it can only be used here in certain situations. We only check for
     these, the most straightforward of the situations:
     (a) pixel values occupy a byte
     (b) the value to be filled is all 0's
     It is when the value to be filled is all 0's and the pixel
     representation is floating point that we take advantage of the
     assumption that 0.0 is represented as all-bits-zero. */
  if (VPixelSize (image) == 1 || value == 0.0)
    memset (first_pixel, (int) value, npixels * VPixelSize (image));

  /* Otherwise, fill by looping over all pixels: */
  else
    switch (VPixelRepn (image)) {
    case VBitRepn:	Fill (VBit);	break;
    case VUByteRepn:	Fill (VUByte);  break;
    case VSByteRepn:	Fill (VSByte);  break;
    case VShortRepn:	Fill (VShort);  break;
    case VLongRepn:	Fill (VLong);   break;
    case VFloatRepn:	Fill (VFloat);  break;
    case VDoubleRepn:	Fill (VDouble); break;
    default: break;
    }

  return TRUE;
}
	
