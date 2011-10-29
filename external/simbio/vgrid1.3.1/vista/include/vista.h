#ifndef VISTA_H
#define VISTA_H 1

/*
 *  Copyright 1994 University of British Columbia
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

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <unistd.h>
#include <limits.h>
#include <stdarg.h>

#ifndef FALSE
#define FALSE 			0
#define TRUE 			1
#endif

#define VMax(a,b)		((a) > (b) ? (a) : (b))
#define VMin(a,b)		((a) < (b) ? (a) : (b))
#define VOffset(type, field) 	((size_t) (((char *) & ((type) 0)->field) - (char *) 0))
#define VOffsetOf(type, field)	VOffset(type *, field)
#define VNumber(array)		((size_t) (sizeof (array) / sizeof ((array)[0])))
#define VNew(type)		((type *) VMalloc (sizeof (type)))
#define VNewString(str) 	((VString) ((str) ? strcpy ((char *) VMalloc (strlen (str) + 1), str) : 0))
#define VFileHeader		"V-data"
#define VFileVersion		2
#define VFileDelimiter		"\f\n"
#define VMaxAttrNameLength	256
#define VRequiredOpt		(& V_RequiredOpt)
#define VOptionalOpt		(& V_OptionalOpt)
#define VAllBands		-1		/* all bands */
#define VAttrListEmpty(l)	((l) == NULL || (l)->next == NULL)
#define VFirstAttr(l,p)		((void) ((p)->list = (l), (p)->ptr = (l)->next))
#define VLastAttr(l,p)		((void) ((p)->list = (l), (p)->ptr = (l)->prev))
#define VAttrExists(p)		((p)->ptr != NULL)
#define VNextAttr(p)		((void) ((p)->ptr = (p)->ptr ? (p)->ptr->next : NULL))
#define VPrevAttr(p)		((void) ((p)->ptr = (p)->ptr ? (p)->ptr->prev : NULL))
#define VGetAttrName(p)		((p)->ptr->name)
#define VGetAttrRepn(p)		((p)->ptr->repn)
#define VListCount(vlist)	((vlist)->count)
#define VListCurr(vlist)	((vlist)->current->item)
#define VListGetCurr(vlist)	((vlist)->current)
#define VListSetCurr(vlist,cur) ((void)((vlist)->current = (cur)))
#define VImageNBands(image)	((image)->nbands)
#define VImageNRows(image)	((image)->nrows)
#define VImageNColumns(image)	((image)->ncolumns)
#define VImageNFrames(image)	((image)->nframes)
#define VImageNViewpoints(image) ((image)->nviewpoints)
#define VImageNColors(image)	((image)->ncolors)
#define VImageNComponents(image) ((image)->ncomponents)
#define VPixelRepn(image)	((image)->pixel_repn)
#define VImageData(image)	((image)->data)
#define VImageAttrList(image)	((image)->attributes)
#define VImageNPixels(image) 	((image)->nbands * (image)->nrows * (image)->ncolumns)
#define VPixelSize(image)	(VRepnSize ((image)->pixel_repn))
#define VPixelPrecision(image)  (VRepnPrecision ((image)->pixel_repn))
#define VPixelRepnName(image)	(VRepnName ((image)->pixel_repn))
#define VPixelMinValue(image)	(VRepnMinValue ((image)->pixel_repn))
#define VPixelMaxValue(image)	(VRepnMaxValue ((image)->pixel_repn))
#define VImageSize(image)	(VImageNPixels(image) * VPixelSize(image))
#define VPixelPtr(image, band, row, column) \
				((VPointer) ((char *) ((image)->band_index[band][row])+(column) * VPixelSize (image)))
#define VPixel(image, band, row, column, type) \
				(* ((type *) (image)->band_index[band][row] + (column)))
#define VPixelArray(image, type) ((type ***) (image)->band_index)
#define VBandIndex(image, frame, viewpoint, color, component) \
				(((((frame) * (image)->nviewpoints + (viewpoint)) * (image)->ncolors + \
					(color)) * (image)->ncomponents) + (component))
#define VSameImageRange(image1, image2)	\
				((image1)->nbands == (image2)->nbands && (image1)->nrows == (image2)->nrows && \
				(image1)->ncolumns == (image2)->ncolumns && (image1)->pixel_repn == (image2)->pixel_repn)
#define VSameImageSize(image1, image2) \
				((image1)->nbands == (image2)->nbands && (image1)->nrows == (image2)->nrows && \
				(image1)->ncolumns == (image2)->ncolumns)
#define VGraphNNodes(graph)	(graph->nnodes)
#define VGraphNFields(graph)	(graph->nfields)
#define VGraphNSize(graph)	(graph->size)
#define VGraphAttrList(graph)	(graph->attributes)
#define VGraphGetNode(graph, nid)	(graph->table[nid-1])
#define VGraphNodeIsFree(graph, nid)	(graph->table[nid-1] == 0)
#define VNodeRepn(graph)	(graph->node_repn)
#define VNodeSize(graph) 	(sizeof(VNodeBaseRec) + (graph->nfields * VRepnPrecision(graph->node_repn)) / 8)
#define VNodeTestVisit(node)	(((VNodeBase)node)->visited == TRUE)
#define VNodeSetVisit(node)	(((VNodeBase)node)->visited = TRUE)
#define VNodeClearVisit(node)	(((VNodeBase)node)->visited = FALSE)


/* Names of generic attributes: */
#define VCommentAttr		"comment"
#define VDataAttr		"data"
#define VHistoryAttr		"history"
#define VLengthAttr		"length"
#define VNameAttr		"name"
#define VNColumnsAttr		"ncolumns"
#define VNRowsAttr		"nrows"
#define VRepnAttr		"repn"

/* Image attribute type names: */
#define VColorInterpAttr	"color_interp"
#define VComponentInterpAttr	"component_interp"
#define VFrameInterpAttr	"frame_interp"
#define VNBandsAttr		"nbands"
#define VNColorsAttr		"ncolors"
#define VNComponentsAttr	"ncomponents"
#define VNFramesAttr		"nframes"
#define VNViewpointsAttr	"nviewpoints"
#define VPixelAspectRatioAttr	"pixel_aspect_ratio"
#define VViewpointInterpAttr	"viewpoint_interp"

/* Graph attribute type names: */
#define VGraphAttr		"Graph"
#define VNGraphNodesAttr	"nnodes"
#define VNGraphSizeAttr		"size"
#define VNNodeFieldsAttr	"nfields"
#define VNNodeWeightsAttr	"useWeights"

/* Macros for generating constants of particular numeric types: */
/* (These definitions may be platform-specific.) */
#define VBitConst(c)		(c)
#define VUByteConst(c)		(c)
#define VSByteConst(c)		(c)
#define VShortConst(c)		(c)
#define VLongConst(c)		(c ## l)
#define VFloatConst(c)		(c ## f)
#define VDoubleConst(c) 	(c)

/* (These definitions may be platform-specific.) */
typedef char VBit;			/* 0 or 1 */
typedef double VDouble;			/* >= 64-bit IEEE floating point */
typedef float VFloat;			/* >= 32-bit IEEE floating point */
typedef int VLong;			/* !! changed, G.L. 19.9.95 !! */
typedef signed char VSByte;		/* integer in [-128,127] */
typedef short VShort;			/* >= 16-bit signed integer */
typedef unsigned char VUByte;		/* integer in [0,255] */
typedef char VBoolean;			/* TRUE or FALSE */
typedef void *VPointer;			/* generic pointer */
typedef const char *VStringConst;	/* null-terminated string constant */
typedef char *VString;			/* null-terminated string */
typedef int VBitPromoted;
typedef int VBooleanPromoted;
typedef double VDoublePromoted;
typedef double VFloatPromoted;
typedef long VLongPromoted;
typedef int VSBytePromoted;
typedef int VShortPromoted;
typedef unsigned int VUBytePromoted;
typedef struct V_ImageRec *VImage;
typedef int VBand;
typedef void VErrorHandler(VStringConst);
typedef void VWarningHandler(VStringConst);
typedef VPointer VCopyMethod(VPointer);
typedef void VDestroyMethod(VPointer);

extern VBoolean V_RequiredOpt, V_OptionalOpt;

/* Codes for referring to representations: */
typedef enum {
    VUnknownRepn,
    VBitRepn,				/* 1-bit integer, [0, 1] */
    VUByteRepn,				/* 8-bit integer, [0, 255] */
    VSByteRepn,				/* 8-bit integer, [-128, 127] */
    VShortRepn,				/* 16-bit integer, [-32768, 32767] */
    VLongRepn,				/* 32-bit integer, [-2**31, 2**31-1] */
    VFloatRepn,				/* 32-bit IEEE floating point */
    VDoubleRepn,			/* 64-bit IEEE floating point */
    VAttrListRepn,			/* attribute list */
    VBooleanRepn,			/* TRUE or FALSE */
    VBundleRepn,			/* object of named type */
    VListRepn,				/* list of opaque objects */
    VPointerRepn,			/* pointer to opaque object */
    VStringRepn,			/* null-terminated string */
    VImageRepn,				/* image */
    VGraphRepn,                         /* graph */
    VNRepnKinds				/* number of predefined types */
} VRepnKind;

/* Values of band interpretation attributes: */
typedef enum {
    VBandInterpNone,			/* no interpretation specified */
    VBandInterpOther,			/* unknown interpretation specified */
    VBandInterpStereoPair,
    VBandInterpRGB,
    VBandInterpComplex,
    VBandInterpGradient,
    VBandInterpIntensity,
    VBandInterpOrientation
} VBandInterp;

/* Dictionary entry: */
typedef struct {
    /* The following are initialized by the dictionary provider: */
    VStringConst keyword;		/* keyword string */
    VLong ivalue;			/* value, if an integer */
    VStringConst svalue;		/* value, if a string */

    /* The following are used only by code in VLookupDictValue: */
    VBoolean icached;			/* whether integer value cached */
    VBoolean fcached;			/* whether float value cached */
    VDouble fvalue;			/* cached floating-point value */
} VDictEntry;

/* Accepted command options are described by a table of these entries: */
typedef struct {
    VStringConst keyword;		/* keyword signalling option */
    VRepnKind repn;			/* type of value supplied by option */
    int number;				/* number of values supplied */
    VPointer value;			/* location for storing value(s) */
    VBoolean *found;			/* whether optional arg found */
    VDictEntry *dict;			/* optional dict of value keywords */
    VStringConst blurb;			/* on-line help blurb */
} VOptionDescRec;

/* If an option takes multiple values, they are represented by a VArgVector: */
typedef struct {
    int number;				/* number of arguments */
    VPointer vector;			/* vector of arguments */
} VArgVector;

/* Each attribute name/value pair is represented by: */
typedef struct V_AttrRec {
    struct V_AttrRec *next;		/* next in list */
    struct V_AttrRec *prev;		/* previous in list */
    VRepnKind repn;			/* rep'n of attribute value */
    VPointer value;			/* pointer to attribute value */
    char name[1];			/* beginning of name string */
} VAttrRec;

typedef VAttrRec *VAttrList;

/* Position within a list of attributes: */
typedef struct {
    VAttrList list;			/* the list */
    struct V_AttrRec *ptr;		/* position within the list */
} VAttrListPosn;

/* Result of trying to retrieve an attribute's value: */
typedef enum {
    VAttrFound,				/* successfully retrieved value */
    VAttrMissing,			/* didn't find attribute */
    VAttrBadValue			/* incompatible value */
} VGetAttrResult;

/* An object whose type is named but not registered: */
typedef struct {
    VAttrList list;			/* object's attribute list value */
    size_t length;			/* length of binary data */
    VPointer data;			/* pointer to binary data */
    char type_name[1];			/* beginning of object's type's name */
} VBundleRec, *VBundle;

typedef VPointer VDecodeMethod(VStringConst, VBundle);
typedef VAttrList VEncodeAttrMethod(VPointer, size_t *);
typedef VPointer VEncodeDataMethod(VPointer, VAttrList, size_t, VBoolean *);

/* Set of methods supporting an object type: */
typedef struct {
    VCopyMethod *copy;
    VDestroyMethod *destroy;
    VDecodeMethod *decode;
    VEncodeAttrMethod *encode_attr;
    VEncodeDataMethod *encode_data;
} VTypeMethods;

/* Information about a representation: */
typedef struct {
    VStringConst name;			/* name string */
    size_t size;			/* size, in bytes */
    int precision;			/* precision, in bits */
    VDouble min_value;			/* min and max representable values */
    VDouble max_value;
    VTypeMethods *methods;		/* associated methods */
} VRepnInfoRec;

/* List element: */
typedef struct V_Node *VNodePtrType;
struct V_Node {
    VPointer item;			/* pointer to data item */
    VNodePtrType prev;			/* pointer to previous node */
    VNodePtrType next;			/* pointer to next node */
};

/* List head: */
typedef struct V_List {
    VNodePtrType current;		/* pointer to current node */
    VNodePtrType head;			/* pointer to head node */
    VNodePtrType tail;			/* pointer to tail node */
    int count;				/* number of nodes in VList */
} *VList;

/* Description of an image: */
typedef struct V_ImageRec {
    int nbands;				/* number of bands */
    int nrows;				/* number of rows */
    int ncolumns;			/* number of columns */
    VRepnKind pixel_repn;		/* representation of pixel values */
    unsigned long flags;		/* various flags */
    VAttrList attributes;		/* list of other image attributes */
    VPointer data;			/* array of image pixel values */
    VPointer *row_index;		/* ptr to first pixel of each row */
    VPointer **band_index;		/* ptr to first row of each band */
    int nframes;			/* number of motion frames */
    int nviewpoints;			/* number of camera viewpoints */
    int ncolors;			/* number of color channels */
    int ncomponents;			/* number of vector components */
} VImageRec;

/* Codes for flags: */
enum {
    VImageSingleAlloc = 0x01		/* one free() releases everything */
};

typedef struct V_GraphRec {
    int nnodes;				/* number of nodes */
    int nfields;			/* size of fields in a node´s private area */
    VRepnKind node_repn;		/* data representation in a node */
    VAttrList attributes;		/* list of other attributes */
    struct VNodestruct **table;		/* node table of Graph */
    int size;				/* number of places in table */
    int lastUsed;			/* last entry used in table */
    int iter;				/* iteration counter in sequential access */
    int useWeights;			/* TRUE iff weights are used */
} VGraphRec, *VGraph;

typedef struct VNodebaseStruct {
    unsigned int hops: 31;		/* numbor of hops in this node */
    unsigned int visited: 1;		/* true iff seen before */
    VFloat weight;			/* weight of this node */
    struct VAdjstruct *head;
} VNodeBaseRec, *VNodeBase;

typedef struct VNodestruct {
    VNodeBaseRec base;
    char data[1];			/* private data area of node starts here */
} VNodeRec, *VNode;

typedef struct VAdjstruct {
    unsigned int id;			/* node reference */
    VFloat weight;			/* weight of this node */
    struct VAdjstruct *next;		/* list of adjacent nodes */
} VAdjRec, *VAdjacency;

extern VRepnInfoRec *VRepnInfo;
extern VDictEntry VBooleanDict[];	/* boolean values */
extern VDictEntry VNumericRepnDict[];	/* numeric representation kinds */
extern VDictEntry VBandInterpDict[];

/* A list of attributes is represented by a header node: */
typedef enum { VLsbFirst, VMsbFirst } VPackOrder;
typedef VBoolean VReadFileFilterProc(VBundle, VRepnKind);

/* Macros for accessing information about representations: */
#define VRepnSize(repn)		(VRepnInfo[repn].size)
#define VRepnPrecision(repn)	(VRepnInfo[repn].precision)
#define VRepnName(repn)		(VRepnInfo[repn].name)
#define VRepnMinValue(repn)	(VRepnInfo[repn].min_value)
#define VRepnMaxValue(repn)	(VRepnInfo[repn].max_value)
#define VRepnMethods(repn)	(VRepnInfo[repn].methods)
#define VIsIntegerRepn(repn)	((repn) >= VBitRepn && (repn) <= VLongRepn)
#define VIsFloatPtRepn(repn)	((repn) == VFloatRepn || (repn) == VDoubleRepn)

/*  Declarations of library routines. */
#ifdef __cplusplus
extern "C" {
#endif
extern VGraph VCreateGraph(int, int, VRepnKind, int);
extern VGraph VCopyGraph(VGraph);
extern void VDestroyGraph(VGraph);
extern int VReadGraphs(FILE *, VAttrList *, VGraph **);
extern VBoolean VWriteGraphs(FILE *, VAttrList, int, VGraph *);
extern int VGraphLookupNode(VGraph, VNode);
extern int VGraphAddNode(VGraph, VNode);
extern int VGraphAddNodeAt(VGraph, VNode, int);
extern int VGraphLinkNodes(VGraph, int, int);
extern int VGraphUnlinkNodes(VGraph, int, int);
extern VPointer VGraphFirstNode(VGraph);
extern VPointer VGraphNextNode(VGraph);
extern void VGraphClearVisit(VGraph);
extern int VGraphResizeFields(VGraph, int);
extern int VGraphNCycles(VGraph);
extern void VGraphToggleNodesFrom(VGraph, int);
extern void VDestroyNode(VGraph, int);
extern void VGraphDestroyNodesFrom(VGraph, int);
extern void VGraphClearHops(VGraph);
extern int VGraphAddAndGrow(VGraph graph, VNode node, int pos);
extern void VNodeRemoveLinks(VGraph graph, int pos);
extern void linkNodes(VGraph graph, VLong a, VLong b);
extern int hasLink (VGraph graph, int a, int b);
extern void unlinkNodes(VGraph graph, VLong a, VLong b);
extern VBoolean VIdentifyFiles(int, VOptionDescRec[], VStringConst, int *, char **, int);
extern VBoolean VParseCommand(int, VOptionDescRec[], int *, char **);
extern void VParseFilterCmd(int, VOptionDescRec[], int, char **, FILE **, FILE **);
extern void VPrintOptions(FILE *, int, VOptionDescRec[]);
extern int VPrintOptionValue(FILE *, VOptionDescRec*);
extern void VReportBadArgs(int, char **);
extern void VReportUsage(VStringConst, int, VOptionDescRec[], VStringConst);
extern void VReportValidOptions(int, VOptionDescRec[]);
extern VBoolean VLoadParameters(int, VOptionDescRec[], VStringConst, VStringConst, VPointer, VBoolean);
extern VBoolean VParseParamDefn(VStringConst, VString, VRepnKind *, VString, VString);
extern VBoolean VParseParamOptions(int, VOptionDescRec[], int *, char **, VPointer);
extern void VPrintParameters(FILE *, int, VOptionDescRec[], VPointer);
extern void VReportValidParamOptions(int, VOptionDescRec[], VPointer);
extern VImage VCreateImage(int, int, int, VRepnKind);
extern VImage VCreateImageLike(VImage);
extern void VDestroyImage(VImage);
extern VDouble VGetPixel(VImage, int, int, int);
extern void VSetPixel(VImage, int, int, int, VDoublePromoted);
extern VImage VCopyImage(VImage, VImage, VBand);
extern VImage VCopyImageAttrs(VImage, VImage);
extern VImage VCopyImagePixels(VImage, VImage, VBand);
extern VBoolean VCopyBand (VImage, VBand, VImage, VBand);
extern VImage VCombineBands(int, VImage[], VBand[], VImage);
extern VImage VCombineBandsVa(VImage, ...);
extern VImage VSelectDestImage(VStringConst, VImage, int, int, int, VRepnKind);
extern VBoolean VSelectBand(VStringConst, VImage, VBand, int *, VPointer *);
extern VBandInterp VImageFrameInterp(VImage);
extern VBandInterp VImageViewpointInterp(VImage);
extern VBandInterp VImageColorInterp(VImage);
extern VBandInterp VImageComponentInterp(VImage);
extern VBoolean VSetBandInterp(VImage, VBandInterp, int, VBandInterp,
		int, VBandInterp, int, VBandInterp, int);
extern int VReadImages(FILE *, VAttrList *, VImage **);
extern VBoolean VWriteImages(FILE *, VAttrList, int, VImage[]);
extern VBoolean VImageStats(VImage src, VBand band, VDouble *pmin,
		VDouble *pmax, VDouble *pmean, VDouble *pvar);
extern FILE *VOpenInputFile(VStringConst, VBoolean);
extern FILE *VOpenOutputFile(VStringConst, VBoolean);
extern int VReadObjects(FILE *, VRepnKind, VAttrList *, VPointer **);
extern VAttrList VReadFile(FILE *, VReadFileFilterProc *);
extern VBoolean VWriteObjects(FILE *, VRepnKind, VAttrList, int, VPointer[]);
extern VBoolean VWriteFile(FILE *, VAttrList);
extern VList VListCreate();
extern VPointer VListFirst(VList);
extern VPointer VListLast(VList);
extern VPointer VListNext(VList);
extern VPointer VListPrev(VList);
extern void VListAdd(VList, VPointer);
extern void VListInsert(VList, VPointer);
extern void VListAppend(VList, VPointer);
extern void VListPrepend(VList, VPointer);
extern VPointer VListRemove(VList);
extern void VListConcat(VList, VList);
extern void VListDestroy(VList, void (*) (VPointer));
extern VPointer VListTrim(VList);
extern VPointer VListSearch(VList, int (*) (), VPointer);
extern VPointer VCalloc(size_t, size_t);
extern void VFree(VPointer);
extern VPointer VMalloc(size_t);
extern VPointer VRealloc(VPointer, size_t);
extern void VAppendAttr(VAttrList, VStringConst, VDictEntry *, VRepnKind, ...);
extern VAttrList VCopyAttrList(VAttrList);
extern VAttrList VCreateAttrList();
extern VBundle VCreateBundle(VStringConst, VAttrList, size_t, VPointer);
extern VBoolean VDecodeAttrValue(VStringConst, VDictEntry *, VRepnKind, VPointer);
extern void VDeleteAttr(VAttrListPosn *);
extern void VDestroyAttrList(VAttrList);
extern void VDestroyBundle(VBundle);
extern VStringConst VEncodeAttrValue(VDictEntry *, VRepnKind, ...);
extern VBoolean VExtractAttr(VAttrList, VStringConst, VDictEntry *, 
		VRepnKind, VPointer, VBooleanPromoted);
extern VGetAttrResult VGetAttr(VAttrList, VStringConst, VDictEntry *, 
		VRepnKind, VPointer);
extern VBoolean VGetAttrValue(VAttrListPosn *, VDictEntry *, VRepnKind, VPointer);
extern void VInsertAttr(VAttrListPosn *, VBooleanPromoted, VStringConst, 
		VDictEntry *, VRepnKind, ...);
extern VBoolean VLookupAttr(VAttrList, VStringConst, VAttrListPosn *);
extern void VPrependAttr(VAttrList, VStringConst, VDictEntry *, VRepnKind, ...);
extern void VSetAttr(VAttrList, VStringConst, VDictEntry *, VRepnKind, ...);
extern void VSetAttrValue(VAttrListPosn *, VDictEntry *, VRepnKind, ...);
extern VDictEntry *VLookupDictKeyword(VDictEntry *, VStringConst);
extern VDictEntry *VLookupDictValue (VDictEntry *, VRepnKind, ...);
extern void VError(VStringConst, ...);
extern void VWarning(VStringConst, ...);
extern void VSystemError(VStringConst, ...);
extern void VSystemWarning(VStringConst, ...);
extern VBoolean VPackData(VRepnKind, size_t, VPointer, VPackOrder, 
		size_t *, VPointer *, VBoolean *);
extern VBoolean VUnpackData(VRepnKind, size_t, VPointer, VPackOrder, 
		size_t *, VPointer *, VBoolean *);
extern void VPackBits(size_t, VPackOrder, VBit *, char *);
extern void VUnpackBits(size_t, VPackOrder, char *, VBit *);
extern VRepnKind VRegisterType(VStringConst, VTypeMethods *);
extern VRepnKind VLookupType(VStringConst);
extern VImage VScaleIntensity(VImage src, double white, double black);
extern VBoolean VFillImage (VImage image, VBand band, VDoublePromoted value);
extern VImage VConvertImageRange (VImage, VImage, VBand, VRepnKind);
extern VImage VBinarizeImage (VImage src, VImage dest, VBand band, double xmin,
			      double xmax);

	
#ifdef __cplusplus
}
#endif

#endif /* VISTA_H */
