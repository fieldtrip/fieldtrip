/*------------------------------------------------------------------------------
 *
 * Header file describing the numerical values used in fif files.
 * 
 * Copyright (c) 2002-2007 Elekta Neuromag Oy.
 * All rights reserved.
 *
 * No part of this program may be photocopied, reproduced,
 * or translated to another program language without the
 * prior written consent of the author.
 *
 * $Header: /sw/tmp/build/libs-fiff-2.1.0/src/RCS/fiff_file.h,v 1.19 2010/04/16 12:32:16 mjk Exp $
 *
 * $Id$
 */

#ifndef _fiff_file_h
#define _fiff_file_h

#define FIFFC_MAJOR_VERSION 1L
#define FIFFC_MINOR_VERSION 2L

#define FIFFC_VERSION (FIFFC_MAJOR_VERSION<<16 | FIFFC_MINOR_VERSION)

/*
 * Tag numbers < 100 are only used during the acquisition
 * they will never appear in a fif file.
 */

/*
 *  Conventions with FIFF_ C macros (mjk 14.01.2000).
 *===========================================================================
 * In order to get some order into the FIFF_ macros, 
 * following macro types are proposed:
 *
 * FIFF_      Fiff tag object identification label.
 * FIFFB_     Fiff block tag value (block type indentification label).
 * FIFFV_     Enumerated value of data having definite meaning in some context.
 *            Would be nice if the context is also present like: 
 *            FIFFV_HPI_ACCEPT_PROGRAM and FIFFV_HPI_ACCEPT_USER.
 * FIFFC_     Independent constant value like FIFFC_MATRIX_MAX_DIM used in
 *            code (having meaning in 'programming context' but not as data).
 * FIFFT_     Fiff type descriptor.
 * FIFFTS_    Fiff type decriptor (FIFFT_) structure definition (see below).
 *
 *
 *
 *  Conventions on type codes: (mjk 14.01.2000)
 *===========================================================================
 * Fiff types are saved using 32 bit numeric indentifiers.
 * The fiff type codes are structured so that they contain two main parts:
 * the 'fundamental structure' and 'type details'. The fundamental structure
 * is coded in the MSB of the 4 byte code. Depending on this code the 
 * interpretation of 'type details' may vary.
 *
 * Current Fundamental structures:
 *----------------------------------------------------------------------
 * Only the MSB is significant. See FIFFTS_FS_MASK
 *
 *    FIFFFS_SCALAR  0x00000000   Scalar type / basic fixed size FIFF record.
 *    FIFFFS_RECORD  0x10000000   <reserved>
 *    FIFFFS_......  0x20000000   <reserved>
 *    FIFFFS_......  0x30000000   <reserved>
 *    FIFFFS_MATRIX  0x40000000   Multidimensional matrix.
 *
 * The lower four bits are reserved for future extensions.
 *
 * Scalar types (FS==0x00):
 *----------------------------------------------------------------------
 * These include the basic scalar types and 'standard' fixed size records
 * used in fiff files. 
 *
 * * It is required that the code is less than 0x0FFF (< 4096). !!!
 * * Fourth byte is currently reserved for future extensions.
 *
 * Current types in this class are:
 *
 *   FIFFT_VOID                  0       Nothing
 *   FIFFT_BYTE                  1       Unsigned? 8 bits
 *   FIFFT_SHORT                 2       Signed 16 bit integer.
 *   FIFFT_INT                   3       Signed 32 bit integer.
 *   FIFFT_FLOAT                 4       Single precision IEEE float (32 bits)
 *   FIFFT_DOUBLE                5       Double precision IEEE float (64 bits)
 *   FIFFT_JULIAN                6	 Julian day. (32 bits)
 *   FIFFT_USHORT                7       Unsigned short (16 bits)
 *   FIFFT_UINT                  8       Unsigned int (32 bits)
 *   FIFFT_ULONG                 9       Unsigned long (64 bits)
 *   FIFFT_STRING               10       Octet, ASCII coding.
 *   FIFFT_ASCII                10
 *   FIFFT_LONG                 11       Long integer (64 bit)
 *   FIFFT_DAU_PACK13           13       13 bit packed format used in HP DAUs.
 *   FIFFT_DAU_PACK14           14       14 bit packed format used in HP DAUs
 *   FIFFT_DAU_PACK16           16       Signed 16 bit integer. (?)
 *   FIFFT_COMPLEX_FLOAT        20       Complex number encoded with floats
 *   FIFFT_COMPLEX_DOUBLE       21       Complex number encoded with doubles
 *   FIFFT_OLD_PACK             23       Neuromag proprietary 16 bit packing.
 *
 * Following are structure types defined in fiff_types.h
 *
 *   FIFFT_CH_INFO_STRUCT        30      Basic info about a measurement chnnel.
 *   FIFFT_ID_STRUCT             31      Unique identifier.
 *   FIFFT_DIR_ENTRY_STRUCT      32      FIFF file directory entry.
 *   FIFFT_DIG_POINT_STRUCT      33      Digitization point.
 *   FIFFT_CH_POS_STRUCT         34      Channel position.
 *   FIFFT_COORD_TRANS_STRUCT    35      Coordinate transformation.
 *   FIFFT_DIG_STRING_STRUCT     36      Digitization string.
 *   FIFFT_STREAM_SEGMENT_STRUCT 37      Data stream segment.
 * 
 * Simple vector of any scalar type can be formed by simply 
 * concatenating into file. The number of of elements is deduced
 * from the size of the data block which should be a multiple of the
 * object size.
 *
 * Futher more, for each scalar type fiff_types.h defines a corresponding
 * C language type. Naming convention is fiff_xxx_t which correspons to 
 * type id FIFFT_xxx.
 * 
 *
 *
 * Matrix types (FS=0x40)
 *----------------------------------------------------------------------
 * 
 * FIFF Matrix type is in principle an arbitrary dimensional rectangular 
 * collection of fixed size elements. However for practical reasons the
 * maximum dimension is currenly restricted arbitrarily to 9
 * (see FIFFC_MATRIX_MAX_DIM) but this restriction may be relaxed in the future.
 *
 * Matrix type codes use following structure:
 *
 *  * Matrix type structure:   0xFFCCyyyy
 *  * Where                     
 *
 *      0xFF......    decribes the 'fundamental structure'.
 *                       FF = 40 denotes a multidimensional matrix.
 *      0x..XX....    describes to basic coding: dense, triangular, sparse etc.
 *      0x....yyyy    describes the element type.
 *
 * Current codings available are:
 *
 *  FIFFTS_MC_DENSE  0x00000000    Dense column oriented matrix
 *  FIFFTS_MC_CCS    0x00100000    Column compressed sparse matrix
 *  FIFFTS_MC_RCS    0x00200000    Row compressed sparse matrix
 * 
 *
 * Dense matrix structure 
 * ----------------------
 * 
 * A(1,1,...),A(2,1...) ... A(1,2,..) ... A(N,M,K...),K,M,N, ...,DIM
 *
 * where DIM is the dimensionality.
 *       N,M,K... are the dimensions.
 *       A(i,j,k...) are the elements.
 * 
 * Note: the 2-dimensional meg_matrix format (routines in fiff_matrix.c) 
 * read and write the matrix in  transposed form, i.e, the storing order 
 * in the FIFF-file is
 * A(1,1),A(1,2) ... A(1,M) ... A(M,M),N,M,DIM
 *
 *
 * Column compressed sparse matrix structure
 *------------------------------------------
 *
 * A(x0,y0,z0...),A(x1,y0,z0),...
 *   ...,x0,x2...,x{NZ-1}, y0,y1,...,y{M-1}, z0,z1...,z{K-1},...,NZ,K,M,N,...,DIM
 *
 * where DIM is the dimensionality.
 *       N,M,K... are the dimensions.
 *       NZ is the number of non zero elements.
 *       A(i,j,k...) are the elements.
 *       x0..x(NZ-1) is row index array (concatenated for all vectors)
 *       y0..y(M-1)  is column start index array for the second dimension.
 *       z1..z(K-1) is slice start index array for the third dimension.
 *         etc.
 * Index arrays are 0 based.
 *
 * Row compressed sparse matrix structure
 *---------------------------------------
 *
 * Similar to column compresses version except that dimensions 1 and 2 are
 * interchanged. Structurally exactly same as ccs_matrix of the transpose.
 *
 */

/* Following definitions are needed only in programs that need to do
 * some 'intelligent' operations depending on arbitrary types.
 * They are rarely needed in user level code.
 *
 *  FIFFTS_FS_MASK     'fundamental structure' bit mask.
 *  FIFFTS_BASE_MASK   'Scalar value' (base value) bit mask.
 *  FIFFTS_MC_MASK     'Basic matrix coding' bit mask.
 *
 * Using type code structure constants directly is depreciated. Use 
 * functions to test type properties.
 */

/*
 * Constants for types
 */

#define FIFFT_VOID                  0
#define FIFFT_BYTE                  1
#define FIFFT_SHORT                 2
#define FIFFT_INT                   3
#define FIFFT_FLOAT                 4
#define FIFFT_DOUBLE                5
#define FIFFT_JULIAN                6	/* Julian day */
#define FIFFT_USHORT                7
#define FIFFT_UINT                  8
#define FIFFT_ULONG                 9
#define FIFFT_STRING               10
#define FIFFT_ASCII                10
#define FIFFT_LONG                 11
#define FIFFT_DAU_PACK13           13
#define FIFFT_DAU_PACK14           14
#define FIFFT_DAU_PACK16           16
#define FIFFT_COMPLEX_FLOAT        20
#define FIFFT_COMPLEX_DOUBLE       21
#define FIFFT_OLD_PACK             23
#define FIFFT_CH_INFO_STRUCT       30
#define FIFFT_ID_STRUCT            31
#define FIFFT_DIR_ENTRY_STRUCT     32
#define FIFFT_DIG_POINT_STRUCT     33
#define FIFFT_CH_POS_STRUCT        34
#define FIFFT_COORD_TRANS_STRUCT   35
#define FIFFT_DIG_STRING_STRUCT    36
#define FIFFT_STREAM_SEGMENT_STRUCT 37

/*
 * These are for matrices of any of the above 
 */
#define FIFFC_MATRIX_MAX_DIM  9

#define FIFFTS_FS_MASK        0xFF000000
#define FIFFTS_BASE_MASK      0x00000FFF
#define FIFFTS_MC_MASK        0x00FF0000

#define FIFFTS_FS_SCALAR      0x00000000
#define FIFFTS_FS_MATRIX      0x40000000	

#define FIFFTS_MC_DENSE       0x00000000
#define FIFFTS_MC_CCS         0x00100000
#define FIFFTS_MC_RCS         0x00200000


#define FIFFT_MATRIX           (FIFFTS_FS_MATRIX | FIFFTS_MC_DENSE)
#define FIFFT_CCS_MATRIX       (FIFFTS_FS_MATRIX | FIFFTS_MC_CCS)
#define FIFFT_RCS_MATRIX       (FIFFTS_FS_MATRIX | FIFFTS_MC_RCS)

#define FIFFT_MATRIX_INT       (FIFFT_MATRIX | FIFFT_INT)
#define FIFFT_MATRIX_FLOAT     (FIFFT_MATRIX | FIFFT_FLOAT)
#define FIFFT_MATRIX_DOUBLE    (FIFFT_MATRIX | FIFFT_DOUBLE)
#define FIFFT_CCS_MATRIX_FLOAT (FIFFT_CCS_MATRIX | FIFFT_FLOAT)
#define FIFFT_RCS_MATRIX_FLOAT (FIFFT_RCS_MATRIX | FIFFT_FLOAT)


#define FIFFM_MESSAGE_GROUP(x) ((x) / 100 == 0)

#define FIFF_NEW_FILE            1
#define FIFF_CLOSE_FILE          2
#define FIFF_DISCARD_FILE        3
#define FIFF_ERROR_MESSAGE       4
#define FIFF_SUSPEND_READING     5
#define FIFF_FATAL_ERROR_MESSAGE 6
#define FIFF_CONNECTION_CHECK    7
#define FIFF_SUSPEND_FILING      8
#define FIFF_RESUME_FILING       9
#define FIFF_RAW_PREBASE        10
#define FIFF_RAW_PICK_LIST      11
#define FIFF_ECHO               12
#define FIFF_RESUME_READING     13
#define FIFF_DACQ_SYSTEM_TYPE   14
#define FIFF_SELECT_RAW_CH      15  /* Instruct rawdisp to select this channel */
#define FIFF_PLAYBACK_MODE      16  /* Tell that we are playing data back from the hard
				     * disks in the data acquisition front end */
#define FIFF_CONTINUE_FILE      17  /* Used to inform that data is saved into a continuation file. */
#define FIFF_JITTER_MAX         18  /* Used to tell the jitter in the timing of data packets */
#define  FIFF_STREAM_SEGMENT    19  /* A segment of data stream */ 
/*
define FIFF_DECIMATION_FACTOR  19  * Collector; not used anywhere? 
*/

#define FIFFV_DACQ_SYSTEM_DAU     0
#define FIFFV_DACQ_SYSTEM_VXI     1
#define FIFFV_DACQ_SYSTEM_RPU     2
#define FIFFV_DACQ_SYSTEM_ORION   3 

#ifdef _DATA_SERVER
#define FIFF_MEM_DATA_BUFFER 10300  /* This is only used by
				     * cm_sender to indicate a data buffer 
				     * which is not logged into the file */
#endif

/*
 * Standard tags used in all blocks
 */

#define FIFF_FILE_ID         100
#define FIFF_DIR_POINTER     101
#define FIFF_DIR             102
#define FIFF_BLOCK_ID        103
#define FIFF_BLOCK_START     104
#define FIFF_BLOCK_END       105
#define FIFF_FREE_LIST       106
#define FIFF_FREE_BLOCK      107            
#define FIFF_NOP             108
#define FIFF_PARENT_FILE_ID  109
#define FIFF_PARENT_BLOCK_ID 110
#define FIFF_BLOCK_NAME      111
#define FIFF_BLOCK_VERSION   112
#define FIFF_CREATOR         113   /* Program that created the file (string)  */
#define FIFF_MODIFIER        114   /* Program that modified the file (string) */

#define FIFF_REF_ROLE        115
#define FIFF_REF_FILE_ID     116
#define FIFF_REF_FILE_NUM    117
#define FIFF_REF_FILE_NAME   118
/* reserverd                 119 */            
#define FIFF_REF_BLOCK_ID    120


/*
 * Megacq saves the parameters in these tags
 */

#define FIFF_DACQ_PARS        150
#define FIFF_DACQ_STIM        151

/* 
 * Structured objects (blocks)
 */

/*
 * MEG/EEG
 */

#define FIFFB_ROOT            999
#define FIFFB_MEAS            100
#define FIFFB_MEAS_INFO       101
#define FIFFB_RAW_DATA        102
#define FIFFB_PROCESSED_DATA  103
#define FIFFB_EVOKED          104
#define FIFFB_MCG_AVE         FIFFB_EVOKED
#define FIFFB_ASPECT          105
#define FIFFB_SUBJECT         106
#define FIFFB_ISOTRAK         107
#define FIFFB_HPI_MEAS        108
#define FIFFB_HPI_RESULT      109
#define FIFFB_HPI_COIL        110
#define FIFFB_PROJECT         111
#define FIFFB_CONTINUOUS_DATA 112
#define FIFFB_VOID            114
#define FIFFB_EVENTS          115
#define FIFFB_INDEX           116
#define FIFFB_DACQ_PARS       117
#define FIFFB_REF             118         /* Fiff referencing mechanism */

#define FIFFB_SMSH_RAW_DATA   119         /* SmartShield raw data */
#define FIFFB_SMSH_ASPECT     120         /* SmartShield averaged data */
#define FIFFB_HPI_SUBSYSTEM   121         /* HPI subsystem */


/*
 * MRI
 */

#define FIFFB_MRI             200         /* MRI/CT data. */
#define FIFFB_MRI_SET         201         /* MRI/CT volume */
#define FIFFB_MRI_SLICE       202         /* MRI/CT slice (image) */
#define FIFFB_MRI_SCENERY     203         /* These are for writing unrelated 'slices' */
#define FIFFB_MRI_SCENE       204	  /* Which are actually 3D scenes... */
#define FIFFB_MRI_SEG         205         /* MRI segmentation data */
#define FIFFB_MRI_SEG_REGION  206         /* One MRI segmentation region */

/*
 * Forward and inverse modelling
 */

#define FIFFB_SPHERE          300	  /* Concentric sphere model related */
#define FIFFB_BEM             310	  /* Boundary-element method */
#define FIFFB_BEM_SURF        311	  /* Boundary-element method surfaces */
#define FIFFB_CONDUCTOR_MODEL 312	  /* One conductor model definition */
#define FIFFB_XFIT_PROJ       313         /* xfit saves the linear projection 
					   * information here */
#define FIFFB_XFIT_PROJ_ITEM  314         /* Each projection item goes here */
#define FIFFB_XFIT_AUX        315         /* xfit saves its auxliary data 
					   * into this block */
/*                            350...
 *                            370          Reserved for MNE estimates (MHa) 
 */

#define FIFFB_BAD_CHANNELS    359         /* Alias of FIFFB_MNE_BAD_CHANNELS */

/*
 * Volume info
 */
#define FIFFB_VOL_INFO        400

/*
 * Sparse matrix, cross-talk correction, and SSS blocks
 */
#define FIFFB_DATA_CORRECTION     500     /* Correction to data */
#define FIFFB_CHANNEL_DECOUPLER   501     /* Cross-talk correction  */
#define FIFFB_SSS_INFO            502     /* SSS processing info */
#define FIFFB_SSS_CAL_ADJUST      503     /* Fine-calibration adjustment data */
#define FIFFB_SSS_ST_INFO         504     /* TSSS info */
#define FIFFB_SSS_BASES           505     /* SSS bases */

#define FIFFB_SMARTSHIELD         510     /* SmartShield data */


/*
 * Different aspects of data
 */

#define FIFFV_ASPECT_AVERAGE       100	  /* Normal average of epochs */
#define FIFFV_ASPECT_STD_ERR       101	  /* Std. error of mean */
#define FIFFV_ASPECT_SINGLE        102	  /* Single epoch cut out from the continuous data */
#define FIFFV_ASPECT_SUBAVERAGE    103	  
#define FIFFV_ASPECT_ALTAVERAGE    104	  /* Alternating subaverage */
#define FIFFV_ASPECT_SAMPLE        105	  /* A sample cut out by graph */
#define FIFFV_ASPECT_POWER_DENSITY 106    /* Power density spectrum */
#define FIFFV_ASPECT_DIPOLE_WAVE   200    /* Dipole amplitude curve */

/*
 * Tags used in data files
 */

#define FIFF_NCHAN           200	  /* Number of channels */
#define FIFF_SFREQ           201	  /* Sampling frequency (Hz) */
#define FIFF_DATA_PACK       202	  /* How the raw data is packed */
#define FIFF_CH_INFO         203          /* Channel descriptor */
#define FIFF_MEAS_DATE       204          /* Measurement date */
#define FIFF_SUBJECT         205	  /* This might be deleted later */
#define FIFF_COMMENT         206          /* This is used in a questionable way... */
#define FIFF_NAVE            207          /* Number of averages */
#define FIFF_FIRST_SAMPLE    208          /* The first sample of an epoch */
#define FIFF_LAST_SAMPLE     209          /* The last sample of an epoch */
#define FIFF_ASPECT_KIND     210          /* Aspect label */
#define FIFF_REF_EVENT       211          /* Reference event */
#define FIFF_EXPERIMENTER    212          /* Experimenter name */
#define FIFF_DIG_POINT       213          /* Digitization point */
#define FIFF_CH_POS_VEC      214          /* Channel positions */
#define FIFF_HPI_SLOPES      215          /* HPI data */
#define FIFF_HPI_NCOIL       216          /* Number of HPI coils */
#define FIFF_REQ_EVENT       217          /* Required event */
#define FIFF_REQ_LIMIT       218          /* Window for required event */
#define FIFF_LOWPASS         219	  /* Analog lowpass */
#define FIFF_BAD_CHS         220          /* List of bad channels */
#define FIFF_ARTEF_REMOVAL   221	  /* Artifact removal */
#define FIFF_COORD_TRANS     222	  /* Coordinate transformation */
#define FIFF_HIGHPASS        223	  /* Analog highpass */
#define FIFF_CH_CALS_VEC     224	  /* This will not occur in new files */
#define FIFF_HPI_BAD_CHS     225          /* List of channels considered to be bad in hpi */
#define FIFF_HPI_CORR_COEFF  226	  /* Hpi curve fit correlations */
#define FIFF_EVENT_COMMENT   227          /* Comment about the events used in averaging */
#define FIFF_NO_SAMPLES      228          /* Number of samples in an epoch */
#define FIFF_FIRST_TIME      229          /* Time scale minimum */
#define FIFF_SUBAVE_SIZE     230	  /* Size of a subaverage */
#define FIFF_SUBAVE_FIRST    231	  /* The first epoch # contained in the
					   * subaverage */
#define FIFF_NAME            233          /* Intended to be a short name. */
#define FIFF_DESCRIPTION     FIFF_COMMENT /* (Textual) Description of an object */
#define FIFF_DIG_STRING      234          /* String of digitized points */

#define FIFF_LINE_FREQ       235          /* Line interference frequency */
#define FIFF_HPI_COIL_FREQ   236          /* HPI coil excitation frequency */
#define FIFF_SIGNAL_CHANNEL  237          /* Signal channel name */

#define FIFFC_HPI_MAX_NCOIL 1000          /* Max value for FIFF_HPI_NCOIL */

/*
 *
 * HPI fitting program tags
 *
 */
#define FIFF_HPI_COIL_MOMENTS       240	  /* Estimated moment vectors for the HPI coil magnetic dipoles */
#define FIFF_HPI_FIT_GOODNESS       241	  /* Three floats indicating the goodness of fit */
#define FIFF_HPI_FIT_ACCEPT         242	  /* Bitmask indicating acceptance (see below) */
#define FIFF_HPI_FIT_GOOD_LIMIT     243	  /* Limit for the goodness-of-fit */
#define FIFF_HPI_FIT_DIST_LIMIT     244	  /* Limit for the coil distance difference */
#define FIFF_HPI_COIL_NO            245	  /* Coil number listed by HPI measurement */
#define FIFF_HPI_COILS_USED         246	  /* List of coils finally used when the transformation
					   * was computed */
#define FIFF_HPI_DIGITIZATION_ORDER 247	  /* Which Isotrak digitization point corresponds to 
					   * each of the coils energized */

#define FIFFV_HPI_ACCEPT_PROGRAM (1<<0)   /* The fit was accepted by the software */
#define FIFFV_HPI_ACCEPT_USER    (1<<1)	  /* The fit was accepted by user action */
#define FIFFV_HPI_ACCEPT_NONE    0

/* Following corresponsds to fields in channel info record. */

#define FIFF_CH_SCAN_NO	           250    /* int32 "Channel scan number. Corresponds to fiffChInfoRec.scanNo field" */
#define FIFF_CH_LOGICAL_NO         251	  /* int32 "Channel logical number. Corresponds to fiffChInfoRec.logNo field" */
#define FIFF_CH_KIND 	           252    /* enum(ch_kind)  "Channel type. Corresponds to fiffChInfoRec.kind field" */
#define FIFF_CH_RANGE              253    /* float          "Conversion from recorded number to (possibly virtual) voltage at the output" */
#define FIFF_CH_CAL                254    /* float          "Calibration coefficients from output voltage to some real units" */
#define FIFF_CH_POS                255    /* ch_pos_rec     "Channel position" */
#define FIFF_CH_UNIT               256    /* enum(unit)     "Unit of the data" */
#define FIFF_CH_UNIT_MUL           257	  /* int            "Unit multiplier exponent. */
#define FIFF_CH_DACQ_NAME          258	  /* string "Name of the channel in the data acquisition system. Same as fiffChInfoRec.name." */

#define FIFF_SSS_FRAME             263    /* SSS coordinate frame */
#define FIFF_SSS_JOB               264    /* SSS job */
#define FIFF_SSS_ORIGIN            265    /* Origin of the SSS inside expansion */
#define FIFF_SSS_ORD_IN            266    /* Order of the SSS inside expansion */
#define FIFF_SSS_ORD_OUT           267    /* Order of the SSS outside expansion */
#define FIFF_SSS_NMAG              268    /* Number of MEG channels */
#define FIFF_SSS_COMPONENTS        269    /* Number of SSS moments */
#define FIFF_SSS_CAL_CHANS         270    /* INT matrix (nmag x 2) for fine-calibrated channel numbers and types */
#define FIFF_SSS_CAL_CORRS         271    /* FLOAT matrix (nmag x 14) for fine-calibration coefficients */
#define FIFF_SSS_ST_CORR           272    /* TSSS subspace correlation */

#define FIFF_GANTRY_TYPE           280    /* enum(gantry_type): Type of the gantry. */
#define FIFF_GANTRY_MODEL          281    /* string: Gantry model. Prefereably the part number. */
#define FIFF_GANTRY_ANGLE          282    /* Tilt angle of the gantry in degrees. */

/* enum(gantry_type) */

#define FIFFV_GANTRY_TYPE_FIXED     1     /* Fixed gantry. */
#define FIFFV_GANTRY_TYPE_UNIAXIAL  2     /* Uni-axial, dewar rotation axis parallel to device y-axis */
#define FIFFV_GANTRY_TYPE_5DEGREES  3     /* Freely adjustable 5 degrees of freedom dewar. */


#define FIFFV_SSS_JOB_NOTHING   0         /* Just copy input to output */
#define FIFFV_SSS_JOB_CTC       1         /* Only cross-talk correction */
#define FIFFV_SSS_JOB_FILTER    2         /* Spatial filtering (default) */
#define FIFFV_SSS_JOB_VIRT      3         /* Reconstruct virtual data */
#define FIFFV_SSS_JOB_HEAD_POS  4         /* Estimate head positions, no movecomp */
#define FIFFV_SSS_JOB_MOVEC_FIT 5         /* Estimate head positions, do movecomp */
#define FIFFV_SSS_JOB_MOVEC_QUA 6         /* Do movecomp from prev estimated positions */
#define FIFFV_SSS_JOB_REC_ALL   7         /* Reconstruct fields from inside and outside moments */
#define FIFFV_SSS_JOB_REC_IN    8         /* Reconstruct fields from inside moments */
#define FIFFV_SSS_JOB_REC_OUT   9         /* Reconstruct fields from outside moments */
#define FIFFV_SSS_JOB_ST       10         /* TSSS */

#define FIFF_SSS_BASE_IN       273        /* DOUBLE matrix for SSS inside basis */
#define FIFF_SSS_BASE_OUT      274        /* DOUBLE matrix for SSS outside basis */
#define FIFF_SSS_BASE_VIRT     275        /* DOUBLE matrix for SSS virtual basis */
#define FIFF_SSS_NORM          276        /* Froebius norm of the inside basis */
#define FIFF_SSS_ITERATE       277        /* Nr of iterations in iterative SSS pseudo-inverse */
#define FIFF_DATA_BUFFER       300        /* Buffer containing measurement data */
#define FIFF_DATA_SKIP         301        /* Data skip in buffers */
#define FIFF_EPOCH             302        /* Buffer containing one epoch and channel */
#define FIFF_DATA_SKIP_SAMP    303        /* Data skip in samples */
#define FIFF_DATA_BUFFER2      304        /* Int_32 data buffer in dacq outgen */
#define FIFF_TIME_STAMP        305        /* Int_32[3] meas time stamp (sec,usec,sampleno) */

#define FIFF_SUBJ_ID           400        /* Subject ID */
#define FIFF_SUBJ_FIRST_NAME   401        /* First name of the subject */
#define FIFF_SUBJ_MIDDLE_NAME  402        /* Middle name of the subject */
#define FIFF_SUBJ_LAST_NAME    403        /* Last name of the subject */
#define FIFF_SUBJ_BIRTH_DAY    404        /* Birthday of the subject */
#define FIFF_SUBJ_SEX          405        /* Sex of the subject */
#define FIFF_SUBJ_HAND         406        /* Handedness of the subject */
#define FIFF_SUBJ_WEIGHT       407        /* Weight of the subject */
#define FIFF_SUBJ_HEIGHT       408        /* Height of the subject */
#define FIFF_SUBJ_COMMENT      409        /* Comment about the subject */
#define FIFF_SUBJ_HIS_ID       410        /* ID used in the Hospital Information System */

#define FIFF_PROJ_ID           500        /* Project ID */
#define FIFF_PROJ_NAME         501        /* Project name */
#define FIFF_PROJ_AIM          502        /* Projct description */
#define FIFF_PROJ_PERSONS      503        /* Persons participating in the project */
#define FIFF_PROJ_COMMENT      504        /* Comment about the project */

/* Special values used in the project tags. */

#define FIFFV_SEX_MALE   1
#define FIFFV_SEX_FEMALE 2

#define FIFFV_HAND_RIGHT 1
#define FIFFV_HAND_LEFT  2

/*
 * Event list saving...
 */
#define FIFF_EVENT_CHANNELS    600	/* Event channel numbers */
#define FIFF_EVENT_LIST        601      /* List of events (integers: 
					 * <sample before after> */
/*
 * Event spec tags
 */
#define FIFF_EVENT_CHANNEL     602	/* Event channel name */
#define FIFF_EVENT_BITS        603      /* Event bits array */

/*
 * Event bitmask constants
 */
#define FIFFC_EVENT_FROMMASK   0
#define FIFFC_EVENT_FROMBITS   1
#define FIFFC_EVENT_TOMASK     2
#define FIFFC_EVENT_TOBITS     3

/* 
 * Tags used in saving SQUID characteristics etc.
 */
#define FIFF_SQUID_BIAS        701
#define FIFF_SQUID_OFFSET      702
#define FIFF_SQUID_GATE        703

/* 
 * Tags for sparse matrices
 */
#define FIFF_DECOUPLER_MATRIX     800                        
#define FIFF_SPARSE_CH_NAME_LIST  FIFF_PROJ_ITEM_CH_NAME_LIST

/*
 * Processing history tags
 */
#define FIFFB_PROCESSING_HISTORY 900     /* Processing history block */
#define FIFFB_PROCESSING_RECORD  901     /*  .. can contain several processing records */

/* 
 * Aspect values used to save characteristic curves of SQUIDs.
 */
#define FIFFV_ASPECT_IFII_LOW  1100
#define FIFFV_ASPECT_IFII_HIGH 1101
#define FIFFV_ASPECT_GATE      1102
/*
 * References
 */
#define FIFF_REF_PATH           1101
/*
 * MRI...
 */
#define FIFF_MRI_SOURCE_PATH               FIFF_REF_PATH
#define FIFF_MRI_SOURCE_FORMAT             2002
#define FIFF_MRI_PIXEL_ENCODING            2003
#define FIFF_MRI_PIXEL_DATA_OFFSET         2004
#define FIFF_MRI_PIXEL_SCALE               2005
#define FIFF_MRI_PIXEL_DATA                2006
#define FIFF_MRI_PIXEL_OVERLAY_ENCODING    2007
#define FIFF_MRI_PIXEL_OVERLAY_DATA        2008

#define FIFF_MRI_BOUNDING_BOX              2009
#define FIFF_MRI_WIDTH                     2010
#define FIFF_MRI_WIDTH_M                   2011
#define FIFF_MRI_HEIGHT                    2012
#define FIFF_MRI_HEIGHT_M                  2013
#define FIFF_MRI_DEPTH                     2014
#define FIFF_MRI_DEPTH_M                   2015
#define FIFF_MRI_THICKNESS                 2016
#define FIFF_MRI_SCENE_AIM                 2017

#define FIFF_MRI_ORIG_SOURCE_PATH          2020
#define FIFF_MRI_ORIG_SOURCE_FORMAT        2021
#define FIFF_MRI_ORIG_PIXEL_ENCODING       2022
#define FIFF_MRI_ORIG_PIXEL_DATA_OFFSET    2023

#define FIFF_MRI_VOXEL_DATA                2030
#define FIFF_MRI_VOXEL_ENCODING            2031

#define FIFF_MRI_MRILAB_SETUP              2100

#define FIFF_MRI_SEG_REGION_ID             2200

/*
 * FIFF_MRI_SOURCE_FORMAT can be one of the following
 *                        A missing FIFF_MRI_SOURCE_FORMAT tag
 *                        indicates that the data is actually in the 
 *                        fiff itself (= FIFF_MRI_FORMAT_FIFF)
 * 
 * If the source format is FIFF_MRI_FORMAT_FIFF 
 * the tags FIFF_MRI_PIXEL_ENCODING and FIFF_MRI_PIXEL_DATA_OFFSET 
 * are missing and should be found by the software loading the data
 * from the FIFF_MRI_PIXEL_DATA tag.
 * 
 */
#define FIFFV_MRI_FORMAT_UNKNOWN           0
#define FIFFV_MRI_FORMAT_MAGNETOM_SHORT    1
#define FIFFV_MRI_FORMAT_MAGNETOM_LONG     2
#define FIFFV_MRI_FORMAT_MERIT             3
#define FIFFV_MRI_FORMAT_SIGNA             4
#define FIFFV_MRI_FORMAT_PIXEL_RAW         5
#define FIFFV_MRI_FORMAT_PIXEL_RAW_PACKED  6
#define FIFFV_MRI_FORMAT_MAGNETOM_NEW      7
#define FIFFV_MRI_FORMAT_FIFF              8
#define FIFFV_MRI_FORMAT_ACR_NEMA          9
#define FIFFV_MRI_FORMAT_DICOM_3          10
#define FIFFV_MRI_FORMAT_VISTA            11
/*
 * FIFF_MRI_PIXEL_ENCODING is one of the following
 */
#define FIFFV_MRI_PIXEL_UNKNOWN             0
#define FIFFV_MRI_PIXEL_BYTE                1
#define FIFFV_MRI_PIXEL_WORD                2
#define FIFFV_MRI_PIXEL_SWAP_WORD           3
#define FIFFV_MRI_PIXEL_FLOAT               4
#define FIFFV_MRI_PIXEL_BYTE_INDEXED_COLOR  5
#define FIFFV_MRI_PIXEL_BYTE_RGB_COLOR      6
#define FIFFV_MRI_PIXEL_BYTE_RLE_RGB_COLOR  7
#define FIFFV_MRI_PIXEL_BIT_RLE             8
/*
 * Forward and inverse modelling...
 */
/*
 * Sphere model     (3000...)
 */
#define FIFF_CONDUCTOR_MODEL_KIND 3000     /* What kind of conductor model */
/*
 * These are the models we support
 */
#define FIFFV_COND_MODEL_UNKNOWN     0      /* Not known */
#define FIFFV_COND_MODEL_SPHERE      1      /* Spherically symmetric */
#define FIFFV_COND_MODEL_BEM_HOMOG   2      /* Homogeneous BEM model */
#define FIFFV_COND_MODEL_BEM         3      /* Multilayer BEM model */

#define FIFF_SPHERE_ORIGIN          3001
#define FIFF_SPHERE_COORD_FRAME     3002   /* Which coordinate frame are we using? */
#define FIFF_SPHERE_LAYERS          3003   /* Array of layer structures */
/*
 * Surfaces for BEM (3100...)
 */
#define FIFF_BEM_SURF_ID            3101   /* int    surface number */
#define FIFF_BEM_SURF_NAME          3102   /* string surface name */
#define FIFF_BEM_SURF_NNODE	    3103   /* int    # of nodes on a surface */
#define FIFF_BEM_SURF_NTRI	    3104   /* int    # number of triangles on a surface */
#define FIFF_BEM_SURF_NODES         3105   /* float  surface nodes (nnode,3) */
#define FIFF_BEM_SURF_TRIANGLES     3106   /* int    surface triangles (ntri,3) */
#define FIFF_BEM_SURF_NORMALS       3107   /* float  surface node normal unit vectors (nnode,3) */
#define FIFF_BEM_SURF_CURVS         3108   /* float  surface node first principal curvature unit 
					    * vectors (nnode,3) */
#define FIFF_BEM_SURF_CURV_VALUES   3109   /* float  the two curvature values (nnode,2) */


#define FIFF_BEM_POT_SOLUTION       3110   /* float ** The solution matrix */
#define FIFF_BEM_APPROX             3111   /* int    approximation method, see below */
#define FIFF_BEM_COORD_FRAME        3112   /* The coordinate frame of the model */
#define FIFF_BEM_SIGMA              3113   /* Conductivity of a compartment */
/*
 * FIFF_BEM_SURF_ID can be one of the following
 */
#define FIFFV_BEM_SURF_ID_UNKNOWN    -1
#define FIFFV_BEM_SURF_ID_BRAIN       1
#define FIFFV_BEM_SURF_ID_CSF         2
#define FIFFV_BEM_SURF_ID_SKULL       3
#define FIFFV_BEM_SURF_ID_HEAD        4

#define FIFFV_BEM_SURF_ID_BLOOD      11
#define FIFFV_BEM_SURF_ID_HEART      12
#define FIFFV_BEM_SURF_ID_LUNGS      13
#define FIFFV_BEM_SURF_ID_TORSO      14

#define FIFFV_BEM_SURF_ID_NM122      21
#define FIFFV_BEM_SURF_UNIT_SPHERE   22
#define FIFFV_BEM_SURF_ID_VV         23
/*
 * FIFF_MRI_SEG_REGION_ID can be one of the following
 */
#define FIFFV_SEG_REGION_ID_UNKNOWN  FIFF_BEM_SURF_ID_UNKNOWN         
#define FIFFV_SEG_REGION_ID_BRAIN    FIFF_BEM_SURF_ID_BRAIN     
#define FIFFV_SEG_REGION_ID_CSF      FIFF_BEM_SURF_ID_CSF       
#define FIFFV_SEG_REGION_ID_SKULL    FIFF_BEM_SURF_ID_SKULL     
#define FIFFV_SEG_REGION_ID_HEAD     FIFF_BEM_SURF_ID_HEAD      
				                               
#define FIFFV_SEG_REGION_ID_BLOOD    FIFF_BEM_SURF_ID_BLOOD     
#define FIFFV_SEG_REGION_ID_HEART    FIFF_BEM_SURF_ID_HEART     
#define FIFFV_SEG_REGION_ID_LUNGS    FIFF_BEM_SURF_ID_LUNGS     
#define FIFFV_SEG_REGION_ID_TORSO    FIFF_BEM_SURF_ID_TORSO     
/*				                               
 * FIFF_BEM_APPROX		     
 */				     
#define FIFFV_BEM_APPROX_CONST        1     /* The constant potential approach */
#define FIFFV_BEM_APPROX_LINEAR       2     /* The linear potential approach */
/*
 * Source descriptions (3200...)
 * The dipole is six floats (position and dipole moment)
 */
#define FIFF_SOURCE_DIPOLE        3201
/*
 * These tags are used by xfit
 */
#define FIFF_XFIT_LEAD_PRODUCTS               3401
#define FIFF_XFIT_MAP_PRODUCTS                3402
#define FIFF_XFIT_GRAD_MAP_PRODUCTS           3403
#define FIFF_XFIT_VOL_INTEGRATION             3404
#define FIFF_XFIT_INTEGRATION_RADIUS          3405
#define FIFF_XFIT_CONDUCTOR_MODEL_NAME        3406
#define FIFF_XFIT_CONDUCTOR_MODEL_TRANS_NAME  3407
#define FIFF_XFIT_CONT_SURF_TYPE              3408   /* Xfit contour surface type */

/*
 * These relate to linear projection
 */
#define FIFF_PROJ_ITEM_KIND          3411
#define FIFF_PROJ_ITEM_TIME          3412
#define FIFF_PROJ_ITEM_DIPOLE        FIFF_SOURCE_DIPOLE
#define FIFF_PROJ_ITEM_IGN_CHS       3413
#define FIFF_PROJ_ITEM_NVEC          3414
#define FIFF_PROJ_ITEM_VECTORS       3415
#define FIFF_PROJ_ITEM_COMMENT       FIFF_COMMENT
#define FIFF_PROJ_ITEM_DESCRIPTION   FIFF_DESCRIPTION
#define FIFF_PROJ_ITEM_DEFINITION    3416
#define FIFF_PROJ_ITEM_CH_NAME_LIST  3417

#define FIFF_XFIT_PROJ_ITEM_KIND     FIFF_PROJ_ITEM_KIND
#define FIFF_XFIT_PROJ_ITEM_TIME     FIFF_PROJ_ITEM_TIME
#define FIFF_XFIT_PROJ_ITEM_DIPOLE   FIFF_PROJ_ITEM_DIPOLE
#define FIFF_XFIT_PROJ_ITEM_IGN_CHS  FIFF_PROJ_ITEM_IGN_CHS
#define FIFF_XFIT_PROJ_ITEM_NVEC     FIFF_PROJ_ITEM_NVEC
#define FIFF_XFIT_PROJ_ITEM_VECTORS  FIFF_PROJ_ITEM_VECTORS
#define FIFF_XFIT_PROJ_ITEM_COMMENT  FIFF_PROJ_ITEM_COMMENT
/*
 * The FIFF_PROJ_ITEM_KIND is an integer,
 * one of the following
 */
#define FIFFV_PROJ_ITEM_NONE        0
#define FIFFV_PROJ_ITEM_FIELD       1
#define FIFFV_PROJ_ITEM_DIP_FIX     2
#define FIFFV_PROJ_ITEM_DIP_ROT     3
#define FIFFV_PROJ_ITEM_HOMOG_GRAD  4
#define FIFFV_PROJ_ITEM_HOMOG_FIELD 5
#define FIFFV_PROJ_ITEM_EEG_AVREF   10         

#define FIFFV_XFIT_PROJ_ITEM_NONE        FIFFV_PROJ_ITEM_NONE
#define FIFFV_XFIT_PROJ_ITEM_FIELD       FIFFV_PROJ_ITEM_FIELD
#define FIFFV_XFIT_PROJ_ITEM_DIP_FIX     FIFFV_PROJ_ITEM_DIP_FIX
#define FIFFV_XFIT_PROJ_ITEM_DIP_ROT     FIFFV_PROJ_ITEM_DIP_ROT
#define FIFFV_XFIT_PROJ_ITEM_HOMOG_GRAD  FIFFV_PROJ_ITEM_HOMOG_GRAD
#define FIFFV_XFIT_PROJ_ITEM_HOMOG_FIELD FIFFV_PROJ_ITEM_HOMOG_FIELD

#define FIFF_XPLOTTER_LAYOUT          3501     /* xplotter layout tag */

/*  FIFF_MNE_xxxx                     3502
 *                                     ...
 *  Reserved for MNE data             3799 
 */
#define FIFF_CH_NAME_LIST             3507     /* Alias of FIFF_MNE_CH_NAME_LIST */

/*
 * These occur in the volume info files
 */
#define FIFF_VOL_ID                  4001
#define FIFF_VOL_NAME                4002
#define FIFF_VOL_OWNER_ID            4003      /* User id of the owner */
#define FIFF_VOL_OWNER_NAME          4004      /* User name of the owner */
#define FIFF_VOL_OWNER_REAL_NAME     4005      /* User name of the owner */
#define FIFF_VOL_TYPE                4006      /* See below... */
#define FIFF_VOL_HOST                4007      /* Where does the volume reside */
#define FIFF_VOL_REAL_ROOT           4008      /* The root of the volume in
						* in the machine where the file 
						* system is mounted */
#define FIFF_VOL_SYMBOLIC_ROOT       4009      /* Symbolic link to the 
						* root of the volume (if any)
						* system is mounted */
#define FIFF_VOL_MOUNT_POINT         4010      /* Last mount point of the volume */
#define FIFF_VOL_BLOCKS              4011      /* Total # of blocks */
#define FIFF_VOL_FREE_BLOCKS         4012      /* # of free blocks */
#define FIFF_VOL_AVAIL_BLOCKS        4013      /* # of free blocks available to non-superuser */
#define FIFF_VOL_BLOCK_SIZE          4014      /* Block size in bytes */
#define FIFF_VOL_DIRECTORY           4015      /* Contents of the volume 
						  in a special format 
						  the data type will be
						  FIFF_VOID */
/*
 * Index
 */
#define FIFF_INDEX_KIND              5001
#define FIFF_INDEX                   5002


/*======================================================================
 * Enumerated types used as tag values.
 *
 *=====================================================================*/

/* Values for FIFF_REF_ROLE. The role of a reference */

#define FIFFV_ROLE_PREV_FILE   1
#define FIFFV_ROLE_NEXT_FILE   2

/*
 * Method by which a projection is defined (FIFF_PROJ_ITEM_DEFINITION).
 * If tag is not present, FIFF_PROJ_BY_COMPLEMENT should be assumed.
 */
#define FIFFV_PROJ_BY_COMPLEMENT     0
#define FIFFV_PROJ_BY_SPACE          1


 /* Volume types used in FIFF_VOL_TYPE */

#define FIFFV_VOL_TYPE_HD            1	       /* Hard disk */
#define FIFFV_VOL_TYPE_MOD           2	       /* Magneto-optical disk */


#endif
