/* 
 * $Id$
 * John Ashburner
 */
 
/* Matlab dependent high level data access and map manipulation routines */

#include "mex.h"
#include "spm_vol_access.h"

#ifndef _SPM_MAPPING_H_
#define _SPM_MAPPING_H_

void free_maps(MAPTYPE *maps, int n);

MAPTYPE *get_maps(const mxArray *ptr, int *n);

void voxdim(MAPTYPE *map, double vdim[3]);

int get_dtype(const mxArray *ptr);

#endif /* _SPM_MAPPING_H_ */
