/*
 * $Id$
 * John Ashburner
 */
 
/* Matlab independent image access routines */

#ifndef _SPM_VOL_ACCESS_H_
#define _SPM_VOL_ACCESS_H_

#include <sys/types.h>
#ifdef SPM_WIN32
typedef char *caddr_t;
#else
#include <sys/mman.h>
#endif

typedef struct maptype
{
    int dim[3];		        /* Dimensions of the volume */
    double *scale, *offset;	/* Scalefactor and offset, such that true_intensity = vox*scale+offset */
    int dtype;		        /* Data-type of volume */
    void **data;	        /* Pointer to data */
    double mat[16];

    caddr_t addr;
    size_t len;
}   MAPTYPE;

int get_datasize(int type);

int resample(int m,MAPTYPE *vol,double *out,double *x, 
	double *y,double *z,int hold, double background);

int resample_d(int m,MAPTYPE *vol,double *out,
	double *gradx, double *grady, double *gradz,
	double *x, double *y, double *z,
	int hold, double background);

int slice(double *mat, double *image, int xdim1,int ydim1, 
	MAPTYPE *vol, int hold,double background);

#endif /* _SPM_VOL_ACCESS_H_ */
