/*
 * $Id$
 * John Ashburner
 */

/* Matlab dependent high level data access and map manipulation routines */

#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#ifdef SPM_WIN32
#include <windows.h>
#include <memory.h>
#else
#include <unistd.h>
#include <sys/mman.h>
#endif

#include "spm_mapping.h"
#include "spm_datatypes.h"

/**************************************************************************/

void free_maps(MAPTYPE *maps, int n)
{
	int j;
	for(j=0; j<n; j++)
	{
		if (maps[j].addr)
		{
#ifdef SPM_WIN32
	(void)UnmapViewOfFile((LPVOID)(maps[j].addr));
#else
	(void)munmap((caddr_t)maps[j].addr, maps[j].len);
#endif
			maps[j].addr=0;
		}
		if (maps[j].data)
		{
			(void)mxFree((char *)maps[j].data);
			maps[j].data=0;
		}
		if (maps[j].scale)
		{
			(void)mxFree((char *)maps[j].scale);
			maps[j].scale=0;
		}
		if (maps[j].offset)
		{
			(void)mxFree((char *)maps[j].offset);
			maps[j].offset=0;
		}
	}
	(void)mxFree((char *)maps);
}

/**************************************************************************/

static void get_map_dat(int i, const mxArray *ptr, MAPTYPE *maps)
{
	mxArray *tmp;
	double *pr;
	int num_dims, j, t, dtype = 0;
	const int *dims;
	unsigned char *dptr;

	tmp=mxGetField(ptr,i,"dat");
	if (tmp == (mxArray *)0)
	{
		free_maps(maps,i);
		mexErrMsgTxt("Cant find dat.");
	}
	if      (mxIsDouble(tmp)) dtype = SPM_DOUBLE;
	else if (mxIsSingle(tmp)) dtype = SPM_FLOAT;
	else if (mxIsInt32 (tmp)) dtype = SPM_SIGNED_INT;
	else if (mxIsUint32(tmp)) dtype = SPM_UNSIGNED_INT;
	else if (mxIsInt16 (tmp)) dtype = SPM_SIGNED_SHORT;
	else if (mxIsUint16(tmp)) dtype = SPM_UNSIGNED_SHORT;
	else if (mxIsInt8  (tmp)) dtype = SPM_SIGNED_CHAR;
	else if (mxIsUint8 (tmp)) dtype = SPM_UNSIGNED_CHAR;
	else
	{
		free_maps(maps,i);
		mexErrMsgTxt("Unknown volume datatype.");
	}
	dptr  = (unsigned char *)mxGetPr(tmp);

	num_dims = mxGetNumberOfDimensions(tmp);
	if (num_dims > 3)
	{
		free_maps(maps,i);
		mexErrMsgTxt("Too many dimensions.");
	}
	dims     = mxGetDimensions(tmp);
	for(j=0; j<num_dims; j++)
		maps[i].dim[j]=dims[j];
	for(j=num_dims; j<3; j++)
		maps[i].dim[j]=1;
	tmp=mxGetField(ptr,i,"dim");
	if (tmp != (mxArray *)0)
	{
		if (mxGetM(tmp)*mxGetN(tmp) != 3)
		{
			free_maps(maps,i);
			mexErrMsgTxt("Wrong sized dim.");
		}
		pr = mxGetPr(tmp);
		if (maps[i].dim[0] != (int)fabs(pr[0]) ||
		    maps[i].dim[1] != (int)fabs(pr[1]) ||
		    maps[i].dim[2] != (int)fabs(pr[2]))
		{
			free_maps(maps,i);
			mexErrMsgTxt("Incompatible volume dimensions in dim.");
		}
	}
	tmp=mxGetField(ptr,i,"dt");
	if (tmp != (mxArray *)0)
	{
		if (mxGetM(tmp)*mxGetN(tmp) != 1 && mxGetM(tmp)*mxGetN(tmp) != 2)
		{
			free_maps(maps,i);
			mexErrMsgTxt("Wrong sized dt.");
		}
		pr = mxGetPr(tmp);
		if (dtype != (int)fabs(pr[0]))
		{
			free_maps(maps,i);
			mexErrMsgTxt("Incompatible datatype in dt.");
		}
	}

	maps[i].addr      = 0;
	maps[i].len       = 0;
	maps[i].dtype  = dtype;
	maps[i].data   = (void  **)mxCalloc(maps[i].dim[2],sizeof(void *));
	maps[i].scale  = (double *)mxCalloc(maps[i].dim[2],sizeof(double));
	maps[i].offset = (double *)mxCalloc(maps[i].dim[2],sizeof(double));

	t   = maps[i].dim[0]*maps[i].dim[1]*get_datasize(maps[i].dtype)/8;
	tmp = mxGetField(ptr,i,"pinfo");
	if (tmp != (mxArray *)0)
	{
		if ((mxGetM(tmp) != 2 && mxGetM(tmp) != 3) || (mxGetN(tmp) != 1 && mxGetN(tmp) != maps[i].dim[2]))
		{
			free_maps(maps,i+1);
			mexErrMsgTxt("Wrong sized pinfo.");
		}
		if (mxGetM(tmp) == 3 && mxGetPr(tmp)[2] != 0)
		{
			free_maps(maps,i+1);
			mexErrMsgTxt("pinfo(3) must equal 0 to read dat field.");
		}
		pr = mxGetPr(tmp);
		if (mxGetN(tmp) == 1)
			for(j=0; j<maps[i].dim[2]; j++)
			{
				maps[i].scale[j]  = pr[0];
				maps[i].offset[j] = pr[1];
				maps[i].data[j]   = &(dptr[j*t]);
			}
		else
			for(j=0; j<maps[i].dim[2]; j++)
			{
				maps[i].scale[j]  = pr[0+j*2];
				maps[i].offset[j] = pr[1+j*2];
				maps[i].data[j]   = &(dptr[j*t]);
			}
	}
	else
		for(j=0; j<maps[i].dim[2]; j++)
		{
			maps[i].scale[j]  = 1.0;
			maps[i].offset[j] = 0.0;
			maps[i].data[j]   = &(dptr[j*t]);
		}
	tmp=mxGetField(ptr,i,"mat");
	if (tmp != (mxArray *)0)
	{
		if (mxGetM(tmp) != 4 || mxGetN(tmp) != 4)
		{
			free_maps(maps,i+1);
			mexErrMsgTxt("Wrong sized mat.");
		}
		pr = mxGetPr(tmp);
		for(j=0; j<16; j++)
			maps[i].mat[j] = pr[j];
	}
	else
	{
		for(j=0; j<16; j++)
			maps[i].mat[j] = 0.0;
		for(j=0; j<4; j++)
			maps[i].mat[j + j*4] = 1.0;
	}
}

/**************************************************************************/

static void get_map_file(int i, const mxArray *ptr, MAPTYPE *maps)
{
	int j;
	mxArray *tmp;
	double *pr;
	int dsize = 0, off;
#ifdef SPM_BIGENDIAN
	int be = 1;
#else
	int be = 0;
#endif

	tmp=mxGetField(ptr,i,"dim");
	if (tmp == (mxArray *)0)
	{
		free_maps(maps,i);
		mexErrMsgTxt("Cant find dim.");
	}
	if (mxGetM(tmp)*mxGetN(tmp) != 3)
	{
		free_maps(maps,i);
		mexErrMsgTxt("Wrong sized dim.");
	}
	pr = mxGetPr(tmp);
	maps[i].dim[0] = (int)fabs(pr[0]);
	maps[i].dim[1] = (int)fabs(pr[1]);
	maps[i].dim[2] = (int)fabs(pr[2]);
	maps[i].dtype  = (int)fabs(pr[3]);
	maps[i].data   = (void  **)mxCalloc(maps[i].dim[2],sizeof(void *));
	maps[i].scale  = (double *)mxCalloc(maps[i].dim[2],sizeof(double));
	maps[i].offset = (double *)mxCalloc(maps[i].dim[2],sizeof(double));

	tmp=mxGetField(ptr,i,"dt");
	if (tmp == (mxArray *)0)
	{
		free_maps(maps,i+1);
		mexErrMsgTxt("Cant find dt.");
	}
	if (mxGetM(tmp)*mxGetN(tmp) != 2)
	{
		free_maps(maps,i+1);
		mexErrMsgTxt("Wrong sized dt.");
	}
	pr = mxGetPr(tmp);
	if ((int)(pr[1])==be)
	{
		if      (pr[0]==    2) {maps[i].dtype = SPM_UNSIGNED_CHAR;  dsize =  8;}
		else if (pr[0]==  256) {maps[i].dtype = SPM_SIGNED_CHAR;    dsize =  8;}
		else if (pr[0]==    4) {maps[i].dtype = SPM_SIGNED_SHORT;   dsize = 16;}
		else if (pr[0]==  512) {maps[i].dtype = SPM_UNSIGNED_SHORT; dsize = 16;}
		else if (pr[0]==    8) {maps[i].dtype = SPM_SIGNED_INT;     dsize = 32;}
		else if (pr[0]==  768) {maps[i].dtype = SPM_UNSIGNED_INT;   dsize = 32;}
		else if (pr[0]==   16) {maps[i].dtype = SPM_FLOAT;          dsize = 32;}
		else if (pr[0]==   64) {maps[i].dtype = SPM_DOUBLE;         dsize = 64;}
		else {free_maps(maps,i+1);mexErrMsgTxt("Unknown datatype.");}
	}
	else
	{
		if      (pr[0]==    2) {maps[i].dtype = SPM_UNSIGNED_CHAR;    dsize =  8;}
		else if (pr[0]==  256) {maps[i].dtype = SPM_SIGNED_CHAR;      dsize =  8;}
		else if (pr[0]==    4) {maps[i].dtype = SPM_SIGNED_SHORT_S;   dsize = 16;}
		else if (pr[0]==  512) {maps[i].dtype = SPM_UNSIGNED_SHORT_S; dsize = 16;}
		else if (pr[0]==    8) {maps[i].dtype = SPM_SIGNED_INT_S;     dsize = 32;}
		else if (pr[0]==  768) {maps[i].dtype = SPM_UNSIGNED_INT_S;   dsize = 32;}
		else if (pr[0]==   16) {maps[i].dtype = SPM_FLOAT_S;          dsize = 32;}
		else if (pr[0]==   64) {maps[i].dtype = SPM_DOUBLE_S;         dsize = 64;}
		else {free_maps(maps,i+1);mexErrMsgTxt("Unknown datatype.");}
	}

	tmp=mxGetField(ptr,i,"fname");
	if (tmp == (mxArray *)0)
	{
		free_maps(maps,i+1);
		mexErrMsgTxt("Cant find fname.");
	}
	if (mxIsChar(tmp))
	{
#ifdef SPM_WIN32
		HANDLE hFile, hMapping;
#endif
		int buflen;
		char *buf;
		int fd;
		struct stat stbuf;
		buflen = mxGetN(tmp)*mxGetM(tmp)+1;
		buf    = mxCalloc(buflen,sizeof(char));
		if (mxGetString(tmp,buf,buflen))
		{
			mxFree(buf);
			free_maps(maps,i+1);
			mexErrMsgTxt("Cant get filename.");
		}
		if ((fd = open(buf, O_RDONLY)) == -1)
		{
			mxFree(buf);
			free_maps(maps,i+1);
			mexErrMsgTxt("Cant open image file.");
		}
		if (fstat(fd, &stbuf) == -1)
		{
			(void)close(fd);
			mxFree(buf);
			free_maps(maps,i+1);
			mexErrMsgTxt("Cant stat image file.");
		}
		maps[i].len = stbuf.st_size;
#ifdef SPM_WIN32
		(void)close(fd);

		/* maps[i].addr = map_file(buf, (caddr_t)0, maps[i].len,
		 *	PROT_READ, MAP_SHARED, (off_t)0); */

		/* http://msdn.microsoft.com/library/default.asp?
		       url=/library/en-us/fileio/base/createfile.asp */
		hFile = CreateFile(buf, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL | FILE_FLAG_RANDOM_ACCESS, NULL);
		mxFree(buf);
		if (hFile == NULL)
			mexErrMsgTxt("Cant open file.  It may be locked by another program.");

		/* http://msdn.microsoft.com/library/default.asp?
		       url=/library/en-us/fileio/base/createfilemapping.asp */
		hMapping = CreateFileMapping(hFile, NULL, PAGE_READONLY, 0, 0, NULL);
		(void)CloseHandle(hFile);
		if (hMapping == NULL)
			mexErrMsgTxt("Cant create file mapping.  It may be locked by another program.");

		/* http://msdn.microsoft.com/library/default.asp?
		       url=/library/en-us/fileio/base/mapviewoffile.asp */
		maps[i].addr    = (caddr_t)MapViewOfFileEx(hMapping, FILE_MAP_READ, 0, 0, maps[i].len, 0);
		(void)CloseHandle(hMapping);
		if (maps[i].addr == NULL)
			mexErrMsgTxt("Cant map view of file.  It may be locked by another program.");

#else
		maps[i].addr = mmap((caddr_t)0, maps[i].len,
			PROT_READ, MAP_SHARED, fd, (off_t)0);
		(void)close(fd);
		if (maps[i].addr == (void *)-1)
		{
			(void)perror("Memory Map");
			mxFree(buf);
			free_maps(maps,i+1);
			mexErrMsgTxt("Cant map image file.");
		}
		mxFree(buf);
#endif
	}


	tmp=mxGetField(ptr,i,"pinfo");
	if (tmp != (mxArray *)0)
	{
		if (mxGetM(tmp) != 3 || (mxGetN(tmp) != 1 && mxGetN(tmp) != maps[i].dim[2]))
		{
			free_maps(maps,i+1);
			mexErrMsgTxt("Wrong sized pinfo.");
		}
		pr = mxGetPr(tmp);
		if (mxGetN(tmp) == 1)
		{
			off = (int)fabs(pr[2]);
			if (off+maps[i].dim[0]*maps[i].dim[1]*maps[i].dim[2]*(dsize/8) > maps[i].len)
			{
				free_maps(maps,i+1);
				mexErrMsgTxt("File too small.");
			}
			for(j=0; j<maps[i].dim[2]; j++)
			{
				maps[i].scale[j]  = pr[0];
				maps[i].offset[j] = pr[1];
				maps[i].data[j]   = maps[i].addr+off+j*maps[i].dim[0]*maps[i].dim[1]*(dsize/8);
			}
		}
		else
		{
			for(j=0; j<maps[i].dim[2]; j++)
			{
				maps[i].scale[j]  = pr[0+j*3];
				maps[i].offset[j] = pr[1+j*3];
				off = (int)fabs(pr[2+j*3]);
				maps[i].data[j]   = maps[i].addr+off;
				if (off+maps[i].dim[0]*maps[i].dim[1]*(dsize/8) > maps[i].len)
				{
					free_maps(maps,i+1);
					mexErrMsgTxt("File too small.");
				}
			}
		}
	}
	else
	{
		if (maps[i].dim[0]*maps[i].dim[1]*maps[i].dim[2]*(dsize/8) > maps[i].len)
		{
			free_maps(maps,i+1);
			mexErrMsgTxt("File too small.");
		}
		for(j=0; j<maps[i].dim[2]; j++)
		{
			maps[i].scale[j]  = 1.0;
			maps[i].offset[j] = 0.0;
			maps[i].data[j]   = maps[i].addr+j*maps[i].dim[0]*maps[i].dim[1]*(dsize/8);
		}
	}

	tmp=mxGetField(ptr,i,"mat");
	if (tmp != (mxArray *)0)
	{
		if (mxGetM(tmp) != 4 || mxGetN(tmp) != 4)
		{
			free_maps(maps,i+1);
			mexErrMsgTxt("Wrong sized mat.");
		}
		pr = mxGetPr(tmp);
		for(j=0; j<16; j++)
			maps[i].mat[j] = pr[j];
	}
	else
	{
		for(j=0; j<16; j++)
			maps[i].mat[j] = 0.0;
		for(j=0; j<4; j++)
			maps[i].mat[j + j*4] = 1.0;
	}
}

/**************************************************************************/

static MAPTYPE *get_maps_struct(const mxArray *ptr, int *n)
{
	MAPTYPE *maps;
	int num_dims, i;
	const int *dims;

	if (!mxIsStruct(ptr))
	{
		mexErrMsgTxt("Not a structure.");
		return((MAPTYPE *)0);
	}
	num_dims  = mxGetNumberOfDimensions(ptr);
	dims      = mxGetDimensions(ptr);
	*n = 1;
	for(i=0; i<num_dims; i++)
		*n *= dims[i];

	maps = (MAPTYPE *)mxCalloc(*n, sizeof(MAPTYPE));

	if (*n > 0)
	{
		if (mxGetField(ptr,0,"dat") == (mxArray *)0)
			for (i=0; i< *n; i++)
				get_map_file(i, ptr, maps);
		else
			for (i=0; i< *n; i++)
				get_map_dat(i, ptr, maps);
	}
	return(maps);
}

/**************************************************************************/

static MAPTYPE *get_maps_3dvol(const mxArray *ptr, int *n)
{
	int num_dims, jj, t, dtype = 0;
	const int *dims;
	MAPTYPE *maps;
	unsigned char *dptr;

	if      (mxIsDouble(ptr)) dtype = SPM_DOUBLE;
	else if (mxIsSingle(ptr)) dtype = SPM_FLOAT;
	else if (mxIsInt32 (ptr)) dtype = SPM_SIGNED_INT;
	else if (mxIsUint32(ptr)) dtype = SPM_UNSIGNED_INT;
	else if (mxIsInt16 (ptr)) dtype = SPM_SIGNED_SHORT;
	else if (mxIsUint16(ptr)) dtype = SPM_UNSIGNED_SHORT;
	else if (mxIsInt8  (ptr)) dtype = SPM_SIGNED_CHAR;
	else if (mxIsUint8 (ptr)) dtype = SPM_UNSIGNED_CHAR;
	else mexErrMsgTxt("Unknown volume datatype.");

	maps = (MAPTYPE *)mxCalloc(1, sizeof(MAPTYPE));

	num_dims = mxGetNumberOfDimensions(ptr);
	if (num_dims > 3)
	{
		mexErrMsgTxt("Too many dimensions.");
	}

	dims     = mxGetDimensions(ptr);
	for(jj=0; jj<num_dims; jj++)
		maps->dim[jj]=dims[jj];
	for(jj=num_dims; jj<3; jj++)
		maps->dim[jj]=1;

	for(jj=0; jj<16; jj++)
		maps->mat[jj] = 0.0;
	for(jj=0; jj<4; jj++)
		maps->mat[jj + jj*4] = 1.0;

	maps->dtype  = dtype;

	maps->data   = (void  **)mxCalloc(maps->dim[2],sizeof(void *));
	maps->scale  = (double *)mxCalloc(maps->dim[2],sizeof(double));
	maps->offset = (double *)mxCalloc(maps->dim[2],sizeof(double));

	t     = maps->dim[0]*maps->dim[1]*get_datasize(maps->dtype)/8;
	dptr   = (unsigned char *)mxGetPr(ptr);

	for(jj=0; jj<maps->dim[2]; jj++)
	{
		maps->scale[jj]  = 1.0;
		maps->offset[jj] = 0.0;
		maps->data[jj]   = &(dptr[jj*t]);
	}

	maps->addr      = 0;
	maps->len       = 0;

	*n = 1;
	return(maps);
}

/**************************************************************************/

MAPTYPE *get_maps(const mxArray *ptr, int *n)
{
	if (mxIsStruct(ptr))
		return(get_maps_struct(ptr, n));
	else if (mxGetNumberOfDimensions(ptr) <= 3 &&
		mxIsNumeric(ptr) && !mxIsComplex(ptr) &&
		!mxIsSparse(ptr))
		return(get_maps_3dvol(ptr, n));
	else
		mexErrMsgTxt("What do I do with this?");
	return((MAPTYPE *)0);
}


/**************************************************************************/

void voxdim(MAPTYPE *map, double vdim[3])
{
	int i, j;
	double t;
	for(j=0; j<3; j++)
	{
		t=0.0;
		for(i=0; i<3; i++)
			t += map->mat[i+j*4]*map->mat[i+j*4];
		vdim[j] = sqrt(t);
	}
}

/**************************************************************************/

int get_dtype(const mxArray *ptr)
{
	mxArray *tmp;
	double *pr;

	if (!mxIsStruct(ptr))
	{
		mexErrMsgTxt("Not a structure.");
		return(0);
	}
	tmp=mxGetField(ptr,0,"dt");
	if (tmp == (mxArray *)0)
	{
		mexErrMsgTxt("Cant find dt.");
	}
	if (mxGetM(tmp)*mxGetN(tmp) != 2)
	{
		mexErrMsgTxt("Wrong sized dt.");
	}
	pr = mxGetPr(tmp);

	return((int)fabs(pr[0]));
}

/**************************************************************************/
