#include <math.h>
#include <mex.h>

#define MAX_SIZE   1200

#define NAN_TO_ZERO(x)   (mxIsFinite(x) ? (x) : 0.0)

/** 1D smoothing filter. 'dest' must have at least 'nx' elements, but 'src' must have nx+2*off elements
    in the form (example off=3)
	    [s0 s0 s0   s0 s1 s2 s3 ... sn_1   sn_1 sn_1 sn_1]
	that is, 'src' contains the true source signal with extra padding (repeated value) at both ends.
	
	Used by smooth_slice. Could be further optimized by creating specialised routines for different 
	kernel sizes (loop unrolling).
*/
void smooth_one(float *dest, int nx, int off, const float *src, const float *kern) {
	int o2 = 2*off + 1;
	int i,j;
	
	for (i=0;i<nx;i++) {
		float h = 0.0;
		for (j=0;j<o2;j++) {
			h += src[i+j]*kern[j];
		}
		dest[i] = h;
	}
}

/* smooth a 2D float array of size (nx,ny) using the kernels 'kx' and 'ky'
   kx must have 2*offx+1 elements and be symmetric + same for ky!
*/
void smooth_slice(int nx, int ny, float *dest, const float *src, int xoff, const float *kx, int yoff, const float *ky) {
	int x, y;
	
	float auxi[MAX_SIZE];
	float auxo[MAX_SIZE];
	
	/* smooth along x direction first */
	for (y=0; y<ny; y++) {
		const float *s_y = src + y * nx;
		float *d_y = dest + y * nx;
		
		float firstY = NAN_TO_ZERO(s_y[0]);
		float lastY = NAN_TO_ZERO(s_y[nx-1]);
		
		for (x=0; x<xoff; x++) {
			auxi[x] = firstY;
			auxi[x+xoff+nx] = lastY;
		}
		for (x=0; x<nx; x++) {
			auxi[x+xoff] = NAN_TO_ZERO(s_y[x]);
		}
		
		smooth_one(auxo, nx, xoff, auxi, kx);
		
		for (x=0;x<nx;x++) {
			d_y[x] = auxo[x];
		}	
	}

	/* now smooth along y direction */
	for (x=0; x<nx; x++) {
		float *dx_ = dest + x;
		
		for (y=0; y<yoff; y++) {
			auxi[y] = dx_[0];
			auxi[y+yoff+ny] = dx_[(ny - 1)*nx];
		}
		for (y=0; y<ny; y++) {
			auxi[y+yoff] = dx_[y*nx];
		}
		
		smooth_one(auxo, ny, yoff, auxi, ky);		
		
		for (y=0;y<ny;y++) {
			dx_[y*nx] = auxo[y];
		}	
	}
	
}

/** Smooth a 3D float array of size (nx,ny,nz) using the kernels 'kx', 'ky' and 'kz'
	kx must have 2*offx+1 elements and be symmetric, same for ky and kz!
	'dest' must have the same size as 'src', source values at the boundaries are 
	repeated.
*/
void smooth_vol(int nx, int ny, int nz, float *dest, const float *src, int offx, const float *kx, int offy, const float *ky, int offz, const float *kz) {
	int z, xy, nxy;
	float auxi[MAX_SIZE];
	float auxo[MAX_SIZE];	
	
	nxy = nx*ny;
	
	for (z=0;z<nz;z++) {
		smooth_slice(nx,ny, dest + z*nxy, src + z*nxy, offx, kx, offy, ky);
	}
	
	/* smooth along z */
	
	for (xy=0; xy<nxy; xy++) {
		float *d = dest + xy;
		
		for (z=0; z<offz; z++) {
			auxi[z] = d[0];
			auxi[z+offz+nz] = d[(nz - 1)*nxy];
		}
		for (z=0; z<nz; z++) {
			auxi[z+offz] = d[z*nxy];
		}
		
		smooth_one(auxo, nz, offz, auxi, kz);
		
		for (z=0;z<nz;z++) {
			d[z*nxy] = auxo[z];
		}	
	}	
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int i, off[3];
	const float *src, *kern[3];
	float *dest;
	char err[400];
	const mwSize *dims;
	
	if (nrhs != 4) {
		mexErrMsgTxt("Usage: Vs = mri_smooth_vol(Vin, kernX, kernY, kernZ)");
	}
	
	for (i=0;i<4;i++) {
		if (!mxIsSingle(prhs[i]) || mxIsComplex(prhs[i])) {
			sprintf(err, "Argument nr. %i is not a real-valued single precision array!", i+1);
			mexErrMsgTxt(err);
		}
	}
	if (mxGetNumberOfDimensions(prhs[0])!=3) mexErrMsgTxt("First argument must be a 3D array");
	dims = mxGetDimensions(prhs[0]);
	src = (const float *) mxGetData(prhs[0]);

	for (i=0;i<3;i++) {
		off[i] = mxGetNumberOfElements(prhs[i+1]);
		
		if ((off[i] & 1) != 1) {
			sprintf(err, "%i-th smoothing kern is of even length - not allowed!", i+1);
			mexErrMsgTxt(err);
		}
		
		off[i] = off[i] >> 1;
		
		if (2*off[i] + dims[i] > MAX_SIZE) {
			sprintf(err, "Too large inputs (source + smoothing) along dimension %i\n", i+1);
			mexErrMsgTxt(err);
		}
		kern[i] = (const float *) mxGetData(prhs[i+1]);
	}
		
	plhs[0] = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);
	dest = (float *) mxGetData(plhs[0]);
	
	smooth_vol(dims[0], dims[1], dims[2], dest, src, off[0], kern[0], off[1], kern[1], off[2], kern[2]);
}
