/*
  © Copyright 2003, C&C Research Laboratories, NEC Europe Ltd.
  

    This file is part of SimBio-Vgrid.

    SimBio-Vgrid is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SimBio-Vgrid is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SimBio-Vgrid.  If not, see <http://www.gnu.org/licenses/>.


*/

/* \file
   \brief  Generate synthetic Vista images for testing purposes

    \author Guntram Berti

     $Id: genimg.C,v 1.12 2004/09/09 15:00:29 berti Exp $
 */

#include <stdio.h> 
#include <vista.h>

#include <cmath>
#include <algorithm>


class plane_eval {
  typedef float c[3];
  c points[4];
  c normals[4];

public:
  plane_eval() { init();}
  void init() {
    static c regtet[4] = {{1,0,-sqrt(2.0)}, {-1,0,-sqrt(2.0)}, {0,1,sqrt(2.0)}, {0,-1,sqrt(2.0)}};
    // regular tet
    for(int i = 0; i < 4; ++i)
      for(int j = 0; j < 3; ++j) {
	points[i][j] = regtet[i][j];
      }
    //    points[0] = {1,0,-sqrt(2.0)};
    //{-1,0,-sqrt(2.0)}, {0,1,sqrt(2.0)}, {0,-1,sqrt(2.0)}};
    typedef int index[3];
    static index idx[4] = { {0, 1, 2}, { 1, 0, 3}, {0, 2, 3}, {1, 3, 2} };
    for(int i = 0; i <= 3; ++i)
      init_plane(i, points[idx[i][0]],points[idx[i][1]],points[idx[i][2]]);
  }
  void init_plane(int i, c p1, c p2, c p3) 
  {
    c d1 = { p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]};
    c d2 = { p3[0]-p1[0], p3[1]-p1[1], p3[2]-p1[2]};
    c ni = { d1[1]*d2[2]-d1[2]*d2[1], d1[2]*d2[0]-d1[0]*d2[2], d1[0]*d2[1]-d1[1]*d2[0]};
    for(int j = 0; j < 3; ++j)
      normals[i][j] = ni[j]; 
    float norm_n = sqrt(dot(normals[i], normals[i]));
    normals[i][0] /= norm_n; 
    normals[i][1] /= norm_n; 
    normals[i][2] /= norm_n; 
  }
  float dot(const c p, const c q) const { return q[0]*p[0]+q[1]*p[1]+q[2]*p[2];}
  int operator()(c x) const {
    float f[4] = { dot(normals[0],x), dot(normals[1],x), dot(normals[2],x), dot(normals[3],x)};
    return std::max_element(f, f+4) -f;
  }
};

inline int mymin(int a, int b, int c) 
{ 
  int m = (a < b ? a : b); 
  return (m < c ? m : c);
}
inline int mymax(int a, int b, int c) 
{ 
  int m = (a > b ? a : b); 
  return (m > c ? m : c);
}

enum geomType { uniform, corner, ball, contball, contellipsoid, shell, two_halves, oblique, oblique2d, brick, checker8, tetra };
static VDictEntry geom_dict[] = {
  { "uniform", uniform}, 
  { "corner", corner},
  { "ball",  ball },
  { "contball",  contball },
  { "contellipsoid",  contellipsoid },
  { "shell", shell },
  { "two-halves", two_halves},
  { "oblique", oblique},
  { "oblique2d", oblique2d},
  { "brick", brick},
  { "checker8", checker8},
  { "tetra", tetra},
  { NULL }
};

int main (int argc, char *argv[])
{
        
        VAttrList list;
        VAttrListPosn posn;
        FILE *outf;
        
        VImage  dst; 
        int     n = 4,  nx=-1, ny=-1, nz=-1;
	int     c = 1;
	int     border = 0;
        FILE    *fp, *datei;
	geomType geom = corner;        
	float ball_radius   = 0.5; // relative radius of ball
        int   two_halves_sep = -1; // separator value for two halfs
	float shell_outer_radius = 0.6;
	float shell_inner_radius = 0.4;
	float ell_x = 1.0;
	float ell_y = 1.0;
	float ell_z = 1.0;
	float scale = -1.0;
	float maxcoord[3] = { -1, -1, -1};
	float mincoord[3] = { 0,   0,  0};
	int   tx = 0, ty = 0, tz = 0;
	int   maxval = 2;

	VBoolean unit_coords = 0;
	VBoolean bg          = 0; // use uniform background 0
	float clip_radius    = -1.0;

        static VOptionDescRec options[] = {
	  { "nx",  VLongRepn, 1, &nx, VOptionalOpt, NULL, "x image dimension" }
	  ,{ "ny",  VLongRepn, 1, &ny, VOptionalOpt, NULL, "y image dimension" }
	  ,{ "nz",  VLongRepn, 1, &nz, VOptionalOpt, NULL, "z image dimension" }
          ,{  "n",   VLongRepn, 1, &n,  VOptionalOpt, NULL, 
	      "x,y,z image dimension. Sets nx,ny,nz if those are unset (= -1)" }
          ,{  "b",   VLongRepn, 1, &border,  VOptionalOpt, NULL, "width of background border"}
	  ,{ "geom", VLongRepn, 1, &geom, VOptionalOpt, geom_dict, "geometry of generated image"}
	  ,{ "r",  VFloatRepn, 1, &ball_radius, VOptionalOpt, NULL, "relative radius of ball"}
	  ,{ "outer_r",  VFloatRepn, 1, &shell_outer_radius, VOptionalOpt, NULL, "relative outer radius of shell"}
	  ,{ "inner_r",  VFloatRepn, 1, &shell_inner_radius, VOptionalOpt, NULL, "relative inner radius of shell"}
	  ,{ "ell_x",  VFloatRepn, 1, &ell_z, VOptionalOpt, NULL, "x scaling of ellipsoid"}
	  ,{ "ell_y",  VFloatRepn, 1, &ell_y, VOptionalOpt, NULL, "y scaling of ellipsoid"}
	  ,{ "ell_z",  VFloatRepn, 1, &ell_z, VOptionalOpt, NULL, "z scaling of ellipsoid"}
	  ,{ "c",  VLongRepn, 1, &c, VOptionalOpt, NULL, "corner dimension" }
	  ,{ "unit_coords",  VBooleanRepn, 1, &unit_coords, VOptionalOpt, NULL, "scale voxel to fit image into unit cube [0,1]^3" }
	  ,{ "bg",           VBooleanRepn, 1, &bg,          VOptionalOpt, NULL, "use uniform background 0 also for non-border voxels"}
	  ,{ "scale",  VFloatRepn, 1, &scale, VOptionalOpt, NULL, "geometry dimension" }
	  ,{ "max",  VLongRepn, 1, &maxval, VOptionalOpt, NULL, "maximal voxel value" }
	  ,{ "x0",  VFloatRepn, 1, &(mincoord[0]), VOptionalOpt, NULL, "minimum x coord"}
	  ,{ "y0",  VFloatRepn, 1, &(mincoord[1]), VOptionalOpt, NULL, "minimum y coord"}
	  ,{ "z0",  VFloatRepn, 1, &(mincoord[2]), VOptionalOpt, NULL, "minimum z coord"}
	  ,{ "x1",  VFloatRepn, 1, &(maxcoord[0]), VOptionalOpt, NULL, "maximum x coord"}
	  ,{ "y1",  VFloatRepn, 1, &(maxcoord[1]), VOptionalOpt, NULL, "maximum y coord"}
	  ,{ "z1",  VFloatRepn, 1, &(maxcoord[2]), VOptionalOpt, NULL, "maximum z coord"}
	  ,{ "tx",  VLongRepn, 1, &tx,            VOptionalOpt, NULL, "x translation of geom (in voxels)"}
	  ,{ "ty",  VLongRepn, 1, &ty,            VOptionalOpt, NULL, "y translation of geom (in voxels)"}
	  ,{ "tz",  VLongRepn, 1, &tz,            VOptionalOpt, NULL, "z translation of geom (in voxels)"}
	  ,{ "cliprad",  VFloatRepn, 1, &clip_radius, VOptionalOpt, NULL, "clipping sphere radius (centered at image center)"}
	};

        VParseFilterCmd(VNumber(options), options, argc, argv, NULL, &outf);
        
	if(nx == -1) nx = n;
	if(ny == -1) ny = n;
	if(nz == -1) nz = n;

	int cx, cy, cz;
	cx = cy = cz = c;
	  
        printf("Core image: nx=%i ny=%i nz=%i\n", nx, ny, nz);

	int ex_nx = nx + 2*border;
	int ex_ny = ny + 2*border;
	int ex_nz = nz + 2*border;
        printf("Extended image: ex_nx=%i ex_ny=%i ex_nz=%i\n", ex_nx, ex_ny, ex_nz);
	  

        dst =  VCreateImage(ex_nz, ex_ny, ex_nx, VUByteRepn);

	float voxelx = 1.0;
	float voxely = 1.0;
	float voxelz = 1.0;
	if(unit_coords) {
	  voxelx = voxely = voxelz = 1.0/mymax(ex_nx, ex_ny, ex_nz);
	}
	if(scale > 0.0) {
	    voxelx *= scale;
	    voxely *= scale;
	    voxelz *= scale;
	}
	if(mincoord[0] < maxcoord[0]) {
	  voxelx = (maxcoord[0]-mincoord[0])/ex_nx;
	  voxely = (maxcoord[1]-mincoord[1])/ex_ny;
	  voxelz = (maxcoord[2]-mincoord[2])/ex_nz;
	}
	char voxel_string[100];
	sprintf(voxel_string, "%f %f %f", voxelx, voxely, voxelz);
	VSetAttr(VImageAttrList(dst), "voxel", NULL, VStringRepn, voxel_string);
	
	fprintf(stdout, "voxel: %f %f %f\n", voxelx, voxely, voxelz);

	int band, row, col;
        /* set destination image entries to one*/
	if(!bg && (geom != contball || geom != contellipsoid))
	  // in principle, don't do this for any geom covering the entire range within border
	  for(band=border; band<nz+border; band++)
	    for(row=border; row<ny+border; row++) 
	      for(col=border; col<nx+border; col++)
		VPixel(dst, band, row, col, VUByte)= maxval/2; //1;
	
	if(geom == uniform) {
	  // nothing to do 
	}       
	else  if(geom == corner) { 
	  // set c x c x c corner voxels to 2
	  for(int z=border; z < cz+border; z++)
	    for(int y=border; y < cy+border; y++) 
	      for(int x=border; x < cx+border; x++)
		VPixel(dst, z, y, x, VUByte)= maxval; //2;
	}
	else if (geom == ball) {
         int cx = ex_nx/2, cy = ex_ny/2, cz = ex_nz/2;
	 float r2 = ball_radius * (mymin(cx,cy,cz)-border);
	 r2 = r2*r2;
	 for(band=border; band<nz+border; band++)
	   for(row=border; row<ny+border; row++) 
	     for(col=border; col<nx+border; col++)
	       if(r2 >= (float)((band-cz)*(band-cz) + (row-cy)*(row-cy) + (col-cx)*(col-cx)))
		 VPixel(dst, band, row, col, VUByte)=2;
	}
	else if (geom == contball || geom == contellipsoid) {
         int cx = ex_nx/2, cy = ex_ny/2, cz = ex_nz/2;
	 int min_n = mymin(nx,ny,nz);
	 float rmax2 = ball_radius*ball_radius*(min_n *min_n); // nx*nx+ny*ny+nz*nz);
	 for(band=border; band<nz+border; band++)
	   for(row=border; row<ny+border; row++) 
	     for(col=border; col<nx+border; col++) {
	       float x_scaled = (col -cx)/ell_x;
	       float y_scaled = (row -cy)/ell_y;
	       float z_scaled = (band-cz)/ell_z;

	       // float rr = (float)( (band-cz)*(band-cz) + (row-cy)*(row-cy) + (col-cx)*(col-cx));
	       float rr = x_scaled*x_scaled + y_scaled*y_scaled + z_scaled*z_scaled;
	       int band_t = band + tz;
	       int row_t  = row  + ty;
	       int col_t  = col  + tx;
	       if( 0 <= band_t && band_t < ex_nz && 0 <= row_t && row_t < ex_ny && 0 < col_t && col_t < ex_nx)
		 VPixel(dst, band_t, row_t, col_t, VUByte)
		   = (unsigned char)( (rmax2 -rr) > 0.0 ? (rmax2 -rr)*255.0/rmax2 : 0.0); 
	     }
	}

	else if (geom == shell) {
         int cx = ex_nx/2, cy = ex_ny/2, cz = ex_nz/2;
	 float r_outer = shell_outer_radius * (mymin(cx,cy,cz)-border);
	 r_outer = r_outer * r_outer;
	 float r_inner = shell_inner_radius * (mymin(cx,cy,cz)-border);
	 r_inner = r_inner * r_inner;
	 for(band=border; band<nz+border; band++)
	   for(row=border; row<ny+border; row++) 
	     for(col=border; col<nx+border; col++) {
	       float rr = ((band-cz)*(band-cz) + (row-cy)*(row-cy) + (col-cx)*(col-cx));
	       if(rr >  r_outer) 
		 VPixel(dst, band, row, col, VUByte)=0;
	       else if( rr >= r_inner)
		 VPixel(dst, band, row, col, VUByte)=maxval/2; //1;
	       else 
		 VPixel(dst, band, row, col, VUByte)=maxval; // 2;
	     }
	}
	else if(geom == two_halves){
	  if(two_halves_sep == -1)
	    two_halves_sep = ex_nx/2;
	  for(int z=border; z<nz+border; z++)
	    for(int y=border; y<ny+border; y++) 
	      for(int x=border; x<nx+border; x++)
		VPixel(dst, z, y, x, VUByte)= (x < two_halves_sep ? 44 : 88);
 
	}
	else if(geom == oblique) {
	  if(two_halves_sep == -1)
	    two_halves_sep = (int) (0.5*(nx+ny+nz+3*border));
	  for(int z=border; z<nz+border; z++)
	    for(int y=border; y<ny+border; y++) 
	      for(int x=border; x<nx+border; x++)
		VPixel(dst, z, y, x, VUByte)= (x+y+z < two_halves_sep ? 44 : 88);
	}
	else if(geom == oblique2d) {
	  if(two_halves_sep == -1)
	    two_halves_sep = 2*border+nx;
	  for(int z=border; z<nz+border; z++)
	    for(int y=border; y<ny+border; y++) 
	      for(int x=border; x<nx+border; x++)
		VPixel(dst, z, y, x, VUByte)= (x+y  < two_halves_sep ? 44 : 88);
	}
	else if(geom == brick) {
	  // set image to 2 inside centered brick [nx/4, nx3/4] x [ny/4, ny3/4] x [nz/4, nz3/4]
	  for(int z = border + nz/4; z < border + (nz*3)/4; ++z)
	    for(int y = border + ny/4; y < border + (ny*3)/4; ++y)
	      for(int x = border + nx/4; x < border + (nx*3)/4; ++x)
		VPixel(dst, z-tz, y-ty, x-tx, VUByte) = maxval; //2;
	}
	else if(geom == checker8) {
	  if(maxval == 2) maxval = 4;
	  for(int z=border; z<nz+border; z++)
	    for(int y=border; y<ny+border; y++) 
	      for(int x=border; x<nx+border; x++) {
		int val = 1 
		  + (x-border < nx/2 ? 0 : maxval/4) 
		  + (y-border < ny/2 ? 0 : maxval/2)
		  + (z-border < nz/2 ? 0 : maxval);
		VPixel(dst, z, y, x, VUByte) = val;
	      }
	}
	else if(geom == tetra) {
	  typedef float c[3];
	  plane_eval tetras;
	  for(int z=border; z<nz+border; z++)
	    for(int y=border; y<ny+border; y++) 
	      for(int x=border; x<nx+border; x++) {
		c X = { (x-border-nx/2.0)*1.0/nx, (y-border-ny/2.0)*1.0/ny, (z-border-nz/2.0)*1.0/nz};
		VPixel(dst, z, y, x, VUByte) = 1 + tetras(X); 
		
	      }
	}
	else {
	  fprintf(stderr, "geom %d unknown!\n",geom);
	  exit(1);
	}
	if(clip_radius > 0) {
	  int cx = ex_nx/2, cy = ex_ny/2, cz = ex_nz/2;
	  float r2 = clip_radius * (mymin(cx,cy,cz)-border);
	  r2 = r2*r2;
	  for(band=border; band<nz+border; band++)
	    for(row=border; row<ny+border; row++) 
	      for(col=border; col<nx+border; col++)
		if(r2 <= (float)((band-cz)*(band-cz) + (row-cy)*(row-cy) + (col-cx)*(col-cx)))
		  VPixel(dst, band, row, col, VUByte)=0;
	}

        /*write out the resulting image*/
        VWriteImages(outf, NULL, 1, &dst) ; 
        return 0;
        
}
