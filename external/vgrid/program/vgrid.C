/*
  © Copyright 2003, C&C Research Laboratories, NEC Europe Ltd.
  © Copyright 2003, Max-Planck-Institute of Cognitive Neuroscience
  

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


/*! \file
  
    Generate FE meshes from labelled volume datasets

    $Id: vgrid.C,v 1.25 2004/09/13 12:37:30 berti Exp $

*/



#include <vista.h>
#include "generateGrid.h"


/*****************************************************************************/

static VDictEntry cons_dict[] = {
   	{ "no",   	no },
	{ "box",   	box },
	{ "plane",   	plane },
   	{ "snodes", 	snodes },
   	{ NULL }
};
static VDictEntry elem_dict[] = {
        { "cube",   	cube },
        { "tetra5",	tetra5 },
        { "tetra6a", 	tetra6a },
        { "tetra6b", 	tetra6b },
        { NULL }
};
static VDictEntry out_dict[] = {
        { "vm",   	vm },
        { "hf",		hf },
        { "ipe", 	ipe },
	{ "ascii",      ascii},
        { "gmv",        gmv},  // GMV format -> 
	{ "gpp",        gpp},  // GeoFEM/GPPView format -> http://geofem.tokyo.rist.or.jp/download_en/gppview.html
	{ "complex2d",  complex2d},
	{ "dx",         dx},
        { NULL }
};
static VDictEntry smooth_dict[] = {
        { "no",		none },
        { "shift",	bshift },
        { "marching",	marching },
        { NULL }
};

static char const* get_elem_type(int t)
{
  switch(t) {
  case cube:    return "cube"; break;
  case tetra5:  return "tetra5"; break;
  case tetra6a: return "tetra6a"; break;
  case tetra6b: return "tetra6b"; break;
  default: return "";
  }
}

/*****************************************************************************/
int main (int argc, char *argv[])
{

	VAttrList list;
	VAttrListPosn posn;
	VImage src;
        FILE *inf, *outf;
		  
	/* default parameter settings */
        static VLong olevel = 1;
        static primitiveType em = cube;
	static smoothType boundary = none;
	
	static VFloat shift = 0.49;
        static VLong min = 4, max = 4;
	static VBoolean ip = false, surf=false;
	static VUByte np=1;
	static VUByte consmode=no;
	static formType form = hf;
	static VBoolean explicit_links = false;	
	static VBoolean use_simbio_mat_ids = false;  
        static VString mat_file_nm = "";
        static VString mat_db_nm   = "";
	static VLong  smooth_iters  = 0;
	static VString surface_smooth_weight_file_nm = "";
	static VString  volume_smooth_weight_file_nm = "";
	static VBoolean nofields = false;

   static VOptionDescRec options[] = {
      { "material", VStringRepn, 1, &mat_file_nm, VOptionalOpt, NULL, "name of material file"},
      { "elem", 	VLongRepn, 1, &em, VOptionalOpt, elem_dict, "Primitive type" },
      { "min", 	VLongRepn, 1, &min, VOptionalOpt, NULL, "Minimal grid resolution" },
      { "max", 	VLongRepn, 1, &max, VOptionalOpt, NULL, "Maximal grid resolution" },
      { "olevel", VLongRepn, 1, &olevel, VOptionalOpt, NULL, "verbosity level" },
      { "smooth", VLongRepn, 1, &boundary, VOptionalOpt, smooth_dict, "Smooth boundaries" },
      { "shift", VFloatRepn, 1, &shift, VOptionalOpt, NULL, "Degree of node shifting (range 0-0.49)" },
      { "format", 	VLongRepn, 1, &form, VOptionalOpt, out_dict, "Output format" },
      { "surface",VBooleanRepn, 1, &surf, VOptionalOpt, NULL, "Surface mesh output" },
      { "np", VUByteRepn, 1, &np, VOptionalOpt, NULL, "# of partitions"},
      { "constraint", VUByteRepn, 1, &consmode, VOptionalOpt, cons_dict, "Constraint Types" },
      { "explicit_links", VBooleanRepn, 1, &explicit_links, VOptionalOpt, NULL, 
	"Write explicit vertex-vertex connectivity"},
      { "simbio_mat_ids", VBooleanRepn, 1, &use_simbio_mat_ids, VOptionalOpt, NULL, 
	"Use the SimBio material IDs defined in D1.2b page 6"},
      { "matdb", VStringRepn, 1, &mat_db_nm, VOptionalOpt, NULL, 
        "name of material DB file, containing data used by -simbio_mat_ids"},
      { "sm_iters",  VLongRepn,  1, &smooth_iters,  VOptionalOpt, NULL, "Smooth marching tet interface if > 0" },
      { "sm_weights", VStringRepn, 1, &surface_smooth_weight_file_nm, VOptionalOpt, NULL, 
        "name of file containing surface smoothing weights for Laplacian smoothing"},
      { "vsm_weights", VStringRepn, 1,&volume_smooth_weight_file_nm, VOptionalOpt, NULL, 
        "name of file containing volume smoothing weights for Laplacian smoothing"},
      { "nofields",  VBooleanRepn, 1, &nofields, VOptionalOpt, NULL, "do not output any fields, except material"},
   };

        VParseFilterCmd(VNumber(options), options, argc, argv, &inf, &outf);

        if (olevel > 0) {
        	fprintf(stdout," --- VGrid --- \n");
        	fflush(stdout);
        }


	// Read source image(s):
	if (! (list = VReadFile(inf, NULL))) return 1;

	// check input parameters:
	if(surf == true)
		boundary = marching;
	if(form == ipe)
		explicit_links=true;
	if (boundary == bshift)  {
		if (em != cube)  {
			fprintf(stderr, "WARNING: vgrid: setting element type to cubes.\n");
			em = cube;
		};
		if (shift < 0 || shift > 0.5)  {
			fprintf(stderr, "WARNING: vgrid: incorrect value for shifting factor!\n");
			return 0;
		};
	};
	if( smooth_iters > 0 && ! explicit_links) {
	  fprintf(stderr, "WARNING: vgrid: setting explicit_links to true for smoothing marching tet interface!\n");
	  explicit_links = true;
	}
	if( boundary == marching && em == cube) {
	  fprintf(stderr, "WARNING: vgrid: setting elem to tetra5 for marching tets!\n");
	  em = tetra5;
	}

	// Operate on each source image:
	int nimg = 0;
	for (VFirstAttr(list, &posn); VAttrExists(&posn); VNextAttr(&posn)) {
	  if (VGetAttrRepn(&posn) == VImageRepn)  {
	    // this is our source image
	    VGetAttrValue(&posn, NULL, VImageRepn, &src);
	    ++nimg;

	    if (olevel > 0) {
	      if(boundary == marching)
		fprintf(stdout," element type = 1 \n");
	      else
		fprintf(stdout," element type = %s \n",get_elem_type(em));
	      fflush(stdout);
	    }

	    std::string mat_db_file        = mat_db_nm;	
	    std::string surface_smooth_weight_file = surface_smooth_weight_file_nm;
	    std::string volume_smooth_weight_file  =  volume_smooth_weight_file_nm;
	    gridGenerator g(src, 
			    em, 
			    min, max, 
			    boundary, 
			    shift, 
			    ip, 
			    surf, 
			    np, 
			    olevel, 
			    consmode,
			    form,
			    mat_file_nm,
			    !explicit_links,
			    use_simbio_mat_ids,
			    mat_db_file,
			    smooth_iters,
			    surface_smooth_weight_file,
			    volume_smooth_weight_file);

	    if(nofields)
	      g.noOutputFields();
	    g.work(outf);

	  }
	}
	if(nimg == 0) {
	  fprintf(stderr, "ERROR: input file does not contain an image!\n");
	  exit(1);
	}


        if (olevel > 0) {
          fprintf(stdout," --- done VGrid  --- \n");
          fflush(stdout);
        }

	return 0;
}

