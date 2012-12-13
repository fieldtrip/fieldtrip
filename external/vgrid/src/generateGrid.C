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


/*
 *  $Id: generateGrid.C,v 1.62 2005/01/26 12:07:02 berti Exp $
 */

#include "generateGrid.h"
#include "vec3.h"
#include "mesh.h"
#include "primitive.h"

#include <vista.h>


#include <math.h>
#include <assert.h>
#include <stdio.h>

#include <string>  
#include <fstream>


using namespace std;  

/***********************************************************************
      
                           helper functions 
           
***********************************************************************/

tetra::edgevec tetra::edgetable[6] 
= { {0,1}, {2,3}, {0,2}, {1,3}, {0,3}, {1,2}}; 


inline double det3(double j[3][3])
{
  return 
    j[0][0]*(j[1][1]*j[2][2]-j[2][1]*j[1][2])
    - j[1][0]*(j[0][1]*j[2][2]-j[2][1]*j[0][2])
    + j[2][0]*(j[0][1]*j[1][2]-j[1][1]*j[0][2]);

}


#define NPLACES	     1			// number of places in vtx cache
#define MAX_LABEL    256


static const int off_x[] = { 0, 1, 0, 1, 0, 1, 0, 1 };
static const int off_y[] = { 0, 0, 1, 1, 0, 0, 1, 1 };
static const int off_z[] = { 0, 0, 0, 0, 1, 1, 1, 1 };


extern unsigned int countNodes(VGraph graph);

static int mapImageToUbyteRepn(VImage& img, int verbosity) 
{
  map<int,unsigned char> newvoxel;
  VRepnKind repn = VPixelRepn(img);
  if(repn == VUByteRepn) {
    return 0;
  }
  else { // map data to reduced repn VUByteRepn
    int nx = VImageNColumns(img); 
    int ny = VImageNRows(img); 
    int nz = VImageNBands(img);
    VImage new_img = VCreateImage(nz,ny,nx, VUByteRepn);
 
    // assuming repn fits into an int.
    map<int, int> hist;
    for(int band=0; band<nz; band++)
      for(int row=0; row<ny; row++) 
	for(int col=0; col<nx; col++) {
	  int p;
	  switch(repn) {
	  case VShortRepn:
	    p = VPixel(img, band, row, col, VShort); break;
	  case VLongRepn:
	    p = VPixel(img, band, row, col, VLong);  break;
	  default:
	    fprintf(stderr, "ERROR: cannot handle representation %d\n", repn);
	    return 1;
	  }
	  if(hist.find(p) == hist.end())
	    hist[p] = 1;
	  else
	    hist[p]++;
	}
    if(hist.size() > 256) 
      fprintf(stderr, "WARNING: to many labels: %d, allowed are %d\n", hist.size(), 256);
    
    unsigned char cnt = 1;
    for(map<int,int>::iterator it = hist.begin(); it != hist.end(); ++it) {
      newvoxel[it->first] = cnt;
      cnt++;
    }
    
    if(hist[0] > 0)
      newvoxel[0] = 0; // preserve special meaning of 0 as background.

    if(verbosity > 0) {
      fprintf(stderr, "Mapping voxel values:\n");
      for(map<int,unsigned char>::const_iterator it = newvoxel.begin(); it != newvoxel.end(); ++it)
	fprintf(stderr, "%d -> %d\n", it->first, it->second);
    }

    for(int band=0; band<nz; band++)
      for(int row=0; row<ny; row++) 
	for(int col=0; col<nx; col++) {
	  switch(repn) {
	  case VShortRepn:
	    VPixel(new_img, band, row, col, VUByte) = newvoxel[VPixel(img, band, row, col, VShort)]; break;
	  case VLongRepn:
	    VPixel(new_img, band, row, col, VUByte) = newvoxel[VPixel(img, band, row, col, VLong)];  break;
	  default:
	    fprintf(stderr, "ERROR: cannot handle representation %d\n", repn);
	    exit(1);
	  }
   	} 

    // copy attributes and swap images
    VImageAttrList(new_img) = VCopyAttrList(VImageAttrList(img));
    VDestroyImage(img);
    assert(VPixelRepn(new_img) == VUByteRepn);
    img = new_img;
  }
  return 0;
}



/***********************************************************************
 ***********************************************************************
      
         class  gridGenerator member functions 

***********************************************************************          
***********************************************************************/

/*
  top level working method: Construct octree, generate mesh, and output mesh to file
*/

bool gridGenerator::work(FILE *outf)
{
  // step 1: allocate grid structures
  // allocate vertex graph
  vertices = VCreateGraph(10, 5, VFloatRepn, false);

  if (vertices == 0)  {
    fprintf(stderr, "ERROR: generateGrid: cannot create vertex graph.\n");
    return false;
  };
	
  // step 2: build an octree representation of the source image
  if (generateCells() == false) return false;

  // step 3: generate a grid from this octree
  if (generateGrid() == false) return false;
	
  // step 4: provide smooth boundaries
  if(surf_arg)
    faces = VCreateGraph(10, 4, VLongRepn, false);

  if (btype == bshift && prim_type == cube) {
    if (shiftBoundaryVertices() == false) return false;
  } else if (btype == marching) {
    if (marchingTetras() == false) return false;
  }
  // step 5: finally add elements to primitives graph
  if(surf_arg == false) {
    if( addPrimitives() == false) return false;
  }
  else
    if ( removeNonSurfaceVertices() == false) return false;
	    
  // step 6: Iterative smoothing after marching tetrahedra algorithm
  if( btype == marching)
    if (smoothMarchingInterfaces() == false) return false;

  // step 7: save results in a file
  return write(outf);
}



/************************************************************************

                   local classes of gridGenerator

***********************************************************************/

string gridGenerator::meshing_type::names[gridGenerator::meshing_type::size_]
= { string("undefined"), string("default"), string("none"), 
    string("uniform_cube"),  string("uniform_tetra5"), string("nonuniform_tets"), 
    string("nonuniform_hybrid"), string("nonuniform_cube") };


void gridGenerator::material_description::read(istream& in)
{ 
  string nm =  "material_description";
  in >> name;              REQUIRE_ALWAYS(in.good(), nm << ": reading 'name' field failed",1);
  in >> maxdim;            REQUIRE_ALWAYS(in.good(), nm << ": reading 'maxdim' field failed",1);
  in >> low_voxel;         REQUIRE_ALWAYS(in.good(), nm << ": reading 'low voxel' field failed",1);
  in >> high_voxel;        REQUIRE_ALWAYS(in.good(), nm << ": reading 'high voxel' field failed",1);
  in >> (int&) mesh_type;  REQUIRE_ALWAYS(in.good(), nm << ": reading 'mesh type' field failed",1);
  in >> weight;            REQUIRE_ALWAYS(in.good(), nm << ": reading 'weight' field failed",1);
}

void gridGenerator::interface_description::read(istream& in)
{ 
  string nm =  "interface_description";
  in >> mat1;   REQUIRE_ALWAYS(in.good(), nm << ": reading 'mat1 name' field failed",1);
  in >> mat2;   REQUIRE_ALWAYS(in.good(), nm << ": reading 'mat2 name' field failed",1);
  in >> maxdim; REQUIRE_ALWAYS(in.good(), nm << ": reading 'maxdim' field failed",1);
}





/************************************************************************

           Constructors, destructor and initalization stuff

***********************************************************************/

gridGenerator::gridGenerator(VImage s, 
			     primitiveType g, 
			     int min, int max,
			     smoothType bt, 
			     double sf, 
			     bool ip, 
			     bool surf, 
			     VUByte np,
			     int olevel,
			     VUByte consmode,
			     formType form,
			     VString mat_file,
			     bool implicit_links_,
			     bool use_simbio_mat_ids_,
			     std::string const& mat_db_file,
			     int   smooth_iters_,
			     std::string const& surface_smooth_weight_file_,
			     std::string const& volume_smooth_weight_file_)
  : mat_file_nm(mat_file),
    MatDB(mat_db_file)
{
  src = s; prim_type = g; 
  if(min <= 0) {
    fprintf(stderr, "WARNING: illegal value %d for min, setting to 1!\n",min);
    min = 1;
  }
  if(max < min) {
    fprintf(stderr, "WARNING: value %d for max is < min, setting to value of min=%d!\n",max,min);
    max = min;
  }
  mindim = min; maxdim = max;
  // now 0 < mindim <= maxdim

  btype = bt; shift = sf; interp = ip; surf_arg=surf; 

  npart=np;  
  if(npart <= 0) {
    fprintf(stderr, "WARNING: npart=%d, setting npart to 1!\n",npart);
    npart = 1;
  }
  out_arg=form;
  nx = VImageNColumns(src); ny = VImageNRows(src); nz = VImageNBands(src);

  verbose = olevel;
  constype=consmode;
  implicit_links = implicit_links_;
  use_simbio_mat_ids = use_simbio_mat_ids_;
  do_output_fields = true;

  // stuff for Laplacian smoothing
  smooth_iters  = smooth_iters_;
  surface_smooth_weight_file_nm =surface_smooth_weight_file_;
  volume_smooth_weight_file_nm =volume_smooth_weight_file_;
  read_smoothing_weights(surface_smooth_weight_file_nm, surface_smooth_weights);
  read_smoothing_weights( volume_smooth_weight_file_nm,  volume_smooth_weights);

  oct = 0; octnodes = 0;
  vertices = 0; primitives = 0; faces = 0;

  bool repn_mismatch = (VPixelRepn(src) != VUByteRepn);
  if(repn_mismatch) {
    fprintf(stderr, "WARNING: Pixel representation will be converted to unsigned byte!\n");
    int err = mapImageToUbyteRepn(src, verbose);
    if(err) 
      exit(err);
  }

  initMaterials();

  if(surf_arg == true)
    if(btype != marching) {
      btype = marching;
      fprintf(stderr, "WARNING: setting boundary smoothing type to marching tets for surface mesh generation!\n");
    }

  if(btype == marching) {
    // currently only for two material
    // either exactly 2 materials with arb. meshing type,
    // or more materials, but only one of which gets meshed
    mat_background  = 0;
    mat_body        = 1; // 
    if ( num_of_materials() > 2) {
      unsigned num_meshed_materials = 0;
      for(size_t m = 0; m < material_list.size(); ++m) {
	if(material_list[m].mesh_type != meshing_type::none) {
	  ++ num_meshed_materials;
	  mat_body = m;
	}
	else mat_background = m;
      }
      if(num_meshed_materials > 1) {
	fprintf(stderr, "WARNING: Marching tets do work only for 2 materials,skipping it.\n");
	btype = none;
      }
    }
  }
	  
  // adapt mindim, maxdim and image sizes 
  // maximal cell size must be 2^k * minimal cell size
  int new_max = min_maxdim;
  while(new_max < max_maxdim)
    new_max *= 2;
  if(new_max != max_maxdim)
    fprintf(stderr, "WARNING: min_maxdim = %d, setting max_maxdim to %d\n", min_maxdim, new_max);
  max_maxdim = new_max;

       
  bool size_mismatch = 
    (  (max_maxdim *(nx/max_maxdim) != nx)
       || (max_maxdim *(ny/max_maxdim) != ny)
       || (max_maxdim *(nz/max_maxdim) != nz));
  // Enlarge image if necessary
  // check if images sizes are multiple of largest possible cell size
  if(size_mismatch) {
    // Enlarge image to new size
    int nx1 = max_maxdim * (nx/max_maxdim + 1);
    int ny1 = max_maxdim * (ny/max_maxdim + 1);
    int nz1 = max_maxdim * (nz/max_maxdim + 1);

    fprintf(stderr,"WARNING: max_maxdim must divide nx,ny,nz!\n");
    fprintf(stderr,"WARNING: enlarging image from size %d, %d, %d to %d,%d,%d\n", nx,ny,nz,nx1,ny1,nz1); 
	  
    VImage nsrc = VCreateImage(nz1, ny1, nx1, VUByteRepn);
    for(int band=0; band<nz1; band++)
      for(int row=0; row<ny1; row++) 
	for(int col=0; col<nx1; col++)
	  VPixel(nsrc, band, row, col, VUByte) = 0; // background
    for(int band=0; band<nz; band++)
      for(int row=0; row<ny; row++) 
	for(int col=0; col<nx; col++)
	  VPixel(nsrc, band, row, col, VUByte) = VPixel(src, band, row, col, VUByte);
    nx = nx1; ny = ny1; nz = nz1;

    VImageAttrList(nsrc) = VCopyAttrList(VImageAttrList(src));
    VDestroyImage(src);

    src = nsrc;
    // if 0 has not occured in the original image, add it.
	  
    if(materials[0] == undefined_material) {      
      materials[0] = material_type(num_of_materials_);
      ++num_of_materials_;
      material_list.push_back(material_description("background",
						   maxdim, 0, 0,
						   meshing_type::none,
						   1.0));
    }
          
  } // enlarge image


  // init vertex cache
  initCache();

  if(prim_type != tetra5 && nonuniform_meshing()) { 
    fprintf(stderr, "WARNING: gridGenerator: Changing grid generation to 'tetra5' \n");
    prim_type = tetra5;
  }

}

/*
  destructor
*/

gridGenerator::~gridGenerator()
{
  deleteCache();
  if (src) 		VDestroyImage(src);        src        = 0;
  if (oct) 		VDestroyImage(oct);        oct        = 0;
  if (octnodes) 		VDestroyImage(octnodes);   octnodes   = 0;
  if (vertices) 		VDestroyGraph(vertices);   vertices   = 0;
  if (primitives) 	VDestroyGraph(primitives); primitives = 0;

  if(cache != 0)
    deleteCache();
}



/*
  Vertex cache
*/

void gridGenerator::initCache()
  // setup a vertex cache
{
  int size = (nx+1)*(ny+1)*(nz+1);
  cache = new vertex* [size]; 
  nv = 0;
  for (int i = 0; i < size; i++) 
    cache[i] = 0;
}



void gridGenerator::deleteCache()
  // delete vertex cache
{
  int size = (nx+1)*(ny+1)*(nz+1);
  for (int i = size-1; i >= 0; i--) 
    delete [] cache[i];
  delete [] cache; cache = 0; 
  nv = 0;
}


void gridGenerator::read_smoothing_weights(std::string const& filenm, std::vector<double> & weights)
{
  double default_weight = 0.5;
  if(filenm == "") {
    weights.push_back(default_weight);
    return;
  }
  ifstream weight_file(filenm.c_str());
  if(! weight_file.is_open()) {
    if(filenm != "") {
      fprintf(stderr, "WARNING: could not open file \"%s\", using default weight!\n", 
	      filenm.c_str());
    }
    weights.push_back(default_weight);
    return;
  }
  double w = 0.0;
  while(weight_file >> w)
    weights.push_back(w);

  if(weights.empty()) {
    fprintf(stderr, "WARNING: no weights found in file \"%s\"\n, using default weight!\n",
	    filenm.c_str());
    weights.push_back(default_weight); 
  }
}

/*****************************************************************************

                            material handling 

****************************************************************************/

gridGenerator::material_type 
gridGenerator::get_material(string const& name) const
{
  for(size_t i = 0; i < material_list.size(); ++i)
    if(material_list[i].name == name)
      return material_type(i);
  fprintf(stderr, "WARNING: undefined material %s detected!\n", name.c_str());
  return undefined_material;
}


string const&
gridGenerator::material_name(gridGenerator::material_type m) const
{
  REQUIRE_ALWAYS( (0 <= m && m < (material_type)num_of_materials()), "m = " << m,1);
  return material_list[m].name;
}


void gridGenerator::print_material_parameters() const 
{
  if(material_list.size() > 0)
    fprintf(stdout, "Materials specified in file %s:\n", mat_file_nm.c_str());
  for(size_t m = 0; m < material_list.size(); ++m) {
    material_description mat = material_list[m];
    fprintf(stdout, "%s \t internal no.: %d \t maxdim: %d \t voxel range: [%d, %d] \t meshing type: %s\n", 
	    mat.name.c_str(),  m,   mat.maxdim, mat.low_voxel, mat.high_voxel, 
	    meshing_type::names[mat.mesh_type].c_str());
  }
  if(maxdim_of_material.size() > material_list.size())
    fprintf(stdout, "Material determined from voxel values:\n");
  for(size_t m = material_list.size(); m < maxdim_of_material.size(); ++m)
    fprintf(stdout,"No. %d \t maxdim: %d \t algorithm: %d\n",
	    m, maxdim_of_material[m], meshing_algorithm[m]);
}






void gridGenerator::initMaterials() 
{
  vector<int> hist(MAX_LABEL, 0);
  materials.resize(MAX_LABEL, undefined_material);
  num_of_materials_ = 0;

  //--- (1) get actual materials present in image --------

  for(int x = 0; x < nx; ++x)
    for(int y = 0; y < ny; ++y)
      for(int z = 0; z < nz; ++z)
	hist[getCellVoxel(x,y,z)]++;

  // (1A) read voxel->material from file
  if(mat_file_nm != "") {
    ifstream material_file(mat_file_nm.c_str());
    if(! material_file.is_open()) {
      fprintf(stderr, "Could not open material file %s\n",mat_file_nm.c_str());
      exit(1);
    }
    if(verbose > 0)  fprintf(stdout, "Opened material file %s\n",mat_file_nm.c_str());
    material_description mat_entry;
    interface_description if_entry;
    string s;
    while (material_file >> s) {
      if(s == "material") {
	material_file >> mat_entry;
	if(mat_entry.maxdim <= 0) {
	  fprintf(stdout, "Setting maxdim of material %s to %d\n", mat_entry.name.c_str(), this->maxdim);
	  mat_entry.maxdim = this->maxdim;
	}
	material_list.push_back(mat_entry);
      }
      else if(s == "interface") {
	material_file >> if_entry;
	if(if_entry.maxdim <= 0) {
	  fprintf(stdout, "WARNING: maxdim of interface (%s,%s) is %d, adjusting to global value %d!\n", 
		  if_entry.mat1.c_str(), if_entry.mat2.c_str(), if_entry.maxdim, this->mindim);
	  if_entry.maxdim = this->mindim;
	}
	material_interfaces.push_back(if_entry);
      }
      else {
	fprintf(stderr, "Invalid key \'%s\' found!\n", s.c_str());
	exit(2);
      }
    }
    num_of_materials_ = material_list.size();
    for(size_t m = 0; m < num_of_materials(); ++m) 
      for(int vox = material_list[m].low_voxel; 
	  vox    <= material_list[m].high_voxel; ++vox)
	materials[vox] = material_type(m);

    // check whether all materials are described in file
    for(int vox = 0; vox < MAX_LABEL; ++vox) {
      if(hist[vox] > 0 && materials[vox] ==  undefined_material) {
	fprintf(stderr, "material for voxel value %d (occuring %d times) ", vox, hist[vox]);
	fprintf(stderr, "not described in material file! Using standard values.\n");
        materials[vox] = material_type(num_of_materials_);
	++num_of_materials_;
	material_list.push_back(material_description("voxel" /*+ as_string(vox)*/, 
						     maxdim, vox, vox, 
						     (vox == 0 ? meshing_type::none : meshing_type::deflt),
						     1.0));
      }
    }
  } // mat_file != ""

  // or (1B): each voxel value defines a different material
  else{
    for(int vox = 0; vox < MAX_LABEL; ++vox) {
      if(hist[vox] > 0) {
	materials[vox] = material_type(num_of_materials_);
	material_list.push_back(material_description("voxel" /*+ as_string(vox)*/, 
						     maxdim, vox, vox, 
						     (vox == 0 ? meshing_type::none : meshing_type::deflt),
						     1.0));
	++num_of_materials_;
      }
      else
	materials[vox] = undefined_material;
    }
  }

  //----  initialize with default values ----  
  maxdim_of_material .resize(num_of_materials(), 0);
  maxdim_of_interface.resize(num_of_materials(), vector<int>(num_of_materials(), 0));
  meshing_algorithm  .resize(num_of_materials(), meshing_type::undefined);


  //--- (2) Set maxdim values for materials and interfaces. ------

  // (2A) use maxdim settings read from file

  for(size_t m = 0; m < material_list.size(); ++m) {
    maxdim_of_material[m] = material_list[m].maxdim;
    meshing_algorithm [m] = material_list[m].mesh_type;
  }
  for(size_t i = 0; i < material_interfaces.size(); ++i) {
    int m1 = get_material(material_interfaces[i].mat1);
    int m2 = get_material(material_interfaces[i].mat2);
    if(m1 == undefined_material || m2 == undefined_material)
      continue;
    maxdim_of_interface[m1][m2] = material_interfaces[i].maxdim;
    maxdim_of_interface[m2][m1] = material_interfaces[i].maxdim;
  }

  // (2B) after reading specific parameters, set defaults for the rest.
  for(size_t m = 0; m < num_of_materials(); ++m) {
    if(maxdim_of_material[m] == 0){
      maxdim_of_material[m] = maxdim;
    }
    for(size_t m1 = 0; m1 < num_of_materials(); ++m1)
      if(maxdim_of_interface[m][m1] == 0 ) {
	maxdim_of_interface[m][m1] = mindim;
      }
  }

  //-------- (4)  set max_maxdim/min_maxdim -----------------

  // init to values of first meshed material 
  for(size_t m = 0; m < num_of_materials(); ++m) {
    if(material_list[m].mesh_type != meshing_type::none) {
      max_maxdim = min_maxdim = maxdim_of_material[m];
      break;
    }
  }
  // Find the max/min of all dimensions of mesh materials and their interfaces
  for(size_t m = 0; m < num_of_materials(); ++m) {
    // consider only materials which are to be meshed
    if(material_list[m].mesh_type != meshing_type::none) {
      max_maxdim = max(max_maxdim, maxdim_of_material[m]);
      min_maxdim = min(min_maxdim, maxdim_of_material[m]);
      // .. and all interfaces
      for(size_t m1 = 0; m1 < num_of_materials(); ++m1) {
	max_maxdim = max(max_maxdim, maxdim_of_interface[m1][m]);
	min_maxdim = min(min_maxdim, maxdim_of_interface[m1][m]);
      }
    }
  }
  if(verbose > 0) {
    fprintf(stdout, "min_maxdim = %d, max_maxdim = %d\n", min_maxdim, max_maxdim);
  }



  if(verbose > 0) {
    for(int vox = 0; vox < MAX_LABEL; ++vox) {
      if(hist[vox] > 0) {
        fprintf(stdout,"Voxel %d occurs %d \t times (%f %s)\n",
		vox,hist[vox], hist[vox]*100.0/(nx*ny*nz),"%");
      }
    }
    print_material_parameters();
  }
} // gridGenerator::initMaterials






/*****************************************************************************

                 member functions for octree creation

*****************************************************************************/



/*  
    Top level octree generation method

*/
bool gridGenerator::generateCells()
{
  if (verbose) {
    fprintf(stdout, "generating cells ...  nx=%d ny=%d nz=%d\n",nx,ny,nz);
    fflush(stdout);
  }

  oct =      VCreateImage(nz,   ny,   nx,   VUByteRepn);
  octnodes = VCreateImage(nz+1, ny+1, nx+1, VBitRepn);
  if (oct == 0 || octnodes == 0)  {
    fprintf(stderr, "generateCells: cannot create octree image.\n");
    return false;
  };

  VFillImage(octnodes, VAllBands, 0);
     

  // Fill octree with cells of size min_maxdim
  VFillImage(oct, VAllBands, 0);
  int d = min_maxdim;
  for (int z = 0; z < nz; z+=d)
    for (int y = 0; y < ny; y+=d)
      for (int x = 0; x < nx; x+=d)
	setCellDim(x,y,z,d);

  // Map voxels to materials
  // FIXME: use different image for materials!
  for (int z = 0; z < nz; z++)
    for (int y = 0; y < ny; y ++)
      for (int x = 0; x < nx; x ++)
	setCellMaterial(x,y,z, get_material(getCellVoxel(x,y,z)));

  // subsample smallest allowed cell size
  subsampleCells(min_maxdim);

  // create maximally admissible cells in octree
  compactCells();

  // print octree cell dim histogram
  vector<int> cellhist(max_maxdim+1,0);
  unsigned cellcnt = 0;
  for (int z = 0; z < nz; z++)
    for (int y = 0; y < ny; y++)
      for (int x = 0; x < nx; x++)
	if(getCellDim(x,y,z) > 0) {
	  cellcnt++;
	  cellhist[getCellDim(x,y,z)]++;
	}
  if(verbose) {
    fprintf(stdout,"Number of octree cells: %d\n",cellcnt);
    for(int d = 1; d <= max_maxdim; ++d)
      if(cellhist[d] > 0) 
	fprintf(stdout," dim %d: %d cells\n",d,cellhist[d]);
  }
  // insert octree cell corners (oct) into octree node table (octnodes)
  projectCorners();

  if (verbose) {
    fprintf(stdout, "\n done generating cells\n");
    fflush(stdout);
  }

  return true;
}



/*
  Get all materials present in the supervoxel [ox,oy,oz,ox+d,oy+d,oz+d]
*/

void gridGenerator::list_materials(int ox, int oy, int oz, int d, 
				   vector<unsigned> & mat) const
{
  assert(mat.size() == num_of_materials());
  assert(getCellDim(ox,oy,oz) > 0);
  assert(getCellDim(ox,oy,oz) <= d);

  // recursivly get materials.
  // the cell of dim d at ox,oy,oz may be subdivided,
  // i.e. d > getCellDim(ox,oy,oz)

  // cell is octree cell, i.e. not subdivided
  if(d == getCellDim(ox,oy,oz)) {
    mat[getCellMaterial(ox,oy,oz)] += 1;
  }
  else { // cell is subdivided: recurse
    for(int x = 0; x <= 1; ++x)
      for(int y = 0; y <= 1; ++y)
	for(int z = 0; z <= 1; ++z) {
	  list_materials(ox+x*d/2, oy+y*d/2, oz+z*d/2, d/2, mat);
	}
  }
}

/*
  get maximal allowed dimension of cell (ox,oy,oz,d)
  get the subsampled material of cell (ox,oy,oz,d)
*/
int gridGenerator::get_maxdim(int ox, int oy, int oz, int d, 
			      material_type& max_mat) const
{
  typedef vector<unsigned> mat_table;
  mat_table mat(num_of_materials(),0); // material -> # of material voxels
  list_materials(ox,oy,oz,d,mat);
  int maxd = max_maxdim; // global maxdim

  max_mat = 0; // material with maximal number of voxels
  for(size_t m =0; m < num_of_materials(); ++m) { 
    // get minimal maxdim for all materials of cell (ox,oy,oz,d) ...
    if(mat[m] > 0) {
      if(mat[m] > mat[max_mat])
	max_mat = m;
      maxd = min(maxd, get_maxdim(material_type(m)));
      // ... and all material interfaces inside cell
      for(size_t m2 = m+1; m2 < num_of_materials(); ++m2)
	if(mat[m2] > 0)
	  maxd = min(maxd, get_maxdim(material_type(m),material_type(m2)));
    }
  }
  // TODO: now consider also geometric definitions of maxdim

  
  return maxd;
}



void gridGenerator::join_cells(int x, int y, int z, int d)
{
  assert(getCellDim(x,y,z) == d/2);

  setCellDim(x, y, z, d);

  // set neighbourhood to 0
  for (int i = 1; i < 8; i++)  {
    int nd = d/2;
    int ex = x+off_x[i]*nd;
    if (ex >= nx) continue;
    int ey = y+off_y[i]*nd;
    if (ey >= ny) continue;
    int ez = z+off_z[i]*nd;
    if (ez >= nz) continue;
    setCellDim(ex, ey, ez, 0);
  }
}



/* 
   Return the value that best represents a cell at (ox, oy, oz) 
   with dimension d 
*/
int gridGenerator::subsampleCell(int ox, int oy, int oz, int d)
{
  vector<unsigned> value(num_of_materials(),0);
  
  for (int z = oz; z < nz && z < oz+d; z++) {
    for (int y = oy; y < ny && y < oy+d; y++) {
      for (int x = ox; x < nx && x < ox+d; x++)  {
	value[getCellMaterial(x,y,z)]++;
      }
    }
  }

  int mx = 0;
  for(unsigned m = 0; m < num_of_materials(); ++m)
    mx = (value[m]*material_list[m].weight  > value[mx]*material_list[mx].weight ? m : mx);
  
  return mx;
}

/*
  Subsampling:
  Generate a new source grid that is subsampled by a factor of mindim
  represent each subsampled region by the most frequent label
*/
void gridGenerator::subsampleCells(int d)
{
  if (verbose) fprintf(stdout, "subsampling cells ... ");
  for (int z = 0; z < nz; z += d) {
    // if (verbose) fprintf(stdout, "subsampling slice %03d\r", z);
    for (int y = 0; y < ny; y += d) {
      for (int x = 0; x < nx; x += d) {
	// represent this cell by a single value
	int v = subsampleCell(x, y, z, d);
	// VSetPixel(src, z, y, x, (double)v);
	setCellMaterial(x,y,z,v);
	setCellDim(x, y, z, d);
      }
    }
  }
  if (verbose) fprintf(stdout, " ... done subsampling cells\n");
}



/*
  Is the supervoxel [(x, y, z), (x+d,y+d,z+d)] homogenous, 
  i.e. composed of 8 cells of dimension d/2 with the same material? 
*/
bool gridGenerator::isHomogenousCell(int x, int y, int z, int d)

{
  // test dimensions of this cell
  int dim = getCellDim(x, y, z);
  if (dim != d/2) return false;

  // then test neighbours
  int v = (int)VPixel(src, z, y, x, VUByte);
  for (int i = 1; i < 8; i++)  {
    // take care of boundaries
    int ex = x+off_x[i]*dim;
    if (ex >= nx) continue;
    int ey = y+off_y[i]*dim;
    if (ey >= ny) continue;
    int ez = z+off_z[i]*dim;
    if (ez >= nz) continue;
		
    // compare dimension with neighbour
    if (getCellDim(ex, ey, ez) != d/2) return false;
		
    // compare value with neighbour
    if (VPixel(src, ez, ey, ex, VUByte) != v) return false;
  };
  return true;
}


/*
  compact regions in the source image up to maximal allowed cell size
  (defined by material and material interface specific resolutions)
  On output, each voxel (x,y,z) in VImage oct represents an octree cell of size
  getCellDim(x,y,z) 
*/
void gridGenerator::compactCells()

{
  // for each cell size d = 2^k * min_maxdim 
  for (int d = 2*min_maxdim; d <= max_maxdim; d *= 2) {
    // for each cell of size d
    for (int z = 0; z < nz; z += d) {
      for (int y = 0; y < ny; y += d) {
	for (int x = 0; x < nx; x += d) {
	  if(getCellDim(x,y,z) == d/2) {
	    material_type new_mat = undefined_material;
	    if(get_maxdim(x,y,z, d, new_mat) >= d) {
	      // cell (x,y,z,d) is admissible
	      join_cells(x,y,z,d);
	      setCellMaterial(x,y,z,new_mat);
	    }
	  }
	}
      }
    }
  }
}

/*
  Mark the corners of each cell in adjacent cells as level 1
  for each cell
*/
void gridGenerator::projectCorners()
{
  for (int z = 0; z < nz; z += 1) {
    for (int y = 0; y < ny; y += 1) {
      for (int x = 0; x < nx; x += 1) {
	// if cell is not marked, continue
	// FIXME: this test  results in missing nodes!
	//if(mustTesselateCell(x,y,z)) {
	int nd = getCellDim(x, y, z);
	if (nd == 0) continue;
	// else mark cell corners in octnodes
	for (int i = 0; i < 8; i++)  {
	  int ex = x+off_x[i]*nd;
	  int ey = y+off_y[i]*nd;
	  int ez = z+off_z[i]*nd;
	  if(ex > nx || ey > ny || ez > nz)
	    fprintf(stderr, 
		    "Warning: (ex,ey,ez) = (%d,%d,%d)\n",
		    ex,ey,ez);
	  else {
	    VPixel(octnodes, ez, ey, ex, VBit) = 1;
	  }
	}
	//} // if(mustTesselateCell(x,y,z))
      }
    }
  }
  // count octree nodes
  unsigned vtxcnt = 0;
  for (int z = 0; z <= nz; z += 1) {
    for (int y = 0; y <= ny; y += 1) {
      for (int x = 0; x <= nx; x += 1) {
	if(isVertex(x,y,z))
	  vtxcnt++;
      }
    }
  }
  if(verbose)
    fprintf(stdout,"Number of octree nodes: %d\n",vtxcnt);

  if(verbose >= 3) {
    // pretty-print octnodes
    fprintf(stdout, "octnodes:\n");
    for (int z = 0; z <= nz; z += mindim) {
      fprintf(stdout,"z = %d\n",z); 
      for (int y = 0; y <= ny; y += mindim) {
	fprintf(stdout, "y = %d\t:",y);
	for (int x = 0; x <= nx; x += mindim) {
	  int b = VPixel(octnodes, z, y, x, VBit);
	  fprintf(stdout,"%d ",b);
	}
	fprintf(stdout,"\n");
      }
      fprintf(stdout,"\n");
    }  
  }
}








/*****************************************************************************

               Member functions for mesh generation from octree

*****************************************************************************/

/*
  Top level mesh generation routine
*/
bool gridGenerator::generateGrid()
{
  if (verbose >= 1) {
    fprintf(stdout,"generating grid ...   (%d slices)\n",nz);
    fflush(stdout);
  }

  // now transform cells into graph
  for (int z = 0; z < nz; z++) {
    if (verbose >= 1 && verbose < 3) {
      fprintf(stdout, "working on slice %03d\r", z);
      fflush(stdout);
    }
    for (int y = 0; y < ny; y++) {
      for (int x = 0; x < nx; x++) {
	if(getCellDim(x, y, z) == 0) continue;
	if(! mustTesselateCell(x,y,z)) continue;
	tesselateCell(x, y, z);
      }
    }
  }

  if (verbose >= 1) {
    fprintf(stdout, "\n ... done generating grid.\n");
    fflush(stdout);
  }
  return true;
}



/*
  Tesselate a general octree cell
*/
void gridGenerator::tesselateCell(int ox, int oy, int oz)
{
  if (verbose > 1) {
    fprintf(stdout,"tesselateCell:  ox=%d oy=%d oz=%d \n",ox,oy,oz);
    fflush(stdout);
  }

  int d = getCellDim(ox, oy, oz);
  octcell c(*this,ox,oy,oz,d);
  //tesselated_cells.push_back(cell(*this,ox,oy,oz,d));

  // prepare point list
  // this is needed only in the case of non-uniform meshing,
  // but is used for detecting errors in the uniform case
  // (we should find at least 8 nodes, and exactly 8 in the uniform case.)
  vector<int> plist;

	  
  // given a high number of small cells it is more efficient
  // to scan the whole volume instead of just the surfaces
  int inc = (nonuniform_meshing() ? 1 : d);
  for (int z = oz; z <= nz && z <= oz+d; z += inc) {
    for (int y = oy; y <= ny && y <= oy+d; y += inc) {
      for (int x = ox; x <= nx && x <= ox+d; x += inc) {
	// is it a vertex?
	if (isVertex(x,y,z)) {
	  vertex v(x, y, z);
	  plist.push_back(addVertex(v));
	}
      }
    }
  }


  if(nonuniform_meshing()) {
    if (plist.size() == 8) {
      // use economic tesselation
      tesselateCube(c);
    }
    else {
      // more than 8 nodes => use midpoint-based tesselation
      vector<triangle> ft; // facet triangulation;
      // triangulate each facet of cell
      get_facet_triangulation(ft, ox  ,oy  ,oz  , ox  , oy+d,oz+d,d,xdir);
      get_facet_triangulation(ft, ox+d,oy  ,oz  , ox+d, oy+d,oz+d,d,xdir);
      get_facet_triangulation(ft, ox  ,oy  ,oz  , ox+d, oy  ,oz+d,d,ydir);
      get_facet_triangulation(ft, ox  ,oy+d,oz  , ox+d, oy+d,oz+d,d,ydir);
      get_facet_triangulation(ft, ox  ,oy  ,oz  , ox+d, oy+d,oz  ,d,zdir);
      get_facet_triangulation(ft, ox  ,oy  ,oz+d, ox+d, oy+d,oz+d,d,zdir);
      // connect each triangle in ft with midpoint m
      int m = addVertex(ox+d/2,oy+d/2,oz+d/2);
      tesselated_cells.push_back(octcell_tesselation(c,m));
      for(vector<triangle>::const_iterator t = ft.begin(); t != ft.end(); ++t) {
	// tesselated_cells.back().addElement(addTetra((*t).V(0), (*t).V(1), (*t).V(2), m));
	tesselated_cells.back().addTetra((*t).V(0), (*t).V(1), (*t).V(2), m, getCellMaterial(c));
	if(verbose >= 3)
	  fprintf(stdout, "[%d,%d,%d] ", (*t).V(0), (*t).V(1), (*t).V(2));
      }
    }
  }
  else { // uniform meshing

    // take the point list from above and perform a tesselation
    if (plist.size() == 8 ) { 
      // bool par = parity(ox,oy,oz,d);
      tesselateCube(c); //plist, par);
    }
    else if (plist.size() > 8) {
      // > 8 nodes should not occur for uniform meshing
      fprintf(stderr, "Cell found with %d > 8 vertices, use non-uniform meshing!\n",
	      plist.size());
    } 
    else { // < 8 nodes should never occur
      fprintf(stderr, 
	      "degenerate cell at (%d,%d,%d), dim= %d, cnt = %d (at boundary?).\n",
	      ox,oy,oz,d,plist.size());
      for(size_t i = 0; i < plist.size(); ++i) {
	vertex* v =  (vertex *)VGraphGetNode(vertices, plist[i]);
	fprintf(stderr, "[%f,%f,%f] ", v->x,v->y,v->z);
      }
      fprintf(stderr, "\n");
    }
  }
}


/*
  Tesselate a regular octree cell
*/
void gridGenerator::tesselateCube(octcell const& c) 
{
  tesselated_cells.push_back(octcell_tesselation(c));
  octcell_tesselation & tc(tesselated_cells.back());
  int p[8];
  int d = c.dim(), ox = c.x(), oy = c.y(), oz = c.z();
  material_type mat = getCellMaterial(c);
  p[0] = addVertex(ox  ,oy  ,oz  ); p[1] = addVertex(ox+d,oy  ,oz  );
  p[2] = addVertex(ox  ,oy+d,oz  ); p[3] = addVertex(ox+d,oy+d,oz  );
  p[4] = addVertex(ox  ,oy  ,oz+d); p[5] = addVertex(ox+d,oy  ,oz+d);
  p[6] = addVertex(ox  ,oy+d,oz+d); p[7] = addVertex(ox+d,oy+d,oz+d);


  bool par = parity(ox,oy,oz,d);
  switch (prim_type)  {
  case tetra5:
    // make tetrahedra for even vs. odd cells
    if (par == false)  {	// even cell
      tc.addTetra(p[0], p[1], p[2], p[4],mat);
      tc.addTetra(p[7], p[1], p[4], p[2],mat);
      tc.addTetra(p[1], p[7], p[3], p[2],mat);
      tc.addTetra(p[1], p[7], p[4], p[5],mat);
      tc.addTetra(p[2], p[6], p[4], p[7],mat);
    } else {		// odd cell
      tc.addTetra(p[4], p[5], p[6], p[0],mat);
      tc.addTetra(p[3], p[5], p[0], p[6],mat);
      tc.addTetra(p[5], p[3], p[7], p[6],mat);
      tc.addTetra(p[5], p[3], p[0], p[1],mat);
      tc.addTetra(p[6], p[2], p[0], p[3],mat);
    };
    break;
  case tetra6a:
    tc.addTetra(p[0], p[1], p[2], p[5],mat);
    tc.addTetra(p[1], p[2], p[3], p[5],mat);
    tc.addTetra(p[0], p[2], p[4], p[5],mat);
    tc.addTetra(p[2], p[3], p[5], p[7],mat);
    tc.addTetra(p[2], p[5], p[6], p[7],mat);
    tc.addTetra(p[2], p[4], p[5], p[6],mat);
    break;
  case tetra6b:
    tc.addTetra(p[0], p[1], p[2], p[4],mat);
    tc.addTetra(p[1], p[2], p[4], p[5],mat);
    tc.addTetra(p[2], p[4], p[5], p[6],mat);
    tc.addTetra(p[1], p[2], p[3], p[5],mat);
    tc.addTetra(p[2], p[3], p[6], p[5],mat);
    tc.addTetra(p[3], p[5], p[6], p[7],mat);
    break;
  case cube:
    tc.addCube(p,mat);
    break;
  case nogrid:
  default:
    fprintf(stderr, "gridGenerator::tesselateCube: no grid mode given.\n");
    break;
  };
}




/*********************************************************************
   helper types for gridGenerator::get_diagonal_triangulation
*********************************************************************/

struct int3 {
  int X[3];
  int3() {}
  int3(int x, int y, int z) { X[0] = x; X[1] = y; X[2] =z;}
  int  operator[](int i) const { return X[i];}
  int& operator[](int i)       { return X[i];}
};

// basis-vectors of facets
typedef int3 fbasis[2];
static fbasis facet_basis[3]  = {{ int3(0,1,0), int3(0,0,1) },   // xdir
				 { int3(1,0,0), int3(0,0,1) },   // ydir
				 { int3(1,0,0), int3(0,1,0) }};  // zdir


bool gridGenerator::isVertex(int x, int y, int z) const
{ return VPixel(octnodes, z,y,x,VBit) != 0;}

int gridGenerator::parity(int x, int y, int z, int d) const
{ 
  assert ( (d*(x/d) == x) && (d*(y/d) == y) && (d * (z/d) == z));
  return ((x/d) ^ (y/d) ^ (z/d)) & 0x1;
}





void gridGenerator::get_diagonal_triangulation(vector<triangle> & triangle_list,
					       int x0, int y0, int z0,
					       int x1, int y1, int z1,
					       int d,  gridGenerator::facet_dir dir) 
{

  if(verbose >= 3) {
    fprintf(stdout, "get_diag_triang (x0=%d,y0=%d,z0=%d,x1=%d,y1=%d,z1=%d,d=%d,dir=%d\n",
	    x0,y0,z0,x1,y1,z1,d,dir);
  }
  int p = parity(x0,y0,z0,d);
  int3 e1 = facet_basis[dir][0];
  int3 e2 = facet_basis[dir][1];
  int3 e3(e1[0]+e2[0], e1[1]+e2[1], e1[2]+e2[2]);
  assert ( x1 == x0 + d*e3[0]);
  assert ( y1 == y0 + d*e3[1]);
  assert ( z1 == z0 + d*e3[2]);


  if(p == 1) {
    triangle_list.push_back(triangle(addVertex(x0,y0,z0), 
				     addVertex(x0+d*e1[0],y0+d*e1[1], z0+d*e1[2]),
				     addVertex(x1,y1,z1)));
    triangle_list.push_back(triangle(addVertex(x0,y0,z0), 
				     addVertex(x0+d*e2[0],y0+d*e2[1], z0+d*e2[2]),
				     addVertex(x1,y1,z1)));
  }
  else {
    triangle_list.push_back(triangle(addVertex(x0,y0,z0), 
				     addVertex(x0+d*e1[0],y0+d*e1[1], z0+d*e1[2]),
				     addVertex(x0+d*e2[0],y0+d*e2[1], z0+d*e2[2])));
    triangle_list.push_back(triangle(addVertex(x0+d*e1[0],y0+d*e1[1], z0+d*e1[2]),
				     addVertex(x0+d*e2[0],y0+d*e2[1], z0+d*e2[2]),
				     addVertex(x1,y1,z1)));
  }
} // gridGenerator::get_diagonal_triangulation




void gridGenerator::get_facet_triangulation(vector<triangle> & triangle_list,
					    int x0, int y0, int z0,
					    int x1, int y1, int z1,
					    int d,  gridGenerator::facet_dir dir) 
{
  if(verbose >= 3) 
    fprintf(stdout, "get_facet_triang(x0=%d,y0=%d,z0=%d,x1=%d,y1=%d,z1=%d,d=%d,dir=%d\n",
	    x0,y0,z0,x1,y1,z1,d,dir);

  int3 e1(facet_basis[dir][0][0], facet_basis[dir][0][1], facet_basis[dir][0][2]);
  int3 e2(facet_basis[dir][1][0], facet_basis[dir][1][1], facet_basis[dir][1][2]);
  int3 e3(e1[0]+e2[0], e1[1]+e2[1], e1[2]+e2[2]);

  assert ( x1 == x0 + d*e3[0]);
  assert ( y1 == y0 + d*e3[1]);
  assert ( z1 == z0 + d*e3[2]);

  // corners must be vertices
  assert( isVertex(x0,y0,z0)); assert( isVertex(x1,y1,z1));
  assert( isVertex(x0+d*e1[0],y0+d*e1[1],z0+d*e1[2]));
  assert( isVertex(x0+d*e2[0],y0+d*e2[1],z0+d*e2[2]));

  if(d == 1) {
    get_diagonal_triangulation(triangle_list,x0,y0,z0,x1,y1,z1,d,dir);      
  }
  else { 
    int mx = (x1+x0)/2;
    int my = (y1+y0)/2;
    int mz = (z1+z0)/2;
    if(isVertex(mx,my,mz)) {
      // midpoint is vertex => facet sudivided => recurse
      get_facet_triangulation(triangle_list,x0,y0,z0, mx,my,mz,d/2,dir);
      get_facet_triangulation(triangle_list,
			      x0+d/2*e1[0], y0+d/2*e1[1], z0+d/2*e1[2],
			      mx+d/2*e1[0], my+d/2*e1[1], mz+d/2*e1[2], d/2,dir);
      get_facet_triangulation(triangle_list,
			      x0+d/2*e2[0], y0+d/2*e2[1], z0+d/2*e2[2],
			      mx+d/2*e2[0], my+d/2*e2[1], mz+d/2*e2[2], d/2,dir);
      get_facet_triangulation(triangle_list,
			      x0+d/2*e3[0], y0+d/2*e3[1], z0+d/2*e3[2],
			      mx+d/2*e3[0], my+d/2*e3[1], mz+d/2*e3[2], d/2,dir);
    }

    else {
      // vertices can only be on the boundary:
      // get them and perform a direct triangulation
      vector<int3> bdv; // boundary vertices of facet
      int3 E[4] = { e1, e2, int3(-e1[0],-e1[1],-e1[2]), int3(-e2[0],-e2[1],-e2[2])};

      // fprintf(stdout, "  bdv: ");
      int x = x0, y = y0, z = z0;
      for(int j = 0; j < 4; ++j) {
	for(int i = 0; i < d; ++i) {
	  x += E[j][0]; y+=E[j][1]; z+=E[j][2];
	  if(isVertex(x,y,z)) {
	    bdv.push_back(int3(x,y,z));
	    // fprintf(stdout, "(%d,%d,%d) ",x,y,z);
	  }
	}
      }
      // fprintf(stdout, "\n");
      assert(bdv.size() >= 4);
      if(bdv.size() == 4) {
	// no subdivided edges: special case
	get_diagonal_triangulation(triangle_list,x0,y0,z0,x1,y1,z1,d,dir);	
      }
      else {
	// more than 4 vertices: connect boundary vertices with midpoint
	int mn = addVertex(mx,my,mz);
	int3 p = bdv[0]; int3 q = bdv.back();
	triangle_list.push_back(triangle(mn, 
					 addVertex(p[0], p[1], p[2]),
					 addVertex(q[0], q[1], q[2])));
	for(unsigned i = 0; i < bdv.size() -1; ++i) {
	  p = bdv[i]; q = bdv[i+1];
	  triangle_list.push_back(triangle(mn, 
					   addVertex(p[0], p[1], p[2]),
					   addVertex(q[0], q[1], q[2])));
	}
      }
    }

  }
} // gridGenerator::get_facet_triangulation






/*! Add elements in tesselated_cells to primitive graph.
 */
bool gridGenerator::addPrimitives() 
{    
  if(verbose >= 1)
    fprintf(stdout, "Generating primitives ...\n");

  primitives = VCreateGraph(10, 9, VLongRepn, false);
  typedef octcell_tesselation::element_iterator elem_iter;
  for(unsigned c = 0; c < tesselated_cells.size(); ++c) {
    for(elem_iter e = tesselated_cells[c].begin_elements();
	e != tesselated_cells[c].end_elements(); ++e)  {
      int id;
      material_type mat = tesselated_cells[c].getMaterial(e);
      // skip elements with unwanted materials
      if(meshing_algorithm[mat] == meshing_type::none)
	continue;
      if((*e).size() == 4) 
	id = addTetra((*e)[0], (*e)[1], (*e)[2], (*e)[3]);
      else if( (*e).size() == 8) {
	int pl[8] = { (*e)[0], (*e)[1], (*e)[2], (*e)[3], (*e)[4], (*e)[5], (*e)[6], (*e)[7]};
	id = addCube(pl);
      }
      else {
	fprintf(stderr, "gridGenerator::addPrimitives(): primitive with %d nodes encountered, skipping!\n", (*e).size());
	continue;
      }
      // keep track of material
      while((int)material_of_primitive.size() <= id) 
	material_of_primitive.push_back(material_type(-1));
      material_of_primitive[id] = mat;
    }
    // could remove tesselation of cell c here
  }

  if(verbose >= 1)
    fprintf(stdout, "done generating primitives.\n");

  return true;
}

/*! Remove vertices not on surface, and update faces graph accordingly.
    
*/
bool gridGenerator::removeNonSurfaceVertices() 
{
  if(verbose >= 1)
    fprintf(stdout, "Remove non-surface vertices ...\n");

  std::vector<int> new_number(vertices->lastUsed+1,-1);
 
  VGraph new_vertices = VCreateGraph(10, 5, VFloatRepn, false);

  for(face* f = (face*) VGraphFirstNode(faces); f; f = (face *) VGraphNextNode(faces)) {
    int v[3];
    for(int i = 0; i < 3; ++i) {
      v[i] = f->vtx[i];
      if(new_number[v[i]] == -1) {
	vertex * old_v = (vertex *)VGraphGetNode(vertices,v[i]);
	vertex * new_v = new vertex(*old_v);
	new_number[v[i]] = VGraphAddAndGrow(new_vertices, (VNode)new_v, new_vertices->lastUsed+1);
      }
      f->vtx[i] = new_number[v[i]];
    }
    for(int i = 0; i < 3; ++i) {
      linkNodes(new_vertices, f->vtx[i], f->vtx[(i+1)%3]);
    }
  }

  VDestroyGraph(vertices);
  vertices = new_vertices;

  if(verbose >= 1)
    fprintf(stdout, "done removing non-surface vertices.\n");
 
  return true; 
}




/*
  For surface mesh generation: Create a new mesh triangle
*/
unsigned int gridGenerator::addTriangle(int i0, int i1, int i2)
{
  linkNodes(vertices, i0, i1);
  linkNodes(vertices, i1, i2);
  linkNodes(vertices, i2, i0);
	
	
  // allocate a new triangle
  assert(i0 != 0); assert(i1 != 0); assert(i2 != 0);
  face f(i0, i1, i2);
  return VGraphAddAndGrow(faces, (VNode)&f, ++faces->lastUsed);
}
 
/*
  Create a new vertex of octree or mesh
*/
unsigned int gridGenerator::addVertex(vertex &v)
{
  int i;
  unsigned int ref = pos(v);

  // get memory at voxel
  vertex *vtx = cache[ref];

  // if necessary, allocate and initialize
  if (vtx == 0) {
    vtx = new vertex [NPLACES];
    memset(vtx, 0, NPLACES*sizeof(vertex));
    cache[ref] = vtx;
  };

  // check if this vertex is already there
  for (i = 0; i < NPLACES; i++)  {
    if (vtx[i].hops == 0) break;		// unused
    if (vtx[i] == v) return vtx[i].hops;	// misuse hops here
  };
  assert(i < NPLACES);

  // else add element
  vtx[i] = v; vtx[i].hops = ++nv; vtx[i].color=0;
  int vpos= VGraphAddAndGrow(vertices, (VNode)&vtx[i], nv);
  if(verbose >= 3)
    fprintf(stdout,"adding vertex %d (%f,%f,%f)\n",nv,v.x,v.y,v.z);         

  return vpos;
}


/*
  Link four points to form a tetrahedron,
  add tetra to primitive list
*/
unsigned int gridGenerator::addTetra(int i0, int i1, int i2, int i3)

{
  // link vertices to form a tetrahedron
  // this could be skipped in the future 
  // (use implicit_connectivity)

  if(! implicit_links) {
    linkNodes(vertices, i0, i1);
    linkNodes(vertices, i0, i2);
    linkNodes(vertices, i0, i3);
    linkNodes(vertices, i1, i2);
    linkNodes(vertices, i1, i3);
    linkNodes(vertices, i2, i3);
  }
  // check orientation of tet i0,i1,i2,i3
  double j[3][3];
  vertex* v[4];
  v[0] = (vertex *)VGraphGetNode(vertices, i0);
  v[1] = (vertex *)VGraphGetNode(vertices, i1);
  v[2] = (vertex *)VGraphGetNode(vertices, i2);
  v[3] = (vertex *)VGraphGetNode(vertices, i3);
  for(int k = 0; k <= 2; ++k) {
    j[0][k] = v[k]->x - v[3]->x; 
    j[1][k] = v[k]->y - v[3]->y; 
    j[2][k] = v[k]->z - v[3]->z; 
  }
  double det = det3(j);
  int ii2 = (det > 0 ? i2 : i3);
  int ii3 = (det > 0 ? i3 : i2);

  // allocate a new tetrahedron
  primitive p;
  assert(i0 != 0); assert(i1 != 0); assert(i2 != 0); assert(i3 != 0);
  p.vcnt = 4; 
  p.id[0] =  i0; p.id[1] =  i1; 
  p.id[2] = ii2; p.id[3] = ii3;
  return VGraphAddAndGrow(primitives, (VNode)&p, ++primitives->lastUsed);		
}


/*
  link eight points to form a cube
  add cube to primitive list
*/
unsigned int gridGenerator::addCube(int *plist)

{
  // link vertices to form a cube
  if(! implicit_links) {

    linkNodes(vertices, plist[0], plist[1]);
    linkNodes(vertices, plist[0], plist[2]);
    linkNodes(vertices, plist[0], plist[4]);
    linkNodes(vertices, plist[3], plist[1]);
    linkNodes(vertices, plist[3], plist[2]);
    linkNodes(vertices, plist[3], plist[7]);

    linkNodes(vertices, plist[5], plist[4]);
    linkNodes(vertices, plist[5], plist[7]);
    linkNodes(vertices, plist[5], plist[1]);
    linkNodes(vertices, plist[6], plist[4]);
    linkNodes(vertices, plist[6], plist[7]);
    linkNodes(vertices, plist[6], plist[2]);
  }
  // allocate a new cube
  primitive p;
  p.vcnt = 8;
  // this ordering has to be changed for ipe visualisation	
  // internal connectivity set to Hughes page 136 convention, JF 20.8.2001
  p.id[0] = plist[0];
  p.id[1] = plist[1];
  p.id[2] = plist[3];
  p.id[3] = plist[2];
  p.id[4] = plist[4];
  p.id[5] = plist[5];
  p.id[6] = plist[7];
  p.id[7] = plist[6];


  if (verbose > 1) {
    fprintf(stdout,"e%d  %i %i %i %i %i %i %i %i\n",primitives->lastUsed,
	    p.id[0],p.id[1],p.id[2],p.id[3],p.id[4],p.id[5],p.id[6],p.id[7]);
    fflush(stdout);
  }

  return VGraphAddAndGrow(primitives, (VNode)&p, ++primitives->lastUsed);
}


/*
  Link eight points forming a cube.
  Only used if vertex connectivity is chosen.
*/
void gridGenerator::linkCube(int id)
{
  primitive *p = (primitive *)VGraphGetNode(primitives, id);
  if (p == 0) return;

  linkNodes(vertices, p->id[0], p->id[1]);
  linkNodes(vertices, p->id[0], p->id[2]);
  linkNodes(vertices, p->id[0], p->id[4]);
  linkNodes(vertices, p->id[3], p->id[1]);
  linkNodes(vertices, p->id[3], p->id[2]);
  linkNodes(vertices, p->id[3], p->id[7]);

  linkNodes(vertices, p->id[5], p->id[4]);
  linkNodes(vertices, p->id[5], p->id[7]);
  linkNodes(vertices, p->id[5], p->id[1]);
  linkNodes(vertices, p->id[6], p->id[4]);
  linkNodes(vertices, p->id[6], p->id[7]);
  linkNodes(vertices, p->id[6], p->id[2]);
}





/*****************************************************************************
 
            Member functions for marching tetrahedra algorithm

*****************************************************************************/

/*
  Trilinear interpolation of binary cell classification on vertices
*/
bool gridGenerator::interpolateOnVertices(VImage f_materials) const 
{
  if(verbose >= 1)
    fprintf(stdout, "interpolateOnVertices  ...\n");

  for(int d = min_maxdim; d <= max_maxdim; d *= 2) {
    // for each cell of size d
    for (int z = 0; z < nz; z += d) {
      for (int y = 0; y < ny; y += d) {
	for (int x = 0; x < nx; x += d) {
	  if(getCellDim(x,y,z) == d) {
	    octcell c(*this, x,y,z,d);
	    material_type m = getCellMaterial(c);
	    float f = (m == mat_body ? -1 : +1);
	    for(int vx = 0; vx <=d; ++vx)
	      for(int vy = 0; vy <=d; ++vy)
		for(int vz = 0; vz <=d; ++vz) {
		  // inverse weight this cell contributes to value of vertex vx,vy,vz 
		  // vertices at corners have weight 1/8, at edges 1/4, at faces 1/2.
		  int w_inner_inv = (0<vx && vx<d? 1:2) * (0<vy && vy<d? 1:2) *(0<vz && vz<d? 1:2);
		  // don't consider inner points
		  if( w_inner_inv > 1) {
		    int gx = x+vx, gy = y+vy, gz=z+vz;
		    // vertices on boundary of grid are incident to fewer cells, so have bigger weight.
		    int w_outer = (0<gx && gx<nx? 1:2)*(0<gy && gy<ny? 1:2)*(0<gz && gz<nz? 1:2);
		    if(VPixel(f_materials,gz,gy,gx,VFloat) == 77)
		      VPixel(f_materials,gz,gy,gx,VFloat) = 0;
		    VPixel(f_materials,gz,gy,gx,VFloat) += f* w_outer*1.0/w_inner_inv;
		  }
		}
	    // what about potential midpoint vertex? -> handled directly in marchingTetras()
	  }
	}
      }
    }
  } 
  if(verbose >= 1)
    fprintf(stdout, "done interpolateOnVertices.\n");

  if(verbose >=3) {
    fprintf(stdout, "Vertex values:\n");
    for(int z = 0; z < nz+1; ++z) {
      fprintf(stdout, "z=%d\n",z);
      for(int y = ny; y >=0; y--) {
	for(int x = 0; x < nx+1; ++x)
	  fprintf(stdout, "%1.2f ", VPixel(f_materials,z,y,x,VFloat));
	fprintf(stdout, "\n");  
      } 
    }  
  }

  return true;
}



/*
  The main marching tetrahedra algorithm
*/
bool gridGenerator::marchingTetras() 
{
  if(verbose >= 1)
    fprintf(stdout, "marching tets ...\n");

  VImage f_materials = VCreateImage(nz+1, ny+1, nx+1,VFloatRepn);
  for(int x = 0; x < nx+1; ++x)
    for(int y = 0; y < ny+1; ++y)
      for(int z = 0; z < nz+1; ++z)
	VPixel(f_materials,z,y,x,VFloat) = 77; // init to dummy value
  
  interpolateOnVertices(f_materials);

  // tesselate untesselated cells at boundary of material to be smoothed by marching tets,
  // e.g. background cells 
  for(int d = min_maxdim; d <= max_maxdim; d *= 2) {
    // for each cell of size d
    for (int z = 0; z < nz; z += d) {
      for (int y = 0; y < ny; y += d) {
	for (int x = 0; x < nx; x += d) {
	  if(getCellDim(x,y,z) == d) {
	    // octcell c(*this, x,y,z,d);
	    // check if at least one  node of c has f(v) != +- 1
	    bool cell_cut = false;
	    for(int vx = 0; vx <= 1; vx++) 
	      for(int vy = 0; vy <= 1; vy++) 
		for(int vz = 0; vz <= 1; vz++) {
		  float f = VPixel(f_materials,z+vz*d, y+vy*d,x+vx*d,VFloat);
		  if(-1 < f && f < 1)
		    cell_cut = true;
		}
	    // cut octcell which has not been tesselated => tesselate
	    // better would be: check if c has been tesselated.
	    if(cell_cut && ! mustTesselateCell(x,y,z)) { 
	      tesselateCell(x,y,z);
	    } 
	  }
	}
      }
    }
  }
  
  // loop over all tesselated cells, cut its tets if necessary
  namespace mt = marching_tets;
  mt::cut_vertex_map cut_edges;
  for(int ic = 0; ic < (int)tesselated_cells.size(); ++ic) {
    typedef octcell_tesselation::element_iterator elem_iter;
    octcell_tesselation  new_tesselation; // tesselation with cut cells
    // possible midpoint has not yet been got a value
    if( tesselated_cells[ic].hasMidpoint()) {
      int x,y,z;
      getVertexIntegerCoords(tesselated_cells[ic].getMidpoint(), x,y,z);
      VPixel(f_materials,z,y,x,VFloat) = ( getCellMaterial(tesselated_cells[ic].TheCell()) == mat_body ? -1 : +1);
    }
    for(elem_iter e = tesselated_cells[ic].begin_elements();
	e != tesselated_cells[ic].end_elements(); ++e)  {
      if( (*e).size() > 4) {
	fprintf(stderr, "WARNING: marchingTetras(): cannot  handle elements with %d > 4 nodes, skipping\n", (*e).size());
	continue;
      }
      // get cutting pattern of tet e
      int pat[4];
      int nplus =0, nminus = 0;
      for(int v = 0; v < 4; ++v) {
	int vx, vy, vz;
	getVertexIntegerCoords((*e)[v], vx,vy,vz);
	float f = VPixel(f_materials,vz,vy,vx, VFloat);
	pat[v] = ( f < -1.0/16 ? -1 : ( f > 1.0/16 ? +1 : 0 ));
	if(pat[v] ==  0) {
	  vertex * vv = (vertex*) VGraphGetNode(vertices, (*e)[v]); assert (vv != 0);
	  vv->color = 1; // for vsmooth
	}
	if(pat[v] == -1) nminus++;
	if(pat[v] == +1) nplus++;
      }
      /*
	if(nminus == 0 && nplus == 0)
	if(verbose > 2) 
	fprintf(stderr, "WARNING: detected tet (%d,%d,%d,%d) with all 4 vertices having value 0!\n",
	(*e)[0], (*e)[1], (*e)[2], (*e)[3]);
      */
      // do we have to cut this tet?
      //if(nminus > 0 && nplus > 0) {
      // subdivide tet e
      makeCuts(*e, cut_edges, f_materials, pat); 
      mt::subdivision sub(pat, mt::tet((*e)[0], (*e)[1], (*e)[2], (*e)[3]), cut_edges);
      if(!surf_arg) { // volume mesh generation
	for(mt::subdivision::tet_iterator ti = sub.FirstTet(); ti != sub.EndTet(); ++ti) {
	  material_type mat = (sub.Sign(ti) == 1 ?  mat_background  :
			       ( sub.Sign(ti) == -1 ? mat_body :
				 getCellMaterial(tesselated_cells[ic].TheCell()) ));
	  if(material_list[mat].mesh_type != meshing_type::none)
	    new_tesselation.addTetra( (*ti)[0], (*ti)[1], (*ti)[2], (*ti)[3] , mat);
	}
      }
      else { // surface mesh generation
	for(mt::subdivision::triangle_iterator ti = sub.FirstSurfaceTriangle(); 
	    ti != sub.EndSurfaceTriangle(); ++ti) {
	  if(sub.Sign(ti) >= 0) {
	    // orient ti such that normal formed by (ti[1]-ti[0]) x (ti[2]-ti[0])
	    // points into positive material.
	    double j[3][3];
	    vertex* v[4];
	    v[0] = (vertex *)VGraphGetNode(vertices, (*ti)[0]);
	    v[1] = (vertex *)VGraphGetNode(vertices, (*ti)[1]);
	    v[2] = (vertex *)VGraphGetNode(vertices, (*ti)[2]);
	    // sub.Sign(ti) >= 0 => there must be at least one positive vertex
	    // ?? In case of 2 positive vertices vp(i) (and 2 triangles), do 
	    // both vp(0) and vp(1) always lie on the same side of each triangle??
	    v[3] = (vertex *)VGraphGetNode(vertices, sub.vp(0));
	    for(int k = 0; k <= 2; ++k) {
	      j[0][k] = v[k]->x - v[3]->x; 
	      j[1][k] = v[k]->y - v[3]->y; 
	      j[2][k] = v[k]->z - v[3]->z; 
	    }
	    double det = det3(j);
	    if(det > 0)
	      addTriangle((*ti)[0], (*ti)[1], (*ti)[2]);
	    else
	      addTriangle((*ti)[0], (*ti)[2], (*ti)[1]);
	  }
	}
      }
      // }
    }
    if(! new_tesselation.empty() || surf_arg)
      tesselated_cells[ic] = new_tesselation; 
  }
  VDestroyImage(f_materials);

  if(verbose >= 1)
    fprintf(stdout, "done marching tets.\n");

  return true;
}




/*
  classify edges for the marching tetrahedra algorithm
*/
bool gridGenerator::makeCuts(std::vector<int> const& t,
			     marching_tets::cut_vertex_map & cuts,
			     VImage f_material,
			     int * pat) 
{
  for(int i = 0; i < 4; ++i)
    for(int j = i+1; j < 4; ++j) {
      if( pat[i]*pat[j] == -1) {
	marching_tets::edge e(t[i],t[j]);
	if(cuts.find(e) == cuts.end()) {
	  int xi,yi,zi, xj, yj,zj;
	  getVertexIntegerCoords(t[i], xi,yi,zi);
	  getVertexIntegerCoords(t[j], xj,yj,zj);
	  float fi = VPixel(f_material,zi,yi,xi,VFloat);
	  float fj = VPixel(f_material,zj,yj,xj,VFloat);
          float lambda = -fi/(fj-fi);
	  float x = (1-lambda)*xi + lambda*xj;
	  float y = (1-lambda)*yi + lambda*yj;
	  float z = (1-lambda)*zi + lambda*zj;
	  // NOTE: we cannot use addVertex here, because that works only for integer coordinates!
	  vertex * v = new vertex(x,y,z);
	  v->color = 1; // for vsmooth
	  v->hops  = ++nv;
	  cuts[e] = VGraphAddAndGrow(vertices, (VNode)v, nv);
	}
      }
    }
  return true;
} 


void gridGenerator::getVertexIntegerCoords(int v, int& vx, int& vy, int& vz)  const
{
  vertex* vv = (vertex*)VGraphGetNode(vertices, v);  assert(vv != 0);
  // if the vertex coords don't have integer values, we have a big problem.
  vx = (int)vv->x;    assert( (float)vx == vv->x);
  vy = (int)vv->y;    assert( (float)vy == vv->y);
  vz = (int)vv->z;    assert( (float)vz == vv->z);
}



/****************************************************************************
 
               member functions for surface smoothing 

****************************************************************************/

/*
  Smooth marked surface vertices obtained from marching tetrahedra algorithm
*/
bool gridGenerator::smoothMarchingInterfaces() 
{
  if(verbose >= 1)
    fprintf(stdout, "smoothing boundary ...\n");
  int nv = vertices->lastUsed;
  std::vector<vec3> new_coords(nv);
  for(int i = 0; i < smooth_iters*surface_smooth_weights.size(); ++i) {
    int numv = 0;
    for (vertex * v = (vertex *) VGraphFirstNode (vertices); v;
	 v = (vertex *) VGraphNextNode (vertices), ++numv){
      if(v->color == 1 ){
	new_coords[numv] = v->smooth_surface_coord(vertices, surface_smooth_weights[i%surface_smooth_weights.size()]); 
      }
      else
	new_coords[numv] = v->smooth_volume_coord (vertices, volume_smooth_weights[i%volume_smooth_weights.size()]);
    }
    numv = 0;
    for (vertex * v = (vertex *) VGraphFirstNode (vertices); v;
	 v = (vertex *) VGraphNextNode (vertices), ++numv){
      //if(v->color == 1 ){
	v->setPoint(new_coords[numv]);
      //}
    }

  }
  if(verbose >= 1)
    fprintf(stdout, "done smoothing boundary.\n");

  return true;
} 


/* 
   Smooth interfaces for uniform brick meshes
*/
bool gridGenerator::shiftBoundaryVertices()
{

  fprintf(stdout, "Shifting boundary vertices, shift=%f\n", shift);
  int mat[8], cnt[8];
	
  std::vector<vec3> displacements(vertices->lastUsed+1, vec3(0,0,0));

  // for each vertex in this mesh
  for (int id = 1; id <= vertices->lastUsed; id++)  {
    vertex *v = (vertex *)VGraphGetNode(vertices, id);
    if (v == 0) continue;
    int ox,oy,oz;
    getVertexIntegerCoords(id,ox,oy,oz);
    // don't shift vertices on boundary
    if(ox == 0 || ox == nx || oy == 0 || oy == ny || oz == 0 || oz == nz)
      continue;
    int dim = getCellDim(ox, oy, oz);
    // find material labels of neighbouring cells
    int j, n = 0;
    for (int i = 0; i < 8; i++)  {
      // take care of boundaries
      // note new convention: vertices are *between* voxels
      int ex = ox-off_x[i]*dim; if (ex < 0 || ex >= nx) continue;
      int ey = oy-off_y[i]*dim; if (ey < 0 || ey >= ny) continue;
      int ez = oz-off_z[i]*dim; if (ez < 0 || ez >= nz) continue;
      int v = (int)VPixel(src, ez, ey, ex, VUByte);
      for (j = 0; j < n; j++)
	if (v == mat[j]) break;
      if (j == n) {	// new material
	mat[n] = v; cnt[n] = 1; n++;
      } else cnt[j]++;
    };
	  
    vec3 d = 0;
    // determine nodal displacement for 2-material boundaries
    if (n == 2) {
      int minor = (cnt[0] < cnt[1])? mat[0]: mat[1];
	    
      n = 0; vec3 c = 0;
      for (int i = 0; i < 8; i++)  {
	int ex = ox-off_x[i]*dim; if (ex < 0 || ex >= nx) continue;
	int ey = oy-off_y[i]*dim; if (ey < 0 || ey >= ny) continue;
	int ez = oz-off_z[i]*dim; if (ez < 0 || ez >= nz) continue;
	if ((int)VPixel(src, ez, ey, ex, VUByte) != minor) continue;
	c.x += (0.5-off_x[i])*dim;
	c.y += (0.5-off_y[i])*dim;
	c.z += (0.5-off_z[i])*dim;
	n++;
      };
      d = c * (shift/n);
    };
    displacements[id] = d;
  }

  // add displacements to vertex coordinates
  for (int id = 1; id <= vertices->lastUsed; id++)  {
    vertex *v = (vertex *)VGraphGetNode(vertices, id);
    if (v == 0) continue;
    v->x += displacements[id].x;
    v->y += displacements[id].y;
    v->z += displacements[id].z;
  }

  return true;
}






/*****************************************************************************

                  member functions for mesh output

*****************************************************************************/

static void printAttributes( VAttrList list)
{
  VAttrRec *a, *a_next;
  if (list) {
    fprintf(stdout," name = <%s>\n",list->name);
    for (a = list->next; a; a = a_next) {
      a_next = a->next;
      fprintf(stdout," name = <%s>\n",a->name);
    }
  }
}

static void trimGraph( VGraph graph, const char* name, int verbose) 
{
  if(verbose > 1) {
    fprintf(stdout,"%s  ->lastUsed  =%d\n",name,graph->lastUsed);
    fprintf(stdout,"%s  ->size      =%d\n",name,graph->size    );
    fprintf(stdout,"%s  ->nnodes    =%d\n",name,graph->nnodes  );
    fprintf(stdout,"countNodes(%s)  =%d\n",name,countNodes(graph));
  }
  assert(graph->lastUsed   == (int)countNodes(graph));
  graph->size     = graph->lastUsed;
  graph->nnodes   = graph->lastUsed;
}

// copy attributes of the src image to the output graph (vertices)

static void copyImageAttributes(VImage img, VGraph graph)
{
  VPointer attr;
  // copy/add compulsory attributes (section 4.2)
  // copy orientation
  if (VGetAttr(VImageAttrList(img), "orientation", 0, VStringRepn, &attr) == VAttrFound) {
    VSetAttr(VGraphAttrList(graph), "orientation", NULL, VStringRepn, attr);
  }

  // copy patient name
  if (VGetAttr(VImageAttrList(img), "patient", 0, VStringRepn, &attr) == VAttrFound) {
    VSetAttr(VGraphAttrList(graph), "patient", NULL, VStringRepn, attr);
  }
  else {
    VSetAttr(VGraphAttrList(graph), "patient", NULL, VStringRepn, "VGrid");
  }
	    
  // copy examination date
  if (VGetAttr(VImageAttrList(img), "date", 0, VStringRepn, &attr) == VAttrFound) {
    VSetAttr(VGraphAttrList(graph), "date", NULL, VStringRepn, attr);
  }
  else {
    VSetAttr(VGraphAttrList(graph), "date", NULL, VStringRepn, "25.7.2002");
  }
	    
  // copy convention
  if (VGetAttr(VImageAttrList(img), "convention", 0, VStringRepn, &attr) == VAttrFound) {
    VSetAttr(VGraphAttrList(graph), "convention", NULL, VStringRepn, attr);
  }
  else {
    VSetAttr(VGraphAttrList(graph), "convention", NULL, VStringRepn, "natural");
  }

  // copy position
  if (VGetAttr(VImageAttrList(img), "position", 0, VStringRepn, &attr) == VAttrFound) {
    VSetAttr(VGraphAttrList(graph), "position", NULL, VStringRepn, attr);
  }

  // copy name (i.e. CT etc.)
  if (VGetAttr(VImageAttrList(img), "name", 0, VStringRepn, &attr) == VAttrFound) {
    VSetAttr(VGraphAttrList(graph), "name", NULL, VStringRepn, attr);
  }

  // copy classes
  if (VGetAttr(VImageAttrList(img), "classes", 0, VStringRepn, &attr) == VAttrFound) {
    VSetAttr(VGraphAttrList(graph), "classes", NULL, VStringRepn, attr);
  }

  // copy voxel dimensions
  if (VGetAttr(VImageAttrList(img), "voxel", 0, VStringRepn, &attr) == VAttrFound) {
    VSetAttr(VGraphAttrList(graph), "voxel", NULL, VStringRepn, attr);
  }
}

bool gridGenerator::write(FILE *outf)
{

  if (vertices == 0 || (primitives == 0 && !surf_arg) || (faces == 0 && surf_arg) )  {
    fprintf(stderr, "gridGenerator::write: incomplete grid.\n");
    return false;
  }
        
  VImage                mp = 0;                     // 1D image storing the material properties
  VImage                iv = 0;                     // 1D image source node
  VImage                pe = 0;                     // 1D image elem partition info
  VImage                pn = 0;                     // 1D image node partition info
  VImage                tI = 0;                     // 1D image nodal constraints
  VImage                fI = 0;                     // 1D image initial nodal forces


  int nv = vertices->lastUsed;		// # of vertices in volume mesh
  int ne = (surf_arg == false ? primitives->lastUsed : faces->lastUsed);
  if (verbose) {
    fprintf(stdout, "final grid contains %d vertices and %d primitives.\n",nv,ne);
  }		


  // allocate attribute images  for volume mesh 	
		
  if(surf_arg==false){ // volume mesh

    mp = VCreateImage(1, 1, ne, VUByteRepn);
    VImageNFrames( mp )     = 1;
    VImageNComponents( mp ) = 1;
    if ( mp == 0 ) {
      fprintf(stderr, "generateCells: cannot create attribute images\n");
      return false;
    }
    if(doOutputFields()) {
      pe = VCreateImage(1, 1, ne, VUByteRepn);
      pn = VCreateImage(1, 1, nv, VUByteRepn);
      iv = VCreateImage(1, 1, nv, VBitRepn  );
      tI = VCreateImage(1, 1, nv, VUByteRepn);
      fI = VCreateImage(3, 1, nv, VFloatRepn);
      if (iv == 0 || pn == 0 || pe == 0 || tI == 0 || fI ==0 )  {
	fprintf(stderr, "generateCells: cannot create attribute images\n");
	return false;
      }

      VImageNFrames( pe )     = 1;
      VImageNComponents( pe ) = 1;
      VImageNFrames( pn )     = 1;
      VImageNComponents( pn ) = 1;
      VImageNFrames( iv )     = 1;
      VImageNComponents( iv ) = 1;
      VImageNFrames( tI )     = 1;
      VImageNComponents( tI ) = 1;
      VImageNFrames( fI )     = 1;
      VImageNComponents( fI ) = 3;
    }
  } // volume mesh

  // modify vertex positions for real-world dimensions (section 3.3)
  VStringConst attr;
  vec3 voxel = 1.0;
	
  if (VGetAttr(VImageAttrList(src), "voxel", 0, VStringRepn, &attr) == VAttrFound){
    sscanf(attr, "%lf%lf%lf", &voxel.x, &voxel.y, &voxel.z);
    for (vertex *v 	= (vertex *)VGraphFirstNode(vertices); v;
	 v= (vertex *)VGraphNextNode(vertices))  {
      vec3 p = v->getPoint() * voxel;
      v->setPoint(p);
    }
  }
  else {
    fprintf(stderr, "WARNING: no voxel attribute found in image, using unit voxel size\n");
  }

  if(verbose)
    fprintf(stdout, "saving grid ...\n");
	
	
  if(surf_arg==false){ // volume mesh
    copyImageAttributes(src, vertices);

    //generate material property images
    if (verbose) fprintf(stdout, "WriteMaterialPid ...\n");
    WriteMaterialPid(mp, ne);
    // add material attribute to primitives graph: list of mat-name mat-id
    if(use_simbio_mat_ids) {
      string d;
      for(int m = 0; m < (int)material_list.size(); ++m) {
	d +=  material_list[m].name;
	char buf[21];  
	snprintf(buf, 20, " %d   ", MatDB.id(material_list[m].name));
	d += buf;
      }
      VSetAttr(VGraphAttrList(primitives), "materials", NULL, VStringRepn, d.c_str());
    }

    if (verbose) fprintf(stdout, "Set attributes   ...\n");
    VAppendAttr(VImageAttrList(mp), "component_repn",   0, VStringRepn, "scalar" );
    VAppendAttr(VImageAttrList(mp), "component_interp", 0, VStringRepn, "element label");
    VAppendAttr(VGraphAttrList(primitives), "matprops",   NULL, VImageRepn, mp);

    // see section 4.2
    VSetAttr(VGraphAttrList(vertices),   "component_interp", NULL, VStringRepn, "vertex");
    VSetAttr(VGraphAttrList(primitives), "component_interp", NULL, VStringRepn, "primitive");
	    
    // see section 4.7
    VSetAttr(VGraphAttrList(primitives), "primitive_interp", NULL, VStringRepn, "volume");
	    
    if(implicit_links)
      VSetAttr(VGraphAttrList(primitives), "implicit_links", NULL, VStringRepn, "true");
    else
      VSetAttr(VGraphAttrList(primitives), "implicit_links", NULL, VStringRepn, "false");
	    

    if(doOutputFields()) {
      // generate the other fields

      if (verbose) fprintf(stdout, "WriteTagImage    ...\n");
      WriteTagImage(tI,nv);
      if (verbose) fprintf(stdout, "WriteForceImage  ...\n");
      WriteForceImage(fI,nv);
      if (verbose) fprintf(stdout, "WriteSourceNodes ...\n");
      WriteSourceNodes(iv,nv);
      if (verbose) fprintf(stdout, "WriteInitialPart ...\n");
      WriteInitialPart(pe, pn, ne, nv); // may adapt npart
      // additional attribute # of partitions
      VSetAttr(VGraphAttrList(vertices), "partition(s)", NULL, VUByteRepn, npart);
	    
	    
      //image attributes
      VAppendAttr(VImageAttrList(fI), "component_repn",   0, VStringRepn, "vector3");
      VAppendAttr(VImageAttrList(fI), "component_interp", 0, VStringRepn, "force N");
	    
      VAppendAttr(VImageAttrList(tI), "component_repn",   0, VStringRepn, "scalar" );
      VAppendAttr(VImageAttrList(tI), "component_interp", 0, VStringRepn, "tag"    );
	    
      VAppendAttr(VImageAttrList(iv), "component_repn",   0, VStringRepn, "scalar" );
      VAppendAttr(VImageAttrList(iv), "component_interp", 0, VStringRepn, "src node tag");
	    
      VAppendAttr(VImageAttrList(pn), "component_repn",   0, VStringRepn, "scalar" );
      VAppendAttr(VImageAttrList(pn), "component_interp", 0, VStringRepn, "nodal partition");
	    
	    
      VAppendAttr(VImageAttrList(pe), "component_repn",   0, VStringRepn, "scalar" );
      VAppendAttr(VImageAttrList(pe), "component_interp", 0, VStringRepn, "element partition");
	    
      //graph attributes
      VAppendAttr(VGraphAttrList(vertices),   "image",      NULL, VImageRepn, tI);
      VAppendAttr(VGraphAttrList(vertices),   "forceimage", NULL, VImageRepn, fI);
      VAppendAttr(VGraphAttrList(vertices),   "srcnodes",   NULL, VImageRepn, iv);
      VAppendAttr(VGraphAttrList(vertices),   "partnode",   NULL, VImageRepn, pn);
      VAppendAttr(VGraphAttrList(primitives), "partelem",   NULL, VImageRepn, pe);
    }// doOutputFields()


    if (verbose > 1) {
      fprintf(stdout,"Attributes of vertices graph:\n");
      printAttributes(vertices->attributes);
      fprintf(stdout,"Attributes of primitives graph:\n");
      printAttributes(primitives->attributes);
    }
   
    VGraph g[2];
    if(out_arg == vm || out_arg==hf || out_arg == gmv || out_arg == gpp || out_arg == ascii || out_arg == dx){
      trimGraph(vertices,   "vertices",   verbose);
      trimGraph(primitives, "primitives", verbose);
      //vm & headfem convention	
      g[0] = vertices;    
      g[1] = primitives;
      if(out_arg == vm || out_arg==hf) 
	VWriteGraphs(outf, NULL, 2 , g); 
      else if(out_arg == gmv) {
	mesh m(vertices, primitives, verbose);
	m.vtogmv(outf);
      }
      else if(out_arg == ascii) {
	mesh m(vertices, primitives, verbose);
	m.vtoa(outf);
      }
      else if(out_arg == gpp) {
	mesh m(vertices, primitives, verbose);
	m.vtogpp(outf);
      }
      else if(out_arg == dx) {
	mesh m(vertices, primitives, verbose);
	m.vtodx(outf);
      }

      else {
	fprintf(stderr, "WARNING: output format %d unknown, using Vista output!\n", out_arg);
	VWriteGraphs(outf, NULL, 2 , g); 
      }
    }
    else{
      // attention this is ipe convention !!!
      VWriteGraphs(outf, NULL, 1, &vertices); 
    }
  } // volume mesh

  else { // surface mesh
    VGraph g[2];
    if(out_arg == vm || out_arg==hf){
      copyImageAttributes(src, vertices);
      VSetAttr(VGraphAttrList(vertices),   "component_interp", NULL, VStringRepn, "vertex");
      VSetAttr(VGraphAttrList(faces),      "component_interp", NULL, VStringRepn, "primitive");
      // see section 4.7
      VSetAttr(VGraphAttrList(faces), "primitive_interp", NULL, VStringRepn, "surface");
	 
      //vm & headfem convention	
      fprintf(stderr, "WARNING: gridGenerator::write : dirty --- cut off the unused tail of VISTA graph!\n");
      trimGraph(vertices,"vertices",verbose);
      trimGraph(faces,   "faces",   verbose);

      g[0] = vertices;
      g[1] = faces;
      VWriteGraphs(outf, NULL, 2, g);
    }
    else if(out_arg == complex2d) {
      fprintf(outf, "%d  %d\n", vertices->lastUsed, faces->lastUsed);
      for (vertex * v = (vertex *) VGraphFirstNode (vertices); v; v = (vertex *) VGraphNextNode (vertices)){
	fprintf(outf, "%f %f %f\n", v->x, v->y, v->z);
      }
      for(face* f = (face*) VGraphFirstNode(faces); f; f = (face *) VGraphNextNode(faces)) {
	fprintf(outf, "3  %d %d %d\n", f->vtx[0]-1, f->vtx[1]-1, f->vtx[2]-1);
      }
    }
    else if(out_arg == gmv) {
      fprintf(outf, "gmvinput ascii\n");
      fprintf(outf, "nodev %d\n",vertices->lastUsed); 
      for (vertex * v = (vertex *) VGraphFirstNode (vertices); v; v = (vertex *) VGraphNextNode (vertices)){
	fprintf(outf, "%f %f %f\n", v->x, v->y, v->z);
      }
      fprintf(outf, "cells 0\n");
      fprintf(outf, "surface %d\n", faces->lastUsed);
      for(face* f = (face*) VGraphFirstNode(faces); f; f = (face *) VGraphNextNode(faces)) {
	fprintf(outf, "3  %d %d %d\n", f->vtx[0], f->vtx[1], f->vtx[2]);
      }
      fprintf(outf,"endgmv\n");
     }
    else if(out_arg == dx) {
      fprintf(outf, "object 1 class array type float rank 1 shape 3  items %d  data follows\n", vertices->lastUsed);
      for (vertex * v = (vertex *) VGraphFirstNode (vertices); v; v = (vertex *) VGraphNextNode (vertices)){
	fprintf(outf, "%f %f %f\n", v->x, v->y, v->z);
      }
      fprintf(outf, "object 2 class array type int rank 1 shape 3 items %d data follows\n", faces->lastUsed);
      for(face* f = (face*) VGraphFirstNode(faces); f; f = (face *) VGraphNextNode(faces)) {
	fprintf(outf, "%d %d %d\n", f->vtx[0]-1, f->vtx[1]-1, f->vtx[2]-1);
      }
      fprintf(outf, "attribute \"element type\" string \"triangles\"\n");
      fprintf(outf, "attribute \"ref\" string \"positions\"\n");
      fprintf(outf, "object \"irregular positions irregular connections\" class field\n");
      fprintf(outf, "component \"positions\" value 1\n");
      fprintf(outf, "component \"connections\" value 2\n");
      fprintf(outf, "end\n");
    }
    else{
      // attention this is ipe convention !!!
      VWriteGraphs(outf, NULL, 1, &vertices); 
    }
  } // surface mesh		

  if(verbose) {
    fprintf(stdout, "... done saving grid\n");
    fflush(stdout);
  }
  return true;	
}




/*! set the mp image to the material ID of primitives

\todo Currently, cells are represented twice: in primitives graph and
the  tesselated_cells list. We probably can do without the primitives graph.
*/
void gridGenerator::WriteMaterialPid(VImage mp, int ne)						
{
  vector<unsigned> hist  (num_of_materials(),0);
  vector<int>      matids(num_of_materials(),-1);

  typedef octcell_tesselation::element_iterator elem_iter;
  //  for(unsigned c = 0; c < tesselated_cells.size(); ++c) {
  //  for(elem_iter e = tesselated_cells[c].begin_elements();
  //	e != tesselated_cells[c].end_elements(); ++e)  {
  for (int i = 1; i <= primitives->lastUsed; i++) {
    primitive *p = (primitive*) VGraphGetNode(primitives, i);
    assert(p != NULL);
    
    material_type m =  material_of_primitive[i]; // getCellMaterial(tesselated_cells[c].TheCell());
    // histogramme of material types      
    hist[m]++;
    
    // make a meaningful material ID out of m. (m is in [0, num_of_materials()[ )
    // Current solution:
    // labels:     0 ...  99  for tetrahedra  (p->vcnt = 4)
    //           100 ... 199  for hexahedra   (p->vcnt = 8)
    //  -> use MTypeof(int m) to get material no. out of m
    //  -> use ETypeof(int m) to get element type out of m
    // (both in headfem-utilities package.)
    
    if(use_simbio_mat_ids)
      m=matids[m] = MatDB.id(material_list[m].name);
    else
      matids[m] = m;

    if (p->vcnt == 4) m += 0;        
    if (p->vcnt == 8) m += 100; 
    
    // set mp image
    VSetPixel(mp, 0, 0, i-1, m);
  }
  
  if(verbose > 0) {
    int sum=0;
    fprintf(stdout, "WriteMaterialPid: histogram:\n");
    for(int m = 0; m < (int)hist.size(); ++m)
      if (hist[m] > 0) {
	fprintf(stdout, " Material No. %d (\"%s\", id %d): %d elements\n",
		m, material_name(m).c_str(), matids[m], hist[m]);
        sum += hist[m];
      }
    fprintf(stdout, " Total number: %d elements\n",sum);
    fflush(stdout);
  }
 
}



void gridGenerator::WriteTagImage(VImage tI, int nv)
{
  for (int i=0; i<nv; i++){
    VUByte u=0;
    VSetPixel(tI, 0, 0, i, u);
  }

  if(constype==box){
    FILE *datei;
    int xmin, ymin, zmin;
    int xdim, ydim, zdim;
    int cid;

    if (verbose > 0) {
      fprintf(stdout,"Setting box constraints.\n");
      fflush(stdout);
    }

    datei=fopen("box.con", "rt");
    if(datei==NULL){
      fprintf(stderr,"ERROR: box.con not provided\n");
      return;
    }
    //read the constraint ID (1=fixed ...)
    fscanf(datei,"%i", &cid);
    if (verbose > 0) {
      fprintf(stdout,"constraint type = %i\n", cid);
      fflush(stdout);
    }

    //read the following box coordinates	
    //xmin, ymin, zmin
    //xdim, ydim, zdim
		
    fscanf(datei,"%i %i %i\n", &xmin, &ymin, &zmin);
    if (verbose > 0) {
      fprintf(stdout,"xmin=%i ymin=%i zmin=%i\n", xmin, ymin, zmin);
      fflush(stdout);
    }

    fscanf(datei,"%i %i %i\n", &xdim, &ydim, &zdim);
    if (verbose > 0) {
      fprintf(stdout,"xdim=%i ydim=%i zdim=%i\n", xdim, ydim, zdim);
      fflush(stdout);
    }

    int flag=0;
    for (int i=1; i<=nv; i++){
           		
      vertex *v = (vertex *)VGraphGetNode(vertices, i);
      assert(v != NULL);
      if(v->x >=xmin && v->x <=xmin+xdim){
	if(v->y >=ymin && v->y <=ymin+ydim){
	  if(v->z >=zmin && v->z <=zmin+zdim){
	    flag++;			
	    VSetPixel(tI, 0, 0, i-1, cid);
	  }
	}
      }
    }
    if(flag == 0){
      fprintf(stderr,"WARNING: no constraint set\n");
      return;
    }
    else {
      if (verbose > 0) {
	fprintf(stdout,"Number of fixed (constraint type %d) nodes = %i \n", cid, flag);
	fflush(stdout);
      }
    }
		
  }
  //if constraints are to be defined on a plane (xy, yz, zx)
  else if(constype==plane){
    FILE *datei;
    int plane_pos;	//x, y or z position of plane
    int plane_dir;	//1:xy  2: yz  3: zx
    int cid;
		
    if (verbose > 0) {
      fprintf(stdout,"Setting plane constraints.");
      fflush(stdout);
    }

    datei=fopen("plane.con", "rt");
    if(datei==NULL){
      fprintf(stderr,"ERROR: plane.con not provided\n");
      return;
    }
    //read the constraint ID (1=fixed ...)
    fscanf(datei,"%i", &cid);
		
    //read the orthogonal plane parameters
    // 1. plane_dir
    // 2. plane_pos
		
    fscanf(datei,"%i\n",  &plane_dir );
    fscanf(datei,"%i\n",  &plane_pos );
		
    //check if specified pos is part of the coordinates
    int flag=0;
    for (int i=1; i<=nv; i++){
           		
      vertex *v = (vertex *)VGraphGetNode(vertices, i);
      assert(v != NULL);
      if(plane_dir == 1){
	if(v->z == plane_pos){
	  flag=1;
	  break;
	}
      }
      if(plane_dir == 2){
	if(v->x == plane_pos){
	  flag=1;
	  break;
	}
      }
      if(plane_dir == 3){
	if(v->y == plane_pos){
	  flag=1;
	  break;
	}
      }
			
    }
    if(!flag){
      if(verbose>0){
	fprintf(stderr, "WARNING: specified plane ");
	fprintf(stderr, "position not part of coordinate values!\n");
	fprintf(stderr, "Searching for next neighbor ...\n");
	fflush(stderr);
      }		
		
      //search for coordinate value closest to plane_pos
      float min=1000;
      float dist, plane=0;
      for (int i=1; i<=nv; i++){
			
	vertex *v = (vertex *)VGraphGetNode(vertices, i);
	assert(v != NULL);				
				
	if(plane_dir == 1){
	  dist=v->z - plane_pos;
	  if(dist<0)
	    dist=dist*(-1);
	  if (dist<min){
	    min=dist;
	    plane=v->z;
	  }
	}
	if(plane_dir == 2){
	  dist=v->x - plane_pos;
	  if(dist<0)
	    dist=dist*(-1);
	  if (dist<min){
	    min=dist;
	    plane=v->x;
	  }	
	}
	if(plane_dir == 3){
	  dist=v->y - plane_pos;
	  if(dist<0)
	    dist=dist*(-1);
	  if (dist<min){
	    min=dist;
	    plane=v->y;
	  }	
	}	
      }
      plane_pos=(int)plane;
    }
    //set attribute values to defined value
    for (int i=1; i<=nv; i++){
			
      vertex *v = (vertex *)VGraphGetNode(vertices, i);
      assert(v != NULL);
			
      if(plane_dir == 1){
	if(v->z == plane_pos){
	  VSetPixel(tI, 0, 0, i-1, cid);
	}
      }
      if(plane_dir == 2){
	if(v->x == plane_pos){
	  VSetPixel(tI, 0, 0, i-1, cid);
	}
      }
      if(plane_dir == 3){
	if(v->y == plane_pos){
	  VSetPixel(tI, 0, 0, i-1, cid);
	}
      }
    }
		
	
  }
  else if(constype==snodes){
		
    FILE *datei;
    int nocn;	//number of constrained nodes
    int cid;		//constrained ID
	
    if (verbose > 0) {
      fprintf(stdout,"Setting constraints on a list of nodes.");
      fflush(stdout);
    }
		
    datei=fopen("snodes.con", "rt");
    if(datei==NULL){
      fprintf(stderr,"ERROR: snodes.con not provided\n");
      return;
    }
    //read the constraint ID (1=fixed ...)
    fscanf(datei,"%i", &cid);
		
    //read the number of constrained nodes
    fscanf(datei,"%i\n",  &nocn );
		
    float *x=(float*)malloc(nocn*sizeof(float));
    float *y=(float*)malloc(nocn*sizeof(float));
    float *z=(float*)malloc(nocn*sizeof(float));
		
    if(x==NULL || y==NULL || z==NULL){
      fprintf(stderr, "ERROR: Probs with x,y,z allocation\n");
      return;
    }
    //read the x, y, z coordinates 	
    for (int i=0; i<nocn; i++){
      fscanf(datei,"%f %f %f\n", &x[i], &y[i], &z[i]);
    }
    int count=0;
    for (int i=1; i<=nv; i++){
           		
      vertex *v = (vertex *)VGraphGetNode(vertices, i);
      assert(v != NULL);
      for (int j=0; j<nocn; j++){
	if(x[j] == v->x && y[j] == v->y && z[j] == v->z){
	  count++;
	  VSetPixel(tI, 0, 0, i-1, cid);
	}
      }
    }
    if(count!=nocn){
      fprintf(stderr, "WARNING: At least one constrained node not in list !!\n");
      fflush(stderr);
      return;
    }
  }
	
  return;
	
}



void gridGenerator::WriteForceImage(VImage fI, int nv)
{

  for (int i=0; i<nv; i++){
    VFloat v=0.0;
    VSetPixel(fI, 0, 0, i, v);
    VSetPixel(fI, 1, 0, i, v);
    VSetPixel(fI, 2, 0, i, v);
  }
}



void gridGenerator::WriteSourceNodes(VImage iv, int nv)
{
	
  for (int i=0; i<nv; i++){	
    VBit v=0;
    VSetPixel(iv, 0, 0, i, v);
  }
}




void gridGenerator::WriteInitialPart(VImage pe, VImage pn, int ne, int nv)
{
  if(npart > nv || npart > ne) {
    fprintf(stderr, "WARNING: number of partitions (%d) to large, adjusting to %d\n", npart, min(ne,nv));
    npart = min(ne,nv);
  } 
 
  if (npart==1) {
    // nodes
    for (int i=0; i<nv; i++){	
      VUByte v=0;
      VSetPixel(pn, 0, 0, i, v);
    }
    // elements
    for (int i=0; i<ne; i++){	
      VUByte v=0;
      VSetPixel(pe, 0, 0, i, v);
    }
  }
  else {
    // initial node partition
    int npf	= nv / npart;	// nodes per file
    VUByte pid  = 0;	// partition ID
		
    for (int i=1; i<=nv; i++){
      if (i%npf != 0)
	VSetPixel(pn, 0, 0, i-1 , pid);
      else{
	VSetPixel(pn, 0, 0, i-1 , pid);
	pid++;
	if (pid==npart) pid--;
      }
    }
		
    // initial element partition	
    int	epf  = ne / npart;	// nodes per file
    VUByte 	eid  = 0;		// partition ID
    int 	ecount=0;
    for (int i=1; i<=ne; i++){
      if(i%epf != 0){
	ecount++;
	VSetPixel(pe, 0, 0, i-1, eid);
      }
      else{
	VSetPixel(pe, 0, 0, i-1, eid);
	eid++;
	if (eid==npart) eid--;
      }
    }
  }				      	
}

