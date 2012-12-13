#ifndef SIMBIO_GENERATE_GRID_H
#define SIMBIO_GENERATE_GRID_H

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


#include <vector>
#include <string>
#include <iostream>

#include <vista.h>

#include "vec3.h"
#include "mesh.h"
#include "generateGrid-util.h"

#include "materialDB.h"
#include "marching-tets.h"

using namespace std; 

/*! \file generateGrid.h
 
   \author Guntram Berti  

   $Id$

   \todo do not overwrite src image (use high resolution for marching tet)
   \todo Preserve material file settings when mapping images with voxel values > 1 byte
   to one byte image
   \todo Adaptive octree generation depending on geometric location
   \todo Enable output of octree as mesh (optionally with constraints for hanging nodes)
*/

enum primitiveType      { nogrid, tetra5, tetra6a, tetra6b, cube };
enum smoothType         { none, bshift, marching };            // boundary smoothing type
enum consType           { no, box, plane, snodes };            // constraints (fixity boundary conditions)  
enum formType           { vm, hf, ipe, gmv, ascii, gpp, complex2d, dx}; // output formats

class tetra_table;

class gridGenerator  {

  VImage                src;                    // source image
  int                   nx, ny, nz;             // source image dimensions
   
  VImage                oct;                    // intermediate octree
  VImage                octnodes;               // octree vertices 
  VGraph                vertices;               // vertex graph
  int                   nv;                     // number of vertices in the 'vertices' VGraph
  VGraph                primitives;             // primitive (cell) graph, for volume mesh generation
  VGraph                faces;                  // faces graph, for surface mesh generation

  bool                  do_output_fields;       // flag whether fields below are written
                                                // (material is always written)


  // mesh generation parameters
  int           mindim, maxdim;           // global minimum and maximum cell size
  int           max_maxdim;               // maximal maxdim over all materials and interfaces
  int           min_maxdim;               // minimal maxdim over all materials and interfaces
  std::string   mat_file_nm;              // file containing material parameters

  primitiveType prim_type;                // target primitive type

  double        shift;                    // shift factor in [0,0.5)
  bool          interp;                   // interpolation flag
  smoothType    btype;                    // type of boundary smoothing
  bool          surf_arg;                 // surface mesh generation (true) or volume m.g. (false)

  vertex     **cache;                     // vertex cache
  void         initCache();
  void         deleteCache();

  // additional mesh attibutes
  VUByte        npart;                    // number of partitions
  int           constype;                 // constraint type (box, ..

  // output parameters
  int           verbose;                  // output level
  bool          implicit_links;           // do not link vertices explicitely
  formType      out_arg;                  // output format
  
  bool          use_simbio_mat_ids;       // use the SimBio material ID convention
  materialDB    MatDB;                    // used for SimBio material strings / IDs



  //------------  material related types, data and functions --------------

  typedef int voxel_type;      // typically grey value [0..MAXLABEL]
  typedef int material_type;   // internal contiguous material ID [0..num_of_materials]
  typedef int simbio_material; // material ID of SimBio

  enum {undefined_material = -1};
  voxel_type    getCellVoxel(int x, int y, int z) const { return VPixel(src, z, y, x, VUByte);}

  // should use a different image for materials!
  // getCellMaterial works only if setCellMaterial called before on src!!!
  void          setCellMaterial(int x, int y, int z, material_type m) const
  { VPixel(src, z, y, x, VUByte) = m;}

  material_type getCellMaterial(octcell const& c) const
  { return getCellMaterial(c.x(), c.y(), c.z());}   

  material_type getCellMaterial(int x, int y, int z) const
  { 
    material_type res = VPixel(src, z, y, x, VUByte); 
    if(res >= (material_type)num_of_materials())
      fprintf(stderr, "WARNING: undefined material %d at %d,%d,%d!\n",res,x,y,z);
    return res;
  }

  material_type get_material(voxel_type vox) const { 
    if(materials[vox] == undefined_material)
      fprintf(stderr, "WARNING: materials[%d] undefined!\n",vox);
    return materials[vox];
  }
  unsigned      num_of_materials() const { return num_of_materials_;}
  material_type get_material(std::string const& mname) const; 
  std::string  const& material_name(material_type m) const;          

  int num_of_materials_;
  vector<material_type> materials; // map voxel values to materials
  vector<material_type> material_of_primitive;

  //! maximal size of cell in material
  vector<int>           maxdim_of_material; // [0 .. num_of_materials]

  //! maximal size of cell at material interface (corr. to global mindim)
  vector<vector<int> >  maxdim_of_interface;

  struct meshing_type {
    enum tag { undefined=0, deflt=1, none=2, 
               uniform_cube=3, uniform_tetra5=4,
               nonuniform_tets=5, nonuniform_hybrid=6, nonuniform_cube=7,
               marching_tet=8, size_=9};
    static std::string names[size_];
  };

  // meshing parameters per material
  struct material_description {
    typedef material_description self;
    string            name;
    int               maxdim;
    voxel_type        low_voxel;
    voxel_type        high_voxel;
    meshing_type::tag mesh_type;
    float             weight; // weighting factor for determining majority  material in cell

    material_description(std::string const& name_ = "", 
                         int maxdim_ = 0, 
                         voxel_type lv = 0, 
                         voxel_type hv = 0, 
                         meshing_type::tag mt = meshing_type::undefined,
                         float w = 1.0) 
      : name(name_),  maxdim(maxdim_), 
        low_voxel(lv), high_voxel(hv), 
        mesh_type(mt), weight(w)
    {}
          
    void read(std::istream & in);
    friend std::istream& operator>>(std::istream& in, self& m) { m.read(in); return in;}
    friend std::ostream& operator<<(std::ostream& out, self const& m)
    { return(out << m.name << ' ' << m.maxdim << ' ' << m.low_voxel 
             << ' ' << m.high_voxel << ' ' << m.mesh_type << ' ' << m.weight);
    }

  };

  // meshing parameters per material interface
  struct interface_description {
    typedef interface_description self;
    string mat1;
    string mat2;
    int maxdim;
    void read(std::istream& in);
    friend std::istream& operator>>(std::istream& in, self& ifd)
    { ifd.read(in); return in;}
    friend std::ostream& operator<<(std::ostream& out, self const& ifd)
    { return (out << ifd.mat1 << ' ' <<  ifd.mat2 << ' ' << ifd.maxdim);}
  };


  //! list of user-defined materials
  std::vector<material_description> material_list;

  //! list of user-defined material interfaces
  std::vector<interface_description> material_interfaces;


  //! preferred algorithm to mesh material
  std::vector<meshing_type::tag>  meshing_algorithm;


  int get_maxdim(material_type m) const { return maxdim_of_material[m];}
  int get_maxdim(material_type m1, material_type m2) const
  { return maxdim_of_interface[m1][m2];}
  //! get maximal dimension of super-voxel
  int  get_maxdim(int ox, int oy, int z, int d, 
                  material_type& max_mat) const;
  bool mustTesselateCell(int x, int y, int z) const
  { 
    material_type m = getCellMaterial(x,y,z);
    return (getCellDim(x,y,z) > 0 
            && meshing_algorithm[m] != meshing_type::none 
            && !(meshing_algorithm[m] == meshing_type::deflt 
                 && this->prim_type == nogrid));
  }

  /*! prepare a list of materials (with frequency) occuring 
    in super-voxel (ox,oy,oz,d) 
  */
  void  list_materials(int ox, int oy, int oz, int d, 
                       vector<unsigned> & mat) const;

  void  initMaterials();
  void  print_material_parameters() const;



  //---------- functions to convert image into octree  ------------------

  bool          generateCells();
  void          subsampleCells(int d);
  int           subsampleCell(int x, int y, int z, int d);
  bool          isHomogenousCell(int x, int y, int z, int d);
  void          compactCells();
  void          projectCorners();




  //--------------- octree functionality ------------

  int           getCellDim(int x, int y, int z) const
  { return (int)VPixel(oct, z, y, x, VUByte); };

  void          setCellDim(int x, int y, int z, int d)
  { VPixel(oct, z, y, x, VUByte) = (unsigned char)d; };

  bool          cellInUse(int x, int y, int z) const
  { return (getCellDim(x, y, z) != 0 
            && VPixel(src, z, y, x, VUByte) != 0)? true: false; }

  void  join_cells(int ox, int oy, int oz, int d);

  //! is x,y,z a corner of an octree cell?
  bool isVertex(int x, int y, int z) const;

  //! linear number of vertex p (in cache)
  unsigned int  pos(vertex &p)
  {
    // this does only work for vertices with integer coordinates!
    int px = (int)p.x; assert ( (float)px == p.x);
    int py = (int)p.y; assert ( (float)py == p.y);
    int pz = (int)p.z; assert ( (float)pz == p.z);
    return pz*(nx+1)*(ny+1) + py*(nx+1) + px; 
  }
  // get location on integer lattice of vertex v.
  // NOTE: this fails for vertices inserted during the marching tets algorithm,
  // which do NOT lie on an integer lattice!
  void  getVertexIntegerCoords(int v, int& x, int& y, int& z) const;




  //----------------- functions to convert octree into grid --------------

  bool          generateGrid();
  void          tesselateCell(int x, int y, int z);
  void          tesselateCube(octcell const& c); // int *plist, bool par);




  //---------- non-uniform mesh generation types, data and functions ----------

  typedef std::vector<octcell_tesselation> cell_tesselation_table;

  //! keep a link between octree cells and elements they generated.
  cell_tesselation_table tesselated_cells;


  //! normal direction of facet 
  enum facet_dir {xdir = 0, ydir = 1, zdir = 2}; 
  //! true iff we do non-uniform meshing
  bool nonuniform_meshing() const { return max_maxdim > min_maxdim;}

  //! 3D checkerbord color (0 or 1) 
  int  parity(int x, int y, int z, int d) const;

  //! tesselate facet into 2 triangles, according to parity of x0,y0,z0
  void get_diagonal_triangulation(vector<triangle> & triangle_list,
                                  int x0, int y0, int z0,
                                  int x1, int y1, int z1,
                                  int d,  facet_dir dir);
  //! tesselate facet into triangles, depending on pattern of vertices
  void get_facet_triangulation(vector<triangle> & triangle_list,
                               int x0, int y0, int z0,
                               int x1, int y1, int z1,
                               int d,  facet_dir dir);



  //--------- handling of cell primitives  -------------

  void          linkCube(int eid);
  unsigned int  addCube(int *plist);
  unsigned int  addTetra(int i0, int i1, int i2, int i3);
  unsigned int  addTriangle(int, int, int);
  unsigned int  addVertex(vertex &v);
  unsigned int  addVertex(int x, int y, int z) 
  { vertex v(x,y,z); return addVertex(v);}
  // add per-cell primitives to primitives graph 
  bool          addPrimitives();


        

  //-------- functions for marching tetrahedra algorithm -------

  material_type mat_background, mat_body;
  // remove vertices not belonging to surface faces (if surface meshing is active)
  bool          removeNonSurfaceVertices(); 
  // interpolate the binary material classifications of cells to vertices
  bool          interpolateOnVertices(VImage f_material) const;
  bool          marchingTetras();
  bool          makeCuts(std::vector<int> const& t,
			 marching_tets::cut_vertex_map & cuts,
			 VImage f_material,
			 int * pat);



  //-------- data and functions for boundary smoothing --------

  int                 smooth_iters;
  std::string         surface_smooth_weight_file_nm;
  std::string         volume_smooth_weight_file_nm;
  std::vector<double> surface_smooth_weights;
  std::vector<double> volume_smooth_weights;

  void read_smoothing_weights(std::string const& filenm, std::vector<double> & weights);

  // smooth boundary of tetrahedral meshes
  bool  smoothMarchingInterfaces();


  // functions for shifting boundary vertices of uniform hexahedral meshes
  bool          shiftBoundaryVertices();



  //---------------  output functions ---------------------

  bool          write(FILE *outf);
  bool          doOutputFields() const { return do_output_fields;}
public:
  void          setOutputFields() { do_output_fields = true;}
  void          noOutputFields()  { do_output_fields = false;}
private:
  //----  write attributes images according to format conventions ----

  void          WriteMaterialPid(VImage mp, int ne);
  void          WriteSourceNodes(VImage iv, int nv);
  void          WriteTagImage   (VImage tI, int nv);
  void          WriteForceImage (VImage fI, int nv);
  void          WriteInitialPart(VImage pe, VImage pn, int ne, int nv);
        
public:
  gridGenerator(VImage s, 
		primitiveType g, 
		int min, int max, 
		smoothType b,
		double d, 
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
		std::string const& volume_smooth_weight_file_);

  ~gridGenerator();     
  bool          work(FILE *outf);
};


#endif
