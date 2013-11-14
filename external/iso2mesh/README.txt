----------------------------------------------------------------------
= iso2mesh: an image-based 3D surface and volumetric mesh generator  =
----------------------------------------------------------------------

Author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>
      Martinos Center for Biomedical Imaging
      Massachusetts General Hospital (Harvard Medical School)
      Bldg. 149, 13th St., Charlestown, MA 02148
Version: 1.0.1 (Mapo Tofu - Update 1)
License: GPL v2 or later (see COPYING) 
      (this license does not cover the binaries under the bin/ 
       directory, see Section III for more details)
URL: http://iso2mesh.sf.net


== Table of Content ==
<toc>

== # Introduction ==

"Iso2mesh" is a Matlab/Octave-based mesh generation toolbox
designed for easy creation of high quality surface and 
tetrahedral meshes from 3D volumetric images. It contains 
a rich set of mesh processing scripts/programs, functioning 
independently or interfacing with external free meshing
utilities. Iso2mesh toolbox can operate directly on 3D
binary, segmented or gray-scale images, such as those
from MRI or CT scans, making it particularly suitable 
for multi-modality medical imaging data analysis or 
multi-physics modeling.

Creating high-quality surface or tetrahedral meshes from 
volumetric images has been a challenging task where 
very limited software and resources could be found. 
Commercial tools, such as Mimics and Amira, are both 
expensive and limited in functionalities. Iso2mesh
was developed as a free alternative to these expensive 
commercial tools and provides researchers a highly
flexible, modularized and streamlined image-based mesh 
generation pipeline. It defines simple interfaces to
perform a wide varieties of meshing tasks, ranging
from 3D volumetric image pre-processing (hole-filling, 
thinning and thickening), surface mesh modeling 
(extraction, remeshing, repairing, and smoothing) to 
volumetric mesh creation. Iso2mesh is cross-platform
and is compatible with both Matlab and GNU Octave 
(a free Matlab clone).

The details of this toolbox can be found in the following
reference:

*Qianqian Fang and David Boas, "Tetrahedral mesh generation from volumetric binary and \
 gray-scale images," Proceedings of IEEE International Symposium on Biomedical Imaging \
 (ISBI 2009), pp. 1142-1145, 2009

== # Function List ==

=== # Streamlined mesh generation - shortcuts ===

==== function [node,elem,face]=v2m(img,isovalues,opt,maxvol,method) ====
 [node,elem,face]=v2m(img,isovalues,opt,maxvol,method)
 volumetric mesh generation from binary or gray-scale volumetric images
 shortcut for vol2mesh
 inputs and outputs are similar to those defined in vol2mesh

==== function [no,el,regions,holes]=v2s(img,isovalues,opt,method) ====
 [no,el,regions,holes]=v2s(img,isovalues,opt,method)
 surface mesh generation from binary or gray-scale volumetric images
 shortcut for vol2surf
 inputs and outputs are similar to those defined in vol2surf

==== function [node,elem,face]=s2m(v,f,keepratio,maxvol) ====
 [node,elem,face]=s2m(v,f,keepratio,maxvol)
 volumetric mesh generation from a closed surface, shortcut for surf2mesh
 inputs and outputs are similar to those defined in surf2mesh

==== function img=s2v(node,face,div) ====
 img=s2v(node,face,div)
 shortcut for surf2vol, coverting a surface to a volumetric image
 input:
	 node: node list of the triangular surface, 3 columns for x/y/z
	 face: triangle node indices, each row is a triangle
	 div:  division number along the shortest edge of the mesh (resolution)
              if not given, div=50
 output:
	 img: a volumetric binary image at position of ndgrid(xi,yi,zi)

==== function newnode=sms(node,face,iter,alpha,method) ====
 newnode=sms(node,face,iter,useralpha,method)
 simplified version of surface mesh smoothing
 input:
    node:  node coordinates of a surface mesh
    face:  face element list of the surface mesh
    iter:  smoothing iteration number
    alpha: scaler, smoothing parameter, v(k+1)=alpha*v(k)+(1-alpha)*mean(neighbors)
    method: same as in smoothsurf, default is 'laplacianhc'
 output:
    newnode: output, the smoothed node coordinates
=== # Streamlined mesh generation ===

==== function [node,elem,face,regions]=vol2mesh(img,ix,iy,iz,opt,maxvol,dofix,method,isovalues) ====
 [node,elem,face,regions]=vol2mesh(img,ix,iy,iz,opt,maxvol,dofix,method,isovalues)
 convert a binary (or multi-valued) volume to tetrahedral mesh
 input:
	 img: a volumetric binary image
	 ix,iy,iz: subvolume selection indices in x,y,z directions
	 opt: as defined in vol2surf.m
	 maxvol: target maximum tetrahedral elem volume
	 dofix: 1: perform mesh validation&repair, 0: skip repairing
	 method: 'cgalsurf' or omit: use CGAL surface mesher
		 'simplify': use binsurface and then simplify
		 'cgalmesh': use CGAL 3.5 3D mesher for direct mesh generation [new]
		 generally speaking, 'cgalmesh' is the most robust path
		 if you want to product meshes from binary or multi-region
		 volumes, however, its limitations include 1) only accept 
		 uint8 volume, and 2) can not extract meshes from gray-scale
		 volumes. If ones goal is to process a gray-scale volume,
		 he/she should use the 'cgalsurf' option. 'simplify' approach
		 is not recommended unless other options failed.
	 isovalues: a list of isovalues where the levelset is defined
 output:
	 node: output, node coordinates of the tetrahedral mesh
	 elem: output, element list of the tetrahedral mesh, the last 
	       column is the region ID
	 face: output, mesh surface element list of the tetrahedral mesh
	       the last column denotes the boundary ID
    region: optional output. if opt.autoregion is set to 1, region
          saves the interior points for each closed surface component

==== function [no,el,regions,holes]=vol2surf(img,ix,iy,iz,opt,dofix,method,isovalues) ====
 [no,el,regions,holes]=vol2surf(img,ix,iy,iz,opt,dofix,method,isovalues)
 converting a 3D volumetric image to surfaces
 input:
	 img: a volumetric binary image; if img is empty, vol2surf will
	      return user defined surfaces via opt.surf if it exists
	 ix,iy,iz: subvolume selection indices in x,y,z directions
	 opt: function parameters
	   if method is 'cgalsurf' or 'cgalpoly':
	     opt=a float number>1: max radius of the Delaunay sphere(element size) 
	     opt.radbound: same as above, max radius of the Delaunay sphere
	     opt.distbound: maximum deviation from the specified isosurfaces
	     opt(1,2,...).radbound: same as above, for each levelset
	   if method is 'simplify':
	     opt=a float number<1: compression rate for surf. simplification
	     opt.keeyratio=a float less than 1: same as above, same for all surf.
	     opt(1,2,..).keeyratio: setting compression rate for each levelset
	   opt(1,2,..).maxsurf: 1 - only use the largest disjointed surface
				0 - use all surfaces for that levelset
          opt(1,2,..).side: - 'upper': threshold at upper interface
                              'lower': threshold at lower interface
	   opt(1,2,..).maxnode: - the maximum number of surface node per levelset
	   opt(1,2,..).holes: user specified holes interior pt list
	   opt(1,2,..).regions: user specified regions interior pt list
	   opt(1,2,..).surf.{node,elem}: add additional surfaces
	   opt(1,2,..).{A,B}: linear transformation for each surface
	   opt.autoregion: if set to 1, vol2surf will try to determine 
              the interior points for each closed surface automatically
	 dofix: 1: perform mesh validation&repair, 0: skip repairing
	 method: - if method is 'simplify', iso2mesh will first call
		   binsurface to generate a voxel-based surface mesh and then
		   use meshresample/meshcheckrepair to create a coarser mesh;
		 - if method is 'cgalsurf', iso2mesh will call the surface
		   extraction program from CGAL to make surface mesh
		 - if method is not specified, 'cgalsurf' is assumed by default
	 isovalues: a list of isovalues where the levelset is defined
 output: 
	 no:  list of nodes on the resulting suface mesh, 3 columns for x,y,z
	 el:  list of trianglular elements on the surface, [n1,n2,n3,region_id]
	 regions: list of interior points for all sub-region, [x,y,z]
	 holes:   list of interior points for all holes, [x,y,z]

==== function [node,elem,face]=surf2mesh(v,f,p0,p1,keepratio,maxvol,regions,holes,forcebox) ====
 [node,elem,face]=surf2mesh(v,f,p0,p1,keepratio,maxvol,regions,holes,forcebox)
 create quality volumetric mesh from isosurface patches
 input parameters:
      v: input, isosurface node list, dimension (nn,3)
         if v has 4 columns, the last column specifies mesh density near each node
      f: input, isosurface face element list, dimension (be,3)
      p0: input, coordinates of one corner of the bounding box, p0=[x0 y0 z0]
      p1: input, coordinates of the other corner of the bounding box, p1=[x1 y1 z1]
      keepratio: input, percentage of elements being kept after the simplification
      maxvol: input, maximum tetrahedra element volume
      regions: list of regions, specifying by an internal point for each region
      holes: list of holes, similar to regions
      forcebox: 1: add bounding box, 0: automatic
 outputs:
      node: output, node coordinates of the tetrahedral mesh
      elem: output, element list of the tetrahedral mesh
      face: output, mesh surface element list of the tetrahedral mesh 
             the last column denotes the boundary ID

==== function img=surf2vol(node,face,xi,yi,zi) ====
 img=surf2vol(node,face,xi,yi,zi)
 convert a triangular surface to a shell of voxels in a 3D image
 input:
	 node: node list of the triangular surface, 3 columns for x/y/z
	 face: triangle node indices, each row is a triangle
	 xi,yi,zi: x/y/z grid for the resulting volume
 output:
	 img: a volumetric binary image at position of ndgrid(xi,yi,zi)
=== # iso2mesh main function backend ===

==== function [node,elem]=binsurface(img,nface) ====
 [node,elem]=binsurface(img,nface)
 fast isosurface extraction from 3D binary images
 input: 
   img:  a 3D binary image
   nface: nface=3 or ignored - for triangular faces, 
          nface=4 - square faces
          nface=0 - return a boundary mask image via node
 output
   elem: integer array with dimensions of NE x nface, each row represents
         a surface mesh face element 
   node: node coordinates, 3 columns for x, y and z respectively
 the outputs of this subroutine can be easily plotted using 
     patch('Vertices',node,'faces',elem,'FaceVertexCData',node(:,3),
           'FaceColor','interp');
 if the surface mesh has triangular faces, one can plot it with
     trisurf(elem,node(:,1),node(:,2),node(:,3))

==== function [node,elem,face]=cgalv2m(vol,opt,maxvol) ====
 [node,elem,face]=cgalv2m(vol,opt,maxvol)
 wrapper for CGAL 3D mesher (CGAL 3.5 or up)
 convert a binary (or multi-valued) volume to tetrahedral mesh
 http://www.cgal.org/Manual/3.5/doc_html/cgal_manual/Mesh_3/Chapter_main.html
 input:
	 vol: a volumetric binary image
	 ix,iy,iz: subvolume selection indices in x,y,z directions
	 opt: parameters for CGAL mesher, if opt is a structure, then
	     opt.radbound: defines the maximum surface element size
	     opt.angbound: defines the miminum angle of a surface triangle
	     opt.distbound: defines the maximum distance between the 
		 center of the surface bounding circle and center of the 
		 element bounding sphere
	     opt.reratio:  maximum radius-edge ratio
	     if opt is a scalar, it only specifies radbound.
	 maxvol: target maximum tetrahedral elem volume
 output:
	 node: output, node coordinates of the tetrahedral mesh
	 elem: output, element list of the tetrahedral mesh, the last 
	      column is the region id
	 face: output, mesh surface element list of the tetrahedral mesh
	      the last column denotes the boundary ID
	      note: each triangle will appear twice in the face list with each
		    one attaches to each side of the interface. one can remove
		    the redundant triangles by unique(face(:,1:3),'rows')

==== function [node,elem,face]=cgals2m(v,f,opt,maxvol) ====
 [node,elem,face]=cgals2m(v,f,opt,maxvol)
 wrapper for CGAL 3D mesher (CGAL 3.5 and newer)
 convert a binary (or multi-valued) volume to tetrahedral mesh
 http://www.cgal.org/Manual/3.5/doc_html/cgal_manual/Mesh_3/Chapter_main.html
 input:
	 v: the node coordinate list of a surface mesh (nn x 3)
	 f: the face element list of a surface mesh (be x 3)
	 opt: parameters for CGAL mesher, if opt is a structure, then
	     opt.radbound: defines the maximum surface element size
	     opt.angbound: defines the miminum angle of a surface triangle
	     opt.distbound: defines the maximum distance between the 
		 center of the surface bounding circle and center of the 
		 element bounding sphere
	     opt.reratio:  maximum radius-edge ratio
	     if opt is a scalar, it only specifies radbound.
	 maxvol: target maximum tetrahedral elem volume
 output:
	 node: output, node coordinates of the tetrahedral mesh
	 elem: output, element list of the tetrahedral mesh, the last 
	      column is the region id
	 face: output, mesh surface element list of the tetrahedral mesh
	      the last column denotes the boundary ID

==== function [node,elem]=vol2restrictedtri(vol,thres,cent,brad,ang,radbound,distbound,maxnode) ====
 [node,elem]=vol2restrictedtri(vol,thres,cent,brad,ang,radbound,distbound,maxnode)
 surface mesh extraction using CGAL mesher
 input:
       vol: a 3D volumetric image
       thres: a scalar as the threshold of of the extraction
       cent: a 3d position (x,y,z) which locates inside the resulting
             mesh, this is automatically computed from vol2surf
       brad: maximum bounding sphere squared of the resulting mesh
       ang: minimum angular constrains of the resulting tranglar elements
            (in degrees)
       radbound: maximum triangle delaunay circle radius
       distbound: maximum delaunay sphere distances
       maxnode: maximum number of surface nodes (even radbound is not reached)
 output:
       node: the list of 3d nodes in the resulting surface (x,y,z)
       elem: the element list of the resulting mesh (3 columns of integers)

==== function img=surf2volz(node,face,xi,yi,zi) ====
 img=surf2volz(node,face,xi,yi,zi)
 convert a triangular surface to a shell of voxels in a 3D image
 along the z-axis
 input:
	 node: node list of the triangular surface, 3 columns for x/y/z
	 face: triangle node indices, each row is a triangle
	 xi,yi,zi: x/y/z grid for the resulting volume
 output:
	 img: a volumetric binary image at position of ndgrid(xi,yi,zi)
=== # iso2mesh primitive meshing functions ===

==== function [node,elem,face]=meshabox(p0,p1,opt,nodesize) ====
 [node,elem,face]=meshabox(p0,p1,opt,maxvol)
 create the surface and tetrahedral mesh of a box geometry
 input: 
   p0:  coordinates (x,y,z) for one end of the box diagnoal
   p1:  coordinates (x,y,z) for the other end of the box diagnoal
   opt: maximum volume of the tetrahedral elements
   nodesize: 1 or a 8x1 array, size of the element near each vertex
 output:
   node: node coordinates, 3 columns for x, y and z respectively
   face: integer array with dimensions of NB x 3, each row represents
         a surface mesh face element 
   elem: integer array with dimensions of NE x 4, each row represents
         a tetrahedron 
 example:
   [node,elem,face]=meshabox([2 3 2],[6 12 15],0.1,1);
   plotmesh(node,elem,'x>4');

==== function [node,face,elem]=meshasphere(c0,r,tsize,maxvol) ====
 [node,face,elem]=meshasphere(c0,r,tsize,maxvol)
 create the surface and tetrahedral mesh of a sphere
 input: 
   c0:  center coordinates (x0,y0,z0) of the sphere
   r:   radius of the sphere
   tsize: maximum surface triangle size on the sphere
   maxvol: maximu volume of the tetrahedral elements
 output
   node: node coordinates, 3 columns for x, y and z respectively
   face: integer array with dimensions of NB x 3, each row represents
         a surface mesh face element 
   elem: integer array with dimensions of NE x 4, each row represents
         a tetrahedron 

==== function [node,face,elem]=meshanellip(c0,rr,tsize,maxvol) ====
 [node,face,elem]=meshanellip(c0,rr,opt)
 create the surface and tetrahedral mesh of an ellipsoid
 input: 
   c0:  center coordinates (x0,y0,z0) of the ellipsoid
   rr:  radii of an ellipsoid, 
        if rr is a scalar, this is a sphere with radius rr
        if rr is a 1x3 or 3x1 vector, it specifies the ellipsoid radii [a,b,c]
        if rr is a 1x5 or 5x1 vector, it specifies [a,b,c,theta,phi]
           where theta and phi are the rotation angles along z and x 
           axes, respectively. Rotation is applied before translation.
   tsize: maximum surface triangle size on the sphere
   maxvol: maximu volume of the tetrahedral elements
 output:
   node: node coordinates, 3 columns for x, y and z respectively
   face: integer array with dimensions of NB x 3, each row represents
         a surface mesh face element 
   elem: integer array with dimensions of NE x 4, each row represents
         a tetrahedron; if ignored, only produces the surface
 example:
   [node,face,elem]=meshanellip([10,10,-10],[30,20,10,pi/4,pi/4],0.5,0.4);
   plotmesh(node,elem,'x>10');axis equal;

==== function [node,face,elem]=meshunitsphere(tsize,maxvol) ====
 [node,face,elem]=meshunitsphere(tsize,maxvol)
 create the surface and/or volumetric mesh of a unit sphere 
 centered at [0 0 0] and radius 1
 input: 
   tsize: maximum size of the surface triangles (from 0 to 1)
   maxvol: maximum volume of the tetrahedron; if one wants to return
           elem without specifying maxvol, maxvol=tsize^3
 output:
   node: node coordinates, 3 columns for x, y and z respectively
   face: integer array with dimensions of NB x 3, each row represents
         a surface mesh face element 
   elem: integer array with dimensions of NE x 4, each row represents
         a tetrahedron. If ignored, this function only produces the surface
 example:
   [node,face]=meshunitsphere(0.05);
   [node,face,elem]=meshunitsphere(0.05,0.01);
   plotmesh(node,elem,'x>0'); axis equal;
=== # Mesh decomposition and query ===

==== function facecell=finddisconnsurf(f) ====
 facecell=finddisconnsurf(f)
 subroutine to extract disconnected surfaces from a 
 cluster of surfaces
 Date: 2008/03/06
 input: 
     f: faces defined by node indices for all surface triangles
 output:
     facecell: separated disconnected surface node indices

==== function openedge=surfedge(f) ====
 openedge=surfedge(f)
 find the edge of an open surface or surface of a volume
 input:
      f: input, surface face element list, dimension (be,3)
 output:
      openedge: list of edges of the specified surface

==== function openface=volface(t) ====
 openface=volface(t)
 find the surface patches of a volume
 input:
      t: input, volumetric element list, dimension (ne,4)
 output:
      openface: list of faces of the specified volume

==== function loops=extractloops(edges) ====
 loops=extractloops(edges)
 extract individual loops from an edge table of a loop
 collection
 input:   
    edges:  two column matrix recording the starting/ending 
             points of all edge segments
 output:
    loops:  output, a single vector separated by NaN, each segment
             is a close-polygon consisted by node IDs

==== function [conn,connnum,count]=meshconn(elem,nn) ====
 [conn,connnum,count]=meshconn(elem,nn)
 create node neighbor list from a mesh
 input:
    elem:  element table of a mesh
    nn  :  total node number of the mesh
 output:
    conn:  output, a cell structure of length nn, conn{n}
           contains a list of all neighboring node ID for node n
    connnum: vector of length nn, denotes the neighbor number of each node
    count: total neighbor numbers

==== function centroid=meshcentroid(v,f) ====
 centroid=meshcentroid(v,f)
 compute the centroids of a mesh defined by nodes and elements
 (surface or tetrahedra) in R^n space
 input:
      v: surface node list, dimension (nn,3)
      f: surface face element list, dimension (be,3)
 output:
      centroid: centroid positions, one row for each element

==== function nodevol=nodevolume(node,elem) ====
 nodevol=nodevolume(node,elem)
 calculate the Voronoi volume of each node in a simplex mesh
 input:
    node:  node coordinates
    elem:  element table of a mesh
 output:
    nodevol:   volume values for all nodes

==== function vol=elemvolume(node,elem,option) ====
 vol=elemvolume(node,elem,option)
 calculate the volume for a list of simplexes
 input:
    node:  node coordinates
    elem:  element table of a mesh
    option: if option='signed', the volume is the raw determinant,
            else, the results will be the absolute values
 output
    vol:   volume values for all elements

==== function [conn,connnum,count]=neighborelem(elem,nn); ====
 [conn,connnum,count]=neighborelem(elem,nn)
 create node neighbor list from a mesh
 input:
    elem:  element table of a mesh
    nn  :  total node number of the mesh
 output:
    conn:  output, a cell structure of length nn, conn{n}
           contains a list of all neighboring elem ID for node n
    connnum: vector of length nn, denotes the neighbor number of each node
    count: total neighbor numbers

==== function facenb=faceneighbors(t,opt) ====
 facenb=faceneighbors(t,opt)
 to find 4 face-neighboring elements of a tetrahedron
 input:
     t: tetrahedron element list, 4 columns of integers
     opt: if opt='surface', return boundary triangle list 
          (should be the same as the face output from v2m)
          otherwise, return the element list for each element:
          each row contains 4 numbers, representing the element
          indices sharing triangular faces [1 2 3],[1 2 4],[1 3 4]
          and [2 3 4] in order, where 1~4 is the node local index.
          if the index is 0, indicating the face has no neighbor
          (i.e. a boundary face)
 output:
     facenb: see opt

==== function edgenb=edgeneighbors(t,opt) ====
 edgenb=edgeneighbors(t,opt)
 to find neighboring triangular elements in a triangule surface
 input:
     t: a triangular surface element list, 3 columns of integers
     opt: if opt='general', return the edge neighbors for a general
          triangular surface: each edge can be shared by more than 2
          triangles; if ignored, we assume all triangles are shared by no
          more than 2 triangles.
 output:
     edgenb: if opt is not supplied, edgenb is a size(t,1) by 3 array with
     each element being the triangle ID of the edge neighbor of that
     triangle. For each row, the order of the neighbors is listed as those
     sharing edges [1 2], [2 3] and [3 1] between the triangle nodes.
     when opt='general', edgenb is a cell array with a length of size(t).
     each member of the cell array is a list of edge neighbors (the order 
     is not defined).

==== function f=maxsurf(facecell) ====
 f=maxsurf(facecell)
 return the surface with the maximum element number (not 
 necessarily in area) from a cell arry of surfaces
 input:
    facecell: a cell array, each element is a face array
 output:
    f: the surface data (node indices) for the surface with most elements

==== function mask=flatsegment(node,edge) ====
 mask=flatsegment(node,edge)
 decompose edge loops into flat segments alone arbitrary planes of the bounding box
 this code is fragile: it can not handle curves with many co-linear
 nodes near the corner point
 input:   
    node:  x,y,z coordinates of each node of the mesh
    edge:  input, a single vector separated by NaN, each segment
           is a close-polygon consisted by node IDs 
 output:
    mask:  output, a cell, each element is a close-polygon 
           on x/y/z plane 

==== function newedge=orderloopedge(edge) ====
 [newedge]=orderloopedge(edge)
 order the node list of a simple loop based on connection sequence
 input: 
        edge: a loop consisted by a sequence of edges, each row 
              is an edge with two integers: start/end node index
 output:
        newedge: reordered edge node list

==== function [X,V,E,F]=mesheuler(face) ====
 [X,V,E,F]=mesheuler(face)
 Euler's charastistics of a mesh
 input: 
   face: a closed surface mesh
 output:
   X: Euler's charastistics
   V: number of vertices 
   E: number of edges
   F: number of faces

==== function seg=bbxflatsegment(node,loop) ====
 seg=bbxflatsegment(node,loop)
 decompose edge loops into flat segments alone x/y/z 
 planes of the bounding box
 input:   
    node:  x,y,z coordinates of each node of the mesh
    loop:  input, a single vector separated by NaN, each segment
             is a close-polygon consisted by node IDs 
 output:
    seg:   output, a single vector separated by NaN, each segment
             is a close-polygon on x/y/z plane 

==== function plane=surfplane(node,face) ====
 plane=surfplane(node,face)
 plane equation coefficients for each face in a surface
 input:
   node: a list of node coordinates (nn x 3)
   face: a surface mesh triangle list (ne x 3)
 output:
   plane: a (ne x 4) array, in each row, it has [a b c d]
        to denote the plane equation as "a*x+b*y+c*z+d=0"

==== function [pt,p0,v0,t,idx]=surfinterior(node,face) ====
 [pt,p0,v0,t,idx]=surfinterior(node,face)
 identify a point that is enclosed by the (closed) surface
 input:
   node: a list of node coordinates (nn x 3)
   face: a surface mesh triangle list (ne x 3)
 output:
   pt: the interior point coordinates [x y z]
   p0: ray origin used to determine the interior point
   v0: the vector used to determine the interior point
   t : ray-tracing intersection distances (with signs) from p0. the
       intersection coordinates can be expressed as p0+t(i)*v0
   idx: index to the face elements that intersect with the ray, order
       match that of t

==== function seeds=surfseeds(node,face) ====
 seeds=surfseeds(node,face)
 calculate a set of interior points with each enclosed by a closed
 component of a surface
 input:
   node: a list of node coordinates (nn x 3)
   face: a surface mesh triangle list (ne x 3)
 output:
   seeds: the interior points coordinates for each closed-surface
          component

==== function quality=meshquality(node,elem) ====
 quality=meshquality(node,elem)
 compute Joe-Liu mesh quality measure of a tetrahedral mesh
 input:
    node:  node coordinates of the mesh (nn x 3)
    elem:  element table of a tetrahedral mesh (ne x 4)
 output
    edge:  edge list; each row is an edge, specified by the starting and
           ending node indices, the total edge number is
           size(elem,1) x nchoosek(size(elem,2),2). All edges are ordered
           by looping through each element first. 

==== function edges=meshedge(elem) ====
 edges=meshedge(elem)
 return all edges in a surface or volumetric mesh
 input:
    elem:  element table of a mesh (support N-d space element)
 output
    edge:  edge list; each row is an edge, specified by the starting and
           ending node indices, the total edge number is
           size(elem,1) x nchoosek(size(elem,2),2). All edges are ordered
           by looping through each element first. 
=== # Mesh processing and reparing ===

==== function [node,elem]=meshcheckrepair(node,elem,opt) ====
 [node,elem]=meshcheckrepair(node,elem,opt)
 check and repair a surface mesh
 input/output:
      node: input/output, surface node list, dimension (nn,3)
      elem: input/output, surface face element list, dimension (be,3)
      opt: options, including
            'duplicated': remove duplicated elements
            'isolated': remove isolated nodes
            'deep': call external jmeshlib to remove non-manifold vertices

==== function newelem=meshreorient(node,elem) ====
 newelem=meshreorient(node,elem)
 reorder nodes in a surface or tetrahedral mesh to ensure all
 elements are oriented consistently
 input:
    node: list of nodes
    elem: list of elements (each row are indices of nodes of each element)
 output:
    newelem: the element list with consistent ordering

==== function elem=removedupelem(elem) ====
 elem=removedupelem(elem)
 remove doubly duplicated elements
 input:
    elem: list of elements (node indices)
 output:
    elem: element list after removing the duplicated elements

==== function [newnode,newelem]=removedupnodes(node,elem) ====
 [newnode,newelem]=removedupnodes(node,elem)
 removing the duplicated nodes from a mesh
 input:
   elem: integer array with dimensions of NE x 4, each row contains
         the indices of all the nodes for each tetrahedron
   node: node coordinates, 3 columns for x, y and z respectively
 output:
   newnode: nodes without duplicates
   newelem: elements with only the unique nodes

==== function [no,el]=removeisolatednode(node,elem) ====
 [no,el]=removeisolatednode(node,elem)
 remove isolated nodes: nodes that are not included in any element
 input:
     node: list of node coordinates
     elem: list of elements of the mesh
 output:
     no: node coordinates after removing the isolated nodes
     el: element list of the resulting mesh

==== function fnew=removeisolatedsurf(v,f,maxdiameter) ====
 fnew=removeisolatedsurf(v,f,maxdiameter)
 remove disjointed surface fragment by using mesh diameter
 input:
    v: list of nodes of the input surface
    f: list of triangles of the input surface
    maxdiameter: maximum bounding box size for surface removal
 ouput:
    fnew: new face list after removing the components smaller than 
          maxdiameter

==== function f=surfaceclean(f,v) ====
 f=surfaceclean(f,v)
 remove surface patches that are located inside 
               the bounding box faces
 input: 
      v: surface node list, dimension (nn,3)
      f: surface face element list, dimension (be,3)  
 output:
      f: faces free of those on the bounding box

==== function eid=getintersecttri(tmppath) ====
 eid=getintersecttri(tmppath)
 get the IDs of self-intersecting elements from tetgen
 call this when tetgen complains about self-intersection
 input: 
   tmppath: working dir, use mwpath('') in most cases
 output:
   eid: an array of all intersecting surface elements, 
     one can read the corresponding node/elem by
     [no,el]=readoff(mwpath('post_vmesh.off'));

==== function elem=delendelem(elem,mask) ====
 elem=delendelem(elem,mask)
 delete elements whose nodes are all edge nodes
 input/output: 
      elem: input/output, surface/volumetric element list
      mask: of length of node number, =0 for internal nodes, =1 for edge nodes
=== # Mesh resampling and optimization ===

==== function [node,elem]=meshresample(v,f,keepratio) ====
 [node,elem]=meshresample(v,f,keepratio)
 resample mesh using CGAL mesh simplification utility
 input:
    v: list of nodes
    f: list of surface elements (each row for each triangle)
    keepratio: decimation rate, a number less than 1, as the percentage
               of the elements after the sampling
 output:
    node: the node coordinates of the sampled surface mesh
    elem: the element list of the sampled surface mesh

==== function [newno,newfc]=remeshsurf(node,face,opt) ====
 [newno,newfc]=remeshsurf(node,face,opt)
 remesh a triangular surface and the output is guaranteed to be
 free of self-intersecting element. This function is similar to 
 meshresample, but it can both downsample or upsample a mesh
 input:
	 node: list of nodes on the input suface mesh, 3 columns for x,y,z
	 face: list of trianglular elements on the surface, [n1,n2,n3,region_id]
	 opt: function parameters
	   opt.gridsize:  resolution for the voxelization of the mesh
	   opt.closesize: if there are openings, set the closing diameter
	   opt.elemsize:  the size of the element of the output surface
	   if opt is a scalar, it defines the elemsize and gridsize=opt/4
 output:
	 newno:  list of nodes on the resulting suface mesh, 3 columns for x,y,z
	 newfc:  list of trianglular elements on the surface, [n1,n2,n3,region_id]

==== function p=smoothsurf(node,mask,conn,iter,useralpha,usermethod,userbeta) ====
 p=smoothsurf(node,mask,conn,iter,useralpha,usermethod,userbeta)
 smoothing a surface mesh
 input:
    node:  node coordinates of a surface mesh
    mask:  flag whether a node is movable: 0 movable, 1 non-movable
           if mask=[], it assumes all nodes are movable
    conn:  input, a cell structure of length size(node), conn{n}
           contains a list of all neighboring node ID for node n,
           this can be computed from meshconn function
    iter:  smoothing iteration number
    useralpha: scaler, smoothing parameter, v(k+1)=alpha*v(k)+(1-alpha)*mean(neighbors)
    usermethod: smoothing method, including 'laplacian','laplacianhc' and 'lowpass'
    userbeta: scaler, smoothing parameter, for 'laplacianhc'
 output:
    p: output, the smoothed node coordinates
 recommendations
    Based on [Bade2006], 'Lowpass' method outperforms 'Laplacian-HC' in volume
    preserving and both are significantly better than the standard Laplacian method
    [Bade2006]  R. Bade, H. Haase, B. Preim, "Comparison of Fundamental Mesh 
                Smoothing Algorithms for Medical Surface Models," 
                Simulation and Visualization, pp. 289-304, 2006. 

==== function [no,el,fc,nodemap]=sortmesh(origin,node,elem,ecol,face,fcol) ====
 [no,el,fc]=sortmesh(origin,node,elem,face)
 sort nodes and elements in a mesh so that the indexed
 nodes and elements are closer to each order
 (this may reduce cache-miss in a calculation)
 input:
    origin: sorting all nodes and elements with the distance and
            angles wrt this location, if origin=[], it will be 
            node(1,:)
    node: list of nodes
    elem: list of elements (each row are indices of nodes of each element)
    ecol: list of columns in elem to participate sorting
    face: list of surface triangles (this can be omitted)
    fcol: list of columns in face to participate sorting
 output:
    no: node coordinates in the sorted order
    el: the element list in the sorted order
    fc: the surface triangle list in the sorted order (can be ignored)
    nodemap: the new node mapping order, no=node(nodemap,:)

==== function [newnode,newelem]=mergemesh(node,elem,varargin) ====
 [newnode,newelem]=mergemesh(node,elem,varargin)
 merge two or more tetrahedral meshes or triangular surfaces
 input: 
      node: node coordinates, dimension (nn,3)
      elem: tetrahedral element or triangle surface (nn,3) to (nn,5)
 output:
      newnode: the node coordinates after merging, dimension (nn,3)
      newelem: tetrahedral element or surfaces after merging (nn,4) or (nhn,5)
 note: you can call meshcheckrepair for the output newnode and
 newelem to remove the duplicated nodes or elements
 example:
   [node1,elem1,face1]=meshabox([0 0 0],[10 10 10],1,1);
   [node2,face2,elem2]=meshasphere([5 5 13.1],3,0.3,3);
   [newnode,newelem]=mergemesh(node1,elem1,node2,elem2);
   plotmesh(newnode,newelem);
   figure;
   [newnode,newface]=mergemesh(node1,face1,node2,face2);
   plotmesh(newnode,newface,'x>5');
=== # File I/O ===

==== function saveasc(v,f,fname) ====
 saveasc(v,f,fname)
 save a surface mesh to FreeSurfer ASC mesh format
 input:
      v: input, surface node list, dimension (nn,3)
      f: input, surface face element list, dimension (be,3)
      fname: output file name

==== function savedxf(node,face,elem,fname) ====
 savedxf(node,face,elem,fname)
 save a surface mesh to DXF format
 input:
      node: input, surface node list, dimension (nn,3)
      face: input, surface face element list, dimension (be,3)
      elem: input, tetrahedral element list, dimension (ne,4)
      fname: output file name

==== function saveinr(vol,fname) ====
 saveinr(vol,fname)
 save a surface mesh to INR Format
 input:
      vol: input, a binary volume
      fname: output file name

==== function saveoff(v,f,fname) ====
 saveoff(v,f,fname)
 save a surface mesh to Geomview Object File Format (OFF)
 input:
      v: input, surface node list, dimension (nn,3)
      f: input, surface face element list, dimension (be,3)
      fname: output file name

==== function savesmf(v,f,fname) ====
 savesmf(v,f,fname)
 save a surface mesh to smf format
 input:
      v: input, surface node list, dimension (nn,3)
      f: input, surface face element list, dimension (be,3)
      fname: output file name

==== function savesurfpoly(v,f,holelist,regionlist,p0,p1,fname,forcebox) ====
 savesurfpoly(v,f,holelist,regionlist,p0,p1,fname)
 save a set of surfaces into poly format (for tetgen)
 input:
      v: input, surface node list, dimension (nn,3)
         if v has 4 columns, the last column specifies mesh density near each node
      f: input, surface face element list, dimension (be,3)
      holelist: list of holes, each hole is represented by an internal point
      regionlist: list of regions, similar to holelist
      p0: coordinate of one of the end of the bounding box
      p1: coordinate for the other end of the bounding box
      fname: output file name
      forcebox: non-empty: add bounding box, []: automatic
                if forcebox is a 8x1 vector, it will be used to 
                specify max-edge size near the bounding box corners

==== function savevrml(node,face,elem,fname) ====
 savevrml(node,face,elem,fname)
 save a surface mesh to VRML 1.0 format
 input:
      node: input, surface node list, dimension (nn,3)
      face: input, surface face element list, dimension (be,3)
      elem: input, tetrahedral element list, dimension (ne,4)
      fname: output file name

==== function [node,elem]=readasc(fname) ====
 [node,elem]=readasc(fname)
 read FreeSurfer ASC mesh format
 input:
      fname: name of the asc file
 output:
      node: node positions of the mesh
      elem: element list of the mesh

==== function dat=readinr(fname) ====
 vol=readinr(fname)
 load a volume from an INR file
 input:
      fname: input file name
 output:
      dat: output, data read from the inr file

==== function [node,elem,face]=readmedit(filename) ====
 [node,elem,face]=readmedit(filename)
 read Medit mesh format
 input:
    fname: name of the medit data file
 output:
    node: node coordinates of the mesh
    elem: list of elements of the mesh	    
    face: list of surface triangles of the mesh	    

==== function [node,elem]=readoff(fname) ====
 [node,elem]=readoff(fname)
 read Geomview Object File Format (OFF)
 input:
    fname: name of the OFF data file
 output:
    node: node coordinates of the mesh
    elem: list of elements of the mesh	    

==== function [node,elem]=readsmf(fname) ====
 [node,elem]=readsmf(fname)
 read simple model format (SMF)
 input: 
    fname: name of the	SMF data file
 output:
    node: node coordinates of the mesh
    elem: list of elements of the mesh

==== function [node,elem,face]=readtetgen(fstub) ====
 [node,elem,face]=readtetgen(fstub)
 read tetgen output files
 input:
    fstub: file name stub
 output:
    node: node coordinates of the tetgen mesh
    elem: tetrahedra element list of the tetgen mesh
    face: surface triangles of the tetgen mesh

==== function flag=deletemeshfile(fname) ====
 flag=deletemeshfile(fname)
 delete a given work mesh file under the working directory
 input: 
     fname: specified file name (without path)
 output:
     flag: not used

==== function binname=mcpath(fname) ====
 binname=mcpath(fname)
 get full executable path by prepending a command directory path
 parameters:
 input:
    fname: input, a file name string
 output:
    binname: output, full file name located in the bin directory
    if global variable ISO2MESH_BIN is set in 'base', it will
    use [ISO2MESH_BIN filesep cmdname] as the command full path,
    otherwise, let matlab pass the cmdname to the shell, which
    will search command in the directories listed in system
    $PATH variable.

==== function tempname=mwpath(fname) ====
 tempname=meshtemppath(fname)
 get full temp-file name by prepend working-directory and current session name
 input:
    fname: input, a file name string
 output:
    tempname: output, full file name located in the working directory
    if global variable ISO2MESH_TEMP is set in 'base', it will use it
    as the working directory; otherwise, will use matlab function tempdir
    to return a working directory.
    if global variable ISO2MESH_SESSION is set in 'base', it will be
    prepended for each file name, otherwise, use supplied file name.

==== function savemedit(node,face,elem,fname) ====
 savedmedit(node,face,elem,fname)
 save a surface or tetrahedral mesh to Medit format
 input:
      node: input, surface node list, dimension (nn,3 or 4)
      face: input, surface face element list, dimension (be,3 or 4)
      elem: input, tetrahedral element list, dimension (ne,4 or 5)
      fname: output file name
=== # Volumetric image pre-processing ===

==== function islands=bwislands(img) ====
 islands=bwislands(img)
 return the indices of non-zero elements in a 2D or 3D image
 grouped by connected regions in a cell array
 input:
	 img: a 2D or 3D array
 output:
	 islands: a cell array, each cell records the indices 
		  of the non-zero elements in img for a connected 
		  region (or an island)

==== function resimg=fillholes3d(img,ballsize) ====
 resimg=fillholes3d(img,ballsize)
 close a 3D image with the speicified gap size and then fill the holes
 input:
    img: a 3D binary image
    ballsize: maximum gap size for image closing
 output:
    resimg: the image free of holes
 this function requires the image processing toolbox for matlab/octave

==== function cleanimg=deislands2d(img,sizelim) ====
 cleanimg=deislands2d(img,sizelim)
 remove isolated islands on a 2D image below speicified size limit
 input:
 	img: a 2D binary image
 	sizelim: a integer as the maximum pixel size of a isolated region
 output:
 	cleanimg: a binary image after removing islands below sizelim

==== function cleanimg=deislands3d(img,sizelim) ====
 cleanimg=deislands3d(img,sizelim)
 remove isolated islands for 3D image (for each slice)
 input:
      img: a 3D volumetric image
      sizelim: maximum island size (in pixels) for each x/y/z slice
 output:
      cleanimg: 3D image after removing the islands

==== function imgdiff=imedge3d(binimg,isdiff) ====
 imgdiff=imedge3d(binimg,isdiff)
 Extract the boundary voxels from a binary image
 author: Aslak Grinsted <ag at glaciology.net>
 input: 
   binimg: a 3D binary image
   isdiff: if isdiff=1, output will be all voxels which 
         is different from its neighbors; if isdiff=0 or 
         ignored, output will be the edge voxels of the 
         non-zero regions in binimg
 output:
   imgdiff: a 3D logical array with the same size as binimg
            with 1 for voxels on the boundary and 0 otherwise 

==== function p=internalpoint(v,aloop) ====
 p=internalpoint(v,aloop)
 imperical function to find an internal point
 of a planar polygon
 input:   
    v:     x,y,z coordinates of each node of the mesh
    aloop:  input, a single vector separated by NaN, each segment
             is a close-polygon consisted by node IDs 
 output:
    p:   output, [x y z] of an internal point of aloop

==== function vol=smoothbinvol(vol,layer) ====
 vol=smoothbinvol(vol,layer)
 perform a memory-limited 3D image smoothing
 input:
     vol: a 3D volumetric image to be smoothed
     layer: number of iterations for the smoothing
 output:
     vol: the volumetric image after smoothing

==== function vol=thickenbinvol(vol,layer) ====
 vol=thickenbinvol(vol,layer)
 thickening a binary volume by a given pixel width
 this is similar to bwmorph(vol,'thicken',3) except 
 this does it in 3d and only does thickening for 
 non-zero elements (and hopefully faster)
 input:
     vol: a volumetric binary image
     layer: number of iterations for the thickenining
 output:
     vol: the volume image after the thickening

==== function vol=thinbinvol(vol,layer) ====
 vol=thinbinvol(vol,layer)
 thinning a binary volume by a given pixel width
 this is similar to bwmorph(vol,'thin',n) except 
 this does it in 3d and only does thinning for 
 non-zero elements (and hopefully faster)
 input:
     vol: a volumetric binary image
     layer: number of iterations for the thickenining
 output:
     vol: the volume image after the thickening
=== # Mesh plotting ===

==== function hm=plotmesh(node,varargin) ====
 hm=plotmesh(node,face,elem,opt)
 plot surface and volumetric meshes
 input: 
      node: a node coordinate list, 3 columns for x/y/z; if node has a 
            4th column, it will be used to set the color at each node.
      face: a triangular surface face list; if face has a 4th column,
            it will be used to separate the surface into 
            sub-surfaces and display them in different colors.
      elem: a tetrahedral element list; if elem has a 5th column,
            it will be used to separate the mesh into 
            sub-domains and display them in different colors.
      opt:  additional options for the plotting
            for simple point plotting, opt can be markers
            or color options, such as 'r.', or opt can be 
            a logic statement to select a subset of the mesh,
            such as 'x>0 & y+z<1'; opt can have more than one
            items to combine these options, for example: 
            plotmesh(...,'x>0','r.'); the range selector must
            appear before the color/marker specifier
 in the event where all of the above inputs have extra settings related to 
 the color of the plot, the priorities are given in the following order:
          opt > node(:,4) > elem(:,5) > face(:,4)
 output:
   hm: handle or handles (vector) to the plotted surfaces
 example:
   h=plotmesh(node,'r.');
   h=plotmesh(node,'x<20','r.');
   h=plotmesh(node,face);
   h=plotmesh(node,face,'y>10');
   h=plotmesh(node,face,'facecolor','r');
   h=plotmesh(node,elem,'x<20');
   h=plotmesh(node,elem,'x<20 & y>0');
   h=plotmesh(node,face,elem);
   h=plotmesh(node,face,elem,'linestyle','--');

==== function hm=plotsurf(node,face,varargin) ====
 hm=plotsurf(node,face,opt)
 plot 3D surface meshes
 input: 
      node: node coordinates, dimension (nn,3); if node has a 
            4th column, it will be used to set the color at each node.
      face: triangular surface face list; if face has a 4th column,
            it will be used to separate the surface into 
            sub-surfaces and display them in different colors.
      opt:  additional options for the plotting, see plotmesh
 output:
   hm: handle or handles (vector) to the plotted surfaces
 example:
   h=plotsurf(node,face);
   h=plotsurf(node,face,'facecolor','r');

==== function hm=plottetra(node,elem,varargin) ====
 hm=plottetra(node,elem,opt)
 plot 3D surface meshes
 input: 
      node: a node coordinate list, 3 columns for x/y/z; if node has a 
            4th column, it will be used to set the color at each node.
      elem: a tetrahedral element list; if elem has a 5th column,
            it will be used to separate the mesh into 
            sub-domains and display them in different colors.
      opt:  additional options for a patch object, see plotmesh
 output:
   hm: handle or handles (vector) to the plotted surfaces
 example:
   h=plottetra(node,elem);
   h=plottetra(node,elem,'facealpha',0.5);

==== function [cutpos,cutvalue,facedata]=qmeshcut(elem,node,value,plane) ====
 [cutpos,cutvalue,facedata]=qmeshcut(elem,node,value,plane)
 fast tetrahedral mesh cross-section plot
 input: 
   elem: integer array with dimensions of NE x 4, each row contains
         the indices of all the nodes for each tetrahedron
   node: node coordinates, 3 columns for x, y and z respectively
   value: a scalar array with the length of node numbers, can have 
          multiple columns 
   plane: defines a plane by 3 points: plane=[x1 y1 z1;x2 y2 z2;x3 y3 z3]
 output:
   cutpos: all the intersections of mesh edges by the plane
           cutpos is similar to node, containing 3 columns for x/y/z
   cutvalue: interpolated values at the intersections, with row number
           being the num. of the intersections, column number being the 
           same as "value".
   facedata: define the intersection polygons in the form of patch "Faces"
 the outputs of this subroutine can be easily plotted using 
     patch('Vertices',cutpos,'Faces',facedata,'FaceVertexCData',cutvalue,...
           'FaceColor','interp');

==== function plottetview(session,method) ====
 plottetview(session,method)
 wrapper for tetview to plot the generated mesh
 input:
	 session: a string to identify the output files for plotting, '' for
	          the default session
    method:  method can be 'cgalsurf' (default), 'simplify', 'cgalpoly'
             'cgalmesh' and 'remesh'
=== # Miscellaneous functions ===

==== function valnew=surfdiffuse(node,tri,val,ddt,iter,type1,opt) ====
 valnew=surfdiffuse(node,tri,val,ddt,iter,type1,opt)
 apply a smoothing/diffusion process on a surface
 input:
     node: list of nodes of the surface mesh
     tri: triangular element list of the surface
     val: vector, scalar value for each node
     ddt: diffusion coefficient multiplied by delta t
     iter: iterations for applying the smoothing
     type1: indices of the nodes which will not be updated
     opt: method, 'grad' for gradient based, and 'simple' for simple average
 output:
     valnew: nodal value vector after the smoothing

==== function [node,elem,face]=volmap2mesh(img,ix,iy,iz,elemnum,maxvol,thickness,Amat,Bvec) ====
 [node,elem,face]=volmap2mesh(img,ix,iy,iz,thickness,elemnum,maxvol,A,B)
 convert a binary volume to tetrahedral mesh followed by an Affine transform
 input: 
        img, ix,iy,iz, elemnum and  maxvol: see vol2mesh.m
        thickness: scale z-dimension of the mesh to specified thickness, 
                   if thickness==0, scaling is bypassed
        Amat: a 3x3 transformation matrix
        Bvec: a 3x1 vector
        Amat and Bvec maps the image index space to real world coordnate system by
                   [x,y,z]_new=Amat*[x,y,z]_old+Bvec

==== function isoctave=isoctavemesh ====
 isoctave=isoctavemesh
 determine whether the code is running in octave
 output:
   isoctave: 1 if in octave, otherwise 0

==== function p=getvarfrom(ws,name) ====
 p=getvarfrom(ws,name)
 get variable value by name from specified work-space
 input:
    ws: name of the work-space, for example, 'base'
    name: name string of the variable
 output:
    p: the value of the specified variable, if the variable does not
       exist, return empty array

==== function [t,u,v]=raytrace(p,v,node,face) ====
 [t,u,v]=raytrace(p,v,node,face)
 perform a Havel-styled ray tracing for a triangular surface
 input:
   p: starting point coordinate of the ray
   v: directional vector of the ray
   node: a list of node coordinates (nn x 3)
   face: a surface mesh triangle list (ne x 3)
 output:
   t: signed distance from p to the intersection point
   u: bary-centric coordinate 1 of the intersection point
   v: bary-centric coordinate 2 of the intersection point
      the final bary-centric triplet is [u,v,1-u-v]
  users can find the IDs of the elements intersecting with the ray by
    idx=find(u>=0 & v>=0 & u+v<=1.0);
 Reference: 
  [1] J. Havel and A. Herout, "Yet faster ray-triangle intersection (using 
          SSE4)," IEEE Trans. on Visualization and Computer Graphics,
          16(3):434-438 (2010)
  [2] Q. Fang, "Comment on 'A study on tetrahedron-based inhomogeneous 
          Monte-Carlo optical simulation'," Biomed. Opt. Express, (in
          press)

==== function [a,b,c,d]=getplanefrom3pt(plane) ====
 [a,b,c,d]=getplanefrom3pt(plane)
 define a plane equation ax+by+cz+d=0 from three 3D points
 input: 
    plane: a 3x3 matrix with each row specifying a 3D point (x,y,z)
 output:
    a,b,c,d: the coefficient for plane equation ax+by+cz+d=0

==== function exesuff=getexeext() ====
 exesuff=getexeext()
 get meshing external tool extension names for the current platform
 output:
     exesuff: file extension for iso2mesh tool binaries

==== function exesuff=fallbackexeext(exesuffix, exename) ====
 exesuff=fallbackexeext(exesuffix, exename)
 get the fallback external tool extension names for the current platform
 input:
     exesuffix: the output executable suffix from getexeext
     exename: the executable name
 output:
     exesuff: file extension for iso2mesh tool binaries

==== function [major,minor,patchnum,extra]=iso2meshver ====
 [major,minor,patchnum,extra]=iso2meshver
      or
 v=iso2meshver
 get the version number of iso2mesh toolbox
 output:
    if you ask for a single output:
      v: a string denotes the current version number; the string is 
       typically in the following format: "major.minor.patch-extra"
       where major/minor/patch are typically integers, and extra can
       be an arbitrary string and is optional
    if you ask for 4 outputs:
     [major,minor,patchnum,extra] are each field of the version string


== # Acknowledgement ==

This toolbox interacts with a number external meshing tools 
to perform the essential functionalities. These tools are listed 
below:

=== bin/tetgen ===

*Summary:tetgen is a compact and fast 3D mesh generator, it performs 
*License: ('''IMPORTANT''') tetgen is free for academic research \
and non-commertial uses only.
*URL: http://tetgen.berlios.de/
*Author: Hang Si <si at wias-berlin.de>
::Research Group: Numerical Mathematics and Scientific Computing
::Weierstrass Institute for Applied Analysis and Stochastics
::Mohrenstr. 39, 10117 Berlin, Germany

=== bin/cgalsurf ===
*Summary: cgalsurf is a utility to extract a surface mesh from \
gray-scale or binary 3D images
*Source: it is a slightly modified version from Surface_mesher \
example from CGAL 3.4
*License: CGAL v3 core library is licensed under QPL (Q Public License) \
other modules are under the Lesser General Public License (LGPL)
*URL: http://www.cgal.org/Manual/3.3/doc_html/cgal_manual/Surface_mesher/Chapter_main.html

=== bin/cgalmesh and bin/cgalpoly ===
*Summary: cgalmesh and cgalpoly are utilities to produce surface \
and volumetric meshes from a multi-valued volumetric image
*Source: it is a slightly modified version from Mesh_3 \
example from CGAL 3.5
*License: CGAL v3 core library is licensed under QPL (Q Public License) \
other modules are under the Lesser General Public License (LGPL)
*URL: http://www.cgal.org/Manual/3.5/doc_html/cgal_manual/Mesh_3/Chapter_main.html

=== bin/cgalsimp2 ===

*Summary: cgalsimp2 performs surface mesh simplification in iso2mesh.
*Source: it is adapted from Surface_mesh_simplification example of CGAL library
*License: CGAL v3 core library is licensed under QPL (Q Public License) \
other modules are under the Lesser General Public License (LGPL)
*URL: http://www.cgal.org/Manual/3.4/doc_html/cgal_manual/Surface_mesh_simplification/Chapter_main.html

=== bin/meshfix ===

*Summary: meshfix is adapted from the sample code of JMeshLib
*License: GPL (GNU General Public License) v2 or later
*URL:http://jmeshlib.sourceforge.net/
*Author:Marco Attene <attene at ge.imati.cnr.it>
::Istituto di Matematica Applicata e Tecnologie Informatiche
::Consiglio Nazionale delle Ricerche
::Via De Marini, 6 (Torre di Francia)
::16149 Genoa - ITALY 

=== bin/tetview ===

*Summary:tetview is a 3D mesh viewer
*License: see bin/tetgen
*URL: http://tetgen.berlios.de/tetview.html
*Author: see bin/tetgen


Note: iso2mesh and the above meshing utilities are considered 
as an "aggregate" rather than "derived work", based on the 
definitions in GPL FAQ (http://www.gnu.org/licenses/gpl-faq.html#MereAggregation)
Therefore, the license of iso2mesh and these utilities are independent.
The iso2mesh license only applies to the scripts and documents/data
in this package and exclude those programs stored in the bin/ directory.
The source codes of the modified meshing utilities are available
separately at iso2mesh's website and retain their upstream licenses.

Your acknowledgement of iso2mesh in your publications or 
presentations would be greatly appreciated by the author of 
this toolbox. The citation information can be found in the
Introduction section.
