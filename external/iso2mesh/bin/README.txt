= External Utilities Used by Iso2Mesh Toolbox =

In this directory, you can find the precompiled binaries
that are used by iso2mesh (http://iso2mesh.sf.net). The 
names and origins of these tools are explained below:


# cgalsurf.*  - Surface Mesh Extraction Utility (built on CGAL)
# cgalmesh.*  - Volume/Surface Mesh Generation Utility (built on CGAL)
# cgalpoly.*  - Volumetric Mesh Generation Utility (built on CGAL)
# cgalsimp2.* - Surface Mesh Simplification Utility (built on CGAL)
# tetgen.*    - A Quality Tetrahedral Mesh Generator and 3D Delaunay Triangulator by Hang Si
# tetview.*   - A Mesh Graphing Utility by Hang Si
# jmeshlib.*  - Mesh Validation and Repairing Utility (build on JMeshLib) by Marco Attene
# meshfix.*   - Mesh Validation and Repairing Utility (build on Meshfix) by Marco Attene
# cork.*      - A robust surface mesh Boolean operation algorithm by Gilbert Bernstein (https://github.com/gilbo/cork)
# PoissonRecon.* - Screened Poisson Surface Reconstruction (Version 8.0) by Michael Kazhdan (http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version8.0/)

To distinguish binaries for different platforms, we added 
matlab mexext output for different platforms as the file
extension (except for Windows). For example, the binary
for cgalsurf on GNU Linux i386 platform is cgalsurf.mexglx,
and that for Windows is cgalsurf.exe.

Please be aware that iso2mesh communicates with these external
tools via pipes and disk files. Thus, iso2mesh and these tools
are considered as "aggregate" and they can be distributed independently
under different licenses.

Nonetheless, all of the above tools are free-software 
(free as in free-beer and freedom). They are distributed under 
one of the OSI-approved open-source licenses. 

For other utilities, cgalsurf, cgalmesh, cgalpoly and 
cgalsimp2 were modified from CGAL 3.x and the binary are
licensed under QPL (The Q Public License v1.0 
http://www.opensource.org/licenses/qtpl). The utility
meshfix was compiled from JMeshLib and is licensed under
GPL v2 or later. The cork binary was compiled using
a modified cork source code, released under the LGPL 
license. The modified source codes of these utilities
are provided on the iso2mesh project management website
(gforge server) and subversion source code repository.
It is COMPLETELY the user's responsibility to use all 
external utilities within the permissions outlined by 
the upstream licenses.

More detailed information regarding these tools can be found
in the last section of the README.txt file under the iso2mesh
main folder.
