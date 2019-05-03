#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <fstream>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Gray_level_image_3.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/Random.h>

// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;

// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

typedef Tr::Geom_traits GT;
typedef CGAL::Gray_level_image_3<GT::FT, GT::Point_3> Gray_level_image;
typedef CGAL::Implicit_surface_3<GT, Gray_level_image> Surface_3;

void usage(char *exename){
	printf("command options:\n\t%s inrfile thres x0 y0 z0 r2 ang br bd maxnode outputoff <randseed> <initpointnum>\nFor example:\n\
\t%s data/skull_2.9.inr 2.9 122. 102. 117. 80000 30 5. 5. 100000 out.off 123456789 100\n",
           exename, exename);
	exit(1);
}

int main(int argc, char** argv) {
  unsigned int randseed, initnum=100;

  // the 'function' is a 3D gray level image
  printf("Surface Mesh Extraction Utility (Based on CGAL %s)\n(modified for iso2mesh by Qianqian Fang)\n\
http://iso2mesh.sf.net\n\n",CGAL_VERSION_STR);

  if(argc!=12&&argc!=13 && argc!=14){
	usage(argv[0]);
  }
  if(argc>=13){
        sscanf(argv[12],"%d",&randseed);
        if(randseed>0){
		printf("RNG seed %d\n",atoi(argv[12]));
		CGAL::Random rd(atoi(argv[12]));
		CGAL::Random::State st;
		rd.save_state(st);
		CGAL::default_random.restore_state(st);
	}
	if(argc>=14){
		sscanf(argv[13],"%d",&initnum);
		printf("set initial mesh size to %d\n",initnum);
	}
  }
  Tr tr;            // 3D-Delaunay triangulation
  C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation
  
  Gray_level_image image(argv[1], atof(argv[2]));
  // Carefully choosen bounding sphere: the center must be inside the
  // surface defined by 'image' and the radius must be high enough so that
  // the sphere actually bounds the whole image.
  GT::Point_3 bounding_sphere_center(atof(argv[3]), atof(argv[4]), atof(argv[5]));

  GT::FT bounding_sphere_squared_radius = atof(argv[6]);

  GT::Sphere_3 bounding_sphere(bounding_sphere_center,
                                   bounding_sphere_squared_radius);

  // definition of the surface, with 10^-5 as relative precision
  Surface_3 surface(image, bounding_sphere, 1e-5);

  // defining meshing criteria
  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(atof(argv[7]),
                                                     atof(argv[8]),
                                                     atof(argv[9]));

  // meshing surface, with the "manifold without boundary" algorithm

  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag(),initnum);

  std::ofstream out(argv[11]);
  CGAL::output_surface_facets_to_off (out, c2t3);
  std::cout << "Final number of points: " << tr.number_of_vertices() << "\n";
}

