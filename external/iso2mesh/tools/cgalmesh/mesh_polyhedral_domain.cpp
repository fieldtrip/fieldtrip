#include <CGAL/AABB_intersections.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

#include <CGAL/Random.h>

// IO
#include <CGAL/IO/Polyhedron_iostream.h>

// Domain
struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria Facet_criteria;
typedef Mesh_criteria::Cell_criteria Cell_criteria;


void usage(char *exename){
	printf("usage:\n\t%s input.off output.mesh <angle|30> <surf-size|6> <approx|4> <rad-edge-ratio|3> <tetra-size|8> <randseed>\n",exename);
	exit(1);
}
int main(int argc,char *argv[])
{
  // Create polyhedron
  Polyhedron polyhedron;

  float angle=30.f,ssize=6.f,approx=4.f,reratio=3.f,vsize=8.f;

  printf("Volumetric Mesh Generation Utility (Based on CGAL %s)\n\
(modified for iso2mesh by Qianqian Fang)\nhttp://iso2mesh.sf.net\n\n",CGAL_VERSION_STR);

  if(argc!=3&&argc!=8&&argc!=9){
  	usage(argv[0]);
  }
  if(argc>=8){
  	angle=atof(argv[3]);
  	ssize=atof(argv[4]);
  	approx=atof(argv[5]);
  	reratio=atof(argv[6]);
  	vsize=atof(argv[7]);
  }
  if(argc==9 && atoi(argv[8])>0){
        printf("RNG seed=%d\n",atoi(argv[8]));
        CGAL::Random rd(atoi(argv[8]));
        CGAL::Random::State st;
        rd.save_state(st);
        CGAL::default_random.restore_state(st);
  }
  std::ifstream input(argv[1]);
  input >> polyhedron;

  // Create domain
  Mesh_domain domain(polyhedron);

  // Mesh criteria
  Facet_criteria facet_criteria(angle, ssize, approx); // angle, size, approximation
  Cell_criteria cell_criteria(reratio, vsize); // radius-edge ratio, size
  Mesh_criteria criteria(facet_criteria, cell_criteria);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

  // Output
  std::ofstream medit_file;
  if(argc>=3)
	medit_file.open(argv[2]);
  else
	medit_file.open("output.mesh");
  c3t3.output_to_medit(medit_file);

  return 0;
}
