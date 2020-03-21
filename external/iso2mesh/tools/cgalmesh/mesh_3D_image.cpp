#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Labeled_image_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>

// start PV //
#include <string>
#include <vector>
#include <CGAL/Mesh_constant_domain_field_3.h>
// end PV //

// Domain
struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Image_3 Image;
typedef CGAL::Labeled_image_mesh_domain_3<Image,K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria    Facet_criteria;
typedef Mesh_criteria::Cell_criteria     Cell_criteria;

// start PV //
typedef CGAL::Mesh_constant_domain_field_3<Mesh_domain::R,
                                           Mesh_domain::Index> Sizing_field_cell;
// end PV //


void usage(char *exename){
	printf("usage:\n\t%s input.inr output.mesh <angle|30> <surf-size|6> <approx|4> <rad-edge-ratio|3> \
<tetra-size|8> <randomseed|-1>\n\
example:\n\t%s input.inr output.mesh 30 6 4 3 8 123456789\n",exename,exename);
	exit(1);
}
int main(int argc,char *argv[])
{
  // Loads image

// start PV //
  //float angle=30.f,ssize=6.f,approx=4.f,reratio=3.f,vsize=8.f;
  float angle=30.f,ssize=6.f,approx=4.f,reratio=3.f,maxvol=0.f;
// end PV //
  int labelid=0, lid;

  printf("Volume/Surface Mesh Generation Utility (Based on CGAL %s)\n\
(modified for iso2mesh by Qianqian Fang and Peter Varga)\nhttp://iso2mesh.sf.net\n\n",CGAL_VERSION_STR);

  if(argc!=3&&argc!=8&&argc!=9){
  	usage(argv[0]);
  }
  if(argc>=8){
  	angle=atof(argv[3]);
  	ssize=atof(argv[4]);
  	approx=atof(argv[5]);
  	reratio=atof(argv[6]);
// start PV //
  	//vsize=atof(argv[7]);
// end PV //
  }
  if(argc==9 && atoi(argv[8])>0){
	printf("RNG seed=%d\n",atoi(argv[8]));
        CGAL::Random rd(atoi(argv[8]));
        CGAL::Random::State st;
        rd.save_state(st);
        CGAL::default_random.restore_state(st);
  }
  Image image;
  image.read(argv[1]);
  Mesh_domain domain(image);

  // Mesh criteria
  Facet_criteria facet_criteria(angle, ssize, approx); // angle, size, approximation

// start PV //
  //Cell_criteria cell_criteria(reratio, vsize); // radius-edge ratio, size

  std::string vsize;
  vsize=argv[7];

  std::vector<double> vsize_vect;
  std::vector<int> labels_vect;

  std::stringstream stream_vsize(vsize);
  std::string word_vsize;

  while( std::getline(stream_vsize, word_vsize, ':') ){
    int len=sscanf(word_vsize.c_str(), "%d=%f", &lid, &maxvol);
    if(len==2){
      labelid=lid;
      if(maxvol<=0.f){
         std::cerr << "cell volume must be positive" << std::endl;
         exit(-2);
      }
      labels_vect.push_back(labelid);
      vsize_vect.push_back(maxvol);
    }else{
      len=sscanf(word_vsize.c_str(), "%f", &maxvol);
      if(len!=1 || maxvol<0 || word_vsize.find_first_of('=')!=std::string::npos){
         std::cerr << "invalid sizing field label '" << word_vsize << "', please check your command" << std::endl;
         exit(-1);
      }
      labels_vect.push_back(++labelid);
      vsize_vect.push_back(maxvol);
    }
  }

  float max_vsize = *max_element(vsize_vect.begin(),vsize_vect.end());

  Sizing_field_cell vsize_cell(max_vsize);
  int volume_dimension = 3;

  for (int i=0; i<vsize_vect.size(); i++)
    vsize_cell.set_size(vsize_vect[i], volume_dimension, 
                domain.index_from_subdomain_index(labels_vect[i]));

  std::cout << "Mesh sizes are (label=size) ";
  for (int i=0; i<vsize_vect.size(); i++){
    std::cout << "(" << labels_vect[i] <<  "=" << vsize_vect[i] << ") ";
  }
  std::cout << std::endl;

  Cell_criteria cell_criteria(reratio, vsize_cell);
// end PV //

  Mesh_criteria criteria(facet_criteria, cell_criteria);

  // Meshing
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

  // Output
  std::ofstream medit_file;
  if(argc>=8)
	medit_file.open(argv[2]);
  else
	medit_file.open("output.mesh");
  c3t3.output_to_medit(medit_file);
  medit_file.close();

  return 0;
}
