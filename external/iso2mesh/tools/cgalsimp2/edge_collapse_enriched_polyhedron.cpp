//-------------------------------------------------------------------------
//===================================================================
//=   CGAL Mesh Simplification Code Modified for iso2mesh toolbox   =
//===================================================================
//
//Modified by Qianqian Fang <fangq at nmr.mgh.harvard.edu>
//patches from 
//  Fernando Cacciola <fernando.cacciola at geometryfactory.com>
//  Andreas Fabri <andreas.fabri at geometryfactory.com>
//
//Date: 2008/03
//
//-------------------------------------------------------------------------

#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

// Adaptor for Polyhedron_3
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>// AF: added
//#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h> // AF: removed

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Extended polyhedron items which include an id() field
#include <CGAL/Polyhedron_items_with_id_3.h>

// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>

// Default cost and placement policies
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h> 

// Non-default cost and placement policies
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_and_length.h> 

typedef CGAL::Simple_cartesian<double> Kernel;

typedef Kernel::Point_3 Point ;

//
// Setup an enriched polyhedron type which stores an id() field in the items
//
typedef CGAL::Polyhedron_3<Kernel,CGAL::Polyhedron_items_with_id_3> Surface; 

typedef Surface::Halfedge_handle Halfedge_handle ;

namespace SMS = CGAL::Surface_mesh_simplification ;

typedef SMS::Edge_profile<Surface> Profile ;


// The following is a Visitor that keeps track of the simplification process.
// In this example the progress is printed real-time and a few statistics are
// recorded (and printed in the end).
//
struct Visitor
{
  Visitor() 
    : collected(0)
    , processed(0)
    , collapsed(0)
    , non_collapsable(0)
    , cost_uncomputable(0) 
    , placement_uncomputable(0) 
  {} 

  // Called on algorithm entry  
  void OnStarted( Surface& ) {} 
  
  // Called on algorithm exit  
  void OnFinished ( Surface& ) { std::cerr << "\n" << std::flush ; } 
  
  // Called when the stop condition returned true
  void OnStopConditionReached( Profile const& ) {} 
  
  // Called during the collecting phase for each edge collected.
  void OnCollected( Profile const&, boost::optional<double> const& )
  {
    ++ collected ;
//    std::cerr << "\rEdges collected: " << collected << std::flush ;
  }                
  
  // Called during the processing phase for each edge selected.
  // If cost is absent the edge won't be collapsed.
  void OnSelected(Profile const&          
                 ,boost::optional<double> cost
                 ,std::size_t             initial
                 ,std::size_t             current
                 )
  {
    ++ processed ;
    if ( !cost )
      ++ cost_uncomputable ;
 /*     
    if ( current == initial )
      std::cerr << "\n" << std::flush ;
    std::cerr << "\r" << current << std::flush ;
 */
  }                
  
  // Called during the processing phase for each edge being collapsed.
  // If placement is absent the edge is left uncollapsed.
  void OnCollapsing(Profile const&          
                   ,boost::optional<Point>  placement
                   )
  {
    if ( placement )
         ++ collapsed;
    else ++ placement_uncomputable ;
  }                
  
  // Called for each edge which failed the so called link-condition,
  // that is, which cannot be collapsed because doing so would
  // turn the surface into a non-manifold.
  void OnNonCollapsable( Profile const& )
  {
    ++ non_collapsable;
  }      

  // AF: added 
  void OnCollapsed(Profile const&, Surface::Vertex_handle)
  {}
  
  std::size_t  collected
             , processed
             , collapsed
             , non_collapsable
             , cost_uncomputable  
             , placement_uncomputable ; 
} ;

template<class GetCost_>
struct Cost_with_fixed_edges : GetCost_
{
 typedef GetCost_ GetCost ;

  // typedef typename GetCost::Profile     Profile ; // AF Profile no longer nested type
  //typedef typename GetCost::Point       Point ; // AF
  //typedef typename GetCost::result_type result_type ; // AF

  template <typename Profile, typename T> // AF: make it a template 
  boost::optional<typename Profile::FT>  
operator()( 
Profile const& aProfile, 
T const& aPlacement ) const // replace result_type
 {
    if ( aProfile.border_edges().size() > 0 )
          return boost::none ;
    else
         return this->GetCost::operator()(aProfile, aPlacement);
 }
} ;

typedef Cost_with_fixed_edges< SMS::Edge_length_cost<Surface> > My_cost ; 
typedef Cost_with_fixed_edges< SMS::LindstromTurk_cost<Surface> > LT_cost ; 


int main( int argc, char** argv ) 
{
  Surface surface; 
  float maxface=0.1;
  int defaultpolicy=0;

  printf("= Surface Mesh Simplification Utility (Based on CGAL %s) =\n(modified for iso2mesh by Qianqian Fang, Fernando Cacciola and Andreas Fabri)\nhttp://iso2mesh.sf.net\n\n",CGAL_VERSION_STR);
  if(argc<2){
  	printf("command options\n\t%s input.off keepratio\nFor example:\n\t%s input.off 0.1\n",argv[0],argv[0]);
	exit(0);
  }
  std::ifstream is(argv[1]) ; is >> surface ;
  if(argc>2) maxface=atof(argv[2]);
  if(maxface<0){
      defaultpolicy=1;
      maxface=-maxface;
  }

  printf("max face ratio=%f\n",maxface);
   
   
  // The items in this polyhedron have an "id()" field 
  // which the default index maps used in the algorithm
  // need to get the index of a vertex/edge.
  // However, the Polyhedron_3 class doesn't assign any value to
  // this id(), so we must do it here:
  int index = 0 ;
  
  for( Surface::Halfedge_iterator eb = surface.halfedges_begin()
     , ee = surface.halfedges_end()
     ; eb != ee
     ; ++ eb
     ) 
    eb->id() = index++;

  printf("edge index number %d\n",index);


  index = 0 ;
  for( Surface::Vertex_iterator vb = surface.vertices_begin()
     , ve = surface.vertices_end()
     ; vb != ve
     ; ++ vb
     ) 
    vb->id() = index++;
    
  // In this example, the simplification stops when the number of undirected edges
  // drops below 10% of the initial count
  SMS::Count_ratio_stop_predicate<Surface> stop(maxface);
     
  Visitor vis ;
  
  // The index maps are not explicitelty passed as in the previous
  // example because the surface items have a proper id() field.
  // On the other hand, we pass here explicit cost and placement
  // function which differ from the default policies, ommited in
  // the previous example.
  
  printf("mesh simplificaton in progress ...\n");

  int r;
  if(defaultpolicy)
      r = SMS::edge_collapse
           (surface
           ,stop
           ,CGAL::parameters::get_cost     (LT_cost())
                 .get_placement   (SMS::LindstromTurk_placement<Surface>())
            .visitor(vis)  // AF: not a pointer
           );
  else
      r = SMS::edge_collapse
           (surface
           ,stop
           ,CGAL::parameters::get_cost     (My_cost())
                 .get_placement   (SMS::Midpoint_placement<Surface>())
                 .visitor(vis)
           );  
  
  std::cout << "\nEdges collected: " << vis.collected
            << "\nEdges proccessed: " << vis.processed
            << "\nEdges collapsed: " << vis.collapsed
            << std::endl
            << "\nEdges not collapsed due to topological constrians: " 
            << vis.non_collapsable
            << "\nEdge not collapsed due to cost computation constrians: " 
            << vis.cost_uncomputable 
            << "\nEdge not collapsed due to placement computation constrians: " 
            << vis.placement_uncomputable 
            << std::endl ; 
            
  std::cout << "\nFinished...\n" << r << " edges removed.\n" 
            << (surface.size_of_halfedges()/2) << " final edges.\n" ;
        
  std::ofstream os( argc > 3 ? argv[3] : "out.off" ) ; os << surface ;
  
  return 0 ;      
}



// EOF //

