function [m]=mesh_empty(varargin)
  % returns an empty mesh
  % USAGE: m=mesh_empty;

  m.nodes=[];
  m.triangles=[];
  m.triangle_regions=[];
  m.tetrahedra=[];
  m.tetrahedron_regions=[];
  m.node_data={};
  m.element_data={};
