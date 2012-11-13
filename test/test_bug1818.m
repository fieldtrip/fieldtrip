function test_bug1818

% to test the reading function of meshes used for FEM headmodel

% TEST test_bug1818
% TEST ft_read_headshape ft_datatype_parcellation

% vista 
 % mesh = ft_read_headshape('~/test_cases/fem/cube2mm3layervorwerk_ns_127_127_127.v');  % error
 
 % path needs to be changed
 mesh = ft_read_headshape('/home/common/matlab/fieldtrip/data/test/bug1818/cube2mm3layervorwerk_ns_127_127_127.v');

 parcellation1 = ft_datatype_parcellation(mesh); 
 parcellation2 = ft_datatype_parcellation(mesh, 'parcellationstyle', 'probabilistic');
 
 if ~(ft_datatype(parcellation1, 'parcellation'))
     error 'vista headshape is not parcellation datatype';
 end
 
 if ~(ft_datatype(parcellation2, 'parcellation'))
     error 'vista headshape is not parcellation datatype';
 end


% tetgen - elements
mesh = ft_read_headshape('/home/common/matlab/fieldtrip/data/test/bug1818/tet_4layer_127_127_127.1.ele');   

parcellation1 = ft_datatype_parcellation(mesh); 
parcellation2 = ft_datatype_parcellation(mesh, 'parcellationstyle', 'probabilistic');
 
 if ~(ft_datatype(parcellation1, 'parcellation'))
     error 'vista headshape is not parcellation datatype';
 end
 
 if ~(ft_datatype(parcellation2, 'parcellation'))
     error 'vista headshape is not parcellation datatype';
 end
 
% tetgen - node 
mesh = ft_read_headshape('/home/common/matlab/fieldtrip/data/test/bug1818/tet_4layer_127_127_127.1.node');   

parcellation1 = ft_datatype_parcellation(mesh); 
parcellation2 = ft_datatype_parcellation(mesh, 'parcellationstyle', 'probabilistic');
 
 if ~(ft_datatype(parcellation1, 'parcellation'))
     error 'vista headshape is not parcellation datatype';
 end
 
 if ~(ft_datatype(parcellation2, 'parcellation'))
     error 'vista headshape is not parcellation datatype';
 end
 
 
 

