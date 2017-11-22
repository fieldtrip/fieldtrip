function savemsh(node,elem,fname,rname)
%
% savemsh(node,elem,fname,rname)
%
% save a tetrahedral mesh to GMSH mesh format
%
% author: Riccardo Scorretti (riccardo.scorretti<at> univ-lyon1.fr)
% date: 2013/07/22
%
% input:
%      node: input, node list, dimension (nn,3)
%      elem: input, tetrahedral mesh element list, dimension (ne,4) or (ne,5) for multi-region meshes
%      fname: output file name
%      rname: name of the regions, cell-array of strings (optional)
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if nargin < 4 , rname = {} ; end
if size(elem,2) < 5 , elem(:,5) = 1 ; end
fid = fopen(fname,'wt');
if(fid==-1)
    error('You do not have permission to save mesh files.');
end
nbNodes = size (node,1);
reg = unique (elem(:,5));
nbRegion = length (reg);
nbElements = size (elem, 1);

% Create the skeleton of the mesh structure
M.Info.version = [];

M.Nodes.nb = 0;
M.Nodes.x = [];
M.Nodes.y = [];
M.Nodes.z = [];

M.Elements.nb = 0;
M.Elements.type = zeros (0, 0, 'uint8');
M.Elements.tableOfNodes = zeros (0, 0, 'uint32');
M.Elements.region = zeros (0, 0, 'uint16');

M.Regions.nb = 0;
M.Regions.name = {};
M.Regions.dimension = [];

% Build the table of nodes

M.Nodes.nb = nbNodes;
M.Nodes.x = node(:,1);
M.Nodes.y = node(:,2);
M.Nodes.z = node(:,3);
clear node

% Build the table of elements

M.Elements.nb = nbElements;
M.Elements.type = uint8(4*ones(nbElements, 1));
M.Elements.tableOfNodes = uint32(elem(:,1:4));
M.Elements.region = uint16(elem(:,5));
clear elem

% Build the table of regions

M.Regions.nb = max(reg);
for k = 1 : nbRegion
   if length(rname) < k , rname{k} = sprintf('region_%d', k) ; end
   M.Regions.name{reg(k)} = sprintf ('%s', rname{k});
   M.Regions.dimension(reg(k)) = 3;
end

% Writhe the header
fprintf (fid, '$MeshFormat\n2.2 0 8\n$EndMeshFormat\n');

% Write the physical names
if M.Regions.nb > 0 
   fprintf (fid, '$PhysicalNames\n');
   fprintf (fid, '%d\n', M.Regions.nb);
   for r = 1 : M.Regions.nb
      name = M.Regions.name{r};
      if isempty (name)
         name = sprintf ('Region_%d', r);
      end
      fprintf (fid, '%d %d "%s"\n', M.Regions.dimension(r), r, name);
   end
   fprintf (fid, '$EndPhysicalNames\n');
end

% Write the nodes
fprintf (fid, '$Nodes\n');
fprintf (fid, '%d\n', size(M.Nodes.x,1));
buffer = [ 1:M.Nodes.nb ; M.Nodes.x' ; M.Nodes.y' ; M.Nodes.z' ];
fprintf (fid, '%d %10.10f %10.10f %10.10f\n', buffer);
fprintf (fid, '$EndNodes\n');

% Write the elements
%
% In order to accelerate the printing, the elements are printed in groups of (blockSize) elements, and
% are grouped by (homogeneous) type. This variable sets the size of each group.
%
blockSize = 100000;

fprintf (fid, '$Elements\n');
fprintf (fid, '%d\n', M.Elements.nb);
for h = 1 : blockSize : ceil(length(M.Elements.type)/blockSize)*blockSize
   e = h : min(length(M.Elements.type) , h+blockSize-1);  % = elements being considered
   type = unique (M.Elements.type(e));                    % = types of elements found in this group
   
   %
   % Process each type of element separately
   %
   for k = 1 : length(type)
      if type(k) == 0, continue; end
      et = e(find(M.Elements.type(e) == type(k)));  % = elements of the group of the same type
      
      %
      % Determine the format for printing the elements
      %
      elementFormat = '%d %d %d %d %d %d\n';
      for n = 1 : 4
         elementFormat = [ elementFormat '%d ' ];
      end
      elementFormat = [ elementFormat '\n' ];
      
      %
      % Collect in a buffer all the data of the elements of index (et)
      %
      buffer = zeros (10 , length(et));
      buffer(1,:) = et;
      buffer(2,:) = type(k);
      buffer(3,:) = 3;
      buffer(4,:) = M.Elements.region(et);
      buffer(5,:) = M.Elements.region(et);
      buffer(6,:) = 0;
      for n = 1 : 4
         buffer(6+n,:) = M.Elements.tableOfNodes(et,n);
      end
      
      %
      % Print all the homogeneous elements in the group with a single instruction
      %
      fprintf (fid, elementFormat, buffer);
   end
end
fprintf (fid, '$EndElements\n');


fclose(fid);
