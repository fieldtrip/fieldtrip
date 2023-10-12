function mesh_save_gmsh4(varargin)
% save a Gmsh mesh file
%
% USAGE:
% mesh_save_gmsh4(m, fn);
% mesh_save_gmsh4(m, fn, 'ascii');  
%
% Saves a gmsh mesh file for use with SimNIBS 2.1. Binary files are written
% as standard (otherwise, specify 'ascii' as third argument)
% 
% Note: Not all features of gmsh meshes are supported (see mesh_load_gmsh4
% for details)
%
% See http://www.geuz.org/gmsh/doc/texinfo/gmsh-full.html#SEC56 for
% documentation of the file format.
% 
% Andre Antunes 17 Oct 2014
% AT 27-Feb-2018 changed to writing of "compact" meshes


%    This program is part of the SimNIBS package.
%    Please check on www.simnibs.org how to cite our work in publications.
%
%    Copyright (C) 2013-2018 Axel Thielscher, Andre Antunes
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.



if nargin < 2
    error('at least two arguments are required: the mesh structure and the filename');
end

if nargin > 3
    disp('You should use mesh_save_gmsh4 as:');
    disp('mesh_save_gmsh4(m, fn)');
    disp('mesh_save_gmsh4(m, fn, ''type''');
    error('Incorrect number of arguments');
end

m  = varargin{1}; % mesh structure
fn = varargin{2}; % file name to write to
if ~ischar(fn)
    error('filename must be a string');
end

version=[2 0 1 8]; % write binary mesh as default
is_binary = true;
% if third argument is 'ascii", then write mesh as ASCII
if nargin == 3
    if strcmpi(varargin{3}, 'ascii')
        version(3) = 0;
        is_binary = false;
        disp('writing ascii');
    else
        error('3rd argument, if existent, must be ''ascii''');
    end
end

[~,~,extHlp] = fileparts(fn);
if ~strcmpi(extHlp,'.msh')
    fn=[fn '.msh'];
end;

fid=fopen(fn, 'W');

%% write initial header
fprintf(fid, '$MeshFormat\n');
fprintf(fid, '%d.%d %d %d\n', version);
if is_binary
    fwrite(fid, 1, 'int');
    fprintf(fid, '\n');
end
fprintf(fid, '$EndMeshFormat\n');

%% write nodes
fprintf(fid, '$Nodes\n');
nr_nodes = size(m.nodes,1);
fprintf(fid, '%d\n', nr_nodes);
if is_binary
    % data is written as:
    % int32 float64 float64 float64
    % I typecast the float64 to integers (this does not change information),
    % so that I can concatenate with the int32 column and write everything with
    % a single fwrite call. 
    % =====
    cx = typecast(m.nodes(:,1), 'uint32');
    cx = reshape(cx, [2 nr_nodes])';
    cy = typecast(m.nodes(:,2), 'uint32');
    cy = reshape(cy, [2 nr_nodes])';
    cz = typecast(m.nodes(:,3), 'uint32');
    cz = reshape(cz, [2 nr_nodes])';
    idx = uint32((1:nr_nodes)');
    fwrite(fid, [idx cx cy cz]', 'uint32');
else
    fprintf(fid, '%d %1.8g %1.8g %1.8g\n', [1:size(m.nodes,1); m.nodes']);
end
fprintf(fid, '$EndNodes\n');


%% write elements
fprintf(fid, '$Elements\n');
nr_triangles = size(m.triangles, 1);
nr_tetr = size(m.tetrahedra,1);
fprintf(fid, '%d\n', nr_triangles+nr_tetr);
if is_binary
    % elm_type, num_elm_follow, num_tags
    if nr_triangles
        fwrite(fid, [2 nr_triangles 2]', 'int');
        idx = int32((1:nr_triangles)');
        % element_number, tag1, tag2, node_number_list
        fwrite(fid, [idx m.triangle_regions m.triangle_regions m.triangles]', 'int');
    end;
    
    if nr_tetr
        fwrite(fid, [4 nr_tetr 2]', 'int');
        idx = int32(nr_triangles + (1:nr_tetr)');
        fwrite(fid, [idx m.tetrahedron_regions m.tetrahedron_regions m.tetrahedra]', 'int');
    end;
else
    % element_number, elm_type, num_tags, tag1, tag2, node_number_list
    if nr_triangles > 0
        fprintf(fid, '%d 2 2 %d %d %d %d %d\n', [1:nr_triangles; m.triangle_regions'; m.triangle_regions'; m.triangles']);
    end
    if size(m.tetrahedra,1) > 0
        fprintf(fid, '%d 4 2 %d %d %d %d %d %d\n', [nr_triangles+(1:size(m.tetrahedra,1)); m.tetrahedron_regions'; m.tetrahedron_regions'; m.tetrahedra']);
    end
end
fprintf(fid, '$EndElements\n');


%% write $NodeData
for i=1:length(m.node_data)
    fprintf(fid, '$NodeData\n');
    nr_nodes = size(m.node_data{i}.data,1);
    if  nr_nodes~=size(m.nodes,1)
        error('node data has to contain one data point per node')
    end;
    nr_comp = size(m.node_data{i}.data,2);
    
    idx = 1:nr_nodes;
    
    fprintf(fid, '1\n"%s"\n1\n0.0\n4\n0\n%d\n%d\n0\n', m.node_data{i}.name, nr_comp, nr_nodes);
    if is_binary
        % data is written as:
        % int32 nr_comp*float64
        % I typecast the float64 to integers (this does not change information),
        % so that I can concatenate with the int32 column and write everything with
        % a single fwrite call. 
        % =====
        fullt = zeros([2*nr_comp+1, nr_nodes], 'uint32');
        fullt(1,:) = idx;
        for j=1:nr_comp
            fullt(2*j:2*j+1,:) = reshape(typecast(m.node_data{i}.data(:,j), 'uint32'), 2, []);
        end
        fwrite(fid, fullt, 'uint32');
    else
        idx = cast(idx, 'double');
        fprt=repmat(' %1.8g', [1 nr_comp]);
        fprintf(fid, ['%d' fprt '\n'], [idx; m.node_data{i}.data']); % AT 24-Mar-2015
    end
    fprintf(fid, '$EndNodeData\n');
end


%% write $ElementData
for i=1:length(m.element_data)
	fprintf(fid, '$ElementData\n');
    if isfield(m.element_data{i}, 'tridata') && ...
            isfield(m.element_data{i}, 'tetdata')
        if  ~isempty(m.element_data{i}.tridata)&&...
            size(m.element_data{i}.tridata,1)~= size(m.triangles,1)
            error('element tridata has to contain one data point per triangle')
        end
        if  ~isempty(m.element_data{i}.tetdata)&&...
            size(m.element_data{i}.tetdata,1)~= size(m.tetrahedra,1)
            error('element tetdata has to contain one data point per tetrahedron')
        end
    
        if  ~isempty(m.element_data{i}.tetdata)&&...
            ~isempty(m.element_data{i}.tridata)&&...    
            size(m.element_data{i}.tetdata,2) ~= size(m.element_data{i}.tridata,2)
            error('element tridata and tetdata has to have the same number of rows')
        end
        data=[m.element_data{i}.tridata; m.element_data{i}.tetdata];
    elseif isfield(m.element_data{1}, 'data')
        if size(m.element_data{i}.data,1)~= size(m.tetrahedra,1) + size(m.triangles,1)
            error('element data has to contain one data point per triangle and tetrahedra')
        end
        data = m.element_data{i}.data;
    else
        error('element data needs a tridata/tetdata pair or a data field')
    end
    nr_elements = size(data,1);
    nr_comp = size(data,2);
    
    idx = 1:nr_elements;
    if isfield(m.element_data{i}, 'tridata') && isempty(m.element_data{i}.tridata)
        idx=idx+size(m.triangles,1);
    end
         
    fprintf(fid, '1\n"%s"\n1\n0.0\n4\n0\n%d\n%d\n0\n', m.element_data{i}.name, nr_comp, nr_elements);
    if is_binary
        % data is written as:
        % int32 nr_comp*float64
        % I typecast the float64 to integers (this does not change information),
        % so that I can concatenate with the int32 column and write everything with
        % a single fwrite call. 
        % =====
        fullt = zeros([2*nr_comp+1, nr_elements], 'uint32');
        fullt(1,:) = idx;
        for j=1:nr_comp
            fullt(2*j:2*j+1,:) = reshape(typecast(data(:,j), 'uint32'), 2, []);
        end
        fwrite(fid, fullt, 'uint32');
    else
        fprt=repmat(' %1.8g', [1 nr_comp]);
        fprintf(fid, ['%d' fprt '\n'], [idx; data']);
    end
    fprintf(fid, '$EndElementData\n');
end


if isfield(m,'element_node_data')&&~isempty(m.element_node_data)
    warning('I am not saving $ElementNodeData information');
end


fclose(fid);