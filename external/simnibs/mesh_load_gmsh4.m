function m = mesh_load_gmsh4( varargin )
% load a Gmsh mesh file (binary or ASCII)
%
% USAGE:
% m=mesh_load_gmsh4(fn);
%
% fn: filename
% m:  mesh-structure with
%       nodes
%       triangles: indices into the nodes, starting at 1
%       triangle_regions: region numbers (gmsh "physical groups")
%       tetrahedra
%       tetrahedron_regions
%
%       node_data{idx}.name: string
%       node_data{idx}.data: vector or matrix
%       
%       element_data{idx}.name
%       element_data{idx}.tridata: vector or matrix (triangle data)
%       element_data{idx}.tetdata: vector or matrix (tetrahedra data)
%
%       element_node_data{idx}.name
%       element_node_data{idx}.data: vector or matrix (only tetrahedra data)
%
%
% Note: This function was written for the meshes created by SimNIBS 2.1 or
% later and does not support all features of Gmsh meshes:
% * only triangles and tetrahedra are supported element types
% * other elements will be skipped; if this occurs, reading of element-data
%   will skipped as well
% * loading binary meshes directly produced by gmsh and SimNIBS versions
%   prior to 2.1 will be very slow
% * node numbers have to be a continuous list of indexes starting at 1
% * triangles are stored prior to tetrahedra
% * element_node_data works only for binary meshes and assumes that 
%   only tetrahedral data is stored
%   
% See http://www.geuz.org/gmsh/doc/texinfo/gmsh-full.html#SEC56 for
% documentation of the file format.
% 
% Andre Antunes 07 Apr 2015
% AT 27-Feb-2018: simplified code, added checks for mesh format

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

%%
tic

if nargin<1
   [fname, pname]= uigetfile('*.msh');
   if isequal(fname,0) || isequal(pname,0); return; end;
   fn = [pname fname];
else
    fn = varargin{1};
end

if ~ischar(fn)
    error(['invalid filename:' fn]);
end

if ~exist(fn, 'file')
    error(['file does not exist:' fn]);
end

m.node_data = {};
m.element_data = {};

fid=fopen(fn);

tline = fgetl(fid);
if ~strcmp('$MeshFormat', tline)
    error('Start tag $MeshFormat expected');
end

% parse 2nd line: version.minor_version is_binary byte_size
tline = fgetl(fid);
version = sscanf(tline, '%f %d %d');
v = int32(version(1));
if v ~= 2 && v ~= 4
    error(['Cant read mesh version ' num2str(version(1))]);
end
if version(3) ~= 8
    error(['expected to read 8 byte precision, but encountered: ' num2str(version(4))]);
end

is_binary = version(2);
if is_binary
    disp('reading binary')
    endianess = fread(fid, 1, 'int');
    if (endianess ~= 1)
        error('Endianess in binary file is different than 1. I cant procceed..\n Was this file created in windows?');
    end
    a = fread(fid, 1, 'char'); % read end of line ( should be ascii 010)
    if a ~= 10; error(['Expected LF (new line), but I read ASCII:' a]); end
end

tline = fgetl(fid);
if ~strcmp('$EndMeshFormat', tline)
    error('End tag $EndMeshFormat expected');
end


%% read nodes
line = fgetl(fid);
while ~strcmp('$Nodes', line)
    line = fgetl(fid);
end

if is_binary
   if v == 2
       m = read_nodes_binary(m, fid);
   else
       m = read_nodes_binary4(m, fid);
   end
else
   if v == 2
       m = read_nodes(m, fid);
   else
       m = read_nodes4(m, fid);
   end
end

%% read elements
line = fgetl(fid);
while ~strcmp('$Elements', line)
    line = fgetl(fid);
end

if is_binary
    if v == 2
        [m, continous_elm_numbers] = read_elements_binary(m, fid);
    else
        [m, continous_elm_numbers] = read_elements_binary4(m, fid);
    end
else
    if v == 2
        [m, continous_elm_numbers] = read_elements(m, fid);
    else
        [m, continous_elm_numbers] = read_elements4(m, fid);
    end
end


%% read data ($NodeData, $ElementData, $ElementNodeData)
tline = fgetl(fid);
while ~feof(fid)
    if strcmp('$NodeData', tline)
        m.node_data{end+1} = read_node_data(fid, size(m.nodes,1), is_binary);
    elseif strcmp('$ElementData', tline)
        m.element_data{end+1} = read_element_data(fid, size(m.triangles,1), size(m.tetrahedra,1), continous_elm_numbers, is_binary);
    elseif strcmp('$ElementNodeData', tline)
        if is_binary
            m.element_node_data = {};
            [data, name, ~] = read_data_binary(fid, tline);
            m.element_node_data{end+1}.data = data;
            m.element_node_data{end}.name = name;
            warning('reading of element_node_data is not tested. Good luck.');
        else
            error('$ElementNodeData is not supported for ASCII meshes.');
        end;
    else
        error(['Unsupported field:' tline]);
    end
    tline = fgetl(fid);
end
    
m.node_data = m.node_data';
m.element_data = m.element_data';
if isfield(m,'element_node_data')
    m.element_node_data = m.element_node_data';
end;

fprintf('Number of Nodes              : %d\n', size(m.nodes,1));
fprintf('Number of Triangles          : %d\n', size(m.triangles,1));
fprintf('Number of Triangle Regions   : %d\n', size(unique(m.triangle_regions),1));
fprintf('Number of Tetrahedra         : %d\n', size(m.tetrahedra,1));
fprintf('Number of Tetrahedron Regions: %d\n', size(unique(m.tetrahedron_regions),1));
 
toc

function m = read_nodes_binary(m, fid)
    % get number of nodes (this line is in ascii)
    tline = fgetl(fid);
    number_of_nodes = sscanf(tline,'%d');
    if ~isnumeric( number_of_nodes )
        error('number of nodes is not a number');
    end
   
    % get number of first node
    node_number(1) = fread(fid, 1 , 'int');
    
    % read 3 double columns, skipping 4 bytes after each read
    m.nodes = fread(fid, [3 number_of_nodes], '3*double',4)';
    
    % get number of last node
    fseek(fid, -(2*4+3*8), 'cof');
    node_number(2) = fread(fid, 1 , 'int', 3*8);   
 
    if (node_number(1) < 1) || ...
       (node_number(2) > number_of_nodes) || ...
       (size(m.nodes,1) ~= number_of_nodes)
        error('node numbers have to be a continuous list of indexes starting at 1');
    end
    
    % sometimes there is an end of line after $Nodes, sometimes there is
    % not... handle that
    read_LF(fid);
   
    % confirm last line
    tline = fgetl(fid);
    if ~strcmp('$EndNodes', tline)
        error('End tag $EndNodes expected'); 
    end


    
function m = read_nodes(m, fid)
    % get number of nodes
    tline = fgetl(fid);
    number_of_nodes = sscanf(tline,'%d');
    if ~isnumeric( number_of_nodes )
        error(['number of nodes is not a number']);
    end
    
    pts = textscan(fid, '%d %f %f %f');
    node_numbers=pts{1};
    m.nodes = [pts{2} pts{3} pts{4}];

    if (min(node_numbers) < 1) || ...
       (max(node_numbers)> number_of_nodes) || ...
       (length(node_numbers) ~= number_of_nodes)
        error('node numbers have to be a continuous list of indexes starting at 1');
    end

    % read end line
    tline = fgetl(fid);
    if ~strcmp('$EndNodes', tline)
        error('End tag $EndNodes expected'); 
    end
    
function m = read_nodes_binary4(m, fid)
    n_blocks = fread(fid,1, 'uint64');
    number_of_nodes = fread(fid,1, 'uint64');
    if ~isnumeric(number_of_nodes)
        error('number of nodes is not a number');
    end
    node_numbers = [];%zeros(number_of_nodes, 0, 'int32');
    m.nodes = [];%zeros(number_of_nodes, 3);
    for b = 1:n_blocks
        fread(fid,2, 'int');
        if fread(fid,1, 'int')
            error('Cant read parametric nodes')
        end
        n_in_block = fread(fid,1, 'uint64');
        node_numbers = [node_numbers; fread(fid, n_in_block, 'int', 3*8)];
        fseek(fid, -n_in_block*(4+3*8)+4, 'cof');
        m.nodes = [m.nodes; fread(fid, [3, n_in_block], '3*double', 4)'];
        fseek(fid, -4, 'cof');
    end
    [node_numbers, order] = sort(node_numbers);
    m.nodes = m.nodes(order, :);
    if (min(node_numbers) < 1) || ...
       (max(node_numbers)> number_of_nodes) || ...
       (length(node_numbers) ~= number_of_nodes)
        error('node numbers have to be a continuous list of indexes starting at 1');
    end
    % sometimes there is an end of line after $Nodes, sometimes there is
    % not... handle that
    read_LF(fid);
   
    % confirm last line
    tline = fgetl(fid);
    if ~strcmp('$EndNodes', tline)
        error('End tag $EndNodes expected'); 
    end
    
function m = read_nodes4(m, fid)
    % get number of nodes
    tline = fgetl(fid);
    l = sscanf(tline,'%d %d');
    n_blocks = l(1);
    number_of_nodes = l(2);
    if ~isnumeric(number_of_nodes)
        error('number of nodes is not a number');
    end
    node_numbers = [];%zeros(number_of_nodes, 0, 'int32');
    m.nodes = [];%zeros(number_of_nodes, 3);
    for b = 1:n_blocks
        l = sscanf(fgetl(fid),'%d %d %d %d');
        if l(3)
            error('Cant read parametric nodes')
        end
        pts = textscan(fid, '%d %f %f %f', l(4));
        node_numbers = [node_numbers; pts{1}];
        m.nodes = [m.nodes; [pts{2} pts{3} pts{4}]];
    end
    [node_numbers, order] = sort(node_numbers);
    m.nodes = m.nodes(order, :);
    if (min(node_numbers) < 1) || ...
       (max(node_numbers)> number_of_nodes) || ...
       (length(node_numbers) ~= number_of_nodes)
        error('node numbers have to be a continuous list of indexes starting at 1');
    end

    % read end line
    tline = fgetl(fid);
    if ~strcmp('$EndNodes', tline)
        error('End tag $EndNodes expected'); 
    end

    

function [m, continous_elm_numbers] = read_elements_binary(m, fid)
    % get number of elements    
    tline = fgetl(fid);
    number_of_elements = sscanf(tline,'%d');
    if ~isnumeric( number_of_elements )
        error('number of elements is not a number');
    end
    
    % loop to read in all elements
    
    m.triangles = [];
    m.triangle_regions = [];
    m.tetrahedra = [];
    m.tetrahedron_regions = [];

    elm_remaining = number_of_elements;
    element_numbers=[];
    element_types = [];
    continous_elm_numbers=true;
    while elm_remaining > 0
    
        type = fread(fid, 1, '*int32');
        nr_elm_follow = fread(fid, 1, '*int32');
        nr_tags = fread(fid, 1, '*int32');
        if (nr_tags ~= 2) 
            error('nr_tags must be equal to 2');
        end

        elm_remaining = elm_remaining - nr_elm_follow;
        % read all the elements in one go
        % idx, tag1, tag2, node1, node2, node3, ... up to nr_nodes_in_element nodes
        buff = fread(fid, [nr_nodes_in_element(type)+3, nr_elm_follow], '*int32'); % +3 is for idx, tag1, tag2
        
        element_numbers = [element_numbers; buff(1,:)'];
        element_types = [element_types; type];
        switch type
            case 2
                if ~isempty(buff)
                    m.triangles = [m.triangles; buff(4:end,:)'];
                    m.triangle_regions = [m.triangle_regions; buff(2,:)']; % this is tag1
                end;
                
            case 4
                if ~isempty(buff)
                    m.tetrahedra = [m.tetrahedra; buff(4:end,:)'];
                    m.tetrahedron_regions = [m.tetrahedron_regions; buff(2,:)']; % this is tag1
                end;    
                
            otherwise
                continous_elm_numbers=false;
                warning(['found element type ' num2str(type) '; skipped during reading! Only reading triangles and tetrahedra, not reading element data']);
        end
    end

    if (min(element_numbers) < 1) || ...
       (max(element_numbers) > number_of_elements)
       
       continous_elm_numbers=false;
       warning('only reading tetrahedra and triangles, not reading element_data')
    end
    
    if any(element_types==2) && any(element_types==4) && ...
       find(element_types==2, 1, 'last')>find(element_types==4, 1)
   
       continous_elm_numbers=false;
       warning('only reading tetrahedra and triangles, not reading element_data')   
    end;
        
    % sometimes there is an end of line after $Elements, sometimes there is
    % not... handle that
    read_LF(fid);
    
    % read end line
    tline = fgetl(fid);
    if ~strcmp('$EndElements', tline)
        error(['End tag $EndElements expected, but I read: ' tline]);
    end

       
function [m, continous_elm_numbers] = read_elements_binary4(m, fid)
    n_blocks = fread(fid,1, 'uint64');
    number_of_elements = fread(fid,1, 'uint64');

    m.triangles = [];
    m.triangle_regions = [];
    m.tetrahedra = [];
    m.tetrahedron_regions = [];
    tr_numbers = [];
    th_numbers = [];
    for b = 1:n_blocks
        tag = fread(fid,1, 'int');
        fread(fid,1, 'int');
        type = fread(fid,1, 'int');
        n_in_block = fread(fid,1, 'uint64');
        buff = fread(fid, [nr_nodes_in_element(type)+1, n_in_block], '*int')';
        switch type
            case 2
                if ~isempty(buff)
                    tr_numbers = [tr_numbers; buff(:, 1)];
                    m.triangles = [m.triangles; buff(:, 2:end)];
                    m.triangle_regions = [m.triangle_regions; tag*ones(n_in_block,1, 'int32')];
                end                
            case 4
                if ~isempty(buff)
                    th_numbers = [th_numbers; buff(:, 1)];
                    m.tetrahedra = [m.tetrahedra; buff(:, 2:end)];
                    m.tetrahedron_regions = [m.tetrahedron_regions; tag*ones(n_in_block,1, 'int32')];
                end
            otherwise
                continous_elm_numbers=false;
                warning(['found element type ' num2str(type) '; skipped during reading! Only reading triangles and tetrahedra, not reading element data']);
        end
    end
    [tr_numbers, order] = sort(tr_numbers);
    m.triangles = m.triangles(order, :);
    m.triangle_regions = m.triangle_regions(order);
    [th_numbers, order] = sort(th_numbers);
    m.tetrahedra = m.tetrahedra(order, :);
    m.tetrahedron_regions = m.tetrahedron_regions(order);
    element_numbers = [tr_numbers; th_numbers];
    if (min(element_numbers) < 1) || ...
       (max(element_numbers) > number_of_elements)
       
       continous_elm_numbers=false;
       warning('only reading tetrahedra and triangles, not reading element_data')
    end
        
    % sometimes there is an end of line after $Elements, sometimes there is
    % not... handle that
    read_LF(fid);
    
    % read end line
    tline = fgetl(fid);
    if ~strcmp('$EndElements', tline)
        error(['End tag $EndElements expected, but I read: ' tline]);
    end        
    
function [m, continous_elm_numbers] = read_elements(m, fid)
    % get number of elements
    tline = fgetl(fid);
    number_of_elements = sscanf(tline,'%d');
    if ~isnumeric( number_of_elements )
        error('number of elements is not a number');
    end
    
    c = textscan(fid,'%d %d %d %d %*d %d %d %d %d');
    element_numbers=c{1};
    number_of_tags=c{3};
    
    if length(element_numbers) ~= number_of_elements
        error('wrong number of elements found');
    end;
    
    if any(number_of_tags ~= 2)
        error('two tag numbers are expected');
    end;
    
    continous_elm_numbers=true;
    if (min(element_numbers) < 1) || ...
       (max(element_numbers)> number_of_elements) || ...
       any(~((c{2} == 2)|(c{2} ==4))) 
       
       continous_elm_numbers=false;
       warning('only reading tetrahedra and triangles, not reading element_data')
    end
    
    if any(c{2}==2) && any(c{2}==4) && ...
       find(c{2}==2, 1, 'last')>find(c{2}==4, 1)
   
       continous_elm_numbers=false;
       warning('only reading tetrahedra and triangles, not reading element_data')   
    end;
    
	% Separate triangles from tetrahedra, throw away rest
    tri = c{2}==2;
    tetr = c{2}==4;
    m.triangles = [c{5}(tri) c{6}(tri) c{7}(tri)];
    m.tetrahedra = [c{5}(tetr) c{6}(tetr) c{7}(tetr) c{8}(tetr)];
    m.triangle_regions = c{4}(tri);
    m.tetrahedron_regions = c{4}(tetr);
    
    % read end line
    tline = fgetl(fid);
    if ~strcmp('$EndElements', tline)
        error(['End tag $EndElements expected, but I read: ' tline]);
    end


    
function [m, continous_elm_numbers] = read_elements4(m, fid)
    % get number of elements
    tline = fgetl(fid);
    l = sscanf(tline,'%d %d');
    n_blocks = l(1);
    number_of_elements = l(2);
    if ~isnumeric( number_of_elements )
        error('number of elements is not a number');
    end
    
    m.triangles = [];
    m.tetrahedra = [];
    m.triangle_regions = [];
    m.tetrahedron_regions = [];
    tr_numbers = [];
    th_numbers = [];
    for b = 1:n_blocks
        l = sscanf(fgetl(fid),'%d %d %d %d');
        if l(3) == 2
            c = textscan(fid,'%d %d %d %d', l(4));
            tr_numbers=[tr_numbers; c{1}];
            m.triangles = [m.triangles; [c{2} c{3} c{4}]];
            m.triangle_regions = [m.triangle_regions; l(1)*ones(l(4),1, 'int32')];
        elseif l(3) == 4
            c = textscan(fid,'%d %d %d %d %d', l(4));
            th_numbers=[th_numbers; c{1}];
            m.tetrahedra = [m.tetrahedra; [c{2} c{3} c{4} c{5}]];
            m.tetrahedron_regions = [m.tetrahedron_regions; l(1)*ones(l(4),1, 'int32')];
        else
            warning('only reading tetrahedra and triangles, not reading element_data')
            c = textscan(fid,'%d %*d %d', l(4));
            continue
        end
    end
    [tr_numbers, order] = sort(tr_numbers);
    m.triangles = m.triangles(order, :);
    m.triangle_regions = m.triangle_regions(order);
    [th_numbers, order] = sort(th_numbers);
    m.tetrahedra = m.tetrahedra(order, :);
    m.tetrahedron_regions = m.tetrahedron_regions(order);
    element_numbers = [tr_numbers; th_numbers];
    continous_elm_numbers=true;
    if (min(element_numbers) < 1) || ...
       (max(element_numbers)> number_of_elements) || ...
       any(~((c{2} == 2)|(c{2} ==4))) 
       
       continous_elm_numbers=false;
       warning('only reading tetrahedra and triangles, not reading element_data')
    end
    
    if any(c{2}==2) && any(c{2}==4) && ...
       find(c{2}==2, 1, 'last')>find(c{2}==4, 1)
   
       continous_elm_numbers=false;
       warning('only reading tetrahedra and triangles, not reading element_data')   
    end
    
    % read end line
    tline = fgetl(fid);
    if ~strcmp('$EndElements', tline)
        error(['End tag $EndElements expected, but I read: ' tline]);
    end



function data = read_node_data(fid, number_of_nodes, isbinary)

    if isbinary
        [dataIn, name, node_element_numbers]= read_data_binary(fid, '$NodeData');
    else
        [dataIn, name, node_element_numbers]= read_data(fid, '$NodeData');
    end;
    
    if (node_element_numbers(1) < 1) || ...
       (node_element_numbers(2) > number_of_nodes) || ...
       (size(dataIn,1) ~= number_of_nodes)
        error('node data has to contain exactly one data point for each node');
    end
        
    % sort data into struct
    data.name = name;
    data.data = dataIn;
       
    
function data = read_element_data(fid, number_of_triangles, number_of_tetrahedra, continous_elm_numbers, isbinary)

    % still reading element_data also for non-continous element number to
    % move file position indicator
    if isbinary 
        [dataIn, name, node_element_numbers]= read_data_binary(fid, '$ElementData');
    else
        [dataIn, name, node_element_numbers]= read_data(fid, '$ElementData');
    end
    
    if continous_elm_numbers
        % sort data into struct
        data.name = name;
        
        if (node_element_numbers(1) == 1) && ...
           (node_element_numbers(2) == number_of_triangles)
            
            if size(dataIn,1) ~= number_of_triangles
                error('element data has to contain one data point for each triangle');
            end;
            
            data.tridata = dataIn;
            data.tetdata = [];
        elseif (node_element_numbers(1) == 1) && ...
               (node_element_numbers(2) == number_of_triangles+number_of_tetrahedra)
           
            if size(dataIn,1) ~= number_of_triangles+number_of_tetrahedra
                error('element data has to contain one data point for each element');
            end;

            data.tridata = dataIn(1:number_of_triangles,:); 
            data.tetdata = dataIn(number_of_triangles+1:end,:);
        elseif (node_element_numbers(1) == number_of_triangles+1) && ...
               (node_element_numbers(2) == number_of_triangles+number_of_tetrahedra)
           
            if size(dataIn,1) ~= number_of_tetrahedra
                error('element data has to contain one data point for each tetrahedron');
            end;
            
            data.tridata = [];
            data.tetdata = dataIn;
        else
           error('reading element data did not succeed');
        end;
    end;

                       
       
function [data, name, node_element_numbers] = read_data_binary(fid, data_type)
    
    % read string tags (including name)
    tline = fgetl(fid);
    if sscanf(tline,'%d') ~= 1; error('nr_string_tags should always be 1'); end
    name = fgetl(fid);

    % read real tags
    tline = fgetl(fid);
    if sscanf(tline,'%d') ~= 1; error('nr_real_tags should always be 1'); end
    tline = fgetl(fid);
    
    % read integer tags (size of data)
    tline = fgetl(fid);
    nlines = sscanf(tline,'%d');
    if nlines ~= 3 && nlines ~= 4
        error('nr_int_tags should always be 3 or 4');
    end
    
    for i=1:nlines
        tline = fgetl(fid);
        int_tags(i) = sscanf(tline, '%d');
    end;
    
    % nr of field components (1, 3, 9) for scalar, vector, tensor
    comp = int_tags(2);
    nr_data = int_tags(3);
    
    % read float data
    node_element_numbers = [0 0]; % node or element numbers of first and last data line
    if strcmp(data_type, '$ElementNodeData')
        % for $ElementNodeData, it is written as:
        % elm-number number-of-nodes-per-element value
        % where value has comp*number_of_nodes_per_element value
        % since we only write to tetrahedra, I will assume there are no triangles
        % with $ElementNodeData

        % read first two columns
        fread(fid, 2, '*uint32');
        data = fread(fid, [4*comp nr_data], [num2str(4*comp) '*double'], 8)';
        fseek(fid, -8, 'cof');
    elseif strcmp(data_type, '$NodeData') || strcmp(data_type, '$ElementData')
        % read the node or element number of the first data line
        node_element_numbers(1) = fread(fid, 1 , 'int');
        
        % read data
        data = fread(fid, [comp nr_data], [num2str(comp) '*double'], 4)';
        
        % read the node or element number of the last data line
        fseek(fid, -(2*4+comp*8), 'cof');
        node_element_numbers(2) = fread(fid, 1 , 'int', comp*8);        
    else
        error('still need to code other data types');
    end
                
    % update name
    if length(name) > 2
        name = name(2:end-1);
    else
        name = '';
    end
     
    % read last line
    tline = fgetl(fid);
    if ~strcmp(['$End' data_type(2:end)], tline)
        error(['End tag $End' data_type ' expected']);
    end



function [data, name, node_element_numbers]= read_data(fid, data_type)

    % read string tags (including name)
    tline = fgetl(fid);
    if sscanf(tline,'%d') ~= 1; error('nr_string_tags should always be 1'); end
    name = fgetl(fid);

    % read real tags
    tline = fgetl(fid);
    if sscanf(tline,'%d') ~= 1; error('nr_real_tags should always be 1'); end
    tline = fgetl(fid);
    
    % read integer tags (size of data)
    tline = fgetl(fid);
    nlines = sscanf(tline,'%d');
    if nlines ~= 3 && nlines ~= 4
        error('nr_int_tags should always be 3 or 4');
    end
    
    for i=1:nlines
        tline = fgetl(fid);
        int_tags(i) = sscanf(tline, '%d');
    end;
    
    comp = int_tags(2);

    % read data
    t_str = repmat(' %f', 1, comp);
    c = textscan(fid, ['%d' t_str '']);
    if size(c{1},1) ~= int_tags(3)
        error(['number of data lines read does not correspond to ' num2str(int_tags(3))]);
    end
    
    if length(name) > 2
        name = name(2:end-1);
    else
        name = '';
    end
    
    node_element_numbers(1)=c{1}(1);
    node_element_numbers(2)=c{1}(end);
    
    data = [];
    for i=1:comp
        data = [data c{i+1}];
    end
        
    % read last line
    tline = fgetl(fid);
    if ~strcmp(['$End' data_type(2:end)], tline)
        error(['End tag $End' data_type ' expected']);
    end

    
function read_LF(fid)
    a = fread(fid, 1, 'char');
    if a == 10
        disp('LF found, but I should be able to handle it');
    elseif a == 36
        % go back 1 byte, I'm reading a "$" character
        fseek(fid, -1, 'cof');
    else
        error(['Dont know what to do with this: ' a]);
    end
    
function n = nr_nodes_in_element(t)
    nne = ...
     [ 2,  3,  4,  4,  8,  6,  5,  3,  6,  9, ...
      10, 27, 18, 14,  1,  8, 20, 15, 13,  9, ...
      10, 12, 15, 15, 21,  4,  5,  6, 20, 35, ...
      56,  0,  0,  0,  0,  0,  0,  0,  0,  0, ...
       0,  0,  0,  0,  0,  0,  0,  0,  0,  0, ...
       0,  0,  0,  0,  0,  0,  0,  0,  0,  0, ...
       0,  0,  0,  0,  0,  0,  0,  0,  0,  0, ...
       0,  0,  0,  0,  0,  0,  0,  0,  0,  0, ...
       0,  0,  0,  0,  0,  0,  0,  0,  0,  0, ...
       0, 64,125];
   n = nne(t);
