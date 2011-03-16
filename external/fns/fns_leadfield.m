function L = fns_leadfield(data,compress,node_sizes,row,col,slice)
% 
% fns_leadfield is a function to extract the lead field matrix from the
% reciprocity data matrix
%
% L = fns_leadfield(recipdata,compress,node_sizes,row,col,slice)
%
% L:            Lead field matrix
% recipdata:    Reciprocity data matrix
% compress:     Compress vector
% node_sizes:   The sizes of the node volume
% row:          Row position of the given node
% col:          Column position of the given node
% slice:        Slice position of the given node


% $Copyright (C) 2010 by Hung Dang$


% Position of the head and tail in the x, y and z direction
ROW = node_sizes(1);
COL = node_sizes(2);
N = length(compress);

% Add 1 to convert between MATLAB index and C index
nnode = row + col * ROW + slice * COL * ROW + 1; 

% Compute the lead field matrix for the given node
if (nnode <= N)
    head_x = compress(nnode + 1);
    tail_x = compress(nnode - 1);
    head_y = compress(nnode + ROW);
    tail_y = compress(nnode - ROW);
    head_z = compress(nnode + COL * ROW);
    tail_z = compress(nnode - COL * ROW);
        
    % Extract the lead field matrix
    if ((head_x > 0) && (tail_x > 0) && (head_y > 0) && (tail_y > 0) && (head_z > 0) && (tail_z > 0))
        % Extract the lead field matrix for x, y and z component
        lx = data(:,head_x) - data(:,tail_x);
        ly = data(:,head_y) - data(:,tail_y);
        lz = data(:,head_z) - data(:,tail_z);
        % Form the final lead field matrix
        L = [lx,ly,lz];
    else
        L = [];
    end
else
    L = [];
end
