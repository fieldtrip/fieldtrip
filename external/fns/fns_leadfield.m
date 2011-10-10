function L = fns_leadfield(forwardSolutions, modelSizes, row, col, slice)
% 
% fns_leadfield is a function that compute the lead field matrix from the
% potential data obtained using Finite Different Method.
%
% L = fns_leadfield(forwardSolutions, row, col, slice)
%
% L:                Lead field matrix
% forwardSolutions: Forward solution data
% row:              Row position of the given node
% col:              Column position of the given node
% slice:            Slice position of the given node

% $Copyright (C) 2010 by Hung Dang$

% Position of the head and tail in the x, y and z direction
ROW = modelSizes(1);
COL = modelSizes(2);
[M, N] = size(forwardSolutions); %#ok<NASGU>

% Add 1 to convert between MATLAB index and C index
nnode = row + col * ROW + slice * COL * ROW + 1; 

% Compute the lead field matrix for the given node
if (nnode <= N)
%     head_x = forwardSolutions.Compress(nnode + 1);
%     tail_x = forwardSolutions.Compress(nnode - 1);
%     head_y = forwardSolutions.Compress(nnode + ROW);
%     tail_y = forwardSolutions.Compress(nnode - ROW);
%     head_z = forwardSolutions.Compress(nnode + COL * ROW);
%     tail_z = forwardSolutions.Compress(nnode - COL * ROW);

    head_x = nnode + 1;
    tail_x = nnode - 1;
    head_y = nnode + ROW;
    tail_y = nnode - ROW;
    head_z = nnode + COL * ROW;
    tail_z = nnode - COL * ROW;


    % Extract the lead field matrix
    if (...
            (head_x > 0) && (tail_x > 0) && ...
            (head_y > 0) && (tail_y > 0) && ...
            (head_z > 0) && (tail_z > 0))
        
        % Extract the lead field matrix for x, y and z component
        lx = forwardSolutions(:,head_x) - forwardSolutions(:,tail_x);
        ly = forwardSolutions(:,head_y) - forwardSolutions(:,tail_y);
        lz = forwardSolutions(:,head_z) - forwardSolutions(:,tail_z);
        
        % Form the final lead field matrix
        L = [lx,ly,lz];
    else
        L = [];
    end
else
    L = [];
end
