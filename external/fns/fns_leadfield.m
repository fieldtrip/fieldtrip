function L = fns_leadfield(forwardSolutions, modelSizes, dipoles)
% 
% fns_leadfield is a function that compute the lead field matrix from the
% potential data obtained using Finite Different Method.
%
% L = fns_leadfield(forwardSolutions, row, col, slice)
%
% L:                Lead field matrix
% forwardSolutions: Forward solution data
% modelSizes:       Size of a volume image
% dipoles:          MRI coordinate locations of dipoles

% $Copyright (C) 2010 by Hung Dang$

% 3D model of a dipole used in FNS
% 
%          0
%          | 0
%          |/
%    0-----x-----0
%         /|
%        0 |
%          0
%
% 

ROW = modelSizes(1) + 1;
COL = modelSizes(2) + 1;

[M, N] = size(forwardSolutions);
[numberOfDipoles, dum] = size(dipoles);

% Allocate memory for lead field matrix
L = zeros(M, numberOfDipoles * 3);

% Compute lead field matrices for dipoles
for dipoleIdx = 1:numberOfDipoles    
    % Add 1 to convert between MATLAB index and C index
    nnode = dipoles(dipoleIdx, 1) + ...
            dipoles(dipoleIdx, 2) * ROW + ...
            dipoles(dipoleIdx, 3) * COL * ROW + 1; 

    % Compute the lead field matrix for the given node
    if ((nnode >= COL * ROW + 1) && (nnode <= (N - COL * ROW)))
        head_x = nnode + 1;
        tail_x = nnode - 1;
        head_y = nnode + ROW;
        tail_y = nnode - ROW;
        head_z = nnode + COL * ROW;
        tail_z = nnode - COL * ROW;        
        
        % Extract the lead field matrix for x, y and z component
        L(:, dipoleIdx * 3 - 2) = forwardSolutions(:,head_x) - forwardSolutions(:,tail_x);
        L(:, dipoleIdx * 3 - 1) = forwardSolutions(:,head_y) - forwardSolutions(:,tail_y);
        L(:, dipoleIdx * 3)     = forwardSolutions(:,head_z) - forwardSolutions(:,tail_z);

    else
        error('FNS:InvalidDipoleLocations', ...
              sprintf(['Dipoles(%d) at (%d, %d, %d) does not belong to the head model. ', ...
                       'You cannot compute lead field for this dipole.'], ...
                      dipoleIdx, ...
                      dipoles(dipoleIdx, 1), ...
                      dipoles(dipoleIdx, 2), ...
                      dipoles(dipoleIdx, 3)));
    end
end