function [X,Y,Z] = adjacency_plot_und(aij,coor)
%ADJACENCY_PLOT     quick visualization tool
%
%   [X,Y,Z] = ADJACENCY_PLOT(AIJ,COOR) takes adjacency matrix AIJ and node
%   spatial coordinates COOR and generates three vectors that can be used
%   for quickly plotting the edges in AIJ. If no coordinates are specified,
%   then each node is assigned a position on a circle. COOR can, in 
%   general, be 2D or 3D.
%
%   Example:
%
%   >> load AIJ;                                % load your adjacency matrix
%   >> load COOR;                               % load 3D coordinates for each node
%   >> [x,y,z] = adjacency_plot_und(AIJ,COOR);  % call function
%   >> plot3(x,y,z);                            % plots network as a single line object
%
%   If COOR were 2D, the PLOT3 command changes to a PLOT command.
%
%   NOTE: This function is similar to MATLAB's GPLOT command.
%
%   Richard Betzel, Indiana University, 2013

n = length(aij);
if nargin < 2
    coor = zeros(n,2);
    for i = 1:n
        coor(i,:) = [cos(2*pi*(i - 1)./n), sin(2*pi*(i - 1)./n)];
    end
end

[i,j] = find(triu(aij,1));
[~, p] = sort(max(i,j));
i = i(p);
j = j(p);

X = [ coor(i,1) coor(j,1)]';
Y = [ coor(i,2) coor(j,2)]';
if size(coor,2) == 3
    Z = [ coor(i,3) coor(j,3)]';
end
if isfloat(coor) || nargout ~= 0
    X = [X; NaN(size(i))'];
    Y = [Y; NaN(size(i))'];
    if size(coor,2) == 3
        Z = [Z; NaN(size(i))'];
    end
end

X = X(:);
Y = Y(:);
if size(coor,2) == 3
    Z = Z(:);
end
