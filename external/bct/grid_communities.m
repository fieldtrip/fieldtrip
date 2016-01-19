function [X,Y,indsort] = grid_communities(c)
%GRID_COMMUNITIES       outline communities along diagonal
%
%   [X Y INDSORT] = GRID_COMMUNITIES(C) takes a vector of community
%   assignments C and returns three output arguments for visualizing the
%   communities. The third is INDSORT, which is an ordering of the vertices
%   so that nodes with the same community assignment are next to one
%   another. The first two arguments are vectors that, when overlaid on the
%   adjacency matrix using the PLOT function, highlight the communities.
%
%   Example:
%
%   >> load AIJ;                                % load adjacency matrix
%   >> [C,Q] = modularity_louvain_und(AIJ);     % get community assignments
%   >> [X,Y,INDSORT] = fcn_grid_communities(C); % call function
%   >> imagesc(AIJ(INDSORT,INDSORT));           % plot ordered adjacency matrix
%   >> hold on;                                 % hold on to overlay community visualization
%   >> plot(X,Y,'r','linewidth',2);             % plot community boundaries
%
%   Inputs:     C,       community assignments
%
%   Outputs:    X,       x coor
%               Y,       y coor
%               INDSORT, indices
%
%   Richard Betzel, Indiana University, 2012
%

%#ok<*AGROW>

nc = max(c);
[c,indsort] = sort(c);

X = [];
Y = [];
for i = 1:nc
    ind = find(c == i);
    if ~isempty(ind)
        mn = min(ind) - 0.5;
        mx = max(ind) + 0.5;
        x = [mn mn mx mx mn NaN];
        y = [mn mx mx mn mn NaN];
        X = [X, x]; 
        Y = [Y, y];
    end
end