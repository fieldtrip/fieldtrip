function [pnt, dhk] = icosahedron();

% ICOSAHEDRON2562 creates a 4-fold refined icosahedron

% Copyright (C) 2003, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

[pnt, dhk] = icosahedron;
[pnt, dhk] = refine(pnt, dhk);
[pnt, dhk] = refine(pnt, dhk);
[pnt, dhk] = refine(pnt, dhk);
[pnt, dhk] = refine(pnt, dhk);

pnt = pnt ./ repmat(sqrt(sum(pnt.^2,2)), 1,3);
