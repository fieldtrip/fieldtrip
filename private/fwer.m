function [h] = fwer(p, q);

% FWER family-wise error rate control using Bonferoni method
%
% Use as
%   h = fwer(p, q)

% Copyright (C) 2005, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% count the number of voxels
V = length(p);

% threshold the probability of all voxels
h = (p<=(q/V));
