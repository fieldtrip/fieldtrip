function [h] = fwer(p, q);

% FWER family-wise error rate control using Bonferoni method
%
% Use as
%   h = fwer(p, q)

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: fwer.m,v $
% Revision 1.1  2005/11/08 16:02:58  roboos
% initial implementation
%

% count the number of voxels
V = length(p);

% threshold the probability of all voxels
h = (p<=(q/V));
