function [s] = nanstd(x, varargin)

% [S] = NANSTD(X, NORMALIZEFLAG, DIM, FLAG) computes the standard deviation
% across all cells in x along the dimension dim. Normalizeflag = 1 normalizes
% by N, normalizeflag = [] or 0 normalizes by N-1. 
% 
% X should be an linear cell-array of matrices for which the size in at 
% least one of the dimensions should be the same for all cells. If flag==1, the mean will
% be subtracted first (default behaviour, but to save time on already demeaned data, it
% can be set to 0).

s = sqrt(nanvar(x, varargin{:}));
