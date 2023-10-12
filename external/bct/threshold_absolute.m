function W = threshold_absolute(W, thr)
% THRESHOLD_ABSOLUTE    Absolute thresholding
% 
%   W_thr = threshold_absolute(W, thr);
%
%   This function thresholds the connectivity matrix by absolute weight
%   magnitude. All weights below the given threshold, and all weights
%   on the main diagonal (self-self connections) are set to 0.
%
%   Inputs: W           weighted or binary connectivity matrix
%           thr         weight treshold
%
%   Output: W_thr       thresholded connectivity matrix
%
%
%   Mika Rubinov, UNSW, 2009-2010

W(1:size(W,1)+1:end)=0;                 %clear diagonal
W(W<thr)=0;                             %apply threshold