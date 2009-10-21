function [h] = fdr(p, q);

% FDR false discovery rate
%
% Use as
%   h = fdr(p, q)
%
% This implements
%   Genovese CR, Lazar NA, Nichols T.
%   Thresholding of statistical maps in functional neuroimaging using the false discovery rate.
%   Neuroimage. 2002 Apr;15(4):870-8.

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: fdr.m,v $
% Revision 1.2  2006/06/07 12:57:09  roboos
% also support n-D input arrays
%
% Revision 1.1  2005/11/08 16:02:58  roboos
% initial implementation
%

% convert the input into a row vector
dim = size(p);
p = reshape(p, 1, numel(p));

% sort the observed uncorrected probabilities
[ps, indx] = sort(p);

% count the number of voxels
V = length(p);

% compute the threshold probability for each voxel
pi = ((1:V)/V)  * q / c(V);

h = (ps<=pi);

% undo the sorting
[dum, unsort] = sort(indx);
h = h(unsort);

% convert the output back into the original format
h = reshape(h, dim);

function s = c(V)
% See Genovese, Lazar and Holmes (2002) page 872, second column, first paragraph
if V<1000
  % compute it exactly
  s = sum(1./(1:V));
else
  % approximate it
  s = log(V) + 0.57721566490153286060651209008240243104215933593992359880576723488486772677766467093694706329174674951463144724980708248096050401448654283622417399764492353625350033374293733773767394279259525824709491600873520394816567;
end

