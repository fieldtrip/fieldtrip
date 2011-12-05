function W = threshold_proportional(W, p)
%THRESHOLD_PROPORTIONAL     Proportional thresholding
%
%   W_thr = threshold_proportional(W, p);
%
%   This function "thresholds" the connectivity matrix by preserving a
%   proportion p (0<p<1) of the strongest weights. All other weights, and
%   all weights on the main diagonal (self-self connections) are set to 0.
%
%   Inputs: W,      weighted or binary connectivity matrix
%           p,      proportion of weights to preserve
%                       range:  p=1 (all weights preserved) to
%                               p=0 (no weights removed)
%
%   Output: W_thr,  thresholded connectivity matrix
%
%
%   Mika Rubinov, UNSW, 2010

n=size(W,1);                                %number of nodes
W(1:n+1:end)=0;                             %clear diagonal
ind=find(W);                                %find all links
E=sortrows([ind W(ind)], -2);               %sort by magnitude
en=round((n^2-n)*p);                        %number of links to be preserved

if rem(en,2) && isequal(W,W.')              %if symmetric matrix
    en=en+1;                                %ensure symmetry is preserved
end

W(E(en+1:end,1))=0;                         %apply threshold