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
%                               p=0 (no weights preserved)
%
%   Output: W_thr,  thresholded connectivity matrix
%
%
%   Mika Rubinov, U Cambridge,
%   Roan LaPlante, Martinos Center, MGH
%   Zitong Zhang, Penn Engineering

%   Modification history:
%   2010: Original (MR)
%   2012: Bug fix for symmetric matrices (RLP)
%   2015: Improved symmetricity test (ZZ)

n=size(W,1);                                %number of nodes
W(1:n+1:end)=0;                             %clear diagonal

if max(max(abs(W-W.'))) < 1e-10;            %if symmetric matrix
    W=triu(W);                              %ensure symmetry is preserved
    ud=2;                                   %halve number of removed links
else
    ud=1;
end

ind=find(W);                                %find all links
E=sortrows([ind W(ind)], -2);               %sort by magnitude
en=round((n^2-n)*p/ud);                     %number of links to be preserved

W(E(en+1:end,1))=0;                         %apply threshold

if ud==2                                    %if symmetric matrix
    W=W+W.';                                %reconstruct symmetry
end
