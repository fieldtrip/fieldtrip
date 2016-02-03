function W = weight_conversion(W, wcm)
% WEIGHT_CONVERSION    Conversion of weights in input matrix
%
%   W_bin = weight_conversion(W, 'binarize');
%   W_nrm = weight_conversion(W, 'normalize');
%   L = weight_conversion(W, 'lengths');
%   W_fix = weight_conversion(W, 'autofix');
%
%   This function may either binarize an input weighted connection matrix,
%   normalize an input weighted connection matrix, convert an input
%   weighted connection matrix to a weighted connection-length matrix, or
%   fix common connection problems in binary or weighted connection matrices.
%
%       Binarization converts all present connection weights to 1.
%
%       Normalization rescales all weight magnitudes to the range [0,1] and
%   should be done prior to computing some weighted measures, such as the
%   weighted clustering coefficient.
%
%       Conversion of connection weights to connection lengths is needed
%   prior to computation of weighted distance-based measures, such as
%   distance and betweenness centrality. In a weighted connection network,
%   higher weights are naturally interpreted as shorter lengths. The
%   connection-lengths matrix here is defined as the inverse of the
%   connection-weights matrix. 
%
%       Autofix removes all Inf and NaN values, remove all self connections 
%   (sets all weights on the main diagonal to 0), ensures that symmetric matrices 
%   are exactly symmetric (by correcting for round-off error), and ensures that 
%   binary matrices are exactly binary (by correcting for round-off error).
%
%   Inputs: W           binary or weighted connectivity matrix
%           wcm         weight-conversion command - possible values:
%                           'binarize'      binarize weights
%                           'normalize'     normalize weights

%                           'lengths'       convert weights to lengths
%                           'autofix'       fixes common weights problems
%
%   Output: W_          output connectivity matrix
%
%
%   Mika Rubinov, U Cambridge, 2012

%   Modification History:
%   Sep 2012: Original
%   Jan 2015: Added autofix feature.

switch wcm
    case 'binarize'
        W=double(W~=0);         % binarize
    case 'normalize'
        W=W./max(abs(W(:)));    % scale by maximal weight
    case 'lengths'
        E=find(W); 
        W(E)=1./W(E);           % invert weights
    case 'autofix'
        % clear diagonal
        n = length(W);
        W(1:n+1:end)=0;

        % remove Infs and NaNs
        idx = isnan(W) | isinf(W);
        if any(any(idx));
            W(idx)=0;
        end

        % ensure exact binariness
        U = unique(W);
        if numel(U)>1
            idx_0 = abs(U)<1e-10;
            idx_1 = abs(U-1)<1e-10;
            if all(idx_0 | idx_1)
                W(idx_0)=0;
                W(idx_1)=1;
            end
        end

        % ensure exact symmetry
        if ~isequal(W,W.');
            if max(max(abs(W-W.'))) < 1e-10;            
                W=(W+W).'/2;
            end
        end
    otherwise
        error('Unknown weight-conversion command.')
end
