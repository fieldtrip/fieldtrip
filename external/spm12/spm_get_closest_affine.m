function [M,R] = spm_get_closest_affine(x,y,w)
% Determine the affine transform mapping x to y
% FORMAT [M,R] = spm_get_closest_affine(X,Y,W)
% X  - n1*n2*n3*3 array of floats representing coordinates.
% Y  - n1*n2*n3*3 array of floats representing coordinates.
% W  - n1*n2*n3   array of floats representing weights.
%
% M  - an affine transform
% R  - a rigid-body transform
%
% The code treats X and Y as reshaped versions (n1*n2*n3) x 3,
% and W as a column vector.
% 
% It generates XX = [X 1]'*diag(W)*[X 1]
% and          XY = [X 1]'*diag(W)*[Y 1]
% 
% These can then be used to compute an affine transform (M),
% by M = (XX\XY)'
% A weighted procrustes decomposition is also performed,
% so that a rigid-body transform matrix (R) is returned.
%
% If W is empty or not passed, then it is assumed to be all ones.
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


XX = zeros(4);
XY = zeros(4);
d  = size(x);
o  = ones(d(1)*d(2),1);
for k=1:size(x,3)
    xk  = reshape(x(:,:,k,:),[d(1)*d(2),3]);
    yk  = reshape(y(:,:,k,:),[d(1)*d(2),3]);
    if nargin>=3 && ~isempty(w)
        wk = reshape(w(:,:,k), [d(1)*d(2),1]);
    else
        wk = o;
    end
    msk = find(all(isfinite(xk),2) & all(isfinite(yk),2) & isfinite(wk));
    X   = [xk(msk,:), o(msk)];
    Y   = [yk(msk,:), o(msk)];
    wk  = wk(msk);
    XX  = XX + double(X'*bsxfun(@times,wk,X));
    XY  = XY + double(X'*bsxfun(@times,wk,Y));
end
M = (XX\XY)';

if nargout>1
    % Procrustes decomposition
    XX1 = XX - XX(:,4)*XX(:,4)'/XX(4,4);
    XY1 = XY - XY(:,4)*XY(4,:) /XY(4,4);
    Z   = (XX1(1:3,1:3)\XY1(1:3,1:3))';
    [U,S,V] = svd(Z);                   % Decompose into rotate, zoom and rotate.
    R   = [U*V' zeros(3,1);0 0 0 1];    % Pure rotation (by taking out the zoom)
    T1  = [eye(4,3) -XY(:,4) /XY(4,4)]; % Initial translation of centre of mass to origin.
    T2  = [eye(4,3) -XY(4,:)'/XY(4,4)]; % Final translation of origin to centre of mass.
    R   = T2 * R * T1;
end
