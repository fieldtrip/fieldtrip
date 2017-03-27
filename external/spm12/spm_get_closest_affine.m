function [M,R] = spm_get_closest_affine(x,y,w1,w2)
% Determine the affine transform mapping x to y
% FORMAT [M,R] = spm_get_closest_affine(X,Y,W1,W2)
% X  - n1*n2*n3*3 array of floats representing coordinates.
% Y  - n1*n2*n3*3 array of floats representing coordinates.
% W1 - n1*n2*n3   array of floats representing weights.
% W2 - n1*n2*n3   array of floats representing weights.
%
% M  - an affine transform
% R  - a rigid-body transform
%
% The code treats X and Y as reshaped versions (n1*n2*n3) x 3,
% and W1 and W2 as column vectors.
% 
% It generates XX = [diag(W1)*X W1]'*diag(W2)*[diag(W1)*X W1]
% and          XY = [diag(W1)*X W1]'*diag(W2)*[Y W1]
% 
% These can then be used to compute an affine transform (M),
% by M = (XX\XY)'
% A weighted procrustes decomposition is also performed,
% so that a rigid-body transform matrix (R) is returned.
%
% If W1 or W2 are empty or not passed, then they are assumed
% to be all ones.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_get_closest_affine.m 6137 2014-08-19 12:43:11Z john $
 
XX = zeros(4);
XY = zeros(4);
d  = size(x);
o  = ones(d(1)*d(2),1);
for k=1:size(x,3),
    xk  = reshape(x(:,:,k,:),[d(1)*d(2),3]);
    if (nargin<3 || isempty(w1)) && (nargin<4 || isempty(w2)),
        ox = o;
        oy = o;
    else
        if nargin>=4 && ~isempty(w1) && ~isempty(w2),
            oy = reshape(w2(:,:,k), [d(1)*d(2),1]);
            ox = reshape(w1(:,:,k), [d(1)*d(2),1]).*oy;
        elseif nargin>=3 && ~isempty(w1),
            ox = reshape(w1(:,:,k), [d(1)*d(2),1]);
            oy = ox;
        elseif nargin>=4 && ~isempty(w2),
            ox = reshape(w2(:,:,k), [d(1)*d(2),1]);
        end
        xk(:,1) = xk(:,1).*ox;
        xk(:,2) = xk(:,2).*ox;
        xk(:,3) = xk(:,3).*ox;
    end
    yk  = reshape(y(:,:,k,:),[d(1)*d(2),3]);
    msk = find(all(isfinite(xk),2) & all(isfinite(yk),2));
    X   = [xk(msk,:), ox(msk)];
    Y   = [yk(msk,:), oy(msk)];
    XX  = XX + double(X'*X);
    XY  = XY + double(X'*Y);
end
M = (XX\XY)';

if nargout>1,
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


