function [A, b, newpos] = regpt2surf(node, elem, p, pmask, A0, b0, cmask, maxiter)
%  [A,b,newpos]=regpt2surf(node,elem,p,pmask,A0,b0,cmask,maxiter)
%  Perform point cloud registration to a triangular surface
%  (surface can be either triangular or cubic), Gauss-Newton method
%  is used for the calculation
%
%  author: Qianqian Fang <q.fang at neu.edu>
%  date: 12/12/2008
%
% parameters:
%      node: node coordinate of the surface mesh (nn x 3)
%      elem: element list of the surface mesh (3 columns for
%            triangular mesh, 4 columns for cubic surface mesh)
%      p: points to be registered, 3 columns for x,y and z respectively
%      pmask: a mask vector with the same length as p, determines the
%         method to handle the point, if pmask(i)=-1, the point is a free
%         node and can be move by the optimization, if pmask(i)=0, the
%         point is fixed; if pmask(i)=n>0, the distance between p(i,:)
%         and node(n,:) will be part of the object function and be optimized
%      A0: a 3x3 matrix, as the initial guess for the affine A matrix (rotation&scaling)
%      b0: a 3x1 vector, as the initial guess for the affine b vector (translation)
%      cmask: a binary 12x1 vector, determines which element of [A(:);b] will be optimized
%          if cmask(i)=0, the corresponding coefficient will not be updated
%      maxiter: a integer, specifying the optimization iterations
%
% outputs:
%      A: 3x3 matrix, the updated affine A matrix
%      b: 3x1 vector, the updated affine b vector
%      newpos: the registered positions for p, newpos=A*p'+b
%
% Please find more information at http://iso2mesh.sf.net/cgi-bin/index.cgi?metch
%
% this function is part of "metch" toobox, see COPYING for license

A = A0;
b = b0(:);

% for simplicity, I wrap A and b into one single vector C
C = [A(:); b];

delta = 1e-4;
newpos = (reshape(C(1:9), 3, 3) * p' + repmat(C(end - 2:end), 1, size(p, 1)))';
nv = nodesurfnorm(node, elem);
clen = length(C);
cuplist = find(cmask == 1);
pfree = find(pmask < 0);
pfix = find(pmask > 0);

% start Gauss-Newton iterations
for iter = 1:maxiter

    % calculate the current residual: the sum of distances to the surface
    dist0 = zeros(length(pfree) + length(pfix), 1);
    if (~isempty(pfree))
        [dist0(pfree), cn0] = dist2surf(node, nv, newpos(pfree, :));
    end
    if (~isempty(pfix))
        fixdist = node(pmask(pfix), :) - newpos(pfix, :);
        dist0(pfix) = sqrt(sum(fixdist .* fixdist, 2));
    end
    fprintf('iter=%d error=%f\n', iter, sum(abs(dist0)));

    % build the Jacobian (sensitivity) matrix
    J = zeros(length(dist0), length(C));
    for i = 1:clen
        if (cmask(i) == 0)
            continue
        end
        dC = C;
        if (C(i))
            dC(i) = C(i) * (1 + delta);
        else
            dC(i) = C(i) + delta;
        end
        newp = (reshape(dC(1:9), 3, 3) * p' + repmat(dC(end - 2:end), 1, size(p, 1)))';

        dist = zeros(length(pfree) + length(pfix), 1);
        if (~isempty(pfree))
            if (length(cn0) == length(pfree))
                dist(pfree) = dist2surf(node, nv, newp(pfree, :), cn0);
            else
                dist(pfree) = dist2surf(node, nv, newp(pfree, :));
            end
        end
        if (~isempty(pfix))
            fixdist = node(pmask(pfix), :) - newp(pfix, :);
            dist(pfix) = sqrt(sum(fixdist .* fixdist, 2));
        end
        % J=dL/dC
        J(:, i) = (dist - dist0) / (dC(i) - C(i));
    end

    % weight the matrix (normalization)
    wj = sqrt(sum(J .* J));
    J(:, cuplist) = J(:, cuplist) ./ repmat(wj(cuplist), length(dist0), 1);

    % calculate the update: J*dC=dL
    dC = (J(:, cuplist) \ dist0) ./ wj(cuplist)';
    C(cuplist) = C(cuplist) - 0.5 * dC;

    % get the updated positions with the calculated A and b
    newpos = (reshape(C(1:9), 3, 3) * p' + repmat(C(end - 2:end), 1, size(p, 1)))';
end

A = reshape(C(1:9), 3, 3);
b = C(end - 2:end);
