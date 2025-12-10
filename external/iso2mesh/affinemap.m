function [A, b] = affinemap(pfrom, pto)
%  [A,b]=affinemap(pfrom,pto)
%
%  calculate an affine transform (A matrix and b vector) to map n
%  vertices from one space to the other using least square solutions
%
%  author: Qianqian Fang <q.fang at neu.edu>
%  date: 12/12/2008
%
% parameters:
%      pfrom: nx3 matrix, each row is a 3d point in original space
%      pto: nx3 matrix, each row is a 3d point in the mapped space
%
% outputs:
%      A: 3x3 matrix, the calculated affine A matrix
%      b: 3x1 vector, the calculated affine b vector
%
% the solution will satisfy the following equation: A*pfrom'+b=pto
%
% Please find more information at http://iso2mesh.sf.net/cgi-bin/index.cgi?metch
%
% this function is part of "metch" toobox, see COPYING for license

bsubmat = eye(3);
ptnum = size(pfrom, 1);
if (size(pto, 1) ~= ptnum)
    error('two inputs should have the same size');
end
amat = zeros(ptnum * 3, 9);
for i = 1:ptnum
    amat(i * 3 - 2:i * 3, :) = kron(bsubmat, pfrom(i, :));
end
amat = [amat, repmat(bsubmat, ptnum, 1)];

bvec = pto';
bvec = bvec(:);

x = amat \ bvec;
A = reshape(x(1:9), 3, 3)';
b = x(end - 2:end);
