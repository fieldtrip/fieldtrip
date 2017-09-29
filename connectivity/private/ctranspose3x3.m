function out = ctranspose3x3(in)

% compute ctranspose of multiple 3x3 matrices, input is 3x3xN

out = conj(in);
out(1,2,:,:) = conj(in(2,1,:,:));
out(2,1,:,:) = conj(in(1,2,:,:));
out(1,3,:,:) = conj(in(3,1,:,:));
out(3,1,:,:) = conj(in(1,3,:,:));
out(2,3,:,:) = conj(in(3,2,:,:));
out(3,2,:,:) = conj(in(2,3,:,:));
