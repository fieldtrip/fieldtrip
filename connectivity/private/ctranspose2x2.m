function out = ctranspose2x2(in)

% compute ctranspose of multiple 2x2 matrices, input is 2x2xN

out = conj(in);
out(1,2,:,:) = conj(in(2,1,:,:));
out(2,1,:,:) = conj(in(1,2,:,:));
