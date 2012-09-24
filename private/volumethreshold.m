function [output] = volumethreshold(input, thresh, str)

% VOLUMETHRESHOLD is a helper function for segmentations

fprintf('thresholding %s at a relative threshold of %0.3f\n', str, thresh);

% mask by taking the negative of the brain, thus ensuring
% that no holes are within the compartment and do a two-pass
% approach to eliminate potential vitamin E capsules etc.

output   = double(input>(thresh*max(input(:))));
[tmp, N] = spm_bwlabel(output, 6);
for k = 1:N
  n(k,1) = sum(tmp(:)==k);
end
output   = double(tmp~=find(n==max(n))); clear tmp;
[tmp, N] = spm_bwlabel(output, 6);
for k = 1:N
  m(k,1) = sum(tmp(:)==k);
end
output   = double(tmp~=find(m==max(m))); clear tmp;