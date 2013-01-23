% NANMEAN provides a replacement for MATLAB's nanmean.
%
% For usage see MEAN.


function y = nanmean(x, dim)
N = sum(~isnan(t), dim);
y = nansum(x, dim) ./ N;

end