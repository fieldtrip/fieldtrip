function y = ft_preproc_slidingrange(dat, width, varargin)

% FT_PREPROC_SLIDINGRANGE computes the range of the data in a sliding time
% window of the width specified. Width should be an odd number (since the
% window needs to be centered on an individual sample).
%
% Use as
%   y = ft_preproc_slidingrange(dat, width, ...)
%
% Optional key-value pair arguments are:
%   'normalize', whether to normalize the range of the data with the square
%                root of the window size

normalize = ft_getopt(varargin, 'normalize', false);

if mod(width+1, 2)
  error('width should be an odd number');
end

% compute half width
h = (width-1)/2;

n = size(dat,2);
minval = zeros(size(dat));
maxval = zeros(size(dat));

for i=1:n
  begsample = i-h;
  endsample = i+h;
  if begsample<1
    begsample = 1;
  end
  if endsample>n
    endsample=n;
  end
  minval(:,i) = min(dat(:,begsample:endsample),[],2);
  maxval(:,i) = max(dat(:,begsample:endsample),[],2);
end

y = maxval - minval;

if istrue(normalize)
  y = y ./ sqrt(width);
end
