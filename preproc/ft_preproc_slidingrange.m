function y = ft_preproc_slidingrange(dat, width, varargin)

% FT_PREPROC_SLIDINGRANGE
%
% Use as
%   y = ft_preproc_slidingrange(dat, width, ...)

normalize = keyval('normalize', varargin); if isempty(normalize), normalize = 'no'; end
normalize = istrue(normalize); % convert yes/no string into boolean

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

if normalize
    y = y ./ sqrt(width);
end
