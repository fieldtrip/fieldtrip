function [dat, baseline] = preproc_baselinecorrect(dat, begsample, endsample)

% PREPROC_BASELINECORRECT performs a baseline correction, e.g. using the
% prestimulus interval of the data or using the complete data
%
% Use as
%   [dat] = preproc_baselinecorrect(dat, begin, end)
% where
%   dat        data matrix (Nchans X Ntime)
%   begsample  index of the begin sample for the baseline estimate
%   endsample  index of the end sample for the baseline estimate
%
% If no begin and end sample are specified for the baseline estimate, it
% will be estimated on the complete data.
%
% See also PREPROC

% Copyright (C) 1998-2008, Robert Oostenveld
%
% $Log: preproc_baselinecorrect.m,v $
% Revision 1.2  2008/05/23 09:13:58  roboos
% cleaned up code and documentation, ensure that all functions are consistent, added proper implementation to the scratch functions
%
% Revision 1.1  2008/05/23 06:54:21  roboos
% created initial scratch version of preprocessing module, to be used in fieldtrip or as stand-alone toolbox (e.g. in spm8 or braingain)
% some functions are copies of existing roboos/misc versions, some just contain some example code for the implementation
%
% Revision 1.3  2003/03/14 10:17:28  roberto
% fixed bug that was introduced by last change, changend from repmat to for-loop
%
% Revision 1.2  2003/03/13 16:44:45  roberto
% fixed bug with multiple epochs and single channel data
%

% determine the size of the data
[Nchans, Nsamples] = size(dat);

% determine the interval to use for baseline correction
if nargin<2 || isempty(begsample)
  begsample = 1;
end
if nargin<3 || isempty(endsample)
  endsample = Nsamples;
end

% estimate the baseline and subtract it
baseline = mean(dat(:,begsample:endsample), 2);

% it is faster to loop over samples than over channels due to the internal memory representation of Matlab
% for chan=1:Nchans
%  dat(chan,:) = dat(chan,:) - baseline(chan);
% end

for sample=1:Nsamples
  dat(:,sample) = dat(:,sample) - baseline;
end


