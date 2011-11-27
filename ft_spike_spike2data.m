function [spikedata] = ft_spike_spike2data(cfg,spike,data)

% FT_SPIKE_SPIKE2DATA converts a point representation SPIKE
% structure to a continuous representation DATA structure.
%
% Use as
%   [spikedata] = ft_spike_spike2data(cfg, spike, data)
% or
%   [spikedata] = ft_spike_spike2data(cfg, spike)
%
% Inputs:
%   Data is a standard continuous structure. DATA is an optionary input. If
%   data is given as input, the output structure will be using the
%   time-axis as given in data, where each spike is assigned a sample on
%   that time-axis. Note that in this case, it is very important that for
%   both data and spike, the trial definitions should be synchronized,
%   i.e., t = 0 has the same meaning, and there are equal amount of trials.
%
% Configurations:
%   cfg.fsample = number in Hz. If DATA.fsample is given as input, cfg.fsample is
%                 not used and set to DATA.fsample. Otherwise, the default is 1000 Hz
%   cfg.sparse    = 'yes' (default) or 'no'
%
%  Outputs:
%   spikedata is a continuous representation standard spike
%   structure. If requested, spikedata.trial contains sparse matrices.
% 
% Appending should be done with FT_APPENDDATA

% Copyright (C) 2010, Martin Vinck
%
% $Id$

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig

% control input spike structure
spike = ft_checkdata(spike,'datatype', 'spike', 'feedback', 'yes');

% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', {'fsample'});
cfg.sparse = ft_getopt(cfg,'sparse', 'yes');

% ensure that the options are valid
cfg = ft_checkopt(cfg,'fsample','doublescalar');
cfg = ft_checkopt(cfg,'sparse', 'char', {'yes', 'no'});

% check whether the data has the right structure
hasAllFields = all(isfield(spike, {'time', 'trial', 'trialtime', 'label'}));
if ~hasAllFields, error('MATLAB:ft_spike_spike2data:wrongStructInput',...
    'input SPIKE should be struct with .time, .trial, .trialtime, .label fields')
end

% check whether all are of right format
correctInp = iscell(spike.time) & iscell(spike.trial) & iscell(spike.label) & size(spike.trialtime,2)==2;
if ~correctInp, error('MATLAB:ft_spike_spike2data:wrongStructInput',...
    '.timestamp, .trial and .label should be cell arrays, time should be nTrials-by-2 matrix')
end

% check whether the data should be appended
doSparse = strcmp(cfg.sparse, 'yes');

% get some sizes
nUnits  = length(spike.label);
nTrials = size(spike.trialtime,1);

% check whether number of trials of DATA matches that of SPIKE: it should
if nargin==3
  nTrialsData = length(data.trial);
  if nTrials~=nTrialsData
    error('MATLAB:ft_spike_spike2data:numberOfTrials',...
      'number trials SPIKE should match number trials of input DATA')
  end
end

% preallocate
spikedata.trial(1:nTrials) = {[]};
spikedata.time(1:nTrials)  = {[]};
for iTrial = 1:nTrials
  
  % create the time-axis that belongs to this trial
  if nargin==3
    timeAx   = data.time{iTrial};
  else
    timeAx   = spike.trialtime(iTrial,1):(1/cfg.fsample):spike.trialtime(iTrial,2);
  end
  
  % convert to continuous
  trialData = zeros(nUnits,length(timeAx));
  for iUnit = 1:nUnits
    
    % get the timestamps and only select those timestamps that are in the trial
    ts       = spike.time{iUnit};
    hasTrial = spike.trial{iUnit}==iTrial;
    ts       = ts(hasTrial);
    
    % get all the samples at once without using loops
    sample   = nearest_nd(timeAx,ts);
    
    % because we have duplicates, simply get our vector by using histc trick
    [N] = histc(sample,1:length(timeAx));
    
    % store it in a matrix
    trialData(iUnit,:) = N;
  end
  
  % put the created data in the original DATA, which is empty if doAppend==0
  if doSparse, trialData = sparse(trialData); end
  spikedata.trial{iTrial} = [spikedata.trial{iTrial};trialData];
  spikedata.time{iTrial} = timeAx;
  
end

% create the associated labels and other aspects of data such as the header
spikedata.label = spike.label;
if nargin~=3
  spikedata.fsample    = cfg.fsample;
end
try
  spikedata.hdr     = spike.hdr;
catch
  spikedata.hdr     = struct([]);
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
if nargin==3
  ft_postamble previous spike data
elseif nargin==2
  ft_postamble previous spike
end
ft_postamble history spikedata


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [indx] = nearest_nd(x,y)

% NEAREST return the index of an n-d matrix to an n-d matrix.
%
% [indx] = nearest_nd(x, y)
%
% Inputs:
%   X can be a n-d matrix of any size (scalar, vector, n-d matrix).
%   Y can be an n-d matrix of any size (scalar, vector, n-d matrix).
%
% If Y is larger than any X, we return the last index that the maximum value of X occurred.
% Otherwise, we return the first occurence of the nearest X.
%
% If Y contains NaNs, we return a NaN for every NaN in Y.
%
% Outputs:
%   INDX is a vector of size Y and contains the indices of the values in X that are
%   closest to the respective value of Y. INDX is a linear index, such that x(INDX) gives
%   the nearest values of X to Y. To convert INDX to subscripts, see IND2SUB.
%
% Example:
% y = 100*rand(5,10);
% x = 1:100;
% indx = nearest_nd(x,y);
% disp(max(max(y-x(indx))))
%
% x = reshape([1:100 100*ones(1,20)],[20 6]);
% y = rand(5,10)*max(x(:));
% y(1) = 100.2;
% indx = nearest_nd(x,y);
% disp(max(max(y-x(indx))))
% Note that indx(1) returns the last occurence and indx is a linear index to 2-d X.
%
% Demonstrate the use with a 3-D x and y array.
% x = reshape([1:100],[10 5 2]);
% y = rand(5,10,2)*max(x(:));
% indx = nearest_nd(x,y);
% d = x(indx)-y; disp(max(d(:)));
% y(1) = NaN;
% indx = nearest_nd(x,y);
% disp(indx(1))
%

% Copyright, Martin Vinck, 2009.

% store the sizes of x and y, this is used to reshape INDX later on
szY = size(y);

% vectorize both x and y
x = x(:);
y = y(:);

% from now on we can treat X and Y as vectors
nY = length(y);
nX = length(x);
hasNan = isnan(y); % indices with nans in y, indx(hasNaN) will be set to NaN later.

if nX==1,
  indx = ones(1,nY);  % only one x value, so nearest is always only element
else
  if nY==1 % in this case only one y value, so use the old NEAREST code
    if y>max(x)
      % return the last occurence of the nearest number
      [dum, indx] = max(flipud(x));
      indx = nX + 1 - indx;
    else
      % return the first occurence of the nearest number
      [mindist, indx] = min(abs(x(:) - y));
    end
  else
    if any(y>max(x))
      % for these return the last occurence of every number as in NEAREST
      indx = zeros(1,nY);
      i = y>max(x);
      [dum,indx] = max(flipud(x));
      indx(i)       = nX + 1 - indx;
      % for the rest return the first occurence of every number
      x = x(:);
      y = y(~i)';
      xRep = x(:,ones(1,length(y)));
      yRep = y(ones(nX,1),:);
      [mindist,indx(~i)] = min(abs(xRep-yRep));
    else
      x = x(:);
      y = y';
      xRep = x(:,ones(1,nY));
      yRep = y(ones(nX,1),:);
      [mindist,indx] = min(abs(xRep-yRep));
    end
  end
end
% return a NaN in INDX for a NaN in Y
indx(hasNan) = NaN;

% reshape the indx back to the y format
if (sum(szY>1)>1 || length(szY)>2) % in this case we are dealing with a matrix
  indx = reshape(indx,[szY]);
end

