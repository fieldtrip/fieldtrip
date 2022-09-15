function [spectrum, freqoi, timeoi] = ft_specest_hilbert(dat, time, varargin)

% FT_SPECEST_HILBERT performs a spectral estimation of data by repeatedly applying a
% bandpass filter and then doing a Hilbert transform.
%
% Use as
%   [spectrum, freqoi, timeoi] = ft_specest_hilbert(dat, time, ...)
% where the input arguments are
%   dat       = matrix of chan*sample
%   time      = vector, containing time in seconds for each sample
% and the output arguments are
%   spectrum  = matrix of nchan*nfreq*ntime of fourier coefficients
%   freqoi    = vector of frequencies in spectrum
%   timeoi    = vector of timebins in spectrum
%
% Optional arguments should be specified in key-value pairs and can include
%   timeoi     = vector, containing time points of interest (in seconds)
%   freqoi     = vector, containing frequencies (in Hz)
%   pad        = number, indicating time-length of data to be padded out to in seconds (split over pre/post; used for spectral interpolation, NOT filtering)
%   padtype    = string, indicating type of padding to be used, can be 'zero', 'mean', 'localmean', 'edge', or 'mirror' (default = 'zero')
%   width      = number or vector, width of band-pass surrounding each element of freqoi
%   bpfilttype = string, filter type, 'but', 'firws', 'fir', 'firls'
%   bpfiltord  = number or vector, filter order
%   bpfiltdir  = string, filter direction, 'onepass', 'onepass-reverse', 'onepass-zerophase', 'onepass-reverse-zerophase', 'onepass-minphase', 'twopass', 'twopass-reverse', 'twopass-average'
%   edgeartnan = 0 (default) or 1, replace edge artifacts due to filtering with NaNs (only applicable for bpfilttype = 'fir'/'firls'/'firws')
%   polyorder  = number, the order of the polynomial to fitted to and removed from the data prior to the fourier transform (default = 0 -> remove DC-component)
%   verbose    = output progress to console (0 or 1, default 1)
%
% See also FT_FREQANALYSIS, FT_SPECEST_MTMFFT, FT_SPECEST_TFR, FT_SPECEST_MTMCONVOL, FT_SPECEST_WAVELET

% Copyright (C) 2010-2022, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% get the optional input arguments
freqoi     = ft_getopt(varargin, 'freqoi');
timeoi     = ft_getopt(varargin, 'timeoi', 'all');
width      = ft_getopt(varargin, 'width', 1);
bpfilttype = ft_getopt(varargin, 'bpfilttype');
bpfiltord  = ft_getopt(varargin, 'bpfiltord');
bpfiltdir  = ft_getopt(varargin, 'bpfiltdir');
bpinstabilityfix = ft_getopt(varargin, 'bpinstabilityfix', 'no');
bpfiltdf         = ft_getopt(varargin, 'bpfiltdf');
bpfiltwintype    = ft_getopt(varargin, 'bpfiltwintype');
bpfiltdev        = ft_getopt(varargin, 'bpfiltdev');

pad        = ft_getopt(varargin, 'pad');
padtype    = ft_getopt(varargin, 'padtype',    'zero');
polyorder  = ft_getopt(varargin, 'polyorder',  0);
fbopt      = ft_getopt(varargin, 'feedback');
verbose    = ft_getopt(varargin, 'verbose',    true);
edgeartnan = ft_getopt(varargin, 'edgeartnan', false);

% these should be Boolean
verbose    = istrue(verbose);
edgeartnan = istrue(edgeartnan);

if isempty(fbopt)
  fbopt.i = 1;
  fbopt.n = 1;
end

% set the default filter type
if isempty(bpfilttype)
  bpfilttype = 'but';
end

% set the default filter direction
if isempty(bpfiltdir)
  if strcmp(bpfilttype, 'firws')
    bpfiltdir = 'onepass-zerophase';
  else
    bpfiltdir = 'twopass';
  end
end

% edge artifacts from filters are only exactly defined for FIR filters
if edgeartnan && ~any(strcmp(bpfilttype, {'firws','fir','firls'}))
  ft_error('edge artifacts are only exactly defined, and can only be replaced by NaNs, for FIR filters');
elseif edgeartnan
  if isempty(bpfiltord)
    ft_error('the filter order for the FIR filters needs to be specified for edge artifact replacement');
  end  
end

% Set n's
[nchan, ndatsample] = size(dat);

% Determine fsample and set total time-length of data
fsample = 1./mean(diff(time));
dattime = ndatsample./fsample; % total time in seconds of input data

% Set a default sampling for the frequencies-of-interest if not defined
if isempty(freqoi)
  freqoi = linspace(2*width, (fsample/3), 50);
end
% check for freqoi = 0 and remove it
if any(freqoi==0)
  freqoi(freqoi==0) = [];
end
nfreqoi = length(freqoi);

% expand width to array if constant width
if numel(width) == 1
  width = ones(1,nfreqoi) * width;
end

% This does not work on integer data
if ~isa(dat, 'double') && ~isa(dat, 'single')
  dat = cast(dat, 'double');
end

% Remove polynomial fit from the data -> default is demeaning
if polyorder >= 0
  dat = ft_preproc_polyremoval(dat, polyorder, 1, ndatsample);
end

% Zero padding
if round(pad * fsample) < ndatsample
  ft_error('the padding that you specified is shorter than the data');
end
if isempty(pad) % if no padding is specified padding is equal to current data length
  pad = dattime;
end
prepad     = floor(((pad - dattime) * fsample) ./ 2);
postpad    = ceil(((pad - dattime) * fsample) ./ 2);

% Set timeboi and timeoi
timeoiinput = timeoi;
offset = round(time(1)*fsample);
if isnumeric(timeoi) % if input is a vector
  timeoi   = unique(round(timeoi .* fsample) ./ fsample);
  timeboi  = round(timeoi .* fsample - offset) + 1;
  ntimeboi = length(timeboi);
elseif strcmp(timeoi,'all') % if input was 'all'
  timeboi  = 1:length(time);
  ntimeboi = length(timeboi);
  timeoi   = time;
end

% throw a warning if input timeoi is different from output timeoi
if isnumeric(timeoiinput)
  if numel(timeoiinput) ~= numel(timeoi) % timeoi will not contain double time-bins when requested
    ft_warning('output time-bins are different from input time-bins, multiples of the same bin were requested but not given');
  else
    if any(abs(timeoiinput-timeoi) >= eps*1e6)
      ft_warning('output time-bins are different from input time-bins');
    end
  end
end

% expand filter order to array if constant filterorder
if numel(bpfiltord) == 1
  bpfiltord = ones(1,nfreqoi) * bpfiltord;
end

% create filter frequencies and check validity
filtfreq = [];
invalidind = [];
for ifreqoi = 1:nfreqoi
  tmpfreq = [freqoi(ifreqoi)-width(ifreqoi) freqoi(ifreqoi)+width(ifreqoi)];
  if all((sign(tmpfreq) == 1))
    filtfreq(end+1,:) = tmpfreq;
  else
    invalidind = [invalidind ifreqoi];
    ft_warning(sprintf('frequency %.2f Hz cannot be estimated with resolution %.2f Hz', freqoi(ifreqoi), width(ifreqoi)));
  end
end

freqoi(invalidind) = [];
nfreqoi = size(filtfreq,1);

% preallocate the result and perform the transform
spectrum = complex(nan(nchan, nfreqoi, ntimeboi), nan(nchan, nfreqoi, ntimeboi));
for ifreqoi = 1:nfreqoi
  str = sprintf('frequency %d (%.2f Hz)', ifreqoi,freqoi(ifreqoi));
  [st, cws] = dbstack;
  if length(st)>1 && strcmp(st(2).name, 'ft_freqanalysis') && verbose
    % specest_convol has been called by ft_freqanalysis, meaning that ft_progress has been initialised
    ft_progress(fbopt.i./fbopt.n, ['trial %d, ',str,'\n'], fbopt.i);
  elseif verbose
    fprintf([str, '\n']);
  end
  
  % filter
  if isempty(bpfiltord)
    ibpfiltord = [];
  else
    ibpfiltord = bpfiltord(ifreqoi);
  end
  flt = ft_preproc_bandpassfilter(dat, fsample, filtfreq(ifreqoi,:), ibpfiltord, bpfilttype, bpfiltdir, bpinstabilityfix, bpfiltdf, bpfiltwintype, bpfiltdev);
  
  % transform
  dum = transpose(hilbert(transpose(ft_preproc_padding(flt, padtype, prepad, postpad))));
  
  if ~edgeartnan % do not replace edge artifacts by NaNs (or keep them as NaNs)
    % select output using timeboi
    spectrum(:,ifreqoi,:) = dum(:,timeboi+prepad);
    
  else % replace edge artifacts by NaNs(i.e. keep them as NaNs
    
    % determine samples to be NaN in output
    % for fir/firws the filter order is always even, so an equal amount of NaNs on either
    % side, but for firls it obtuse in the code. To be safe, ceil(ord/2) is used for NaNing purpose
    % when applicable
    switch bpfiltdir
      case {'onepass','onepass-minphase'} % artifact is at the start as long as full order
        edgeartind = 1:bpfiltord(ifreqoi); % the low level filtering functions enforce the order to be smaller than the data length
      case 'onepass-reverse' % artifact is at the end as long as full order
        edgeartind = (ndatsample-bpfiltord(ifreqoi)+1):ndatsample;
      case {'twopass', 'twopass-reverse', 'twopass-average'} % artifact is at the start and end as long as full order
        edgeartind = [1:bpfiltord(ifreqoi) (ndatsample-bpfiltord(ifreqoi)+1):ndatsample];
      case {'onepass-zerophase','onepass-reverse-zerophase'} % artifact is at the start and end as long as half order
        halforder = ceil(bpfiltord(ifreqoi)/2);
        edgeartind = [1:halforder (ndatsample-halforder+1):ndatsample];
    end
    
    % select output using timeboi and ignore those bins that are part of the edge artifacts (i.e. identically as to specest_mtmconvol)
    reqtimeboi = setdiff(timeboi,edgeartind);
    spectrum(:,ifreqoi,reqtimeboi) = dum(:,reqtimeboi+prepad);
  end  
end







