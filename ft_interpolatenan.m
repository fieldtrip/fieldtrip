function [dataout] = ft_interpolatenan(cfg, datain)

% FT_INTERPOLATENAN interpolates time series that contains segments of nans obtained
% by replacing artifactual data with nans using, for example, FT_REJECTARTIFACT, or
% by redefining trials with FT_REDEFINETRIAL resulting in trials with gaps.
%
% Use as
%   outdata = ft_interpolatenan(cfg, indata)
% where cfg is a configuration structure and the input data is obtained from FT_PREPROCESSING.
%
% The configuration should contain
%   cfg.method      = string, interpolation method, see INTERP1 (default = 'linear')
%   cfg.prewindow   = value, length of data prior to interpolation window, in seconds (default = 1)
%   cfg.postwindow  = value, length of data after interpolation window, in seconds (default = 1)
%   cfg.feedback    = string, 'no', 'text', 'textbar', 'gui' (default = 'text')
%
% This function only interpolates over time, not over space. If you want to
% interpolate using spatial information, e.g. using neighbouring channels, you should
% use FT_CHANNELREPAIR.
%
% To facilitate data-handling and distributed computing, you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_REJECTARTIFACT, FT_REDEFINETRIAL, FT_CHANNELREPAIR

% Copyright (C) 2003-2020, Jan-Mathijs Schoffelen & Robert Oostenveld
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the initial part deals with parsing the input options and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    datain
ft_preamble provenance datain

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
datain = ft_checkdata(datain, 'datatype', {'raw+comp', 'raw'}, 'feedback', 'yes', 'hassampleinfo', 'yes');

% get the options
cfg.method      = ft_getopt(cfg, 'method',    'linear'); % default is linear, can be 'nearest', 'linear', 'spline', 'pchip', 'cubic', 'v5cubic', 'makima'
cfg.prewindow   = ft_getopt(cfg, 'prewindow',  1);       % default is 1 second
cfg.postwindow  = ft_getopt(cfg, 'postwindow', 1);       % default is 1 seconds
cfg.feedback    = ft_getopt(cfg, 'feedback', 'etf');

% check if the input is valid
cfg = ft_checkopt(cfg, 'prewindow', 'numericscalar');
cfg = ft_checkopt(cfg, 'postwindow', 'numericscalar');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prewindow  = round(cfg.prewindow * datain.fsample);  % Express window in samples
postwindow = round(cfg.postwindow * datain.fsample); % Express window in samples

% Start with a copy of the input data
dataout = datain;

% Let users know that the interpolation will start and initialize the progress indicator
ntrl = numel(datain.trial);
nchan = numel(datain.label);
fprintf('Initializing %s interpolation of %d trials\n', cfg.method, ntrl);
ft_progress('init', cfg.feedback, 'Processing trial...');

for i=1:ntrl
  ft_progress(i/ntrl, 'Processing trial %d from %d', i, ntrl);
  tim = datain.time{i};
  for j=1:nchan
    dat = datain.trial{i}(j,:);
    replace = isnan(dat); % Find samples that have been replaced by nans
    if ~any(replace)
      continue
    end
    onset  = find(diff([0 replace])>0);
    offset = find(diff([replace 0])<0);
    
    for k=1:numel(onset)
      begsample = onset(k)-prewindow;
      endsample = offset(k)+postwindow;
      if begsample<1
        ft_warning('not enough samples for prewindow')
        begsample = 1;
      end
      if endsample>numel(dat)
        ft_warning('not enough samples for postwindow')
        endsample = numel(dat);
      end
      
      x = tim(begsample:endsample);
      y = dat(begsample:endsample);
      xx = x; % this is where we want to know the interpolated values
      x = x(~replace(begsample:endsample)); % remove the part that needs to be interpolated
      y = y(~replace(begsample:endsample)); % remove the part that needs to be interpolated
      
      yy = interp1(x, y, xx, cfg.method); % this may contain nans
      
      % The default extrapolation behavior of INTERP1 with four input arguments is to
      % extrapolate for 'spline', 'pchip' and 'makima', and to use nan for other
      % methods.
      
      if begsample==1
        % there may be nans at the beginning, replace the data with mean of the values that are not nan
        f = find(~isnan(yy), 1, 'first');
        yy(1:f-1) = nanmean(yy);
      elseif endsample==numel(dat)
        % there may be nans at the end, replace the data with mean of the values that are not nan
        f = find(~isnan(yy), 1, 'last');
        yy(f+1:end) = nanmean(yy);
      end
      
      % insert the interpolated data
      dataout.trial{i}(j,begsample:endsample) = yy;
      
    end % for all nan-segments
  end % for all channels
end % for all trials
ft_progress('close');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cleanup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ft_postamble debug
ft_postamble previous   datain
ft_postamble provenance dataout
ft_postamble history    dataout
ft_postamble savevar    dataout
