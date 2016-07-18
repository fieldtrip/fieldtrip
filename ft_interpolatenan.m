function dataout = ft_interpolatenan(cfg, datain)

% FT_INTERPOLATENAN interpolates data that contains segments of nans
% obtained by replacing artifactual data with nans using, for example, FT_REJECTARTIFACT,
% or by redefining trials with FT_REDEFINETRIAL resulting in trials with
% gaps.
%
% Use as
%  outdata = ft_interpolatenan(cfg, indata)
% where indata is data as obtained from FT_PREPROCESSING
% and cfg is a configuratioun structure that should contain
%
%  cfg.method      = string, interpolation method, see HELP INTERP1 (default = 'linear')
%  cfg.prewindow   = value, length of data prior to interpolation window, in seconds (default = 1)
%  cfg.postwindow  = value, length of data after interpolation window, in seconds (default = 1)
%
%
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following options:
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_REJECTARTIFACT, FT_REDEFINETRIAL

% Copyright (C) 2003-2011, Jan-Mathijs Schoffelen & Robert Oostenveld
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
ft_preamble loadvar datain
ft_preamble provenance datain
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
datain = ft_checkdata(datain, 'datatype', {'raw+comp', 'raw'}, 'feedback', 'yes', 'hassampleinfo', 'yes');

% check if the input is valid
cfg             = ft_checkconfig(cfg, 'allowedval',{'method','nearest','linear','spline','pchip','cubic','v5cubic'});
cfg             = ft_checkopt(cfg, 'prewindow', 'numericscalar');
cfg             = ft_checkopt(cfg, 'postwindow', 'numericscalar');

% get the options
method          = ft_getopt(cfg, 'method', 'linear');        % default is linear
prewindow       = ft_getopt(cfg, 'prewindow', 1);        % default is 1 second
postwindow      = ft_getopt(cfg, 'postwindow', 1);        % default is 1 seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataout = datain;
ntrl = numel(datain.trial);
prewindow = round(prewindow * datain.fsample); % Express window in samples
postwindow = round(postwindow * datain.fsample); % Express window in samples

% Let users know that the interpolation will start and initialize the
% progress indicator
fprintf('Initializing %s interpolation of %d trials\n',method,ntrl);
ft_progress('init', 'etf');

for i=1:ntrl
  ft_progress(i/ntrl, 'Processing trial %d out of %d',i,ntrl);
  nchan = size(datain.trial{i},1);
  replace = isnan(datain.trial{i}); % Find samples that have been replaced by nans
  if any(replace(:)) % Check whether any values should be interpolated
    [idx_start_r, idx_start_c] = find(diff(replace,[],2)==1); % Determine onset of nan-chunk
    [dum, idx_end_c] = find(diff(replace,[],2)==-1); % Determine offset of nan-chunk
    idx_start_c = idx_start_c + 1; % correct for shift due to using diff

    % Loop across found nan-chunks
    for j=1:size(idx_start_c,1)
     sample_window = [idx_start_c(j)-prewindow:idx_start_c(j)-1 idx_end_c(j)+1:idx_end_c(j)+postwindow]; % Indices of time-points used for interpolation
     if any(sample_window<1) || any(sample_window>size(datain.trial{i},2)) % Check whether sampling window falls within data range
       ft_warning('Sample window partially outside of data-range, using less samples');
       sample_window(sample_window<1|sample_window>size(datain.trial{i},2))=[];
     elseif any(isnan(datain.trial{i}(idx_start_r(j),sample_window))) % Check whether sampling window overlaps with other chunk of nans
      error('Sample window overlaps with other chunk of nans');
     end;
     fill_window = idx_start_c(j):idx_end_c(j); % Indices of time-points that will be interpolated
     fill = interp1(sample_window, datain.trial{i}(idx_start_r(j),sample_window),fill_window,method); % Interpolation
     dataout.trial{i}(idx_start_r(j),fill_window)=fill; % Add interpolated segments to dataout
    end;
  end;
end; % i=1:ntrl
ft_progress('close');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cleanup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   datain
ft_postamble provenance dataout
ft_postamble history    dataout
ft_postamble savevar    dataout
