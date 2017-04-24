function [data] = ft_appenddata(cfg, varargin)

% FT_APPENDDATA concatenates multiple raw data structures that have been preprocessed
% separately into a single raw data structure.
%
% Use as
%   data = ft_appenddata(cfg, data1, data2, data3, ...)
% where the configuration can be empty.
%
% If the input datasets all have the same channels, the trials will be concatenated.
% This is useful for example if you have different experimental conditions, which,
% besides analyzing them separately, for some reason you also want to analyze
% together. The function will check for consistency in the order of the channels. If
% the order is inconsistent the channel order of the output will be according to the
% channel order of the first data structure in the input.
%
% If the input datasets have different channels, but the same number of trials, the
% channels will be concatenated within each trial. This is useful for example if the
% data that you want to analyze contains both MEG and EMG channels which require
% different preprocessing options.
%
% Occasionally, the data needs to be concatenated in the trial dimension while
% there's a slight discrepancy in the channels in the input data (e.g. missing
% channels in one of the data structures). The function will then return a data
% structure containing only the channels which are present in all inputs.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure. The data structure in the input file should be a
% cell array for this particular function.
%
% See also FT_PREPROCESSING, FT_DATAYPE_RAW, FT_APPENDTIMELOCK, FT_APPENDFREQ,
% FT_APPENDSENS, FT_APPENDSOURCE

% Copyright (C) 2005-2008, Robert Oostenveld
% Copyright (C) 2009-2011, Jan-Mathijs Schoffelen
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    varargin
ft_preamble provenance varargin
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
for i=1:length(varargin)
  % FIXME: raw+comp is not always dealt with correctly
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', {'raw+comp', 'raw'}, 'feedback', 'no');
end

% set the defaults
cfg.appendsens = ft_getopt(cfg, 'appendsens', 'no');
cfg.appenddim  = ft_getopt(cfg, 'appenddim', []);
cfg.parameter  = {'trial'}; % this is hard-coded, it is used for consistency with ft_appendtimelock and ft_appendfreq

if isempty(cfg.appenddim) || strcmp(cfg.appenddim, 'auto')
  if checkchan(varargin{:}, 'identical') && checktime(varargin{:}, 'identical', cfg.tolerance)
    cfg.appenddim = 'rpt';
  elseif checkchan(varargin{:}, 'unique')
    cfg.appenddim = 'chan';
%   elseif checktime(varargin{:}, 'unique', cfg.tolerance)
%     cfg.appenddim = 'time';
  else
    error('cfg.appenddim should be specified');
  end
end
fprintf('concatenating over the "%s" dimension\n', cfg.appenddim);

% use a low-level function that is shared with the other ft_appendxxx functions
data = append_common(cfg, varargin{:});

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous varargin
ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data
