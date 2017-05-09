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
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', {'raw', 'raw+comp'}, 'feedback', 'no');
end

% set the defaults
cfg.appendsens = ft_getopt(cfg, 'appendsens', 'no');
cfg.appenddim  = ft_getopt(cfg, 'appenddim', []);
cfg.tolerance  = ft_getopt(cfg, 'tolerance', 1e-5);

isequaltime  = true;
isequallabel = true;
issamelabel  = true;
for i=2:numel(varargin)
  isequaltime  = isequaltime  && isequal(varargin{i}.time , varargin{1}.time );
  isequallabel = isequallabel && isequal(varargin{i}.label, varargin{1}.label);
  issamelabel  = issamelabel  && numel(intersect(varargin{i}.label, varargin{1}.label))==numel(varargin{1}.label);
end

if isempty(cfg.appenddim) || strcmp(cfg.appenddim, 'auto')
  if isequallabel || issamelabel
    cfg.appenddim = 'rpt';
  elseif isequaltime && ~isequallabel && ~issamelabel
    cfg.appenddim = 'chan';
  else
    error('cannot append this data');
  end
end
fprintf('concatenating over the "%s" dimension\n', cfg.appenddim);

% ft_selectdata cannot create the union of the data contained in cell-arrays
% make a dummy without the actual data, but keep trialinfo/sampleinfo/grad/elec/opto
dummy = cell(size(varargin));
for i=1:numel(varargin)
  dummy{i} = removefields(varargin{i}, {'trial', 'time'});
  % add a dummy data field, this cause the datatype to become 'chan'
  dummy{i}.dummy       = ones(numel(dummy{i}.label),1);
  dummy{i}.dummydimord = 'chan'; 
end

% don't do any data appending inside the common function
cfg.parameter = {};
% use a low-level function that is shared with the other ft_appendxxx functions
% this will handle the trialinfo/sampleinfo/grad/elec/opto
data = append_common(cfg, dummy{:});
% this is the actual data field that will be appended further down
cfg.parameter = {'trial'};

switch cfg.appenddim
  case 'chan'
    ntrl = length(varargin{1}.trial);
    dat = varargin{1}.trial;
    lab = varargin{1}.label(:);
    for i=2:numel(varargin)
      for j=1:ntrl
        dat{j} = cat(1, dat{j}, varargin{i}.trial{j});
      end
      lab = cat(1, lab, varargin{i}.label(:));
    end
    data.trial = dat;
    data.label = lab;
    data.time  = varargin{1}.time;
    
  case 'rpt'
    if isequallabel
      % the channels are the same and sorted in the same order
      dat = varargin{1}.trial;
      tim = varargin{1}.time;
      for i=2:numel(varargin)
        dat = cat(2, dat, varargin{i}.trial);
        tim = cat(2, tim, varargin{i}.time);
      end
      data.trial = dat;
      data.time  = tim;
      data.label = varargin{1}.label;
    else
      % the channels are the same, but they do not appear in the same order
      % sort them according to the first dataset
      dat = varargin{1}.trial;
      tim = varargin{1}.time;
      for i=2:numel(varargin)
        % find the indices of the channels in each dataset, sorted according to the first dataset
        [dum, indx] = match_str(varargin{1}.label, varargin{i}.label);
        for j=1:numel(varargin{i}.trial)
          dat = cat(2, dat, varargin{i}.trial{j}(indx,:));
          tim = cat(2, tim, varargin{i}.time{j});
        end
      end
      data.trial = dat;
      data.time  = tim;
      data.label = varargin{1}.label;
    end
    
  otherwise
    error('unsupported cfg.appenddim');
end % switch

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous varargin
ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data
