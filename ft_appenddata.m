function [data] = ft_appenddata(cfg, varargin)

% FT_APPENDDATA concatenates multiple raw data structures that have been preprocessed
% separately into a single raw data structure.
%
% Use as
%   data = ft_appenddata(cfg, data1, data2, data3, ...)
%
% The following configuration options are supported:
%   cfg.keepsampleinfo  = 'yes', 'no', 'ifmakessense' (default = 'ifmakessense')
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
% If you concatenate trials and the data originates from the same original datafile,
% the sampleinfo is consistent and you can specify cfg.keepsampleinfo='yes'. If the
% data originates from different datafiles, the sampleinfo is inconsistent and does
% not point to the same recording, hence you should specify cfg.keepsampleinfo='no'.
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
% cell-array for this particular function.
%
% See also FT_PREPROCESSING, FT_DATAYPE_RAW, FT_APPENDTIMELOCK, FT_APPENDFREQ,
% FT_APPENDSOURCE, FT_APPENDSENS

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

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% set the defaults
cfg.keepsampleinfo = ft_getopt(cfg, 'keepsampleinfo', 'ifmakessense');

try
  % although not 100% robust, this could make some users becoming aware of the issue of overlapping trials
  for i=1:numel(varargin)
    dataset{i}       = ft_findcfg(varargin{i}.cfg, 'dataset');
    hassampleinfo(i) = isfield(varargin{i}, 'sampleinfo');
  end
  if ~all(strcmp(dataset, dataset{1})) && ~strcmp(cfg.keepsampleinfo, 'no')
    ft_warning('the data has overlapping segments or originates from different recordings on disk');
    ft_warning('please consider specifying cfg.keepsampleinfo=''no''')
  end
end % try

haschantype = false;
haschanunit = false;
for i=1:length(varargin)
  % if one of them has chantype or chanunit, we want it for the others as well
  haschantype = haschantype || isfield(varargin{i}, 'chantype');
  haschanunit = haschanunit || isfield(varargin{i}, 'chanunit');
end

% ensure that the input data is valid for this function
for i=1:length(varargin)
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', {'raw', 'raw+comp'}, 'feedback', 'no', 'haschantype', haschantype, 'haschanunit', haschanunit, 'hassampleinfo', cfg.keepsampleinfo);
end

% set the defaults
cfg.appendsens = ft_getopt(cfg, 'appendsens', 'no');
cfg.appenddim  = ft_getopt(cfg, 'appenddim', []);
cfg.tolerance  = ft_getopt(cfg, 'tolerance', 1e-5); % this is passed to append_common, which passes it to ft_selectdata

isequaltime  = true;
isequallabel = true;
issamelabel  = true;
isequalfsample = true;
for i=2:numel(varargin)
  isequaltime  = isequaltime  && isequal(varargin{i}.time , varargin{1}.time );
  isequallabel = isequallabel && isequal(varargin{i}.label, varargin{1}.label);
  issamelabel  = issamelabel  && isempty(setxor(varargin{i}.label, varargin{1}.label));
  isequalfsample = isequalfsample && isfield(varargin{i},'fsample') && isfield(varargin{1},'fsample') && isequal(varargin{i}.fsample, varargin{1}.fsample);
end

if isempty(cfg.appenddim) || strcmp(cfg.appenddim, 'auto')
  if isequallabel || issamelabel
    cfg.appenddim = 'rpt';
  elseif isequaltime && ~isequallabel && ~issamelabel
    cfg.appenddim = 'chan';
  else
    ft_error('cannot append this data');
  end
end
ft_info('concatenating over the "%s" dimension\n', cfg.appenddim);

% ft_selectdata cannot create the union of the data contained in cell-arrays
% make a dummy without the actual data, but keep trialinfo/sampleinfo/grad/elec/opto
% also remove the topo/unmixing/topolabel, if present, otherwise it is not
% possible to concatenate raw and comp data. Note that in append_common the
% topo etc. is removed anyhow, when appenddim = 'chan'
dummy = cell(size(varargin));
hastopo     = false;
hasunmixing = false;
for i=1:numel(varargin)
  dummy{i}    = removefields(varargin{i}, {'trial', 'time'});
  hastopo     = isfield(dummy{i}, 'topo');
  hasunmixing = isfield(dummy{i}, 'unmixing');
  if strcmp(cfg.appenddim, 'chan')
    % this avoids a bunch of confusing warnings in append_common
    dummy{i} = removefields(dummy{i}, {'topo', 'unmixing', 'topolabel'});
  end
  % add a dummy data field, this causes the datatype to become 'chan'
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
    if isequaltime
      ntrl = length(varargin{1}.trial);
      dat = varargin{1}.trial;
      lab = varargin{1}.label(:);
      if hastopo && isfield(varargin{1}, 'topo')
        topo    = varargin{1}.topo;
        topolab = varargin{1}.topolabel(:);
      else
        topo    = eye(numel(varargin{1}.label));
        topolab = varargin{1}.label(:);
      end
      if hasunmixing && isfield(varargin{1}, 'unmixing')
        unmixing = varargin{1}.unmixing;
      else
        unmixing = eye(numel(varargin{1}.label));
      end
      
      for i=2:numel(varargin)
        for j=1:ntrl
          dat{j} = cat(1, dat{j}, varargin{i}.trial{j});
        end
        lab = cat(1, lab, varargin{i}.label(:));
        
        if hastopo && isfield(varargin{i}, 'topo')
          topo    = blkdiag(topo, varargin{i}.topo);
          topolab = cat(1, topolab, varargin{i}.topolabel);
        else
          topo    = blkdiag(topo, eye(numel(varargin{i}.label)));
          topolab = cat(1, topolab, varargin{i}.label(:));
        end
        if hasunmixing && isfield(varargin{i}, 'unmixing')
          unmixing = blkdiag(unmixing, varargin{i}.unmixing);
        else
          unmixing = blkdiag(unmixing, eye(numel(varargin{i}.label)));
        end
        
      end
      data.label = lab; % replace the one from append_common
      data.trial = dat;
      data.time  = varargin{1}.time;
      if haschantype
        chantype = varargin{1}.chantype(:);
        for i=2:numel(varargin)
          chantype = cat(1, chantype, varargin{i}.chantype(:));
        end
        data.chantype = chantype;
      end
      if haschanunit
        chanunit = varargin{1}.chanunit(:);
        for i=2:numel(varargin)
          chanunit = cat(1, chanunit, varargin{i}.chanunit(:));
        end
        data.chanunit = chanunit;
      end
      
      if isfield(varargin{1}, 'topodimord'),     data.topodimord     = varargin{1}.topodimord;     end
      if isfield(varargin{1}, 'unmixingdimord'), data.unmixingdimord = varargin{1}.unmixingdimord; end
      if hastopo,     data.topo      = topo;    end
      if hastopo,     data.topolabel = topolab; end
      if hasunmixing, data.unmixing  = unmixing; end
      
    else
      ft_error('data has different time, cannot append over channels');
    end
    
  case 'rpt'
    if isequallabel
      % the channels are the same and sorted in the same order
      dat = cell(1,0);
      tim = cell(1,0);
      for i=1:numel(varargin)
        dat = cat(2, dat, varargin{i}.trial);
        tim = cat(2, tim, varargin{i}.time);
      end
      data.trial = dat;
      data.time  = tim;
    else
      % the channels are not the same, or do not appear in the same order
      % only return the intersection of the channels
      dat = cell(1,0);
      tim = cell(1,0);
      for i=1:numel(varargin)
        % find the indices of the channels in each dataset, sorted according to the first input
        [dum, indx] = match_str(data.label, varargin{i}.label);
        for j=1:numel(varargin{i}.trial)
          dat = cat(2, dat, {varargin{i}.trial{j}(indx,:)});
          tim = cat(2, tim, {varargin{i}.time{j}});
        end
      end
      data.trial = dat;
      data.time  = tim;
      data.label;       % keep it as determined by append_common
    end
    
  otherwise
    ft_error('unsupported cfg.appenddim');
end % switch

if isequalfsample
  data.fsample = varargin{1}.fsample;
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous varargin
ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data
