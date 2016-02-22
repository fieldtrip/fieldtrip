function [freq] = ft_appendfreq(cfg, varargin)

% FT_APPENDFREQ concatenates multiple frequency or time-frequency data
% structures that have been processed separately. If the input data
% structures contain different channels, it will be concatenated along the
% channel direction. If the channels are identical in the input data
% structures, the data will be concatenated along the repetition dimension.
%
% Use as
%  combined = ft_appendfreq(cfg, freq1, freq2, ...)
%
%  cfg.parameter  = String. Specifies the name of the field to concatenate.
%                   For example, to concatenate freq1.powspctrm,
%                   freq2.powspctrm etc, use cft.parameter = 'powspctrm'.
%
% The configuration can optionally contain
%  cfg.appenddim  = String. The dimension to concatenate over (default is 'auto').
%  cfg.tolerance  = Double. Tolerance determines how different the units of
%                   frequency structures are allowed to be to be considered
%                   compatible (default: 1e-5).
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a
% *.mat file on disk and/or the output data will be written to a *.mat file.
% These mat files should contain only a single variable, corresponding with
% the input/output structure.
%
% See also FT_FREQANALYSIS, FT_APPENDDATA, FT_APPENDTIMELOCK, FT_APPENDSOURCE

% Copyright (C) 2011, Robert Oostenveld
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
ft_preamble loadvar varargin
ft_preamble provenance varargin
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
for i=1:length(varargin)
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'freq');
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'required', 'parameter');

% set the defaults
cfg.channel    = ft_getopt(cfg, 'channel', 'all');
cfg.appenddim  = ft_getopt(cfg, 'appenddim',  'auto');
cfg.tolerance  = ft_getopt(cfg, 'tolerance',  1e-5);

% do a basic check to see whether the dimords match
Ndata = length(varargin);
dimord = cell(1,Ndata);
for i=1:Ndata
  dimord{i} = varargin{i}.dimord;
end
dimordmatch = all(strcmp(dimord{1}, dimord));

if ~dimordmatch
  error('the dimords of the input data structures are not equal');
end

% create the output structure from scratch
freq   = [];

tol    = cfg.tolerance;
dimtok = tokenize(dimord{1}, '_');
switch cfg.appenddim
  case 'auto'

    % only allow to append across observations if these are present in the data
    if any(strcmp(dimtok, 'rpt'))
      cfg.appenddim = 'rpt';
    elseif any(strcmp(dimtok, 'rpttap'))
      cfg.appenddim = 'rpttap';
    elseif any(strcmp(dimtok, 'subj'))
      cfg.appenddim = 'subj';
    else
      % we need to check whether the other dimensions are the same.
      % if not, consider some tolerance.
      boolval1 = checkchan(varargin{:}, 'identical');
      boolval2 = checkfreq(varargin{:}, 'identical', tol);

      if isfield(varargin{1}, 'time'),
        boolval3 = checktime(varargin{:}, 'identical', tol);
        if boolval1 && boolval2 && boolval3
          % each of the input datasets contains a single repetition (perhaps an average), these can be concatenated
          cfg.appenddim = 'rpt';
        elseif ~boolval1 && boolval2 && boolval3
          cfg.appenddim = 'chan';
        elseif boolval1 && ~boolval2 && boolval3
          cfg.appenddim = 'freq';
        elseif boolval1 && boolval2 && ~boolval3
          cfg.appenddim = 'time';
        else
          error('the input datasets have multiple non-identical dimensions, this function only appends one dimension at a time');
        end
      else
        if boolval1 && boolval2
          % each of the input datasets contains a single repetition (perhaps an average), these can be concatenated
          cfg.appenddim = 'rpt';
        elseif ~boolval1 && boolval2
          cfg.appenddim = 'chan';
        elseif boolval1 && ~boolval2
          cfg.appenddim = 'freq';
        end
      end

    end % determine the dimension for appending
end

switch cfg.appenddim
  case {'rpt' 'rpttap' 'subj'}
    catdim = find(ismember(dimtok, {'rpt' 'rpttap' 'subj'}));
    if numel(catdim)==0
      catdim = 0;
    elseif numel(catdim)==1
      % this is OK
    elseif numel(catdim)>1
      error('ambiguous dimord for concatenation');
    end

    % if any of these are present, concatenate
    % if not prepend the dimord with rpt (and thus shift the dimensions)

    % here we need to check whether the other dimensions are the same. if
    % not, consider some tolerance.
    boolval1 = checkchan(varargin{:}, 'identical');
    boolval2 = checkfreq(varargin{:}, 'identical', tol);
    if isfield(varargin{1}, 'time'),
      boolval3 = checktime(varargin{:}, 'identical', tol);
    else
      boolval3 = true;
    end

    if any([boolval2 boolval3]==false)
      error('appending across observations is not possible, because the frequency and/or temporal dimensions are incompatible');
    end

    % select and reorder the channels that are in every dataset
    tmpcfg           = [];
    tmpcfg.channel   = cfg.channel;
    tmpcfg.tolerance = cfg.tolerance;
    [varargin{:}]    = ft_selectdata(tmpcfg, varargin{:});
    for i=1:Ndata
      [cfg_rolledback, varargin{i}] = rollback_provenance(cfg, varargin{i});
    end
    cfg = cfg_rolledback;

    % update the dimord
    if catdim==0
      freq.dimord = ['rpt_',varargin{1}.dimord];
      % FIXME append dof
    else
      freq.dimord = varargin{1}.dimord;
      % FIXME append dof
      % before append cumtapcnt cumsumcnt trialinfo check if there's a
      % subfield in each dataset. Append fields of different dataset might
      % lead in empty and/or non-existing fields in a particular dataset
      hascumsumcnt = [];
      hascumtapcnt = [];
      hastrialinfo = [];
      for i=1:Ndata
        if isfield(varargin{i},'cumsumcnt');
          hascumsumcnt(end+1) = 1;
        else
          hascumsumcnt(end+1) = 0;
        end
        if isfield(varargin{i},'cumtapcnt');
          hascumtapcnt(end+1) = 1;
        else
          hascumtapcnt(end+1) = 0;
        end
        if isfield(varargin{i},'trialinfo');
          hastrialinfo(end+1) = 1;
        else
          hastrialinfo(end+1) = 0;
        end
      end

      % screen concatenable fields
      if ~checkfreq(varargin{:}, 'identical', tol)
        error('the freq fields of the input data structures are not equal');
      else
        freq.freq=varargin{1}.freq;
      end
      if ~sum(hascumsumcnt)==0 && ~(sum(hascumsumcnt)==Ndata);
        error('the cumsumcnt fields of the input data structures are not equal');
      else
        iscumsumcnt=unique(hascumsumcnt);
      end
      if ~sum(hascumtapcnt)==0 && ~(sum(hascumtapcnt)==Ndata);
        error('the cumtapcnt fields of the input data structures are not equal');
      else
        iscumtapcnt=unique(hascumtapcnt);
      end
      if ~sum(hastrialinfo)==0 && ~(sum(hastrialinfo)==Ndata);
        error('the trialinfo fields of the input data structures are not equal');
      else
        istrialinfo=unique(hastrialinfo);
      end

      % concatenating fields
      for i=1:Ndata;
        if iscumsumcnt;
          cumsumcnt{i}=varargin{i}.cumsumcnt;
        end
        if iscumtapcnt;
          cumtapcnt{i}=varargin{i}.cumtapcnt;
        end
        if istrialinfo;
          trialinfo{i}=varargin{i}.trialinfo;
        end
      end

      % fill in the rest of the descriptive fields
      if iscumsumcnt;
        freq.cumsumcnt = cat(catdim,cumsumcnt{:});
        clear cumsumcnt;
      end
      if iscumtapcnt;
        freq.cumtapcnt = cat(catdim,cumtapcnt{:});
        clear cumtapcnt;
      end
      if istrialinfo;
        freq.trialinfo = cat(catdim,trialinfo{:});
        clear trialinfo;
      end
    end

    freq.label = varargin{1}.label;
    freq.freq  = varargin{1}.freq;
    if isfield(varargin{1}, 'time'), freq.time = varargin{1}.time; end

  case 'chan'
    catdim = find(strcmp('chan', dimtok));
    if isempty(catdim)
      % try chancmb
      catdim = find(strcmp('chancmb', dimtok));
    elseif numel(catdim)>1
      error('ambiguous dimord for concatenation');
    end

    % check whether all channels are unique and throw an error if not
    [boolval, list] = checkchan(varargin{:}, 'unique');
    if ~boolval
      error('the input data structures have non-unique channels, concatenation across channel is not possible');
    end

    if isfield(varargin{1}, 'time')
      if ~checktime(varargin{:}, 'identical', tol)
        error('the input data structures have non-identical time bins, concatenation across channels not possible');
      end
    end

    if ~checkfreq(varargin{:}, 'identical', tol)
      error('the input data structures have non-identical frequency bins, concatenation across channels not possible');
    end

    % update the channel description
    freq.label = list;

    % fill in the rest of the descriptive fields
    freq.freq  = varargin{1}.freq;
    if isfield(varargin{1}, 'time'), freq.time = varargin{1}.time; end
    freq.dimord = varargin{1}.dimord;

  case 'freq'
    catdim = find(strcmp('freq', dimtok));

    % check whether all frequencies are unique and throw an error if not
    [boolval, list] = checkfreq(varargin{:}, 'unique', tol);
    if ~boolval
      error('the input data structures have non-unique frequency bins, concatenation across frequency is not possible');
    end

    if ~checkchan(varargin{:}, 'identical')
      error('the input data structures have non-identical channels, concatenation across frequency not possible');
    end
    if isfield(varargin{1}, 'time')
      if ~checktime(varargin{:}, 'identical', tol)
        error('the input data structures have non-identical time bins, concatenation across channels not possible');
      end
    end

    % update the frequency description
    freq.freq = list(:)';

    % fill in the rest of the descriptive fields
    freq.label  = varargin{1}.label;
    freq.dimord = varargin{1}.dimord;
    if isfield(varargin{1}, 'time'), freq.time = varargin{1}.time; end

  case 'time'
    catdim = find(strcmp('time', dimtok));

    % check whether all time points are unique and throw an error if not
    [boolval, list] = checktime(varargin{:}, 'unique', tol);
    if ~boolval
      error('the input data structures have non-unique time bins, concatenation across time is not possible');
    end

    if ~checkchan(varargin{:}, 'identical')
      error('the input data structures have non-identical channels, concatenation across time not possible');
    end
    if ~checkfreq(varargin{:}, 'identical', tol)
      error('the input data structures have non-identical frequency bins, concatenation across time not possible');
    end

    % update the time description
    freq.time = list(:)';

    % fill in the rest of the descriptive fields
    freq.label  = varargin{1}.label;
    freq.freq   = varargin{1}.freq;
    freq.dimord = varargin{1}.dimord;

  otherwise
    error('it is not allowed to concatenate across dimension %s',cfg.appenddim);
end

param = cfg.parameter;
if ~iscell(param), param = {param}; end

% are we appending along the channel dimension?
catchan = strcmp(cfg.appenddim, 'chan');
chandim = find(strcmp('chan', dimtok));

% concatenate the numeric data
for k = 1:numel(param)
  tmp = cell(1,Ndata);

  % get the numeric data 'param{k}' if present
  for m = 1:Ndata

    if ~isfield(varargin{m}, param{k})
      error('parameter %s is not present in all data sets', param{k});
    end
    tmp{m} = varargin{m}.(param{k});

    % if we are not appending along the channel dimension, make sure we
    % reorder the channel dimension across the different data sets. At this
    % point we can be sure that all data sets have identical channels.
    if ~catchan && m > 1
      [a,b] = match_str(varargin{1}.label, varargin{m}.label);
      if ~all(a==b)
        tmp{m} = reorderdim(tmp{m}, chandim, b);
      end
    end
  end

  if catdim==0,
    ndim    = length(size(tmp{1}));
    freq.(param{k}) = permute(cat(ndim+1,tmp{:}),[ndim+1 1:ndim]);
  else
    freq.(param{k}) = cat(catdim,tmp{:});
  end
end % for k = 1:numel(param)

% deal with the sensor information, if present
if isfield(varargin{1}, 'grad') || isfield(varargin{1}, 'elec')
  keepsensinfo = true;

  if isfield(varargin{1}, 'grad'), sensfield = 'grad'; end
  if isfield(varargin{1}, 'elec'), sensfield = 'elec'; end

  for k = 2:Ndata
    keepsensinfo = keepsensinfo && isequaln(varargin{1}.(sensfield), varargin{k}.(sensfield));
  end

  if keepsensinfo,
    freq.(sensfield) = varargin{1}.(sensfield);
  end
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous varargin
ft_postamble provenance freq
ft_postamble history freq
ft_postamble savevar freq


%---------------------------------------------
% subfunction to check uniqueness of freq bins
function [boolval, faxis] = checkfreq(varargin)

% last input is always the required string
tol      = varargin{end};
required = varargin{end-1};
varargin = varargin(1:end-2);

Ndata = numel(varargin);
Nfreq = zeros(1,Ndata);
faxis = zeros(1,0);
for i=1:Ndata
  Nfreq(i) = numel(varargin{i}.freq);
  faxis    = [faxis;varargin{i}.freq(:)];
end

if strcmp(required, 'unique')
  boolval = numel(unique(faxis))==numel(faxis) && ~all(isnan(faxis));
  % the second condition is included when the freq is set to dummy nan
elseif strcmp(required, 'identical')
  % the number of frequency bins needs at least to be the same across
  % inputs
  boolval = all(Nfreq==Nfreq(1));
  if boolval
    % then check whether the axes are equal
    faxis   = reshape(faxis, Nfreq(1), []);
    boolval = all(all(abs(faxis - repmat(faxis(:,1), 1, Ndata))<tol)==1);
    faxis   = faxis(:,1);
  end
end
