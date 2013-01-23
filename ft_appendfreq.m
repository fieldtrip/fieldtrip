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
%                   freq2.powspecrum etc, use cft.parameter = 'powspctrm'.
%
% The configuration can optionally contain
%  cfg.appenddim  = String. The dimension to concatenate over (default:
%                   'auto').
%  cfg.tolerance  = Double. Tolerance determines how different the units of
%                   frequency structures are allowed to be to be considered
%                   compatible (default: 1e-5).
%
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following options:
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
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble provenance
ft_preamble trackconfig
ft_preamble debug
ft_preamble loadvar varargin

% check if the input data is valid for this function
for i=1:length(varargin)
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'freq', 'feedback', 'yes');
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'required', 'parameter');

% set the defaults
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
    % determine the appenddim and recursively call ft_appendfreq
    tmpcfg = cfg;
    
    % only allow to append across observations if these are present in the data
    if any(strcmp(dimtok, 'rpt'))
      tmpcfg.appenddim = 'rpt';
      
    elseif any(strcmp(dimtok, 'rpttap'))
      tmpcfg.appenddim = 'rpttap';
      
    elseif any(strcmp(dimtok, 'subj'))
      tmpcfg.appenddim = 'subj';
      
    else
      % we need to check whether the other dimensions are the same.
      % if not, consider some tolerance.
      boolval1 = checkchan(varargin{:}, 'identical');
      boolval2 = checkfreq(varargin{:}, 'identical', tol);
      
      if isfield(varargin{1}, 'time'),
        boolval3 = checktime(varargin{:}, 'identical', tol);
        if boolval1 && boolval2 && boolval3
          % each of the input datasets contains a single repetition (perhaps an average), these can be concatenated
          tmpcfg.appenddim = 'rpt';
        elseif ~boolval1 && boolval2 && boolval3
          tmpcfg.appenddim = 'chan';
        elseif boolval1 && ~boolval2 && boolval3
          tmpcfg.appenddim = 'freq';
        elseif boolval1 && boolval2 && ~boolval3
          tmpcfg.appenddim = 'time';
        end
      else
        if boolval1 && boolval2
          % each of the input datasets contains a single repetition (perhaps an average), these can be concatenated
          tmpcfg.appenddim = 'rpt';
        elseif ~boolval1 && boolval2
          tmpcfg.appenddim = 'chan';
        elseif boolval1 && ~boolval2
          tmpcfg.appenddim = 'freq';
        end
      end
      
      freq = ft_appendfreq(tmpcfg, varargin{:});
      return;
    end % determining the dimension for appending
    
    % otherwise we need to check whether the other dimensions are the same. if
    % not, consider some tolerance.
    boolval1 = checkchan(varargin{:}, 'identical');
    boolval2 = checkfreq(varargin{:}, 'identical', tol);
    if isfield(varargin{1}, 'time'),
      boolval3 = checktime(varargin{:}, 'identical', tol);
      if boolval1 && boolval2 && boolval3
        tmpcfg.appenddim = 'rpt';
      elseif ~boolval1 && boolval2 && boolval3
        tmpcfg.appenddim = 'chan';
      elseif boolval1 && ~boolval2 && boolval3
        tmpcfg.appenddim = 'freq';
      elseif boolval1 && boolval2 && ~boolval3
        tmpcfg.appenddim = 'time';
      end
    else
      if boolval1 && boolval2
        tmpcfg.appenddim = 'rpt';
      elseif ~boolval1 && boolval2
        tmpcfg.appenddim = 'chan';
      elseif boolval1 && ~boolval2
        tmpcfg.appenddim = 'freq';
      end
    end
    freq = ft_appendfreq(tmpcfg, varargin{:});
    return;
    
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
    
    if any([boolval1 boolval2 boolval3]==false)
      error('appending across observations is not possible, because the dimensions are incompatible');
    end
    
    % update the dimord
    if catdim==0
      freq.dimord = ['rpt_',varargin{1}.dimord];
      % FIXME append dof
    else
      freq.dimord = varargin{1}.dimord;
      % FIXME append cumtapcnt cumsumcnt trialinfo dof
    end
    
    % fill in the rest of the descriptive fields
    freq.label = varargin{1}.label;
    freq.freq  = varargin{1}.freq;
    if isfield(varargin{1}, 'time'), freq.time = varargin{1}.time; end
    
  case 'chan'
    catdim = strmatch('chan', dimtok);
    if isempty(catdim)
      % try chancmb
      catdim = strmatch('chancmb', dimtok);
    elseif numel(catdim)>1
      error('ambiguous dimord for concatenation');
    end
    
    % check whether all channels are unique and throw an error if not
    [boolval, list] = checkchan(varargin{:}, 'unique');
    if ~boolval
      error('the input data structures have non-unique channels, concatenation across channel is not possible');
    end
    
    % update the channel description
    freq.label = list;
    
    % fill in the rest of the descriptive fields
    freq.freq  = varargin{1}.freq;
    if isfield(varargin{1}, 'time'), freq.time = varargin{1}.time; end
    freq.dimord = varargin{1}.dimord;
    
  case 'freq'
    catdim = strmatch('freq', dimtok);
    
    % check whether all frequencies are unique and throw an error if not
    [boolval, list] = checkfreq(varargin{:}, 'unique', tol);
    if ~boolval
      error('the input data structures have non-unique frequency bins, concatenation across frequency is not possible');
    end
    
    % update the frequency description
    freq.freq = list(:)';
    
    % fill in the rest of the descriptive fields
    freq.label  = varargin{1}.label;
    freq.dimord = varargin{1}.dimord;
    if isfield(varargin{1}, 'time'), freq.time = varargin{1}.time; end
    
  case 'time'
    catdim = strmatch('time', dimtok);
    
    % check whether all time points are unique and throw an error if not
    [boolval, list] = checktime(varargin{:}, 'unique', tol);
    if ~boolval
      error('the input data structures have non-unique time bins, concatenation across time is not possible');
    end
    
    % update the time description
    freq.time = list(:)';
    
    % fill in the rest of the descriptive fields
    freq.label  = varargin{1}.label;
    freq.freq   = varargin{1}.freq;
    freq.dimord = varargin{1}.dimord;
    
  otherwise
end

% FIXME do a check on whether the parameters are present in all datasets
param = cfg.parameter;
if ~iscell(param), param = {param}; end

% concatenate the numeric data
for k = 1:numel(param)
  tmp = cell(1,Ndata);
  % get the numeric data 'param{k}' if present
  for m = 1:Ndata
    tmp{m} = varargin{m}.(param{k});
  end
  
  if catdim==0,
    ndim    = length(size(tmp{1}));
    freq.(param{k}) = permute(cat(ndim+1,tmp{:}),[ndim+1 1:ndim]);
  else
    freq.(param{k}) = cat(catdim,tmp{:});
  end
end % for k = 1:numel(param)

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
ft_postamble previous varargin
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

%---------------------------------------------
% subfunction to check uniqueness of time bins
function [boolval, taxis] = checktime(varargin)

% last input is always the required string
tol      = varargin{end};
required = varargin{end-1};
varargin = varargin(1:end-2);

Ndata = numel(varargin);
Ntime = zeros(1,Ndata);
taxis = zeros(1,0);
for i=1:Ndata
  Ntime(i) = numel(varargin{i}.time);
  taxis    = [taxis;varargin{i}.time(:)];
end

if strcmp(required, 'unique')
  boolval = numel(unique(taxis))==numel(taxis) && ~all(isnan(taxis));
  % the second condition is included when the time is set to dummy nan
elseif strcmp(required, 'identical')
  % the number of time bins needs at least to be the same across
  % inputs
  boolval = all(Ntime==Ntime(1));
  if boolval
    % then check whether the axes are equal
    taxis   = reshape(taxis, Ntime(1), []);
    boolval = all(all(abs(taxis - repmat(taxis(:,1), 1, Ndata))<tol)==1);
    taxis   = taxis(:,1);
  end
end

%--------------------------------------------------
% subfunction to check uniqueness of channel labels
function [boolval, list] = checkchan(varargin)

% last input is always the required string
required = varargin{end};
varargin = varargin(1:end-1);

Ndata = numel(varargin);
Nchan = zeros(1,Ndata);
list  = cell(0,1);
for i=1:Ndata
  Nchan(i) = numel(varargin{i}.label);
  list     = [list;varargin{i}.label(:)];
end

if strcmp(required, 'unique')
  boolval = numel(unique(list))==numel(list);
elseif strcmp(required, 'identical')
  boolval = 1;
end

