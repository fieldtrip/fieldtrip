function [source] = ft_appendsource(cfg, varargin)

% FT_APPENDSOURCE concatenates multiple volumetric source reconstruction data
% structures that have been processed separately.
%
% Use as
%   combined = ft_appendsource(cfg, source1, source2, ...)
%
% If the source reconstructions were computed for different ROIs or different slabs
% of a regular 3D grid (as indicated by the source positions), the data will be
% concatenated along the spatial dimension.
%
% If the source reconstructions were computed on the same source positions, but for
% different frequencies and/or latencies, e.g. for time-frequency spectrally
% decomposed data, the data will be concatenated along the frequency and/or time
% dimension, but only of the frequency or time axes are well-behaved, i.e. all data
% points along the dimension of interest should be sortable across data objects;
% interleaving across data objects is not possible.
%
% See also FT_SOURCEANALYSIS, FT_DATATYPE_SOURCE, FT_APPENDDATA, FT_APPENDTIMELOCK,
% FT_APPENDFREQ

% Copyright (C) 2011-2021, Robert Oostenveld
% Copyright (C) 2022-, Jan-Mathijs Schoffelen and Robert Oostenveld
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

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
for i=1:length(varargin)
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'source', 'feedback', 'yes', 'hassampleinfo', 'ifmakessense');
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'required', 'parameter');

% set the defaults
cfg.appenddim  = ft_getopt(cfg, 'appenddim',  'auto');
cfg.tolerance  = ft_getopt(cfg, 'tolerance',  1e-5); % this is passed to ft_selectdata

% do a basic check to see whether the dimords match
Ndata = length(varargin);
dimord = cell(1,Ndata);
for i=1:Ndata
  dimord{i} = getdimord(varargin{i}, cfg.parameter);
end
dimordmatch = all(strcmp(dimord{1}, dimord));

if ~dimordmatch
  ft_error('the dimords of the input data structures are not equal');
end

% create the output structure from scratch
source = [];

tol    = cfg.tolerance;
dimtok = tokenize(dimord{1}, '_');

matchpos  = checkpos(varargin{:}, 'identical', tol);
matchfreq = true;
matchtime = true;
if isfield(varargin{1}, 'freq'), matchfreq = checkfreq(varargin{:}, 'identical', tol); end
if isfield(varargin{1}, 'time'), matchtime = checktime(varargin{:}, 'identical', tol); end

if strcmp(cfg.appenddim, 'auto')
  % if there are observations in the data, only allow to append across observations
  if any(strcmp(dimtok, 'rpt'))
    cfg.appenddim = 'rpt';
  elseif any(strcmp(dimtok, 'rpttap'))
    cfg.appenddim = 'rpttap';
  elseif any(strcmp(dimtok, 'subj'))
    cfg.appenddim = 'subj';
  else
    % try to guess the appenddim
    if matchpos
      % append across either time or freq, provided one of the other match
      if matchfreq && matchtime
        cfg.appenddim = 'rpt';
      elseif ~matchfreq && matchtime
        cfg.appenddim = 'freq';
      elseif matchfreq && ~matchtime
        cfg.appenddim = 'time';
      elseif ~matchfreq && ~matchtimee
        ft_error('the input datasets have multiple non-identical dimensions, this function only appends one dimension at a time');
      end
    elseif ~matchpos
      % append across pos, provided time/freq (if present) are matched
      if matchfreq && matchtime
        cfg.appenddim = 'pos';
      else
        ft_error('the input datasets have multiple non-identical dimensions, this function only appends one dimension at a time');
      end
    end
  end % determine the dimension for appending
end % if auto

catdim = find(ismember(dimtok, cfg.appenddim));
switch cfg.appenddim
  case {'rpt' 'rpttap' 'subj'}
    if numel(catdim)==0
      catdim = 0;
    elseif numel(catdim)==1
      % this is OK
    elseif numel(catdim)>1
      ft_error('ambiguous dimord for concatenation');
    end
    
    % if any of these are present, concatenate
    % if not prepend the dimord with rpt (and thus shift the dimensions)
    
    % here we need to check whether the other dimensions are the same. if
    % not, consider some tolerance.
    if any([matchpos matchfreq matchtime]==false)
      ft_error('appending across observations is not possible, because the spatial, spectral and/or temporal dimensions are incompatible');
    end
    
    % determine the union of all input data
    tmpcfg = keepfields(cfg, {'tolerance', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
    [varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
    for i=1:Ndata
      [cfg_rolledback, varargin{i}] = rollback_provenance(cfg, varargin{i});
    end
    cfg = cfg_rolledback;
    
    % update the dimord
    if catdim==0
      source.dimord = ['rpt_',dimord{1}];
      % FIXME append dof
    else
      source.dimord = dimord{1};
      % before append cumtapcnt cumsumcnt trialinfo check if there's a
      % subfield in each dataset. Append fields of different dataset might
      % lead in empty and/or non-existing fields in a particular dataset
      hascumsumcnt = false(1, Ndata);
      hascumtapcnt = false(1, Ndata);
      hastrialinfo = false(1, Ndata);
      for i=1:Ndata
        if isfield(varargin{i}, 'cumsumcnt'), hascumsumcnt(i) = true; end
        if isfield(varargin{i}, 'cumtapcnt'), hascumtapcnt(i) = true; end
        if isfield(varargin{i}, 'trialinfo'), hastrialinfo(i) = true; end
      end
      
      % screen concatenable fields
      if ~matchfreq
        ft_error('the freq fields of the input data structures are not equal');
      else
        source.freq=varargin{1}.freq;
      end
      if ~sum(hascumsumcnt)==0 && ~(sum(hascumsumcnt)==Ndata)
        ft_error('the cumsumcnt fields of the input data structures are not equal');
      else
        iscumsumcnt=unique(hascumsumcnt);
      end
      if ~sum(hascumtapcnt)==0 && ~(sum(hascumtapcnt)==Ndata)
        ft_error('the cumtapcnt fields of the input data structures are not equal');
      else
        iscumtapcnt=unique(hascumtapcnt);
      end
      if ~sum(hastrialinfo)==0 && ~(sum(hastrialinfo)==Ndata)
        ft_error('the trialinfo fields of the input data structures are not equal');
      else
        istrialinfo=unique(hastrialinfo);
      end
      
      % concatenating fields
      for i=1:Ndata
        if iscumsumcnt
          cumsumcnt{i}=varargin{i}.cumsumcnt;
        end
        if iscumtapcnt
          cumtapcnt{i}=varargin{i}.cumtapcnt;
        end
        if istrialinfo
          trialinfo{i}=varargin{i}.trialinfo;
        end
      end
      
      % fill in the rest of the descriptive fields
      if iscumsumcnt
        source.cumsumcnt = cat(catdim,cumsumcnt{:});
        clear cumsumcnt;
      end
      if iscumtapcnt
        source.cumtapcnt = cat(catdim,cumtapcnt{:});
        clear cumtapcnt;
      end
      if istrialinfo
        source.trialinfo = cat(catdim,trialinfo{:});
        clear trialinfo;
      end
    end
    
    source.pos = varargin{1}.pos;
    if isfield(varargin{1}, 'freq'), source.freq = varargin{1}.freq; end
    if isfield(varargin{1}, 'time'), source.time = varargin{1}.time; end
    if isfield(varargin{1}, 'tri'),  source.tri  = varargin{1}.tri;  end % FIXME assumed equal
    if isfield(varargin{1}, 'dim'),  source.dim  = varargin{1}.dim;  end % FIXME assumed equal
    
  case 'pos'
    % ensure that the positions are in the same units
    if ~isfield(varargin{1}, 'unit')
      varargin{1} = ft_determine_units(varargin{1});
    end

    if ~matchfreq || ~matchtime
      ft_error('when appending across positions, all other dimensions should match');
    end
    
    pos = cell(Ndata,1);
    tri = cell(0,1);
    Ntotal = 0;
    for i=1:Ndata
      varargin{i} = ft_convert_units(varargin{i}, varargin{1}.unit);
      pos{i} = varargin{i}.pos;
      if isfield(varargin{i}, 'tri')
        tri{end+1} = varargin{i}.tri + Ntotal;
      end
      Ntotal = Ntotal + size(varargin{i}.pos,1);
    end
    source = keepfields(varargin{1}, {'time' 'freq'});
    source.pos = cat(1, pos{:});
    if ~isempty(tri)
      source.tri = cat(1, tri{:});
    end

  case 'freq'
    
    if ~matchtime || ~matchpos
      ft_error('when appending across time, all other dimensions should match');
    end

    % ensure that contatenation along the freq dimension can be done by a
    % simple cat operation along the specified dimension, i.e. the
    % individual frequency axes should not be interleaved, and unique
    freq = cell(Ndata, 1);
    freqbounds = zeros(Ndata, 2);
    for i=1:Ndata
      freq{i} = varargin{i}.freq;
      freqbounds(i,:) = [freq{i}(1) freq{i}(end)];
    end
    ok = double(freqbounds(:,1)>freqbounds(:,2)') + double(freqbounds(:,2)<freqbounds(:,1)');

    if any(ok(:)==1)
      ft_error('when appending across freq, there should be no overlap in the individual frequency axes');
    end

    [srt, ix] = sort(freqbounds(:,1));
    if ~isequal(ix(:), (1:Ndata)')
      ft_warning('reordering the input data arguments to get a well-behaved frequency axis');
      varargin = varargin(ix);
      freq     = freq(ix);
    end

    source = keepfields(varargin{1}, {'pos' 'dim' 'time' 'tri'});
    source.freq = cat(2, freq{:});

  case 'time'
   
    if ~matchfreq || ~matchpos
      ft_error('when appending across time, all other dimensions should match');
    end
    
    % ensure that contatenation along the time dimension can be done by a
    % simple cat operation along the specified dimension, i.e. the
    % individual time axes should not be interleaved, and unique
    time = cell(1,Ndata);
    timebounds = zeros(Ndata, 2);
    for i=1:Ndata
      time{i} = varargin{i}.time;
      timebounds(i,:) = [time{i}(1) time{i}(end)];
    end
    ok = double(timebounds(:,1)>timebounds(:,2)') + double(timebounds(:,2)<timebounds(:,1)');
    
    if any(ok(:)==2) || all(ok(:)==0)
      ft_error('when appending across time, there should be no overlap in the individual time axes');
    end

    [srt, ix] = sort(timebounds(:,1));
    if ~isequal(ix(:), (1:Ndata)')
      ft_warning('reordering the input data arguments to get a well-behaved time axis');
      varargin = varargin(ix);
      time     = time(ix);
    end

    source = keepfields(varargin{1}, {'pos' 'dim' 'freq' 'tri'});
    source.time = cat(2, time{:});

  otherwise
    ft_error('it is not allowed to concatenate across dimension %s',cfg.appenddim);
end

param = cfg.parameter;
if iscell(param) && numel(param)>1
  ft_error('it is not possible yet to append multiple parameters in a single call');
elseif iscell(param)
  param = param{1};
end

% concatenate the numeric data
tmp = cell(1,Ndata);

% get the numeric data 'param{k}' if present
for m = 1:Ndata
  
  if ~isfield(varargin{m}, param)
    ft_error('parameter %s is not present in all data sets', param);
  end
  tmp{m} = varargin{m}.(param);
  
end

if catdim==0
  ndim    = length(size(tmp{1}));
  source.(param) = permute(cat(ndim+1,tmp{:}),[ndim+1 1:ndim]);
else
  source.(param) = cat(catdim,tmp{:});
end


% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous varargin
ft_postamble provenance source
ft_postamble history source
ft_postamble savevar source
