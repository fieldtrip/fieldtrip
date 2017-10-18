function [source] = ft_appendsource(cfg, varargin)

% FT_APPENDSOURCE concatenates multiple volumetric source reconstruction
% data structures that have been processed seperately.
%
% If the source reconstructions were computed for different ROIs or
% different slabs of a regular 3D grid (as indicated by the source
% positions), the data will be concatenated along the spatial dimension.
%
% If the source reconstructions were computed on the same source
% positions, but for different frequencies and/or latencies, e.g. for
% time-frequency spectrally decomposed data, the data will be concatenared
% along the frequency and/or time dimension.
%
% Use as
%   combined = ft_appendsource(cfg, source1, source2, ...)
%
% See also FT_SOURCEANALYSIS, FT_APPENDDATA, FT_APPENDFREQ, FT_APPENDSOURCE

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
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'source', 'feedback', 'yes', 'hassampleinfo', 'ifmakessense');
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
      boolval1 = checkpos(varargin{:}, 'identical');
      if isfield(varargin{1}, 'freq'),
        boolval2 = checkfreq(varargin{:}, 'identical', tol);
      else
        boolval2 = true;
      end

      if isfield(varargin{1}, 'time'),
        boolval3 = checktime(varargin{:}, 'identical', tol);
        if boolval1 && boolval2 && boolval3
          % each of the input datasets contains a single repetition (perhaps an average), these can be concatenated
          cfg.appenddim = 'rpt';
        elseif ~boolval1 && boolval2 && boolval3
          cfg.appenddim = 'pos';
        elseif boolval1 && ~boolval2 && boolval3
          cfg.appenddim = 'freq';
        elseif boolval1 && boolval2 && ~boolval3
          cfg.appenddim = 'time';
        else
          ft_error('the input datasets have multiple non-identical dimensions, this function only appends one dimension at a time');
        end
      else
        if boolval1 && boolval2
          % each of the input datasets contains a single repetition (perhaps an average), these can be concatenated
          cfg.appenddim = 'rpt';
        elseif ~boolval1 && boolval2
          cfg.appenddim = 'pos';
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
      ft_error('ambiguous dimord for concatenation');
    end

    % if any of these are present, concatenate
    % if not prepend the dimord with rpt (and thus shift the dimensions)

    % here we need to check whether the other dimensions are the same. if
    % not, consider some tolerance.
    boolval1 = checkpos(varargin{:}, 'identical', tol);
    if isfield(varargin{1}, 'freq'),
      boolval2 = checkfreq(varargin{:}, 'identical', tol);
    else
      boolval2 = true;
    end

    if isfield(varargin{1}, 'time'),
      boolval3 = checktime(varargin{:}, 'identical', tol);
    else
      boolval3 = true;
    end

    if any([boolval1 boolval2 boolval3]==false)
      ft_error('appending across observations is not possible, because the spatial, spectral and/or temporal dimensions are incompatible');
    end

    % select and reorder the channels that are in every dataset
    tmpcfg           = [];
    %tmpcfg.channel   = cfg.channel;
    tmpcfg.tolerance = cfg.tolerance;
    [varargin{:}]    = ft_selectdata(tmpcfg, varargin{:});
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
      hascumsumcnt = [];
      hascumtapcnt = [];
      hastrialinfo = [];
      for i=1:Ndata
        if isfield(varargin{i}, 'cumsumcnt');
          hascumsumcnt(end+1) = 1;
        else
          hascumsumcnt(end+1) = 0;
        end
        if isfield(varargin{i}, 'cumtapcnt');
          hascumtapcnt(end+1) = 1;
        else
          hascumtapcnt(end+1) = 0;
        end
        if isfield(varargin{i}, 'trialinfo');
          hastrialinfo(end+1) = 1;
        else
          hastrialinfo(end+1) = 0;
        end
      end

      % screen concatenable fields
      if ~checkfreq(varargin{:}, 'identical', tol)
        ft_error('the freq fields of the input data structures are not equal');
      else
        source.freq=varargin{1}.freq;
      end
      if ~sum(hascumsumcnt)==0 && ~(sum(hascumsumcnt)==Ndata);
        ft_error('the cumsumcnt fields of the input data structures are not equal');
      else
        iscumsumcnt=unique(hascumsumcnt);
      end
      if ~sum(hascumtapcnt)==0 && ~(sum(hascumtapcnt)==Ndata);
        ft_error('the cumtapcnt fields of the input data structures are not equal');
      else
        iscumtapcnt=unique(hascumtapcnt);
      end
      if ~sum(hastrialinfo)==0 && ~(sum(hastrialinfo)==Ndata);
        ft_error('the trialinfo fields of the input data structures are not equal');
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
        source.cumsumcnt = cat(catdim,cumsumcnt{:});
        clear cumsumcnt;
      end
      if iscumtapcnt;
        source.cumtapcnt = cat(catdim,cumtapcnt{:});
        clear cumtapcnt;
      end
      if istrialinfo;
        source.trialinfo = cat(catdim,trialinfo{:});
        clear trialinfo;
      end
    end

    source.pos   = varargin{1}.pos;
    if isfield(varargin{1}, 'freq'), source.freq = varargin{1}.freq; end
    if isfield(varargin{1}, 'time'), source.time = varargin{1}.time; end
    if isfield(varargin{1}, 'tri'),  source.tri  = varargin{1}.tri;  end

  case 'pos'
    % FIXME
    ft_error('this functionality does not work.....yet');


  case 'freq'
    % FIXME
    ft_error('this functionality does not work.....yet');

  case 'time'
    % FIXME
    ft_error('this functionality does not work.....yet');

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

if catdim==0,
  ndim    = length(size(tmp{1}));
  source.(param) = permute(cat(ndim+1,tmp{:}),[ndim+1 1:ndim]);
else
  source.(param) = cat(catdim,tmp{:});
end


% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous varargin
ft_postamble provenance source
ft_postamble history source
ft_postamble savevar source
