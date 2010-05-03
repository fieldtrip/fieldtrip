function [data] = checkdata(data, varargin)

% CHECKDATA checks the input data of the main FieldTrip functions, e.g. whether
% the type of data strucure corresponds with the required data. If neccessary
% and possible, this function will adjust the data structure to the input
% requirements (e.g. change dimord, average over trials, convert inside from
% index into logical).
%
% If the input data does NOT correspond to the requirements, this function
% is supposed to give a elaborate warning message and if applicable point
% the user to external documentation (link to website).
%
% Use as
%   [data] = checkdata(data, ...)
%
% Optional input arguments should be specified as key-value pairs and can include
%   feedback           = yes, no
%   datatype           = raw, freq, timelock, comp, spike, source, volume, dip
%   dimord             = any combination of time, freq, chan, refchan, rpt, subj, chancmb, rpttap, pos
%   senstype           = ctf151, ctf275, ctf151_planar, ctf275_planar, neuromag122, neuromag306, bti148, bti248, bti248_planar, magnetometer, electrode
%   inside             = logical, index
%   ismeg              = yes, no
%   hastrials          = yes, no
%   hasunits           = yes, no
%   hastrialdef        = yes, no
%   hasoffset          = yes, no (only applies to raw data)
%   hascumtapcnt       = yes, no (only applies to freq data)
%   hasdof             = yes, no
%   cmbrepresentation  = sparse, full (applies to covariance and cross-spectral density)
%
% For some options you can specify multiple values, e.g.
%   [data] = checkdata(data, 'senstype', {'ctf151', 'ctf275'}), e.g. in megrealign
%   [data] = checkdata(data, 'datatype', {'timelock', 'freq'}), e.g. in sourceanalysis

% Copyright (C) 2007-2009, Robert Oostenveld
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

% in case of an error this function could use dbstack for more detailled
% user feedback
%
% this function should replace/encapsulate
%   fixdimord
%   fixinside
%   fixprecision
%   fixvolume
%   data2raw
%   raw2data
%   grid2transform
%   transform2grid
%   fourier2crsspctrm
%   freq2cumtapcnt
%   sensortype
%   time2offset
%   offset2time
%
% other potential uses for this function:
%   time -> offset in freqanalysis
%   average over trials
%   csd as matrix

% get the optional input arguments
feedback      = keyval('feedback',      varargin); if isempty(feedback), feedback = 'no'; end
dtype         = keyval('datatype',      varargin); % should not conflict with the datatype function
dimord        = keyval('dimord',        varargin);
stype         = keyval('senstype',      varargin); % senstype is a function name which should not be masked
ismeg         = keyval('ismeg',         varargin);
inside        = keyval('inside',        varargin); % can be logical or index
hastrials     = keyval('hastrials',     varargin);
hasunits      = keyval('hasunits',      varargin);
hastrialdef   = keyval('hastrialdef',   varargin); if isempty(hastrialdef), hastrialdef = 'no'; end
hasoffset     = keyval('hasoffset',     varargin); if isempty(hasoffset), hasoffset = 'no'; end
hasdimord     = keyval('hasdimord',     varargin); if isempty(hasdimord), hasdimord = 'no'; end
hascumtapcnt  = keyval('hascumtapcnt',  varargin);
hasdof        = keyval('hasdof',        varargin); if isempty(hasdof), hasdof = 'no'; end
haspow        = keyval('haspow',        varargin); if isempty(haspow), haspow = 'no'; end
cmbrepresentation = keyval('cmbrepresentation',  varargin);
channelcmb    = keyval('channelcmb',   varargin);
sourcedimord  = keyval('sourcedimord', varargin);
sourcerepresentation = keyval('sourcerepresentation', varargin);
keepoutside   = keyval('keepoutside',  varargin);

% determine the type of input data
% this can be raw, freq, timelock, comp, spike, source, volume, dip
israw      = datatype(data, 'raw');
isfreq     = datatype(data, 'freq');
istimelock = datatype(data, 'timelock');
iscomp     = datatype(data, 'comp');
isspike    = datatype(data, 'spike');
isvolume   = datatype(data, 'volume');
issource   = datatype(data, 'source');
isdip      = datatype(data, 'dip');
ismvar     = datatype(data, 'mvar');
isfreqmvar = datatype(data, 'freqmvar');

% FIXME use the istrue function on ismeg and hasxxx options

if ~isequal(feedback, 'no')
  if israw
    nchan = length(data.label);
    ntrial = length(data.trial);
    fprintf('the input is raw data with %d channels and %d trials\n', nchan, ntrial);
  elseif isfreq
    nchan = length(data.label);
    nfreq = length(data.freq);
    if isfield(data, 'time'), ntime = num2str(length(data.time)); else ntime = 'no'; end
    fprintf('the input is freq data with %d channels, %d frequencybins and %s timebins\n', nchan, nfreq, ntime);
  elseif istimelock
    nchan = length(data.label);
    ntime = length(data.time);
    fprintf('the input is timelock data with %d channels and %d timebins\n', nchan, ntime);
  elseif iscomp
    ncomp = length(data.label);
    nchan = length(data.topolabel);
    fprintf('the input is component data with %d components and %d original channels\n', ncomp, nchan);
  elseif isspike
    nchan = length(data.label);
    fprintf('the input is spike data\n');
  elseif isvolume
    fprintf('the input is volume data with dimensions [%d %d %d]\n', data.dim(1), data.dim(2), data.dim(3));
  elseif issource
    nsource = size(data.pos, 1);
    fprintf('the input is source data with %d positions\n', nsource);
  elseif isdip
    fprintf('the input is dipole data\n');
  elseif ismvar
    fprintf('the input is mvar data\n');
  elseif isfreqmvar
    fprintf('the input is freqmvar data\n');
  end
end % give feedback

if isfreq || istimelock || iscomp || issource || isvolume
  % ensure consistency between the dimord string and the axes that describe the data dimensions
  data = fixdimord(data);
end

if istimelock
  % remove the unwanted fields
  if isfield(data, 'numsamples'),       data = rmfield(data, 'numsamples');       end
  if isfield(data, 'numcovsamples'),    data = rmfield(data, 'numcovsamples');    end
  if isfield(data, 'numblcovsamples'),  data = rmfield(data, 'numblcovsamples');  end
end

if issource && isvolume
  % it should be either one or the other
  % the choice here is to represent it as volume description since that is simpler to handle
  % remove the unwanted fields
  if isfield(data, 'pos'),    data = rmfield(data, 'pos');    end
  if isfield(data, 'xgrid'),  data = rmfield(data, 'xgrid');  end
  if isfield(data, 'ygrid'),  data = rmfield(data, 'ygrid');  end
  if isfield(data, 'zgrid'),  data = rmfield(data, 'zgrid');  end
  issource = false;
end

if issource || isvolume
  % these don't contain a dimord but in the frequency domain case they could 
  % contain a .frequency field rather than a .freq field
  if isfield(data, 'frequency'), 
    data.freq = data.frequency;
    data      = rmfield(data, 'frequency');
  end
end

if ~isempty(dtype)
  if ~isa(dtype, 'cell')
    dtype = {dtype};
  end

  okflag = 0;
  for i=1:length(dtype)
    % check that the data matches with one or more of the required datatypes
    switch dtype{i}
      case 'raw'
        okflag = okflag + israw;
      case 'freq'
        okflag = okflag + isfreq;
      case 'timelock'
        okflag = okflag + istimelock;
      case 'comp'
        okflag = okflag + iscomp;
      case 'spike'
        okflag = okflag + isspike;
      case 'volume'
        okflag = okflag + isvolume;
      case 'source'
        okflag = okflag + issource;
      case 'dip'
        okflag = okflag + isdip;
      case 'mvar'
        okflag = okflag + ismvar;
      case 'freqmvar'
        okflag = okflag + isfreqmvar;
    end % switch dtype
  end % for dtype
  
  if ~okflag
    % try to convert the data
    for iCell = 1:length(dtype)
      if isequal(dtype(iCell), {'source'}) && isvolume
        data = volume2source(data);
        isvolume = 0;
        issource = 1;
        okflag = 1;
      elseif isequal(dtype(iCell), {'volume'}) && issource
        data = source2volume(data);
        isvolume = 1;
        issource = 0;
        okflag = 1;
      elseif isequal(dtype(iCell), {'raw'}) && issource
        data = data2raw(data);
        issource = 0;
        israw = 1;
        okflag = 1;
      elseif isequal(dtype(iCell), {'raw'}) && istimelock
        data = timelock2raw(data);
        istimelock = 0;
        israw = 1;
        okflag = 1;
      elseif isequal(dtype(iCell), {'timelock'}) && israw
        data = raw2timelock(data);
        israw = 0;
        istimelock = 1;
        okflag = 1;
      elseif isequal(dtype(iCell), {'raw'}) && isfreq
        data = freq2raw(data);
        isfreq = 0;
        israw = 1;
        okflag = 1;
      elseif isequal(dtype(iCell), {'raw'}) && iscomp
        data = comp2raw(data);
        iscomp = 0;
        israw = 1;
        okflag = 1;
      end
    end % for iCell
  end % if okflag

  if ~okflag
    % construct an error message
    if length(dtype)>1
      str = sprintf('%s, ', dtype{1:(end-2)});
      str = sprintf('%s%s or %s', str, dtype{end-1}, dtype{end});
    else
      str = dtype{1};
    end
    str = sprintf('This function requires %s data as input.', str);
    error(str);
  end % if okflag
end

if ~isempty(dimord)
  if ~isa(dimord, 'cell')
    dimord = {dimord};
  end

  if isfield(data, 'dimord')
    okflag = any(strcmp(data.dimord, dimord));
  else
    okflag = 0;
  end

  if ~okflag
    % construct an error message
    if length(dimord)>1
      str = sprintf('%s, ', dimord{1:(end-2)});
      str = sprintf('%s%s or %s', str, dimord{end-1}, dimord{end});
    else
      str = dimord{1};
    end
    str = sprintf('This function requires data with a dimord of %s.', str);
    error(str);
  end % if okflag
end

if ~isempty(stype)
  if ~isa(stype, 'cell')
    stype = {stype};
  end

  if isfield(data, 'grad') || isfield(data, 'elec')
    if any(strcmp(ft_senstype(data), stype));
      okflag = 1;
    else
      okflag = 0;
    end
  else
    okflag = 0;
  end

  if ~okflag
    % construct an error message
    if length(stype)>1
      str = sprintf('%s, ', stype{1:(end-2)});
      str = sprintf('%s%s or %s', str, stype{end-1}, stype{end});
    else
      str = stype{1};
    end
    str = sprintf('This function requires %s data as input, but you are giving %s data.', str, ft_senstype(data));
    error(str);
  end % if okflag
end

if ~isempty(ismeg)
  if isequal(ismeg, 'yes')
    okflag = isfield(data, 'grad');
  elseif isequal(ismeg, 'no')
    okflag = ~isfield(data, 'grad');
  end

  if ~okflag && isequal(ismeg, 'yes')
    error('This function requires MEG data with a ''grad'' field');
  elseif ~okflag && isequal(ismeg, 'no')
    error('This function should not be given MEG data with a ''grad'' field');
  end % if okflag
end

if ~isempty(inside)
  % TODO absorb the fixinside function into this code
  data   = fixinside(data, inside);
  okflag = isfield(data, 'inside');

  if ~okflag
    % construct an error message
    error('This function requires data with an ''inside'' field.');
  end % if okflag
end

%if isvolume
%  % ensure consistent dimensions of the volumetric data
%  % reshape each of the volumes that is found into a 3D array
%  param = parameterselection('all', data);
%  dim   = data.dim;
%  for i=1:length(param)
%    tmp  = getsubfield(data, param{i});
%    tmp  = reshape(tmp, dim);
%    data = setsubfield(data, param{i}, tmp);
%  end
%end

if issource || isvolume,
  % these are not used any more
  if isfield(data, 'xgrid'),  data = rmfield(data, 'xgrid');  end
  if isfield(data, 'ygrid'),  data = rmfield(data, 'ygrid');  end
  if isfield(data, 'zgrid'),  data = rmfield(data, 'zgrid');  end

  if isequal(hasunits, 'yes') && ~isfield(data, 'units')
    % calling convert_units with only the input data adds the units without converting
    data = ft_convert_units(data);
  end
  
  % the following section is to make a dimord-consistent representation of
  % volume and source data, taking trials, time and frequency into account
  if isequal(hasdimord, 'yes') && (~isfield(data, 'dimord') || ~strcmp(data.dimord,sourcedimord))

    % determine the size of the data
    if isfield(data, 'dimord'),
      dimtok = tokenize(data.dimord, '_');
      if ~isempty(strmatch('time', dimtok)), Ntime = length(data.time); else Ntime = 1; end
      if ~isempty(strmatch('freq', dimtok)), Nfreq = length(data.freq); else Nfreq = 1; end
    else
      Nfreq = 1;
      Ntime = 1;
    end

    %convert old style source representation into new style
    if isfield(data, 'avg') && isfield(data.avg, 'mom') && (isfield(data, 'freq') || isfield(data, 'frequency')) && strcmp(sourcedimord, 'rpt_pos'),
      %frequency domain source representation convert to single trial power
      Npos   = size(data.pos,1);
      Nrpt   = length(data.cumtapcnt);
      tmpmom = zeros(Npos, size(data.avg.mom{data.inside(1)},2));
      tmpmom(data.inside,:) = cat(1,data.avg.mom{data.inside});
      tmppow = zeros(Npos, Nrpt);
      tapcnt = [0;cumsum(data.cumtapcnt)];
      for k = 1:Nrpt
        Ntap = tapcnt(k+1)-tapcnt(k);
        tmppow(data.inside,k) = sum(abs(tmpmom(data.inside,(tapcnt(k)+1):tapcnt(k+1))).^2,2)./Ntap;
      end
      data.pow = tmppow';
      data     = rmfield(data, 'avg');
      if strcmp(inside, 'logical'),
        data     = fixinside(data, 'logical');
        data.inside = repmat(data.inside(:)',[Nrpt 1]);
      end
    elseif isfield(data, 'avg') && isfield(data.avg, 'mom') && (isfield(data, 'freq') || isfield(data, 'frequency')) && strcmp(sourcedimord, 'rpttap_pos'),
      %frequency domain source representation convert to single taper fourier coefficients
      Npos   = size(data.pos,1);
      Nrpt   = sum(data.cumtapcnt);
      data.fourierspctrm = complex(zeros(Nrpt, Npos), zeros(Nrpt, Npos));
      data.fourierspctrm(:, data.inside) = transpose(cat(1, data.avg.mom{data.inside}));
      data   = rmfield(data, 'avg');
    elseif isfield(data, 'avg') && isfield(data.avg, 'mom') && isfield(data, 'time') && strcmp(sourcedimord, 'pos_time'),
      Npos   = size(data.pos,1);
      Nrpt   = 1;
      tmpmom = zeros(Npos, size(data.avg.mom{data.inside(1)},2));
      tmpmom(data.inside,:) = cat(1,data.avg.mom{data.inside});
      data.mom = tmpmom;
      if isfield(data.avg, 'noise'),
        tmpnoise = data.avg.noise(:);
        data.noise = tmpnoise(:,ones(1,size(tmpmom,2)));
      end
      data = rmfield(data, 'avg');
      Ntime = length(data.time);
    elseif isfield(data, 'trial') && isfield(data.trial(1), 'mom') && isfield(data, 'time') && strcmp(sourcedimord, 'rpt_pos_time'),
      Npos   = size(data.pos,1);
      Nrpt   = length(data.trial);
      Ntime  = length(data.time);
      tmpmom = zeros(Nrpt, Npos, Ntime);
      for k = 1:Nrpt
        tmpmom(k,data.inside,:) = cat(1,data.trial(k).mom{data.inside});
      end
      data     = rmfield(data, 'trial');
      data.mom = tmpmom;
    elseif isfield(data, 'trial') && isstruct(data.trial)
      Nrpt = length(data.trial);
    else
      Nrpt = 1;
    end

    % start with an initial specification of the dimord and dim
    if (~isfield(data, 'dim') || ~isfield(data, 'dimord'))
      if issource
        % at least it should have a Nx3 pos
        data.dim    = size(data.pos, 1);
        data.dimord = 'pos';
      elseif isvolume
        % at least it should have a 1x3 dim
        data.dim    = data.dim;
        data.dimord = 'dim1_dim2_dim3';
      end
    end

    % add the additional dimensions
    if Nfreq>1
      data.dimord = [data.dimord '_freq'];
      data.dim    = [data.dim     Nfreq];
    end
    if Ntime>1
      data.dimord = [data.dimord '_time'];
      data.dim    = [data.dim     Ntime];
    end
    if Nrpt>1 && strcmp(sourcedimord, 'rpt_pos'),
      data.dimord = ['rpt_' data.dimord];
      data.dim    = [Nrpt   data.dim ];
    elseif Nrpt>1 && strcmp(sourcedimord, 'rpttap_pos'),
      data.dimord = ['rpttap_' data.dimord];
      data.dim    = [Nrpt   data.dim ];
    end

    % the nested trial structure is not compatible with dimord
    if isfield(data, 'trial') && isstruct(data.trial)
      param = fieldnames(data.trial);
      for i=1:length(param)
        if isa(data.trial(1).(param{i}), 'cell')
          concat = cell(data.dim(1), prod(data.dim(2:end)));
        else
          concat = zeros(data.dim(1), prod(data.dim(2:end)));
        end
        for j=1:length(data.trial)
          tmp = data.trial(j).(param{i});
          concat(j,:) = tmp(:);
        end % for each trial
        data.trial = rmfield(data.trial, param{i});
        data.(param{i}) = reshape(concat, data.dim);
      end % for each param
      data = rmfield(data, 'trial');
    end
  end

  % ensure consistent dimensions of the source reconstructed data
  % reshape each of the source reconstructed parameters
  if issource && prod(data.dim)==size(data.pos,1)
    dim = [prod(data.dim) 1];
  elseif isfield(data, 'dim'),
    dim = [data.dim 1];
  else
    %HACK
    dimtok = tokenize(data.dimord, '_');
    for i=1:length(dimtok)
      if strcmp(dimtok(i), 'pos')
        dim(1,i) = size(getsubfield(data,dimtok{i}),1);
      elseif strcmp(dimtok(i), 'rpt')
        dim(1,i) = nan;
      else
        dim(1,i) = length(getsubfield(data,dimtok{i}));
      end
    end
    i = find(isnan(dim));
    if ~isempty(i)
      n = fieldnames(data);
      for ii=1:length(n)
        numels(1,ii) = numel(getfield(data,n{ii}));
      end
      nrpt = numels./prod(dim(setdiff(1:length(dim),i)));
      nrpt = nrpt(nrpt==round(nrpt));
      dim(i) = max(nrpt);
    end
    if numel(dim)==1, dim(1,2) = 1; end;
  end

  % these fields should not be reshaped
  exclude = {'cfg' 'fwhm' 'leadfield' 'q' 'rough'};
  if ~strcmp(inside, 'logical')
    % also exclude the inside/outside from being reshaped
    exclude = cat(2, exclude, {'inside' 'outside'});
  end

  param = setdiff(parameterselection('all', data), exclude);
  for i=1:length(param)
    if any(param{i}=='.')
      % the parameter is nested in a substructure, which can have multiple elements (e.g. source.trial(1).pow, source.trial(2).pow, ...)
      % loop over the substructure array and reshape for every element
      tok  = tokenize(param{i}, '.');
      sub1 = tok{1};  % i.e. this would be 'trial'
      sub2 = tok{2};  % i.e. this would be 'pow'
      tmp1 = getfield(data, sub1);
      for j=1:numel(tmp1)
        tmp2 = getfield(tmp1(j), sub2);
        tmp2 = reshape(tmp2, dim);
        tmp1(j) = setfield(tmp1(j), sub2, tmp2);
      end
      data = setfield(data, sub1, tmp1);
    else
      tmp  = getfield(data, param{i});
      tmp  = reshape(tmp, dim);
      data = setfield(data, param{i}, tmp);
    end
  end

end

if isequal(hastrials, 'yes')
  okflag = isfield(data, 'trial');
  if ~okflag
    error('This function requires data with a ''trial'' field');
  end % if okflag
end

if isequal(hastrialdef, 'yes')
  data = fixtrialdef(data);
end

if isequal(hasoffset, 'yes')
  okflag = isfield(data, 'offset');

  if ~okflag && isfield(data, 'time') && isa(data.time, 'cell')
    if ~isfield(data, 'fsample')
      data.fsample = 1/(data.time{1}(2)-data.time{1}(1));
    end
    for i=1:length(data.time);
      data.offset(i) = time2offset(data.time{i}, data.fsample);
    end
    okflag = 1;
  elseif ~okflag && datatype(data, 'mvar')
    data.offset = 0;
    okflag = 1;
  end

  if ~okflag
    error('This function requires data with an ''offset'' field');
  end % if okflag

elseif isequal(hasoffset, 'no') && isfield(data, 'offset')
  data = rmfield(data, 'offset');
end % if hasoffset

if isequal(hascumtapcnt, 'yes') && ~isfield(data, 'cumtapcnt')
  error('This function requires data with a ''cumtapcnt'' field');
elseif isequal(hascumtapcnt, 'no') && isfield(data, 'cumtapcnt')
  data = rmfield(data, 'cumtapcnt');
end % if hascumtapcnt

if isequal(hasdof, 'yes') && ~isfield(data, 'hasdof')
  error('This function requires data with a ''dof'' field');
elseif isequal(hasdof, 'no') && isfield(data, 'hasdof')
  data = rmfield(data, 'cumtapcnt');
end % if hasdof

if ~isempty(cmbrepresentation)
  if istimelock
    data = fixcov(data, cmbrepresentation);
  elseif isfreq
    data = fixcsd(data, cmbrepresentation, channelcmb);
  elseif isfreqmvar
    data = fixcsd(data, cmbrepresentation, channelcmb);
  else
    error('This function requires data with a covariance, coherence or cross-spectrum');
  end
end % cmbrepresentation

if issource && strcmp(keepoutside, 'no'),
  % remove all grid points that are marked as outside
  data = source2sparse(data);
end

if issource && ~isempty(sourcerepresentation)
  data = fixsource(data, 'type', sourcerepresentation);
end

if issource && ~isempty(haspow)
 data = fixsource(data, 'type', sourcerepresentation, 'haspow', haspow);
end 

if isfield(data, 'grad')
  % ensure that the gradiometer balancing is specified
  if ~isfield(data.grad, 'balance') || ~isfield(data.grad.balance, 'current')
    data.grad.balance.current = 'none';
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% represent the covariance matrix in a particular manner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = fixcov(data, desired)
if isfield(data, 'cov')     && ~isfield(data, 'labelcmb')
  current = 'full';
elseif isfield(data, 'cov') &&  isfield(data, 'labelcmb')
  current = 'sparse';
else
  error('Could not determine the current representation of the covariance matrix');
end
if isequal(current, desired)
  % nothing to do
elseif strcmp(current, 'full') && strcmp(desired, 'sparse')
  % FIXME should be implemented
  error('not yet implemented');
elseif strcmp(current, 'sparse') && strcmp(desired, 'full')
  % FIXME should be implemented
  error('not yet implemented');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = comp2raw(data)
% just remove the component topographies
data = rmfield(data, 'topo');
data = rmfield(data, 'topolabel');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = volume2source(data)
if isfield(data, 'dimord')
  % it is a modern source description
else
  % it is an old-fashioned source description
  xgrid = 1:data.dim(1);
  ygrid = 1:data.dim(2);
  zgrid = 1:data.dim(3);
  [x y z] = ndgrid(xgrid, ygrid, zgrid);
  data.pos = warp_apply(data.transform, [x(:) y(:) z(:)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = source2volume(data)

if isfield(data, 'dimord')
  % it is a modern source description

  %this part depends on the assumption that the list of positions is describing a full 3D volume in
  %an ordered way which allows for the extraction of a transformation matrix
  %i.e. slice by slice
  try,
    if isfield(data, 'dim'),
      data.dim = pos2dim3d(data.pos, data.dim);
    else
      data.dim = pos2dim3d(data);
    end
  catch
  end
end

if isfield(data, 'dim') && length(data.dim)>=3,
  % it is an old-fashioned source description, or the source describes a regular 3D volume in pos
  xgrid = 1:data.dim(1);
  ygrid = 1:data.dim(2);
  zgrid = 1:data.dim(3);
  [x y z] = ndgrid(xgrid, ygrid, zgrid);
  ind = [x(:) y(:) z(:)];    % these are the positions expressed in voxel indices along each of the three axes
  pos = data.pos;            % these are the positions expressed in head coordinates
  % represent the positions in a manner that is compatible with the homogeneous matrix multiplication,
  % i.e. pos = H * ind
  ind = ind'; ind(4,:) = 1;
  pos = pos'; pos(4,:) = 1;
  % recompute the homogeneous transformation matrix
  data.transform = pos / ind;
end

% remove the unwanted fields
if isfield(data, 'pos'),    data = rmfield(data, 'pos');    end
if isfield(data, 'xgrid'),  data = rmfield(data, 'xgrid');  end
if isfield(data, 'ygrid'),  data = rmfield(data, 'ygrid');  end
if isfield(data, 'zgrid'),  data = rmfield(data, 'zgrid');  end

% make inside a volume
data = fixinside(data, 'logical');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = freq2raw(freq)

if strcmp(freq.dimord, 'rpt_chan_freq_time')
  dat = freq.powspctrm;
elseif strcmp(freq.dimord, 'rpttap_chan_freq_time')
  warning('converting fourier representation into raw data format. this is experimental code');
  dat = freq.fourierspctrm;
else
  error('this only works for dimord=''rpt_chan_freq_time''');
end

nrpt  = size(dat,1);
nchan = size(dat,2);
nfreq = size(dat,3);
ntime = size(dat,4);
data = [];
% create the channel labels like "MLP11@12Hz""
k = 0;
for i=1:nfreq
  for j=1:nchan
    k = k+1;
    data.label{k} = sprintf('%s@%dHz', freq.label{j}, freq.freq(i));
  end
end
% reshape and copy the data as if it were timecourses only
for i=1:nrpt
  data.time{i}  = freq.time;
  data.trial{i} = reshape(dat(i,:,:,:), nchan*nfreq, ntime);
  if any(isnan(data.trial{i}(1,:))),
    tmp = data.trial{i}(1,:);
    begsmp = find(isfinite(tmp),1, 'first');
    endsmp = find(isfinite(tmp),1, 'last' );
    data.trial{i} = data.trial{i}(:, begsmp:endsmp);
    data.time{i}  = data.time{i}(begsmp:endsmp);
  end
end
nsmp = cellfun('size',data.time,2);
seln = find(nsmp>1,1, 'first');
data.fsample = 1/(data.time{seln}(2)-data.time{seln}(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = raw2timelock(data)
ntrial = numel(data.trial);
nchan  = numel(data.label);
if ntrial==1
  data.time   = data.time{1};
  data.avg    = data.trial{1};
  data        = rmfield(data, 'trial');
  data.dimord = 'chan_time';
else
  % determine the location of the trials relative to the resulting combined time axis
  begtime  = cellfun(@min,  data.time);
  endtime  = cellfun(@max,  data.time);
  begsmp   = round((begtime - min(begtime)) * data.fsample) + 1;
  endsmp   = round((endtime - min(begtime)) * data.fsample) + 1;

  % create a combined time axis and concatenate all trials
  tmptime  = min(begtime):(1/data.fsample):max(endtime);
  tmptrial = zeros(ntrial, nchan, length(tmptime)) + nan;
  for i=1:ntrial
    tmptrial(i,:,begsmp(i):endsmp(i)) = data.trial{i};
  end

  % update the trial definition
  trl = findcfg(data.cfg, 'trl');
  if ~isempty(trl)
    begpad = begsmp - min(begsmp);   % number of nan-samples added to the begin
    endpad = max(endsmp) - endsmp;   % number of nan-samples added to the end
    data.cfg.trlold = trl;
    trl(:,1) = trl(:,1) - begpad(:);
    trl(:,2) = trl(:,2) + endpad(:);
    trl(:,3) = trl(:,3) - begpad(:);
    data.cfg.trl = trl;
  end

  % construct the output timelocked data
  % data.avg     = reshape(nanmean(tmptrial,     1), nchan, length(tmptime));
  % data.var     = reshape(nanvar (tmptrial, [], 1), nchan, length(tmptime))
  % data.dof     = reshape(sum(~isnan(tmptrial), 1), nchan, length(tmptime));
  data.trial   = tmptrial;
  data.time    = tmptime;
  data.dimord = 'rpt_chan_time';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = timelock2raw(data)
switch data.dimord
  case 'chan_time'
    data.trial{1} = data.avg;
    data.time     = {data.time};
    data          = rmfield(data, 'avg');
  case 'rpt_chan_time'
    tmptrial = {};
    tmptime  = {};
    ntrial = size(data.trial,1);
    nchan  = size(data.trial,2);
    ntime  = size(data.trial,3);
    for i=1:ntrial
      tmptrial{i} = reshape(data.trial(i,:,:), [nchan, ntime]);
      tmptime{i}  = data.time;
    end
    data       = rmfield(data, 'trial');
    data.trial = tmptrial;
    data.time  = tmptime;
  case 'subj_chan_time'
    tmptrial = {};
    tmptime  = {};
    ntrial = size(data.individual,1);
    nchan  = size(data.individual,2);
    ntime  = size(data.individual,3);
    for i=1:ntrial
      tmptrial{i} = reshape(data.individual(i,:,:), [nchan, ntime]);
      tmptime{i}  = data.time;
    end
    data       = rmfield(data, 'individual');
    data.trial = tmptrial;
    data.time  = tmptime;
  otherwise
    error('unsupported dimord');
end
% remove the unwanted fields
if isfield(data, 'avg'),        data = rmfield(data, 'avg'); end
if isfield(data, 'var'),        data = rmfield(data, 'var'); end
if isfield(data, 'cov'),        data = rmfield(data, 'cov'); end
if isfield(data, 'dimord'),     data = rmfield(data, 'dimord'); end
if isfield(data, 'numsamples'), data = rmfield(data, 'numsamples'); end
if isfield(data, 'dof'),        data = rmfield(data, 'dof'); end
