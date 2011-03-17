function [data] = ft_checkdata(data, varargin)

% FT_CHECKDATA checks the input data of the main FieldTrip functions, e.g. whether
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
%   [data] = ft_checkdata(data, ...)
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
%   [data] = ft_checkdata(data, 'senstype', {'ctf151', 'ctf275'}), e.g. in megrealign
%   [data] = ft_checkdata(data, 'datatype', {'timelock', 'freq'}), e.g. in sourceanalysis

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
%    You should have received a copy of the GNU General Publhasoffsetic License
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
dtype         = keyval('datatype',      varargin); % should not conflict with the ft_datatype function
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

% determine the type of input data
% this can be raw, freq, timelock, comp, spike, source, volume, dip
israw      = ft_datatype(data, 'raw');
isfreq     = ft_datatype(data, 'freq');
istimelock = ft_datatype(data, 'timelock');
iscomp     = ft_datatype(data, 'comp');
isspike    = ft_datatype(data, 'spike');
isvolume   = ft_datatype(data, 'volume');
issource   = ft_datatype(data, 'source');
isdip      = ft_datatype(data, 'dip');
ismvar     = ft_datatype(data, 'mvar');
isfreqmvar = ft_datatype(data, 'freqmvar');

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

if issource && isvolume
  % it should be either one or the other: the choice here is to
  % represent it as volume description since that is simpler to handle
  % the conversion is done by remove the grid positions
  data = rmfield(data, 'pos');
  issource = false;
end

% the ft_datatype_XXX functions ensures the consistency of the XXX datatype
% and provides a detailled description of the dataformat and its history
if     israw
  data = ft_datatype_raw(data);
elseif isfreq
  data = ft_datatype_freq(data);
elseif istimelock 
  data = ft_datatype_timelock(data);
elseif iscomp
  data = ft_datatype_comp(data);
elseif isspike
  data = ft_datatype_spike(data);
elseif isvolume
  data = ft_datatype_vol(data);
elseif issource
  data = ft_datatype_source(data);
elseif isdip
  data = ft_datatype_dip(data);
elseif ismvar || isfreqmvar
  data = ft_datatype_mvar(data);
end

if ~isempty(dtype)
  if ~isa(dtype, 'cell')
    dtype = {dtype};
  end

  okflag = 0;
  for i=1:length(dtype)
    % check that the data matches with one or more of the required ft_datatypes
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

if isequal(hasunits, 'yes') && ~isfield(data, 'units')
  % calling convert_units with only the input data adds the units without converting
  data = ft_convert_units(data);
end

if issource || isvolume,
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
  elseif issource && any(~cellfun('isempty',strfind(fieldnames(data), 'dimord')))
    dim = [size(data.pos,1) 1]; %sparsely represented source structure new style
  elseif isfield(data, 'dim'),
    dim = [data.dim 1];
  elseif isfield(data, 'dimord'),
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
    data.offset = data.offset(:); % ensure that it is a column vector
    okflag = 1;
  elseif ~okflag && ft_datatype(data, 'mvar')
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

if issource && ~isempty(sourcerepresentation)
  data = fixsource(data, 'type', sourcerepresentation);
end

if issource && ~strcmp(haspow, 'no')
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% represent the cross-spectral density matrix in a particular manner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = fixcsd(data, desired, channelcmb)

% FIXCSD converts univariate frequency domain data (fourierspctrm) into a bivariate
% representation (crsspctrm), or changes the representation of bivariate frequency
% domain data (sparse/full/sparsewithpow, sparsewithpow only works for crsspctrm or
% fourierspctrm)

% Copyright (C) 2010, Jan-Mathijs Schoffelen, Robert Oostenveld 

if isfield(data, 'crsspctrm') && isfield(data, 'powspctrm')
  current = 'sparsewithpow';
elseif isfield(data, 'powspctrm')
  current = 'sparsewithpow';
elseif isfield(data, 'fourierspctrm') && ~isfield(data, 'labelcmb')
  current = 'fourier';
elseif ~isfield(data, 'labelcmb')
  current = 'full';
elseif isfield(data, 'labelcmb')
  current = 'sparse';
else
  error('Could not determine the current representation of the %s matrix', param);
end

% first go from univariate fourier to the required bivariate representation
if strcmp(current, 'fourier') && strcmp(desired, 'fourier')
  % nothing to do
elseif strcmp(current, 'fourier') && strcmp(desired, 'sparsewithpow')
  dimtok = tokenize(data.dimord, '_');
  if ~isempty(strmatch('rpttap',   dimtok)),
    nrpt = length(data.cumtapcnt);
    flag = 0;
  else
    nrpt = 1;
    flag = 1;
  end
  if ~isempty(strmatch('freq',  dimtok)), nfrq=length(data.freq);      else nfrq = 1; end
  if ~isempty(strmatch('time',  dimtok)), ntim=length(data.time);      else ntim = 1; end
  
  fastflag = all(data.cumtapcnt(:)==data.cumtapcnt(1));

  %create auto-spectra
  nchan     = length(data.label);
  if fastflag
    % all trials have the same amount of tapers
    powspctrm = zeros(nrpt,nchan,nfrq,ntim);
    ntap      = data.cumtapcnt(1);
    for p = 1:ntap
      powspctrm = powspctrm + abs(data.fourierspctrm(p:ntap:end,:,:,:,:)).^2;
    end
    powspctrm = powspctrm./ntap;
  else
    % different amount of tapers
    powspctrm = zeros(nrpt,nchan,nfrq,ntim)+i.*zeros(nrpt,nchan,nfrq,ntim);
    sumtapcnt = [0;cumsum(data.cumtapcnt(:))];
    for p = 1:nrpt
      indx   = (sumtapcnt(p)+1):sumtapcnt(p+1);
      tmpdat = data.fourierspctrm(indx,:,:,:);
      powspctrm(p,:,:,:) = (sum(tmpdat.*conj(tmpdat),1))./data.cumtapcnt(p);
    end
  end

  %create cross-spectra
  if ~isempty(channelcmb),
    ncmb      = size(channelcmb,1);
    cmbindx   = zeros(ncmb,2);
    labelcmb  = cell(ncmb,2);
    for k = 1:ncmb
      ch1 = find(strcmp(data.label, channelcmb(k,1)));
      ch2 = find(strcmp(data.label, channelcmb(k,2)));
      if ~isempty(ch1) && ~isempty(ch2),
        cmbindx(k,:)  = [ch1 ch2];
        labelcmb(k,:) = data.label([ch1 ch2])';
      end
    end

    crsspctrm = zeros(nrpt,ncmb,nfrq,ntim)+i.*zeros(nrpt,ncmb,nfrq,ntim);
    if fastflag
      for p = 1:ntap
        tmpdat1   = data.fourierspctrm(p:ntap:end,cmbindx(:,1),:,:,:);
        tmpdat2   = data.fourierspctrm(p:ntap:end,cmbindx(:,2),:,:,:);
        crsspctrm = crsspctrm + tmpdat1.*conj(tmpdat2);
      end
      crsspctrm = crsspctrm./ntap;
    else
      for p = 1:nrpt
        indx    = (sumtapcnt(p)+1):sumtapcnt(p+1);
        tmpdat1 = data.fourierspctrm(indx,cmbindx(:,1),:,:);
        tmpdat2 = data.fourierspctrm(indx,cmbindx(:,2),:,:);
        crsspctrm(p,:,:,:) = (sum(tmpdat1.*conj(tmpdat2),1))./data.cumtapcnt(p);
      end
    end
    data.crsspctrm = crsspctrm;
    data.labelcmb  = labelcmb;
  end
  data.powspctrm = powspctrm;
  data           = rmfield(data, 'fourierspctrm');
  if ntim>1,
    data.dimord = 'chan_freq_time';
  else
    data.dimord = 'chan_freq';
  end
  
  if nrpt>1,
    data.dimord = ['rpt_',data.dimord];
  end
  
  if flag, siz = size(data.crsspctrm); data.crsspctrm = reshape(data.crsspctrm, siz(2:end)); end
elseif strcmp(current, 'fourier') && strcmp(desired, 'sparse')

  if isempty(channelcmb), error('no channel combinations are specified'); end
  dimtok = tokenize(data.dimord, '_');
  if ~isempty(strmatch('rpttap',   dimtok)),
    nrpt = length(data.cumtapcnt);
    flag = 0;
  else
    nrpt = 1;
    flag = 1;
  end
  if ~isempty(strmatch('freq',  dimtok)), nfrq=length(data.freq); else nfrq = 1; end
  if ~isempty(strmatch('time',  dimtok)), ntim=length(data.time); else ntim = 1; end
  
  ncmb      = size(channelcmb,1);
  cmbindx   = zeros(ncmb,2);
  labelcmb  = cell(ncmb,2);
  for k = 1:ncmb
    ch1 = find(strcmp(data.label, channelcmb(k,1)));
    ch2 = find(strcmp(data.label, channelcmb(k,2)));
    if ~isempty(ch1) && ~isempty(ch2),
      cmbindx(k,:)  = [ch1 ch2];
      labelcmb(k,:) = data.label([ch1 ch2])';
    end
  end

  sumtapcnt = [0;cumsum(data.cumtapcnt(:))];
  fastflag  = all(data.cumtapcnt(:)==data.cumtapcnt(1));
  
  if fastflag && nrpt>1
    ntap = data.cumtapcnt(1);
    
    % compute running sum across tapers
    siz = [size(data.fourierspctrm) 1];
    
    for p = 1:ntap
      indx      = p:ntap:nrpt*ntap;
      
      if p==1.

        tmpc = zeros(numel(indx), size(cmbindx,1), siz(3), siz(4)) + ... 
           1i.*zeros(numel(indx), size(cmbindx,1), siz(3), siz(4)); 
      end
      
      for k = 1:size(cmbindx,1)
        tmpc(:,k,:,:) = data.fourierspctrm(indx,cmbindx(k,1),:,:).*  ...
                   conj(data.fourierspctrm(indx,cmbindx(k,2),:,:));
      end

      if p==1
        crsspctrm = tmpc;
      else
        crsspctrm = tmpc + crsspctrm;
      end
    end
    crsspctrm = crsspctrm./ntap;
  else
    crsspctrm = zeros(nrpt, ncmb, nfrq, ntim);
    for p = 1:nrpt
      indx    = (sumtapcnt(p)+1):sumtapcnt(p+1);
      tmpdat1 = data.fourierspctrm(indx,cmbindx(:,1),:,:);
      tmpdat2 = data.fourierspctrm(indx,cmbindx(:,2),:,:);
      crsspctrm(p,:,:,:) = (sum(tmpdat1.*conj(tmpdat2),1))./data.cumtapcnt(p);
    end
  end
  data.crsspctrm = crsspctrm;
  data.labelcmb  = labelcmb;
  data           = rmfield(data, 'fourierspctrm');
  data           = rmfield(data, 'label');
  if ntim>1,
    data.dimord = 'chan_freq_time';
  else
    data.dimord = 'chan_freq';
  end
  
  if nrpt>1,
    data.dimord = ['rpt_',data.dimord];
  end

  if flag, siz = size(data.crsspctrm); data.crsspctrm = reshape(data.crsspctrm, siz(2:end)); end
elseif strcmp(current, 'fourier') && strcmp(desired, 'full')

  % this is how it is currently and the desired functionality of prepare_freq_matrices
  dimtok = tokenize(data.dimord, '_');
  if ~isempty(strmatch('rpttap',   dimtok)),
    nrpt = size(data.cumtapcnt, 1);
    flag = 0;
  else
    nrpt = 1;
    flag = 1;
  end
  if ~isempty(strmatch('rpttap',dimtok)), nrpt=size(data.cumtapcnt, 1); else nrpt = 1; end
  if ~isempty(strmatch('freq',  dimtok)), nfrq=length(data.freq);       else nfrq = 1; end
  if ~isempty(strmatch('time',  dimtok)), ntim=length(data.time);       else ntim = 1; end
  if any(data.cumtapcnt(1,:) ~= data.cumtapcnt(1,1)), error('this only works when all frequencies have the same number of tapers'); end
  nchan     = length(data.label);
  crsspctrm = zeros(nrpt,nchan,nchan,nfrq,ntim);
  sumtapcnt = [0;cumsum(data.cumtapcnt(:,1))];
  for k = 1:ntim
    for m = 1:nfrq
      for p = 1:nrpt
        %FIXME speed this up in the case that all trials have equal number of tapers
        indx   = (sumtapcnt(p)+1):sumtapcnt(p+1);
        tmpdat = transpose(data.fourierspctrm(indx,:,m,k));
        crsspctrm(p,:,:,m,k) = (tmpdat*tmpdat')./data.cumtapcnt(p);
        clear tmpdat;
      end
    end
  end
  data.crsspctrm = crsspctrm;
  data           = rmfield(data, 'fourierspctrm');
  
  if ntim>1,
    data.dimord = 'chan_chan_freq_time';
  else
    data.dimord = 'chan_chan_freq';
  end
  
  if nrpt>1,
    data.dimord = ['rpt_',data.dimord];
  end   
  
  % remove first singleton dimension
  if flag || nrpt==1, siz = size(data.crsspctrm); data.crsspctrm = reshape(data.crsspctrm, siz(2:end)); end

elseif strcmp(current, 'fourier') && strcmp(desired, 'fullfast'),

  dimtok = tokenize(data.dimord, '_');
  nrpt = size(data.fourierspctrm, 1);    
  nchn = numel(data.label);    
  nfrq = numel(data.freq);  
  if ~isempty(strmatch('time',  dimtok)), ntim=numel(data.time); else ntim = 1; end
  
  data.fourierspctrm = reshape(data.fourierspctrm, [nrpt nchn nfrq*ntim]);
  data.fourierspctrm(~isfinite(data.fourierspctrm)) = 0;
  crsspctrm = complex(zeros(nchn,nchn,nfrq*ntim));
  for k = 1:nfrq*ntim
    tmp = transpose(data.fourierspctrm(:,:,k));
    n   = sum(tmp~=0,2);
    crsspctrm(:,:,k) = tmp*tmp'./n(1);
  end
  data           = rmfield(data, 'fourierspctrm');
  data.crsspctrm = reshape(crsspctrm, [nchn nchn nfrq ntim]);
  if isfield(data, 'time'),
    data.dimord = 'chan_chan_freq_time';
  else
    data.dimord = 'chan_chan_freq';
  end

end % convert to the requested bivariate representation

% from one bivariate representation to another
if isequal(current, desired)
  % nothing to do

elseif (strcmp(current, 'full')       && strcmp(desired, 'fourier')) || ...
    (strcmp(current, 'sparse')        && strcmp(desired, 'fourier')) || ...
    (strcmp(current, 'sparsewithpow') && strcmp(desired, 'fourier'))
  % this is not possible
  error('converting the cross-spectrum into a Fourier representation is not possible');

elseif strcmp(current, 'full') && strcmp(desired, 'sparsewithpow')
  error('not yet implemented');
elseif strcmp(current, 'sparse') && strcmp(desired, 'sparsewithpow')
  % convert back to crsspctrm/powspctrm representation: useful for plotting functions etc
  indx     = labelcmb2indx(data.labelcmb);
  autoindx = indx(indx(:,1)==indx(:,2), 1);
  cmbindx  = setdiff([1:size(indx,1)]', autoindx);
  
  if strcmp(data.dimord(1:3), 'rpt')
    data.powspctrm = data.crsspctrm(:, autoindx, :, :);
    data.crsspctrm = data.crsspctrm(:, cmbindx,  :, :);
  else
    data.powspctrm = data.crsspctrm(autoindx, :, :);
    data.crsspctrm = data.crsspctrm(cmbindx,  :, :);
  end 
  data.label    = data.labelcmb(autoindx,1);
  data.labelcmb = data.labelcmb(cmbindx, :);
  
  if isempty(cmbindx)
    data = rmfield(data, 'crsspctrm');
    data = rmfield(data, 'labelcmb');
  end
  
elseif strcmp(current, 'full') && strcmp(desired, 'sparse')
  dimtok = tokenize(data.dimord, '_');
  if ~isempty(strmatch('rpt',   dimtok)), nrpt=numel(data.cumtapcnt); else nrpt = 1; end
  if ~isempty(strmatch('freq',  dimtok)), nfrq=numel(data.freq);      else nfrq = 1; end
  if ~isempty(strmatch('time',  dimtok)), ntim=numel(data.time);      else ntim = 1; end
  nchan    = length(data.label);
  ncmb     = nchan*nchan;
  labelcmb = cell(ncmb, 2);
  cmbindx  = zeros(nchan, nchan);
  k = 1;
  for j=1:nchan
    for m=1:nchan
      labelcmb{k, 1} = data.label{m};
      labelcmb{k, 2} = data.label{j};
      cmbindx(m,j)   = k;
      k = k+1;
    end
  end
  
  % reshape all possible fields
  fn = fieldnames(data);
  for ii=1:numel(fn)
    if numel(data.(fn{ii})) == nrpt*ncmb*nfrq*ntim;
      if nrpt>1,
        data.(fn{ii}) = reshape(data.(fn{ii}), nrpt, ncmb, nfrq, ntim);
      else
        data.(fn{ii}) = reshape(data.(fn{ii}), ncmb, nfrq, ntim);
      end
    end
  end
  % remove obsolete fields
  data           = rmfield(data, 'label');
  try, data      = rmfield(data, 'dof'); end
  % replace updated fields
  data.labelcmb  = labelcmb;
  if ntim>1,
    data.dimord = 'chancmb_freq_time';
  else
    data.dimord = 'chancmb_freq';
  end

  if nrpt>1,
    data.dimord = ['rpt_',data.dimord];
  end

elseif strcmp(current, 'sparsewithpow') && strcmp(desired, 'sparse')

  % this representation for sparse data contains autospectra
  % as e.g. {'A' 'A'} in labelcmb
  if isfield(data, 'crsspctrm'),
    dimtok         = tokenize(data.dimord, '_');
    catdim         = match_str(dimtok, {'chan' 'chancmb'});
    data.crsspctrm = cat(catdim, data.powspctrm, data.crsspctrm);
    data.labelcmb  = [data.label(:) data.label(:); data.labelcmb];
    data           = rmfield(data, 'powspctrm');
  else
    data.crsspctrm = data.powspctrm;
    data.labelcmb  = [data.label(:) data.label(:)];
    data           = rmfield(data, 'powspctrm');
  end
  data = rmfield(data, 'label');

elseif strcmp(current, 'sparse') && strcmp(desired, 'full')
  dimtok = tokenize(data.dimord, '_');
  if ~isempty(strmatch('rpt',   dimtok)), nrpt=numel(data.cumtapcnt); else nrpt = 1; end
  if ~isempty(strmatch('freq',  dimtok)), nfrq=numel(data.freq);      else nfrq = 1; end
  if ~isempty(strmatch('time',  dimtok)), ntim=numel(data.time);      else ntim = 1; end
  
  if ~isfield(data, 'label')
    data.label = unique(data.labelcmb(:));
  end

  nchan     = length(data.label);
  ncmb      = size(data.labelcmb,1);
  cmbindx   = zeros(nchan,nchan);

  for k = 1:size(data.labelcmb,1)
    ch1 = find(strcmp(data.label, data.labelcmb(k,1)));
    ch2 = find(strcmp(data.label, data.labelcmb(k,2)));
    if ~isempty(ch1) && ~isempty(ch2),
      cmbindx(ch1,ch2) = k;
    end
  end

  complete = all(cmbindx(:)~=0);

  fn = fieldnames(data);
  for ii=1:numel(fn)
    if numel(data.(fn{ii})) == nrpt*ncmb*nfrq*ntim;
      if nrpt==1,
        data.(fn{ii}) = reshape(data.(fn{ii}), [nrpt ncmb nfrq ntim]);
      end

      tmpall = nan(nrpt,nchan,nchan,nfrq,ntim);

      for j = 1:nrpt
        for k = 1:ntim
          for m = 1:nfrq
            tmpdat = nan(nchan,nchan);
            indx   = find(cmbindx);
            if ~complete
              % this realizes the missing combinations to be represented as the
              % conjugate of the corresponding combination across the diagonal
              tmpdat(indx) = reshape(data.(fn{ii})(j,cmbindx(indx),m,k),[numel(indx) 1]);
              tmpdat       = ctranspose(tmpdat);
            end
            tmpdat(indx)    = reshape(data.(fn{ii})(j,cmbindx(indx),m,k),[numel(indx) 1]);
            tmpall(j,:,:,m,k) = tmpdat;
          end % for m
        end % for k
      end % for j

      % replace the data in the old representation with the new representation
      if nrpt>1,
        data.(fn{ii}) = tmpall;
      else
        data.(fn{ii}) = reshape(tmpall, [nchan nchan nfrq ntim]);
      end
    end % if numel
  end % for ii

  % remove obsolete fields
  try, data      = rmfield(data, 'powspctrm');  end
  try, data      = rmfield(data, 'labelcmb');   end
  try, data      = rmfield(data, 'dof');        end

  if ntim>1,
    data.dimord = 'chan_chan_freq_time';
  else
    data.dimord = 'chan_chan_freq';
  end

  if nrpt>1,
    data.dimord = ['rpt_',data.dimord];
  end

elseif strcmp(current, 'sparse') && strcmp(desired, 'fullfast')
  dimtok = tokenize(data.dimord, '_');
  if ~isempty(strmatch('rpt',   dimtok)), nrpt=numel(data.cumtapcnt); else nrpt = 1; end
  if ~isempty(strmatch('freq',  dimtok)), nfrq=numel(data.freq);      else nfrq = 1; end
  if ~isempty(strmatch('time',  dimtok)), ntim=numel(data.time);      else ntim = 1; end
  
  if ~isfield(data, 'label')
    data.label = unique(data.labelcmb(:));
  end

  nchan     = length(data.label);
  ncmb      = size(data.labelcmb,1);
  cmbindx   = zeros(nchan,nchan);

  for k = 1:size(data.labelcmb,1)
    ch1 = find(strcmp(data.label, data.labelcmb(k,1)));
    ch2 = find(strcmp(data.label, data.labelcmb(k,2)));
    if ~isempty(ch1) && ~isempty(ch2),
      cmbindx(ch1,ch2) = k;
    end
  end

  complete = all(cmbindx(:)~=0);

  fn = fieldnames(data);
  for ii=1:numel(fn)
    if numel(data.(fn{ii})) == nrpt*ncmb*nfrq*ntim;
      if nrpt==1,
        data.(fn{ii}) = reshape(data.(fn{ii}), [nrpt ncmb nfrq ntim]);
      end

      tmpall = nan(nchan,nchan,nfrq,ntim);

      for k = 1:ntim
        for m = 1:nfrq
          tmpdat = nan(nchan,nchan);
          indx   = find(cmbindx);
          if ~complete
            % this realizes the missing combinations to be represented as the
            % conjugate of the corresponding combination across the diagonal
            tmpdat(indx) = reshape(nanmean(data.(fn{ii})(:,cmbindx(indx),m,k)),[numel(indx) 1]);
            tmpdat       = ctranspose(tmpdat);
          end
          tmpdat(indx)    = reshape(nanmean(data.(fn{ii})(:,cmbindx(indx),m,k)),[numel(indx) 1]);
          tmpall(:,:,m,k) = tmpdat;
        end % for m
      end % for k

      % replace the data in the old representation with the new representation
      if nrpt>1,
        data.(fn{ii}) = tmpall;
      else
        data.(fn{ii}) = reshape(tmpall, [nchan nchan nfrq ntim]);
      end
    end % if numel
  end % for ii

  % remove obsolete fields
  try, data      = rmfield(data, 'powspctrm');  end
  try, data      = rmfield(data, 'labelcmb');   end
  try, data      = rmfield(data, 'dof');        end

  if ntim>1,
    data.dimord = 'chan_chan_freq_time';
  else
    data.dimord = 'chan_chan_freq';
  end

elseif strcmp(current, 'sparsewithpow') && strcmp(desired, 'full')
  % this is how is currently done in prepare_freq_matrices
  data = ft_checkdata(data, 'cmbrepresentation', 'sparse');
  data = ft_checkdata(data, 'cmbrepresentation', 'full');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert to new source representation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output] = fixsource(input, varargin)

% FIXSOURCE converts old style source structures into new style source structures and the
% other way around
%
% Use as:
%   fixsource(input, type)
%    where input is a source structure,
%
% Typically, old style source structures contain
%   avg.XXX or trial.XXX fields
%
% The new style source structure contains:
%   source.pos
%   source.dim (optional, if the list of positions describes a 3D volume
%   source.XXX the old style subfields in avg/trial
%   source.XXXdimord string how to interpret the respective XXX field:
%     e.g. source.leadfield = cell(1,Npos), source.leadfielddimord = '{pos}_chan_ori'
%          source.mom       = cell(1,Npos), source.momdimord       = '{pos}_ori_rpttap'

type   = keyval('type',   varargin{:});
haspow = keyval('haspow', varargin{:});

if isempty(type),   type   = 'old'; end
if isempty(haspow), haspow = 'no';  end

fnames = fieldnames(input);
tmp    = cell2mat(strfind(fnames, 'dimord')); %get dimord like fields
if any(tmp>1),
  current = 'new';
elseif any(tmp==1),
  %don't know what to do yet data is JM's own invention
  current = 'old';
else
  current = 'old';
end

if strcmp(current, type),
  %do nothing
  output = input;
  
  %return
elseif strcmp(current, 'old') && strcmp(type, 'new'),
  %go from old to new

  if isfield(input, 'avg'),
    stuff  = getfield(input, 'avg');
    output = rmfield(input,  'avg');
  elseif isfield(input, 'trial'),
    stuff  = getfield(input, 'trial');
    output = rmfield(input,  'trial');
  else
    %this could occur later in the pipeline, e.g. when doing group statistics using individual subject
    %descriptive statistics
    error('the input does not contain an avg or trial field');
  end

  %-------------------------------------------------  
  %remove and rename the specified fields if present
  removefields = {'xgrid';'ygrid';'zgrid';'method'};
  renamefields = {'frequency' 'freq'; 'csdlabel' 'orilabel'};
  fnames       = fieldnames(output);
  for k = 1:numel(fnames)
    ix = strmatch(fnames{k}, removefields);
    if ~isempty(ix),
      output = rmfield(output, fnames{k});
    end
    ix = strmatch(fnames{k}, renamefields(:,1), 'exact');
    if ~isempty(ix),
      output = setfield(output, renamefields{ix,2}, ...
                        getfield(output, renamefields{ix,1}));
      output = rmfield(output, fnames{k});
    end
  end

  %----------------------------------------------------------------------
  %put the stuff originally in avg or trial one level up in the structure
  fnames       = fieldnames(stuff(1));
  npos         = size(input.pos,1);
  nrpt         = numel(stuff);
  for k = 1:numel(fnames)
    if nrpt>1,
      %multiple trials
      %(or subjects FIXME not yet implemented, nor tested)
      tmp  = getfield(stuff(1), fnames{k});
      siz  = size(tmp);
      if isfield(input, 'cumtapcnt') && strcmp(fnames{k}, 'mom')
        %pcc based mom is orixrpttap
        %tranpose to keep manageable
        for kk = 1:numel(input.inside)
          indx = input.inside(kk);
          tmp{indx} = permute(tmp{indx}, [2 1 3]); 
        end
        nrpttap = sum(input.cumtapcnt);
        sizvox  = [size(tmp{input.inside(1)}) 1];
        sizvox  = [nrpttap sizvox(2:end)];
      elseif strcmp(fnames{k}, 'mom'),
        %this is then probably not a frequency based mom
        nrpttap = numel(stuff);
        sizvox  = [size(tmp{input.inside(1)}) 1];
        sizvox  = [nrpttap sizvox];
      elseif iscell(tmp)
        nrpttap = numel(stuff);
        sizvox  = [size(tmp{input.inside(1)}) 1];
        sizvox  = [nrpttap sizvox];
      end
      
      if siz(1) ~= npos && siz(2) ==npos,
        tmp = transpose(tmp);
      end
      
      if iscell(tmp)
        %allocate memory for cell-array
        tmpall = cell(npos,1);
        for n = 1:numel(input.inside)
          tmpall{input.inside(n)} = zeros(sizvox);
        end
      else
        %allocate memory for matrix
        tmpall = zeros([npos nrpt siz(2:end)]);
      end
      
      cnt = 0;
      for m = 1:nrpt
        tmp = getfield(stuff(m), fnames{k});
        siz = size(tmp);
        if siz(1) ~= npos && siz(2) ==npos,
          tmp = transpose(tmp);
        end
        
        if ~iscell(tmp),
          tmpall(:,m,:,:,:) = tmp;
        else
          for n = 1:numel(input.inside)
            indx   = input.inside(n);
            tmpdat = tmp{indx};
            if isfield(input, 'cumtapcnt') && strcmp(fnames{k}, 'mom'),
              if n==1, siz1 = size(tmpdat,2); end
            else
              if n==1, siz1 = 1; end
            end
            tmpall{indx}(cnt+[1:siz1],:,:,:,:) = tmpdat;
            if n==numel(input.inside), cnt  = cnt + siz1;     end
          end
        end
      end
      output    = setfield(output, fnames{k}, tmpall);
      newdimord = createdimord(output, fnames{k}, 1);
      if ~isempty(newdimord)
        output    = setfield(output, [fnames{k},'dimord'], newdimord);
      end
    
    else
      tmp = getfield(stuff, fnames{k});
      siz = size(tmp);
      if isfield(input, 'cumtapcnt') && strcmp(fnames{k}, 'mom')
        %pcc based mom is orixrpttap
        %tranpose to keep manageable
        for kk = 1:numel(input.inside)
          indx = input.inside(kk);
          tmp{indx} = permute(tmp{indx}, [2 1 3]); 
        end
      end
      if siz(1) ~= npos && siz(2) ==npos,
        tmp = transpose(tmp);
      end
      output    = setfield(output, fnames{k}, tmp);
      newdimord = createdimord(output, fnames{k}); 
      if ~isempty(newdimord)
        output    = setfield(output, [fnames{k},'dimord'], newdimord);
      end
    end
  end
  
  if isfield(output, 'csdlabel')
    output = setfield(output, 'orilabel', getfield(output, 'csdlabel'));
    output = rmfield(output,  'csdlabel');
  end

  if isfield(output, 'leadfield')
    % add dimord to leadfield as well. since the leadfield is not in
    % the original .avg or .trial field it has not yet been taken care of
    output.leadfielddimord = createdimord(output, 'leadfield');  
  end
  
  if isfield(output, 'ori')
    % convert cell-array ori into matrix
    ori = zeros(3,npos) + nan;
    try,
      ori(:,output.inside) = cat(2, output.ori{output.inside});
    catch
      %when oris are in wrong orientation (row rather than column)
      for k = 1:numel(output.inside)
        ori(:,output.inside(k)) = output.ori{output.inside(k)}';
      end
     end
    output.ori = ori;
  end
  current = 'new';
 
elseif strcmp(current, 'new') && strcmp(type, 'old')
  %go from new to old
  error('not implemented yet');
end

if strcmp(current, 'new') && strcmp(haspow, 'yes'), 

  %----------------------------------------------
  %convert mom into pow if requested and possible
  convert = 0;
  if isfield(output, 'mom') && size(output.mom{output.inside(1)},2)==1,
    convert = 1;
  else
    warning('conversion from mom to pow is not possible, either because there is no mom in the data, or because the dimension of mom>1. in that case call ft_sourcedescriptives first with cfg.projectmom');
  end
   
  if isfield(output, 'cumtapcnt')
    convert = 1 & convert;
  else
    warning('conversion from mom to pow will not be done, because cumtapcnt is missing');
  end
    
  if convert,  
    npos = size(output.pos,1);
    nrpt = numel(output.cumtapcnt);
    tmpmom = cat(2,output.mom{output.inside});
    tmppow = zeros(npos, nrpt);
    tapcnt = [0;cumsum(output.cumtapcnt(:))];
    for k = 1:nrpt
      ntap = tapcnt(k+1)-tapcnt(k);  
      tmppow(output.inside,k) = sum(abs(tmpmom((tapcnt(k)+1):tapcnt(k+1),:)).^2,1)./ntap;
    end
    output.pow       = tmppow;
    output.powdimord = ['pos_rpt_freq']; 
  end  

elseif strcmp(current, 'old') && strcmp(haspow, 'yes')
  warning('construction of single trial power estimates is not implemented here using old style source representation');

end
 

%--------------------------------------------------------
function [dimord] = createdimord(output, fname, rptflag);

if nargin==2, rptflag = 0; end

tmp = getfield(output, fname);

dimord = '';
dimnum = 1;
hasori = isfield(output, 'ori'); %if not, this is probably singleton and not relevant at the end

if iscell(tmp) && (size(output.pos,1)==size(tmp,dimnum) || size(output.pos,1)==size(tmp,2))
  dimord = [dimord,'{pos}'];
  dimnum = dimnum + 1;
elseif ~iscell(tmp) && size(output.pos,1)==size(tmp,dimnum)
  dimord = [dimord,'pos'];
  dimnum = dimnum + 1;
end

switch fname
  case 'cov'
    if hasori, dimord = [dimord,'_ori_ori']; end;
  case 'csd'
    if hasori, dimord = [dimord,'_ori_ori']; end;
  case 'csdlabel'
    dimord = dimord;
  case 'filter'
    dimord = [dimord,'_ori_chan']; 
  case 'leadfield'
    %if hasori,
      dimord = [dimord,'_chan_ori'];
    %else
    %  dimord = [dimord,'_chan'];
    %end
  case 'mom'
    if isfield(output, 'cumtapcnt') && sum(output.cumtapcnt)==size(tmp{output.inside(1)},1)
      if hasori,
        dimord = [dimord,'_rpttap_ori'];
      else
        dimord = [dimord,'_rpttap'];
      end
    elseif isfield(output, 'time')
      if rptflag,
        dimord = [dimord,'_rpt'];
        dimnum = dimnum + 1;
      end
      if numel(output.time)==size(tmp{output.inside(1)},dimnum)
        dimord = [dimord,'_ori_time'];
      end
    end
    
    if isfield(output, 'freq') && numel(output.freq)>1,
      dimord = [dimord,'_freq'];
    end    
  case 'nai'
    if isfield(output, 'freq') && numel(output.freq)==size(tmp,dimnum)
      dimord = [dimord,'_freq'];
    end
  case 'noise'
    if isfield(output, 'freq') && numel(output.freq)==size(tmp,dimnum)
      dimord = [dimord,'_freq'];
    end
  case 'noisecsd'
    if hasori, dimord = [dimord,'_ori_ori']; end
  case 'ori'
    dimord = '';
  case 'pow'
    if isfield(output, 'cumtapcnt') && numel(output.cumtapcnt)==size(tmp,dimnum)
      dimord = [dimord,'_rpt'];
      dimnum = dimnum + 1;
    end
    
    if isfield(output, 'freq') && numel(output.freq)>1 && numel(output.freq)==size(tmp,dimnum)
      dimord = [dimord,'_freq'];
      dimnum = dimnum+1;
    end
    
    if isfield(output, 'time') && numel(output.time)>1 && numel(output.time)==size(tmp,dimnum)
      dimord = [dimord,'_time'];
      dimnum = dimnum+1;
    end
    
    otherwise
      warning('skipping unknown fieldname %s', fname);
      %error(sprintf('unknown fieldname %s', fname));
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

if isfield(freq, 'trialinfo'), data.trialinfo = freq.trialinfo; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = raw2timelock(data)

nsmp = cellfun('size',data.time,2);
data   = ft_checkdata(data, 'hastrialdef', 'yes');
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

  % update the sampleinfo
  begpad = begsmp - min(begsmp);
  endpad = max(endsmp) - endsmp;
  if isfield(data, 'sampleinfo')
    data.sampleinfo = data.sampleinfo + [-begpad(:) endpad(:)];
  end
  
  % construct the output timelocked data
  % data.avg     = reshape(nanmean(tmptrial,     1), nchan, length(tmptime));
  % data.var     = reshape(nanvar (tmptrial, [], 1), nchan, length(tmptime))
  % data.dof     = reshape(sum(~isnan(tmptrial), 1), nchan, length(tmptime));
  data.trial   = tmptrial;
  data.time    = tmptime;
  data.dimord = 'rpt_chan_time';
  data = rmfield(data, 'fsample'); % fsample in timelock data is obsolete
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = timelock2raw(data)
try
  nsmp = cellfun('size',data.time,2);
catch
  nsmp = size(data.time,2);
end
switch data.dimord
  case 'chan_time'
    data.trial{1} = data.avg;
    data.time     = {data.time};
    data          = rmfield(data, 'avg');
    seln = find(nsmp>1,1, 'first');
    data.fsample = 1/(data.time{seln}(2)-data.time{seln}(1));
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
    seln = find(nsmp>1,1, 'first');
    data.fsample = 1/(data.time{seln}(2)-data.time{seln}(1));
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
    seln = find(nsmp>1,1, 'first');
    data.fsample = 1/(data.time{seln}(2)-data.time{seln}(1));    
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
