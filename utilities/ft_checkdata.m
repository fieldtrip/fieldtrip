function [data] = ft_checkdata(data, varargin)

% FT_CHECKDATA checks the input data of the main FieldTrip functions, e.g. whether the
% type of data strucure corresponds with the required data. If necessary and possible,
% this function will adjust the data structure to the input requirements (e.g. change
% dimord, average over trials, convert inside from index into logical).
%
% If the input data does NOT correspond to the requirements, this function will give a
% warning message and if applicable point the user to external documentation (link to
% website).
%
% Use as
%   [data] = ft_checkdata(data, ...)
%
% Optional input arguments should be specified as key-value pairs and can include
%   feedback           = yes, no
%   datatype           = raw, freq, timelock, comp, spike, source, mesh, dip, volume, segmentation, parcellation
%   dimord             = any combination of time, freq, chan, refchan, rpt, subj, chancmb, rpttap, pos
%   senstype           = ctf151, ctf275, ctf151_planar, ctf275_planar, neuromag122, neuromag306, bti148, bti248, bti248_planar, magnetometer, electrode
%   inside             = logical, index
%   ismeg              = yes, no
%   iseeg              = yes, no
%   isnirs             = yes, no
%   hasunit            = yes, no
%   hascoordsys        = yes, no
%   haschantype        = yes, no
%   haschanunit        = yes, no
%   hassampleinfo      = yes, no, ifmakessense (applies to raw and timelock data)
%   hascumtapcnt       = yes, no (only applies to freq data)
%   hasdim             = yes, no
%   hasdof             = yes, no
%   cmbrepresentation  = sparse, full (applies to covariance and cross-spectral density)
%   fsample            = sampling frequency to use to go from SPIKE to RAW representation
%   segmentationstyle  = indexed, probabilistic (only applies to segmentation)
%   parcellationstyle  = indexed, probabilistic (only applies to parcellation)
%   hasbrain           = yes, no (only applies to segmentation)
%   trialinfostyle     = matrix, table or empty
%
% For some options you can specify multiple values, e.g.
%   [data] = ft_checkdata(data, 'senstype', {'ctf151', 'ctf275'}), e.g. in megrealign
%   [data] = ft_checkdata(data, 'datatype', {'timelock', 'freq'}), e.g. in sourceanalysis
%
% See also FT_DATATYPE_XXX for each of the respective data types.

% Copyright (C) 2007-2015, Robert Oostenveld
% Copyright (C) 2010-2012, Martin Vinck
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

% in case of an error this function could use dbstack for more detailled
% user feedback
%
% this function should replace/encapsulate
%   fixdimord
%   fixinside
%   fixprecision
%   fixvolume
%   fixpos
%   data2raw
%   raw2data
%   grid2transform
%   transform2grid
%   fourier2crsspctrm
%   freq2cumtapcnt
%   sensortype
%   time2offset
%   offset2time
%   fixsens -> this is kept a separate function because it should also be
%              called from other modules
%
% other potential uses for this function:
%   time -> offset in freqanalysis
%   average over trials
%   csd as matrix

% FIXME the following is difficult, if not impossible, to support without knowing the parameter
% FIXME it is presently (dec 2014) not being used anywhere in FT, so can be removed
%   hastrials          = yes, no

% get the optional input arguments
feedback             = ft_getopt(varargin, 'feedback', 'no');
dtype                = ft_getopt(varargin, 'datatype'); % should not conflict with the ft_datatype function
dimord               = ft_getopt(varargin, 'dimord');
stype                = ft_getopt(varargin, 'senstype'); % senstype is a function name which should not be masked
ismeg                = ft_getopt(varargin, 'ismeg');
iseeg                = ft_getopt(varargin, 'iseeg');
isnirs               = ft_getopt(varargin, 'isnirs');
inside               = ft_getopt(varargin, 'inside'); % can be 'logical' or 'index'
hastrials            = ft_getopt(varargin, 'hastrials');
hasunit              = ft_getopt(varargin, 'hasunit', 'no');
hascoordsys          = ft_getopt(varargin, 'hascoordsys', 'no');
haschantype          = ft_getopt(varargin, 'haschantype', 'no');
haschanunit          = ft_getopt(varargin, 'haschanunit', 'no');
hassampleinfo        = ft_getopt(varargin, 'hassampleinfo', 'ifmakessense');
hasdim               = ft_getopt(varargin, 'hasdim');
hascumtapcnt         = ft_getopt(varargin, 'hascumtapcnt');
hasdof               = ft_getopt(varargin, 'hasdof');
cmbrepresentation    = ft_getopt(varargin, 'cmbrepresentation');
channelcmb           = ft_getopt(varargin, 'channelcmb');
fsample              = ft_getopt(varargin, 'fsample');
segmentationstyle    = ft_getopt(varargin, 'segmentationstyle'); % this will be passed on to the corresponding ft_datatype_xxx function
parcellationstyle    = ft_getopt(varargin, 'parcellationstyle'); % this will be passed on to the corresponding ft_datatype_xxx function
hasbrain             = ft_getopt(varargin, 'hasbrain');
trialinfostyle       = ft_getopt(varargin, 'trialinfostyle');

% check whether people are using deprecated stuff
depHastrialdef = ft_getopt(varargin, 'hastrialdef');
if (~isempty(depHastrialdef))
  ft_warning('ft_checkdata option ''hastrialdef'' is deprecated; please use ''hassampleinfo'' instead');
  hassampleinfo = depHastrialdef;
end

% determine the type of input data
israw           = ft_datatype(data, 'raw');
isfreq          = ft_datatype(data, 'freq');
istimelock      = ft_datatype(data, 'timelock');
iscomp          = ft_datatype(data, 'comp');
isspike         = ft_datatype(data, 'spike');
isvolume        = ft_datatype(data, 'volume');
issegmentation  = ft_datatype(data, 'segmentation');
isparcellation  = ft_datatype(data, 'parcellation');
issource        = ft_datatype(data, 'source');
isdip           = ft_datatype(data, 'dip');
ismvar          = ft_datatype(data, 'mvar');
isfreqmvar      = ft_datatype(data, 'freqmvar');
ischan          = ft_datatype(data, 'chan');
ismesh          = ft_datatype(data, 'mesh');
% FIXME use the istrue function on ismeg and hasxxx options

if ~isequal(feedback, 'no') % can be 'yes' or 'text'
  if iscomp
    % it can be comp and raw/timelock/freq at the same time, therefore this has to go first
    nchan = size(data.topo,1);
    ncomp = size(data.topo,2);
    ft_info('the input is component data with %d components and %d original channels\n', ncomp, nchan);
  end  % if iscomp
  
  if ismesh
    % it can be comp and source at the same time, therefore this has to go first
    data = fixpos(data);
    npos = 0;
    ntri = 0;
    nhex = 0;
    ntet = 0;
    % the data can contain multiple surfaces
    for i=1:numel(data)
      npos = npos+size(data.pos,1);
      if isfield(data, 'tri'), ntri = ntri+size(data.tri,1); end
      if isfield(data, 'hex'), nhex = nhex+size(data.hex,1); end
      if isfield(data, 'tet'), ntet = ntet+size(data.tet,1); end
    end
    if isfield(data,'tri')
      ft_info('the input is mesh data with %d vertices and %d triangles\n', npos, ntri);
    elseif isfield(data,'hex')
      ft_info('the input is mesh data with %d vertices and %d hexahedrons\n', npos, nhex);
    elseif isfield(data,'tet')
      ft_info('the input is mesh data with %d vertices and %d tetrahedrons\n', npos, ntet);
    else
      ft_info('the input is mesh data with %d vertices', npos);
    end
  end % if ismesh
  
  if israw
    nchan = length(data.label);
    ntrial = length(data.trial);
    ft_info('the input is raw data with %d channels and %d trials\n', nchan, ntrial);
  elseif istimelock
    nchan = length(data.label);
    ntime = length(data.time);
    ft_info('the input is timelock data with %d channels and %d timebins\n', nchan, ntime);
  elseif isfreq
    if isfield(data, 'label')
      nchan = length(data.label);
      nfreq = length(data.freq);
      if isfield(data, 'time'), ntime = num2str(length(data.time)); else ntime = 'no'; end
      ft_info('the input is freq data with %d channels, %d frequencybins and %s timebins\n', nchan, nfreq, ntime);
    elseif isfield(data, 'labelcmb')
      nchan = length(data.labelcmb);
      nfreq = length(data.freq);
      if isfield(data, 'time'), ntime = num2str(length(data.time)); else ntime = 'no'; end
      ft_info('the input is freq data with %d channel combinations, %d frequencybins and %s timebins\n', nchan, nfreq, ntime);
    else
      ft_error('cannot infer freq dimensions');
    end
  elseif isspike
    nchan  = length(data.label);
    ft_info('the input is spike data with %d channels\n', nchan);
  elseif isvolume
    if issegmentation
      subtype = 'segmented volume';
    else
      subtype = 'volume';
    end
    ft_info('the input is %s data with dimensions [%d %d %d]\n', subtype, data.dim(1), data.dim(2), data.dim(3));
    clear subtype
  elseif issource
    data = fixpos(data); % ensure that positions are in pos, not in pnt
    nsource = size(data.pos, 1);
    if isparcellation
      subtype = 'parcellated source';
    else
      subtype = 'source';
    end
    if isfield(data, 'dim')
      ft_info('the input is %s data with %d brainordinates on a [%d %d %d] grid\n', subtype, nsource, data.dim(1), data.dim(2), data.dim(3));
    else
      ft_info('the input is %s data with %d brainordinates\n', subtype, nsource);
    end
    clear subtype
  elseif isdip
    ft_info('the input is dipole data\n');
  elseif ismvar
    ft_info('the input is mvar data\n');
  elseif isfreqmvar
    ft_info('the input is freqmvar data\n');
  elseif ischan
    nchan = length(data.label);
    if isfield(data, 'brainordinate')
      ft_info('the input is parcellated data with %d parcels\n', nchan);
    else
      ft_info('the input is chan data with %d channels\n', nchan);
    end
  end % if israw etc.
end % give feedback

if issource && isvolume
  % it should be either one or the other: the choice here is to represent it as volume description since that is simpler to handle
  % the conversion is done by removing the grid positions
  data = rmfield(data, 'pos');
  issource = false;
end

if isfield(data, 'trialinfo')
  if strcmp(trialinfostyle, 'table')
    if ismatrix(data.trialinfo)
      data.trialinfo = array2table(data.trialinfo);
    end
  elseif strcmp(trialinfostyle, 'matrix')
    if istable(data.trialinfo)
      data.trialinfo = table2array(data.trialinfo);
    end
  else
    % no conversion is needed
  end
end

% the ft_datatype_XXX functions ensures the consistency of the XXX datatype
% and provides a detailed description of the dataformat and its history
if iscomp % this should go before israw/istimelock/isfreq
  data = ft_datatype_comp(data, 'hassampleinfo', hassampleinfo);
elseif israw
  data = ft_datatype_raw(data, 'hassampleinfo', hassampleinfo);
elseif istimelock
  data = ft_datatype_timelock(data, 'hassampleinfo', hassampleinfo);
elseif isfreq
  data = ft_datatype_freq(data);
elseif isspike
  data = ft_datatype_spike(data);
elseif issegmentation % this should go before isvolume
  data = ft_datatype_segmentation(data, 'segmentationstyle', segmentationstyle, 'hasbrain', hasbrain);
elseif isvolume
  data = ft_datatype_volume(data);
elseif isparcellation % this should go before issource
  data = ft_datatype_parcellation(data, 'parcellationstyle', parcellationstyle);
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
      case 'raw+comp'
        okflag = okflag + (israw & iscomp);
      case 'freq+comp'
        okflag = okflag + (isfreq & iscomp);
      case 'timelock+comp'
        okflag = okflag + (istimelock & iscomp);
      case 'source+mesh'
        okflag = okflag + (issource & ismesh);
      case 'raw'
        okflag = okflag + (israw & ~iscomp);
      case 'freq'
        okflag = okflag + (isfreq & ~iscomp);
      case 'timelock'
        okflag = okflag + (istimelock & ~iscomp);
      case 'comp'
        okflag = okflag + (iscomp & ~(israw | istimelock | isfreq));
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
      case 'chan'
        okflag = okflag + ischan;
      case 'segmentation'
        okflag = okflag + issegmentation;
      case 'parcellation'
        okflag = okflag + isparcellation;
      case 'mesh'
        okflag = okflag + ismesh;
    end % switch dtype
  end % for dtype
  
  % try to convert the data if needed
  for iCell = 1:length(dtype)
    if okflag
      % the requested datatype is specified in descending order of
      % preference (if there is a preference at all), so don't bother
      % checking the rest of the requested data types if we already
      % succeeded in converting
      break;
    end
    if isequal(dtype(iCell), {'parcellation'}) && issegmentation
      data = volume2source(data); % segmentation=volume, parcellation=source
      data = ft_datatype_parcellation(data);
      issegmentation = 0;
      isvolume = 0;
      isparcellation = 1;
      issource = 1;
      okflag = 1;
    elseif isequal(dtype(iCell), {'segmentation'}) && isparcellation
      data = source2volume(data); % segmentation=volume, parcellation=source
      data = ft_datatype_segmentation(data);
      isparcellation = 0;
      issource = 0;
      issegmentation = 1;
      isvolume = 1;
      okflag = 1;
    elseif isequal(dtype(iCell), {'source'}) && isvolume
      data = volume2source(data);
      data = ft_datatype_source(data);
      isvolume = 0;
      issource = 1;
      okflag = 1;
    elseif isequal(dtype(iCell), {'volume'}) && (ischan || istimelock || isfreq)
      if isfield(data, 'brainordinate')
        data = parcellated2source(data);
        data = ft_datatype_volume(data);
      else
        ft_error('cannot convert channel-level data to volumetric representation');
      end
      ischan = 0; istimelock = 0; isfreq = 0;
      isvolume = 1;
      okflag = 1;
    elseif isequal(dtype(iCell), {'source'}) && (ischan || istimelock || isfreq)
      if isfield(data, 'brainordinate')
        data = parcellated2source(data);
        data = ft_datatype_source(data);
      else
        data = chan2source(data);
        data = ft_datatype_source(data);
      end % converting channel data
      ischan = 0; istimelock = 0; isfreq = 0;
      issource = 1;
      okflag = 1;
    elseif isequal(dtype(iCell), {'volume'}) && issource
      data = source2volume(data);
      data = ft_datatype_volume(data);
      isvolume = 1;
      issource = 0;
      okflag = 1;
    elseif isequal(dtype(iCell), {'raw'}) && issource
      data = source2raw(data);
      data = ft_datatype_raw(data, 'hassampleinfo', hassampleinfo);
      issource = 0;
      israw = 1;
      okflag = 1;
    elseif isequal(dtype(iCell), {'raw'}) && istimelock
      if iscomp
        data = removefields(data, {'topo', 'topolabel', 'topodimord', 'unmixing', 'unmixingdimord'}); % these fields are not desired
        iscomp = 0;
      end
      data = timelock2raw(data);
      data = ft_datatype_raw(data, 'hassampleinfo', hassampleinfo);
      istimelock = 0;
      israw = 1;
      okflag = 1;
    elseif isequal(dtype(iCell), {'raw+comp'}) && istimelock && iscomp
      data = timelock2raw(data);
      data = ft_datatype_raw(data, 'hassampleinfo', hassampleinfo);
      istimelock = 0;
      iscomp = 1;
      israw = 1;
      okflag = 1;
    elseif isequal(dtype(iCell), {'timelock+comp'}) && israw && iscomp
      data = raw2timelock(data);
      data = ft_datatype_timelock(data, 'hassampleinfo', hassampleinfo);
      istimelock = 1;
      iscomp = 1;
      israw = 0;
      okflag = 1;
    elseif isequal(dtype(iCell), {'comp'}) && israw  && iscomp
      data = keepfields(data, {'label', 'topo', 'topolabel', 'unmixing', 'elec', 'grad', 'cfg'}); % these are the only relevant fields
      data = ft_datatype_comp(data, 'hassampleinfo', hassampleinfo);
      israw = 0;
      iscomp = 1;
      okflag = 1;
    elseif isequal(dtype(iCell), {'comp'}) && istimelock && iscomp
      data = keepfields(data, {'label', 'topo', 'topolabel', 'unmixing', 'elec', 'grad', 'cfg'}); % these are the only relevant fields
      data = ft_datatype_comp(data, 'hassampleinfo', hassampleinfo);
      istimelock = 0;
      iscomp = 1;
      okflag = 1;
    elseif isequal(dtype(iCell), {'comp'}) && isfreq && iscomp
      data = keepfields(data, {'label', 'topo', 'topolabel', 'unmixing', 'elec', 'grad', 'cfg'}); % these are the only relevant fields
      data = ft_datatype_comp(data, 'hassampleinfo', 'no'); % freq data does not have sampleinfo
      isfreq = 0;
      iscomp = 1;
      okflag = 1;
    elseif isequal(dtype(iCell), {'raw'}) && israw
      if iscomp
        data = removefields(data, {'topo', 'topolabel', 'topodimord', 'unmixing', 'unmixingdimord'}); % these fields are not desired
        iscomp = 0;
      end
      data = ft_datatype_raw(data, 'hassampleinfo', hassampleinfo);
      okflag = 1;
    elseif isequal(dtype(iCell), {'timelock'}) && istimelock
      if iscomp
        data = removefields(data, {'topo', 'topolabel', 'topodimord', 'unmixing', 'unmixingdimord'}); % these fields are not desired
        iscomp = 0;
      end
      data = ft_datatype_timelock(data, 'hassampleinfo', hassampleinfo);
      okflag = 1;
    elseif isequal(dtype(iCell), {'freq'}) && isfreq
      if iscomp
        data = removefields(data, {'topo', 'topolabel', 'topodimord', 'unmixing', 'unmixingdimord'}); % these fields are not desired
        iscomp = 0;
      end
      data = ft_datatype_freq(data);
      okflag = 1;
    elseif isequal(dtype(iCell), {'timelock'}) && israw
      if iscomp
        data = removefields(data, {'topo', 'topolabel', 'topodimord', 'unmixing', 'unmixingdimord'}); % these fields are not desired
        iscomp = 0;
      end
      data = raw2timelock(data);
      data = ft_datatype_timelock(data, 'hassampleinfo', hassampleinfo);
      israw = 0;
      istimelock = 1;
      okflag = 1;
    elseif isequal(dtype(iCell), {'raw'}) && isfreq
      if iscomp
        data = removefields(data, {'topo', 'topolabel', 'topodimord', 'unmixing', 'unmixingdimord'}); % these fields are not desired
        iscomp = 0;
      end
      data = freq2raw(data);
      data = ft_datatype_raw(data, 'hassampleinfo', hassampleinfo);
      isfreq = 0;
      israw = 1;
      okflag = 1;
      
    elseif isequal(dtype(iCell), {'raw'}) && ischan
      data = chan2timelock(data);
      data = timelock2raw(data);
      data = ft_datatype_raw(data, 'hassampleinfo', hassampleinfo);
      ischan = 0;
      israw = 1;
      okflag = 1;
    elseif isequal(dtype(iCell), {'timelock'}) && ischan
      data = chan2timelock(data);
      data = ft_datatype_timelock(data, 'hassampleinfo', hassampleinfo);
      ischan = 0;
      istimelock = 1;
      okflag = 1;
    elseif isequal(dtype(iCell), {'freq'}) && ischan
      data = chan2freq(data);
      data = ft_datatype_freq(data);
      ischan = 0;
      isfreq = 1;
      okflag = 1;
    elseif isequal(dtype(iCell), {'spike'}) && israw
      data = raw2spike(data);
      data = ft_datatype_spike(data);
      israw = 0;
      isspike = 1;
      okflag = 1;
    elseif isequal(dtype(iCell), {'raw'}) && isspike
      data = spike2raw(data,fsample);
      data = ft_datatype_raw(data, 'hassampleinfo', hassampleinfo);
      isspike = 0;
      israw   = 1;
      okflag  = 1;
    end
  end % for iCell
  
  if ~okflag
    % construct an error message
    typestr = printor(dtype, true);
    helpfun = cell(size(dtype));
    for i=1:numel(dtype)
      helpfun{i} = sprintf('ft_datatype_%s', dtype{i});
    end
    helpfun = helpfun(cellfun(@exist, helpfun)>0);
    if ~isempty(helpfun)
      helpstr = printor(helpfun);
      ft_error('This function requires %s data as input, see %s.', typestr, helpstr);
    else
      ft_error('This function requires %s data as input.', typestr);
    end
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
    ft_error('This function requires data with a dimord of %s.', printor(dimord, true));
  end % if okflag
end

if ~isempty(stype)
  if ~isa(stype, 'cell')
    stype = {stype};
  end
  
  if isfield(data, 'grad') || isfield(data, 'elec') || isfield(data, 'opto')
    if any(strcmp(ft_senstype(data), stype))
      okflag = 1;
    elseif any(cellfun(@ft_senstype, repmat({data}, size(stype)), stype))
      % this is required to detect more general types, such as "meg" or "ctf" rather than "ctf275"
      okflag = 1;
    else
      okflag = 0;
    end
  else
    % the data does not contain a sensor array
    okflag = 0;
  end
  
  if ~okflag
    % construct an error message
    ft_error('This function requires data with an %s sensor array.', printor(stype, true));
  end % if okflag
end

if ~isempty(ismeg)
  if isequal(ismeg, 'yes')
    okflag = isfield(data, 'grad');
  elseif isequal(ismeg, 'no')
    okflag = ~isfield(data, 'grad');
  end
  
  if ~okflag && isequal(ismeg, 'yes')
    ft_error('This function requires MEG data with a ''grad'' field');
  elseif ~okflag && isequal(ismeg, 'no')
    ft_error('This function should not be given MEG data with a ''grad'' field');
  end % if okflag
end

if ~isempty(iseeg)
  if isequal(iseeg, 'yes')
    okflag = isfield(data, 'elec');
  elseif isequal(iseeg, 'no')
    okflag = ~isfield(data, 'elec');
  end
  
  if ~okflag && isequal(iseeg, 'yes')
    ft_error('This function requires EEG data with an ''elec'' field');
  elseif ~okflag && isequal(iseeg, 'no')
    ft_error('This function should not be given EEG data with an ''elec'' field');
  end % if okflag
end

if ~isempty(isnirs)
  if isequal(isnirs, 'yes')
    okflag = isfield(data, 'opto');
  elseif isequal(isnirs, 'no')
    okflag = ~isfield(data, 'opto');
  end
  
  if ~okflag && isequal(isnirs, 'yes')
    ft_error('This function requires NIRS data with an ''opto'' field');
  elseif ~okflag && isequal(isnirs, 'no')
    ft_error('This function should not be given NIRS data with an ''opto'' field');
  end % if okflag
end

if ~isempty(inside)
  if strcmp(inside, 'index')
    ft_warning('the indexed representation of inside/outside source locations is deprecated');
  end
  % TODO absorb the fixinside function into this code
  data   = fixinside(data, inside);
  okflag = isfield(data, 'inside');
  
  if ~okflag
    % construct an error message
    ft_error('This function requires data with an ''inside'' field.');
  end % if okflag
end

if istrue(hasunit) && ~isfield(data, 'unit')
  % calling convert_units with only the input data adds the units without converting
  data = ft_determine_units(data);
end % if hasunit

if istrue(hascoordsys) && ~isfield(data, 'coordsys')
  data = ft_determine_coordsys(data);
end % if hascoordsys

if istrue(haschantype) && ~isfield(data, 'chantype')
  data.chantype = ft_chantype(data);
end % if haschantype

if istrue(haschanunit) && ~isfield(data, 'chanunit')
  data.chanunit = ft_chanunit(data);
end % if haschanunit

if isequal(hastrials, 'yes')
  hasrpt = isfield(data, 'trial');
  if ~hasrpt && isfield(data, 'dimord')
    % instead look in the dimord for rpt or subj
    hasrpt = ~isempty(strfind(data.dimord, 'rpt')) || ...
      ~isempty(strfind(data.dimord, 'rpttap')) || ...
      ~isempty(strfind(data.dimord, 'subj'));
  end
  if ~hasrpt
    ft_error('This function requires data with a ''trial'' field');
  end % if hasrpt
elseif isequal(hastrials, 'no') && istimelock
  if ~isfield(data, 'avg') && (isfield(data, 'trial') || isfield(data, 'individual'))
    % average on the fly
    tmpcfg = [];
    tmpcfg.keeptrials = 'no';
    data = ft_timelockanalysis(tmpcfg, data);
  end
end

if strcmp(hasdim, 'yes') && ~isfield(data, 'dim')
  data.dim = pos2dim(data.pos);
elseif strcmp(hasdim, 'no') && isfield(data, 'dim')
  data = rmfield(data, 'dim');
end % if hasdim

if strcmp(hascumtapcnt, 'yes') && ~isfield(data, 'cumtapcnt')
  ft_error('This function requires data with a ''cumtapcnt'' field');
elseif strcmp(hascumtapcnt, 'no') && isfield(data, 'cumtapcnt')
  data = rmfield(data, 'cumtapcnt');
end % if hascumtapcnt

if strcmp(hasdof, 'yes') && ~isfield(data, 'dof')
  ft_error('This function requires data with a ''dof'' field');
elseif strcmp(hasdof, 'no') && isfield(data, 'dof')
  data = rmfield(data, 'dof');
end % if hasdof

if ~isempty(cmbrepresentation)
  if istimelock
    data = fixcov(data, cmbrepresentation);
  elseif isfreq
    data = fixcsd(data, cmbrepresentation, channelcmb);
  elseif isfreqmvar
    data = fixcsd(data, cmbrepresentation, channelcmb);
  else
    ft_error('This function requires data with a covariance, coherence or cross-spectrum');
  end
end % cmbrepresentation

if isfield(data, 'grad')
  % ensure that the gradiometer structure is up to date
  data.grad = ft_datatype_sens(data.grad);
end

if isfield(data, 'elec')
  % ensure that the electrode structure is up to date
  data.elec = ft_datatype_sens(data.elec);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% represent the covariance matrix in a particular manner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = fixcov(data, desired)
if any(isfield(data, {'cov', 'corr'}))
  if ~isfield(data, 'labelcmb')
    current = 'full';
  else
    current = 'sparse';
  end
else
  ft_error('Could not determine the current representation of the covariance matrix');
end
if isequal(current, desired)
  % nothing to do
elseif strcmp(current, 'full') && strcmp(desired, 'sparse')
  % FIXME should be implemented
  ft_error('not yet implemented');
elseif strcmp(current, 'sparse') && strcmp(desired, 'full')
  % FIXME should be implemented
  ft_error('not yet implemented');
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
  ft_error('Could not determine the current representation of the %s matrix', param);
end

% first go from univariate fourier to the required bivariate representation
if isequal(current, desired)
  % nothing to do
  
elseif strcmp(current, 'fourier') && strcmp(desired, 'sparsewithpow')
  dimtok = tokenize(data.dimord, '_');
  if ~isempty(strmatch('rpttap',   dimtok))
    nrpt = size(data.cumtapcnt,1);
    flag = 0;
  else
    nrpt = 1;
  end
  if ~isempty(strmatch('freq',  dimtok)), nfrq=length(data.freq);      else nfrq = 1; end
  if ~isempty(strmatch('time',  dimtok)), ntim=length(data.time);      else ntim = 1; end
  
  fastflag = all(data.cumtapcnt(:)==data.cumtapcnt(1));
  flag     = nrpt==1; % needed to truncate the singleton dimension upfront
  
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
    powspctrm = zeros(nrpt,nchan,nfrq,ntim) + zeros(nrpt,nchan,nfrq,ntim)*1i;
    sumtapcnt = [0;cumsum(data.cumtapcnt(:))];
    for p = 1:nrpt
      indx   = (sumtapcnt(p)+1):sumtapcnt(p+1);
      tmpdat = data.fourierspctrm(indx,:,:,:);
      powspctrm(p,:,:,:) = (sum(tmpdat.*conj(tmpdat),1))./data.cumtapcnt(p);
    end
  end
  
  %create cross-spectra
  if ~isempty(channelcmb)
    ncmb      = size(channelcmb,1);
    cmbindx   = zeros(ncmb,2);
    labelcmb  = cell(ncmb,2);
    for k = 1:ncmb
      ch1 = find(strcmp(data.label, channelcmb(k,1)));
      ch2 = find(strcmp(data.label, channelcmb(k,2)));
      if ~isempty(ch1) && ~isempty(ch2)
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
  if ntim>1
    data.dimord = 'chan_freq_time';
  else
    data.dimord = 'chan_freq';
  end
  
  if nrpt>1
    data.dimord = ['rpt_',data.dimord];
  end
  
  if flag
    siz = size(data.powspctrm);
    data.powspctrm = reshape(data.powspctrm, [siz(2:end) 1]);
    if isfield(data, 'crsspctrm')
      siz = size(data.crsspctrm);
      data.crsspctrm = reshape(data.crsspctrm, [siz(2:end) 1]);
    end
  end
elseif strcmp(current, 'fourier') && strcmp(desired, 'sparse')
  
  if isempty(channelcmb), ft_error('no channel combinations are specified'); end
  dimtok = tokenize(data.dimord, '_');
  if ~isempty(strmatch('rpttap',   dimtok))
    nrpt = size(data.cumtapcnt,1);
    flag = 0;
  else
    nrpt = 1;
  end
  if ~isempty(strmatch('freq',  dimtok)), nfrq=length(data.freq); else nfrq = 1; end
  if ~isempty(strmatch('time',  dimtok)), ntim=length(data.time); else ntim = 1; end
  
  flag      = nrpt==1; % flag needed to squeeze first dimension if singleton
  ncmb      = size(channelcmb,1);
  cmbindx   = zeros(ncmb,2);
  labelcmb  = cell(ncmb,2);
  for k = 1:ncmb
    ch1 = find(strcmp(data.label, channelcmb(k,1)));
    ch2 = find(strcmp(data.label, channelcmb(k,2)));
    if ~isempty(ch1) && ~isempty(ch2)
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
  if ntim>1
    data.dimord = 'chancmb_freq_time';
  else
    data.dimord = 'chancmb_freq';
  end
  
  if nrpt>1
    data.dimord = ['rpt_',data.dimord];
  end
  
  if flag
    if isfield(data,'powspctrm')
      % deal with the singleton 'rpt', i.e. remove it
      siz = size(data.powspctrm);
      data.powspctrm = reshape(data.powspctrm, [siz(2:end) 1]);
    end
    if isfield(data,'crsspctrm')
      % this conditional statement is needed in case there's a single channel
      siz            = size(data.crsspctrm);
      data.crsspctrm = reshape(data.crsspctrm, [siz(2:end) 1]);
    end
  end
elseif strcmp(current, 'fourier') && strcmp(desired, 'full')
  
  % this is how it is currently and the desired functionality of prepare_freq_matrices
  dimtok = tokenize(data.dimord, '_');
  if ~isempty(strmatch('rpttap',   dimtok))
    nrpt = size(data.cumtapcnt, 1);
    flag = 0;
  else
    nrpt = 1;
    flag = 1;
  end
  if ~isempty(strmatch('rpttap',dimtok)), nrpt=size(data.cumtapcnt, 1); else nrpt = 1; end
  if ~isempty(strmatch('freq',  dimtok)), nfrq=length(data.freq);       else nfrq = 1; end
  if ~isempty(strmatch('time',  dimtok)), ntim=length(data.time);       else ntim = 1; end
  if any(data.cumtapcnt(1,:) ~= data.cumtapcnt(1,1)), ft_error('this only works when all frequencies have the same number of tapers'); end
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
  
  if isfield(data, 'trialinfo'),  data = rmfield(data, 'trialinfo'); end
  if isfield(data, 'sampleinfo'), data = rmfield(data, 'sampleinfo'); end
  if isfield(data, 'cumsumcnt'),  data = rmfield(data, 'cumsumcnt');  end
  if isfield(data, 'cumtapcnt'),  data = rmfield(data, 'cumtapcnt');  end
  
end % convert to the requested bivariate representation

% from one bivariate representation to another
if isequal(current, desired)
  % nothing to do
  
elseif (strcmp(current, 'full')       && strcmp(desired, 'fourier')) || ...
    (strcmp(current, 'sparse')        && strcmp(desired, 'fourier')) || ...
    (strcmp(current, 'sparsewithpow') && strcmp(desired, 'fourier'))
  % this is not possible
  ft_error('converting the cross-spectrum into a Fourier representation is not possible');
  
elseif strcmp(current, 'full') && strcmp(desired, 'sparsewithpow')
  ft_error('not yet implemented');
  
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
  if ~isempty(strmatch('rpt',   dimtok)), nrpt=size(data.cumtapcnt,1); else nrpt = 1; end
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
    if numel(data.(fn{ii})) == nrpt*ncmb*nfrq*ntim
      if nrpt>1
        data.(fn{ii}) = reshape(data.(fn{ii}), nrpt, ncmb, nfrq, ntim);
      else
        data.(fn{ii}) = reshape(data.(fn{ii}), ncmb, nfrq, ntim);
      end
    end
  end
  % remove obsolete fields
  data           = rmfield(data, 'label');
  try data      = rmfield(data, 'dof'); end
  % replace updated fields
  data.labelcmb  = labelcmb;
  if ntim>1
    data.dimord = 'chancmb_freq_time';
  else
    data.dimord = 'chancmb_freq';
  end
  
  if nrpt>1
    data.dimord = ['rpt_',data.dimord];
  end
  
elseif strcmp(current, 'sparsewithpow') && strcmp(desired, 'sparse')
  % this representation for sparse data contains autospectra as e.g. {'A' 'A'} in labelcmb
  if isfield(data, 'crsspctrm')
    dimtok         = tokenize(data.dimord, '_');
    catdim         = match_str(dimtok, {'chan' 'chancmb'});
    data.crsspctrm = cat(catdim, data.powspctrm, data.crsspctrm);
    data.labelcmb  = [data.label(:) data.label(:); data.labelcmb];
    data           = rmfield(data, 'powspctrm');
    data.dimord    = strrep(data.dimord, 'chan_', 'chancmb_');
  else
    data.crsspctrm = data.powspctrm;
    data.labelcmb  = [data.label(:) data.label(:)];
    data           = rmfield(data, 'powspctrm');
    data.dimord    = strrep(data.dimord, 'chan_', 'chancmb_');
  end
  data = rmfield(data, 'label');
  
elseif strcmp(current, 'sparse') && strcmp(desired, 'full')
  dimtok = tokenize(data.dimord, '_');
  if ~isempty(strmatch('rpt',   dimtok)), nrpt=size(data.cumtapcnt,1); else nrpt = 1; end
  if ~isempty(strmatch('freq',  dimtok)), nfrq=numel(data.freq);      else nfrq = 1; end
  if ~isempty(strmatch('time',  dimtok)), ntim=numel(data.time);      else ntim = 1; end
  
  if ~isfield(data, 'label')
    % ensure that the bivariate spectral factorization results can be
    % processed. FIXME this is experimental and will not work if the user
    % did something weird before
    for k = 1:numel(data.labelcmb)
      tmp = tokenize(data.labelcmb{k}, '[');
      data.labelcmb{k} = tmp{1};
    end
    data.label = unique(data.labelcmb(:));
  end
  
  nchan     = length(data.label);
  ncmb      = size(data.labelcmb,1);
  cmbindx   = zeros(nchan,nchan);
  
  for k = 1:size(data.labelcmb,1)
    ch1 = find(strcmp(data.label, data.labelcmb(k,1)));
    ch2 = find(strcmp(data.label, data.labelcmb(k,2)));
    if ~isempty(ch1) && ~isempty(ch2)
      cmbindx(ch1,ch2) = k;
    end
  end
  
  complete = all(cmbindx(:)~=0);
  
  % remove obsolete fields
  try data      = rmfield(data, 'powspctrm');  end
  try data      = rmfield(data, 'labelcmb');   end
  try data      = rmfield(data, 'dof');        end
  
  fn = fieldnames(data);
  for ii=1:numel(fn)
    if numel(data.(fn{ii})) == nrpt*ncmb*nfrq*ntim
      if nrpt==1
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
      if nrpt>1
        data.(fn{ii}) = tmpall;
      else
        data.(fn{ii}) = reshape(tmpall, [nchan nchan nfrq ntim]);
      end
    end % if numel
  end % for ii
  
  if ntim>1
    data.dimord = 'chan_chan_freq_time';
  else
    data.dimord = 'chan_chan_freq';
  end
  
  if nrpt>1
    data.dimord = ['rpt_',data.dimord];
  end
  
elseif strcmp(current, 'sparse') && strcmp(desired, 'fullfast')
  dimtok = tokenize(data.dimord, '_');
  if ~isempty(strmatch('rpt',   dimtok)), nrpt=size(data.cumtapcnt,1); else nrpt = 1; end
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
    if ~isempty(ch1) && ~isempty(ch2)
      cmbindx(ch1,ch2) = k;
    end
  end
  
  complete = all(cmbindx(:)~=0);
  
  fn = fieldnames(data);
  for ii=1:numel(fn)
    if numel(data.(fn{ii})) == nrpt*ncmb*nfrq*ntim
      if nrpt==1
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
      if nrpt>1
        data.(fn{ii}) = tmpall;
      else
        data.(fn{ii}) = reshape(tmpall, [nchan nchan nfrq ntim]);
      end
    end % if numel
  end % for ii
  
  % remove obsolete fields
  try data      = rmfield(data, 'powspctrm');  end
  try data      = rmfield(data, 'labelcmb');   end
  try data      = rmfield(data, 'dof');        end
  
  if ntim>1
    data.dimord = 'chan_chan_freq_time';
  else
    data.dimord = 'chan_chan_freq';
  end
  
elseif strcmp(current, 'sparsewithpow') && any(strcmp(desired, {'full', 'fullfast'}))
  % recursively call ft_checkdata, but ensure channel order to be the same as the original input.
  origlabelorder = data.label; % keep track of the original order of the channels
  data       = ft_checkdata(data, 'cmbrepresentation', 'sparse');
  data.label = origlabelorder; % this avoids the labels to be alphabetized in the next call
  data       = ft_checkdata(data, 'cmbrepresentation', 'full');
  
end % convert from one to another bivariate representation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [source] = chan2source(data)
chanpos = zeros(0,3);
chanlab = cell(0,1);
posunit = [];
if isfield(data, 'elec')
  chanpos = cat(1, chanpos, data.elec.chanpos);
  chanlab = cat(1, chanlab, data.elec.label);
  if isfield(data.elec, 'unit')
    posunit = data.elec.unit;
  end
end
if isfield(data, 'grad')
  chanpos = cat(1, chanpos, data.grad.chanpos);
  chanlab = cat(1, chanlab, data.grad.label);
  if isfield(data.grad, 'unit')
    posunit = data.grad.unit;
  end
end
if isfield(data, 'opto')
  chanpos = cat(1, chanpos, data.opto.chanpos);
  chanlab = cat(1, chanlab, data.opto.label);
  if isfield(data.opto, 'unit')
    posunit = data.opto.unit;
  end
end

fn = fieldnames(data);
fn = setdiff(fn, {'label', 'time', 'freq', 'hdr', 'cfg', 'grad', 'elec', 'dimord', 'unit'}); % remove irrelevant fields
fn(~cellfun(@isempty, regexp(fn, 'dimord$'))) = []; % remove irrelevant (dimord) fields
sel = false(size(fn));
for i=1:numel(fn)
  try
    sel(i) = ismember(getdimord(data, fn{i}), {'chan', 'chan_time', 'chan_freq', 'chan_freq_time', 'chan_chan'});
  end
end
parameter = fn(sel);

% determine the channel indices for which the position is known
[datsel, possel] = match_str(data.label, chanlab);

source = [];
source.pos = chanpos(possel, :);
if ~isempty(posunit)
  source.unit = posunit;
end
for i=1:numel(parameter)
  dat = data.(parameter{i});
  dimord = getdimord(data, parameter{i});
  dimtok = tokenize(dimord, '_');
  for dim=1:numel(dimtok)
    if strcmp(dimtok{dim}, 'chan')
      dat = dimindex(dat, dim, {datsel});
      dimtok{dim} = 'pos';
    end
  end
  dimord = sprintf('%s_', dimtok{:});
  dimord = dimord(1:end-1); % remove the last '_'
  % copy the data to the source representation
  source.(parameter{i})            = dat;
  source.([parameter{i} 'dimord']) = dimord;
end
% copy the descriptive fields, these are necessary for visualising the data in ft_sourceplot
source = copyfields(data, source, {'time', 'freq'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [source] = parcellated2source(data)
if ~isfield(data, 'brainordinate')
  ft_error('projecting parcellated data onto the full brain model geometry requires the specification of brainordinates');
end
% the main structure contains the functional data on the parcels
% the brainordinate sub-structure contains the original geometrical model
source = ft_checkdata(data.brainordinate, 'datatype', 'source');
data   = rmfield(data, 'brainordinate');
if isfield(data, 'cfg')
  source.cfg = data.cfg;
end

fn = fieldnames(data);
fn = setdiff(fn, {'label', 'time', 'freq', 'hdr', 'cfg', 'grad', 'elec', 'dimord', 'unit'}); % remove irrelevant fields
fn(~cellfun(@isempty, regexp(fn, 'dimord$'))) = []; % remove irrelevant (dimord) fields
sel = false(size(fn));
for i=1:numel(fn)
  try
    sel(i) = ismember(getdimord(data, fn{i}), {'chan', 'chan_time', 'chan_freq', 'chan_freq_time', 'chan_chan'});
  end
end
parameter = fn(sel);

fn = fieldnames(source);
sel = false(size(fn));
for i=1:numel(fn)
  tmp = source.(fn{i});
  sel(i) = iscell(tmp) && isequal(tmp(:), data.label(:));
end
parcelparam = fn(sel);
if numel(parcelparam)~=1
  ft_error('cannot determine which parcellation to use');
else
  parcelparam = parcelparam{1}(1:(end-5)); % minus the 'label'
end

for i=1:numel(parameter)
  source.(parameter{i}) = unparcellate(data, source, parameter{i}, parcelparam);
end

% copy the descriptive fields, these are necessary for visualising the data in ft_sourceplot
source = copyfields(data, source, {'time', 'freq'});


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
  data.pos = ft_warp_apply(data.transform, [x(:) y(:) z(:)]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = source2volume(data)

fn = fieldnames(data);
fd = nan(size(fn));
for i=1:numel(fn)
  fd(i) = ndims(data.(fn{i}));
end

if ~isfield(data, 'dim')
  % this part depends on the assumption that the list of positions is describing a full 3D volume in
  % an ordered way which allows for the extraction of a transformation matrix, i.e. slice by slice
  data.dim = pos2dim(data.pos);
  try
    % if the dim is correct, it should be possible to obtain the transform
    ws = warning('off', 'MATLAB:rankDeficientMatrix');
    pos2transform(data.pos, data.dim);
    warning(ws);
  catch
    % remove the incorrect dim
    data = rmfield(data, 'dim');
  end
end

if isfield(data, 'dim')
  data.transform = pos2transform(data.pos, data.dim);
end

% remove the unwanted fields
data = removefields(data, {'pos', 'xgrid', 'ygrid', 'zgrid', 'tri', 'tet', 'hex'});

% make inside a volume
data = fixinside(data, 'logical');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = freq2raw(freq)

if isfield(freq, 'powspctrm')
  param = 'powspctrm';
elseif isfield(freq, 'fourierspctrm')
  param = 'fourierspctrm';
else
  ft_error('not supported for this data representation');
end

if strcmp(freq.dimord, 'rpt_chan_freq_time') || strcmp(freq.dimord, 'rpttap_chan_freq_time')
  dat = freq.(param);
elseif strcmp(freq.dimord, 'chan_freq_time')
  dat = freq.(param);
  dat = reshape(dat, [1 size(dat)]); % add a singleton dimension
else
  ft_error('not supported for dimord %s', freq.dimord);
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
  if any(sum(isnan(data.trial{i}),1)==size(data.trial{i},1))
    tmp = sum(~isfinite(data.trial{i}),1)==size(data.trial{i},1);
    begsmp = find(~tmp,1, 'first');
    endsmp = find(~tmp,1, 'last' );
    data.trial{i} = data.trial{i}(:, begsmp:endsmp);
    data.time{i}  = data.time{i}(begsmp:endsmp);
  end
end

if isfield(freq, 'trialinfo'), data.trialinfo = freq.trialinfo; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tlck] = raw2timelock(data)

data   = ft_checkdata(data, 'hassampleinfo', 'yes');
ntrial = numel(data.trial);
nchan  = numel(data.label);

if ntrial==1
  tlck.time   = data.time{1};
  tlck.avg    = data.trial{1};
  tlck.label  = data.label;
  tlck.dimord = 'chan_time';
  tlck        = copyfields(data, tlck, {'grad', 'elec', 'opto', 'cfg', 'trialinfo', 'topo', 'topodimord', 'topolabel', 'unmixing', 'unmixingdimord'});
  
else
  % the code below tries to construct a general time-axis where samples of all trials can fall on
  % find the earliest beginning and latest ending
  begtime = min(cellfun(@min, data.time));
  endtime = max(cellfun(@max, data.time));
  % find 'common' sampling rate
  fsample = 1./nanmean(cellfun(@mean, cellfun(@diff,data.time, 'uniformoutput', false)));
  % estimate number of samples
  nsmp = round((endtime-begtime)*fsample) + 1; % numerical round-off issues should be dealt with by this round, as they will/should never cause an extra sample to appear
  % construct general time-axis
  time = linspace(begtime,endtime,nsmp);
  
  % concatenate all trials
  tmptrial = nan(ntrial, nchan, length(time));
  
  begsmp = nan(ntrial, 1);
  endsmp = nan(ntrial, 1);
  for i=1:ntrial
    begsmp(i) = nearest(time, data.time{i}(1));
    endsmp(i) = nearest(time, data.time{i}(end));
    tmptrial(i,:,begsmp(i):endsmp(i)) = data.trial{i};
  end
  
  % construct the output timelocked data
  tlck.trial   = tmptrial;
  tlck.time    = time;
  tlck.dimord  = 'rpt_chan_time';
  tlck.label   = data.label;
  tlck         = copyfields(data, tlck, {'grad', 'elec', 'opto', 'cfg', 'trialinfo', 'sampleinfo', 'topo', 'topodimord', 'topolabel', 'unmixing', 'unmixingdimord'});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = timelock2raw(data)
fn = getdatfield(data);
if any(ismember(fn, {'trial', 'individual', 'avg'}))
  % trial, individual and avg (in that order) should be preferred over all other data fields
  % see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2965#c12
  fn = fn(ismember(fn, {'trial', 'individual', 'avg'}));
end
dimord = cell(size(fn));
for i=1:numel(fn)
  % determine the dimensions of each of the data fields
  dimord{i} = getdimord(data, fn{i});
end
% the fields trial, individual and avg (with their corresponding default dimord) are preferred
if sum(strcmp(dimord, 'rpt_chan_time'))==1
  fn = fn{strcmp(dimord, 'rpt_chan_time')};
  ft_info('constructing trials from "%s"\n', fn);
  dimsiz = getdimsiz(data, fn);
  ntrial = dimsiz(1);
  nchan  = dimsiz(2);
  ntime  = dimsiz(3);
  tmptrial = {};
  tmptime  = {};
  for j=1:ntrial
    tmptrial{j} = reshape(data.(fn)(j,:,:), [nchan, ntime]);
    tmptime{j}  = data.time;
  end
  data       = rmfield(data, fn);
  data.trial = tmptrial;
  data.time  = tmptime;
elseif sum(strcmp(dimord, 'subj_chan_time'))==1
  fn = fn{strcmp(dimord, 'subj_chan_time')};
  ft_info('constructing trials from "%s"\n', fn);
  dimsiz = getdimsiz(data, fn);
  nsubj = dimsiz(1);
  nchan  = dimsiz(2);
  ntime  = dimsiz(3);
  tmptrial = {};
  tmptime  = {};
  for j=1:nsubj
    tmptrial{j} = reshape(data.(fn)(j,:,:), [nchan, ntime]);
    tmptime{j}  = data.time;
  end
  data       = rmfield(data, fn);
  data.trial = tmptrial;
  data.time  = tmptime;
elseif sum(strcmp(dimord, 'chan_time'))==1
  fn = fn{strcmp(dimord, 'chan_time')};
  ft_info('constructing single trial from "%s"\n', fn);
  data.time  = {data.time};
  data.trial = {data.(fn)};
  data = rmfield(data, fn);
else
  ft_error('unsupported data structure');
end
% remove unwanted fields
data = removefields(data, {'avg', 'var', 'cov', 'dimord', 'numsamples' ,'dof'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = chan2freq(data)
data.dimord = [data.dimord '_freq'];
data.freq   = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = chan2timelock(data)
data.dimord = [data.dimord '_time'];
data.time   = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spike] = raw2spike(data)
ft_info('converting raw data into spike data\n');
nTrials 	 = length(data.trial);
[spikelabel] = detectspikechan(data);
spikesel     = match_str(data.label, spikelabel);
nUnits       = length(spikesel);
if nUnits==0
  ft_error('cannot convert raw data to spike format since the raw data structure does not contain spike channels');
end

trialTimes  = zeros(nTrials,2);
for iUnit = 1:nUnits
  unitIndx = spikesel(iUnit);
  spikeTimes  = []; % we dont know how large it will be, so use concatenation inside loop
  trialInds   = [];
  for iTrial = 1:nTrials
    
    % read in the spike times
    [spikeTimesTrial]    = getspiketimes(data, iTrial, unitIndx);
    nSpikes              = length(spikeTimesTrial);
    spikeTimes           = [spikeTimes; spikeTimesTrial(:)];
    trialInds            = [trialInds; ones(nSpikes,1)*iTrial];
    
    % get the begs and ends of trials
    hasNum = find(~isnan(data.time{iTrial}));
    if iUnit==1, trialTimes(iTrial,:) = data.time{iTrial}([hasNum(1) hasNum(end)]); end
  end
  
  spike.label{iUnit}     = data.label{unitIndx};
  spike.waveform{iUnit}  = [];
  spike.time{iUnit}      = spikeTimes(:)';
  spike.trial{iUnit}     = trialInds(:)';
  
  if iUnit==1, spike.trialtime             = trialTimes; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = spike2raw(spike, fsample)

if nargin<2 || isempty(fsample)
  timeDiff = abs(diff(sort([spike.time{:}])));
  fsample  = 1/min(timeDiff(timeDiff>0));
  ft_warning('Desired sampling rate for spike data not specified, automatically resampled to %f', fsample);
end

% get some sizes
nUnits  = length(spike.label);
nTrials = size(spike.trialtime,1);

% preallocate
data.trial(1:nTrials) = {[]};
data.time(1:nTrials)  = {[]};
for iTrial = 1:nTrials
  
  % make bins: note that the spike.time is already within spike.trialtime
  x = [spike.trialtime(iTrial,1):(1/fsample):spike.trialtime(iTrial,2)];
  timeBins   = [x x(end)+1/fsample] - (0.5/fsample);
  time       = (spike.trialtime(iTrial,1):(1/fsample):spike.trialtime(iTrial,2));
  
  % convert to continuous
  trialData = zeros(nUnits,length(time));
  for iUnit = 1:nUnits
    
    % get the timestamps and only select those timestamps that are in the trial
    ts       = spike.time{iUnit};
    hasTrial = spike.trial{iUnit}==iTrial;
    ts       = ts(hasTrial);
    
    N = histc(ts,timeBins);
    if isempty(N)
      N = zeros(1,length(timeBins)-1);
    else
      N(end) = [];
    end
    
    % store it in a matrix
    trialData(iUnit,:) = N;
  end
  
  data.trial{iTrial} = trialData;
  data.time{iTrial}  = time;
  
end % for all trials

% create the associated labels and other aspects of data such as the header
data.label = spike.label;
data.fsample = fsample;
if isfield(spike,'hdr'), data.hdr = spike.hdr; end
if isfield(spike,'cfg'), data.cfg = spike.cfg; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = source2raw(source)

fn = fieldnames(source);
fn = setdiff(fn, {'pos', 'dim', 'transform', 'time', 'freq', 'cfg'});
for i=1:length(fn)
  dimord{i} = getdimord(source, fn{i});
end
sel = strcmp(dimord, 'pos_time');
assert(sum(sel)>0, 'the source structure does not contain a suitable field to represent as raw channel-level data');
assert(sum(sel)<2, 'the source structure contains multiple fields that can be represented as raw channel-level data');
fn     = fn{sel};
dimord = dimord{sel};

switch dimord
  case 'pos_time'
    % add fake raw channel data to the original data structure
    data.trial{1} = source.(fn);
    data.time{1}  = source.time;
    % add fake channel labels
    data.label = {};
    for i=1:size(source.pos,1)
      data.label{i} = sprintf('source%d', i);
    end
    data.label = data.label(:);
    data.cfg = source.cfg;
  otherwise
    % FIXME other formats could be implemented as well
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for detection of channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spikelabel, eeglabel] = detectspikechan(data)

maxRate = 2000; % default on what we still consider a neuronal signal: this firing rate should never be exceeded

% autodetect the spike channels
ntrial = length(data.trial);
nchans  = length(data.label);
spikechan = zeros(nchans,1);
for i=1:ntrial
  for j=1:nchans
    hasAllInts    = all(isnan(data.trial{i}(j,:)) | data.trial{i}(j,:) == round(data.trial{i}(j,:)));
    hasAllPosInts = all(isnan(data.trial{i}(j,:)) | data.trial{i}(j,:)>=0);
    T = nansum(diff(data.time{i}),2); % total time
    fr            = nansum(data.trial{i}(j,:),2) ./ T;
    spikechan(j)  = spikechan(j) + double(hasAllInts & hasAllPosInts & fr<=maxRate);
  end
end
spikechan = (spikechan==ntrial);

spikelabel = data.label(spikechan);
eeglabel   = data.label(~spikechan);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spikeTimes] = getspiketimes(data, trial, unit)
spikeIndx       = logical(data.trial{trial}(unit,:));
spikeCount      = data.trial{trial}(unit,spikeIndx);
spikeTimes      = data.time{trial}(spikeIndx);
if isempty(spikeTimes), return; end
multiSpikes     = find(spikeCount>1);
% get the additional samples and spike times, we need only loop through the bins
[addSamples, addTimes]   = deal([]);
for iBin = multiSpikes(:)' % looping over row vector
  addTimes     = [addTimes ones(1,spikeCount(iBin))*spikeTimes(iBin)];
  addSamples   = [addSamples ones(1,spikeCount(iBin))*spikeIndx(iBin)];
end
% before adding these times, first remove the old ones
spikeTimes(multiSpikes) = [];
spikeTimes              = sort([spikeTimes(:); addTimes(:)]);
