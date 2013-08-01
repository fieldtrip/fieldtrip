function [sens] = ft_datatype_sens(sens, varargin)

% FT_DATATYPE_SENS describes the FieldTrip structure that represents
% an EEG, ECoG, or MEG sensor array. This structure is commonly called
% "elec" for EEG and "grad" for MEG, or more general "sens" for either
% one.
%
% The structure for MEG gradiometers and/or magnetometers contains
%    sens.label    = Mx1 cell-array with channel labels
%    sens.chanpos  = Mx3 matrix with channel positions
%    sens.chanori  = Mx3 matrix with channel orientations, used for synthetic planar gradient computation
%    sens.tra      = MxN matrix to combine coils into channels
%    sens.coilpos  = Nx3 matrix with coil positions
%    sens.coilori  = Nx3 matrix with coil orientations
%    sens.balance  = structure containing info about the balancing, See FT_APPLY_MONTAGE
%
% The structure for EEG or ECoG channels contains
%    sens.label    = Mx1 cell-array with channel labels
%    sens.chanpos  = Mx3 matrix with channel positions
%    sens.tra      = MxN matrix to combine electrodes into channels
%    sens.elecpos  = Nx3 matrix with electrode positions
% In case sens.tra is not present in the EEG sensor array, the channels
% are assumed to be average referenced.
%
% The following fields are optional
%    sens.type     = string with the MEG or EEG acquisition system, see FT_SENSTYPE
%    sens.chantype = Mx1 cell-array with the type of the channel, see FT_CHANTYPE
%    sens.chanunit = Mx1 cell-array with the units of the channel signal, e.g. 'T', 'fT' or 'fT/cm'
%    sens.fid      = structure with fiducial information
%
% Revision history:
%
% (upcoming) The default units for all sensor descriptions will be changed to
%  SI units, i.e. T, V and m. It will be possible to convert the amplitude and
%  distance units (e.g. from T to fT and from m to mm).
%  It will be possible to express planar and axial gradiometer channels in units
%  of amplitude or in units of amplitude/distance (i.e. gradient).
%
% (2011v2/latest) The chantype and chanunit have been added for MEG.
%
% (2011v1) To facilitate determining the position of channels (e.g. for plotting)
%  in case of balanced MEG or bipolar EEG, an explicit distinction has been made
%  between chanpos+chanori and coilpos+coilori (for MEG) and chanpos and elecpos
%  (for EEG). The pnt and ori fields are removed
%
% (2010) Added support for bipolar or otherwise more complex linear combinations
%  of EEG electrodes using sens.tra, similar to MEG.
%
% (2009) Noice reduction has been added for MEG systems in the balance field.
%
% (2006) The optional fields sens.type and sens.unit were added.
%
% (2003) The initial version was defined, which looked like this for EEG
%    sens.pnt     = Mx3 matrix with electrode positions
%    sens.label   = Mx1 cell-array with channel labels
% and like this for MEG
%    sens.pnt     = Nx3 matrix with coil positions
%    sens.ori     = Nx3 matrix with coil orientations
%    sens.tra     = MxN matrix to combine coils into channels
%    sens.label   = Mx1 cell-array with channel labels
%
% See also FT_READ_SENS, FT_SENSTYPE, FT_CHANTYPE, FT_APPLY_MONTAGE, CTF2GRAD, FIF2GRAD,
% BTI2GRAD, YOKOGAWA2GRAD, ITAB2GRAD

% Copyright (C) 2011-2013, Robert Oostenveld & Jan-Mathijs Schoffelen
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

% undocumented options for the upcoming (2013?) format
%   amplitude     = string, can be 'T' or 'fT'
%   distance      = string, can be 'm', 'cm' or 'mm'
%   scaling       = string, can be 'amplitude' or 'amplitude/distance'

% these are for remembering the type on subsequent calls with the same input arguments
persistent previous_argin previous_argout

current_argin = [{sens} varargin];
if isequal(current_argin, previous_argin)
  % don't do the whole cheking again, but return the previous output from cache
  sens = previous_argout{1};
  return
end

% get the optional input arguments, which should be specified as key-value pairs
version   = ft_getopt(varargin, 'version', 'latest');
amplitude = ft_getopt(varargin, 'amplitude');
distance  = ft_getopt(varargin, 'distance');
scaling   = ft_getopt(varargin, 'scaling');

if strcmp(version, 'latest')
  % NOTE TO SELF: once finalized, the code for 2011v2 and 2013 has to be merged into a single latest version
  version = '2011v2';
end

if isempty(sens)
  return;
end

switch version
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case '2011v2'
    if ~isempty(amplitude) || ~isempty(distance) || ~isempty(scaling)
      warning('amplitude, distance and scaling are not supported for version "%s"', version);
    end
    
    % This speeds up subsequent calls to ft_senstype and channelposition.
    % However, if it is not more precise than MEG or EEG, don't keep it in
    % the output (see further down).
    if ~isfield(sens, 'type')
      sens.type = ft_senstype(sens);
    end
    
    % there are many cases which deal with either eeg or meg
    ismeg = ft_senstype(sens, 'meg');
    
    if isfield(sens, 'pnt')
      if ismeg
        % sensor description is a MEG sensor-array, containing oriented coils
        sens.coilpos = sens.pnt; sens = rmfield(sens, 'pnt');
        sens.coilori = sens.ori; sens = rmfield(sens, 'ori');
      else
        % sensor description is something else, EEG/ECoG etc
        sens.elecpos = sens.pnt; sens = rmfield(sens, 'pnt');
      end
    end
    
    if ~isfield(sens, 'chanpos')
      if ismeg
        % sensor description is a MEG sensor-array, containing oriented coils
        [chanpos, chanori, lab] = channelposition(sens, 'channel', 'all');
        % the channel order can be different in the two representations
        [selsens, selpos] = match_str(sens.label, lab);
        sens.chanpos = nan(length(sens.label), 3);
        sens.chanori = nan(length(sens.label), 3);
        % insert the determined position/orientation on the appropriate rows
        sens.chanpos(selsens,:) = chanpos(selpos,:);
        sens.chanori(selsens,:) = chanori(selpos,:);
        if length(selsens)~=length(sens.label)
          warning('cannot determine the position and orientation for all channels');
        end
      else
        % sensor description is something else, EEG/ECoG etc
        % note that chanori will be all NaNs
        [chanpos, chanori, lab] = channelposition(sens, 'channel', 'all');
        % the channel order can be different in the two representations
        [selsens, selpos] = match_str(sens.label, lab);
        sens.chanpos = nan(length(sens.label), 3);
        % insert the determined position/orientation on the appropriate rows
        sens.chanpos(selsens,:) = chanpos(selpos,:);
        if length(selsens)~=length(sens.label)
          warning('cannot determine the position and orientation for all channels');
        end
      end
    end
    
    if ~isfield(sens, 'chantype') || all(strcmp(sens.chantype, 'unknown'))
      if ismeg
        sens.chantype = ft_chantype(sens);
      else
        % FIXME for EEG we have not yet figured out how to deal with this
      end
    end
    
    if ~isfield(sens, 'chanunit') || all(strcmp(sens.chanunit, 'unknown'))
      if ismeg
        sens.chanunit = ft_chanunit(sens);
      else
        % FIXME for EEG we have not yet figured out how to deal with this
      end
    end
    
    if ~isfield(sens, 'unit')
      sens = ft_convert_units(sens);
    end
    
    if any(strcmp(sens.type, {'meg', 'eeg', 'magnetometer', 'electrode', 'unknown'}))
      % this is not sufficiently informative, so better remove it
      % see also http://bugzilla.fcdonders.nl/show_bug.cgi?id=1806
      sens = rmfield(sens, 'type');
    end
    
    if size(sens.chanpos,1)~=length(sens.label) || ...
        isfield(sens, 'tra') && size(sens.tra,1)~=length(sens.label) || ...
        isfield(sens, 'tra') && isfield(sens, 'elecpos') && size(sens.tra,2)~=size(sens.elecpos,1) || ...
        isfield(sens, 'tra') && isfield(sens, 'coilpos') && size(sens.tra,2)~=size(sens.coilpos,1) || ...
        isfield(sens, 'tra') && isfield(sens, 'coilori') && size(sens.tra,2)~=size(sens.coilori,1) || ...
        isfield(sens, 'chanpos') && size(sens.chanpos,1)~=length(sens.label) || ...
        isfield(sens, 'chanori') && size(sens.chanori,1)~=length(sens.label)
      error('inconsistent number of channels in sensor description');
    end
    
    if ismeg
      % ensure that the magnetometer/gradiometer balancing is specified
      if ~isfield(sens, 'balance') || ~isfield(sens.balance, 'current')
        sens.balance.current = 'none';
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'upcoming' % this is under development and expected to become the standard in 2013
    
    % update it to the latest standard version
    sens = ft_datatype_sens(sens, 'version', '2011v2');
    
    % perform some sanity checks
    if ~ft_senstype(sens, 'meg')
      error('unsupported type of gradiometer array "%s"', ft_senstype(sens));
    end
    
    if ~any(strcmp(amplitude, {'T' 'mT' 'uT' 'nT' 'pT' 'fT'}))
      error('unsupported unit of amplitude "%s"', amplitude);
    end
    
    if ~any(strcmp(distance, {'m' 'dm' 'cm' 'mm'}))
      error('unsupported unit of distance "%s"', distance);
    end
    
    sel_m  = ~cellfun(@isempty, regexp(sens.chanunit, '/m$'));
    sel_dm = ~cellfun(@isempty, regexp(sens.chanunit, '/dm$'));
    sel_cm = ~cellfun(@isempty, regexp(sens.chanunit, '/cm$'));
    sel_mm = ~cellfun(@isempty, regexp(sens.chanunit, '/mm$'));
    
    if     strcmp(sens.unit, 'm') && (any(sel_dm) || any(sel_cm) || any(sel_mm))
      error('inconsistent units in input gradiometer');
    elseif strcmp(sens.unit, 'dm') && (any(sel_m) || any(sel_cm) || any(sel_mm))
      error('inconsistent units in input gradiometer');
    elseif strcmp(sens.unit, 'cm') && (any(sel_m) || any(sel_dm) || any(sel_mm))
      error('inconsistent units in input gradiometer');
    elseif strcmp(sens.unit, 'mm') && (any(sel_m) || any(sel_dm) || any(sel_cm))
      error('inconsistent units in input gradiometer');
    end
    
    % update the units of distance
    sens = ft_convert_units(sens, distance);
    
    % update the units of amplitude
    nchan = length(sens.label);
    for i=1:nchan
      if ~isempty(regexp(sens.chanunit{i}, ['/' distance '$'], 'once'))
        % this channel is expressed as amplitude per distance
        sens.tra(i,:)    = sens.tra(i,:) * scalingfactor(sens.chanunit{i}, [amplitude '/' distance]);
        sens.chanunit{i} = [amplitude '/' distance];
      elseif ~isempty(regexp(sens.chanunit{i}, 'T$', 'once'))
        % this channel is expressed as amplitude
        sens.tra(i,:)    = sens.tra(i,:) * scalingfactor(sens.chanunit{i}, amplitude);
        sens.chanunit{i} = amplitude;
      else
        error('unexpected channel unit "%s" in channel %d', i, sens.chanunit{i});
      end
    end
    
    % update the gradiometer scaling
    if strcmp(scaling, 'amplitude')
      for i=1:nchan
        if strcmp(sens.chanunit{i}, [amplitude '/' distance])
          % this channel is expressed as amplitude per distance
          coil = find(abs(sens.tra(i,:))~=0);
          if length(coil)~=2
            error('unexpected number of coils contributing to channel %d', i);
          end
          baseline = norm(sens.coilpos(coil(1),:) - sens.coilpos(coil(2),:));
          sens.tra(i,:)    = sens.tra(i,:)*baseline;  % scale with the baseline distance
          sens.chanunit{i} = amplitude;
        elseif strcmp(sens.chanunit{i}, amplitude)
          % no conversion needed
        else
          error('unexpected channel unit "%s" in channel %d', i, sens.chanunit{i});
        end % if
      end % for
      
    elseif strcmp(scaling, 'amplitude/distance')
      for i=1:nchan
        if strcmp(sens.chanunit{i}, amplitude)
          % this channel is expressed as amplitude
          coil = find(abs(sens.tra(i,:))~=0);
          if length(coil)==1
            % this is a magnetometer channel, no conversion needed
            continue
          elseif length(coil)~=2
            error('unexpected number of coils (%d) contributing to channel %s (%d)', length(coil), sens.label{i}, i);
          end
          baseline = norm(sens.coilpos(coil(1),:) - sens.coilpos(coil(2),:));
          sens.tra(i,:)    = sens.tra(i,:)/baseline; % scale with the baseline distance
          sens.chanunit{i} = [amplitude '/' distance];
        elseif strcmp(sens.chanunit{i}, [amplitude '/' distance])
          % no conversion needed
        else
          error('unexpected channel unit "%s" in channel %d', i, sens.chanunit{i});
        end % if
      end % for
      
    else
      error('incorrect specification of gradiometer scaling');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  otherwise
    error('converting to version %s is not supported', version);
    
end % switch

% this makes the display with the "disp" command look better
sens = sortfieldnames(sens);

% remember the current input and output arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
current_argout = {sens};
previous_argin  = current_argin;
previous_argout = current_argout;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b = sortfieldnames(a)
fn = sort(fieldnames(a));
for i=1:numel(fn)
  b.(fn{i}) = a.(fn{i});
end
