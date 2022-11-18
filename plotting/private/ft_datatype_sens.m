function [sens] = ft_datatype_sens(sens, varargin)

% FT_DATATYPE_SENS describes the FieldTrip structure that represents an MEG, EEG,
% sEEG, ECoG, or NIRS sensor array. This structure is commonly called "grad" for MEG,
% "elec" for EEG and intranial EEG, "opto" for NIRS, or in general "sens" if it could
% be any one.
%
% For all sensor types a distinction should be made between the channel (i.e. the
% output of the transducer that is A/D converted) and the sensor, which may have some
% spatial extent. For example in MEG gradiometers are comprised of multiple coils and
% with EEG you can have a bipolar channel, where the position of the channel can be
% represented as in between the position of the two electrodes.
%
% The structure for MEG gradiometers and/or magnetometers contains
%    sens.label      = Mx1 cell-array with channel labels
%    sens.chanpos    = Mx3 matrix with channel positions
%    sens.chanori    = Mx3 matrix with channel orientations, used for synthetic planar gradient computation
%    sens.coilpos    = Nx3 matrix with coil positions
%    sens.coilori    = Nx3 matrix with coil orientations
%    sens.tra        = MxN matrix to combine coils into channels
%    sens.balance    = structure containing info about the balancing, See FT_APPLY_MONTAGE
% and optionally
%    sens.chanposold = Mx3 matrix with original channel positions (in case sens.chanpos has been updated to contain NaNs, e.g. after FT_COMPONENTANALYSIS)
%    sens.chanoriold = Mx3 matrix with original channel orientations
%    sens.labelold   = Mx1 cell-array with original channel labels
%
% The structure for EEG, sEEG or ECoG channels contains
%    sens.label      = Mx1 cell-array with channel labels
%    sens.chanpos    = Mx3 matrix with channel positions (often the same as electrode positions)
%    sens.elecpos    = Nx3 matrix with electrode positions
%    sens.tra        = MxN matrix to combine electrodes into channels
% In case sens.tra is not present in the EEG sensor array, the channels
% are assumed to be average referenced.
%
% The structure for NIRS channels contains
%    sens.label         = Mx1 cell-array with channel labels
%    sens.chanpos       = Mx3 matrix with position of the channels (usually halfway the transmitter and receiver)
%    sens.optopos       = Nx3 matrix with the position of individual optodes
%    sens.optotype      = Nx1 cell-array with information about the type of optode (receiver or transmitter)
%    sens.optolabel     = Nx1 cell-array with optode labels
%    sens.wavelength    = 1xK vector of all wavelengths that were used
%    sens.tra           = MxN matrix that specifies for each of the M channels which of the N optodes transmits at which wavelength (positive integer from 1 to K), or receives (negative ingeger from 1 to K)
%
% The following fields apply to MEG, EEG, sEEG and ECoG
%    sens.chantype = Mx1 cell-array with the type of the channel, see FT_CHANTYPE
%    sens.chanunit = Mx1 cell-array with the units of the channel signal, e.g. 'V', 'fT' or 'T/cm', see FT_CHANUNIT
%
% Optional fields:
%    type, unit, fid, chantype, chanunit, coordsys
%
% Historical fields:
%    pnt, pos, ori, pnt1, pnt2, fiberpos, fibertype, fiberlabel, transceiver, transmits, laserstrength
%
% Revision history:
% (2020/latest) Updated the specification of the NIRS sensor definition.
%   Dropped the laserstrength and renamed transmits into tra for consistency.
%
% (2019/latest) Updated the specification of the NIRS sensor definition.
%   Use "opto" instead of "fibers", see http://bit.ly/33WaqWU for details.
%
% (2016) The chantype and chanunit have become required fields.
%  Original channel details are specified with the suffix "old" rather than "org".
%  All numeric values are represented in double precision.
%  It is possible to convert the amplitude and distance units (e.g. from T to fT and
%  from m to mm) and it is possible to express planar and axial gradiometer channels
%  either in units of amplitude or in units of amplitude/distance (i.e. proper
%  gradient).
%
% (2011v2) The chantype and chanunit have been added for MEG.
%
% (2011v1) To facilitate determining the position of channels (e.g. for plotting)
%  in case of balanced MEG or bipolar EEG, an explicit distinction has been made
%  between chanpos+chanori and coilpos+coilori (for MEG) and chanpos and elecpos
%  (for EEG). The pnt and ori fields are removed.
%
% (2010) Added support for bipolar or otherwise more complex linear combinations
%  of EEG electrodes using sens.tra, similar to MEG.
%
% (2009) Noise reduction has been added for MEG systems in the balance field.
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

% Copyright (C) 2011-2020, Robert Oostenveld & Jan-Mathijs Schoffelen
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

% undocumented options for the 2016 format and up
%   amplitude     = string, can be 'T' or 'fT'
%   distance      = string, can be 'm', 'cm' or 'mm'
%   scaling       = string, can be 'amplitude' or 'amplitude/distance'

% these are for speeding up subsequent calls with the same input arguments
persistent previous_argin previous_argout

current_argin = [{sens} varargin];
if isequal(current_argin, previous_argin)
  % don't do the whole checking again, but return the previous output from cache
  sens = previous_argout{1};
  return
end

% get the optional input arguments, which should be specified as key-value pairs
version   = ft_getopt(varargin, 'version', 'latest');
amplitude = ft_getopt(varargin, 'amplitude'); % should be 'V' 'uV' 'T' 'mT' 'uT' 'nT' 'pT' 'fT'
distance  = ft_getopt(varargin, 'distance');  % should be 'm' 'dm' 'cm' 'mm'
scaling   = ft_getopt(varargin, 'scaling');   % should be 'amplitude' or 'amplitude/distance', the default depends on the senstype

if ~isempty(amplitude) && ~any(strcmp(amplitude, {'V' 'uV' 'T' 'mT' 'uT' 'nT' 'pT' 'fT'}))
  ft_error('unsupported unit of amplitude "%s"', amplitude);
end

if ~isempty(distance) && ~any(strcmp(distance, {'m' 'dm' 'cm' 'mm'}))
  ft_error('unsupported unit of distance "%s"', distance);
end

if strcmp(version, 'latest')
  version = '2020';
end

if isempty(sens)
  return;
end

% this is needed further down
nchan = length(sens.label);

% these are used at multiple places, therefore we determine them only once
if isfield(sens, 'coilpos')
  ismeg = true;
  iseeg = true;
  isnirs = false;
elseif isfield(sens, 'elecpos')
  ismeg = false;
  iseeg = true;
  isnirs = false;
elseif isfield(sens, 'optopos')
  ismeg = false;
  iseeg = false;
  isnirs = true;
else
  % doing it this way takes a lot more CPU time
  ismeg  = ft_senstype(sens, 'meg');
  iseeg  = ft_senstype(sens, 'eeg');
  isnirs = ft_senstype(sens, 'nirs');
end

switch version
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case '2020'
    % update it to the previous standard version
    new_argin = ft_setopt(varargin, 'version', '2019');
    sens      = ft_datatype_sens(sens, new_argin{:});
    if isfield(sens, 'coordsys')
      sens = fixcoordsys(sens);
    end

    if isnirs
      sens = renamefields(sens, 'transmits', 'tra'); % this makes it more consistent with EEG and MEG
      sens = removefields(sens, {'laserstrength'});
      
      if isfield(sens, 'optotype') && isfield(sens, 'tra')
        % the read_artinis_oxy3 file returns the wrong sign for the receivers and transmitters
        % but since it is p-code, it cannot be fixed
        tmp = sens.tra(:, strcmp(sens.optotype, 'transmitter')); correctT = all(tmp(:)>=0);
        tmp = sens.tra(:, strcmp(sens.optotype, 'receiver'   )); correctR = all(tmp(:)<=0);
        if ~correctT
          ft_warning('flipping sign for transmitters');
          sens.tra(:, strcmp(sens.optotype, 'transmitter')) = -sens.tra(:, strcmp(sens.optotype, 'transmitter'));
        end
        if ~correctR
          ft_warning('flipping sign for receivers');
          sens.tra(:, strcmp(sens.optotype, 'receiver')) = -sens.tra(:, strcmp(sens.optotype, 'receiver'));
        end
      end
      
    end % ifnirs
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case '2019'
    % update it to the previous standard version
    new_argin = ft_setopt(varargin, 'version', '2016');
    sens      = ft_datatype_sens(sens, new_argin{:});
    
    if isnirs
      % rename some NIRS field names
      sens = renamefields(sens, 'fiberpos', 'optopos');
      sens = renamefields(sens, 'fibertype', 'optotype');
      sens = renamefields(sens, 'fiberlabel', 'optolabel');
      sens = renamefields(sens, 'transceiver', 'transmits');
      
      % these might be present due to the reading functions but are not part of the user/technical documentation yet, so better not include them for now
      % especially the chanunit field needs some careful thought when converting between optical densities and chromophore concentrations.
      sens = removefields(sens, {'chantype', 'chanunit'});
    end
    
    % ensure all positions are represented as 3D, not 2D
    fn = {'chanpos', 'optopos', 'elecpos', 'coilpos'};
    for i=1:numel(fn)
      if isfield(sens, fn{i}) && size(sens.(fn{i}),2)==2
        ft_notice('converting %s from 2D to 3D', fn{i});
        sens.(fn{i})(:,3) = 0;
      end
    end % for
    
    if isfield(sens, 'optolabel')
      if length(sens.optolabel)~=length(unique(sens.optolabel))
        ft_warning('non-unique optode labels detected');
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case '2016'
    % update it to the previous standard version
    new_argin = ft_setopt(varargin, 'version', '2011v2');
    sens      = ft_datatype_sens(sens, new_argin{:});
    
    % rename from org to old (reverse = false)
    sens = fixoldorg(sens, false);
    
    % ensure that all numbers are represented in double precision
    % this only affects the top-level fields and does not recurse
    sens = ft_struct2double(sens, 1);
    
    % in version 2011v2 this was optional, now it is required
    if ~isfield(sens, 'chantype') || all(strcmp(sens.chantype, 'unknown'))
      sens.chantype = ft_chantype(sens);
    end
    
    % in version 2011v2 this was optional, now it is required
    if ~isfield(sens, 'chanunit') || all(strcmp(sens.chanunit, 'unknown'))
      sens.chanunit = ft_chanunit(sens);
    end
    
    if ~isempty(distance)
      % update the units of distance, this also updates the tra matrix
      sens = ft_convert_units(sens, distance);
    else
      % determine the default, this may be needed to set the scaling
      distance = sens.unit;
    end
    
    if ~isempty(amplitude) && isfield(sens, 'tra')
      % update the tra matrix for the units of amplitude, this ensures that
      % the leadfield values remain consistent with the units
      for i=1:nchan
        if ~isempty(regexp(sens.chanunit{i}, 'm$', 'once'))
          % this channel is expressed as amplitude per distance
          sens.tra(i,:)    = sens.tra(i,:) * ft_scalingfactor(sens.chanunit{i}, [amplitude '/' distance]);
          sens.chanunit{i} = [amplitude '/' distance];
        elseif ~isempty(regexp(sens.chanunit{i}, '[T|V]$', 'once'))
          % this channel is expressed as amplitude
          sens.tra(i,:)    = sens.tra(i,:) * ft_scalingfactor(sens.chanunit{i}, amplitude);
          sens.chanunit{i} = amplitude;
        else
          ft_error('unexpected channel unit "%s" in channel %d', sens.chanunit{i}, i);
        end
      end
    else
      % determine the default amplityde, this may be needed to set the scaling
      if any(~cellfun(@isempty, regexp(sens.chanunit, '^T')))
        % one of the channel units starts with T
        amplitude = 'T';
      elseif any(~cellfun(@isempty, regexp(sens.chanunit, '^fT')))
        % one of the channel units starts with fT
        amplitude = 'fT';
      elseif any(~cellfun(@isempty, regexp(sens.chanunit, '^V')))
        % one of the channel units starts with V
        amplitude = 'V';
      elseif any(~cellfun(@isempty, regexp(sens.chanunit, '^uV')))
        % one of the channel units starts with uV
        amplitude = 'uV';
      else
        % this unknown amplitude will cause a problem if the scaling needs to be changed between amplitude and amplitude/distance
        amplitude = 'unknown';
      end
    end
    
    % perform some sanity checks
    if ismeg
      sel_m  = ~cellfun(@isempty, regexp(sens.chanunit, '/m$'));
      sel_dm = ~cellfun(@isempty, regexp(sens.chanunit, '/dm$'));
      sel_cm = ~cellfun(@isempty, regexp(sens.chanunit, '/cm$'));
      sel_mm = ~cellfun(@isempty, regexp(sens.chanunit, '/mm$'));
      
      if     strcmp(sens.unit, 'm') && (any(sel_dm) || any(sel_cm) || any(sel_mm))
        ft_error('inconsistent units in input gradiometer');
      elseif strcmp(sens.unit, 'dm') && (any(sel_m) || any(sel_cm) || any(sel_mm))
        ft_error('inconsistent units in input gradiometer');
      elseif strcmp(sens.unit, 'cm') && (any(sel_m) || any(sel_dm) || any(sel_mm))
        ft_error('inconsistent units in input gradiometer');
      elseif strcmp(sens.unit, 'mm') && (any(sel_m) || any(sel_dm) || any(sel_cm))
        ft_error('inconsistent units in input gradiometer');
      end
      
      if ~strcmp(amplitude, 'unknown') && ~strcmp(distance, 'unknown')
        
        % the default should be "amplitude/distance" for neuromag and "amplitude" for all others
        if isempty(scaling)
          if ft_senstype(sens, 'neuromag') && ~any(contains(sens.chanunit, '/'))
            ft_warning('assuming that the default scaling should be amplitude/distance rather than amplitude');
            scaling = 'amplitude/distance';
          elseif ft_senstype(sens, 'yokogawa440') && any(contains(sens.chanunit, '/'))
            ft_warning('assuming that the default scaling should be amplitude rather than amplitude/distance');
            scaling = 'amplitude';
          end
        end
        
        % update the gradiometer scaling
        if strcmp(scaling, 'amplitude') && isfield(sens, 'tra')
          for i=1:nchan
            if strcmp(sens.chanunit{i}, [amplitude '/' distance])
              % this channel is expressed as amplitude per distance
              coil = find(abs(sens.tra(i,:))~=0);
              if length(coil)~=2
                ft_error('unexpected number of coils contributing to channel %d', i);
              end
              baseline         = norm(sens.coilpos(coil(1),:) - sens.coilpos(coil(2),:));
              sens.tra(i,:)    = sens.tra(i,:)*baseline;  % scale with the baseline distance
              sens.chanunit{i} = amplitude;
            elseif strcmp(sens.chanunit{i}, amplitude)
              % no conversion needed
            elseif isfield(sens, 'balance') && strcmp(sens.balance.current, 'comp')
              % no conversion needed
            elseif isfield(sens, 'balance') && strcmp(sens.balance.current, 'planar')
              % no conversion needed
            else
              % see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=3144
              ft_warning('unexpected channel unit "%s" in channel %d', sens.chanunit{i}, i);
            end % if
          end % for nchan
          
        elseif strcmp(scaling, 'amplitude/distance') && isfield(sens, 'tra')
          for i=1:nchan
            if strcmp(sens.chanunit{i}, amplitude)
              % this channel is expressed as amplitude
              coil = find(abs(sens.tra(i,:))~=0);
              if length(coil)==1 || strcmp(sens.chantype{i}, 'megmag')
                % this is a magnetometer channel, no conversion needed
                continue
              elseif length(coil)~=2
                ft_error('unexpected number of coils (%d) contributing to channel %s (%d)', length(coil), sens.label{i}, i);
              end
              baseline         = norm(sens.coilpos(coil(1),:) - sens.coilpos(coil(2),:));
              sens.tra(i,:)    = sens.tra(i,:)/baseline; % scale with the baseline distance
              sens.chanunit{i} = [amplitude '/' distance];
            elseif strcmp(sens.chanunit{i}, [amplitude '/' distance])
              % no conversion needed
            else
              % see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=3144
              ft_warning('unexpected channel unit "%s" in channel %d', sens.chanunit{i}, i);
            end % if
          end % for nchan
          
        end % if strcmp scaling
      end % if amplitude and scaling are not unknown
      
    else
      sel_m  = ~cellfun(@isempty, regexp(sens.chanunit, '/m$'));
      sel_dm = ~cellfun(@isempty, regexp(sens.chanunit, '/dm$'));
      sel_cm = ~cellfun(@isempty, regexp(sens.chanunit, '/cm$'));
      sel_mm = ~cellfun(@isempty, regexp(sens.chanunit, '/mm$'));
      if any(sel_m | sel_dm | sel_cm | sel_mm)
        ft_error('scaling of amplitude/distance has not been considered yet for EEG');
      end
      
    end % if iseeg or ismeg
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case '2011v2'
    
    % rename from old to org (reverse = true)
    sens = fixoldorg(sens, true);
    
    if ~isempty(amplitude) || ~isempty(distance) || ~isempty(scaling)
      ft_warning('amplitude, distance and scaling are not supported for version "%s"', version);
    end
    
    % This speeds up subsequent calls to ft_senstype and channelposition.
    % However, if it is not more precise than MEG or EEG, don't keep it in
    % the output (see further down).
    if ~isfield(sens, 'type')
      sens.type = ft_senstype(sens);
    end
    
    if isfield(sens, 'pnt')
      if ismeg
        % sensor description is a MEG sensor-array, containing oriented coils
        sens.coilpos = sens.pnt; sens = rmfield(sens, 'pnt');
        sens.coilori = sens.ori; sens = rmfield(sens, 'ori');
      else
        % sensor description is something else, EEG/ECoG/sEEG, etc
        sens.elecpos = sens.pnt; sens = rmfield(sens, 'pnt');
      end
    end

    if isfield(sens, 'pos')
      if ismeg
        % sensor description is a MEG sensor-array, containing oriented coils
        sens.coilpos = sens.pos; sens = rmfield(sens, 'pos');
        sens.coilori = sens.ori; sens = rmfield(sens, 'ori');
      else
        % sensor description is something else, EEG/ECoG/sEEG, etc
        sens.elecpos = sens.pos; sens = rmfield(sens, 'pos');
      end
    end

    if ~isfield(sens, 'chanpos')
      if ismeg
        % sensor description is a MEG sensor-array, containing oriented coils
        [chanpos, chanori, lab] = channelposition(sens);
        % the channel order can be different in the two representations
        [selsens, selpos] = match_str(sens.label, lab);
        sens.chanpos = nan(length(sens.label), 3);
        sens.chanori = nan(length(sens.label), 3);
        % insert the determined position/orientation on the appropriate rows
        sens.chanpos(selsens,:) = chanpos(selpos,:);
        sens.chanori(selsens,:) = chanori(selpos,:);
        if length(selsens)~=length(sens.label)
          ft_warning('cannot determine the position and orientation for all channels');
        end
      else
        % sensor description is something else, EEG/ECoG etc
        % note that chanori will be all NaNs
        [chanpos, chanori, lab] = channelposition(sens);
        % the channel order can be different in the two representations
        [selsens, selpos] = match_str(sens.label, lab);
        sens.chanpos = nan(length(sens.label), 3);
        % insert the determined position/orientation on the appropriate rows
        sens.chanpos(selsens,:) = chanpos(selpos,:);
        if length(selsens)~=length(sens.label)
          ft_warning('cannot determine the position and orientation for all channels');
        end
      end
    end
    
    if ~isfield(sens, 'chantype') || all(strcmp(sens.chantype, 'unknown'))
      if ismeg
        sens.chantype = ft_chantype(sens);
      else
        % for EEG it is not required
      end
    end
    
    if ~isfield(sens, 'unit')
      % this should be done prior to calling ft_chanunit, since ft_chanunit uses this for planar neuromag channels
      sens = ft_determine_units(sens);
    end
    
    if ~isfield(sens, 'chanunit') || all(strcmp(sens.chanunit, 'unknown'))
      if ismeg
        sens.chanunit = ft_chanunit(sens);
      else
        % for EEG it is not required
      end
    end
    
    if any(strcmp(sens.type, {'meg', 'eeg', 'magnetometer', 'electrode', 'unknown'}))
      % this is not sufficiently informative, so better remove it
      % see also http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=1806
      sens = rmfield(sens, 'type');
    end
    
    if size(sens.chanpos,1)~=length(sens.label) || ...
        isfield(sens, 'tra') && size(sens.tra,1)~=length(sens.label) || ...
        isfield(sens, 'tra') && isfield(sens, 'elecpos') && size(sens.tra,2)~=size(sens.elecpos,1) || ...
        isfield(sens, 'tra') && isfield(sens, 'coilpos') && size(sens.tra,2)~=size(sens.coilpos,1) || ...
        isfield(sens, 'tra') && isfield(sens, 'coilori') && size(sens.tra,2)~=size(sens.coilori,1) || ...
        isfield(sens, 'chanpos') && size(sens.chanpos,1)~=length(sens.label) || ...
        isfield(sens, 'chanori') && size(sens.chanori,1)~=length(sens.label)
      ft_error('inconsistent number of channels in sensor description');
    end
    
    if ismeg
      % ensure that the magnetometer/gradiometer balancing is specified
      if ~isfield(sens, 'balance') || ~isfield(sens.balance, 'current')
        sens.balance.current = 'none';
      end
      
      % try to add the chantype and chanunit to the CTF G1BR montage
      if isfield(sens, 'balance') && isfield(sens.balance, 'G1BR') && ~isfield(sens.balance.G1BR, 'chantype')
        sens.balance.G1BR.chantypeorg = repmat({'unknown'}, size(sens.balance.G1BR.labelorg));
        sens.balance.G1BR.chanunitorg = repmat({'unknown'}, size(sens.balance.G1BR.labelorg));
        sens.balance.G1BR.chantypenew = repmat({'unknown'}, size(sens.balance.G1BR.labelnew));
        sens.balance.G1BR.chanunitnew = repmat({'unknown'}, size(sens.balance.G1BR.labelnew));
        % the synthetic gradient montage does not fundamentally change the chantype or chanunit
        [sel1, sel2] = match_str(sens.balance.G1BR.labelorg, sens.label);
        sens.balance.G1BR.chantypeorg(sel1) = sens.chantype(sel2);
        sens.balance.G1BR.chanunitorg(sel1) = sens.chanunit(sel2);
        [sel1, sel2] = match_str(sens.balance.G1BR.labelnew, sens.label);
        sens.balance.G1BR.chantypenew(sel1) = sens.chantype(sel2);
        sens.balance.G1BR.chanunitnew(sel1) = sens.chanunit(sel2);
      end
      
      % idem for G2BR
      if isfield(sens, 'balance') && isfield(sens.balance, 'G2BR') && ~isfield(sens.balance.G2BR, 'chantype')
        sens.balance.G2BR.chantypeorg = repmat({'unknown'}, size(sens.balance.G2BR.labelorg));
        sens.balance.G2BR.chanunitorg = repmat({'unknown'}, size(sens.balance.G2BR.labelorg));
        sens.balance.G2BR.chantypenew = repmat({'unknown'}, size(sens.balance.G2BR.labelnew));
        sens.balance.G2BR.chanunitnew = repmat({'unknown'}, size(sens.balance.G2BR.labelnew));
        % the synthetic gradient montage does not fundamentally change the chantype or chanunit
        [sel1, sel2] = match_str(sens.balance.G2BR.labelorg, sens.label);
        sens.balance.G2BR.chantypeorg(sel1) = sens.chantype(sel2);
        sens.balance.G2BR.chanunitorg(sel1) = sens.chanunit(sel2);
        [sel1, sel2] = match_str(sens.balance.G2BR.labelnew, sens.label);
        sens.balance.G2BR.chantypenew(sel1) = sens.chantype(sel2);
        sens.balance.G2BR.chanunitnew(sel1) = sens.chanunit(sel2);
      end
      
      % idem for G3BR
      if isfield(sens, 'balance') && isfield(sens.balance, 'G3BR') && ~isfield(sens.balance.G3BR, 'chantype')
        sens.balance.G3BR.chantypeorg = repmat({'unknown'}, size(sens.balance.G3BR.labelorg));
        sens.balance.G3BR.chanunitorg = repmat({'unknown'}, size(sens.balance.G3BR.labelorg));
        sens.balance.G3BR.chantypenew = repmat({'unknown'}, size(sens.balance.G3BR.labelnew));
        sens.balance.G3BR.chanunitnew = repmat({'unknown'}, size(sens.balance.G3BR.labelnew));
        % the synthetic gradient montage does not fundamentally change the chantype or chanunit
        [sel1, sel2] = match_str(sens.balance.G3BR.labelorg, sens.label);
        sens.balance.G3BR.chantypeorg(sel1) = sens.chantype(sel2);
        sens.balance.G3BR.chanunitorg(sel1) = sens.chanunit(sel2);
        [sel1, sel2] = match_str(sens.balance.G3BR.labelnew, sens.label);
        sens.balance.G3BR.chantypenew(sel1) = sens.chantype(sel2);
        sens.balance.G3BR.chanunitnew(sel1) = sens.chanunit(sel2);
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  otherwise
    ft_error('converting to version %s is not supported', version);
    
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
