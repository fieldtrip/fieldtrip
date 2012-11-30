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
%    sens.elecpos  = Nx3 matrix with coil positions
% In case sens.tra is not present in the EEG sensor array, the channels
% are assumed to be average referenced.
%
% The following fields are optional
%    sens.type     = string with the MEG or EEG acquisition system, see FT_SENSTYPE
%    sens.chantype = Mx1 cell-array with the type of the channel, see FT_CHANTYPE
%    sens.chanunit = Mx1 cell-array with the units of the channel signal, e.g. 'T', 'fT' or 'fT/cm'
%
% Revision history:
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

% Copyright (C) 2011, Robert Oostenveld & Jan-Mathijs Schoffelen
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

% get the optional input arguments, which should be specified as key-value pairs
version = ft_getopt(varargin, 'version', 'latest');

if strcmp(version, 'latest')
  version = '2011v2';
end

if isempty(sens)
  return;
end

switch version
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case '2011v2'
    
    % This speeds up subsequent calls to ft_senstype and channelposition.
    % However, if it is not more precise than MEG or EEG, don't keep it in
    % the output (see further down).
    if ~isfield(sens, 'type')
      sens.type = ft_senstype(sens);
    end
    
    % there are many cases which deal with either eeg or meg
    isgrad = ft_senstype(sens, 'meg');
    
    if isfield(sens, 'pnt') || isfield(sens, 'ori')
      if isgrad
        % sensor description is a MEG sensor-array, containing oriented coils
        [chanpos, chanori, chanlab] = channelposition(sens, 'channel', 'all');
        sens.coilori = sens.ori; sens = rmfield(sens, 'ori');
        sens.coilpos = sens.pnt; sens = rmfield(sens, 'pnt');
        sens.chanpos = chanpos;
        sens.chanori = chanori;
      else
        % sensor description is something else, EEG/ECoG etc
        % note that chanori will be all NaNs
        [chanpos, chanori, chanlab] = channelposition(sens, 'channel', 'all');
        sens.elecpos = chanpos; sens = rmfield(sens, 'pnt');
        sens.chanpos = chanpos;
      end
    end
    
    if ~isfield(sens, 'chanpos')
      if isgrad
        % sensor description is a MEG sensor-array, containing oriented coils
        [chanpos, chanori, lab] = channelposition(sens, 'channel', 'all');
        sens.chanpos = chanpos;
        sens.chanori = chanori;
      else
        % sensor description is something else, EEG/ECoG etc
        [chanpos, lab] = channelposition(sens, 'channel', 'all');
        sens.chanpos   = chanpos;
      end
    end
    
    if ~isfield(sens, 'chantype') || all(strcmp(sens.chantype, 'unknown'))
      if isgrad
        sens.chantype = ft_chantype(sens);
      else
        % FIXME for EEG we have not yet figured out how to deal with this
      end
    end
    
    if ~isfield(sens, 'chanunit') || all(strcmp(sens.chanunit, 'unknown'))
      if isgrad
        sens.chanunit = ft_chanunit(sens);
      else
        % FIXME for EEG we have not yet figured out how to deal with this
      end
    end
    
    if ~isfield(sens, 'unit')
      sens = ft_convert_units(sens);
    end
    
    if any(strcmp(sens.type, {'meg', 'eeg', 'magnetometer', 'electrode'}))
      % this is not sufficiently informative, so better remove it
      % see also http://bugzilla.fcdonders.nl/show_bug.cgi?id=1806
      sens = rmfield(sens, 'type');
    end
    
    if size(sens.chanpos,1)~=length(sens.label) || ...
       isfield(sens, 'tra') && size(sens.tra,1)~=length(sens.label) || ...
       isfield(sens, 'tra') && isfield(sens, 'elecpos') && size(sens.tra,2)~=size(sens.elecpos,1) || ...
       isfield(sens, 'tra') && isfield(sens, 'coilpos') && size(sens.tra,2)~=size(sens.coilpos,1) || ...
       isfield(sens, 'tra') && isfield(sens, 'coilori') && size(sens.tra,2)~=size(sens.coilori,1)
     error('inconsistent number of channels in sensor description');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  otherwise
    error('converting to version %s is not supported', version);
    
end % switch

% this makes the display with the "disp" command look better
sens = sortfieldnames(sens);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b = sortfieldnames(a)
fn = sort(fieldnames(a));
for i=1:numel(fn)
  b.(fn{i}) = a.(fn{i});
end

