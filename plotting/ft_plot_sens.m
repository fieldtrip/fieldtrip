function hs = ft_plot_sens(sens, varargin)

% FT_PLOT_SENS plots the position of the channels in the EEG or MEG sensor array
%
% Use as
%   ft_plot_sens(sens, ...)
% where the first argument is the sensor array as returned by READ_SENS
% or PREPARE_VOL_SENS.
%
% Optional input arguments should come in key-value pairs and can include
%   style    = plotting style for the points representing the channels, see plot3 (default = 'k.')
%   coil     = true/false, plot each individual coil or the channelposition (default = false)
%   label    = show the label, can be 'off', 'label', 'number' (default = 'off')
%   chantype = string or cell-array with strings, for example {'meg', 'megref'} (default is all)
%
% Example
%   sens = ft_read_sens('Subject01.ds');
%   ft_plot_sens(sens, 'style', 'r*')

% Copyright (C) 2009, Robert Oostenveld
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

ws = warning('on', 'MATLAB:divideByZero');

% ensure that the sensor description is up-to-date (Aug 2011)
sens = ft_datatype_sens(sens);

% get the optional input arguments
style    = ft_getopt(varargin, 'style',  'k.');
coil     = ft_getopt(varargin, 'coil',   false);
label    = ft_getopt(varargin, 'label',  'off');
chantype = ft_getopt(varargin, 'chantype');

% select a subset of channels to be plotted
if ~isempty(chantype)
  
  if ischar(chantype)
    chantype = {chantype};
  end
  chansel = match_str(sens.chantype, chantype);
  
  % remove the balancing from the sensor definition, e.g. 3rd order gradients, PCA-cleaned data or ICA projections
  sens = undobalancing(sens);
  
  % remove the channels that are not selected
  sens.chanpos = sens.chanpos(chansel,:);
  sens.label   = sens.label(chansel);
  
  % remove the magnetometer and gradiometer coils that are not in one of the selected channels
  if isfield(sens, 'tra') && isfield(sens, 'coilpos')
    sens.tra     = sens.tra(chansel,:);
    coilsel      = any(sens.tra~=0,1);
    sens.coilpos = sens.coilpos(coilsel,:);
    sens.coilori = sens.coilori(coilsel,:);
    sens.tra     = sens.tra(:,coilsel);
  end
  
  % FIXME note that I have not tested this on any complicated electrode definitions
  % remove the electrodes that are not in one of the selected channels
  if isfield(sens, 'tra') && isfield(sens, 'elecpos')
    sens.tra     = sens.tra(chansel,:);
    elecsel      = any(sens.tra~=0,1);
    sens.elecpos = sens.elecpos(elecsel,:);
    sens.tra     = sens.tra(:,elecsel);
  end
  
end % selecting channels and coils

% everything is added to the current figure
holdflag = ishold;
if ~holdflag
  hold on
end

if istrue(coil)
  % simply plot the position of all coils or electrodes
  if isfield(sens, 'coilpos')
    pnt = sens.coilpos;
  elseif isfield(sens, 'elecpos')
    pnt = sens.elecpos;
  end
  
  hs = plot3(pnt(:,1), pnt(:,2), pnt(:,3), style);
  
else
  % determine the position of each channel, which is for example the mean of
  % two bipolar electrodes, or the bottom coil of a axial gradiometer
  hs = plot3(sens.chanpos(:,1), sens.chanpos(:,2), sens.chanpos(:,3), style);
  
  if ~isempty(label) && ~any(strcmp(label, {'off', 'no'}))
    for i=1:length(sens.label)
      switch label
        case {'on', 'yes'}
          str = sens.label{i};
        case {'label' 'labels'}
          str = sens.label{i};
        case {'number' 'numbers'}
          str = num2str(i);
        otherwise
          error('unsupported value for option ''label''');
      end % switch
      text(sens.chanpos(i,1), sens.chanpos(i,2), sens.chanpos(i,3), str);
    end % for
  end % if empty or off/no
  
end

axis vis3d
axis equal

if ~nargout
  clear hs
end
if ~holdflag
  hold off
end

warning(ws); % revert to original state
