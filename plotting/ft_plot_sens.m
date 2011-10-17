function hs = ft_plot_sens(sens, varargin)

% FT_PLOT_SENS plots the position of the channels in the EEG or MEG sensor array
%
% Use as
%   ft_plot_sens(sens, ...)
% where the first argument is the sensor array as returned by READ_SENS
% or PREPARE_VOL_SENS.
%
% Optional input arguments should come in key-value pairs and can include
%   'style'    plotting style for the points representing the channels, see plot3 (default = 'k.')
%   'coil'     true/false, plot each individual coil or the channelposition (default = false)
%   'label'    show the label, can be 'off', 'label', 'number' (default = 'off')
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
sens = fixsens(sens);

% get the optional input arguments
style = ft_getopt(varargin, 'style',  'k.');
coil  = ft_getopt(varargin, 'coil',   false);
label = ft_getopt(varargin, 'label',  'off');

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
  
  if ~isempty(label)
    for i=1:length(sens.label)
      switch label
        case {'on', 'yes'}
          str = sens.label{i};
        case {'off', 'no'}
          str = '';
        case {'label' 'labels'}
          str = sens.label{i};
        case {'number' 'numbers'}
          str = num2str(i);
        otherwise
          error('unsupported value for option ''label''');
      end % switch
      text(sens.chanpos(i,1), sens.chanpos(i,2), sens.chanpos(i,3), str);
    end % for
  end % if
  
end

axis vis3d
axis equal

if ~nargout
  clear hs
end
if ~holdflag
  hold off
end

warning(ws); %revert to original state
