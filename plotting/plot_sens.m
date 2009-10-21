function hs = plot_sens(sens, varargin)

% PLOT_SENS plots the position of the channels in the EEG or MEG sensor array
%
% Use as
%   plot_sens(sens, ...)
% where the first argument is the sensor array as returned by READ_SENS
% or PREPARE_VOL_SENS. 
%
% Optional input arguments should come in key-value pairs and can include
%   'style'    plotting style for the points representing the channels, see plot3 (default = 'k.')
%   'coil'     true/false, plot each individual coil or the channelposition (default = false)
%
% Example
%   sens = read_sens('Subject01.ds');
%   plot_sens(sens, 'style', 'r*')

% Copyright (C) 2009, Robert Oostenveld
%
% $Log: plot_sens.m,v $
% Revision 1.6  2009/07/29 07:07:02  roboos
% added option for plotting all coils
%
% Revision 1.5  2009/06/04 07:30:07  roboos
% added example
%
% Revision 1.4  2009/05/12 18:12:26  roboos
% added handling of hold on/off
%
% Revision 1.3  2009/04/14 19:48:28  roboos
% added keyvalcheck
%
% Revision 1.2  2009/04/08 06:35:05  roboos
% first implementation, covers the use in headmodelplot
%

% get the optional input arguments
keyvalcheck(varargin, 'optional', {'style', 'coil'});
style = keyval('style', varargin); if isempty(style), style = 'k.'; end
coil  = keyval('coil',  varargin); if isempty(coil), coil = false; end

% convert yes/no string into boolean value
coil = istrue(coil);

% everything is added to the current figure
holdflag = ishold;
hold on

if coil
  % simply plot the position of all coils
  hs = plot3(sens.pnt(:,1), sens.pnt(:,2), sens.pnt(:,3), style);
else
  % determine the position of each channel, which is for example the mean of
  % two bipolar electrodes, or the bottom coil of a axial gradiometer
  [chan.pnt, chan.label] = channelposition(sens);
  hs = plot3(chan.pnt(:,1), chan.pnt(:,2), chan.pnt(:,3), style);
end

axis vis3d
axis equal

if ~nargout
  clear hs
end
if ~holdflag
  hold off
end


