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
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

warning('on', 'MATLAB:divideByZero');

% get the optional input arguments
keyvalcheck(varargin, 'optional', {'style', 'coil', 'label'});
style = keyval('style', varargin); if isempty(style), style = 'k.'; end
coil  = keyval('coil',  varargin); if isempty(coil), coil = false; end
label = keyval('label', varargin); if isempty(coil), coil = false; end

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
  
  if ~isempty(label)
    for i=1:length(chan.label)
      switch label
        case {'on', 'yes'}
          str = chan.label{i};
        case {'off', 'no'}
          str = '';
        case {'label' 'labels'}
          str = chan.label{i};
        case {'number' 'numbers'}
          str = num2str(i);
        otherwise
          error('unsupported value for option ''label''');
      end % switch
      text(chan.pnt(i,1), chan.pnt(i,2), chan.pnt(i,3), str);
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


