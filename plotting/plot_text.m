function [varargout] = plot_text(X, Y, str, varargin)

% PLOT_TEXT helper function for plotting text, which can also be used in
% combination with the multiple channel layout display in FieldTrip.
%
% Use as
%   plot_text(X, Y, ...)
% where optional input arguments should come in key-value pairs and may
% include
%   hpos
%   vpos
%   width
%   height
%   hlim
%   vlim
%   Color
%   FontSize
%   FontName
%   HorizontalAlignment

% Copyrights (C) 2009, Robert Oostenveld
%
% $Log: plot_text.m,v $
% Revision 1.6  2009/10/07 13:52:03  roboos
% updated documentation
%
% Revision 1.5  2009/08/05 08:52:09  roboos
% added HorizontalAlignment option
% use CamelCase for options that are passed on to the default set() function
%
% Revision 1.4  2009/06/02 15:40:36  giopia
% added varargout to pass handle
%
% Revision 1.3  2009/04/14 19:48:28  roboos
% added keyvalcheck
%
% Revision 1.2  2009/04/14 14:31:08  roboos
% many small changes to make it fully functional
%

% get the optional input arguments
keyvalcheck(varargin, 'optional', {'hpos', 'vpos', 'width', 'height', 'hlim', 'vlim', 'Color', 'FontSize', 'FontName', 'HorizontalAlignment'});
hpos        = keyval('hpos',      varargin);
vpos        = keyval('vpos',      varargin);
width       = keyval('width',     varargin);
height      = keyval('height',    varargin);
hlim        = keyval('hlim',      varargin);
vlim        = keyval('vlim',      varargin);
Color       = keyval('Color',     varargin);  if isempty(Color), Color = 'k'; end
FontSize    = keyval('FontSize',  varargin);
FontName    = keyval('FontName',  varargin);
HorizontalAlignment = keyval('HorizontalAlignment',  varargin); if isempty(HorizontalAlignment), HorizontalAlignment = 'center'; end

if isempty(hlim) && isempty(vlim) && isempty(hpos) && isempty(vpos) && isempty(height) && isempty(width)
  % no scaling is needed, the input X and Y are already fine
  % use a shortcut to speed up the plotting
  
else
  % use the full implementation
  abc = axis;
  if isempty(hlim)
    hlim = abc([1 2]);
  end
  
  if isempty(vlim)
    vlim = abc([3 4]);
  end
  
  if isempty(hpos);
    hpos = (hlim(1)+hlim(2))/2;
  end
  
  if isempty(vpos);
    vpos = (vlim(1)+vlim(2))/2;
  end
  
  if isempty(width),
    width = hlim(2)-hlim(1);
  end
  
  if isempty(height),
    height = vlim(2)-vlim(1);
  end
  
  % first shift the horizontal axis to zero
  X = X - (hlim(1)+hlim(2))/2;
  % then scale to length 1
  X = X ./ (hlim(2)-hlim(1));
  % then scale to the new width
  X = X .* width;
  % then shift to the new horizontal position
  X = X + hpos;
  
  % first shift the vertical axis to zero
  Y = Y - (vlim(1)+vlim(2))/2;
  % then scale to length 1
  Y = Y ./ (vlim(2)-vlim(1));
  % then scale to the new width
  Y = Y .* height;
  % then shift to the new vertical position
  Y = Y + vpos;
  
end % shortcut

h = text(X, Y, str);
set(h, 'HorizontalAlignment', HorizontalAlignment);
set(h, 'Color', Color);
if ~isempty(FontSize), set(h, 'FontSize', FontSize); end
if ~isempty(FontName), set(h, 'FontName', FontName); end

% the (optional) output is the handle
if nargout == 1;
  varargout{1} = h;
end
