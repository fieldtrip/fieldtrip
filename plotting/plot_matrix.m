function plot_matrix(varargin)

% PLOT_MATRIX
%
% Use as
%   plot_matrix(C, ...)
% where C is a 2 dimensional MxN matrix, or
%   plot_matrix(X, Y, C, ...)
% where X and Y describe the 1xN horizontal and 1xM vertical axes 
% respectively.
%
% Additional options should be specified in key-value pairs and can be
%   'hpos'
%   'vpos'
%   'width'
%   'height'
%   'hlim'
%   'vlim'
%   'clim'
%   'box'     can be 'yes' or 'no'
%
% Example use
%   plot_matrix(randn(30,50), 'width', 1, 'height', 1, 'hpos', 0, 'vpos', 0)

% Copyrights (C) 2009, Robert Oostenveld
%
% $Log: plot_matrix.m,v $
% Revision 1.6  2009/07/14 16:13:53  roboos
% added implementation for clim and highlight
%
% Revision 1.5  2009/06/16 07:51:51  crimic
% small change
%
% Revision 1.4  2009/04/15 20:02:17  roboos
% changed input parsing, fixed 1-pixel offset, added box option
%
% Revision 1.3  2009/04/14 14:31:08  roboos
% many small changes to make it fully functional
%

if nargin>2 && all(cellfun(@isnumeric, varargin(1:3)))
  % the function was called like imagesc(x, y, c, ...)
  hdat = varargin{1};
  vdat = varargin{2};
  cdat = varargin{3};
  varargin = varargin(4:end);
else
  % the function was called like plot(c, ...)
  cdat = varargin{1};
  vdat = 1:size(cdat,1);
  hdat = 1:size(cdat,2);
  varargin = varargin(2:end);
end

% get the optional input arguments
keyvalcheck(varargin, 'optional', {'hpos', 'vpos', 'width', 'height', 'hlim', 'vlim', 'clim', 'box','highlight','highlightstyle'});
hpos   = keyval('hpos',   varargin);
vpos   = keyval('vpos',   varargin);
width  = keyval('width',  varargin);
height = keyval('height', varargin);
hlim   = keyval('hlim',   varargin);
vlim   = keyval('vlim',   varargin);
clim   = keyval('clim',   varargin);
box    = keyval('box',    varargin); if isempty(box), box = false; end
highlight      = keyval('highlight',       varargin);
highlightstyle = keyval('highlightstyle',  varargin); if isempty(highlightstyle), highlightstyle = 'opacity'; end
% axis   = keyval('axis',   varargin); if isempty(axis), axis = false; end
% label  = keyval('label',  varargin); % FIXME
% style  = keyval('style',  varargin); % FIXME

% convert the yes/no strings into boolean values
box  = istrue(box);

if isempty(hlim)
  hlim = 'maxmin';
end

if isempty(vlim)
  vlim = 'maxmin';
end

if isempty(clim)
  clim = 'maxmin';
end

if ischar(hlim)
  switch hlim
    case 'maxmin'
      hlim = [min(hdat) max(hdat)];
    case 'absmax'
      hlim = max(abs(hdat));
      hlim = [-hlim hlim];
    otherwise
      error('unsupported option for hlim')
  end % switch
end % if ischar

if ischar(vlim)
  switch vlim
    case 'maxmin'
      vlim = [min(vdat) max(vdat)];
    case 'absmax'
      vlim = max(abs(vdat));
      vlim = [-vlim vlim];
    otherwise
      error('unsupported option for vlim')
  end % switch
end % if ischar

if ischar(clim)
  switch clim
    case 'maxmin'
      clim = [min(cdat(:)) max(cdat(:))];
    case 'absmax'
      clim = max(abs(cdat(:)));
      clim = [-clim clim];
    otherwise
      error('unsupported option for clim')
  end % switch
end % if ischar

if isempty(hpos);
  hpos = (hlim(1)+hlim(2))/2;
end

if isempty(vpos);
  vpos = (vlim(1)+vlim(2))/2;
end

if isempty(width),
  width = hlim(2)-hlim(1);
  width = width * length(hdat)/(length(hdat)-1);
  autowidth = true;
else
  autowidth = false;
end

if isempty(height),
  height = vlim(2)-vlim(1);
  height = height * length(vdat)/(length(vdat)-1);
  autoheight = true;
else
  autoheight = false;
end

% hlim
% vlim

% first shift the horizontal axis to zero
hdat = hdat - (hlim(1)+hlim(2))/2;
% then scale to length 1
hdat = hdat ./ (hlim(2)-hlim(1));
% then scale to compensate for the patch size
hdat = hdat * (length(hdat)-1)/length(hdat);
% then scale to the new width
hdat = hdat .* width;
% then shift to the new horizontal position
hdat = hdat + hpos;

% first shift the vertical axis to zero
vdat = vdat - (vlim(1)+vlim(2))/2;
% then scale to length 1
vdat = vdat ./ (vlim(2)-vlim(1));
% then scale to compensate for the patch size
vdat = vdat * (length(vdat)-1)/length(vdat);
% then scale to the new width
vdat = vdat .* height;
% then shift to the new vertical position
vdat = vdat + vpos;

h = uimagesc(hdat, vdat, cdat, clim);

if ~isempty(highlight)
  switch highlightstyle
    case 'opacity'
      set(h,'AlphaData',highlight);
      set(h, 'AlphaDataMapping', 'scaled');
      alim([0 1]);
    case 'outline'
      % the significant voxels could be outlined with a black contour
      error('unsupported highlightstyle')
    otherwise
      error('unsupported highlightstyle')
  end % switch highlightstyle
end

if box
  boxposition = zeros(1,4);
  % this plots a box around the original hpos/vpos with appropriate width/height
  boxposition(1) = hpos - width/2;
  boxposition(2) = hpos + width/2;
  boxposition(3) = vpos - height/2;
  boxposition(4) = vpos + height/2;
  plot_box(boxposition);
end
