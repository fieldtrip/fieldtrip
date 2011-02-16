function h = uimagesc(varargin)
%UIMAGESC  Display scaled image with uneven axis.
%   UIMAGESC(...) is the same as UIMAGE(...) except the data is scaled
%   to use the full colormap. See UIMAGE for details.
%
%   Note: UIMAGESC is based on Matlab's original IMAGESC, Revision 5.11.4.5.
%   UIMAGESC simply calls UIMAGE with a scaled colormap.
% 
%   F. Moisy - adapted from TMW
%   Revision: 1.01,  Date: 2006/06/13.
%
%   See also IMAGE, IMAGESC, UIMAGE.
%
% This function is downloaded on Oct 24th 2008 from www.mathworks.com/matlabcentral/fileexchange/11368

% History:
%   2006/06/12: v1.00, first version.
%
% $Id$

clim = [];
switch (nargin),
  case 0,
    hh = uimage('CDataMapping','scaled');
  case 1,
    hh = uimage(varargin{1},'CDataMapping','scaled');
  case 3,
    hh = uimage(varargin{:},'CDataMapping','scaled');
  otherwise,

    % Determine if last input is clim
    if isequal(size(varargin{end}),[1 2])
      str = false(length(varargin),1);
      for n=1:length(varargin)
        str(n) = ischar(varargin{n});
      end
      str = find(str);
      if isempty(str) || (rem(length(varargin)-min(str),2)==0),
        clim = varargin{end};
        varargin(end) = []; % Remove last cell
      else
        clim = [];
      end
    else
      clim = [];
    end
    hh = uimage(varargin{:},'CDataMapping','scaled');
end

% Get the parent Axes of the image
cax = ancestor(hh,'axes');

if ~isempty(clim),
  set(cax,'CLim',clim)
elseif ~ishold(cax),
  set(cax,'CLimMode','auto')
end

if nargout > 0
    h = hh;
end
