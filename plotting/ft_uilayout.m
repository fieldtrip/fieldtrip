function ft_uilayout(h, varargin)

% FT_UILAYOUT is a helper function to make a consistent graphical user interafce with
% multiple control elements. This function will find all elements with a specific tag
% and style, and update or position them consistently.
%
% Use as
%   ft_uilayout(h, 'tag', '...', 'style', '...', ...)
% where h is the figure handle and 'tag' and 'style' are used to specifying which
% user control elements in the figure should be selected.
%
% You can pass most options from UICONTROL as key-value pair, such as
% 'BackgroundColor', 'CallBack', 'Clipping', 'Enable', 'FontAngle', 'FontName',
% 'FontSize', 'FontUnits', 'FontWeight', 'ForegroundColor', 'HorizontalAlignment',
% 'Max', 'Min', 'Position', 'Selected', 'String', 'Units', 'Value', 'Visible'.
%
% In addition to the options from UICONTROL, you can use the following key-value
% pairs for a consistent placement of multiple GUI elements relative to each other:
%   'hpos'         = 'auto'       puts elements in horizontal adjacent order with a fixed distance of 0.01
%                    'align'      adjusts the horizontal position of all elements to the first element
%                    'distribute' puts elements in horizontal adjacent order such that they distribute evenly
%                    scalar       sets the horizontal position of elements to the specified scalar
%   'vpos'         = 'auto'       puts elements in vertical adjacent order with a fixed distance of 0.01
%                    'align'      adjusts the vertical position of all elements to the first element
%                    'distribute' puts elements in vertical adjacent order such that they distribute evenly
%                    scalar       sets the vertical position of elements to the specified scalar
%   'width'        = scalar       sets the width of elements to the specified scalar
%   'height'       = scalar       sets the height of elements to the specified scalar
%   'halign'       = 'left'       aligns the horizontal position of elements to the left
%                    'right'      aligns the horizontal position of elements to the right
%   'valign'       = 'top'        aligns the vertical position of elements to the top
%                    'bottom'     aligns the vertical position of elements to the bottom
%   'halign'       = 'left'       aligns the horizontal position of elements to the left
%                    'right'      aligns the horizontal position of elements to the right
%   'hshift'       = scalar       shift the elements in horizontal direction
%   'vshift'       = scalar       shift the elements in vertical direction
%
% Here is an example that positions a number of buttons in a 2x3 grid. It makes use
% of regular expressions to match the tags to the rows and columns.
%
%   h = figure;
%   uicontrol('style', 'pushbutton', 'string', '11', 'tag', 'row1_column1');
%   uicontrol('style', 'pushbutton', 'string', '12', 'tag', 'row1_column2');
%   uicontrol('style', 'pushbutton', 'string', '13', 'tag', 'row1_column3');
%   uicontrol('style', 'pushbutton', 'string', '21', 'tag', 'row2_column1');
%   uicontrol('style', 'pushbutton', 'string', '22', 'tag', 'row2_column2');
%   uicontrol('style', 'pushbutton', 'string', '23', 'tag', 'row2_column3');
%
%   ft_uilayout(h, 'tag', '^row1', 'vpos', 100);
%   ft_uilayout(h, 'tag', '^row2', 'vpos', 200);
%
%   ft_uilayout(h, 'tag', 'column1$', 'hpos', 100);
%   ft_uilayout(h, 'tag', 'column2$', 'hpos', 200);
%   ft_uilayout(h, 'tag', 'column3$', 'hpos', 300);
%
%   ft_uilayout(h, 'tag', '.*', 'BackGroundColor', [1 0 0]);
%
% See also UICONTROL

% Copyright (C) 2009-2015, Robert Oostenveld
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

% these are used to select the uicontrol elements that will be updated
tag    = ft_getopt(varargin, 'tag');
style  = ft_getopt(varargin, 'style');

if any(tag=='.' | tag=='*' | tag=='?' | tag=='^' | tag=='$')
  % use regular expressions
  tagopt = {'-regexp'};
else
  tagopt = {};
end

% determine all children
if ~isempty(tag) && ~isempty(style)
  c = findobj(h, tagopt{:}, 'tag', tag, 'style', style);
elseif ~isempty(tag)
  c = findobj(h, tagopt{:}, 'tag', tag);
elseif ~isempty(style)
  c = findobj(h, 'style', style);
else
  c = findobj(h);
end
c = flipud(c);
% fprintf('selected %d elements\n', numel(c));

% these are the normal features of an uicontrol
BackgroundColor         = ft_getopt(varargin, 'BackgroundColor');
CallBack                = ft_getopt(varargin, 'CallBack');
Clipping                = ft_getopt(varargin, 'Clipping');
Enable                  = ft_getopt(varargin, 'Enable');
FontAngle               = ft_getopt(varargin, 'FontAngle');
FontName                = ft_getopt(varargin, 'FontName');
FontSize                = ft_getopt(varargin, 'FontSize');
FontUnits               = ft_getopt(varargin, 'FontUnits');
FontWeight              = ft_getopt(varargin, 'FontWeight');
ForegroundColor         = ft_getopt(varargin, 'ForegroundColor');
HorizontalAlignment     = ft_getopt(varargin, 'HorizontalAlignment');
Max                     = ft_getopt(varargin, 'Max');
Min                     = ft_getopt(varargin, 'Min');
Position                = ft_getopt(varargin, 'Position');
Selected                = ft_getopt(varargin, 'Selected');
String                  = ft_getopt(varargin, 'String');
Units                   = ft_getopt(varargin, 'Units');
Value                   = ft_getopt(varargin, 'Value');
Visible                 = ft_getopt(varargin, 'Visible');
Tag                     = ft_getopt(varargin, 'retag');   % this is to change the tag on the selected items
Style                   = ft_getopt(varargin, 'restyle'); % this is to change the style on the selected items

feature = {
  'BackgroundColor'
  'CallBack'
  'Clipping'
  'Enable'
  'FontAngle'
  'FontName'
  'FontSize'
  'FontUnits'
  'FontWeight'
  'ForegroundColor'
  'HorizontalAlignment'
  'Max'
  'Min'
  'Position'
  'Selected'
  'String'
  'Units'
  'Value'
  'Visible'
  'Tag'   % this is used for retag
  'Style' % this is used for restyle
  };

for i=1:length(feature)
  val = eval(feature{i});
  if ~isempty(val)
    % fprintf('setting %s\n', feature{i});
    for j=1:length(c)
      set(c(j), feature{i}, val)
    end
  end
end

% these are special features to help with the positioning of the elements
hpos   = ft_getopt(varargin, 'hpos');
vpos   = ft_getopt(varargin, 'vpos');
halign = ft_getopt(varargin, 'halign', 'left');
valign = ft_getopt(varargin, 'valign', 'bottom');
width  = ft_getopt(varargin, 'width');
height = ft_getopt(varargin, 'height');
hshift = ft_getopt(varargin, 'hshift');
vshift = ft_getopt(varargin, 'vshift');

if isempty(hpos) && isempty(vpos) && isempty(width) && isempty(height) && isempty(hshift) && isempty(vshift)
  % re-positioning of the elements is not needed
  return
end

% this will contain the horizontal position, vertical position, with and height
pos = zeros(length(c), 4);

for i=1:length(c)
  % this only works if the units are normalized
  pos(i,:) = get(c(i), 'position');
end

if ~isempty(width)
  pos(:,3) = width;
end
width = pos(:,3);

if ~isempty(height)
  pos(:,4) = height;
end
height = pos(:,4);

if ~isempty(hpos)
  if isequal(hpos, 'auto')
    scale = (1 - 0.01 - 0.01*length(c)) / sum(width);
    if scale>0 && scale<1
      % fprintf('adjusting width with %f\n', scale);
      width = width*scale;
      pos(:,3) = width;
    end
    hpos = cumsum([0.01; width+0.01]);
    if isequal(halign, 'right')
      hpos = 1-hpos(end) + hpos - 0.01;
    end
    hpos = hpos(1:end-1);
  elseif isequal(hpos, 'align')
    if isequal(halign, 'right')
      hpos = pos(end,1); % the position of the last element
    else % default behavior
      hpos = pos(1,1); % the position of the first element
    end
  elseif isequal(hpos, 'distribute')
    minpos = min(pos(:,1));
    maxpos = max(pos(:,1));
    hpos = linspace(minpos, maxpos, length(c));
  end
  pos(:,1) = hpos;
  pos(:,3) = width;
end % hpos

if ~isempty(vpos)
  if isequal(vpos, 'auto')
    scale = (1 - 0.01 - 0.01*length(c)) / sum(height);
    if scale>0 && scale<1
      % fprintf('adjusting height with %f\n', scale);
      height = height*scale;
      pos(:,4) = height;
    end
    vpos = cumsum([0.01; height+0.01]);
    if isequal(valign, 'bottom') % default
      vpos = 1-vpos(end) + vpos - 0.01;
    end
    vpos = vpos(end-1:-1:1);
  elseif isequal(vpos, 'align')
    if isequal('valign', 'bottom')
      vpos = pos(end,2); % the position of the last element
    else % default behavior
      vpos = pos(1,2); % the position of the first element
    end
  elseif isequal(vpos, 'distribute')
    minpos = min(pos(:,2));
    maxpos = max(pos(:,2));
    vpos = linspace(minpos, maxpos, length(c));
  end
  pos(:,2) = vpos;
end % vpos

if ~isempty(hshift)
  % apply horizontal shift
  pos(:,1) = pos(:,1) + hshift;
end

if ~isempty(vshift)
  % apply vertical shift
  pos(:,2) = pos(:,2) + vshift;
end

% assign the new/automatic position to each of the elements
for i=1:length(c)
  set(c(i), 'position', pos(i,:));
end
