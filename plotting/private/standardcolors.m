function y = standardcolors(x)

% STANDARDCOLORS looks up the RGB value for a named color that is specified as
% a string, or looks up the name given the RGB value.
%
% Use as
%   rgb = standardcolors(name)
% or
%   name = standardcolors(rgb)
% or
%   list = standardcolors
%
% This returns a predefined color as [red green blue] values, according to
% the following mapping:
%   red               = [255   0   0]/255;
%   green             = [  0 192   0]/255;
%   blue              = [  0   0 255]/255;
%   magenta           = [255 255   0]/255;
%   cyan              = [  0 255 255]/255;
%   yellow            = [255 255   0]/255;
%   white             = [255 255 255]/255;
%   black             = [  0   0   0]/255;
%   brain             = [202 100 100]/255;
%   skull             = [140  85  85]/255
%   cortex            = [255 213 119]/255;
%   cortex_light      = [199 194 169]/255;
%   cortex_dark       = [100  97  85]/255;
%   skin              = [249 223 192]/255;
%   skin_light        = [249 223 192]/255;
%   skin_medium_light = [225 194 158]/255;
%   skin_medium       = [188 142 106]/255;
%   skin_medium_dark  = [155 102	65]/255;
%   skin_dark         = [ 91  71  61]/255;
%
% The different skin-based colors follow the Fitzpatrick scale with type I and II
% combined, and return RGB values that approximate those used by Apple in the emoji
% skin tones. See also https://emojipedia.org/emoji-modifier-sequence/
%
% If no specific skin tone is specified, this function returns a light skin color.
% This corresponds with that of one of the developers who approximated his own skin
% color more than 15 years ago upon the first implementation of this function.
%
% See also HTMLCOLORS, COLORSPEC2RGB, FT_COLORMAP, COLORMAP, COLORMAPEDITOR, BREWERMAP, MATPLOTLIB, CMOCEAM

% Copyright (C) 2009-2025, Robert Oostenveld
% Copyright (C) 2025, Jan-Mathijs Schoffelen
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

mapping = {
  'red',               [255   0   0]
  'green',             [  0 192   0]
  'blue',              [  0   0 255]
  'magenta',           [255 255   0]
  'cyan',              [  0 255 255]
  'yellow',            [255 255   0]
  'white',             [255 255 255]
  'black',             [  0   0   0]
  'brain'              [202 100 100]
  'skull',             [140  85  85]
  'cortex',            [255 213 119]
  'cortex_light',      [199 194 169]
  'cortex_dark',       [100  97  85]
  'skin',              [249 223 192]
  'skin_light',        [249 223 192]
  'skin_medium_light', [225 194 158]
  'skin_medium',       [188 142 106]
  'skin_medium_dark',  [155 102	 65]
  'skin_dark',         [ 91  71  61]
  };

if nargin==0
  % return the list with names
  y = mapping(:,1);

elseif ischar(x)
  if isequal(x, 'skin')
    msgId = 'FieldTrip:plotting:private:standardcolors';
    ft_notice('once', msgId);
    ft_notice(msgId, 'The color ''skin'' by default results in a light skin, you can also explicitly specify ''skin_light'','' skin_medium_light'', ''skin_medium'', ''skin_medium_dark'', or ''skin_dark''');
  end
  % look up the corresponding RGB values
  sel = strcmp(mapping(:,1), x);
  y = mapping{sel, 2};
  y = y/255;

elseif isnumeric(x)
  x = x*255;
  % look up the corresponding name
  rgb = vertcat(mapping{:,2});
  rgb(:,1) = rgb(:,1) - x(1);
  rgb(:,2) = rgb(:,2) - x(2);
  rgb(:,3) = rgb(:,3) - x(3);
  [dum, sel] = min(sum(rgb.^2, 2));
  y = mapping{sel,1};
end
