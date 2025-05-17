function ft_headlight(~,~)

% FT_HEADLIGHT places a light along the direction of the camera and updates the light
% position as you rotate the scene and the camera view point changes.
%
% Use as
%   surf(peaks);
%   ft_headlight;
%
% See https://stackoverflow.com/questions/30921003/matlab-how-to-make-camera-light-follow-3d-rotation
%
% See also CAMLIGHT, ROTATE3D

% Copyrights (C) 2025, Robert Oostenveld
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

% call the function as it would be called from the ActionPostCallback
h = gcf;
e.Axes = gca;
update(h, e);

r = rotate3d;
r.ActionPostCallback = @update;
r.Enable = 'on';

function update(h, e)
c = findall(e.Axes, 'type', 'light');
if isempty(c)
  % place the first light
  camlight('headlight');
else
  % update the first light
  camlight(c(1), 'headlight');
end