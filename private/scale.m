function [H] = scale(S);

% SCALE returns the homogenous coordinate transformation matrix
% corresponding to a scaling along the x, y and z-axis
% 
% Use as
%   [H] = translate(S)
% where
%   S		[sx, sy, sz] scaling along each of the axes
%   H 	corresponding homogenous transformation matrix

% Copyright (C) 2000-2005, Robert Oostenveld
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

H = [
  S(1) 0    0    0 
  0    S(2) 0    0
  0    0    S(3) 0
  0    0    0    1
  ];

  
