function file = dccnpath(filename)

% DCCNPATH updates the path in the filename of a test file located on DCCN central
% storage. It helps to read test file from Linux, Windows or macOS computers in the
% DCCN.
%
% Use as
%  filename = dccnpath(filename)
%
% The function assumes that on a Windows machine the DCCN home network drive is
% mapped to the H: drive on Windows, or that it is mounted on /Volume/home for macOS.
%
% Similar functionality as this function could be realized with an anonymous function
% like this:
%
% if ispc
%   dccnpath = @(filename) strrep(strrep(filename,'/home','H:'),'/','\');
% else
%   dccnpath = @(filename) strrep(strrep(filename,'H:','/home'),'\','/');
% end

% Copyright (C) 2012-2019, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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

if ispc
  file = strrep(filename,'/home','H:');
  file = strrep(file,'/','\');
else
  file = strrep(filename,'H:','/home');
  file = strrep(file,'\','/');
end

[p, f, x] = fileparts(file);

alternative1 = [f x]; % this should not be used when it is only "test", since that is too generic
alternative2 = strrep(filename, '/home/common/matlab/fieldtrip', '/Volumes/home/common/matlab/fieldtrip');
alternative3 = strrep(filename, '/home/common/matlab/fieldtrip', '/Volumes/128GB/data/test');

% exist(alternative1, 'file') returns 7 in case the present working directory matches alternative1
if endsWith(pwd, alternative1)
  % it is not a valid match
  alternative1 = '';
end

if ~exist(file, 'file') && exist(alternative1, 'file') && ~strcmp(alternative1, 'test')
  warning('using local copy %s instead of %s', alternative1, file);
  file = alternative1;
elseif ~exist(file, 'file') && exist(alternative2, 'file')
  warning('using local copy %s instead of %s', alternative2, file);
  file = alternative2;
elseif ~exist(file, 'file') && exist(alternative3, 'file')
  warning('using local copy %s instead of %s', alternative3, file);
  file = alternative3;
end
