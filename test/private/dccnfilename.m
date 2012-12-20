function file = dccnfilename(filename)

% DCCNFILENAME applies correct OS path formatting to any given dccn home
% path. 
% The function assumes that on a Windows machine the dccn home
% network drive is mapped to H:. The function basically converts between
% H:\ and /home and uses the appropriate file seperator depending on the
% OS.
%
% Copyright (C) 2012, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

if IsWindows
  file = strrep(filename,'/home','H:');
  file = strrep(file,'/','\');
else
  file = strrep(filename,'H:','/home');
  file = strrep(file,'\','/');
end;

end