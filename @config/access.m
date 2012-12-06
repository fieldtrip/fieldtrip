function y = access(x, cmd);

% ACCESS Return the number of accesses (assignments and references) to a CONFIGURATION object.

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

if nargin==1
  fprintf('----------- values     -----------\n');
  disp(x.value);

  fprintf('----------- original   -----------\n');
  disp(x.original);

  fprintf('----------- assignment -----------\n');
  disp(x.assign);

  fprintf('----------- reference -----------\n');
  disp(x.reference);
else
  switch cmd(1)
    case 'v'
      y = x.value;
    case 'r'
      y = x.reference;
    case 'a'
      y = x.assign;
    case 'o'
      y = x.original;
    otherwise
      error('Incorrect command for accessing a config object.')
  end
end
