function y = access(x, cmd, key, val)

% ACCESS Return the number of accesses (assignments and references) to a CONFIGURATION object.
%

% Copyright (C) 2012-2015, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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

if nargin==1
  fprintf('----------- values     -----------\n');
  disp(x.value);

  fprintf('----------- original   -----------\n');
  disp(x.original);

  fprintf('----------- assignment -----------\n');
  disp(x.assign);

  fprintf('----------- reference -----------\n');
  disp(x.reference);

  fprintf('----------- hidden -----------\n');
  disp(x.hidden);
  
else
  switch cmd(1)
    case 'v' % values
      y = x.value;
    case 'o' % original
      y = x.original;
    case 'a' % assignment
      y = x.assign;
    case 'r' % reference
      y = x.reference;
    case 'g' % get hidden
      y = x.hidden.(key);
    case 's' % set hidden
      x.hidden.(key) = val;
      y = x;
    otherwise
      error('Incorrect command for accessing a config object.')
  end
end
