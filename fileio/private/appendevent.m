function c = appendevent(a, b)

% APPENDEVENT

% Copyright (C) 2008, Robert Oostenveld
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

if isempty(a)
  c = b(:);
elseif isempty(b)
  c = a(:);
else
  c = a(:);
  for i=1:numel(b)
    c(end+1).type     = b(i).type;
    c(end  ).value    = b(i).value;
    c(end  ).sample   = b(i).sample;
    if isfield(b, 'timestamp')
      c(end  ).timestamp = b(i).timestamp; % optional
    end
    if isfield(b, 'offset')
      c(end  ).offset    = b(i).offset;    % optional
    end
    if isfield(b, 'duration')
      c(end  ).duration  = b(i).duration;  % optional
    end
  end
end
