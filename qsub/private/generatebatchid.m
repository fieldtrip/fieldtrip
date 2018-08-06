function id = generatebatchid(batch)

% GENERATEBATCHID generates a unique string that can be used as batch identifier.
%
% See also GENERATEJOBID, GENERATESESSIONID

% Copyright (C) 2011-2012, Robert Oostenveld
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

if nargin~=1
  error('incorrect number of input arguments');
end

id = sprintf('%s_%s_p%d_b%d', getusername(), gethostname(), getpid(), batch);
id = fixname(id);

