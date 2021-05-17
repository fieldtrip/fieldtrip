function tsv = ft_read_tsv(filename)

% FT_READ_TSV reads information from a tab-separated-values file and represents it as a MATLAB table
%
% Use as
%   tsv = ft_read_tsv(filename)
%
% See also FT_WRITE_TSV, FT_READ_JSON, FT_WRITE_JSON, READTABLE, WRITETABLE

% Copyright (C) 2018-2021, Robert Oostenveld
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

ft_info('reading ''%s''\n', filename);
tsv = readtable(filename, 'Delimiter', 'tab', 'FileType', 'text', 'TreatAsEmpty', 'n/a', 'ReadVariableNames', true);
