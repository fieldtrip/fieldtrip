function fieldtrip2spss(filename, labels, data)

% FIELDTRIP2SPSS compiles data and correpsonding labels into a textfile,
% suitable for importing into SPSS or JASP (jasp-stats.org).
%
% Use as
%   fieldtrip2spss(filename, labels, data)
%
% When exporting from MATLAB, set:
%   - filename; should be string (e.g. 'counts.txt')
%   - labels; should be a cell array (e.g. {'ones', 'twos', 'threes'})
%   - data; should be either a vector or matrix (e.g. [1 2 3; 1 2 3; 1 2 3])
%
% When importing to SPSS, set;
%   - variables included at top of file: 'yes'
%   - first case of data on line number: '2' (default)
%   - delimiter appearing between variables: 'tab' (default)
%
% In case the columns that make up the data matrix have unequal lengths
% (e.g. because of different number of subjects per group), use:
%   data         = ones(30,2)*9999
%   data(1:30,1) = 1 (30 subj in Group 1)
%   data(1:20,2) = 2 (20 subj in Group 2)
% After importing to SPSS, click the Missing cell in the Variable View
% window and enter 9999 as the missing value definition.

% Copyright (C) 2011-2017, Arjen Stolk
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


% check whether data and labels have the same lengths
if ~isequal(size(data,2),size(labels,2))
    ft_warning('the data and labels have unequal number of columns');
end

% print labels and append data
fprintf('exporting to %s \n', filename);
txt = sprintf('%s\t',labels{:});
txt(end) = '';
dlmwrite(filename, txt, '');
dlmwrite(filename, data, '-append', 'delimiter', '\t', 'precision', 4);
