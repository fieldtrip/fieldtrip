function fieldtrip2spss(filename, labels, data)

% FIELDTRIP2SPSS compiles a dataset and its labels into a textfile, 
% suitable for import to SPSS. 
%
% When exporting from MATLAB, use:
%   - filename; should be string (e.g. 'counts.txt')
%   - labels; should be a cell array (e.g. {'ones','twos','threes'})
%   - data; should be either a vector or matrix (e.g. [1 2 3; 1 2 3; 1 2 3])
%
% When importing to SPSS, use;
%   - variables included at top of file: 'yes'
%   - first case of data on line number: '2' (default)
%   - delimiter appearing between variables: 'tab' (default)
%
% In case the columns that make up the data matrix have unequal lengths
% (e.g. because of different number of subjects per group), use:
%   data         = ones(30,2)*9999
%   data(:,1)    = Group1 (30 subj) 
%   data(1:20,2) = Group2 (only 20 subj)
% After importing to SPSS, click the Missing cell in the Variable View
% window and enter 9999 as the missing value definition.
%
% See also NUTMEG2FIELDTRIP, SPASS2FIELDTRIP, LORETA2FIELDTRIP

% Copyright (C) 2011-2013, Arjen Stolk
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble provenance

% check whether data and labels have the same lengths
if ~isequal(size(data,2),size(labels,2))
    warning('the data and labels have unequal number of columns');
end

% print labels and append data
fprintf('exporting to %s \n', filename);
txt = sprintf('%s\t',labels{:});
txt(end) = '';
dlmwrite(filename, txt, '');
dlmwrite(filename, data, '-append', 'delimiter', '\t', 'precision', 4);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble provenance