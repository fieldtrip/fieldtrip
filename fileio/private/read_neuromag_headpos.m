function [data] = read_neuromag_headpos(filename)

% READ_NEUROMAG_HEADPOS reads head position information from file. The file
% contains information about Time, Quaternions (q1-q6), goodness of
% fit (g-value) and error.
% Time       q1       q2       q3       q4       q5       q6       g-value  error    
%
% data = read_neuromag_headpos(filename)
%
% where the returned structure data has the fields
%   data.data      Contains the numeric values
%   data.textdata  Contains the Column name
%   data.coldata   Contains the Column name

% Copyright (C) 2017, Simon Homoelle
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

data = importdata(filename);
data.data = data.data';
data.colheaders{2} = 'QUAT001';
data.colheaders{3} = 'QUAT002';
data.colheaders{4} = 'QUAT003';
data.colheaders{5} = 'QUAT004';
data.colheaders{6} = 'QUAT005';
data.colheaders{7} = 'QUAT006';

data.textdata{2} = 'QUAT001';
data.textdata{3} = 'QUAT002';
data.textdata{4} = 'QUAT003';
data.textdata{5} = 'QUAT004';
data.textdata{6} = 'QUAT005';
data.textdata{7} = 'QUAT006';


