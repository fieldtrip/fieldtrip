function nirs = fieldtrip2homer(data)

% FIELDTRIP2HOMER converts a continuous raw data structure from FieldTrip format to
% Homer format.
%
% Use as
%   nirs = fieldtrip2homer(data)
% where the input data structure is formatted according to the output of
% FT_PREPROCESSING and the output nirs structure is according to Homer.
%
% See https://www.nitrc.org/plugins/mwiki/index.php/homer2:Homer_Input_Files#NIRS_data_file_format
% for a description of the Homer data structure.
%
% See also HOMER2FIELDTRIP, FT_PREPROCESSING, FT_DATATYPE_RAW

% Copyright (C) 2020, Robert Oostenveld
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

% ensure that the input is according to the required format
data = ft_checkdata(data, 'datatype', 'raw', 'feedback', true);

hdr = ft_fetch_header(data);
dat = ft_fetch_data(data);

nirs.t = data.time{1};
nirs.d = dat(strcmp(hdr.chantype, 'nirs'), :)';
nirs.s = dat(strcmp(hdr.chantype, 'stimulus'), :)';
nirs.aux = dat(strcmp(hdr.chantype, 'aux'), :)';
nirs.CondNames = hdr.label(strcmp(hdr.chantype, 'stimulus'));
nirs.SD = opto2homer(hdr.opto);
