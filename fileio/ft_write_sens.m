function ft_write_sens(filename, sens, varargin)

% FT_WRITE_SENS writes electrode information to an external file for further processing in external software.
%
% Use as
%  ft_write_sens(filename, sens, ...)
%
% The specified filename can already contain the filename extention,
% but that is not required since it will be added automatically.
%
% Additional options should be specified in key-value pairs and can be
%   'format'     string, see below
%
% The supported file formats are
%   bioimage_mgrid
%
% See also FT_READ_SENS, FT_DATATYPE_SENS

% Copyright (C) 2017, Robert Oostenveld
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

% get the options
format = ft_getopt(varargin, 'format');

% ensure the input is according to the latest standards
sens = ft_datatype_sens(sens);

switch format
  case 'bioimage_mgrid'
    [p, f, x] = fileparts(filename);
    filename = fullfile(p, [f '.mgrid']);
    write_bioimage_mgrid(filename, sens);
    
  otherwise
    ft_error('unsupported format "%s"', format);
end % switch format

