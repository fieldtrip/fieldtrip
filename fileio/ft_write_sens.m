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
%   besa_sfp
%   polhemus_pos
%   matlab
%
% See also FT_READ_SENS, FT_DATATYPE_SENS

% Copyright (C) 2017-2019, Robert Oostenveld
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

% ensure that the directory exists
isdir_or_mkdir(fileparts(filename));

switch format
  case 'bioimage_mgrid'
    [p, f, x] = fileparts(filename);
    filename = fullfile(p, [f '.mgrid']);
    write_bioimage_mgrid(filename, sens);
    
  case 'besa_sfp'
    % this always seems to be in mm
    sens = ft_convert_units(sens, 'mm');
    [p, f, x] = fileparts(filename);
    filename = fullfile(p, [f '.sfp']);
    fid = fopen_or_error(filename, 'wt');
    for i=1:numel(sens.label)
      fprintf(fid, '%s\t%g\t%g\t%g\n', sens.label{i}, sens.elecpos(i,1), sens.elecpos(i,2), sens.elecpos(i,3));
    end
    fclose(fid);
    
  case 'matlab'
    [p, f, x] = fileparts(filename);
    filename = fullfile(p, [f '.mat']);
    save(filename, sens);
    
  case 'polhemus_pos'
    % this always seems to be in mm
    sens = ft_convert_units(sens, 'mm');
    % the files that I know only have the electrode number, not any name
    % they also end with nasion, left, right, but those are not contained in the elec structure
    % the correct amount of whitespace might be slightly different
    [p, f, x] = fileparts(filename);
    filename = fullfile(p, [f '.pos']);
    fid = fopen_or_error(filename, 'wt');
    fprintf(fid, ' %d\n', numel(sens.label));
    for i=1:numel(sens.label)
      fprintf(fid, ' %-3d%27.14g%27.14g%27.14g\n', i, sens.elecpos(i,1), sens.elecpos(i,2), sens.elecpos(i,3));
    end
    fclose(fid);
    
  otherwise
    ft_error('unsupported format "%s"', format);
end % switch format

