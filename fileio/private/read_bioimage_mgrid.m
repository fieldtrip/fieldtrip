function [elec] = read_bioimage_mgrid(mgridfile)

% READ_BIOIMAGE_MGRID reads BioImage Suite *.mgrid files and converts them
% into a FieldTrip-compatible elec datatype structure with electrode
% positions in xyz coordinates (equals voxel coordinates in mm)
%
% Use as
%   elec = read_bioimage_mgrid(filename)
% where the filename has the .mgrid file extension
%
% See also FT_READ_SENS, FT_DATATYPE_SENS

% Copyright (C) 2017, Arjen Stolk & Sandon Griffin
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

% define output
elec.label = {}; % N x 1 cell array
elec.elecpos = []; % N x 3 matrix

% conditional variables
WaitForPos = 0;
WaitForDescript = 0;
WaitForDim = 0;

% open and read ascii-file line by line
fileline = 0;
fid = fopen(mgridfile,'r'); % open ascii-file
while fileline >= 0 % read line by line
  
  fileline = fgets(fid); % read a line
  if fileline > 0
    
    % grid number
    if ~isempty(findstr(fileline,'# Electrode Grid '))
      GridNr = sscanf(fileline(findstr(fileline, 'Grid '):end),'Grid %i');
      WaitForDescript = 1;
      WaitForDim = 1;
    end
    
    % grid description
    if ~isempty(findstr(fileline,'#Description')) && WaitForDescript
      nextfileline = fgets(fid);
      GridDescript = sscanf(nextfileline,'%s');
      WaitForDescript = 0;
    end
    
    % grid dimension
    if ~isempty(findstr(fileline,'#Dimensions')) && WaitForDim
      nextfileline = fgets(fid);
      GridDim = sscanf(nextfileline,'%i %i')'; % row & column
      WaitForDim = 0;
    end
    
    % electrode number
    if ~isempty(findstr(fileline,'# Electrode '))
      type = sscanf(fileline(findstr(fileline,'Electrode '):end),'Electrode %s');
      if ~strcmp(type, 'Grid')
        ElecNr = sscanf(fileline(findstr(fileline,'Electrode '):end),'Electrode %i %i')';
        WaitForPos = 1;
      end
    end
    
    % electrode position
    if ~isempty(findstr(fileline,'#Position')) && WaitForPos
      nextfileline = fgets(fid);
      ElecPos = sscanf(nextfileline,'%f %f %f')'; % x, y, and z
      WaitForPos = 0;
    end
    
    % electrode present
    if ~isempty(findstr(fileline,'#Electrode Present'))
      nextfileline = fgets(fid);
      ElecPres = logical(sscanf(nextfileline,'%f'))'; % 1 = present; 0 = not present
      if ~ElecPres
        ElecPos = NaN(1,3);
      end
    end
    
    % store
    if ~isempty(findstr(fileline,'#Value'))
      elec.label{end+1,1} = [GridDescript num2str(GridDim(2) - ElecNr(2) + GridDim(2)*ElecNr(1))];
      elec.elecpos(end+1,:) = ElecPos;
    end
    
  end % if fileline
end % end of while loop
fclose(fid);

% assume a unipolar montage, where the reference is implicit and the electrode positions are identical to the channel positions
elec.chanpos = elec.elecpos;
elec.tra = eye(size(elec.elecpos,1));

% for now it remains unclear what units and coordinate system this is expressed, see https://github.com/fieldtrip/fieldtrip/pull/342
% elec.coordsys = ERROR
% elec.unit = ERROR
