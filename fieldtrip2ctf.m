function fieldtrip2ctf(filename, data, varargin)

% FIELDTRIP2CTF saves a FieldTrip data structures to a corresponding CTF file. The
% file to which the data is exported depends on the input data structure that you
% provide.
%
% Use as
%   fieldtrip2ctf(filename, data, ...)
% where "filename" is a string, "data" is a FieldTrip data structure, and
% additional options can be specified as key-value pairs.
%
% The FieldTrip "montage" structure (see FT_APPLY_MONTAGE and the cfg.montage
% option in FT_PREPROCESSING) can be exported to a CTF "Virtual Channels" file.
%
% At this moment the support for other FieldTrip structures and CTF fileformats is
% still limited, but this function serves as a placeholder for future improvements.
%
% See also FT_VOLUMEWRITE, FT_SOURCEWRITE, FT_WRITE_DATA

% Copyright (C) 2015, Robert Oostenveld
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

type = ft_datatype(data);
switch type
  case 'montage'
    % Virtual channels are weighted linear combinations of real channels collected
    % by the CTF MEG System
    
    % rename it for convenience
    montage = data; clear data;
    
    fid = fopen(filename, 'wt');
    assert(fid>0, 'could not open file for writing');
    fprintf(fid, '// Virtual channel configuration\n');
    for i=1:length(montage.labelnew)
      sel = find(montage.tra(i,:) ~= 0);
      fprintf(fid, '\n');
      fprintf(fid, 'VirtualChannel\n');
      fprintf(fid, '{\n');
      fprintf(fid, 'Name:\t%s\n', montage.labelnew{i});
      fprintf(fid, 'Unit:\n');
      for j=1:numel(sel)
        fprintf(fid, 'Ref:\t%s,%f\n', montage.labelorg{sel(j)}, montage.tra(i,sel(j)));
      end
      fprintf(fid, '}\n');
    end
    
    % case 'event'
    %   if ~ft_filetype(filename, 'ctf_ds')
    %     error('you should specify the directory name of a CTF dataset to which the MarkerFile.mrk will be added');
    %   end
    %   
    % The MarkerFile.mrk file requires exact line spacing, otherwise the software
    % will fail to read the it properly. Two blank lines must be present between
    % sections and three blank lines must be present at the end of the file.
    
    % case 'raw'
    % this should be written with FT_WRITE_DATA

    % case 'timelock'
    % this should be written as a single trial with FT_WRITE_DATA

    % case 'volume'
    % this could be supported for singlesphere and localspheres

    % for the following representations of processed data there is no suitable CTF file format 
    % case 'freq'
    % case 'source'
    % case 'comp'
    % case 'spike'
    % case 'source'
    % case 'dip'

  otherwise
    error('unsuported data structure "%s" for exporting to CTF', type);
end

