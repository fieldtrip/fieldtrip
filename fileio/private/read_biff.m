function [this] = read_biff(filename, opt)

% READ_BIFF reads data and header information from a BIFF file
%
% This is a attemt for a reference implementation to read the BIFF
% file format as defined by the Clinical Neurophysiology department of
% the University Medical Centre, Nijmegen.
%
% read all data and information
%   [data]  = read_biff(filename)
% or read a selected top-level chunk
%   [chunk] = read_biff(filename, chunkID)
%
% known top-level chunk id's are
%   data    : measured data         (matrix)
%   dati    : information on data       (struct)
%   expi    : information on experiment (struct)
%   pati    : information on patient    (struct)
%   evnt    : event markers         (struct)

% Copyright (C) 2000, Robert Oostenveld
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

define_biff;
this = [];

fid = fopen_or_error(filename, 'r');
fseek(fid,0,'bof');

[id, siz] = chunk_header(fid);

switch id
  case 'SEMG'
    child  = subtree(BIFF, id);
    this   = read_biff_chunk(fid, id, siz, child);

  case 'LIST'
    fprintf('skipping unimplemented chunk id="%s" size=%4d\n', id, siz);

  case 'CAT '
    fprintf('skipping unimplemented chunk id="%s" size=%4d\n', id, siz);

  otherwise
    fprintf('skipping unrecognized chunk id="%s" size=%4d\n', id, siz);
    fseek(fid, siz, 'cof');

end                         % switch
fclose(fid);                        % close file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION read_biff_chunk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function this = read_biff_chunk(fid, id, siz, chunk)

% start with empty structure
this = [];

if strcmp(id, 'null')               % this is an empty chunk
  fprintf('skipping empty chunk id="%s" size=%4d\n', id, siz);
  assert(~feof(fid));
  fseek(fid, siz, 'cof');

elseif isempty(chunk)               % this is an unrecognized chunk
  fprintf('skipping unrecognized chunk id="%s" size=%4d\n', id, siz);
  assert(~feof(fid));
  fseek(fid, siz, 'cof');

else
  eoc   = ftell(fid) + siz;
  name  = char(chunk.desc(2));
  type  = char(chunk.desc(3));

  fprintf('reading chunk id= "%s" size=%4d name="%s"\n', id, siz, name);

  switch type
    case 'group'

      while ~feof(fid) && ftell(fid)<eoc
        % read all subchunks
        [id, siz] = chunk_header(fid);
        child     = subtree(chunk, id);

        if ~isempty(child)
          % read data and add subchunk data to chunk structure
          name  = char(child.desc(2));
          val   = read_biff_chunk(fid, id, siz, child);
          this  = setfield(this, name, val);
        else
          fprintf('skipping unrecognized chunk id="%s" size=%4d\n', id, siz);
          fseek(fid, siz, 'cof');
        end
      end                     % while

    case 'string'
      this = char(fread(fid, siz, 'uchar')');

    case {'char', 'uchar', 'int8', 'int16', 'int32', 'int64', 'uint8', 'uint16', 'uint32', 'float32', 'float64'}
      this = fread(fid, 1, type);

    case {'int8vec', 'int16vec', 'int32vec', 'int64vec', 'uint8vec', 'uint16vec', 'uint32vec', 'float32vec', 'float64vec'}
      ncol = fread(fid, 1, 'uint32');
      this = fread(fid, ncol, type(1:(length(type)-3)));

    case {'int8mat', 'int16mat', 'int32mat', 'int64mat', 'uint8mat', 'uint16mat', 'uint32mat', 'float32mat', 'float64mat'}
      nrow = fread(fid, 1, 'uint32');
      ncol = fread(fid, 1, 'uint32');
      this = fread(fid, [nrow, ncol], type(1:(length(type)-3)));

    otherwise
      fseek(fid, siz, 'cof');         % skip this chunk
      sprintf('unimplemented data type "%s" in chunk "%s"', type, id);
      % ft_warning(ans);
  end                     % switch chunk type
end                     % else

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION subtree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function child = subtree(parent, id)
blank = strfind(id, ' ');
while ~isempty(blank)
  id(blank) = '_';
  blank = strfind(id, ' ');
end
elem  = fieldnames(parent);                   % list of all subitems
num   = find(strcmp(elem, id));               % number in parent tree
if numel(num) == 1
  child = getfield(parent, char(elem(num)));  % child subtree
else
  child = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION chunk_header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [id, siz] = chunk_header(fid)
id  = char(fread(fid, 4, 'uchar')');          % read chunk ID
siz = fread(fid, 1, 'uint32');            % read chunk size
if strcmp(id, 'GRP ') || strcmp(id, 'BIFF')
  id  = char(fread(fid, 4, 'uchar')');        % read real chunk ID
  siz = siz - 4;                  % reduce size by 4
end
