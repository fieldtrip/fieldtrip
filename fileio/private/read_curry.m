function s = read_curry(filename)

% READ_CURRY reads and parses Curry V2 and V4 ascii files and returns the
% content in a structure that is similar to the block-structured layout of
% the file. This function does not interpret the content of the file, but
% is intended as a helper function for READ_CURRY_XXX functions (where XXX
% is the extension of the file).
%
% Use as
%   s = read_curry(filename)

% Copyright (C) 2005, Robert Oostenveld
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


fid = fopen_or_error(filename, 'rt');

s    = [];
line = [];

while ~feof(fid)
  line = fgetl(fid);
  blockname = '';

  % search for the beginning of a block
  if ~isempty(strfind(line, 'START'))
    [blockname, blocktype] = strtok(line);
    blockname = strtrim(blockname);
    blocktype = strtrim(strtok(blocktype));
    val = {};

    % fprintf('reading block %s\n', blockname);
    while isempty(strfind(line, 'END'))

      % read the next line
      line = fgetl(fid);
      % process the data in the line
      if strcmp(blocktype, 'START') && ~isempty(strfind(line, '='))
        % the line looks like "lhs = rhs"
        [lhs, rhs] = strtok(line, '=');
        lhs = strtrim(lhs);
        rhs(1) = ' ';
        rhs = strtrim(rhs);
        val = str2double(rhs);
        if isnan(val)
          val = rhs;
        end
        s = setfield(s, blockname, lhs, val);
      elseif strcmp(blocktype, 'START_LIST')
        % this contains unstructured data, process this line as being part of a list
        val{end+1} = line;
      else
        % process this line differently
        % FIXME
      end

      if isnumeric(line) || feof(fid)
        ft_error('unexpected end of file');
      end

    end % while ~END
  end % if START of block

  % store the collected cell-array if it was part of an unstructured list
  if strcmp(blocktype, 'START_LIST')
    val = val(1:(end-1));   % remove the last item, which is the end of the block
    val = val(:);           % make it a column array
    s = setfield(s, blockname, val);
  end
  
  % clear them to prevent confusion on the next loop
  blockname = '';
  blocktype = '';
  
end  % while ~feof

fclose(fid);
