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
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information


fid = fopen(filename, 'rt');
if fid<0
  error(sprintf('could not open file: %s', filename));
end

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
      if strcmp(blocktype, 'START') & ~isempty(strfind(line, '='))
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

      if isnumeric(line) | feof(fid)
        error('unexpected end of file');
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
