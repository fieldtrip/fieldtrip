function [hdr] = neuralynx_getheader(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for reading the 16384 byte header from any Neuralynx file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid     = fopen(filename, 'rb', 'ieee-le');
buf     = fread(fid, 16*1024, 'uint8=>char');
fclose(fid);

buf     = buf(:)';
nl      = find(buf==10);    % determine the new-lines
cr      = find(buf==13);    % determine the carriage-returns
begchar = [1 nl(1:(end-1))];
endchar = nl - 1;
num     = length(nl);

hdr        = [];
hdr.Header = buf; % remember the full header in its original format

for i=1:num
  line = fliplr(deblank(fliplr(deblank(char(buf(begchar(i):endchar(i)))))));
  if numel(line)==0
    % line is empty
    continue
  elseif line(1)=='#'
    % line contains a comment
    continue
  else
    % strip the '-' sign
    while line(1)=='-'
      line = line(2:end);
    end
    % replace tabs with spaces
    line(find(line==9)) = ' ';
    % cut into pieces
    item = strread(line, '%s');
    if length(item)==2
      key = item{1};
      val = item{2};
      if any(val(1)=='-01234567989')
        % try to convert to number
        val = str2num(val);
        if isempty(val)
          % revert to the original text
          val = item{2};
        end
      end
      % remove unuseable characters from the variable name (key)
      key = key(key~=':');
      % assign the value to the header structure
      hdr = setfield(hdr, key, val);
    end
  end
end

