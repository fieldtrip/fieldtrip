function write_neuralynx_nts(filename, nts);

% WRITE_NEURALYNX_NTS writes spike timestamps to a NTS file
%
% Use as
%   write_neuralynx_nts(filename, nts)
%
% See also READ_NEURALYNX_NTS

% Copyright (C) 2007, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

if ~isa(nts.TimeStamp, 'uint64')
  error('timestamps should be uint64');
end

fid  = fopen(filename, 'wb', 'ieee-le');

% construct the header
buf  = [];
buf  = [buf sprintf('######## Neuralynx Data File Header\r\n')];
buf  = [buf sprintf('## File Name: %s\r\n', filename)];
buf  = [buf sprintf('## Time Opened: (m/d/y): %s  At Time: %s\r\n', datestr(clock, 'mm/dd/yy'), datestr(clock, 'HH:MM:SS'))];
f = fieldnames(nts.hdr);
for i=1:length(f)
  v = getfield(nts.hdr, f{i});
  switch class(v)
    case 'char'
      buf = [buf sprintf('-%s\t%s\r\n', f{i}, v)];
    case 'double'
      buf = [buf sprintf('-%s\t%g\r\n', f{i}, v)];
    otherwise
      error('unknown class in writing header');
  end
end

% pad the rest of the header with zeros and write it to file
buf((end+1):16384) = 0;
fwrite(fid, buf);

% The format of a clustered electrode record is
%   int64 TimeStamp

fwrite(fid, nts.TimeStamp, 'uint64');

% close the file
fclose(fid);



