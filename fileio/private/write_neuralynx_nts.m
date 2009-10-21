function write_neuralynx_nts(filename, nts);

% WRITE_NEURALYNX_NTS writes spike timestamps to a NTS file
%
% Use as
%   write_neuralynx_nts(filename, nts)
%
% See also READ_NEURALYNX_NTS

% Copyright (C) 2007, Robert Oostenveld
%
% $Log: write_neuralynx_nts.m,v $
% Revision 1.1  2009/01/14 09:12:16  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.1  2007/03/21 17:12:04  roboos
% renamed NTE in NTS (filenames and function names)
%
% Revision 1.2  2007/03/21 12:53:51  roboos
% updated the documentation
%
% Revision 1.1  2007/03/14 16:08:09  roboos
% new implementation
%

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



