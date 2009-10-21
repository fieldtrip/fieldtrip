function write_neuralynx_ncs(filename, ncs);

% WRITE_NEURALYNX_NCS writes continuous data to a NCS file
% The input data should be scaled in uV.
%
% Use as
%   write_neuralynx_ncs(filename, ncs)
%
% See also READ_NEURALYNX_NCS

% Copyright (C) 2005-2007, Robert Oostenveld
%
% $Log: write_neuralynx_ncs.m,v $
% Revision 1.1  2009/01/14 09:12:16  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.7  2007/03/21 12:51:20  roboos
% included the scaling to int16 AD values into this function, i.e. the input data should be uV (double)
%
% Revision 1.6  2007/03/19 16:54:52  roboos
% changed numeric representation of something in the ascii header
%
% Revision 1.5  2007/02/21 09:52:27  roboos
% changed the writing of the data to reflect t he changed representation (2D numeric array instead of 1D cell array)
%
% Revision 1.4  2007/01/09 09:40:03  roboos
% write timestamps as unsigned (uint64 instead of int64)
%
% Revision 1.3  2005/12/02 09:01:42  roboos
% fixed a comment
%
% Revision 1.2  2005/09/09 12:27:02  roboos
% cleaned up the code, changed the looping over records, added the number of records to the output
%
% Revision 1.1  2005/08/05 13:41:39  roboos
% new implementation
%

if ~isa(ncs.TimeStamp, 'uint64')
  error('timestamps should be uint64');
end

% convert the data from uV into V
ncs.dat = ncs.dat * 1e-6;
% scale the data and convert to 16 bits,
% this has to be done prior to writing the header
ADMaxValue = double(intmax('int16'));
ADMaxVolts = max(abs(ncs.dat(:)));
if ADMaxVolts>0
  ADBitVolts = ADMaxVolts / ADMaxValue;
else
  ADBitVolts = 1;
end
ncs.dat = int16(ncs.dat / ADBitVolts);
% update the header with the calibration values
ncs.hdr.ADBitVolts = ADBitVolts;
ncs.hdr.ADMaxValue = ADMaxValue;

% construct the header
buf  = [];
buf  = [buf sprintf('######## Neuralynx Data File Header\r\n')];
buf  = [buf sprintf('## File Name: %s\r\n', filename)];
buf  = [buf sprintf('## Time Opened: (m/d/y): %s  At Time: %s\r\n', datestr(clock, 'mm/dd/yy'), datestr(clock, 'HH:MM:SS'))];
f = fieldnames(ncs.hdr);
for i=1:length(f)
  v = getfield(ncs.hdr, f{i});
  switch class(v)
    case 'char'
      buf = [buf sprintf('-%s\t%s\r\n', f{i}, v)];
    case 'double'
      buf = [buf sprintf('-%s\t%s\r\n', f{i}, num2str(v))];
    otherwise
      error('unknown class in writing header');
  end
end

% pad the rest of the header with zeros
buf((end+1):16384) = 0;

% open the file and write the header
fid  = fopen(filename, 'wb', 'ieee-le');
fwrite(fid, buf);

% The format of a continuous sampled record is
%   int64 TimeStamp
%   int32 ChanNumber
%   int32 SampFreq
%   int32 NumValidSamp
%   int16 Samp[0] ... int16 Samp[511]
% Note that if NumValidSamp < 512, Samp[n], where n >= NumValidSamp, will
% contain random data.

for i=1:size(ncs.dat,2)
  % write a single continuous data record
  fwrite(fid, ncs.TimeStamp(i)   , 'uint64');
  fwrite(fid, ncs.ChanNumber(i)  , 'int32');
  fwrite(fid, ncs.SampFreq(i)    , 'int32');
  fwrite(fid, ncs.NumValidSamp(i), 'int32');
  fwrite(fid, ncs.dat(:,i)       , 'int16');
end

% close the file
fclose(fid);
