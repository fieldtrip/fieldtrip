function [output] = read_mff_bin(filename, begblock, endblock, chanindx)

% READ_MFF_BIN
%
% Use as
%   [hdr] = read_mff_bin(filename)
% or
%   [dat] = read_mff_bin(filename, begblock, endblock);

fid = fopen(filename,'r');

if fid == -1
  error('wrong filename') % could not find signal(n)
end

needhdr = (nargin==1);
needdat = (nargin==4);


if needhdr
  hdr = read_mff_block(fid, [], [], 'skip');
  prevhdr = hdr(end);
  for i=2:hdr.opthdr.nblocks
    hdr(i) = read_mff_block(fid, prevhdr, [], 'skip');
    prevhdr = hdr(end);
  end
  % assign the output variable
  output = hdr;
  
elseif needdat
  prevhdr = [];
  dat     = {};
  block   = 0;
  while true
    block = block+1;
    
    if block<begblock
      [hdr] = read_mff_block(fid, prevhdr, chanindx, 'skip');
      prevhdr = hdr;
      continue % with the next block
    end
    
    if block==begblock
      [hdr, dat] = read_mff_block(fid, prevhdr, chanindx, 'read');
      prevhdr = hdr;
      continue % with the next block
    end
    
    if block<=endblock
      [hdr, dat(:,end+1)] = read_mff_block(fid, prevhdr, chanindx, 'read');
      prevhdr = hdr;
      continue % with the next block
    end
    
    break % out of the while loop
    
  end % while
  
  % assign the output variable
  output = dat;
  
end % need header or data

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [hdr, dat] = read_mff_block(fid, prevhdr, chanindx, action)

if nargin<3
  prevhdr = [];
end

if nargin<4
  action = 'none';
end

%-Endianness
endian = 'ieee-le';

% General information
hdr.version = fread(fid, 1, 'int32', endian);

if hdr.version==0
  % the header did not change compared to the previous one
  % no additional information is present in the file
  % the file continues with the actual data
  hdr = prevhdr;
  hdr.version = 0;
else
  hdr.headersize = fread(fid, 1, 'int32', endian);
  hdr.datasize   = fread(fid, 1, 'int32', endian); % the documentation specified that this includes the data, but that seems not to be the case
  hdr.nsignals   = fread(fid, 1, 'int32', endian);
  
  % channel-specific information
  hdr.offset = fread(fid, hdr.nsignals, 'int32', endian);
  
  % signal depth and frequency for each channel
  for i = 1:hdr.nsignals
    hdr.depth(i)   = fread(fid, 1, 'int8', endian);
    hdr.fsample(i) = fread(fid, 1, 'bit24', endian); %ingnie: is bit24 the same as int24?
  end
  
  %-Optional header length
  hdr.optlength   = fread(fid, 1, 'int32', endian);
  
  if hdr.optlength
    hdr.opthdr.EGItype  = fread(fid, 1, 'int32', endian);
    hdr.opthdr.nblocks  = fread(fid, 1, 'int64', endian);
    hdr.opthdr.nsamples = fread(fid, 1, 'int64', endian);
    hdr.opthdr.nsignals = fread(fid, 1, 'int32', endian);
  else
    hdr.opthdr = [];
  end
  
  % determine the number of samples for each channel
  hdr.nsamples = diff(hdr.offset);
  % the last one has to be determined by looking at the total data block length
  hdr.nsamples(end+1) = hdr.datasize - hdr.offset(end);
  % divide by the number of bytes in each channel
  hdr.nsamples = hdr.nsamples(:) ./ (hdr.depth(:)./8);
  
end % reading the rest of the header

switch action
  case 'read'
    dat = cell(length(chanindx), 1);
    currchan = 1;
    for i = 1:hdr.nsignals
      switch hdr.depth(i) % length in bit
        case 16
          slen = 2; % sample length, for fseek
          datatype = 'int16';
        case 32
          slen = 4; % sample length, for fseek
          datatype = 'single';
        case 64
          slen = 8; % sample length, for fseek
          datatype = 'double';
      end % case
      
      tmp = fread(fid, [1 hdr.nsamples(i)], datatype);
      if ~isempty(intersect(chanindx,i)) %only keep channels that are requested
        dat{currchan} = tmp;
        currchan = currchan + 1;
      end
    end % for
    
  case 'skip'
    dat = {};
    fseek(fid, hdr.datasize, 'cof');
    
  case 'none'
    dat = {};
    return;
end % switch action

