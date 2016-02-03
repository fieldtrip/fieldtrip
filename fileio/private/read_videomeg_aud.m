function out = read_videomeg_aud(filename, varargin)

% READ_VIDEOMEG_AUD
%
% Use as
%   hdr = read_videomeg_aud(filename)
% or
%   dat = read_videomeg_aud(filename, hdr, begsample, endsample)
%
% See also READ_VIDEOMEG_VID, LOAD_AUDIO0123, LOAD_VIDEO123

needhdr = numel(varargin)==0;
needdat = numel(varargin)==3;

% the following code is largely copied from LOAD_AUDIO0123

if needhdr
  fid    = fopen(filename, 'rb');
  magic  = fread(fid, [1 length('ELEKTA_AUDIO_FILE')], 'uchar=>char');
  ver    = fread(fid, 1, 'uint32');
  srate  = fread(fid, 1, 'uint32');
  nchans = fread(fid, 1, 'uint32');
  assert(strcmp(magic, 'ELEKTA_AUDIO_FILE'));
  assert(any(ver==[0 1 2 3]));
  
  dt_start = ftell(fid); % this is where the data starts
  
  if(ver==2 || ver==3)
    attrib_sz = 8 + 8 + 4;  % timestamp, id, blklen
  else
    attrib_sz = 8 + 4;      % timestamp, blklen
  end
  
  ts = fread(fid, 1, 'uint64=>uint64');
  
  if(ver==2 || ver==3)
    id = fread(fid, 1, 'uint64=>uint64');
  else
    id = -1;
  end
  
  blklen = fread(fid, 1, 'uint32'); % the size of the data in each block in bytes
  blksmp = blklen/(nchans*2);       % the size of the data in each block in samples
  
  fseek(fid, dt_start, 'bof');
  ts = fread(fid, inf, '1*uint64=>uint64', blklen+attrib_sz-8);
  
  fseek(fid, 0 ,'eof');
  dt_end = ftell(fid); % this is where the data ends
  fclose(fid);
  
  numblk = (dt_end-dt_start)/(blklen+attrib_sz);
  
  % convert it into a FieldTrip header
  hdr         = [];
  hdr.Fs      = srate;
  hdr.nChans  = nchans;
  for i=1:nchans
    hdr.label{i}    = num2str(i);
    hdr.chantype{i} = 'audio';
    hdr.chanunit{i} = 'unknown';
  end
  hdr.nTrials     = 1;
  hdr.nSamples    = numblk*blksmp;
  hdr.nSamplesPre = 0;
  
  % compute the relation between timstamps and audio samples
  % the timestamp is specified at the start of each block
  hdr.FirstTimeStamp      = ts(1);
  hdr.TimeStampPerSample  = double(ts(end)-ts(1))/((numblk-1)*blksmp);
  
  % also keep some original file details
  hdr.orig.magic      = magic;
  hdr.orig.ver        = ver;
  hdr.orig.srate      = srate;
  hdr.orig.nchans     = nchans;
  hdr.orig.dt_start   = dt_start;
  hdr.orig.attrib_sz  = attrib_sz;
  hdr.orig.ts         = ts;
  hdr.orig.id         = id;
  hdr.orig.blklen     = blklen;
  hdr.orig.blksmp     = blksmp;
  hdr.orig.numblk     = numblk;
  
  % return the header as output
  out = hdr;
  
elseif needdat
  
  hdr       = varargin{1};
  begsample = varargin{2};
  endsample = varargin{3};
  
  % get the header details, as if the code above was just executed
  magic      = hdr.orig.magic;
  ver        = hdr.orig.ver;
  srate      = hdr.orig.srate;
  nchans      = hdr.orig.nchans;
  attrib_sz  = hdr.orig.attrib_sz;
  dt_start   = hdr.orig.dt_start;
  ts         = hdr.orig.ts;
  id         = hdr.orig.id;
  blklen     = hdr.orig.blklen;
  
  % data is organized in blocks of "blklen" bytes, with each block having a header of "attrib_sz" bytes plus the data.
  
  % size of the data in each block, expressed as nchans*nsamples
  blksiz = blklen/(nchans*2);
  
  % determine the first and last block that span the requested data
  begblk = floor((begsample-1)/blksiz + 1);
  endblk = floor((endsample-1)/blksiz + 1);
  
  % determine the start and end pointer in the file
  begptr = dt_start + (begblk-1) * (blklen+attrib_sz);
  endptr = dt_start + (endblk  ) * (blklen+attrib_sz);
  
  % read the required blocks of data from the file
  fid = fopen(filename, 'rb');
  fseek(fid, begptr, 'bof');
  blk = fread(fid, [1 endptr-begptr], 'uint8=>uint8');
  fclose(fid);
  
  blk = reshape(blk, (blklen+attrib_sz), []);
  
  if(ver==2 || ver==3)
    ts = blk(1:8,:);
    ts = typecast(ts(:), 'uint64');
    id = blk(9:16,:);
    id = typecast(id(:), 'uint64');
    bs = blk(17:20,:);
    bs = typecast(bs(:), 'uint32');
    dat = blk(21:end,:);
    dat = typecast(dat(:), 'int16');
  else
    ts = blk(1:8,:);
    ts = typecast(ts(:), 'uint64');
    id = -1*ones(size(ts));
    bs = blk(9:12,:);
    bs = typecast(bs(:), 'uint32');
    dat = blk(13:end,:);
    dat = typecast(dat(:), 'int16');
    
  end
  
  dat = reshape(dat, 2, []); % FIXME needs to be checked, it might have to be transposed
  
  % cut out the desired section of data based on sample numbers
  begsample = begsample - (begblk-1)*blksiz;
  endsample = endsample - (begblk-1)*blksiz;
  dat = dat(:, begsample:endsample);
  
  % return the data as output
  out = dat;
end
