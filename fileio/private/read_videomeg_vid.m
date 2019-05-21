function out = read_videomeg_vid(filename, varargin)

% READ_VIDEOMEG_VID
%
% Use as
%   hdr = read_videomeg_vid(filename)
% or
%   dat = read_videomeg_vid(filename, hdr, begsample, endsample)
%
% See also READ_VIDEOMEG_AUD, LOAD_AUDIO0123, LOAD_VIDEO123

needhdr = numel(varargin)==0;
needdat = numel(varargin)==3;

% the following code is largely copied from LOAD_VIDEO0123

if needhdr
  fid   = fopen_or_error(filename, 'rb');
  magic = fread(fid, [1 length('ELEKTA_VIDEO_FILE')], 'uchar=>char');
  ver   = fread(fid, 1, 'uint32');
  
  assert(strcmp(magic, 'ELEKTA_VIDEO_FILE'));
  assert(any(ver==[1 2 3]));
  
  if ver==3
    site_id     = fread(fid, 1, 'uint8');
    sender_flag = fread(fid, 1, 'uint8');
    is_sender   = (sender_flag==1);
  else
    site_id     = -1;
    sender_flag = -1;
    is_sender   = -1;
  end
  
  ts         = [];
  chunk_id   = [];
  frame_size = [];
  frame_ptr  = [];
  frame      = [];
  
  while ~feof(fid)
    try
      ts(end+1) = fread(fid, 1, 'uint64=>uint64');
      if(ver > 1)
        chunk_id(end+1) = fread(fid, 1, 'uint64=>uint64');
      else
        chunk_id(end+1) = -1;
      end
      
      frame_size(end+1)   = fread(fid, 1, 'uint32');
      frame_ptr (end+1)   = ftell(fid);
      if isempty(frame)
        % read a single frame to determine the size
        frame = fread(fid, frame_size(end), 'uchar=>uchar');
      else
        fseek(fid, frame_size(end), 'cof');
      end
    end
  end
  fclose(fid);
  
  % we need to decompress a single JPEG image here to determine its size
  % see http://stackoverflow.com/questions/18659586/from-raw-bits-to-jpeg-without-writing-into-a-file
  
  % decode the JPEG frame using Java
  jImg = javax.imageio.ImageIO.read(java.io.ByteArrayInputStream(frame));
  dim(1) = jImg.getHeight;
  dim(2) = jImg.getWidth;
  
  % convert Java image to MATLAB image
  % note that it is a monochrome image, hence not RGB
  img = reshape(typecast(jImg.getData.getDataStorage, 'uint8'), dim(2), dim(1))';
  
  hdr = [];
  hdr.Fs          = 30;
  hdr.nChans      = prod(dim);
  hdr.nSamples    = numel(frame_ptr);
  hdr.nSamplesPre = 0;
  hdr.nTrials     = 1;
  for i=1:hdr.nChans
    hdr.label{i}  = num2str(i);
  end

  % compute the relation between timstamps and video frames
  hdr.FirstTimeStamp      = ts(1);
  hdr.TimeStampPerSample  = double((ts(end))-ts(1))./(numel(ts)-1);

  % also keep some original file details
  hdr.orig.site_id    = site_id;
  hdr.orig.ts         = ts;
  hdr.orig.chunk_id   = chunk_id;
  hdr.orig.frame_size = frame_size;
  hdr.orig.frame_ptr  = frame_ptr;
  hdr.orig.dim        = dim;
  hdr.orig.frame      = img;
  
  out = hdr;
  
elseif needdat
  
  hdr       = varargin{1};
  begsample = varargin{2};
  endsample = varargin{3};
  
  % get the header details, as if the code above was just executed
  dim         = hdr.orig.dim;
  frame       = hdr.orig.frame;
  frame_size  = hdr.orig.frame_size;
  frame_ptr   = hdr.orig.frame_ptr;
  
  nframes = endsample-begsample+1;
  dat     = zeros(prod(dim), nframes);
  
  fid = fopen_or_error(filename, 'rb');
  for i=1:nframes
    fseek(fid, frame_ptr(begsample-1+i), 'bof');
    frame = fread(fid, frame_size(begsample-1+i), 'uchar=>uchar');
    jImg  = javax.imageio.ImageIO.read(java.io.ByteArrayInputStream(frame));
    img   = reshape(typecast(jImg.getData.getDataStorage, 'uint8'), dim(2), dim(1))';
    dat(:,i) = img(:);
  end
  fclose(fid);
  
  out = dat;
  
end
