function [su_hdr,ex_hdr,se_hdr,im_hdr,pix_hdr,im_offset] = GE_readHeader(IFileName)
%
% [su_hdr,ex_hdr,se_hdr,im_hdr,pix_hdr,im_offset] = GE_readHeader(IFileName)
% Reads the header info from a GE lx2 or 5.X file 
%
%
% Souheil J. Inati
% Dartmouth College
% May 2000
% souheil.inati@dartmouth.edu
%

%%%% Open the IFile %%%%%
[fid,message] = fopen(IFileName,'r','b');         %note: 'b' = Big-Endian format
if fid == -1
  error(message)
end

% Check for the magic number at the first word
magic = fread(fid,1,'int32');
if magic == 1229801286
  % This is a 5.X format image
  byte_align = 0;
  pix_hdr_offset = 0;
  su_hdr_size = 114;
  ex_hdr_size = 1024;
  se_hdr_size = 1020;
  im_hdr_size = 1022;
else
  % Magic no. not at the 1st word, try at the LX2 position
  fseek(fid,3228,-1);
  magic = fread(fid,1,'int32');
  if magic == 1229801286
    % This is an LX2 format image
    byte_align = 1;
    pix_hdr_offset = 3228;
    su_hdr_size = 116;
    ex_hdr_size = 1040;
    se_hdr_size = 1028;
    im_hdr_size = 1044;
  else
    msg = 'This is not a 5.X or LX2 format image.  No Magic Number.';
    error(msg);
  end
end

% Load the pixel header
fseek(fid,pix_hdr_offset,-1);
pix_hdr = GE_readHeaderPixel(fid, byte_align);

% Compute the offsets
su_hdr_offset = pix_hdr.img_p_dbHdr;
ex_hdr_offset = su_hdr_offset + su_hdr_size;
se_hdr_offset = ex_hdr_offset + ex_hdr_size;
im_hdr_offset = se_hdr_offset + se_hdr_size;
im_offset = pix_hdr_offset + pix_hdr.img_hdr_length;

% Check for epirecon images
if (pix_hdr.img_l_dbHdr == 0 & byte_align==0)
  error('This is a epirecon image. No header')
end

% Load the suite header
fseek(fid,su_hdr_offset,-1);
su_hdr = GE_readHeaderSuite(fid, byte_align);

% Load the exam header
fseek(fid,ex_hdr_offset,-1);
ex_hdr = GE_readHeaderExam(fid, byte_align);

% Load the series header
fseek(fid,se_hdr_offset,-1);
se_hdr = GE_readHeaderSeries(fid, byte_align);

% Load the image header
fseek(fid,im_hdr_offset,-1);
im_hdr = GE_readHeaderImage(fid, byte_align);

% Close the file
fclose(fid);

return





