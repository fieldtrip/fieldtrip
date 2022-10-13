function pix_hdr = GE_readHeaderPixel(fid, byte_align)
%
% pix_hdr = read_pixel_header(fid, byte_align)
%
% Loads the pixel header from the file with filed id fid
% and returns it as a structure.
% if byte_align = 1 then 32-bit alignment (SGI, LX2 format)
% if byte_align = 0 then 16-bit alignment (Sun, 5.X format)
%
% Souheil J. Inati
% Dartmouth College
% May 2000
% souheil.inati@dartmouth.edu
%

% define the structure and read
pix_hdr = struct('img_magic', fread(fid,1,'int32'));
pix_hdr = setfield(pix_hdr, 'img_hdr_length', fread(fid,1,'int32'));% length of total header in bytes and
                                                                    % a byte displacement to the 'pixel data area'
pix_hdr = setfield(pix_hdr, 'img_width', fread(fid,1,'int32'));  % width (pixels) of image
pix_hdr = setfield(pix_hdr, 'img_height', fread(fid,1,'int32')); % height (pixels) of image
pix_hdr = setfield(pix_hdr, 'img_depth', fread(fid,1,'int32'));  % depth (1, 8, 16, or 24 bits) of pixel
pix_hdr = setfield(pix_hdr, 'img_compress', fread(fid,1,'int32'));   % type of compression; see IC_* below
pix_hdr = setfield(pix_hdr, 'img_dwindow', fread(fid,1,'int32'));    % default window setting
pix_hdr = setfield(pix_hdr, 'img_dlevel', fread(fid,1,'int32')); % default level setting
pix_hdr = setfield(pix_hdr, 'img_bgshade', fread(fid,1,'int32'));    % background shade to use for non-image
pix_hdr = setfield(pix_hdr, 'img_ovrflow', fread(fid,1,'int32'));    % overflow value
pix_hdr = setfield(pix_hdr, 'img_undflow', fread(fid,1,'int32'));    % underflow value
pix_hdr = setfield(pix_hdr, 'img_top_offset', fread(fid,1,'int32')); % number of blank lines at image top
pix_hdr = setfield(pix_hdr, 'img_bot_offset', fread(fid,1,'int32')); % number of blank lines at image bottom
pix_hdr = setfield(pix_hdr, 'img_version', fread(fid,1,'int16'));    % version of the header structure
if byte_align, fseek(fid,2,0); end % 32-bit alignment
pix_hdr = setfield(pix_hdr, 'img_checksum', fread(fid,1,'uint16')); % 16 bit end_around_carry sum of pixels
pix_hdr = setfield(pix_hdr, 'img_p_id', fread(fid,1,'int32'));   % a byte disp to unique image identifier
pix_hdr = setfield(pix_hdr, 'img_l_id', fread(fid,1,'int32'));   % byte length of unique image identifier
pix_hdr = setfield(pix_hdr, 'img_p_unpack', fread(fid,1,'int32'));   % a byte disp to 'unpack control'
pix_hdr = setfield(pix_hdr, 'img_l_unpack', fread(fid,1,'int32'));   % byte length of 'unpack control'
pix_hdr = setfield(pix_hdr, 'img_p_compress', fread(fid,1,'int32')); % a byte disp to 'compression control'
pix_hdr = setfield(pix_hdr, 'img_l_compress', fread(fid,1,'int32')); % byte length of 'compression control'
pix_hdr = setfield(pix_hdr, 'img_p_histo', fread(fid,1,'int32'));    % a byte disp to 'histogram control'
pix_hdr = setfield(pix_hdr, 'img_l_histo', fread(fid,1,'int32'));    % byte length of 'histogram control'
pix_hdr = setfield(pix_hdr, 'img_p_text', fread(fid,1,'int32')); % a byte disp to 'text plane data'
pix_hdr = setfield(pix_hdr, 'img_l_text', fread(fid,1,'int32')); % byte length of 'text plane data'
pix_hdr = setfield(pix_hdr, 'img_p_graphics', fread(fid,1,'int32')); % a byte disp to 'graphics plane data'
pix_hdr = setfield(pix_hdr, 'img_l_graphics', fread(fid,1,'int32')); % byte length of 'graphics plane data'
pix_hdr = setfield(pix_hdr, 'img_p_dbHdr', fread(fid,1,'int32'));    % a byte disp to 'data base header data'
pix_hdr = setfield(pix_hdr, 'img_l_dbHdr', fread(fid,1,'int32'));    % byte length of 'data base header data'
pix_hdr = setfield(pix_hdr, 'img_levelOffset', fread(fid,1,'int32'));% value to add to stored Pixel Data values
                                                                     % to get the correct presentation value
pix_hdr = setfield(pix_hdr, 'img_p_user', fread(fid,1,'int32')); % byte displacement to user defined data
pix_hdr = setfield(pix_hdr, 'img_l_user', fread(fid,1,'int32')); % byte length of user defined data

return
