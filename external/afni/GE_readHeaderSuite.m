function su_hdr = GE_readHeaderSuite(fid, byte_align)
%
% su_hdr = GE_readHeaderSuite(fid, byte_align)
%
% Loads the suite header from the file with filed id fid
% and returns it as a structure. 
% if byte_align = 1 then 32-bit alignment (SGI, LX2 format)
% if byte_align = 0 then 16-bit alignment (Sun, 5.X format)
%
%
% Souheil J. Inati
% Dartmouth College
% May 2000
% souheil.inati@dartmouth.edu
%

% define the structure and read
su_hdr = struct('su_id', fread(fid,4,'char'));                   % Suite ID
su_hdr = setfield(su_hdr, 'su_uniq', fread(fid,1,'int16'));      % The Make-Unique Flag
su_hdr = setfield(su_hdr, 'su_diskid', fread(fid,1,'char'));     % Disk ID
su_hdr = setfield(su_hdr, 'prodid', fread(fid,13,'char'));       % Product ID
su_hdr = setfield(su_hdr, 'su_verscre', fread(fid,2,'char'));    % Genesis Version
su_hdr = setfield(su_hdr, 'su_verscur', fread(fid,2,'char'));    % Genesis Version
su_hdr = setfield(su_hdr, 'su_checksum', fread(fid,1,'uint32')); % Suite  Record Checksum
su_hdr = setfield(su_hdr, 'su_padding', fread(fid,85,'char'));   % Spare Space
fseek(fid,1,0); % 16-bit alignment
if byte_align, fseek(fid,2,0); end % 32-bit alignment

return
