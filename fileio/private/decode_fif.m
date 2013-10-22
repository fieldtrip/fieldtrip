function hdr = decode_fif(chunk)

% DECODE_FIF is a helper function for real-time processing of Neuromag data. This
% function is used to decode the content of the optional neuromag_fif chunck.
%
% See also DECODE_RES4, DECODE_NIFTI1, SAP2MATLAB

% the binary blob was created on the little-endian Intel Linux acquisition
% computer, whereas the default for fiff files is that they are stored in
% big-endian byte order. MATLAB is able to swap the bytes on the fly by specifying
% 'le" or "be" to fopen. The normal MNE fiff_open function assumes that it is big
% endian, hence here we have to open it as little endian.

% decode the fif chunck content, first write the binary blob to disk, then use the standard
% reading functions to decode the blob
filename = tempname;

% wtite the binary blob to disk, byte-by-byte to avoid any swapping between little and big-endian content
F = fopen(filename, 'w');
fwrite(F, chunk, 'uint8');
fclose(F);

% calling ft_read_header rather than the low-level MNE functions ensures that mne2grad
% and other auxiliary functions will be called
hdr = ft_read_header(filename, 'headerformat', 'neuromag_mne', 'endian', 'L');

% clean up the temporary file
delete(filename);
