function [t] = neuralynx_timestamp(filename, num)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for reading a single timestamp of a single channel Neuralynx file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

headersize = 16384;
switch filetype(filename)
  case 'neuralynx_ncs'
    recordsize = 1044;  % in bytes
  case 'neuralynx_nse'
    recordsize = 112;   % in bytes
  case 'neuralynx_nst'
    recordsize = 304;   % in bytes
  case 'neuralynx_nts'
    recordsize = 8;     % in bytes
  case 'neuralynx_ntt'
    recordsize = 304;   % in bytes
end

fid = fopen(filename, 'rb', 'ieee-le');

if (ispc)
  % this is to fix a bug in the windwos version which does not want to do uint64=>uint64
  % however this code will fail if the MSB is set (only likely in very long recordings)
  if ~isinf(num)
    % read the timestamp of the indicated record
    fseek(fid, headersize + (num-1)*recordsize, 'bof');
    t = fread(fid, 1, 'uint64=>integer*8');
  else
    % read the timestamp of the last record
    fseek(fid, -recordsize, 'eof');
    t = fread(fid, 1, 'uint64=>integer*8');
  end
  t = uint64(t);
else
  if ~isinf(num)
    % read the timestamp of the indicated record
    fseek(fid, headersize + (num-1)*recordsize, 'bof');
    t = fread(fid, 1, 'uint64=>uint64');
  else
    % read the timestamp of the last record
    fseek(fid, -recordsize, 'eof');
    t = fread(fid, 1, 'uint64=>uint64');
  end
end

fclose(fid);
