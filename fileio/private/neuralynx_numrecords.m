function [t] = neuralynx_numrecords(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for determining the number of records in a single channel Neualynx file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

headersize = 16384;
switch filetype(filename)
  case 'neuralynx_ncs'
    recordsize = 1044;  % in bytes
  case 'neuralynx_nse'
    recordsize = 112;   % in bytes
  case 'neuralynx_nts'
    recordsize = 8;     % in bytes
  otherwise
    error('unsupported filetype');
end

fid = fopen(filename, 'rb', 'ieee-le');
fseek(fid, 0, 'eof');
t = floor((ftell(fid) - headersize)/recordsize);
fclose(fid);
