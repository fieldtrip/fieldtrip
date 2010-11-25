function sev = read_tdt_sev(filename, dtype, begsample, endsample)

% READ_TDT_SEV
%
% Use as
%   sev = read_tdt_sev(filename, dtype, begsample, endsample)
%
% Note: sev files contain raw broadband data that is streamed to the RS4


fid = fopen(filename, 'rb', 'ieee-le');

if nargin<3
  begsample = 1;
end

if nargin<4
  endsample = inf;
end

switch dtype
  case 0
    fmt   = 'float32';
    wsize = 4;
  case 1
    fmt   = 'int32';
    wsize = 4;
  case 2
    fmt   = 'int16';
    wsize = 2;
  case 3
    fmt   = 'int8';
    wsize = 1;
  case 4
    fmt   = 'double';
    wsize = 8;
  case 5
    error('don''t know what a DFORM_QWORD is');
  otherwise
    error('unknown dtype');
end

fseek(fid, begsample*wsize, 'cof');
sev = fread(fid, endsample-begsample, fmt);
fclose(fid);

