function [ts] = timestamp_plexon(tsl, tsh);

% TIMESTAMP_PLEXON merge the low and high part of the timestamps
% into a single uint64 value

if ~isa(tsl, 'uint32')
  error('invalid input');
elseif ~isa(tsh, 'uint16')
  error('invalid input');
end

% convert the 16 bit high timestamp into a 32 bit integer
dum = zeros(2, length(tsh), 'uint16');
if littleendian
  dum(1,:) = tsh(:);
else
  dum(2,:) = tsh(:);
end
tsh = typecast(dum(:), 'uint32');

% convert the 32 bit low and 32 bit high timestamp into a 64 bit integer
dum = zeros(2, length(tsh), 'uint32');
if littleendian
  dum(1,:) = tsl;
  dum(2,:) = tsh;
else
  dum(1,:) = tsh;
  dum(2,:) = tsl;
end

ts = typecast(dum(:), 'uint64');
ts = reshape(ts, size(tsl));

