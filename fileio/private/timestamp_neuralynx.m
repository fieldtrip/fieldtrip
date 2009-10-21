function [ts] = timestamp_neuralynx(tsl, tsh);

% TIMESTAMP_NEURALYNX merge the low and high part of Neuralynx timestamps
% into a single uint64 value

if ~isa(tsl, 'uint32')
  error('invalid input');
elseif ~isa(tsh, 'uint32')
  error('invalid input');
end

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
