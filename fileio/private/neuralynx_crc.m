function crc = neuralynx_crc(dat, dim)

% NEURALYNX_CRC computes a cyclic redundancy check
%
% Use as
%   crc = neuralynx_crc(dat)
%
% Note that the CRC is computed along the first dimension.

% Copyright (C) 2007, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information


nchans   = size(dat,1);
nsamples = size(dat,2);

crc = zeros(1,nsamples,class(dat));

for i=1:nchans
  crc = bitxor(crc, dat(i,:));
end

