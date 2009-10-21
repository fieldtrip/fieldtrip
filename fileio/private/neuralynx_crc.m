function crc = neuralynx_crc(dat, dim)

% NEURALYNX_CRC computes a cyclic redundancy check
%
% Use as
%   crc = neuralynx_crc(dat)
%
% Note that the CRC is computed along the first dimension.

% Copyright (C) 2007, Robert Oostenveld
%
% $Log: neuralynx_crc.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.1  2007/02/20 08:52:04  roboos
% new implementation
%


nchans   = size(dat,1);
nsamples = size(dat,2);

crc = zeros(1,nsamples,class(dat));

for i=1:nchans
  crc = bitxor(crc, dat(i,:));
end

