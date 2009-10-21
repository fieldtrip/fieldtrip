function [vol] = ama2vol(ama)

% AMA2VOL
%
% Use as
%   vol = ama2vol(ama)

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: ama2vol.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.1  2008/01/28 19:56:13  roboos
% moved code from subfunction in prepare_bemmodel and prepare_vol_sens into seperate function
% renamed from convert_ama2vol
%

vol  = [];
ngeo = length(ama.geo);
for i=1:ngeo
  vol.bnd(i).pnt = ama.geo(i).pnt;
  vol.bnd(i).tri = ama.geo(i).dhk;
  vol.cond(i) = ama.geo(i).sigmam;
end
vol.mat = ama.bi;
npnt = size(vol.mat,2);
if size(vol.mat,1)<npnt
  vol.mat(npnt, npnt) = 0;    % it should be a square matrix
end
vol.mat  = sparse(vol.mat);   % convert to sparse for faster multiplications
vol.type = 'dipoli';

