function [vol] = ama2vol(ama)

% AMA2VOL
%
% Use as
%   vol = ama2vol(ama)

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: ama2vol.m,v $
% Revision 1.1  2008/12/24 10:25:41  roboos
% cleaned up the dipoli wrapper, replaced the binary by a better one and added a copy of the helper functions (from fileio)
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

