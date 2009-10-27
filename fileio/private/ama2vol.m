function [vol] = ama2vol(ama)

% AMA2VOL
%
% Use as
%   vol = ama2vol(ama)

% Copyright (C) 2008, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

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

