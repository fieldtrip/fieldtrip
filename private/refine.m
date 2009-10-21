function [pntr, dhkr] = refine(pnt, dhk, method, varargin);

% REFINE a 3D surface that is described by a triangulation
%
% Use as
%   [pnt, tri] = refine(pnt, tri)
%   [pnt, tri] = refine(pnt, tri, 'updown', numtri)
%
% The default method is to refine the mesh globally by inserting a vertex at
% each edge according to the algorithm described in Banks, 1983.
%
% The alternative 'updown' method refines the mesh a couple of times
% using Banks' algorithm, followed by a downsampling using the REDUCEPATCH
% function.

% The Banks method is a memory efficient implementation which remembers
% the previously inserted vertices. The refinement algorithm executes in
% linear time with the number of triangles.

% Copyright (C) 2002-2005, Robert Oostenveld
%
% $Log: refine.m,v $
% Revision 1.4  2005/05/06 08:54:31  roboos
% added 'updown method, using matlab reducepatch function
%
% Revision 1.3  2003/03/11 15:35:20  roberto
% converted all files from DOS to UNIX
%
% Revision 1.2  2003/03/04 21:46:19  roberto
% added CVS log entry and synchronized all copyright labels
%

if nargin<3
  method = 'banks';
end

switch lower(method)
case 'banks'
  npnt   = size(pnt,1);
  ndhk   = size(dhk,1);
  insert = spalloc(3*npnt,3*npnt,3*ndhk);

  dhkr  = zeros(4*ndhk,3);		% allocate memory for the new triangles
  pntr  = zeros(npnt+3*ndhk,3);		% allocate memory for the maximum number of new vertices
  pntr(1:npnt,:) = pnt;			% insert the original vertices
  current = npnt;

  for i=1:ndhk

    if ~insert(dhk(i,1),dhk(i,2))
      current = current + 1;
      pntr(current,:) = (pnt(dhk(i,1),:) + pnt(dhk(i,2),:))/2;
      insert(dhk(i,1),dhk(i,2)) = current;
      insert(dhk(i,2),dhk(i,1)) = current;
      v12 = current;
    else
      v12 = insert(dhk(i,1),dhk(i,2));
    end

    if ~insert(dhk(i,2),dhk(i,3))
      current = current + 1;
      pntr(current,:) = (pnt(dhk(i,2),:) + pnt(dhk(i,3),:))/2;
      insert(dhk(i,2),dhk(i,3)) = current;
      insert(dhk(i,3),dhk(i,2)) = current;
      v23 = current;
    else
      v23 = insert(dhk(i,2),dhk(i,3));
    end

    if ~insert(dhk(i,3),dhk(i,1))
      current = current + 1;
      pntr(current,:) = (pnt(dhk(i,3),:) + pnt(dhk(i,1),:))/2;
      insert(dhk(i,3),dhk(i,1)) = current;
      insert(dhk(i,1),dhk(i,3)) = current;
      v31 = current;
    else
      v31 = insert(dhk(i,3),dhk(i,1));
    end

    % add the 4 new triangles with the correct indices
    dhkr(4*(i-1)+1, :) = [dhk(i,1) v12 v31];
    dhkr(4*(i-1)+2, :) = [dhk(i,2) v23 v12];
    dhkr(4*(i-1)+3, :) = [dhk(i,3) v31 v23];
    dhkr(4*(i-1)+4, :) = [v12 v23 v31];

  end

  % remove the space for the vertices that was not used
  pntr = pntr(1:current, :);

case 'updown'
  ndhk = size(dhk,1);
  while ndhk<varargin{1}
    % increase the number of triangles by a factor of 4
    [pnt, dhk] = refine(pnt, dhk, 'banks');
    ndhk = size(dhk,1);
  end
  % reduce number of triangles using Matlab function
  [dhkr, pntr] = reducepatch(dhk, pnt, varargin{1});

otherwise
  error(['unsupported method: ' method]);
end

