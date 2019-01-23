function [posr, trir, texr] = refine(pos, tri, method, varargin)

% REFINE a 3D surface that is described by a triangulation
%
% Use as
%   [pos, tri]          = refine(pos, tri)
%   [pos, tri]          = refine(pos, tri, 'banks')
%   [pos, tri, texture] = refine(pos, tri, 'banks', texture)
%   [pos, tri]          = refine(pos, tri, 'updown', numtri)
%
% If no method is specified, the default is to refine the mesh globally by bisecting
% each edge according to the algorithm described in Banks, 1983.
%
% The Banks method allows the specification of a subset of triangles to be refined
% according to Banks' algorithm. Adjacent triangles will be gracefully dealt with.
%
% The alternative 'updown' method refines the mesh a couple of times
% using Banks' algorithm, followed by a downsampling using the REDUCEPATCH
% function.
%
% If the textures of the vertices are specified, the textures for the new
% vertices are computed
%
% The Banks method is a memory efficient implementation which remembers the
% previously inserted vertices. The refinement algorithm executes in linear
% time with the number of triangles. It is mentioned in
% http://www.cs.rpi.edu/~flaherje/pdf/fea8.pdf, which also contains the original
% reference.

% Copyright (C) 2002-2014, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%

if nargin<3
  method = 'banks';
end

texture = [];
numtri  = [];

if nargin>3
  switch lower(method)
    case 'banks'
      texture = varargin{1};
    case 'updown'
      numtri = varargin{1};
  end % switch
end

switch lower(method)
  case 'banks'
    if ~isempty(texture)
      npnt   = size(pos,1);
      ntri   = size(tri,1);
      ntex   = size(texture,1);
      
      assert(ntex==ntri, 'invalid size of texture');
      
      insert = spalloc(3*npnt,3*npnt,3*ntri);
      trir  = zeros(4*ntri,3);      % allocate memory for the new triangles
      posr  = zeros(npnt+3*ntri,3); % allocate memory for the maximum number of new vertices
      texr  = zeros(ntex+3*ntri,2);
      posr(1:npnt,:) = pos;         % insert the original vertices
      texr(1:ntex,:) = texture;
      current = npnt;
      
      for i=1:ntri
        
        if ~insert(tri(i,1),tri(i,2))
          current = current + 1;
          posr(current,:) = (pos(tri(i,1),:) + pos(tri(i,2),:))/2;
          texr(current,:) = (texture(tri(i,1),:) + texture(tri(i,2),:))/2;
          insert(tri(i,1),tri(i,2)) = current;
          insert(tri(i,2),tri(i,1)) = current;
          v12 = current;
        else
          v12 = insert(tri(i,1),tri(i,2));
        end
        
        if ~insert(tri(i,2),tri(i,3))
          current = current + 1;
          posr(current,:) = (pos(tri(i,2),:) + pos(tri(i,3),:))/2;
          texr(current,:) = (texture(tri(i,2),:) + texture(tri(i,3),:))/2;
          insert(tri(i,2),tri(i,3)) = current;
          insert(tri(i,3),tri(i,2)) = current;
          v23 = current;
        else
          v23 = insert(tri(i,2),tri(i,3));
        end
        
        if ~insert(tri(i,3),tri(i,1))
          current = current + 1;
          posr(current,:) = (pos(tri(i,3),:) + pos(tri(i,1),:))/2;
          texr(current,:) = (texture(tri(i,3),:) + texture(tri(i,1),:))/2;
          insert(tri(i,3),tri(i,1)) = current;
          insert(tri(i,1),tri(i,3)) = current;
          v31 = current;
        else
          v31 = insert(tri(i,3),tri(i,1));
        end
        
        % add the 4 new triangles with the correct indices
        trir(4*(i-1)+1, :) = [tri(i,1) v12 v31];
        trir(4*(i-1)+2, :) = [tri(i,2) v23 v12];
        trir(4*(i-1)+3, :) = [tri(i,3) v31 v23];
        trir(4*(i-1)+4, :) = [v12 v23 v31];
        
      end
      posr = posr(1:current, :);
      texr = texr(1:current, :);
      
    else
      % there is no texture
      
      npnt   = size(pos,1);
      ntri   = size(tri,1);
      insert = spalloc(3*npnt,3*npnt,3*ntri);
      
      trir  = zeros(4*ntri,3);      % allocate memory for the new triangles
      posr  = zeros(npnt+3*ntri,3); % allocate memory for the maximum number of new vertices
      posr(1:npnt,:) = pos;         % insert the original vertices
      current = npnt;
      
      for i=1:ntri
        
        if ~insert(tri(i,1),tri(i,2))
          current = current + 1;
          posr(current,:) = (pos(tri(i,1),:) + pos(tri(i,2),:))/2;
          insert(tri(i,1),tri(i,2)) = current;
          insert(tri(i,2),tri(i,1)) = current;
          v12 = current;
        else
          v12 = insert(tri(i,1),tri(i,2));
        end
        
        if ~insert(tri(i,2),tri(i,3))
          current = current + 1;
          posr(current,:) = (pos(tri(i,2),:) + pos(tri(i,3),:))/2;
          insert(tri(i,2),tri(i,3)) = current;
          insert(tri(i,3),tri(i,2)) = current;
          v23 = current;
        else
          v23 = insert(tri(i,2),tri(i,3));
        end
        
        if ~insert(tri(i,3),tri(i,1))
          current = current + 1;
          posr(current,:) = (pos(tri(i,3),:) + pos(tri(i,1),:))/2;
          insert(tri(i,3),tri(i,1)) = current;
          insert(tri(i,1),tri(i,3)) = current;
          v31 = current;
        else
          v31 = insert(tri(i,3),tri(i,1));
        end
        
        % add the 4 new triangles with the correct indices
        trir(4*(i-1)+1, :) = [tri(i,1) v12 v31];
        trir(4*(i-1)+2, :) = [tri(i,2) v23 v12];
        trir(4*(i-1)+3, :) = [tri(i,3) v31 v23];
        trir(4*(i-1)+4, :) = [v12 v23 v31];
      end
      % remove the space for the vertices that was not used
      posr = posr(1:current, :);
    end
    
  case 'updown'
    ntri = size(tri,1);
    while ntri<numtri
      % increase the number of triangles by a factor of 4
      [pos, tri] = refine(pos, tri, 'banks');
      ntri = size(tri,1);
    end
    % reduce number of triangles using MATLAB function
    [trir, posr] = reducepatch(tri, pos, numtri);
    
  otherwise
    ft_error('unsupported method "%s"', method);
end

