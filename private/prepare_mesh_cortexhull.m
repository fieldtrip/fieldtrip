function headshape = prepare_mesh_cortexhull(cfg)

% PREPARE_MESH_CORTEXHULL creates a mesh representing the cortex hull, i.e.
% the smoothed envelope around the pial surface created by FreeSurfer
%
% This function relies on the FreeSurfer and iso2mesh software packages
%
% Configuration options:
%   cfg.headshape    = a filename containing the pial surface computed by
%                      FreeSurfer recon-all ('/path/to/surf/lh.pial')
%   cfg.fshome       = FreeSurfer folder location
%                      (default: '/Applications/freesurfer')
%   cfg.resolution   = resolution of the volume delimited by headshape being
%                      floodfilled by mris_fill (default: 1)
%   cfg.outer_surface_sphere = diameter of the sphere used by make_outer_surface
%                      to close the sulci using morphological operations (default: 15)
%   cfg.smooth_steps = number of standard smoothing iterations (default: 0)
%   cfg.laplace_steps = number of Laplacian (non-shrinking) smoothing
%                      iterations (default: 2000)
%   cfg.fixshrinkage = reduce possible shrinkage due to smoothing (default: 'no')
%   cfg.expansion_mm = amount in mm with which the hull is re-expanded, applies
%                      when cfg.fixshrinkage = 'yes' (default: 'auto')
%
% See also FT_PREPARE_MESH

% Copyright (C) 2012-2018, Arjen Stolk, Gio Piantoni, Andrew Dykstra
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
% $Id$


% get the default options
surf                 = ft_getopt(cfg, 'headshape');
fshome               = ft_getopt(cfg, 'fshome', '/Applications/freesurfer');
resolution           = ft_getopt(cfg, 'resolution', 1);
outer_surface_sphere = ft_getopt(cfg, 'outer_surface_sphere', 15);
smooth_steps         = ft_getopt(cfg, 'smooth_steps', 0); % previous default was 60
laplace_steps        = ft_getopt(cfg, 'laplace_steps', 2000); % replaces smoothing using mris_smooth
fixshrinkage         = ft_getopt(cfg, 'fixshrinkage', 'no');
expansion_mm         = ft_getopt(cfg, 'expansion_mm', 'auto'); % applies when fixshrinkage is 'yes'

% add the FreeSurfer environment
fprintf('adding the FreeSurfer environment\n')
addpath([fshome '/matlab']); % where make_outer_surface is located
setenv('FREESURFER_HOME', fshome);
PATH = getenv('PATH');
setenv('PATH', [PATH ':' fshome '/bin']); % where mris_fill is located

% temporary files
surf_filled   = [tempname() '_pial.filled.mgz'];
surf_outer    = [tempname() '_pial_outer'];
surf_smooth   = [tempname() '_pial_smooth'];

% fill the pial volume
cmd = sprintf('mris_fill -c -r %d %s %s', resolution, surf, surf_filled);
system(['source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' cmd]);

% make outer surface by closing the gaps (modified)
make_outer_surface(surf_filled, outer_surface_sphere, surf_outer)

% smooth the surface using FreeSurfer's mris_smooth (not applied by default)
cmd = sprintf('mris_smooth -nw -n %d %s %s', smooth_steps, surf_outer, surf_smooth);
system(['source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' cmd]);
headshape = ft_read_headshape(surf_smooth);

% smooth the surface using iso2mesh's smoothsurf
if laplace_steps >= 1
  ft_hastoolbox('iso2mesh',1);
  fprintf('Laplacian smoothing for %d iterations\n', laplace_steps)
  conn = meshconn(headshape.tri, size(headshape.pos,1)); % determine neighbors
  headshape.pos = smoothsurf(headshape.pos, [], conn, laplace_steps, 0, 'laplacianhc', .2);
end

% fix shrinkage if needed 
% FIXME: might be too generous of a hull and ideally involves local rather than global expansion
if strcmp(fixshrinkage, 'yes')
  try
    pial = ft_read_headshape(cfg.headshape);
    inside = intriangulation(pial.pos, pial.tri, headshape.pos); % determine inside/outside points
    if numel(find(inside==0)) > numel(find(inside==1)) % assuming the majority of points are inside
      inside = ~inside;
    end
    if any(inside==0)
      fprintf('expanding the hull to include more outside points\n');
      expansion = zeros(size(headshape.pos,1),1);
      idx = find(inside==0); % outside point indices
      for p = 1:numel(idx) % determine pial-hull distance for outside points
        expansion(idx(p)) = min(sqrt( (headshape.pos(idx(p),1)-pial.pos(:,1)).^2 + ...
          (headshape.pos(idx(p),2)-pial.pos(:,2)).^2 + (headshape.pos(idx(p),3)-pial.pos(:,3)).^2 ));
      end
      if strcmp(expansion_mm, 'auto')
        expansion_mm = mean(expansion(idx)); % mean outside distance
      else
        expansion_mm = cfg.expansion_mm;
      end
      surf_expanded = [tempname() '_pial_expanded'];
      write_surf(surf_expanded, headshape.pos, headshape.tri);
      cmd = sprintf('mris_expand %s %d %s', surf_expanded, expansion_mm, surf_expanded); % global expansion
      system(['source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' cmd]);
      headshape = ft_read_headshape(surf_expanded);
      %     % fix shrinkage locally - FIXME: output showing expansion at unexpected locations
      %     thicknessfile = [fileparts(tempname()) filesep 'expansion'];
      %     write_curv(thicknessfile, expansion, size(headshape.pos,1)); % write thickness file used for expansion
      %     pialfile = [fileparts(tempname()) filesep 'rh.pial']; % unclear why this file is also needed
      %     write_surf(pialfile, headshape.pos, headshape.tri);
      %     spherefile = [fileparts(tempname()) filesep 'rh.sphere']; % unclear why this file is also needed
      %     write_surf(spherefile, headshape.pos, headshape.tri);
      %     cmd = sprintf('mris_expand -thickness -thickness_name %s %s %d %s', thicknessfile, surf_smooth, -1, surf_expanded);
      %     system(['source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' cmd]);
      %     headshape = ft_read_headshape(surf_expanded);
    else
      fprintf('no outside hull points found\n');
    end
  catch
    fprintf('hull expansion failed\n');
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_outer_surface (filled_volume, se_diameter, output_surface)

% Original Author: Marie Schaer
% Date: 2007/11/14
%
% This function takes as an input the binary volume resulting of the
% filling of the any surface (usually the pial one) using mris_fill, and
% will close the sulci using morphological operation, with a sphere as the
% structural element.
%
% Parameters:
% se_diameter is the diameter of the sphere (in mm), use 15mm by default
% to close the sulci.
%
% Utilities to write the surface to the freesurfer format are modified from
% "freesurfer_write_surf" (from Darren Weber's bioelectromagnetism toolbox),
% according to a suggestion by Don Hagler in FreeSurfer's mailing list on
% August 3, 2007.
%
% Example: make_outer_surface('lh.pial.mgz',15,'lh.outer-pial')


fprintf('reading filled volume...\n');
vol=MRIread(filled_volume);
volume=vol.vol;
volume(volume==1)=255;
fprintf('closing volume...\n');

% first apply a very soft gaussian filter, with sigma = 1mm, in order to
% facilitate the closing
Gaussian = fspecial('gaussian',[2 2],1);
image_f=zeros(256,256,256);
for slice=1:256
  temp = double(reshape(volume(:,:,slice),256,256));
  image_f(:,:,slice) = conv2(temp,Gaussian,'same');
end
image2=zeros(size(image_f));
image2(image_f<=25)=0;
image2(image_f>25)=255;

se=strel('sphere',se_diameter); % changed from 'ball' into 'sphere' (Stolk, Dec 2018)
BW2=imclose(image2,se);
thresh = max(BW2(:))/2;
i=find(BW2<=thresh);
BW2(i)=0;
i=find(BW2>thresh);
BW2(i)=255;

[f,v] = isosurface(BW2,100);

v2=[129-v(:,1) v(:,3)-129 129-v(:,2)]; % in order to cope with the different orientation
v=v2;

fprintf('morphological closing done.\n');
fprintf('writing outer surface...\n');

fname=output_surface;
vert = v;
face = f - 1;
vnum = size(vert,1);
fnum = size(face,1);

% open it as a big-endian file
fid = fopen(fname, 'wb', 'b');
TRIANGLE_FILE_MAGIC_NUMBER = 16777214;
fwrite3(fid, TRIANGLE_FILE_MAGIC_NUMBER);

% Ouput a couple of text lines with creation date
str = sprintf('created from matlab on %s\n',datestr(now));
fwrite(fid, str,'char');
fwrite(fid, vnum,'int32');
fwrite(fid, fnum,'int32');

% reshape vert into column array and write
vert = reshape(vert',size(vert,1)*size(vert,2),1);
fwrite(fid, vert,'float32');

% reshape face into column array and write
face = reshape(face',size(face,1)*size(face,2),1);
fwrite(fid, face,'int32');
fclose(fid) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function in = intriangulation(vertices,faces,testp,heavytest)
% intriangulation: Test points in 3d wether inside or outside a (closed) triangulation
% usage: in = intriangulation(vertices,faces,testp,heavytest)
%
% arguments: (input)
%  vertices   - points in 3d as matrix with three columns
%
%  faces      - description of triangles as matrix with three columns.
%               Each row contains three indices into the matrix of vertices
%               which gives the three cornerpoints of the triangle.
%
%  testp      - points in 3d as matrix with three columns
%
%  heavytest  - int n >= 0. Perform n additional randomized rotation tests.
%
% IMPORTANT: the set of vertices and faces has to form a watertight surface!
%
% arguments: (output)
%  in - a vector of length size(testp,1), containing 0 and 1.
%       in(nr) =  0: testp(nr,:) is outside the triangulation
%       in(nr) =  1: testp(nr,:) is inside the triangulation
%       in(nr) = -1: unable to decide for testp(nr,:)
%
% Thanks to Adam A for providing the FEX submission voxelise. The
% algorithms of voxelise form the algorithmic kernel of intriangulation.
%
% Thanks to Sven to discussions about speed and avoiding problems in
% special cases.
%
% Example usage:
%
%      n = 10;
%      vertices = rand(n, 3)-0.5; % Generate random points
%      tetra = delaunayn(vertices); % Generate delaunay triangulization
%      faces = freeBoundary(TriRep(tetra,vertices)); % use free boundary as triangulation
%      n = 1000;
%      testp = 2*rand(n,3)-1; % Generate random testpoints
%      in = intriangulation(vertices,faces,testp);
%      % Plot results
%      h = trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3));
%      set(h,'FaceColor','black','FaceAlpha',1/3,'EdgeColor','none');
%      hold on;
%      plot3(testp(:,1),testp(:,2),testp(:,3),'b.');
%      plot3(testp(in==1,1),testp(in==1,2),testp(in==1,3),'ro');
%
% See also: intetrahedron, tsearchn, inpolygon
%
% Author: Johannes Korsawe, heavily based on voxelise from Adam A.
% E-mail: johannes.korsawe@volkswagen.de
% Release: 1.3
% Release date: 25/09/2013

% check number of inputs
if nargin<3,
  fprintf('??? Error using ==> intriangulation\nThree input matrices are needed.\n');in=[];return;
end
if nargin==3,
  heavytest = 0;
end
% check size of inputs
if size(vertices,2)~=3 || size(faces,2)~=3 || size(testp,2)~=3,
  fprintf('??? Error using ==> intriagulation\nAll input matrices must have three columns.\n');in=[];return;
end
ipmax = max(faces(:));zerofound = ~isempty(find(faces(:)==0, 1));
if ipmax>size(vertices,1) || zerofound,
  fprintf('??? Error using ==> intriangulation\nThe triangulation data is defect. use trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3)) for test of deficiency.\n');return;
end

% loop for heavytest
inreturn = zeros(size(testp,1),1);VER = vertices;TESTP = testp;

for n = 1:heavytest+1,
  
  % Randomize
  if n>1,
    v=rand(1,3);D=rotmatrix(v/norm(v),rand*180/pi);vertices=VER*D;testp = TESTP*D;
  else,
    vertices=VER;
  end
  
  % Preprocessing data
  meshXYZ = zeros(size(faces,1),3,3);
  for loop = 1:3,
    meshXYZ(:,:,loop) = vertices(faces(:,loop),:);
  end
  
  % Basic idea (ingenious from FeX-submission voxelise):
  % If point is inside, it will cross the triangulation an uneven number of times in each direction (x, -x, y, -y, z, -z).
  
  % The function VOXELISEinternal is about 98% identical to its version inside voxelise.m.
  % This includes the elaborate comments. Thanks to Adam A!
  
  % z-direction:
  % intialization of results and correction list
  [in,cl] = VOXELISEinternal(testp(:,1),testp(:,2),testp(:,3),meshXYZ);
  
  % x-direction:
  % has only to be done for those points, that were not determinable in the first step --> cl
  [in2,cl2] = VOXELISEinternal(testp(cl,2),testp(cl,3),testp(cl,1),meshXYZ(:,[2,3,1],:));
  % Use results of x-direction that determined "inside"
  in(cl(in2==1)) = 1;
  % remaining indices with unclear result
  cl = cl(cl2);
  
  % y-direction:
  % has only to be done for those points, that were not determinable in the first and second step --> cl
  [in3,cl3] = VOXELISEinternal(testp(cl,3),testp(cl,1),testp(cl,2),meshXYZ(:,[3,1,2],:));
  
  % Use results of y-direction that determined "inside"
  in(cl(in3==1)) = 1;
  % remaining indices with unclear result
  cl = cl(cl3);
  
  % mark those indices, where all three tests have failed
  in(cl) = -1;
  
  if n==1,
    inreturn = in;  % Starting guess
  else,
    % if ALWAYS inside, use as inside!
    %        I = find(inreturn ~= in);
    %        inreturn(I(in(I)==0)) = 0;
    
    % if AT LEAST ONCE inside, use as inside!
    I = find(inreturn ~= in);
    inreturn(I(in(I)==1)) = 1;
    
  end
  
end

in = inreturn;

%==========================================================================
function [OUTPUT,correctionLIST] = VOXELISEinternal(testx,testy,testz,meshXYZ)

% Prepare logical array to hold the logical data:
OUTPUT = false(size(testx,1),1);

%Identify the min and max x,y coordinates of the mesh:
meshZmin = min(min(meshXYZ(:,3,:)));meshZmax = max(max(meshXYZ(:,3,:)));

%Identify the min and max x,y,z coordinates of each facet:
meshXYZmin = min(meshXYZ,[],3);meshXYZmax = max(meshXYZ,[],3);

%======================================================
% TURN OFF DIVIDE-BY-ZERO WARNINGS
%======================================================
%This prevents the Y1predicted, Y2predicted, Y3predicted and YRpredicted
%calculations creating divide-by-zero warnings.  Suppressing these warnings
%doesn't affect the code, because only the sign of the result is important.
%That is, 'Inf' and '-Inf' results are ok.
%The warning will be returned to its original state at the end of the code.
warningrestorestate = warning('query', 'MATLAB:divideByZero');
%warning off MATLAB:divideByZero

%======================================================
% START COMPUTATION
%======================================================

correctionLIST = [];   %Prepare to record all rays that fail the voxelisation.  This array is built on-the-fly, but since
%it ought to be relatively small should not incur too much of a speed penalty.

% Loop through each testpoint.
% The testpoint-array will be tested by passing rays in the z-direction through
% each x,y coordinate of the testpoints, and finding the locations where the rays cross the mesh.
facetCROSSLIST = zeros(1,1e3);  % uses countindex: nf
nm = size(meshXYZmin,1);
for loop = 1:length(OUTPUT),
  
  nf = 0;
  %    % - 1a - Find which mesh facets could possibly be crossed by the ray:
  %    possibleCROSSLISTy = find( meshXYZmin(:,2)<=testy(loop) & meshXYZmax(:,2)>=testy(loop) );
  
  %    % - 1b - Find which mesh facets could possibly be crossed by the ray:
  %    possibleCROSSLIST = possibleCROSSLISTy( meshXYZmin(possibleCROSSLISTy,1)<=testx(loop) & meshXYZmax(possibleCROSSLISTy,1)>=testx(loop) );
  
  % Do - 1a - and - 1b - faster
  possibleCROSSLISTy = find((testy(loop)-meshXYZmin(:,2)).*(meshXYZmax(:,2)-testy(loop))>0);
  possibleCROSSLISTx = (testx(loop)-meshXYZmin(possibleCROSSLISTy,1)).*(meshXYZmax(possibleCROSSLISTy,1)-testx(loop))>0;
  possibleCROSSLIST = possibleCROSSLISTy(possibleCROSSLISTx);
  
  if isempty(possibleCROSSLIST)==0  %Only continue the analysis if some nearby facets were actually identified
    
    % - 2 - For each facet, check if the ray really does cross the facet rather than just passing it close-by:
    
    % GENERAL METHOD:
    % 1. Take each edge of the facet in turn.
    % 2. Find the position of the opposing vertex to that edge.
    % 3. Find the position of the ray relative to that edge.
    % 4. Check if ray is on the same side of the edge as the opposing vertex.
    % 5. If this is true for all three edges, then the ray definitely passes through the facet.
    %
    % NOTES:
    % 1. If the ray crosses exactly on an edge, this is counted as crossing the facet.
    % 2. If a ray crosses exactly on a vertex, this is also taken into account.
    
    for loopCHECKFACET = possibleCROSSLIST'
      
      %Check if ray crosses the facet.  This method is much (>>10 times) faster than using the built-in function 'inpolygon'.
      %Taking each edge of the facet in turn, check if the ray is on the same side as the opposing vertex.  If so, let testVn=1
      
      Y1predicted = meshXYZ(loopCHECKFACET,2,2) - ((meshXYZ(loopCHECKFACET,2,2)-meshXYZ(loopCHECKFACET,2,3)) * (meshXYZ(loopCHECKFACET,1,2)-meshXYZ(loopCHECKFACET,1,1))/(meshXYZ(loopCHECKFACET,1,2)-meshXYZ(loopCHECKFACET,1,3)));
      YRpredicted = meshXYZ(loopCHECKFACET,2,2) - ((meshXYZ(loopCHECKFACET,2,2)-meshXYZ(loopCHECKFACET,2,3)) * (meshXYZ(loopCHECKFACET,1,2)-testx(loop))/(meshXYZ(loopCHECKFACET,1,2)-meshXYZ(loopCHECKFACET,1,3)));
      
      if (Y1predicted > meshXYZ(loopCHECKFACET,2,1) && YRpredicted > testy(loop)) || (Y1predicted < meshXYZ(loopCHECKFACET,2,1) && YRpredicted < testy(loop)) || (meshXYZ(loopCHECKFACET,2,2)-meshXYZ(loopCHECKFACET,2,3)) * (meshXYZ(loopCHECKFACET,1,2)-testx(loop)) == 0
        %                testV1 = 1;   %The ray is on the same side of the 2-3 edge as the 1st vertex.
      else
        %                testV1 = 0;   %The ray is on the opposite side of the 2-3 edge to the 1st vertex.
        % As the check is for ALL three checks to be true, we can continue here, if only one check fails
        continue;
      end %if
      
      Y2predicted = meshXYZ(loopCHECKFACET,2,3) - ((meshXYZ(loopCHECKFACET,2,3)-meshXYZ(loopCHECKFACET,2,1)) * (meshXYZ(loopCHECKFACET,1,3)-meshXYZ(loopCHECKFACET,1,2))/(meshXYZ(loopCHECKFACET,1,3)-meshXYZ(loopCHECKFACET,1,1)));
      YRpredicted = meshXYZ(loopCHECKFACET,2,3) - ((meshXYZ(loopCHECKFACET,2,3)-meshXYZ(loopCHECKFACET,2,1)) * (meshXYZ(loopCHECKFACET,1,3)-testx(loop))/(meshXYZ(loopCHECKFACET,1,3)-meshXYZ(loopCHECKFACET,1,1)));
      if (Y2predicted > meshXYZ(loopCHECKFACET,2,2) && YRpredicted > testy(loop)) || (Y2predicted < meshXYZ(loopCHECKFACET,2,2) && YRpredicted < testy(loop)) || (meshXYZ(loopCHECKFACET,2,3)-meshXYZ(loopCHECKFACET,2,1)) * (meshXYZ(loopCHECKFACET,1,3)-testx(loop)) == 0
        %                testV2 = 1;   %The ray is on the same side of the 3-1 edge as the 2nd vertex.
      else
        %                testV2 = 0;   %The ray is on the opposite side of the 3-1 edge to the 2nd vertex.
        % As the check is for ALL three checks to be true, we can continue here, if only one check fails
        continue;
      end %if
      
      Y3predicted = meshXYZ(loopCHECKFACET,2,1) - ((meshXYZ(loopCHECKFACET,2,1)-meshXYZ(loopCHECKFACET,2,2)) * (meshXYZ(loopCHECKFACET,1,1)-meshXYZ(loopCHECKFACET,1,3))/(meshXYZ(loopCHECKFACET,1,1)-meshXYZ(loopCHECKFACET,1,2)));
      YRpredicted = meshXYZ(loopCHECKFACET,2,1) - ((meshXYZ(loopCHECKFACET,2,1)-meshXYZ(loopCHECKFACET,2,2)) * (meshXYZ(loopCHECKFACET,1,1)-testx(loop))/(meshXYZ(loopCHECKFACET,1,1)-meshXYZ(loopCHECKFACET,1,2)));
      if (Y3predicted > meshXYZ(loopCHECKFACET,2,3) && YRpredicted > testy(loop)) || (Y3predicted < meshXYZ(loopCHECKFACET,2,3) && YRpredicted < testy(loop)) || (meshXYZ(loopCHECKFACET,2,1)-meshXYZ(loopCHECKFACET,2,2)) * (meshXYZ(loopCHECKFACET,1,1)-testx(loop)) == 0
        %                testV3 = 1;   %The ray is on the same side of the 1-2 edge as the 3rd vertex.
      else
        %                testV3 = 0;   %The ray is on the opposite side of the 1-2 edge to the 3rd vertex.
        % As the check is for ALL three checks to be true, we can continue here, if only one check fails
        continue;
      end %if
      
      nf=nf+1;facetCROSSLIST(nf)=loopCHECKFACET;
      
    end %for
    
    % Use only values ~=0
    facetCROSSLIST = facetCROSSLIST(1:nf);
    
    % - 3 - Find the z coordinate of the locations where the ray crosses each facet:
    gridCOzCROSS = zeros(1,nf);
    for loopFINDZ = facetCROSSLIST
      
      % METHOD:
      % 1. Define the equation describing the plane of the facet.  For a
      % more detailed outline of the maths, see:
      % http://local.wasp.uwa.edu.au/~pbourke/geometry/planeeq/
      %    Ax + By + Cz + D = 0
      %    where  A = y1 (z2 - z3) + y2 (z3 - z1) + y3 (z1 - z2)
      %           B = z1 (x2 - x3) + z2 (x3 - x1) + z3 (x1 - x2)
      %           C = x1 (y2 - y3) + x2 (y3 - y1) + x3 (y1 - y2)
      %           D = - x1 (y2 z3 - y3 z2) - x2 (y3 z1 - y1 z3) - x3 (y1 z2 - y2 z1)
      % 2. For the x and y coordinates of the ray, solve these equations to find the z coordinate in this plane.
      
      planecoA = meshXYZ(loopFINDZ,2,1)*(meshXYZ(loopFINDZ,3,2)-meshXYZ(loopFINDZ,3,3)) + meshXYZ(loopFINDZ,2,2)*(meshXYZ(loopFINDZ,3,3)-meshXYZ(loopFINDZ,3,1)) + meshXYZ(loopFINDZ,2,3)*(meshXYZ(loopFINDZ,3,1)-meshXYZ(loopFINDZ,3,2));
      planecoB = meshXYZ(loopFINDZ,3,1)*(meshXYZ(loopFINDZ,1,2)-meshXYZ(loopFINDZ,1,3)) + meshXYZ(loopFINDZ,3,2)*(meshXYZ(loopFINDZ,1,3)-meshXYZ(loopFINDZ,1,1)) + meshXYZ(loopFINDZ,3,3)*(meshXYZ(loopFINDZ,1,1)-meshXYZ(loopFINDZ,1,2));
      planecoC = meshXYZ(loopFINDZ,1,1)*(meshXYZ(loopFINDZ,2,2)-meshXYZ(loopFINDZ,2,3)) + meshXYZ(loopFINDZ,1,2)*(meshXYZ(loopFINDZ,2,3)-meshXYZ(loopFINDZ,2,1)) + meshXYZ(loopFINDZ,1,3)*(meshXYZ(loopFINDZ,2,1)-meshXYZ(loopFINDZ,2,2));
      planecoD = - meshXYZ(loopFINDZ,1,1)*(meshXYZ(loopFINDZ,2,2)*meshXYZ(loopFINDZ,3,3)-meshXYZ(loopFINDZ,2,3)*meshXYZ(loopFINDZ,3,2)) - meshXYZ(loopFINDZ,1,2)*(meshXYZ(loopFINDZ,2,3)*meshXYZ(loopFINDZ,3,1)-meshXYZ(loopFINDZ,2,1)*meshXYZ(loopFINDZ,3,3)) - meshXYZ(loopFINDZ,1,3)*(meshXYZ(loopFINDZ,2,1)*meshXYZ(loopFINDZ,3,2)-meshXYZ(loopFINDZ,2,2)*meshXYZ(loopFINDZ,3,1));
      
      if abs(planecoC) < 1e-14
        planecoC=0;
      end
      
      gridCOzCROSS(facetCROSSLIST==loopFINDZ) = (- planecoD - planecoA*testx(loop) - planecoB*testy(loop)) / planecoC;
      
    end %for
    
    if isempty(gridCOzCROSS),continue;end
    
    %Remove values of gridCOzCROSS which are outside of the mesh limits (including a 1e-12 margin for error).
    gridCOzCROSS = gridCOzCROSS( gridCOzCROSS>=meshZmin-1e-12 & gridCOzCROSS<=meshZmax+1e-12 );
    
    %Round gridCOzCROSS to remove any rounding errors, and take only the unique values:
    gridCOzCROSS = round(gridCOzCROSS*1e10)/1e10;
    
    % Replacement of the call to unique (gridCOzCROSS = unique(gridCOzCROSS);) by the following line:
    tmp = sort(gridCOzCROSS);I=[0,tmp(2:end)-tmp(1:end-1)]~=0;gridCOzCROSS = [tmp(1),tmp(I)];
    
    % - 4 - Label as being inside the mesh all the voxels that the ray passes through after crossing one facet before crossing another facet:
    
    if rem(numel(gridCOzCROSS),2)==0  % Only rays which cross an even number of facets are voxelised
      
      for loopASSIGN = 1:(numel(gridCOzCROSS)/2)
        voxelsINSIDE = (testz(loop)>gridCOzCROSS(2*loopASSIGN-1) & testz(loop)<gridCOzCROSS(2*loopASSIGN));
        OUTPUT(loop) = voxelsINSIDE;
        if voxelsINSIDE,break;end
      end %for
      
      
    elseif numel(gridCOzCROSS)~=0    % Remaining rays which meet the mesh in some way are not voxelised, but are labelled for correction later.
      correctionLIST = [ correctionLIST; loop ];
    end %if
    
  end %if
  
end %for

%======================================================
% RESTORE DIVIDE-BY-ZERO WARNINGS TO THE ORIGINAL STATE
%======================================================

warning(warningrestorestate)

% J.Korsawe: A correction is not possible as the testpoints need not to be
%            ordered in any way.
%            voxelise contains a correction algorithm which is appended here
%            without changes in syntax.
return

%======================================================
% USE INTERPOLATION TO FILL IN THE RAYS WHICH COULD NOT BE VOXELISED
%======================================================
%For rays where the voxelisation did not give a clear result, the ray is
%computed by interpolating from the surrounding rays.
countCORRECTIONLIST = size(correctionLIST,1);

if countCORRECTIONLIST>0
  
  %If necessary, add a one-pixel border around the x and y edges of the
  %array.  This prevents an error if the code tries to interpolate a ray at
  %the edge of the x,y grid.
  if min(correctionLIST(:,1))==1 || max(correctionLIST(:,1))==numel(gridCOx) || min(correctionLIST(:,2))==1 || max(correctionLIST(:,2))==numel(gridCOy)
    gridOUTPUT     = [zeros(1,voxcountY+2,voxcountZ);zeros(voxcountX,1,voxcountZ),gridOUTPUT,zeros(voxcountX,1,voxcountZ);zeros(1,voxcountY+2,voxcountZ)];
    correctionLIST = correctionLIST + 1;
  end
  
  for loopC = 1:countCORRECTIONLIST
    voxelsforcorrection = squeeze( sum( [ gridOUTPUT(correctionLIST(loopC,1)-1,correctionLIST(loopC,2)-1,:) ,...
      gridOUTPUT(correctionLIST(loopC,1)-1,correctionLIST(loopC,2),:)   ,...
      gridOUTPUT(correctionLIST(loopC,1)-1,correctionLIST(loopC,2)+1,:) ,...
      gridOUTPUT(correctionLIST(loopC,1),correctionLIST(loopC,2)-1,:)   ,...
      gridOUTPUT(correctionLIST(loopC,1),correctionLIST(loopC,2)+1,:)   ,...
      gridOUTPUT(correctionLIST(loopC,1)+1,correctionLIST(loopC,2)-1,:) ,...
      gridOUTPUT(correctionLIST(loopC,1)+1,correctionLIST(loopC,2),:)   ,...
      gridOUTPUT(correctionLIST(loopC,1)+1,correctionLIST(loopC,2)+1,:) ,...
      ] ) );
    voxelsforcorrection = (voxelsforcorrection>=4);
    gridOUTPUT(correctionLIST(loopC,1),correctionLIST(loopC,2),voxelsforcorrection) = 1;
  end %for
  
  %Remove the one-pixel border surrounding the array, if this was added
  %previously.
  if size(gridOUTPUT,1)>numel(gridCOx) || size(gridOUTPUT,2)>numel(gridCOy)
    gridOUTPUT = gridOUTPUT(2:end-1,2:end-1,:);
  end
  
end %if

%disp([' Ray tracing result: ',num2str(countCORRECTIONLIST),' rays (',num2str(countCORRECTIONLIST/(voxcountX*voxcountY)*100,'%5.1f'),'% of all rays) exactly crossed a facet edge and had to be computed by interpolation.'])

%==========================================================================

function D = rotmatrix(v,deg)
% calculate the rotation matrix about v by deg degrees

deg=deg/180*pi;
if deg~=0,
  v=v/norm(v);
  v1=v(1);v2=v(2);v3=v(3);ca=cos(deg);sa=sin(deg);
  D=[ca+v1*v1*(1-ca),v1*v2*(1-ca)-v3*sa,v1*v3*(1-ca)+v2*sa;
    v2*v1*(1-ca)+v3*sa,ca+v2*v2*(1-ca),v2*v3*(1-ca)-v1*sa;
    v3*v1*(1-ca)-v2*sa,v3*v2*(1-ca)+v1*sa,ca+v3*v3*(1-ca)];
else,
  D=eye(3,3);
end

