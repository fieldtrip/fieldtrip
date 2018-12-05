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
%   cfg.smooth_steps = number of standard smoothing iterations (default: 60)
%   cfg.fixshrinkage = reduce possible shrinkage due to standard smoothing
%                      (default: 'no')
%   cfg.expansion_mm = amount in mm used to globally expand hull, applies
%                      when cfg.fixshrinkage = 'yes' (default: 'auto')
%   cfg.laplace_steps = number of alternative (non-shrinking) smoothing
%                      iterations (default: 0, recommended: 1000 when used alone)
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
smooth_steps         = ft_getopt(cfg, 'smooth_steps', 60);
fixshrinkage         = ft_getopt(cfg, 'fixshrinkage', 'no');
expansion_mm         = ft_getopt(cfg, 'expansion_mm', 'auto'); % applies when fixshrinkage is 'yes'
laplace_steps        = ft_getopt(cfg, 'laplace_steps', 0);

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

% fill the surface
cmd = sprintf('mris_fill -c -r %d %s %s', resolution, surf, surf_filled);
system(['source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' cmd]);

% make outer surface
make_outer_surface(surf_filled, outer_surface_sphere, surf_outer)

% smooth the surface using mris_smooth
cmd = sprintf('mris_smooth -nw -n %d %s %s', smooth_steps, surf_outer, surf_smooth);
system(['source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' cmd]);
headshape = ft_read_headshape(surf_smooth);

% fix shrinkage if needed
if strcmp(fixshrinkage, 'yes')
  pial = ft_read_headshape(cfg.headshape);
  inside = inpolyhedron(pial.tri, pial.pos, headshape.pos); % determine inside/outside points
  if numel(find(inside==0)) > numel(find(inside==1)) % assuming the majority of points are inside
    inside = ~inside;
  end
  if any(inside==0)
    fprintf('expanding the hull to include more outside points\n');
    surf_expanded = [tempname() '_pial_expand'];
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
    cmd = sprintf('mris_expand %s %d %s', surf_smooth, expansion_mm, surf_expanded); % global expansion
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
end

% smooth using iso2mesh (non-shrinking)
if laplace_steps >= 1
  ft_hastoolbox('iso2mesh',1);
  fprintf('non-shrinking smoothing for %d iterations\n', laplace_steps)
  conn = meshconn(headshape.tri, size(headshape.pos,1)); % determine neighbors
  headshape.pos = smoothsurf(headshape.pos, [], conn, laplace_steps, 0, 'laplacianhc', .2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IN = inpolyhedron(varargin)
%INPOLYHEDRON  Tests if points are inside a 3D triangulated (faces/vertices) surface
%
%   IN = INPOLYHEDRON(FV,QPTS) tests if the query points (QPTS) are inside
%   the patch/surface/polyhedron defined by FV (a structure with fields
%   'vertices' and 'faces'). QPTS is an N-by-3 set of XYZ coordinates. IN
%   is an N-by-1 logical vector which will be TRUE for each query point
%   inside the surface. By convention, surface normals point OUT from the
%   object (see FLIPNORMALS option below if to reverse this convention).
%
%   INPOLYHEDRON(FACES,VERTICES,...) takes faces/vertices separately, rather than in
%   an FV structure.
%
%   IN = INPOLYHEDRON(..., X, Y, Z) voxelises a mask of 3D gridded query points
%   rather than an N-by-3 array of points. X, Y, and Z coordinates of the grid
%   supplied in XVEC, YVEC, and ZVEC respectively. IN will return as a 3D logical
%   volume with SIZE(IN) = [LENGTH(YVEC) LENGTH(XVEC) LENGTH(ZVEC)], equivalent to
%   syntax used by MESHGRID. INPOLYHEDRON handles this input faster and with a lower
%   memory footprint than using MESHGRID to make full X, Y, Z query points matrices.
%
%   INPOLYHEDRON(...,'PropertyName',VALUE,'PropertyName',VALUE,...) tests query
%   points using the following optional property values:
%
%   TOL           - Tolerance on the tests for "inside" the surface. You can think of
%   tol as the distance a point may possibly lie above/below the surface, and still
%   be perceived as on the surface. Due to numerical rounding nothing can ever be
%   done exactly here. Defaults to ZERO. Note that in the current implementation TOL
%   only affects points lying above/below a surface triangle (in the Z-direction).
%   Points coincident with a vertex in the XY plane are considered INside the surface.
%   More formal rules can be implemented with input/feedback from users.
%
%   GRIDSIZE      - Internally, INPOLYHEDRON uses a divide-and-conquer algorithm to
%   split all faces into a chessboard-like grid of GRIDSIZE-by-GRIDSIZE regions.
%   Performance will be a tradeoff between a small GRIDSIZE (few iterations, more
%   data per iteration) and a large GRIDSIZE (many iterations of small data
%   calculations). The sweet-spot has been experimentally determined (on a win64
%   system) to be correlated with the number of faces/vertices. You can overwrite
%   this automatically computed choice by specifying a GRIDSIZE parameter.
%
%   FACENORMALS   - By default, the normals to the FACE triangles are computed as the
%   cross-product of the first two triangle edges. You may optionally specify face
%   normals here if they have been pre-computed.
%
%   FLIPNORMALS   - (Defaults FALSE). To match a wider convention, triangle
%   face normals are presumed to point OUT from the object's surface. If
%   your surface normals are defined pointing IN, then you should set the
%   FLIPNORMALS option to TRUE to use the reverse of this convention.
%
%   Example:
%       tmpvol = zeros(20,20,20);       % Empty voxel volume
%       tmpvol(5:15,8:12,8:12) = 1;     % Turn some voxels on
%       tmpvol(8:12,5:15,8:12) = 1;
%       tmpvol(8:12,8:12,5:15) = 1;
%       fv = isosurface(tmpvol, 0.99);  % Create the patch object
%       fv.faces = fliplr(fv.faces);    % Ensure normals point OUT
%       % Test SCATTERED query points
%       pts = rand(200,3)*12 + 4;       % Make some query points
%       in = inpolyhedron(fv, pts);     % Test which are inside the patch
%       figure, hold on, view(3)        % Display the result
%       patch(fv,'FaceColor','g','FaceAlpha',0.2)
%       plot3(pts(in,1),pts(in,2),pts(in,3),'bo','MarkerFaceColor','b')
%       plot3(pts(~in,1),pts(~in,2),pts(~in,3),'ro'), axis image
%       % Test STRUCTURED GRID of query points
%       gridLocs = 3:2.1:19;
%       [x,y,z] = meshgrid(gridLocs,gridLocs,gridLocs);
%       in = inpolyhedron(fv, gridLocs,gridLocs,gridLocs);
%       figure, hold on, view(3)        % Display the result
%       patch(fv,'FaceColor','g','FaceAlpha',0.2)
%       plot3(x(in), y(in), z(in),'bo','MarkerFaceColor','b')
%       plot3(x(~in),y(~in),z(~in),'ro'), axis image
%
% See also: UNIFYMESHNORMALS (on the <a href="http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=43013">file exchange</a>)

% TODO-list
% - Optmise overall memory footprint. (need examples with MEM errors)
% - Implement an "ignore these" step to speed up calculations for:
%     * Query points outside the convex hull of the faces/vertices input
% - Get a better/best gridSize calculation. User feedback?
% - Detect cases where X-rays or Y-rays would be better than Z-rays?

%
%   Author: Sven Holcombe
% - 10 Jun 2012: Version 1.0
% - 28 Aug 2012: Version 1.1 - Speedup using accumarray
% - 07 Nov 2012: Version 2.0 - BEHAVIOUR CHANGE
%                Query points coincident with a VERTEX are now IN an XY triangle
% - 18 Aug 2013: Version 2.1 - Gridded query point handling with low memory footprint.
% - 10 Sep 2013: Version 3.0 - BEHAVIOUR CHANGE
%                NEW CONVENTION ADOPTED to expect face normals pointing IN
%                Vertically oriented faces are now ignored. Speeds up
%                computation and fixes bug where presence of vertical faces
%                produced NaN distance from a query pt to facet, making all
%                query points under facet erroneously NOT IN polyhedron.
% - 25 Sep 2013: Version 3.1 - Dropped nested unique call which was made
%                mostly redundant via v2.1 gridded point handling. Also
%                refreshed grid size selection via optimisation.
% - 25 Feb 2014: Version 3.2 - Fixed indeterminate behaviour for query
%                points *exactly* in line with an "overhanging" vertex.
% - 11 Nov 2016: Version 3.3 - Used quoted semicolons ':' inside function
%                handle calls to conform with new 2015b interpreter
%%

% FACETS is an unpacked arrangement of faces/vertices. It is [3-by-3-by-N],
% with 3 1-by-3 XYZ coordinates of N faces.
[facets, qPts, options] = parseInputs(varargin{:});
numFaces = size(facets,3);
if ~options.griddedInput            % SCATTERED QUERY POINTS
  numQPoints = size(qPts,1);
else                                % STRUCTURED QUERY POINTS
  numQPoints = prod(cellfun(@numel,qPts(1:2)));
end

% Precompute 3d normals to all facets (triangles). Do this via the cross
% product of the first edge vector with the second. Normalise the result.
allEdgeVecs = facets([2 3 1],:,:) - facets(:,:,:);
if isempty(options.facenormals)
  allFacetNormals =  bsxfun(@times, allEdgeVecs(1,[2 3 1],:), allEdgeVecs(2,[3 1 2],:)) - ...
    bsxfun(@times, allEdgeVecs(2,[2 3 1],:), allEdgeVecs(1,[3 1 2],:));
  allFacetNormals = bsxfun(@rdivide, allFacetNormals, sqrt(sum(allFacetNormals.^2,2)));
else
  allFacetNormals = permute(options.facenormals,[3 2 1]);
end
if options.flipnormals
  allFacetNormals = -allFacetNormals;
end
% We use a Z-ray intersection so we don't even need to consider facets that
% are purely vertically oriented (have zero Z-component).
isFacetUseful = allFacetNormals(:,3,:) ~= 0;

%% Setup grid referencing system
% Function speed can be thought of as a function of grid size. A small number of grid
% squares means iterating over fewer regions (good) but with more faces/qPts to
% consider each time (bad). For any given mesh/queryPt configuration, there will be a
% sweet spot that minimises computation time. There will also be a constraint from
% memory available - low grid sizes means considering many queryPt/faces at once,
% which will require a larger memory footprint. Here we will let the user specify
% gridsize directly, or we will estimate the optimum size based on prior testing.
if ~isempty(options.gridsize)
  gridSize = options.gridsize;
else
  % Coefficients (with 95% confidence bounds):
  p00 =         -47;    p10 =       12.83;    p01 =       20.89;
  p20 =      0.7578;    p11 =      -6.511;    p02 =      -2.586;
  p30 =     -0.1802;    p21 =      0.2085;    p12 =      0.7521;
  p03 =     0.09984;    p40 =    0.005815;    p31 =    0.007775;
  p22 =    -0.02129;    p13 =    -0.02309;
  GSfit = @(x,y)p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y + p12*x*y^2 + p03*y^3 + p40*x^4 + p31*x^3*y + p22*x^2*y^2 + p13*x*y^3;
  gridSize = min(150 ,max(1, ceil(GSfit(log(numQPoints),log(numFaces)))));
  if isnan(gridSize), gridSize = 1; end
end

%% Find candidate qPts -> triangles pairs
% We have a large set of query points. For each query point, find potential
% triangles that would be pierced by vertical rays through the qPt. First,
% a simple filter by XY bounding box

% Calculate the bounding box of each facet
minFacetCoords = permute(min(facets(:,1:2,:),[],1),[3 2 1]);
maxFacetCoords = permute(max(facets(:,1:2,:),[],1),[3 2 1]);

% Set rescale values to rescale all vertices between 0(-eps) and 1(+eps)
scalingOffsetsXY = min(minFacetCoords,[],1) - eps;
scalingRangeXY = max(maxFacetCoords,[],1) - scalingOffsetsXY + 2*eps;

% Based on scaled min/max facet coords, get the [lowX lowY highX highY] "grid" index
% of all faces
lowToHighGridIdxs = floor(bsxfun(@rdivide, ...
  bsxfun(@minus, ... % Use min/max coordinates of each facet (+/- the tolerance)
  [minFacetCoords-options.tol maxFacetCoords+options.tol],...
  [scalingOffsetsXY scalingOffsetsXY]),...
  [scalingRangeXY scalingRangeXY]) * gridSize) + 1;

% Build a grid of cells. In each cell, place the facet indices that encroach into
% that grid region. Similarly, each query point will be assigned to a grid region.
% Note that query points will be assigned only one grid region, facets can cover many
% regions. Furthermore, we will add a tolerance to facet region assignment to ensure
% a query point will be compared to facets even if it falls only on the edge of a
% facet's bounding box, rather than inside it.
cells = cell(gridSize);
[unqLHgrids,~,facetInds] = unique(lowToHighGridIdxs,'rows');
tmpInds = accumarray(facetInds(isFacetUseful),find(isFacetUseful),[size(unqLHgrids,1),1],@(x){x});
for xi = 1:gridSize
  xyMinMask = xi >= unqLHgrids(:,1) & xi <= unqLHgrids(:,3);
  for yi = 1:gridSize
    cells{yi,xi} = cat(1,tmpInds{xyMinMask & yi >= unqLHgrids(:,2) & yi <= unqLHgrids(:,4)});
    % The above line (with accumarray) is faster with equiv results than:
    % % cells{yi,xi} = find(ismember(facetInds, xyInds));
  end
end
% With large number of facets, memory may be important:
clear lowToHightGridIdxs LHgrids facetInds tmpInds xyMinMask minFacetCoords maxFacetCoords

%% Compute edge unit vectors and dot products

% Precompute the 2d unit vectors making up each facet's edges in the XY plane.
allEdgeUVecs = bsxfun(@rdivide, allEdgeVecs(:,1:2,:), sqrt(sum(allEdgeVecs(:,1:2,:).^2,2)));

% Precompute the inner product between edgeA.edgeC, edgeB.edgeA, edgeC.edgeB
allEdgeEdgeDotPs = sum(allEdgeUVecs .* -allEdgeUVecs([3 1 2],:,:),2) - 1e-9;

%% Gather XY query locations
% Since query points are most likely given as a (3D) grid of query locations, we only
% need to consider the unique XY locations when asking which facets a vertical ray
% through an XY location would pierce.
if ~options.griddedInput            % SCATTERED QUERY POINTS
  qPtsXY = @(varargin)qPts(:,1:2);
  qPtsXYZViaUnqIndice = @(ind)qPts(ind,:);
  outPxIndsViaUnqIndiceMask = @(ind,mask)ind(mask);
  outputSize = [size(qPts,1),1];
  reshapeINfcn = @(INMASK)INMASK;
  minFacetDistanceFcn = @minFacetToQptDistance;
else                                % STRUCTURED QUERY POINTS
  [xmat,ymat] = meshgrid(qPts{1:2});
  qPtsXY = [xmat(:) ymat(:)];
  % A standard set of Z locations will be shifted around by different
  % unqQpts XY coordinates.
  zCoords = qPts{3}(:) * [0 0 1];
  qPtsXYZViaUnqIndice = @(ind)bsxfun(@plus, zCoords, [qPtsXY(ind,:) 0]);
  % From a given indice and mask, we will turn on/off the IN points under
  % that indice based on the mask. The easiest calculation is to setup
  % the IN matrix as a numZpts-by-numUnqPts mask. At the end, we must
  % unpack/reshape this 2D mask to a full 3D logical mask
  numZpts = size(zCoords,1);
  baseZinds = 1:numZpts;
  outPxIndsViaUnqIndiceMask = @(ind,mask)(ind-1)*numZpts + baseZinds(mask);
  outputSize = [numZpts, size(qPtsXY,1)];
  reshapeINfcn = @(INMASK)reshape(INMASK', cellfun(@numel, qPts([2 1 3])));
  minFacetDistanceFcn = @minFacetToQptsDistance;
end

% Start with every query point NOT inside the polyhedron. We will
% iteratively find those query points that ARE inside.
IN = false(outputSize);
% Determine with grids each query point falls into.
qPtGridXY = floor(bsxfun(@rdivide, bsxfun(@minus, qPtsXY(':',':'), scalingOffsetsXY),...
  scalingRangeXY) * gridSize) + 1;
[unqQgridXY,~,qPtGridInds] = unique(qPtGridXY,'rows');
% We need only consider grid indices within those already set up
ptsToConsidMask = ~any(qPtGridXY<1 | qPtGridXY>gridSize, 2);
if ~any(ptsToConsidMask)
  IN = reshapeINfcn(IN);
  return;
end
% Build the reference list
cellQptContents = accumarray(qPtGridInds(ptsToConsidMask),find(ptsToConsidMask), [],@(x){x});
gridsToCheck = unqQgridXY(~any(unqQgridXY<1 | unqQgridXY>gridSize, 2),:);
cellQptContents(cellfun('isempty',cellQptContents)) = [];
gridIndsToCheck = sub2ind(size(cells), gridsToCheck(:,2), gridsToCheck(:,1));

% For ease of multiplication, reshape qPt XY coords to [1-by-2-by-1-by-N]
qPtsXY = permute(qPtsXY(':',':'),[4 2 3 1]);

% There will be some grid indices with query points but without facets.
emptyMask = cellfun('isempty',cells(gridIndsToCheck))';
for i = find(~emptyMask)
  % We get all the facet coordinates (ie, triangle vertices) of triangles
  % that intrude into this grid location. The size is [3-by-2-by-N], for
  % the [3vertices-by-XY-by-Ntriangles]
  allFacetInds = cells{gridIndsToCheck(i)};
  candVerts = facets(:,1:2,allFacetInds);
  % We need the XY coordinates of query points falling into this grid.
  allqPtInds = cellQptContents{i};
  queryPtsXY = qPtsXY(:,:,:,allqPtInds);
  
  % Get unit vectors pointing from each triangle vertex to my query point(s)
  vert2ptVecs = bsxfun(@minus, queryPtsXY, candVerts);
  vert2ptUVecs = bsxfun(@rdivide, vert2ptVecs, sqrt(sum(vert2ptVecs.^2,2)));
  % Get unit vectors pointing around each triangle (along edge A, edge B, edge C)
  edgeUVecs = allEdgeUVecs(:,:,allFacetInds);
  % Get the inner product between edgeA.edgeC, edgeB.edgeA, edgeC.edgeB
  edgeEdgeDotPs = allEdgeEdgeDotPs(:,:,allFacetInds);
  % Get inner products between each edge unit vec and the UVs from qPt to vertex
  edgeQPntDotPs = sum(bsxfun(@times, edgeUVecs, vert2ptUVecs),2);
  qPntEdgeDotPs = sum(bsxfun(@times,vert2ptUVecs, -edgeUVecs([3 1 2],:,:)),2);
  % If both inner products 2 edges to the query point are greater than the inner
  % product between the two edges themselves, the query point is between the V
  % shape made by the two edges. If this is true for all 3 edge pair, the query
  % point is inside the triangle.
  resultIN = all(bsxfun(@gt, edgeQPntDotPs, edgeEdgeDotPs) & bsxfun(@gt, qPntEdgeDotPs, edgeEdgeDotPs),1);
  resultONVERTEX = any(any(isnan(vert2ptUVecs),2),1);
  result = resultIN | resultONVERTEX;
  qPtHitsTriangles = any(result,3);
  % If NONE of the query points pierce ANY triangles, we can skip forward
  if ~any(qPtHitsTriangles), continue, end
  
  % In the next step, we'll need to know the indices of ALL the query points at
  % each of the distinct XY coordinates. Let's get their indices into "qPts" as a
  % cell of length M, where M is the number of unique XY points we had found.
  for ptNo = find(qPtHitsTriangles(:))'
    % Which facets does it pierce?
    piercedFacetInds = allFacetInds(result(1,1,:,ptNo));
    
    % Get the 1-by-3-by-N set of triangle normals that this qPt pierces
    piercedTriNorms = allFacetNormals(:,:,piercedFacetInds);
    
    % Pick the first vertex as the "origin" of a plane through the facet. Get the
    % vectors from each query point to each facet origin
    facetToQptVectors = bsxfun(@minus, ...
      qPtsXYZViaUnqIndice(allqPtInds(ptNo)),...
      facets(1,:,piercedFacetInds));
    
    % Calculate how far you need to go up/down to pierce the facet's plane.
    % Positive direction means "inside" the facet, negative direction means
    % outside.
    facetToQptDists = bsxfun(@rdivide, ...
      sum(bsxfun(@times,piercedTriNorms,facetToQptVectors),2), ...
      abs(piercedTriNorms(:,3,:)));
    
    % Since it's possible for two triangles sharing the same vertex to
    % be the same distance away, I want to sum up all the distances of
    % triangles that are closest to the query point. Simple case: The
    % closest triangle is unique Edge case: The closest triangle is one
    % of many the same distance and direction away. Tricky case: The
    % closes triangle has another triangle the equivalent distance
    % but facing the opposite direction
    IN( outPxIndsViaUnqIndiceMask(allqPtInds(ptNo), ...
      minFacetDistanceFcn(facetToQptDists)<options.tol...
      )) = true;
  end
end

% If they provided X,Y,Z vectors of query points, our output is currently a
% 2D mask and must be reshaped to [LEN(Y) LEN(X) LEN(Z)].
IN = reshapeINfcn(IN);

%% Called subfunctions

% vertices = [
%     0.9046    0.1355   -0.0900
%     0.8999    0.3836   -0.0914
%     1.0572    0.2964   -0.0907
%     0.8735    0.1423   -0.1166
%     0.8685    0.4027   -0.1180
%     1.0337    0.3112   -0.1173
%     0.9358    0.1287   -0.0634
%     0.9313    0.3644   -0.0647
%     1.0808    0.2816   -0.0641
% ];
% faces = [
%      1     2     5
%      1     5     4
%      2     3     6
%      2     6     5
%      3     1     4
%      3     4     6
%      6     4     5
%      2     1     8
%      8     1     7
%      3     2     9
%      9     2     8
%      1     3     7
%      7     3     9
%      7     9     8
% ];
% point = [vertices(3,1),vertices(3,2),1.5];


function closestTriDistance = minFacetToQptDistance(facetToQptDists)
% FacetToQptDists is a 1pt-by-1-by-Nfacets array of how far you need to go
% up/down to pierce each facet's plane. If the Qpt was directly over an
% "overhang" vertex, then two facets with opposite orientation will be
% equally distant from the Qpt, with one distance positive and one
% negative. In such cases, it is impossible for the Qpt to actually be
% "inside" this pair of facets, so their distance is updated to Inf.

[~,minInd] = min(abs(facetToQptDists),[],3);
while any( abs(facetToQptDists + facetToQptDists(minInd)) < 1e-15 )
  % Since the above comparison is made every time, but the below variable
  % setting is done only in the rare case that a query point coincides
  % with an overhang vertex, it is more efficient to re-compute the
  % equality when it's true, rather than store the result every time.
  facetToQptDists( abs(facetToQptDists) - abs(facetToQptDists(minInd)) < 1e-15) = inf;
  if ~any(isfinite(facetToQptDists))
    break;
  end
  [~,minInd] = min(abs(facetToQptDists),[],3);
end
closestTriDistance = facetToQptDists(minInd);

function closestTriDistance = minFacetToQptsDistance(facetToQptDists)
% As above, but facetToQptDists is an Mpts-by-1-by-Nfacets array.

% The multi-point version is a little more tricky. While below is quite a
% bit slower when the while loop is entered, it is very rarely entered and
% very fast to make just the initial comparison.
[minVals,minInds] = min(abs(facetToQptDists),[],3);
while any(...
    any(abs(bsxfun(@plus,minVals,facetToQptDists))<1e-15,3) & ...
    any(abs(bsxfun(@minus,minVals,facetToQptDists))<1e-15,3))
  maskP = abs(bsxfun(@plus,minVals,facetToQptDists))<1e-15;
  maskN = abs(bsxfun(@minus,minVals,facetToQptDists))<1e-15;
  mustAlterMask = any(maskP,3) & any(maskN,3);
  for i = find(mustAlterMask)'
    facetToQptDists(i,:,maskP(i,:,:) | maskN(i,:,:)) = inf;
  end
  [newMv,newMinInds] = min(abs(facetToQptDists(mustAlterMask,:,:)),[],3);
  minInds(mustAlterMask) = newMinInds(:);
  minVals(mustAlterMask) = newMv(:);
end
% Below is a tiny speedup on basically a sub2ind call.
closestTriDistance = facetToQptDists((minInds-1)*size(facetToQptDists,1) + (1:size(facetToQptDists,1))');

%% Input handling subfunctions
function [facets, qPts, options] = parseInputs(varargin)

% Gather FACES and VERTICES
if isstruct(varargin{1})                        % inpolyhedron(FVstruct, ...)
  if ~all(isfield(varargin{1},{'vertices','faces'}))
    error( 'Structure FV must have "faces" and "vertices" fields' );
  end
  faces = varargin{1}.faces;
  vertices = varargin{1}.vertices;
  varargin(1) = []; % Chomp off the faces/vertices
  
else                                            % inpolyhedron(FACES, VERTICES, ...)
  faces = varargin{1};
  vertices = varargin{2};
  varargin(1:2) = []; % Chomp off the faces/vertices
end

% Unpack the faces/vertices into [3-by-3-by-N] facets. It's better to
% perform this now and have FACETS only in memory in the main program,
% rather than FACETS, FACES and VERTICES
facets = vertices';
facets = permute(reshape(facets(:,faces'), 3, 3, []),[2 1 3]);

% Extract query points
if length(varargin)<2 || ischar(varargin{2})    % inpolyhedron(F, V, [x(:) y(:) z(:)], ...)
  qPts = varargin{1};
  varargin(1) = []; % Chomp off the query points
else                                            % inpolyhedron(F, V, xVec, yVec, zVec, ...)
  qPts = varargin(1:3);
  % Chomp off the query points and tell the world that it's gridded input.
  varargin(1:3) = [];
  varargin = [varargin {'griddedInput',true}];
end

% Extract configurable options
options = parseOptions(varargin{:});

% Check if face normals are unified
if options.testNormals
  options.normalsAreUnified = checkNormalUnification(faces);
end

function options = parseOptions(varargin)
IP = inputParser;
if verLessThan('matlab', 'R2013b')
  fcn = 'addParamValue';
else
  fcn = 'addParameter';
end
IP.(fcn)('gridsize',[], @(x)isscalar(x) && isnumeric(x))
IP.(fcn)('tol', 0, @(x)isscalar(x) && isnumeric(x))
IP.(fcn)('tol_ang', 1e-5, @(x)isscalar(x) && isnumeric(x))
IP.(fcn)('facenormals',[]);
IP.(fcn)('flipnormals',false);
IP.(fcn)('griddedInput',false);
IP.(fcn)('testNormals',false);
IP.parse(varargin{:});
options = IP.Results;
