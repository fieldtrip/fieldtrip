function [h,Coord,Color,height] = som_dendrogram(Z,varargin)

%SOM_DENDROGRAM Visualize a dendrogram.
%
% [h,Coord,Color,height] = som_dendrogram(Z, [[argID,] value, ...])
%
%  Z = som_linkage(sM); 
%  som_dendrogram(Z); 
%  som_dendrogram(Z,sM); 
%  som_dendrogram(Z,'coord',co); 
%
%  Input and output arguments ([]'s are optional):
%   h        (vector) handle to the arc lines
%   Z        (matrix) size n-1 x 1, the hierarchical cluster matrix
%                     returned by functions like LINKAGE and SOM_LINKAGE
%                     n is the number of original data samples.
%   [argID,  (string) See below. The values which are unambiguous can 
%    value]  (varies) be given without the preceeding argID.
%   Coord    (matrix) size 2*n-1 x {1,2}, the coordinates of the
%                     original data samples and cluster nodes used 
%                     in the visualization
%   Color    (matrix) size 2*n-1 x 3, the colors of ...
%   height   (vector) size 2*n-1 x 1, the heights of ...
%   
% Here are the valid argument IDs and corresponding values. The values 
% which are unambiguous (marked with '*') can be given without the
% preceeding argID.
%   'data'  *(struct) map or data struct: many other optional 
%                     arguments require this 
%            (matrix) data matrix
%   'coord'  (matrix) size n x 1 or n x 2, the coordinates of 
%                     the original data samples either in 1D or 2D
%            (matrix) size 2*n-1 x {1,2}, the coordinates of both
%                     original data samples and each cluster
%           *(string) 'SOM', 'pca#', 'sammon#', or 'cca#': the coordinates
%                     are calculated using the given data and the 
%                     required projection algorithm. The '#' at the
%                     end of projection algorithms refers to the 
%                     desired output dimension and can be either 1 or 2
%                     (2 by default). In case of 'SOM', the unit
%                     coordinates (given by SOM_VIS_COORDS) are used.
%   'color'  (matrix) size n x 3, the color of the original data samples
%            (matrix) size 2*n-1 x 3, the colors of both original 
%                     data samples and each cluster
%            (string) color specification, e.g. 'r.', used for each node
%   'height' (vector) size n-1 x 1, the heights used for each cluster
%            (vector) size 2*n-1 x 1, the heights used for both original
%                     data samples and each cluster
%           *(string) 'order', the order of combination determines height
%                     'depth', the depth at which the combination
%                     happens determines height
%   'linecolor' (string) color specification for the arc color, 'k' by default
%               (vector) size 1 x 3
%
% See also SOM_LINKAGE, DENDROGRAM.

% Copyright (c) 2000 by Juha Vesanto
% Contributed to SOM Toolbox on June 16th, 2000 by Juha Vesanto
% http://www.cis.hut.fi/projects/somtoolbox/
 
% Version 2.0beta juuso 160600

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read the arguments

% Z 
nd = size(Z,1)+1; 
nc = size(Z,1); 

% varargin
Coordtype = 'natural'; Coord = []; codim = 1; 
Colortype = 'none'; Color = [];
height = [zeros(nd,1); Z(:,3)]; 
M = []; 
linecol = 'k'; 

i=1; 
while i<=length(varargin), 
  argok = 1; 
  if ischar(varargin{i}), 
    switch varargin{i}, 
     case 'data', i = i + 1; M = varargin{i}; 
     case 'coord', 
      i=i+1; 
      if isnumeric(varargin{i}), Coord = varargin{i}; Coordtype = 'given'; 
      else 
	if strcmp(varargin{i},'SOM'), Coordtype = 'SOM'; 
	else Coordtype = 'projection'; Coord = varargin{i}; 
	end
      end
     case 'color', 
      i=i+1; 
      if isempty(varargin{i}), Colortype = 'none'; 
      elseif ischar(varargin{i}), Colortype = 'colorspec'; Color = varargin{i};
      else Colortype = 'given'; Color = varargin{i}; 
      end
     case 'height', i=i+1; height = varargin{i}; 
     case 'linecolor', i=i+1; linecol = varargin{i}; 
     case 'SOM', 
      Coordtype = 'SOM'; 
     case {'pca','pca1','pca2','sammon','sammon1','sammon2','cca','cca1','cca2'}, 
      Coordtype = 'projection'; Coord = varargin{i}; 
     case {'order','depth'}, height = varargin{i};  
    end
  elseif isstruct(varargin{i}), M = varargin{i}; 
  else
    argok = 0; 
  end
  if ~argok, 
    disp(['(som_dendrogram) Ignoring invalid argument #' num2str(i+1)]); 
  end
  i = i+1; 
end

switch Coordtype, 
 case 'SOM', 
  if isempty(M) | ~any(strcmp(M.type,{'som_map','som_topol'})) , 
    error('Cannot determine SOM coordinates without a SOM.'); 
  end
  if strcmp(M.type,'som_map'), M = M.topol; end
 case 'projection', 
  if isempty(M), error('Cannot do projection without the data.'); end
  if isstruct(M), 
    if strcmp(M.type,'som_data'), M = M.data; 
    elseif strcmp(M.type,'som_map'), M = M.codebook; 
    end
  end  
  if size(M,1) ~= nd, 
    error('Given data must be equal in length to the number of original data samples.')
  end    
 case 'given', 
  if size(Coord,1) ~= nd & size(Coord,1) ~= nd+nc, 
    error('Size of given coordinate matrix does not match the cluster hierarchy.');
  end  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialization

% Coordinates
switch Coordtype, 
 case 'natural', o = leavesorder(Z)'; [dummy,Coord] = sort(o); codim = 1; 
 case 'SOM', Coord = som_vis_coords(M.lattice,M.msize); codim = 2;   
 case 'projection', 
  switch Coord, 
   case {'pca','pca2'},       Coord = pcaproj(M,2);   codim = 2; 
   case 'pca1',               Coord = pcaproj(M,1);   codim = 1; 
   case {'cca','cca2'},       Coord = cca(M,2,20);    codim = 2; 
   case 'cca1',               Coord = cca(M,1,20);    codim = 1; 
   case {'sammon','sammon2'}, Coord = sammon(M,2,50); codim = 2; 
   case 'sammon1',            Coord = sammon(M,1,50); codim = 1; 
  end
 case 'given', codim = min(size(Coord,2),2); % nill
end

if size(Coord,1) == nd, 
  Coord = [Coord; zeros(nc,size(Coord,2))]; 
  for i=(nd+1):(nd+nc), 
    leaves = leafnodes(Z,i,nd);
    if any(leaves), Coord(i,:) = mean(Coord(leaves,:),1); else Coord(i,:) = Inf; end
  end
end

% Colors
switch Colortype, 
 case 'colorspec', % nill
 case 'none', Color = ''; 
 case 'given',
  if size(Color,1) == nd, 
    Color = [Color; zeros(nc,3)]; 
    for i=(nd+1):(nd+nc), 
      leaves = leafnodes(Z,i,nd);
      if any(leaves), Color(i,:) = mean(Color(leaves,:),1); 
      else Color(i,:) = 0.8; 
      end
    end
  end
end

% height
if ischar(height), 
  switch height, 
   case 'order', height = [zeros(nd,1); [1:nc]']; 
   case 'depth', height = nodedepth(Z); height = max(height) - height; 
  end
else
  if length(height)==nc, height = [zeros(nd,1); height]; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% draw

% the arcs
lfrom = []; lto = []; 
for i=1:nd+nc, 
  if i<=nd, ch = [];  
  elseif ~isfinite(Z(i-nd,3)), ch = []; 
  else ch = Z(i-nd,1:2)'; 
  end
  if any(ch), 
    lfrom = [lfrom; i*ones(length(ch),1)]; 
    lto = [lto; ch]; 
  end
end

% the coordinates of the arcs
if codim == 1, 
  Lx = [Coord(lfrom), Coord(lto), Coord(lto)];
  Ly = [height(lfrom), height(lfrom), height(lto)];
  Lz = []; 
else
  Lx = [Coord(lfrom,1), Coord(lto,1), Coord(lto,1)];
  Ly = [Coord(lfrom,2), Coord(lto,2), Coord(lto,2)];
  Lz = [height(lfrom), height(lfrom), height(lto)];
end

washold = ishold; 
if ~washold, cla; end

% plot the lines
if isempty(Lz), 
  h = line(Lx',Ly','color',linecol); 
else 
  h = line(Lx',Ly',Lz','color',linecol); 
  if ~washold, view(3); end
  rotate3d on
end

% plot the nodes
hold on
switch Colortype, 
 case 'none', % nill
 case 'colorspec', 
  if codim == 1, plot(Coord,height,Color); 
  else plot3(Coord(:,1), Coord(:,2), height, Color); 
  end
 case 'given', 
  som_grid('rect',[nd+nc 1],'line','none','Coord',[Coord, height],...
	   'Markersize',10,'Markercolor',Color);
end
if ~washold, hold off, end

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% subfunctions

function depth = nodedepth(Z)

  nd = size(Z,1)+1; 
  nc = size(Z,1);  
  depth = zeros(nd+nc,1); 
  ch = nc+nd-1; 
  while any(ch), 
    c  = ch(1); ch = ch(2:end);
    if c>nd & isfinite(Z(c-nd,3)), 
      chc = Z(c-nd,1:2); 
      depth(chc) = depth(c) + 1; 
      ch = [ch, chc]; 
    end
  end
  return;

function inds = leafnodes(Z,i,nd)

  inds = []; 
  ch = i; 
  while any(ch), 
    c  = ch(1); ch = ch(2:end);
    if c>nd & isfinite(Z(c-nd,3)), ch = [ch, Z(c-nd,1:2)]; end
    if c<=nd, inds(end+1) = c; end 
  end
  return;

function order = leavesorder(Z)

  nd = size(Z,1)+1;
  order = 2*nd-1; 
  nonleaves = 1; 
  while any(nonleaves), 
    j = nonleaves(1); 
    ch = Z(order(j)-nd,1:2);
    if j==1, oleft = []; else oleft = order(1:(j-1)); end
    if j==length(order), oright = []; else oright = order((j+1):length(order)); end
    order = [oleft, ch, oright];
    nonleaves = find(order>nd); 
  end
  return; 


  