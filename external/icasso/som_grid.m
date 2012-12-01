function [S,m,l,t,s]=som_grid(varargin)

%SOM_GRID Visualization of a SOM grid
%
% [sGrid,m,l,t,s]=som_grid(sGrid, ['argID', value, ...])
% [sGrid,m,l,t,s]=som_grid(topol, ['argID', value, ...])
% [sGrid,m,l,t,s]=som_grid(lattice, msize, ['argID', value, ...])
%
% Input and output arguments ([]'s are optional)
%  sGrid    (struct) som_grid struct (see output arguments)
%  topol    (struct) map or topol struct for giving the topology
%           (cell array) of form {'lattice', msize, ['shape']}. 
%                    Default value for 'shape' is 'sheet'.
%  lattice  (string) 'hexa', 'rect' 
%           (matrix) size M x M, defines topological connections             
%  msize    (vector) 1x2 vector defines the grid size, M=msize(1)*msize(2)
%  ['argID',(string) Other arguments can be given as 'argID', value   
%   value]  (varies) pairs. See list below for valid values.
%
%  sGrid    (struct) with fields S.msize, S.shape, S.lattice, S.coord, S.marker, 
%                    S.markersize, S.markercolor, S.line, S.linewidth, S.linecolor,
%                    S.surf, S.label, S.labelsize, S.labelcolor
%  m        (matrix) handels to LINE objects (unit markers) 
%  l        (matrix) handles to LINE objects (lines connecting the units)
%  t        (matrix) handles to TEXT objects (labels)
%  s        (scalar) handle  to SURF object  (surface between units)
%
%  Here are the valid argument IDs (case insensitive) and
%  associated values: 
%  'Coord'       Mx2 or Mx3 matrix of coordinates 
%                (default: according to lattice as in som_cplane)
%  'Marker'      string 'o','+','x','*','v','^','<','>','h','s','d','p','.', 
%                'none' or Mx1 cell or char array of these strings 
%                Default: 'o'.  
%  'MarkerSize'  scalar or Mx1 matrix of double. Default: 6.
%  'MarkerColor' ColorSpec or Mx3 matrix of RGB triples. Default: 'k'.
%  'Line'        string '-',':','--' or '-.' or 'none'. Default: '-'.
%  'Surf'        [], Mx1 or Mx3 matrix of RGB triples 
%                to define surface values. Default: [] = no surf. 
%                Note: shading is turned to 'interp'.
%  'LineWidth'   scalar or MxM matrix, default: 0.5
%  'LineColor'   ColorSepc, MxMx3 matrix of RGB triples or a cell array 
%                of form {r g b} where r,g, and b are MxM  
%                (sparse) matrices of R,G, and B values
%  'Label'       Mx1 char array, cell array of strings size MxL 
%                or [] to indicate no labels, default: [] = no labels.
%  'LabelSize'   scalar
%  'LabelColor'  ColorSpec or string 'none', default: 'g'.
%
% For more help, try 'type som_grid' or check out online documentation.
% See also SOM_CONNECTION, SOM_SHOW, SOM_CPLANE, SOM_SET, SCATTER, SCATTER3.

%%%%%%%%%%%%%%% DETAILED DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% som_grid
%
% PURPOSE 
% 
%  To visualize the SOM grid in various ways 
% 
% SYNTAX
%
%  [sGrid,m,l,t,s]=som_grid(sGrid)
%  [sGrid,m,l,t,s]=som_grid(sTopol)
%  [sGrid,m,l,t,s]=som_grid(sMap)
%  [sGrid,m,l,t,s]=som_grid({lattice, msize, [shape]}) 
%  [sGrid,m,l,t,s]=som_grid(lattice, msize)
%  [sGrid,m,l,t,s]=som_grid(..., ['argID', value, ...])
%
% DESCRIPTION 
%
% The SOM can be defined as a set of units (neurons) and their
% topological relations.  This function is used to visualize these in
% various ways. The units may be drawn using different markers and
% colors, in different sizes and in different locations in 2D or
% 3D. However the topological neighborhood is limited to be
% 2-dimensional. The connections between these units may be drawn using
% lines having different thicknesses and colors. Labeling text may be
% plotted on the units. It is possible also to draw a surface between
% the units. The surface coloring is either indexed (one value per
% unit) or fixed RGB (a 1x3 RGB triple per unit).
%
% REQUIRED INPUT ARGUMENTS 
%
% Note: M is the number of map units.
% 
% The first (or first two) argument may have various different types of values
%
%   1. sGrid   (struct) som_grid struct (the output of this function)
%  
%     The struct initiates the visualization. The argID-value -pairs
%     are used to alter the initiation.
%
%   Following argument types may be used to give the topology for the grid
% 
%   2. sTopol  (struct) som_topol struct 
%   3. sMap    (struct) som_map struct (only topology matters)
%   4. {lattice, msize} or {lattice, msize, sheet} (cell array) 
%       - lattice must be 'hexa' or 'rect'
%       - msize must be a 1x2 vector
%       - shape (if specified) must be string 'sheet', 'cyl' or 'toroid'
%         If shape is not given it is 'sheet' by default.
%   5. lattice (string or matrix) AND msize (1x2 vector) as two separate arguments
%       - lattice may be string 'rect' or 'hexa' or a connection matrix
%         (see SOM_CONNECTION) to define a free topology. This connection
%         matrix is of size MxM and its element i,j (i<j) is set
%         to 1 if there is a connection between units i and j, otherwise to 
%         zero. Shape is set to 'sheet' by default. Shape does not have any 
%         meaning if a free topology is specified, anyway.
%       - msize must be a 1x2 vector
%
%  In cases 2...5 the sGrid structure is initiated by default values
%  which are set in SOM_SET. These include black markers 'o' (6pt),
%  light gray conncection lines (graph edges), unit coordinates
%  according to the lattice ('hexa','rect'), no labels, and no
%  surface.
%
%  OPTIONAL INPUT ARGUMENTS 
%   
%  Note: M is the number of map units.
%   
%  Here is a list of the valid arguments IDs and the associated
%  values (identifiers are case insensitive):
%
%  'Coord'     Unit coordinates
%              This defines the coordinates of the units. Default: the
%              topological coordinates (calculated as in function
%              SOM_VIS_COORDS and SOM_CPLANE). If the topology is free
%              (lattice is a connection matrix) this argument is obligatory!
%     (matrix) size Mx2 of 2D coordinates for each unit
%     (matrix) size Mx3 of 3D coordinates for each unit
% 
% 'Marker'       Unit markers, default is 'o'.
%     (string) 'o','+','x','*','v','^','<','>','h','s','d', 'p','.', or 'none'
%              give the same marker for each unit. 
%     (cell array) of size Mx1 of previous strings gives individual 
%              markers for each unit.
%     
% 'MarkerSize'   Size (pt) of unit markers, default is 6 (pt).
%     (scalar) gives the same size for every unit.  
%     (matrix) Mx1 gives an individual size for each unit marker.  
%
%  'MarkerColor' Unit marker colors, default is 'k'
%     (ColorSpec) gives the same color each unit. 
%     (matrix) Mx3 of RGB triples gives individual color for each unit
%              Note that indexed coloring - like in SOM_CPLANE - is
%              not possible. If indexed coloring is needed, you can
%              use SOM_NORMCOLOR to calculate RGB colors that
%              emulate indexed coloring. However, the colors for the
%              units are fixed, so changing colormap will not
%              change the colors.
%  
%  'Line'        Line type, default is '-'.
%     (string) '-',':','--' or '-.' or 'none'. Only one linetype in
%              grid is allowed.  
%
%  'LineWidth'   Width of the topological connection lines (edges)  
%     (scalar) gives the same width for each line. Default is 0.5.
%     (matrix) MxM sparse (or full) matrix gives individual width for 
%              each connection. The element (i,j), i<j, gives the line width for 
%              connection between nodes i and j. (The sparse form is
%              recommended for saving memory, a full matrix works as well,
%              of course). Note that only the elements satisfying i<j
%              matter - as the elememts for which j >= i are ignored in
%              order to avoid ambiguous situations if the matrix would be 
%              non-symmetric. The "connections to oneself" is not drawn. 
%
%              Line width zero is valid and causes the line to disappear.
%           
%  'LineColor'   Color of connection lines, default is [0.9 0.9 0.9].
%     (ColorSpec) gives the same color for each line
%     (matrix) MxMx3 matrix of RGB triples gives individual width for 
%              each connection. The element (i,j,:), i<j, gives the RGB triple for 
%              line between nodes i and j.     
%     (cell array) of three sparse (or full) matrices {r,g,b} where 
%              r(i,j), g(i,j) and b(i,j) gives the R,G, and B values in the RGB
%              triple for the line between nodes i and j. (The motivation for this
%              form is the fact that a 3D arrays can't use sparse format in 
%              Matlab version 5.1.)
%  
%              Note that only the elements satisfying i<j matter - the elememts 
%              for which j >= i are ignored in order to avoid ambiguous situations 
%              if the matrix was non-symmetric. The "connections to oneself"
%              is not drawn. 
%    
%
%  'Label'       Labels for units, default is [].
%     (empty)  [] means no labels.
%     (char array) of size Mx1. Element (i,:) has the label for unit i.     
%     (cell array) of size MxL consisting of sets of labels. Element {i,:} 
%              contains the labeling for unit i. 
%               In case of multiple labels, the labels for one unit are shown 
%               in one column centered at that unit.
%    
%   'LabelSize'   Text size of labels (points), default is 10.
%     (scalar) Default is 10.
%
%   'LabelColor'  Color of labels, default is 'c' (cyan).
%     (ColorSpec) gives the same color for each label string 'xor'
%                 sets the colors automatically so that they differ
%                 from the background (using Matlab's built-in xor-color feature.)
%    
%  'Surf'         Surface between nodes, default is [].
%     (empty)  [] gives no surface
%     (vector) Mx1 gives an indexed interpolated color surface between 
%              units using the actual colormap.
%     (matrix) Mx3 matrix of RGB triples gives a interpolated color surface 
%              between units using fixed RGB colors.
%    
%              Note that the interpolation is done using Matlab's built-in
%              color interpolation for SURF objects.
%
% OUTPUT ARGUMENTS
%
%  sGrid    (struct) with fields S.msize, S.shape, S.lattice, S.coord, S.marker, 
%                    S.markersize, S.markercolor, S.line, S.linewidth, S.linecolor,
%                    S.surf, S.label, S.labelsize, S.labelcolor
%
%  m        (matrix) handels to LINE objects (unit markers) 
%
%  l        (matrix) handles to LINE objects (lines connecting the units)
% 
%  t        (matrix) handles to TEXT objects (labels)
%
%  s        (scalar) handle  to SURF object  (surface between units)
%
% EXAMPLES
%
% % Make map of size 15x10 on random data:
% 
%    map=som_make(rand(1000,4),'msize',[15 10], 'lattice', 'hexa');
%
% % Draw the grid using two frist varable values as coordinates
% % and store the sGrid struct in varable S:
%
%    S=som_grid(map, 'coord', map.codebook(:,[1 2]))
%
% %Define some things: 
% %
% % Create a cell array of size 150x1 that divides map in to two label classes
% % 'circles' and 'squares'
%   
%    L(1:75,1)='o'; L(76:150,1)='s'; L = cellstr(L);
%    
% % Create a coloring according to the 3rd variable according to current
% % colormap: 
%    
%    C = som_normcolor(map.codebook(:,3));
% 
% % Change the visualization: use black lines, unit markers in M and unit
% % color in C, and set unit size to 10:
%
%    S=som_grid(S, 'linecolor', 'k', 'marker', L, 'MarkerColor',C, ...
%     'MarkerSize', 10);
%
% % Do a new visualization, use indexed color surface calcualted from the
% % first variable, coordinates according to the lattice (default) but 
% % no markers nor lines:
% 
%    S=som_grid(map,'line','none','marker','none','surf',map.codebook(:,1));
%
% % Set coordinates according to three last varables
%
%    som_grid(S,'coord',map.codebook(:,2:4));
%
% % Create a random connection matrix R1 and the usual hexagonal
% % neighborhood connection matrix R2: 
% 
%    R1=sparse(rand(150,150)>0.9); 
%    R2=som_connection(map);
%
% % Show these connections. Note that coordinates _must_ now be given
% % explicitly: we form default topological coordinates using
% % som_unit_coords.
%
%    som_grid(R1,map.topol.msize,'coord',som_unit_coords(map));
%    som_grid(R2,map.topol.msize,'coord',som_unit_coords(map));
%    
% % Show connections (R1 AND R2)
%    som_grid(R2.*R2,map.topol.msize,'coord',som_unit_coords(map));
% 
% OBJECT TAGS
%
%  No tags are set.
%
% SEE ALSO
%
%  som_show        The basic map visualization routine
%  som_cplane      The basic component plane visualization
%  som_connection  The basic topological connections
%  scatter         Scatter plots
%  scatter3        3-dimensional scatter plots 

% Copyright (c) 1999-2000 by the SOM toolbox programming team.
% http://www.cis.hut.fi/projects/somtoolbox/             

% Version 2.0beta Johan 061099 juuso 151199 310300

%% Init %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

True=1; False=0;                % const.
m=[]; l=[]; t=[]; s=[];         % default values for outputs
Ref=som_set('som_grid');        % reference struct

num_of_args=length(varargin);   % numb. of varargins

if num_of_args==0,
  S=som_set('som_grid');
  return;
end

switch class(varargin{1})
case 'struct'
  S=varargin{1};
  first_identifier=2;
  if ~isfield(S,'type'),
    error('Input struct is invalid: field ''type'' is missing.');
  end
  switch  S.type
  case 'som_grid'
    S=varargin{1};
    first_identifier=2;
  case 'som_map'
    Ref.lattice=S.topol.lattice;
    Ref.msize=S.topol.msize;
    Ref.shape=S.topol.shape;
    S=Ref;
    first_identifier=2;
  case 'som_topol'
    Ref.lattice=S.lattice;
    Ref.msize=S.msize;
    Ref.shape=S.shape;
    S=Ref;
    first_identifier=2;
  otherwise
    error('Input struct has to be of type som_grid, som_map or som_topol.');
  end
case 'cell'
  S=varargin{1};
  first_identifier=2;
  if vis_valuetype(S,{'topol_cell_no_shape'}),
    Ref.lattice=S{1};
    Ref.msize=S{2};
  elseif vis_valuetype(S,{'topol_cell'}),
    Ref.lattice=S{1};
    Ref.msize=S{2};
    Ref.shape=S{3};
  else
    error(['The cell value for 1st argument has to be {lattice, msize}' ...
	  'or {lattice, msize, shape}.']); 
  end
  S=Ref; 
case{'double','sparse','char'} 
  % Set defaults
  S=Ref;
  first_identifier=3; 
  if num_of_args<2,
    error('Not enough input arguments.');
  end
  S.lattice=varargin{1};
  S.msize=varargin{2};
otherwise
  error('Invalid input arguments!');
end  

%% Check input args %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=first_identifier:2:num_of_args,
  if ischar(varargin{i}) & isfield(Ref,lower(varargin{i})),
    if i+1>num_of_args,
      error('Invalid identifier-value pairs or wrong argument order.');
    else
      S=setfield(S,lower(varargin{i}),varargin{i+1});
    end
  elseif ischar(varargin{i}), 
    error(['Identifier ''' varargin{i} ''' is unknown.']);
  else
    error('Invalid identifier-value pairs or wrong argument order.');
  end
end 

% msize

if ~vis_valuetype(S.msize,{'1x2'}),
  error('msize has to be a 1x2 vector.');
end
munits=prod(S.msize);

% Default coordinates according to negihborhood

if isempty(S.coord),
  if ischar(S.lattice),
    switch S.lattice,
    case{'hexa','rect'}
      S.coord=som_vis_coords(S.lattice,S.msize);
    otherwise
      error('String value for lattice must be ''hexa'' or ''rect''.');
    end
  else
    error('Lattice is not ''hexa'' or ''rect'': coordinates must be given.');
  end
end

% connections

type=class(S.lattice);
switch type
case {'sparse','double'}   % free topology
  fixedline=False;   
case 'char'                % default topologies (hexa,char)
  switch S.lattice
  case 'hexa'
    hexa=True;
  case 'rect'
    hexa=False;
  otherwise
    error('Unknown lattice or neighborhood.');
  end

  % If topology is hexa/rect but linetype, color etc. is 
  % not constant, the topology is set to free

  if size(S.linewidth,1)>1 | size(S.linecolor,1)>1 | ...
    iscell(S.linecolor) % matrix or cell = not constant 
    fixedline=False;
    S.lattice=som_connection({S.lattice,S.msize,S.shape});
  else
    fixedline=True;
  end
end

% Check coordinate matrix size and set dummy zeros to z-axis
% if 2D coordinates (always 3D plots!)

if ~vis_valuetype(S.coord,{[munits 2],[munits 3]}),
   error('Coordinate matrix has wrong size.');
elseif size(S.coord,2)==2,
   S.coord(:,3)=0;
end

% Fixed marker size, color, type? 

if size(S.markersize,1)>1 | size(S.markercolor,1)>1 | size(S.marker,1)>1
   fixedmarker=False;
else
   fixedmarker=True;
end

% Check labels

if ~vis_valuetype(S.label,{'chararray','2Dcellarray_of_char'}) ...
      & ~isempty(S.label),
  error('Labels should be in a char array or cell array of strings.');
elseif ischar(S.label)
  S.label=cellstr(S.label);
end

if size(S.label,1) ~= munits & ~isempty(S.label),
  error('Number of labels and map size do not match.');
end

% Check line width, marker size, marker color,
% label size label color and surf sizes&types:

if ~vis_valuetype(S.linewidth,{[munits munits] [1 1]}),
  error('LineWidth matrix value has wrong size or dimension.');
elseif any(S.linewidth(:)<0),
  error('All elements of LineWidth must be non-negative.');
elseif ~vis_valuetype(S.markersize,{[munits 1] [1 1]}), 
  error('MarkerSize matrix value has wrong size or dimension.');
elseif any(S.markersize(:)<0), 
  error('All elements of MarkerSize must be non-negative.');
elseif ~vis_valuetype(S.markercolor,{'1x3rgb','colorstyle'}) & ...
      ~vis_valuetype(S.markercolor,{[munits 3],'nx3rgb'},'all'),
  error('MarkerColor should be a ColorSpec or Mx3 matrix of RGB triples.');
elseif ~vis_valuetype(S.labelcolor,{'1x3rgb','colorstyle','xor'}),
  error('LabelColor shoud be a ColorSpec or ''xor'' or ''none''.')
elseif ~vis_valuetype(S.labelsize,{'1x1'})
  error('LabelSize should be a scalar.');
elseif ~isempty(S.surf) & ~vis_valuetype(S.surf,{[munits 1] [munits 3]});
  error('Surf matrix value has wrong size or dimension.');
end

% Check marker type & size

if vis_valuetype(S.marker,{'cellcolumn_of_char'}) 
  % Don't bother to check the mareker strings in this case
  % let the plot3 handle them; it returns quite understandable
  % error messages, anyway
  
  if ~size(S.marker) == [munits 1],
    error(['Marker should be one of Matlab''s valid marker type,' ...
	   ' string ''none'' or a Mx1 cell array of these.']); 
  end
elseif ~vis_valuetype(S.marker,{'markerstyle','none'}),
      error(['Marker should be one of Matlab''s valid marker type,' ...
	   ' string ''none'' or a Mx1 cell array of these.']); 
end

% Check line type & size: only one line style allowed

if ~vis_valuetype(S.line,{'linestyle','none'}) 
  error(['Line should be a valid Matlab''s line style string or' ...
	 ' string ''none''.']);
end
	
% Check line color

if iscell(S.linecolor),
  if ndims(S.linecolor) ~= 2 | any(size(S.linecolor) ~= [1 3]),
    error('Cell input for LineColor should be of form {r,g,b}.')
  elseif ~vis_valuetype(S.linecolor{1},{[munits munits],'nxn[0,1]'},'all')| ...
	~vis_valuetype(S.linecolor{2},{[munits munits],'nxn[0,1]'},'all')| ...
	~vis_valuetype(S.linecolor{3},{[munits munits],'nxn[0,1]'},'all'),
    error(['In cell input {r,g,b} some matrix r,g or b is invalid: ' ...
	   'Size must be MxM and values in interval [0,1].']);
  end 
elseif ~vis_valuetype(S.linecolor,{'colorstyle','1x3rgb'}) & ...
      ~vis_valuetype(S.linecolor,{'nxnx3rgb', [munits munits 3]},'all'),
  error('Invalid LineColor: see help text for valid values.'),
elseif vis_valuetype(S.linecolor, {'none'}),
  error('LineColor ''none'' not allowed: set Line to ''none'' instead.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Action

memhold=ishold; % take hold state
if ~memhold
   cla;
end
hold on;

% Set surf if it exist

if ~isempty(S.surf),
   for i=1:3,
      s(:,:,i)=reshape(S.coord(:,i),S.msize);
   end
   s(:,:,4:3+size(S.surf,2))=reshape(S.surf,[S.msize size(S.surf,2)]);
   s=surf(s(:,:,1),s(:,:,2),s(:,:,3),s(:,:,4:end));
   set(s,'EdgeColor','none','Marker','none','FaceColor','interp');
end


if fixedline,
  % Line properties are fixed: draw fast, but
  % if line is set to 'none' set empty handle ans skip
  if strcmp(S.line,'none')
    l={};
  else
    p1=reshape(S.coord, [S.msize 3]);
    p2=zeros(size(p1)-[0 1 0]);
    p2(1:2:end,:,:)=p1(1:2:end,2:end,:);
    p2(2:2:end,:,:)=p1(2:2:end,1:end-1,:);
    
    l{1}=plot3(p1(:,:,1), p1(:,:,2), p1(:,:,3), ...
	       'Color', S.linecolor(1,:), ...
	       'LineWidth', S.linewidth(1), ...
	       'LineStyle', S.line);
    l{2}=plot3(p1(:,:,1)', p1(:,:,2)', p1(:,:,3)', ...
	       'Color', S.linecolor(1,:), ...
	       'LineWidth', S.linewidth(1), ...
	       'LineStyle', S.line);
    if hexa,
      l{3}=plot3(p2(:,:,1), p2(:,:,2), p2(:,:,3), ...
		 'Color', S.linecolor(1,:), ...
		 'LineWidth', S.linewidth(1), ...
		 'LineStyle', S.line);
    end
  end
  l=cat(1,l{:});
else
   % Variable properties: draw connection by connection
   
   [I,J,lw]=find(S.lattice); 
   x=[S.coord(I,1)'; S.coord(J,1)']; 
   y=[S.coord(I,2)'; S.coord(J,2)'];
   z=[S.coord(I,3)'; S.coord(J,3)'];
   if S.linewidth(1)==0,
      linewidth=0.5;
   else
      linewidth=S.linewidth(1);
   end
   if ndims(S.linecolor) ~=  3
     if isstr(S.linecolor)  
       l=plot3(x, y, z, ...
	       'Color', S.linecolor, ...
	       'LineWidth', linewidth, ...
	       'LineStyle',S.line);
     else 
       if iscell(S.linecolor)
         lcolor=[S.linecolor{1}(1,1) S.linecolor{2}(1,1) S.linecolor{3}(1,1)];
         l=plot3(x, y, z, ...
		 'Color', lcolor, ...
		 'LineWidth', linewidth, ...
		 'LineStyle',S.line);
       else
         l=plot3(x, y, z, ...
                 'Color', S.linecolor(1,:), ...
                 'LineWidth', linewidth, ...
                 'LineStyle',S.line);
       end
     end
   else
     l=plot3(x, y, z, ...
	     'Color', S.linecolor(1,1,:), ...
	     'LineWidth', linewidth, ...
	     'LineStyle',S.line);
   end
end

if fixedmarker,

  % If marker is set to 'none' skip and set empty handle 
  if strcmp(S.marker,'none')
    m=[]; 
  else
    % Fixed markers: draw all in one command

    m=plot3(S.coord(:,1), S.coord(:,2), S.coord(:,3), ... 
	    'LineStyle', 'none', ...
	    'Marker', S.marker, ...
	    'MarkerSize', S.markersize(1), ...
	    'MarkerFaceColor', S.markercolor(1,:), ...
	    'MarkerEdgeColor', S.markercolor(1,:));
  end
else
  % Variable marker properties: draw marker by marker

  x=[S.coord(:,1)'; S.coord(:,1)']; 
  y=[S.coord(:,2)'; S.coord(:,2)'];
  z=[S.coord(:,3)'; S.coord(:,3)'];
  if iscell(S.marker)
    marker=S.marker{1};
  else
    marker=S.marker(1);
  end
  sz=max(S.markersize(1),0.1);
  m=plot3(x, y, z, ... 
	  'LineStyle', 'none', ...
	  'Marker', marker, ...
	  'MarkerSize', sz, ... 
	  'MarkerFaceColor', S.markercolor(1,:), ...
	  'MarkerEdgeColor', S.markercolor(1,:));
end

L=length(l); 
n=munits;

%%% Set variable properties %%%

% Line width

if length(S.linewidth)>1 
   lwidth=diag(S.linewidth(I,J)); 

   % Handle zero width
   iszero=(lwidth == 0);lwidth(iszero)=0.5;
   for i=1:length(l),
     set(l(i),'LineWidth', lwidth(i));
   end
   if ~isempty(iszero), % zero width
      set(l(iszero),'Visible','off');
   end
end

% Line color

if size(S.linecolor,1)>1 | iscell(S.linecolor)
   if length(size(S.linecolor)) == 3 | iscell(S.linecolor) 
     if ~iscell(S.linecolor)
       for i=1:L
         set(l(i),'Color',S.linecolor(I(i),J(i),:));
       end
     else
       for i=1:L
         lcolor=[S.linecolor{1}(I(i),J(i)),...
                 S.linecolor{2}(I(i),J(i)),...
                 S.linecolor{3}(I(i),J(i))];
         set(l(i),'Color',lcolor);
       end
     end
   else
     for i=1:L,
       set(l(i),'Color', S.linecolor(I(i),:));
     end
   end
end

% Marker size

if length(S.markersize)>1
   % handle zero size
   iszero=find(~S.markersize);
   S.markersize(iszero)=1;
   for i=1:n,
      set(m(i),'MarkerSize', S.markersize(i));
   end
   if ~isempty(iszero), % zero size
      set(m(iszero),'Visible','off');
   end
end

% Marker type

if size(S.marker,1)>1
   S.marker=char(S.marker);
   for i=1:n,
      set(m(i),'Marker', S.marker(i));
   end
end

% Marker color

if size(S.markercolor,1)>1
   for i=1:n,
     set(m(i),'MarkerFaceColor', S.markercolor(i,:), ...
	      'MarkerEdgeColor', S.markercolor(i,:));
   end
end

% Set labels if they exist

if ~isempty(S.label)
  if vis_valuetype(S.labelcolor,{'xor'}),
    S.labelcolor='g';
    XOR=1;
  else
    XOR=0;
  end
  if vis_valuetype(S.labelcolor,{'none'}),
    S.labelcolor='g';
    VIS = 1;
  else
    VIS = 0;
  end
  for i=1:size(S.label,1),
    L=cat(1,S.label(i,:)); 
    for j=length(L):-1:1,
      if isempty(L{j}),
	L=L(1:end-1); 
      end
    end
    
    if isempty(L),
      L='';
    end
    t(i)=text(S.coord(i,1), S.coord(i,2), S.coord(i,3), L,...
	'FontSize', S.labelsize, 'Color',S.labelcolor, ...
	'HorizontalAlignment', 'center');
  end
  if XOR
    set(t,'EraseMode','xor');
  end
  if VIS
    set(t,'Visible','off');
  end 
else
   t=[];
end

%% Set hold state

if ~memhold,
   hold off; 
end

if nargout==0,
  clear S m l t s;
end