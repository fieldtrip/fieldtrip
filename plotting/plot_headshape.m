function hs = plot_headshape(headshape,varargin)

% PLOT_HEADSHAPE visualizes the shape of a head generated from a variety of files 
% (like CTF and Polhemus). The headshape and fiducials can for example be used for coregistration.
%
% Use as
%   hs = plot_headshape(shape, varargin)
%
% Graphic facilities are available for vertices and fiducials (Nasion, Left, Right ...). A list of
% the arguments is given below with the correspondent admitted choices.
%
%     'vertexcolor'   ['brain', 'cortex', 'skin', 'black', 'red', 'r', ..., [0.5 1 0], ...]
%     'fidcolor'      ['brain', 'cortex', 'skin', 'black', 'red', 'r', ..., [0.5 1 0], ...]
%     'fidmarker'     ['.', '*', '+',  ...]
%     'fidlabel'      ['yes', 'no', 1, 0, 'true', 'false']
%     'transform'     transformation matrix for the fiducials, converts MRI
%                       voxels into head shape coordinates
%
% Example
%  [shape] = read_headshape(filename);
%   plot_headshape(shape)

% Copyright (C) 2009, Cristiano Micheli
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

warning('on', 'MATLAB:divideByZero');

% get the optional input arguments
vertexcolor = keyval('vertexcolor', varargin); if isempty(vertexcolor), vertexcolor='r'; end
fidcolor    = keyval('fidcolor',    varargin); if isempty(fidcolor), fidcolor='g'; end
fidmarker   = keyval('fidmarker',   varargin); if isempty(fidmarker), fidmarker='.'; end
fidlabel    = keyval('fidlabel',    varargin); if isempty(fidlabel), fidlabel='no'; end
transform   = keyval('transform',    varargin); if isempty(transform), transform=[]; end

% start with empty return values
hs      = [];

% everything is added to the current figure
holdflag = ishold;
hold on

pnt = headshape.pnt;
bnd.pnt = pnt;
bnd.tri = [];

hs  = plot_mesh(bnd, 'vertices', 'yes', 'vertexcolor',vertexcolor,'vertexsize',10);

if isfield(headshape, 'fid')
  fid = headshape.fid;
  if ~isempty(transform)
    % plot the fiducials
    fidc = fid.pnt;
    try
      fidc = warp_apply(transform, fidc);
    end
    hs = plot3(fidc(:,1), fidc(:,2), fidc(:,3), 'Marker',fidmarker,'MarkerEdgeColor',fidcolor);
    % show the fiducial labels
    if isfield(fid,'label') && istrue(fidlabel)
      for node_indx=1:size(fidc,1)
        str = sprintf('%s', fid.label{node_indx});
        h   = text(fidc(node_indx, 1), fidc(node_indx, 2), fidc(node_indx, 3), str, ...
                  'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','Interpreter','none');
        hs  = [hs; h];
      end
    end    
  end
   
end
if nargout==0
  clear hs
end
if ~holdflag
  hold off
end
