function hs = ft_plot_headshape(headshape,varargin)

% FT_PLOT_HEADSHAPE visualizes the shape of a head from a variety of
% acquisition system. Usually the head shape is measured with a
% Polhemus tracker and some proprietary software (e.g. from CTF, BTi
% or Yokogawa). The headshape and fiducials can be used for coregistration.
% If present in the headshape, the location of the fiducials will also
% be shown.
%
% Use as
%   hs = ft_plot_headshape(shape, ...)
% where the shape is a structure obtained from FT_READ_HEADSHAPE.
%
% Optional arguments should come in key-value pairs and can include
%   'vertexcolor' = color specification as [r g b] values or a string, for example 'brain', 'cortex', 'skin', 'red', 'r'
%   'vertexsize'  = scalar value specifying the size of the vertices (default = 10)
%   'fidcolor'    = color specification as [r g b] values or a string, for example 'brain', 'cortex', 'skin', 'red', 'r'
%   'fidmarker'   = ['.', '*', '+',  ...]
%   'fidlabel'    = ['yes', 'no', 1, 0, 'true', 'false']
%   'transform'   = transformation matrix for the fiducials, converts MRI
%                   voxels into head shape coordinates
%
% Example
%   shape = ft_read_headshape(filename);
%   ft_plot_headshape(shape)

% Copyright (C) 2009, Cristiano Micheli
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

ws = warning('on', 'MATLAB:divideByZero');

if ~isstruct(headshape) && isnumeric(headshape) && size(headshape,2)==3
  % the input seems like a list of points, convert into something that resembles a headshape
  warning('off', 'MATLAB:warn_r14_stucture_assignment');
  headshape.pnt = headshape;
end

% get the optional input arguments
vertexcolor = ft_getopt(varargin, 'vertexcolor',  'r');
vertexsize  = ft_getopt(varargin, 'vertexsize',   10);
facecolor   = ft_getopt(varargin, 'facecolor',    'none');
edgecolor   = ft_getopt(varargin, 'edgecolor',    'none');
fidcolor    = ft_getopt(varargin, 'fidcolor',     'g');
fidmarker   = ft_getopt(varargin, 'fidmarker',    '*');
fidlabel    = ft_getopt(varargin, 'fidlabel',     true);
transform   = ft_getopt(varargin, 'transform');

% start with empty return values
hs      = [];

% everything is added to the current figure
holdflag = ishold;
if ~holdflag
  hold on
end

pnt = headshape.pnt;
bnd.pnt = pnt;
bnd.tri = [];

ft_plot_mesh(bnd, 'vertexcolor',vertexcolor,'vertexsize',vertexsize);

if isfield(headshape, 'fid')
  fid = headshape.fid;
  if ~isempty(transform)
    % spatially transform the fiducials
    % FIXME what is the reason for this?
    fid.pnt = warp_apply(transform, fid.pnt);
  end
  
  % show the fiducial labels
  for i=1:size(fid.pnt,1)
    hs = plot3(fid.pnt(i,1), fid.pnt(i,2), fid.pnt(i,3), 'Marker',fidmarker,'MarkerEdgeColor',fidcolor);
    if isfield(fid,'label') && istrue(fidlabel)
      str = sprintf('%s', fid.label{i});
      h   = text(fid.pnt(i, 1), fid.pnt(i, 2), fid.pnt(i, 3), str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','Interpreter','none');
      hs  = [hs; h];
    end
  end
end

if nargout==0
  clear hs
end
if ~holdflag
  hold off
end

warning(ws); %revert to original state

