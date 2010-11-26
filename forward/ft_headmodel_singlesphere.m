function vol = ft_headmodel_singlesphere(pnt, varargin)

% FT_HEADMODEL_SINGLESPHERE fits a spherical model to a set of points that
% describe the head surface.
%
% Use as
%   vol = ft_headmodel_singlesphere(pnt, ...)
%
% Optional input arguments should be specified in key-value pairs and can
% include
%   headshape        = string, filename with headshape
%   conductivity     = number, conductivity of the sphere

% get the optional arguments
headshape    = keyval('headshape', varargin);
conductivity = keyval('conductivity', varargin);

% start with an empty volume conductor
vol = [];

if ~isempty(headshape)
  headshape = ft_read_headshape(headshape);
  % use the headshape points instead of the user-specified ones
  pnt = headshape.pnt;
  if isfield(headshape, 'unit')
    % copy the geometrical units into he volume conductor
    vol.unit = headshape.unit;
  end
end

% fit a single sphere to all headshape points
[single_o, single_r] = fitsphere(pnt);

vol.r = single_r;
vol.o = single_o;
vol.c = conductivity;
vol.type = 'singlesphere';
