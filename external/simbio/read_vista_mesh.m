function [varargout] = read_vista_mesh(varargin)

% READ_VISTA_MESH is implemented as mex file
%
% Use as
%   [nodes,elements,labels] = read_vista_mesh(filename);
% where
%   filename = the name of the Vista mesh (with extension .v)
%
% $Id$

error('The mex file %s is missing', [mfilename '.' mexext]);
