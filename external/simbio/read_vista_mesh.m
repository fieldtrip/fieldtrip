function [nodes,elements,labels] = read_vista_mesh(filename)

% Reads a vista format mesh using correspondently named compiled mex file.
%
% Use as
%   [nodes,elements,labels] = read_vista_mesh(filename);
% where
%   filename = the name of the Vista mesh (with extension .v)

error('The mex file %s is missing', [mfilename '.' mexext]);
