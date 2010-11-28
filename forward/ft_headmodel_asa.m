function vol = ft_headmodel_asa(filename, varargin)

% FT_HEADMODEL_ASA reads a volume conduction model from an ASA *.vol
% file
%
% ASA is commercial software (http://www.ant-neuro.com) that supports
% among others the boundary element method (BEM) for EEG. This function
% allows you to read an EEG BEM volume conduction model from an ASA
% format file (*.vol) and use that for leadfield computations in
% MATLAB. Constructing the geometry of the head model from an anatomical
% MRI and the computation of the BEM system are both handled by ASA.
% 
% Use as
%   vol = ft_headmodel_asa(filename)
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% read the headmodel from file
% this works for ASA version 2.x, perhaps also for ASA 3.x 
vol = ft_read_vol(filename);

