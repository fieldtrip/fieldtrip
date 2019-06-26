function headmodel = ft_headmodel_asa(filename)

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
%   headmodel = ft_headmodel_asa(filename)
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% Copyright (C) 2012, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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

% read the headmodel from file
% this works for ASA version 2.x, perhaps also for ASA 3.x 
headmodel = ft_read_headmodel(filename);

