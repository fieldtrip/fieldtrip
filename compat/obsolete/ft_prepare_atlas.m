function [atlas, cfg] = ft_prepare_atlas(varargin)

% FT_PREPARE_ATLAS reads in a specified atlas with coordinates and
% anatomical labels. It either uses the AFNI brik file that is available
% from http://afni.nimh.nih.gov/afni/doc/misc/ttatlas_tlrc, or it
% uses one of the WFU atlasses available from http://fmri.wfubmc.edu.
%
% This function is called by other FieldTrip functions that make
% use of an atlas, for example for plotting or for selection of an
% anatomical region of interest.
%
% Use as
%   [atlas] = ft_prepare_atlas(cfg)
%
% where the configuration should contain
%   cfg.atlas      string, filename of the atlas to use
%
% See also FT_VOLUMELOOKUP, FT_SOURCEPLOT, FT_SOURCESTATISTICS

% Copyright (C) 2005-2011, Robert Oostenveld, Ingrid Nieuwenhuis
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

warning('ft_prepare_atlas is only a compatibility wrapper, which will soon be removed. Please instead call ft_read_atlas.');

% assume first input to be a cfg structure, containing the atlas filename
filename = varargin{1}.atlas;
atlas    = ft_read_atlas(filename);

