function [source] = ft_appendsource(cfg, varargin)

% FT_APPENDSOURCE concatenates multiple volumetric source reconstruction
% data structures that have been processed seperately. 
%
% If the source reconstructions were computed for different ROIs or
% different slabs of a regular 3D grid (as indicated by the source
% positions), the data will be concatenated along the spatial dimension.
%
% If the source reconstructions were computed on the same source
% positions, but for different frequencies and/or latencies, e.g. for
% time-frequency spectrally decomposed data, the data will be concatenared
% along the frequency and/or time dimension.
%
% Use as
%   combined = ft_appendsource(cfg, source1, source2, ...)
%
% See also FT_SOURCEANALYSIS, FT_APPENDDATA, FT_APPENDFREQ, FT_APPENDSOURCE

% Copyright (C) 2011, Robert Oostenveld
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble provenance
ft_preamble trackconfig
ft_preamble debug
ft_preamble loadvar varargin

% check if the input data is valid for this function
for i=1:length(varargin)
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'source', 'feedback', 'yes', 'hassampleinfo', 'ifmakessense');
end

% use a helper function to select the consistent parts of the data and to concatenate it
source = ft_selectdata(varargin{:}, 'param', {'pow'});

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
ft_postamble previous varargin
ft_postamble history source
ft_postamble savevar source

