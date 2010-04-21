function ft_sourcewrite(cfg, volume)

% FT_VOLUMEWRITE exports source analysis results to an Analyze MRI file
% that can subsequently be read into BrainVoyager or MRIcro
%
% Warning: FT_NORMALISEVOLUME has been renamed to FT_VOLUMENORMALISE
% Warning: backward compatibility will be removed in the future

% Copyright (C) 2005-2006, F.C. Donders Centre
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

warning('SOURCEWRITE has been renamed to VOLUMEWRITE');
warning('backward compatibility will be removed in the future');

ft_volumewrite(cfg, volume);
