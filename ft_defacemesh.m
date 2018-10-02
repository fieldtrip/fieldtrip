function mesh = ft_defacemesh(cfg, mesh)

% FT_DEFACEMESH allows you to de-identify a scalp surface mesh by erasing specific
% regions, such as the face and ears. The graphical user interface allows you to
% position a box over the anatomical data inside which all vertices will be removed.
% You might have to call this function multiple times when both face and ears need to
% be removed. Following defacing, you should check the result with FT_PLOT_MESH.
%
% Use as
%   mesh = ft_defacevolume(cfg, mesh)
%
% The configuration can contain the following options
%   cfg.translate  = initial position of the center of the box (default = [0 0 0])
%   cfg.scale      = initial size of the box along each dimension (default is automatic)
%   cfg.translate  = initial rotation of the box (default = [0 0 0])
%   cfg.selection  = which voxels to keep, can be 'inside' or 'outside' (default = 'outside')
%
% See also FT_ANONIMIZEDATA, FT_DEFACEVCOLUME, FT_ANALYSISPIPELINE, FT_PLOT_MESH

% Copyright (C) 2015-2016, Robert Oostenveld
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    mesh
ft_preamble provenance mesh
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% the actual work is done by FT_DEFACEVOLUME
tmpcfg = cfg;
tmpcfg.showcallinfo = 'no';
mesh = ft_defacevolume(tmpcfg, mesh);
% restore provenance information and put back cfg.callinfo
tmpcallinfo = cfg.showcallinfo;
[cfg, mesh] = rollback_provenance(cfg, mesh);
cfg.showcallinfo = tmpcallinfo;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous mesh
ft_postamble provenance mesh
ft_postamble history mesh
ft_postamble savevar mesh
