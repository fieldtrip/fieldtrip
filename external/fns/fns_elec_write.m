function fns_elec_write(pnt,vsize,dimpos,elecfile)
% 
% This function write the electrodes locations on disk and other
% information
% The rgn structure contains the following fields:
%   locations is the location of the electrodes
%   gridlocs is the location of the electrodes again
%   voxel_sizes is the dimension of a voxel in the cartesian coordinates in mm [3X1]
%   nodes_sizes is the dimension of the dipoles grid [NX NY NZ]
%   status is the electrodes status
%   values is spare
%   info is spare

% $Copyright (C) 2010 by Hung Dang$

% initialize a region structure rgn
rgn = struct('locations',pnt, ...
             'gridlocs',pnt, ...
             'voxel_sizes',vsize, ... 
             'node_sizes',dimpos, ...
             'status',ones(size(pnt,1),1), ...
             'values',[],'info',[]);

fns_region_write(elecfile,rgn);