function [inside] = inside_vol(pos, vol)

% INSIDE_VOL locates dipole locations inside/outside the source
% compartment of a volume conductor model.
%
% [inside] = inside_vol(pos, vol, ...)
%
% where the input should be
%   pos      Nx3 matrix with dipole positions
%   vol      structure with volume conductor model
% and the output is
%   inside   list of dipoles inside the brain compartment
%            (1=inside, 0=outisde)
%
% Additional optional input arguments shoudl be given in key value pairs
% and can include
%   <none>

% Copyright (C) 2003-2007, Robert Oostenveld
%
% $Log: inside_vol.m,v $
% Revision 1.5  2009/02/05 10:20:35  roboos
% added bemcp as volume type
%
% Revision 1.4  2009/01/21 11:15:57  roboos
% moved function back from forwinv/private into public section, because it is part of the public API
%
% Revision 1.1  2009/01/21 10:46:10  roboos
% moved from forwinv/* and forwinv/mex/* directory to forwinv/private/* to make the CVS layout consistent with the release version
%
% Revision 1.2  2008/09/29 12:04:41  roboos
% use logical (built-in) instead of boolean (simulink)
%
% Revision 1.1  2008/09/20 13:41:35  roboos
% moved content of find_inside_vol to new inside_vol function with slightly different interface
% added wrapper for spm
%
%
% %%%
% Switch from find_inside_vol.m to inside_vol.m
% 2008/09/09 chrisp
% %%%
%
% Revision 1.7  2008/04/15 20:36:21  roboos
% added explicit handling of various BEM implementations, i.e. for all
% voltype variants
%
% Revision 1.6  2007/07/25 08:34:05  roboos
% switched to using voltype helper function
% also support N-shell concentric sphere model
% changed detection for single and concentric sphere model, now explicitely
% using distance and source compartment
%
% Revision 1.5  2004/09/06 08:46:27  roboos
% moved reading of neuromag BEM boundary into fieldtrip's prepare_vol_sens
%
% Revision 1.4  2004/09/03 09:07:17  roboos
% added support for finding dipoles in neuromag BEM model
%
% Revision 1.3  2003/08/04 09:19:30  roberto
% fixed bug for dipole on BEM volume boundary
%
% Revision 1.2  2003/03/24 12:30:06  roberto
% added support for multi-sphere volume conductor model
%

% determine the type of volume conduction model
switch voltype(vol)

  % single-sphere or multiple concentric spheres
  case {'singlesphere' 'concentric'}
    if ~isfield(vol, 'source')
      % locate the innermost compartment and remember it
      [dum, vol.source] = min(vol.r);
    end
    if isfield(vol, 'o')
      % shift dipole positions toward origin of sphere
      tmp = pos - repmat(vol.o, size(pos,1), 1);
    else
      tmp = pos;
    end
    distance = sqrt(sum(tmp.^2, 2))-vol.r(vol.source);
    % positive if outside, negative if inside
    inside   = distance<0;

    % multi-sphere volume conductor model
  case 'multisphere'

    % nspheres = size(vol.r,1);
    % ndipoles = size(pos,1);
    % inside = zeros(ndipoles,1);
    % for sph=1:nspheres
    % for dip=1:ndipoles
    %   if inside(dip)
    %     % the dipole has already been detected in one of the other spheres
    %     continue
    %   end
    %   inside(dip) = (norm(pos(dip,:) - vol.o(sph,:)) <= vol.r(sph));
    % end
    % end
    % outside = find(inside==0);
    % inside  = find(inside==1);

    % this is a much faster implementation
    nspheres = size(vol.r,1);
    ndipoles = size(pos,1);
    inside = zeros(ndipoles,1);
    for sph=1:nspheres
      % temporary shift dipole positions toward origin
      if isfield(vol, 'o')
        tmp = pos - repmat(vol.o(sph,:), [ndipoles 1]);
      else
        tmp = pos;
      end
      flag = (sqrt(sum(tmp.^2,2)) <= vol.r(sph));
      inside = inside + flag;
    end
    inside  = inside>0;

    % realistic BEM volume conductor model
  case {'bem', 'dipoli', 'bemcp', 'asa', 'avo', 'nolte', 'neuromag'}
    if ~isfield(vol, 'source')
      % locate the innermost compartment and remember it
      vol.source = find_innermost_boundary(vol.bnd);
    end
    % use the specified source compartment
    pnt = vol.bnd(vol.source).pnt;
    tri = vol.bnd(vol.source).tri;
    % determine the dipole positions that are inside the brain compartment
    inside  = bounding_mesh(pos, pnt, tri);

    % unrecognized volume conductor model
  otherwise
    error('unrecognized volume conductor model');
end

% ensure that these are column vectors
inside  = logical(inside(:));
