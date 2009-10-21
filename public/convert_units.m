function [obj] = convert_units(obj, target);

% CONVERT_UNITS changes the geometrical dimension to the specified SI unit.
% The units of the input object is determined from the structure field
% object.unit, or is estimated based on the spatial extend of the structure, 
% e.g. a volume conduction model of the head should be approximately 20 cm large.
%
% Use as
%   [object] = convert_units(object, target)
%
% The following input objects are supported
%   simple dipole position
%   electrode definition
%   gradiometer array definition
%   volume conductor definition
%   dipole grid definition
%   anatomical mri
%
% Possible target units are 'm', 'dm', 'cm ' or 'mm'.
%
% See READ_VOL, READ_SENS

% Copyright (C) 2005-2008, Robert Oostenveld
%
% $Log: convert_units.m,v $
% Revision 1.8  2009/03/26 14:57:04  roboos
% moved the unit estimation to a seperate function
%
% Revision 1.7  2009/03/11 11:28:24  roboos
% detect as mm for very wide anatomical MRIs (happens at the fcdc with the ctf MRIs)
%
% Revision 1.6  2009/02/05 10:20:35  roboos
% added bemcp as volume type
%
% Revision 1.5  2009/01/08 17:18:45  roboos
% fixed bug mom->pos for source structures
% added unit detection for volumes with a transform and a dim
%
% Revision 1.4  2008/07/21 20:29:22  roboos
% small change in output on screen
%
% Revision 1.3  2008/04/14 20:53:58  roboos
% added detection for headshape and/or fiducials
% fixed bug in scaling of fiducials
%
% Revision 1.2  2008/04/14 19:29:58  roboos
% cleanded up code and added autodetection based on geometrical size of object
% changed interface, no forced input type is possible
% changed from dos to unix
%
% Revision 1.1  2005/03/03 11:01:39  roboos
% already old (and unused) implementation, but sofar this function was not included in CVS
%

% This function consists of three parts:
%   1) determine the input units
%   2) determine the requested scaling factor to obtain the output units
%   3) try to apply the scaling to the known geometrical elements in the input object

% determine the unit-of-dimension of the input object
if isfield(obj, 'unit')
  % use the units specified in the object
  unit = obj.unit;
else
  % try to estimate the units from the object
  type = voltype(obj);
  if ~strcmp(type, 'unknown')
    switch type
      case 'infinite'
        % there is nothing to do to convert the units
        unit = target;

      case 'singlesphere'
        size = obj.r;

      case 'multisphere'
        size = median(obj.r);

      case 'concentric'
        size = max(obj.r);

      case 'nolte'
        size = norm(range(obj.bnd.pnt));

      case {'bem' 'dipoli' 'bemcp' 'asa' 'avo'}
        size = norm(range(obj.bnd(1).pnt));

      otherwise
        error('cannot determine geometrical units of volume conduction model');
    end % switch
    
    % determine the units by looking at the size
    unit = estimate_units(size);

  elseif senstype(obj, 'meg')
    size = norm(range(obj.pnt));
    unit = estimate_units(size);

  elseif senstype(obj, 'eeg')
    size = norm(range(obj.pnt));
    unit = estimate_units(size);

  elseif isfield(obj, 'pnt') && ~isempty(obj.pnt)
    size = norm(range(obj.pnt));
    unit = estimate_units(size);
    
  elseif isfield(obj, 'pos') && ~isempty(obj.pos)
    size = norm(range(obj.pos));
    unit = estimate_units(size);

  elseif isfield(obj, 'transform') && ~isempty(obj.transform)
    % construct the corner points of the voxel grid in head coordinates
    xi = 1:obj.dim(1);
    yi = 1:obj.dim(2);
    zi = 1:obj.dim(3);
    pos = [
      xi(  1) yi(  1) zi(  1)
      xi(  1) yi(  1) zi(end)
      xi(  1) yi(end) zi(  1)
      xi(  1) yi(end) zi(end)
      xi(end) yi(  1) zi(  1)
      xi(end) yi(  1) zi(end)
      xi(end) yi(end) zi(  1)
      xi(end) yi(end) zi(end)
      ];
    pos = warp_apply(obj.transform, pos);
    size = norm(range(pos));
    unit = estimate_units(size);
    
  elseif isfield(obj, 'fid') && isfield(obj.fid, 'pnt') && ~isempty(obj.fid.pnt)
    size = norm(range(obj.fid.pnt));
    unit = estimate_units(size);
    
  else
    error('cannot determine geometrical units');
    
  end % recognized type of volume conduction model or sensor array
end % determine input units

if nargin<2
  % just remember the units in the output and return
  obj.unit = unit;
  return
elseif strcmp(unit, target)
  % no conversion is needed
  obj.unit = unit;
  return
end

% give some information about the conversion
fprintf('converting units from ''%s'' to ''%s''\n', unit, target)

if strcmp(unit, 'm')
  unit2meter = 1;
elseif strcmp(unit, 'dm')
  unit2meter = 0.1;
elseif strcmp(unit, 'cm')
  unit2meter = 0.01;
elseif strcmp(unit, 'mm')
  unit2meter = 0.001;
end

% determine the unit-of-dimension of the output object
if strcmp(target, 'm')
  meter2target = 1;
elseif strcmp(target, 'dm')
  meter2target = 10;
elseif strcmp(target, 'cm')
  meter2target = 100;
elseif strcmp(target, 'mm')
  meter2target = 1000;
end

% combine the units into one scaling factor
scale = unit2meter * meter2target;

% volume conductor model
if isfield(obj, 'r'), obj.r = scale * obj.r; end
if isfield(obj, 'o'), obj.o = scale * obj.o; end
if isfield(obj, 'bnd'), for i=1:length(obj.bnd), obj.bnd(i).pnt = scale * obj.bnd(i).pnt; end, end

% gradiometer array
if isfield(obj, 'pnt1'), obj.pnt1 = scale * obj.pnt1; end
if isfield(obj, 'pnt2'), obj.pnt2 = scale * obj.pnt2; end
if isfield(obj, 'prj'),  obj.prj  = scale * obj.prj;  end

% gradiometer array, electrode array, head shape or dipole grid
if isfield(obj, 'pnt'), obj.pnt = scale * obj.pnt; end

% fiducials
if isfield(obj, 'fid') && isfield(obj.fid, 'pnt'), obj.fid.pnt = scale * obj.fid.pnt; end

% dipole grid
if isfield(obj, 'pos'), obj.pos = scale * obj.pos; end

% anatomical MRI or functional volume
if isfield(obj, 'transform'),
  H = diag([scale scale scale 1]);
  obj.transform = H * obj.transform;
end

% remember the unit
obj.unit = target;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%   Use as
%     r = range(x)
%   or you can also specify the dimension along which to look by
%     r = range(x, dim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = range(x, dim)
if nargin==1
  r = max(x) - min(x);
else
  r = max(x, [], dim) - min(x, [], dim);
end
